import argparse
import json
import struct
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np

STEP  = 13.517831473
MAGIC = 0x554C414D
LAM   = 2.57144749847630
WAVE  = 2.44344296778474
CHUNK = 2_000_000

def _meta(path: str) -> tuple[int, int]:
    with open(path, "rb") as f:
        magic, _ = struct.unpack("<II", f.read(8))
        if magic != MAGIC:
            raise ValueError(f"Not a valid saved.bin: {path!r}")
        max_computed, count = struct.unpack("<QQ", f.read(16))
    return int(count), int(max_computed)


def _ulam_chunks(path: str, limit: int | None = None):
    with open(path, "rb") as f:
        magic, _ = struct.unpack("<II", f.read(8))
        if magic != MAGIC:
            raise ValueError(f"Not a valid saved.bin: {path!r}")
        _, total = struct.unpack("<QQ", f.read(16))
        if limit is not None:
            total = min(int(total), int(limit))
        pos = 0
        while pos < total:
            n     = min(CHUNK, total - pos)
            a_off = np.frombuffer(f.read(n * 4), dtype=np.int32)
            idx   = np.arange(pos + 1, pos + 1 + n, dtype=np.float64)
            yield a_off.astype(np.int64) + (idx * STEP).astype(np.int64)
            pos  += n


def _load_ulams_up_to(path: str, max_val: int) -> np.ndarray:
    parts: list[np.ndarray] = []
    for chunk in _ulam_chunks(path):
        mask = chunk <= max_val
        if mask.any():
            parts.append(chunk[mask])
        if not mask.all():
            break
    return np.concatenate(parts) if parts else np.array([], dtype=np.int64)


def _progress(iterable, total: int, desc: str):
    try:
        from tqdm import tqdm
        return tqdm(iterable, total=total, desc=desc)
    except ImportError:
        return iterable

def cmd_gaps(args: argparse.Namespace) -> None:
    from_json: str | None = getattr(args, "from_json", None)

    if from_json:
        with open(from_json) as fj:
            raw = json.load(fj)
        gaps: Counter[int] = Counter({int(k): int(v) for k, v in raw.items()})
        total_pairs = sum(gaps.values())
        subtitle = f"source: {Path(from_json).name}"
        print(f"Loaded {len(gaps):,} distinct gap values from {from_json!r}")
    else:
        count, max_computed = _meta(args.state)
        print(f"Computing gaps for {count:,} Ulam numbers  (max={max_computed:,}) ...")
        gaps = Counter()
        prev: int | None = None
        n_chunks = (count + CHUNK - 1) // CHUNK
        for chunk in _progress(_ulam_chunks(args.state), n_chunks, "Reading"):
            if prev is not None:
                gaps[int(chunk[0]) - prev] += 1
            diffs = np.diff(chunk)
            vals, cnts = np.unique(diffs, return_counts=True)
            for v, c in zip(vals.tolist(), cnts.tolist()):
                gaps[int(v)] += int(c)
            prev = int(chunk[-1])
        total_pairs = sum(gaps.values())
        subtitle = f"{count:,} Ulam numbers  (max = {max_computed:,})"

        json_path = args.json or str(Path(args.output).with_suffix(".json"))
        with open(json_path, "w") as fj:
            json.dump({str(k): v for k, v in sorted(gaps.items())}, fj, indent=2)
        print(f"Gap counts  ->  {json_path!r}")

    sorted_gaps = sorted(gaps.items())
    xs_all = [g for g, _ in sorted_gaps]
    ys_all = [c for _, c in sorted_gaps]

    SPLIT1, SPLIT2 = 50, 300

    def _band(lo, hi):
        pairs = [(x, y) for x, y in sorted_gaps if lo <= x <= hi]
        if not pairs:
            return [], []
        return zip(*pairs)  # type: ignore[return-value]

    xs1, ys1 = _band(1, SPLIT1)
    xs2, ys2 = _band(SPLIT1 + 1, SPLIT2)

    fig, axes = plt.subplots(
        3, 1, figsize=(16, 14),
        gridspec_kw={"height_ratios": [3, 2.2, 1.5]},
    )
    fig.suptitle(
        f"Consecutive-gap distribution -- {subtitle}\n"
        f"Total pairs: {total_pairs:,}   Distinct gaps: {len(gaps):,}",
        fontsize=13, y=0.99,
    )

    ax1 = axes[0]
    if xs1:
        ax1.bar(list(xs1), list(ys1), width=0.7,
                color="steelblue", edgecolor="white", linewidth=0.4)
    if args.log:
        ax1.set_yscale("log")
        ax1.set_ylabel("Count (log)", fontsize=11)
    else:
        ax1.set_ylabel("Count", fontsize=11)
        ax1.yaxis.set_major_formatter(
            FuncFormatter(lambda v, _: f"{v/1e6:.0f}M" if v >= 1e6 else f"{v:,.0f}")
        )
    ax1.set_title(f"Gaps 1 - {SPLIT1}  (dominant peaks)", fontsize=11, loc="left", pad=4)
    ax1.set_xlabel("Gap  (U_{n+1} - U_{n})", fontsize=10)
    ax1.set_xlim(0, SPLIT1 + 1)
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    ax1.grid(axis="y", alpha=0.3, linewidth=0.5)

    ax2 = axes[1]
    if xs2:
        ax2.bar(list(xs2), list(ys2), width=1.6,
                color="darkorange", edgecolor="white", linewidth=0.3)
    ax2.set_yscale("log")
    ax2.set_ylabel("Count (log)", fontsize=11)
    ax2.set_xlabel("Gap  (U_{n+1} - U_{n})", fontsize=10)
    ax2.set_title(f"Gaps {SPLIT1 + 1} - {SPLIT2}  (secondary structure, log y)", fontsize=11, loc="left", pad=4)
    ax2.set_xlim(SPLIT1, SPLIT2 + 5)
    ax2.xaxis.set_major_locator(MultipleLocator(25))
    ax2.xaxis.set_minor_locator(MultipleLocator(5))
    ax2.grid(axis="y", alpha=0.3, linewidth=0.5)

    ax3 = axes[2]
    ax3.scatter(xs_all, ys_all, s=8, color="forestgreen", alpha=0.75, linewidths=0)
    ax3.set_yscale("log")
    ax3.set_xscale("log")
    ax3.set_ylabel("Count (log)", fontsize=10)
    ax3.set_xlabel("Gap value (log)", fontsize=10)
    ax3.set_title("Full range -- log-log overview", fontsize=11, loc="left", pad=4)
    ax3.grid(alpha=0.25, linewidth=0.5)

    plt.tight_layout(rect=(0, 0, 1, 0.96))
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Histogram   ->  {args.output!r}")

def cmd_cosine(args: argparse.Namespace) -> None:
    count, max_computed = _meta(args.state)
    lam: float = args.lam
    max_val: int | None = args.max_val if args.max_val > 0 else None

    try:
        import mpmath  # type: ignore[import]
        mpmath.mp.dps = 50
        _lam_mp = mpmath.mpf(repr(lam))
        _confirm = lambda v: bool(mpmath.cos(_lam_mp * v) >= 0)
        precision_note = "mpmath 30-digit precision"
    except ImportError:
        lam_f = np.float64(lam)
        _confirm = lambda v: bool(np.cos(lam_f * np.float64(v)) >= 0.0)
        precision_note = "double precision"

    range_desc = f"<= {max_val:,}" if max_val else f"<= {max_computed:,}"
    print(f"Checking  cos({lam} * U_n) < 0  for Ulam numbers {range_desc}  [{precision_note}] ...")

    lam_f64 = np.float64(lam)
    n_chunks = (count + CHUNK - 1) // CHUNK
    exception_count = 0
    known        = {2, 3, 47, 69}
    known_found  : set[int] = set()
    unexpected_examples: list[int] = []

    out = args.output
    with open(out, "w") as f:
        f.write(f"# cos({lam} * U_n) >= 0  -  inequality does NOT hold\n")
        f.write(f"# Ulam numbers {range_desc}, precision: {precision_note}\n")
        f.write("# Expected per Wikipedia (U_n <= 10^9): {2, 3, 47, 69}\n\n")
        for chunk in _progress(_ulam_chunks(args.state), n_chunks, "Checking"):
            if max_val is not None:
                chunk = chunk[chunk <= max_val]
                if len(chunk) == 0:
                    break
            candidates = chunk[np.cos(lam_f64 * chunk.astype(np.float64)) >= 0.0]
            for v in candidates.tolist():
                vi = int(v)
                if _confirm(vi):
                    f.write(f"{vi}\n")
                    exception_count += 1
                    if vi in known:
                        known_found.add(vi)
                    elif len(unexpected_examples) < 10:
                        unexpected_examples.append(vi)

    print(f"Exceptions  ->  {out!r}  ({exception_count:,} found)")
    unexpected_count = exception_count - len(known_found)
    missing = known - known_found
    if unexpected_count == 0 and not missing:
        print("  OK  Matches known exceptions exactly: {2, 3, 47, 69}")
    if unexpected_count > 0:
        print(f"  !!  {unexpected_count:,} unexpected exceptions"
              f"  (first few): {unexpected_examples[:10]}")
    if missing:
        print(f"  !!  Known exceptions not found"
              f" (outside computed range?): {sorted(missing)}")

def cmd_line(args: argparse.Namespace) -> None:
    count, _ = _meta(args.state)
    n = min(args.n, count)
    print(f"Loading first {n:,} Ulam numbers ...")

    parts: list[np.ndarray] = []
    for chunk in _ulam_chunks(args.state, limit=n):
        parts.append(chunk)
    arr = np.concatenate(parts)[:n]
    max_val = int(arr[-1])

    presence_mb = (max_val + 1) / 1_000_000
    if presence_mb > 2_000:
        print(
            f"Warning: presence array would be {presence_mb:,.0f} MB "
            f"(max_val={max_val:,}). Skipping heatmap strip."
        )
        has_heatmap = False
    else:
        has_heatmap = True
    presence = np.zeros(max_val + 1, dtype=np.uint8) if has_heatmap else None
    if has_heatmap:
        assert presence is not None
        presence[arr] = 1

    nrows = 2 if has_heatmap else 1
    height_ratios = [4, 1] if has_heatmap else [1]
    fig, axes_all = plt.subplots(
        nrows, 1, figsize=(18, 4),
        gridspec_kw={"height_ratios": height_ratios},
        squeeze=False,
    )
    axes_row = axes_all[:, 0]
    fig.suptitle(
        f"First {n:,} Ulam numbers on the number line  "
        f"(U1 = {int(arr[0])},  U_{n:,} = {int(arr[-1]):,})",
        fontsize=13,
    )
    ax = axes_row[0]
    ax.axhline(0, color="black", lw=1.2)

    colors = matplotlib.colormaps["plasma"](np.linspace(0.1, 0.9, n))
    segs = [[(float(u), -0.45), (float(u), 0.45)] for u in arr]
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=0.8))

    label_n = min(25, n)
    for i, u in enumerate(arr[:label_n]):
        y = 0.68 if i % 2 == 0 else -0.68
        ax.text(float(u), y, str(int(u)),
                ha="center", va="center", fontsize=7)

    ax.set_xlim(-1, max_val + 2)
    ax.set_ylim(-1.1, 1.1)
    ax.set_yticks([])
    ax.set_xticks([])

    if has_heatmap:
        assert presence is not None
        ax2 = axes_row[1]
        ax2.imshow(
            presence.reshape(1, -1),
            cmap="plasma", aspect="auto",
            extent=[0, max_val + 1, 0, 1],
            interpolation="nearest",
        )
        ax2.set_yticks([])
        ax2.set_xlabel("Value", fontsize=11)
        ax2.set_xlim(-1, max_val + 2)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Number line ->  {args.output!r}")


def cmd_spiral(args: argparse.Namespace) -> None:
    size    = args.size
    max_val = size * size
    print(f"Building {size}x{size} Ulam spiral  (covers 1 - {max_val:,}) ...")

    ulams    = _load_ulams_up_to(args.state, max_val)
    ulam_set = set(ulams.tolist())
    print(f"  {len(ulam_set):,} Ulam numbers in range.")

    grid    = np.zeros((size, size), dtype=np.uint8)
    r, c    = size // 2, size // 2
    if 1 in ulam_set:
        grid[r, c] = 1

    # Spiral walk: right 1, up 1, left 2, down 2, right 3, up 3, ...
    n      = 2
    dr, dc = 0, 1
    step   = 1

    while n <= max_val:
        for _ in range(2):
            for _ in range(step):
                r += dr
                c += dc
                if 0 <= r < size and 0 <= c < size and n in ulam_set:
                    grid[r, c] = 1
                n += 1
                if n > max_val:
                    break
            dr, dc = -dc, dr
            if n > max_val:
                break
        step += 1

    _, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(grid, cmap="binary_r", interpolation="nearest")
    ax.set_title(
        f"Ulam spiral - {size}x{size}  ({len(ulam_set):,} Ulam numbers)",
        fontsize=14,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Spiral      ->  {args.output!r}")

def cmd_basis(args: argparse.Namespace) -> None:
    n: int = args.n
    _, max_computed = _meta(args.state)

    if max_computed < n:
        print(
            f"Warning: state covers up to {max_computed:,} but n = {n:,}. "
            "Results above max_computed may be incomplete."
        )

    print(f"Allocating DP array for 1..{n:,} ...")
    reachable = np.zeros(n + 1, dtype=bool)
    reachable[0] = True
    n_loaded = 0

    est_chunks = max(1, int(n / STEP) // CHUNK + 1)
    for chunk in _progress(_ulam_chunks(args.state), est_chunks, "DP"):
        in_range = chunk[chunk <= n]
        for u in in_range.tolist():
            reachable[u : n + 1] |= reachable[: n + 1 - u]
        n_loaded += len(in_range)
        if len(in_range) < len(chunk):
            break

    unreachable = (np.flatnonzero(~reachable[1:]) + 1).tolist()
    n_unreachable = len(unreachable)

    out = args.output
    with open(out, "w") as f:
        f.write(f"# Integers in [1, {n:,}] not expressible as a sum of distinct Ulam numbers\n")
        f.write(f"# Ulam basis size: {n_loaded:,}  (all Ulam numbers <= {n:,})\n")
        f.write(f"# Total unreachable: {n_unreachable:,}\n\n")
        f.writelines(f"{v}\n" for v in unreachable)

    print(f"  Ulam basis size : {n_loaded:,}")
    print(f"  Reachable       : {n - n_unreachable:,} / {n:,}  ({100*(n - n_unreachable)/n:.4f} %)")
    print(f"  Unreachable     : {n_unreachable:,}")
    if 0 < n_unreachable <= 20:
        print(f"  Full list       : {unreachable}")
    elif n_unreachable > 20:
        print(f"  First 10        : {unreachable[:10]}")
    print(f"Output          ->  {out!r}")

def cmd_residue(args: argparse.Namespace) -> None:
    n_bins: int = args.bins
    count, max_computed = _meta(args.state)
    print(f"Computing residue distribution for {count:,} Ulam numbers  ({n_bins:,} bins) ...")

    wave = np.float64(WAVE)
    bin_counts = np.zeros(n_bins, dtype=np.int64)
    n_chunks = (count + CHUNK - 1) // CHUNK

    for chunk in _progress(_ulam_chunks(args.state), n_chunks, "Reading"):
        r = np.mod(chunk.astype(np.float64), wave) / wave
        idx = np.clip((r * n_bins).astype(np.int64), 0, n_bins - 1)
        bin_counts += np.bincount(idx, minlength=n_bins)

    probs = bin_counts.astype(np.float64) * n_bins / max_computed
    bin_centers = (np.arange(n_bins) + 0.5) / n_bins
    bw = 1.0 / n_bins

    _, ax = plt.subplots(figsize=(12, 6))
    ax.bar(bin_centers, probs, width=bw * 0.98, color="#555555", edgecolor="none")

    for x in (1/3, 2/3):
        ax.axvline(x, color="#cc2222", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.axvspan(1/3, 2/3, alpha=0.05, color="green", zorder=0)

    density = count / max_computed
    ax.axhline(density, color="steelblue", linestyle=":", linewidth=1.2,
               label=f"Mean density ~={density:.5f}  (1 in {1/density:.2f})")
    ax.legend(fontsize=10, loc="upper left")

    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(
        f"Residue modulo lambda  (lambda = {WAVE},  normalized to [0, 1))",
        fontsize=11,
    )
    ax.set_ylabel("P(Ulam | residue)", fontsize=11)
    ax.set_title(
        f"Probability of being Ulam as a function of residue modulo lambda\n"
        f"{count:,} Ulam numbers  (max = {max_computed:,})   {n_bins:,} bins",
        fontsize=12,
    )
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Residue plot ->  {args.output!r}")


def cmd_growth(args: argparse.Namespace) -> None:
    count, max_computed = _meta(args.state)
    n_points: int = args.points
    step = max(1, count // n_points)
    actual_points = (count + step - 1) // step
    print(
        f"Sampling every {step:,}-th of {count:,} Ulam numbers "
        f"-> {actual_points:,} points ..."
    )

    idx_parts: list[np.ndarray] = []
    val_parts: list[np.ndarray] = []
    pos = 0
    n_chunks = (count + CHUNK - 1) // CHUNK

    for chunk in _progress(_ulam_chunks(args.state), n_chunks, "Reading"):
        n = len(chunk)
        first_local = int((-pos) % step)
        local_sel = np.arange(first_local, n, step, dtype=np.int64)
        if len(local_sel):
            idx_parts.append(pos + local_sel + 1)
            val_parts.append(chunk[local_sel])
        pos += n

    ns = np.concatenate(idx_parts).astype(np.float64)
    us = np.concatenate(val_parts).astype(np.float64)
    trend = ns * STEP

    _, ax = plt.subplots(figsize=(12, 6))
    ax.plot(ns, us, color="steelblue", linewidth=0.6, label=r"$U_n$")
    ax.plot(ns, trend, color="tomato", linewidth=1.0, linestyle="--", alpha=0.8,
            label=f"Linear trend  $U_n \\approx {STEP} \\cdot n$")

    if args.linear:
        ax.set_xscale("linear")
        ax.set_yscale("linear")
        scale_note = "linear"
    elif args.log:
        ax.set_xscale("linear")
        ax.set_yscale("log")
        scale_note = "log y"
    else:
        ax.set_xscale("log")
        ax.set_yscale("log")
        scale_note = "log-log"

    ax.set_xlabel("Index $n$", fontsize=11)
    ax.set_ylabel("$U_n$", fontsize=11)
    ax.set_title(
        f"Ulam number growth  ({actual_points:,} sample points  |  {count:,} total)\n"
        f"max = {max_computed:,}   scale: {scale_note}",
        fontsize=12,
    )
    ax.legend(fontsize=10)
    ax.grid(alpha=0.25, linewidth=0.5, which="both")
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Growth plot ->  {args.output!r}")

def main() -> None:
    p = argparse.ArgumentParser(
        prog="visualization_script.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-s", "--state", default="saved.bin", metavar="FILE",
        help="path to saved.bin  (default: saved.bin)",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pg = sub.add_parser("gaps", help="histogram of consecutive gaps + JSON")
    pg.add_argument("-o", "--output", default="ulam_gaps.png",  metavar="FILE")
    pg.add_argument("--json", default=None, metavar="FILE",
                    help="JSON output path  (default: <output>.json)")
    pg.add_argument("--from-json", default=None, metavar="FILE",
                    help="load pre-computed gap counts from FILE; skips state file")
    pg.add_argument("--log", action="store_true",
                    help="log-scale y axis on panel 1 (panels 2-3 always log)")
    pg.set_defaults(func=cmd_gaps)

    pc = sub.add_parser("cosine", help="find exceptions to the cosine inequality")
    pc.add_argument("-s",default="saved.bin", metavar="FILE")
    pc.add_argument("-o", "--output", default="ulam_cosine_exceptions.txt", metavar="FILE")
    pc.add_argument("--lam", type=float, default=LAM, metavar="FLOAT",
                    help=f"lam in cos(lam*Un) < 0  (default: {LAM})")
    pc.add_argument("--max-val", type=int, default=0, metavar="INT",
                    help="only check Ulam numbers <= this value  (default: all)")
    pc.set_defaults(func=cmd_cosine)

    pl = sub.add_parser("line", help="number-line visualization")
    pl.add_argument("-o", "--output", default="ulam_line.png",  metavar="FILE")
    pl.add_argument("-n", type=int, default=200, metavar="INT",
                    help="how many Ulam numbers to show  (default: 200)")
    pl.set_defaults(func=cmd_line)

    ps = sub.add_parser("spiral", help="Ulam-spiral grid image")
    ps.add_argument("-o", "--output", default="ulam_spiral.png", metavar="FILE")
    ps.add_argument("--size", type=int, default=500, metavar="INT",
                    help="grid side length in cells  (default: 500)")
    ps.set_defaults(func=cmd_spiral)

    pr = sub.add_parser(
        "residue",
        help="P(Ulam | residue mod lambda) distribution -- replicates Gibbs (2015) Fig. 1",
    )
    pr.add_argument("-o", "--output", default="ulam_residue.png", metavar="FILE")
    pr.add_argument("--bins", type=int, default=1200, metavar="INT",
                    help="number of histogram bins  (default: 1200, same as Gibbs paper)")
    pr.set_defaults(func=cmd_residue)

    pg2 = sub.add_parser("growth", help="plot U_n vs n growth curve")
    pg2.add_argument("-o", "--output", default="ulam_growth.png", metavar="FILE")
    pg2.add_argument("--points", type=int, default=10_000, metavar="INT",
                     help="number of sample points to plot  (default: 10000)")
    scale_grp = pg2.add_mutually_exclusive_group()
    scale_grp.add_argument("--log", action="store_true",
                           help="log y-axis only")
    scale_grp.add_argument("--linear", action="store_true",
                           help="fully linear axes  (default is log-log)")
    pg2.set_defaults(func=cmd_growth)

    pb = sub.add_parser(
        "basis",
        help="find integers in [1, N] not expressible as sums of distinct Ulam numbers",
    )
    pb.add_argument("n", type=int, metavar="N",
                    help="check all integers in [1, N]")
    pb.add_argument("-o", "--output", default="ulam_basis_unreachable.txt", metavar="FILE")
    pb.set_defaults(func=cmd_basis)

    args = p.parse_args()
    skip_state = args.cmd == "gaps" and getattr(args, "from_json", None)
    if not skip_state and not Path(args.state).exists():
        p.error(f"State file not found: {args.state!r}")
    args.func(args)


if __name__ == "__main__":
    main()
