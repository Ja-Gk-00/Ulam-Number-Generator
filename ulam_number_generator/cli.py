import argparse
import sys

from .core import UlamComputer, calculate_all_ulams, check_if_ulam, peek_state


def cmd_compute(args: argparse.Namespace) -> None:
    result = calculate_all_ulams(
        args.n,
        save_file=args.save,
        load_file=args.load,
        track_pairs=args.pairs,
        pairs_file=args.pairs_out,
        adjacent_pairs_file=args.adj_pairs_out,
        verbose=not args.quiet,
    )

    pct = len(result) / args.n * 100 if args.n else 0.0
    summary = f"{len(result):,} / {args.n:,} Ulam numbers ({pct:.4f} %)"

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("\n".join(map(str, result)) + "\n")
        if not args.quiet:
            print(f"Written {summary} to {args.output!r}.")
    elif not args.quiet:
        if len(result) <= 50:
            print(" ".join(map(str, result)))
            print(summary)
        else:
            print(summary)
            print(f"First 10: {result[:10]}")
            print(f"Last  10: {result[-10:]}")


def cmd_check(args: argparse.Namespace) -> None:
    result = check_if_ulam(args.n, state_file=args.state)
    label = "IS a Ulam number" if result else "is NOT a Ulam number"
    print(f"{args.n} {label}.")
    sys.exit(0 if result else 1)


def cmd_list(args: argparse.Namespace) -> None:
    comp = UlamComputer.from_file(args.state)
    cnt = comp.count
    limit = min(args.limit, cnt) if args.limit else cnt
    for i in range(limit):
        print(comp.get_ulam(i))
    if not args.quiet:
        remaining = cnt - limit
        if remaining:
            print(f"… ({remaining:,} more, use --limit 0 for all)")


def cmd_pairs(args: argparse.Namespace) -> None:
    comp = UlamComputer.from_file(args.state)
    if not comp.track_pairs:
        print(
            "Error: the state file was not computed with --pairs enabled.",
            file=sys.stderr,
        )
        sys.exit(1)
    pairs = comp.get_pairs()
    limit = min(args.limit, len(pairs)) if args.limit else len(pairs)
    for a, b in pairs[:limit]:
        print(f"{a} {b}")
    if not args.quiet:
        remaining = len(pairs) - limit
        if remaining:
            print(f"… ({remaining:,} more pairs, use --limit 0 for all)")


def cmd_info(args: argparse.Namespace) -> None:
    comp = UlamComputer.from_file(args.state)
    print(f"State file      : {args.state}")
    print(f"max_computed    : {comp.max_computed:,}")
    print(f"Ulam count      : {comp.count:,}")
    print(f"track_pairs     : {comp.track_pairs}")
    if comp.track_pairs:
        print(f"Pairs stored    : {int(comp._lib.ulam_pairs_count(comp._state)):,}")
    if comp.count >= 2:
        print(f"First Ulam      : {comp.get_ulam(0)}")
        print(f"Last  Ulam      : {comp.get_ulam(comp.count - 1)}")


def cmd_adj_pairs(args: argparse.Namespace) -> None:
    comp = UlamComputer.from_file(args.state)
    adj = comp.get_adjacent_pairs()
    limit = min(args.limit, len(adj)) if args.limit else len(adj)
    for a, b in adj[:limit]:
        print(f"{a} {b}")
    if not args.quiet:
        remaining = len(adj) - limit
        if remaining:
            print(f"… ({remaining:,} more pairs, use --limit 0 for all)")


def cmd_status(args: argparse.Namespace) -> None:
    try:
        info = peek_state(args.file)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"max_computed : {info.max_computed:,}")
    print(f"ulam_count   : {info.ulam_count:,}")
    print(f"track_pairs  : {info.track_pairs}")


def cmd_build(_args: argparse.Namespace) -> None:
    from ._loader import build

    try:
        build(verbose=True)
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ulam",
        description="Ulam number generator - fast C-backed Python module",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  ulam compute 10000 --save ck.bin --pairs --pairs-out pairs.csv\n"
            "  ulam compute 50000 --load ck.bin --save ck.bin   # resume\n"
            "  ulam check   47   --state ck.bin\n"
            "  ulam list         --state ck.bin --limit 20\n"
            "  ulam pairs        --state ck.bin --limit 10\n"
            "  ulam info         --state ck.bin\n"
            "  ulam build\n"
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("compute", help="Compute Ulam numbers up to N")
    p.add_argument("n", type=int, metavar="N", help="compute all Ulam numbers ≤ N")
    p.add_argument(
        "--save", metavar="FILE", help="save state to FILE on completion or interruption"
    )
    p.add_argument("--load", metavar="FILE", help="resume computation from FILE (if it exists)")
    p.add_argument(
        "--output", metavar="FILE", help="write resulting Ulam numbers (one per line) to FILE"
    )
    p.add_argument(
        "--pairs", action="store_true", help="track and store every consecutive pair (U_k, U_{k+1})"
    )
    p.add_argument(
        "--pairs-out",
        metavar="FILE",
        dest="pairs_out",
        help="write sequence-consecutive pairs (U_k,U_{k+1}) to CSV (requires --pairs)",
    )
    p.add_argument(
        "--adj-pairs-out",
        metavar="FILE",
        dest="adj_pairs_out",
        help="write adjacent-integer pairs (n, n+1 both Ulam) to CSV FILE",
    )
    p.add_argument("--quiet", "-q", action="store_true", help="suppress progress output")
    p.set_defaults(func=cmd_compute, adj_pairs_out=None)

    p = sub.add_parser("check", help="Check if N is a Ulam number")
    p.add_argument("n", type=int, metavar="N")
    p.add_argument(
        "--state", metavar="FILE", help="use saved state for O(1) lookup when N is in scope"
    )
    p.set_defaults(func=cmd_check, quiet=False)

    p = sub.add_parser("list", help="List Ulam numbers from a saved state")
    p.add_argument("--state", metavar="FILE", required=True)
    p.add_argument(
        "--limit",
        type=int,
        metavar="K",
        default=50,
        help="max numbers to print (0 = all, default 50)",
    )
    p.add_argument("--quiet", "-q", action="store_true")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("pairs", help="List consecutive Ulam pairs from a saved state")
    p.add_argument("--state", metavar="FILE", required=True)
    p.add_argument(
        "--limit",
        type=int,
        metavar="K",
        default=50,
        help="max pairs to print (0 = all, default 50)",
    )
    p.add_argument("--quiet", "-q", action="store_true")
    p.set_defaults(func=cmd_pairs)

    p = sub.add_parser(
        "adj-pairs",
        help="List adjacent-integer Ulam pairs (n, n+1) from a saved state",
    )
    p.add_argument("--state", metavar="FILE", required=True)
    p.add_argument(
        "--limit",
        type=int,
        metavar="K",
        default=50,
        help="max pairs to print (0 = all, default 50)",
    )
    p.add_argument("--quiet", "-q", action="store_true")
    p.set_defaults(func=cmd_adj_pairs)

    p = sub.add_parser("status", help="Quickly show progress of a state file")
    p.add_argument("file", metavar="FILE", help="state file to inspect")
    p.set_defaults(func=cmd_status, quiet=False)

    p = sub.add_parser("info", help="Show summary of a saved state file")
    p.add_argument("--state", metavar="FILE", required=True)
    p.set_defaults(func=cmd_info, quiet=False)

    p = sub.add_parser("build", help="(Re-)compile the C shared library")
    p.set_defaults(func=cmd_build, quiet=False)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
