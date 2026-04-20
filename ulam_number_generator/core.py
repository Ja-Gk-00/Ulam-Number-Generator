import ctypes
import os
import signal
import struct
from collections.abc import Iterator
from types import FrameType
from typing import Self

from tqdm import tqdm

from ._loader import get_lib


class UlamComputer:
    def __init__(self, track_pairs: bool = False) -> None:
        self._lib = get_lib()
        self._state = self._lib.ulam_init(int(track_pairs))
        if not self._state:
            raise MemoryError("ulam_init returned NULL — out of memory")

    def __del__(self) -> None:
        if getattr(self, "_state", None) and getattr(self, "_lib", None):
            self._lib.ulam_free(self._state)
            self._state = None

    @classmethod
    def from_file(cls, path: str | os.PathLike[str]) -> Self:
        lib = get_lib()
        state = lib.ulam_load(str(path).encode())
        if not state:
            raise ValueError(
                f"Failed to load Ulam state from {path!r}. "
                "File may be missing, corrupt, or from a different version."
            )
        obj = cls.__new__(cls)
        obj._lib = lib
        obj._state = state
        return obj

    def compute_up_to(self, n: int) -> None:
        if not self._lib.ulam_compute_up_to(self._state, int(n)):
            raise MemoryError(f"Out of memory while computing Ulam numbers up to {n:,}")

    def is_ulam(self, n: int) -> bool | None:
        r = self._lib.ulam_is_ulam(self._state, int(n))
        return None if r == -1 else bool(r)

    @property
    def max_computed(self) -> int:
        return int(self._lib.ulam_max_computed(self._state))

    @property
    def count(self) -> int:
        return int(self._lib.ulam_count(self._state))

    @property
    def track_pairs(self) -> bool:
        return bool(self._lib.ulam_track_pairs(self._state))

    def get_ulam(self, index: int) -> int:
        return int(self._lib.ulam_get(self._state, int(index)))

    def get_all_ulams(self) -> list[int]:
        cnt = self.count
        if cnt == 0:
            return []
        return [int(self._lib.ulam_get(self._state, i)) for i in range(cnt)]

    def get_adjacent_pairs(self) -> list[tuple[int, int]]:
        out_count = ctypes.c_uint64(0)
        raw = self._lib.ulam_find_adjacent_pairs(self._state, ctypes.byref(out_count))
        n = int(out_count.value)
        if n == 0 or not raw:
            return []
        addr = ctypes.cast(raw, ctypes.c_void_p).value
        assert addr is not None
        arr = (ctypes.c_uint64 * (n * 2)).from_address(addr)
        result = [(int(arr[i * 2]), int(arr[i * 2 + 1])) for i in range(n)]
        self._lib.ulam_free_ptr(raw)
        return result

    def get_pairs(self) -> list[tuple[int, int]]:
        if not self.track_pairs:
            raise ValueError(
                "Consecutive pairs are not tracked. Create the UlamComputer with track_pairs=True."
            )
        cnt = int(self._lib.ulam_pairs_count(self._state))
        if cnt == 0:
            return []
        raw = self._lib.ulam_pairs_ptr(self._state)
        addr = ctypes.cast(raw, ctypes.c_void_p).value
        assert addr is not None, "ulam_pairs_ptr returned NULL"
        arr = (ctypes.c_uint64 * (cnt * 2)).from_address(addr)
        return [(int(arr[i * 2]), int(arr[i * 2 + 1])) for i in range(cnt)]

    def save(self, path: str | os.PathLike[str]) -> None:
        """Persist full state (including sum-count table) to *path*."""
        if not self._lib.ulam_save(self._state, str(path).encode()):
            raise OSError(f"Failed to save Ulam state to {path!r}")

    def __repr__(self) -> str:
        return (
            f"UlamComputer(max_computed={self.max_computed:,}, "
            f"count={self.count:,}, track_pairs={self.track_pairs})"
        )

    def __contains__(self, n: int) -> bool:
        r = self.is_ulam(n)
        if r is None:
            raise ValueError(
                f"{n} has not been evaluated yet (max_computed={self.max_computed:,}). "
                "Call compute_up_to first."
            )
        return r

    def __iter__(self) -> Iterator[int]:
        for i in range(self.count):
            yield int(self._lib.ulam_get(self._state, i))


def _compute_with_progress(comp: UlamComputer, n: int, *, verbose: bool) -> None:

    start = comp.max_computed
    chunk = max(1, (n - start) // 200)
    cur = start
    with tqdm(
        total=n,
        initial=cur,
        unit="",
        unit_scale=True,
        desc="Computing",
        disable=not verbose,
        dynamic_ncols=True,
        bar_format=(
            "{desc}: {percentage:5.1f}% |{bar}| {n_fmt}/{total_fmt}  [{elapsed}<{remaining}]"
        ),
    ) as bar:
        while cur < n:
            nxt = min(cur + chunk, n)
            comp.compute_up_to(nxt)
            bar.update(nxt - cur)
            bar.set_postfix_str(f"Ulam: {comp.count:,}", refresh=False)
            cur = nxt


def calculate_all_ulams(
    n: int,
    *,
    save_file: str | None = None,
    load_file: str | None = None,
    track_pairs: bool = False,
    pairs_file: str | None = None,
    adjacent_pairs_file: str | None = None,
    verbose: bool = False,
) -> list[int]:

    comp: UlamComputer | None = None
    if load_file and os.path.exists(load_file):
        try:
            if verbose:
                print(f"Loading state from {load_file!r} …")
            comp = UlamComputer.from_file(load_file)
            if verbose:
                print(
                    f"  Resumed: max_computed={comp.max_computed:,}, "
                    f"Ulams found so far={comp.count:,}"
                )
        except (ValueError, OSError) as exc:
            if verbose:
                print(f"  Could not load state ({exc}); starting fresh.")
            comp = None
    if comp is None:
        comp = UlamComputer(track_pairs=track_pairs)

    def _save_and_exit(signum: int, frame: FrameType | None) -> None:
        if save_file:
            print(
                f"\nInterrupted (signal {signum}) — saving state to {save_file!r} …",
                flush=True,
            )
            try:
                comp.save(save_file)  # type: ignore[union-attr]
                print("  Saved.", flush=True)
            except Exception as exc:
                print(f"  Save failed: {exc}", flush=True)
        raise SystemExit(1)

    old_sigint = signal.signal(signal.SIGINT, _save_and_exit)
    old_sigterm = signal.signal(signal.SIGTERM, _save_and_exit)

    try:
        already = comp.max_computed
        if adjacent_pairs_file:
            ok = comp._lib.ulam_open_adj_stream(comp._state, adjacent_pairs_file.encode())
            if not ok:
                raise OSError(f"Could not open adjacent-pairs stream: {adjacent_pairs_file!r}")

        try:
            if already >= n:
                if verbose:
                    pct = comp.count / n * 100
                    print(
                        f"State already covers up to {already:,} — no computation needed.\n"
                        f"  {comp.count:,} / {n:,} Ulam numbers ({pct:.4f} %)"
                    )
            else:
                _compute_with_progress(comp, n, verbose=verbose)
        finally:
            if adjacent_pairs_file:
                comp._lib.ulam_close_adj_stream(comp._state)

        pct = comp.count / n * 100
        if verbose:
            print(f"Done.  {comp.count:,} / {n:,} Ulam numbers ({pct:.4f} %)")

        if save_file:
            comp.save(save_file)
            if verbose:
                print(f"State saved to {save_file!r}.")

        if track_pairs and pairs_file:
            pairs = comp.get_pairs()
            with open(pairs_file, "w", newline="") as fh:
                fh.write("a,b\n")
                for a, b in pairs:
                    fh.write(f"{a},{b}\n")
            if verbose:
                print(f"Sequence-consecutive pairs saved to {pairs_file!r} ({len(pairs):,} pairs).")

        if adjacent_pairs_file and verbose:
            print(f"Adjacent-integer pairs streamed to {adjacent_pairs_file!r}.")

        return comp.get_all_ulams()

    finally:
        signal.signal(signal.SIGINT, old_sigint)
        signal.signal(signal.SIGTERM, old_sigterm)


def check_if_ulam(
    n: int,
    *,
    state_file: str | None = None,
) -> bool:
    if n <= 0:
        return False

    if state_file and os.path.exists(state_file):
        comp = UlamComputer.from_file(state_file)
        result = comp.is_ulam(n)
        if result is not None:
            return result
        comp.compute_up_to(n)
        return bool(comp.is_ulam(n))

    comp = UlamComputer()
    comp.compute_up_to(n)
    return bool(comp.is_ulam(n))


_MAGIC = 0x554C414D
_VERSION = 4  # compressed int32 a_off format
_VERSION3 = 3  # Gibbs uint64 ulams + bitset (still loadable)
_VERSION2 = 2  # legacy sc[] format (still loadable)
_HEADER_SIZE = 24  # magic(4) + version(4) + max_computed(8) + count(8)


class StateInfo:
    __slots__ = ("max_computed", "path", "track_pairs", "ulam_count")

    def __init__(
        self,
        path: str,
        max_computed: int,
        ulam_count: int,
        track_pairs: bool,
    ) -> None:
        self.path = path
        self.max_computed = max_computed
        self.ulam_count = ulam_count
        self.track_pairs = track_pairs

    def __repr__(self) -> str:
        return (
            f"StateInfo(path={self.path!r}, "
            f"max_computed={self.max_computed:,}, "
            f"ulam_count={self.ulam_count:,}, "
            f"track_pairs={self.track_pairs})"
        )


def peek_state(path: str | os.PathLike[str]) -> StateInfo:
    spath = str(path)
    try:
        with open(spath, "rb") as fh:
            raw = fh.read(_HEADER_SIZE + 4)
    except OSError as exc:
        raise ValueError(f"Cannot open state file {spath!r}: {exc}") from exc

    if len(raw) < _HEADER_SIZE:
        raise ValueError(f"State file too small (got {len(raw)} bytes): {spath!r}")

    magic, version = struct.unpack_from("<II", raw, 0)
    if magic != _MAGIC:
        raise ValueError(f"Not a ulam state file (bad magic): {spath!r}")
    if version not in (_VERSION, _VERSION3, _VERSION2):
        raise ValueError(
            f"Unsupported state file version {version} "
            f"(expected {_VERSION2}, {_VERSION3}, or {_VERSION}): {spath!r}"
        )

    max_computed, ulam_count = struct.unpack_from("<QQ", raw, 8)

    track_pairs = _peek_track_pairs(spath, ulam_count, version)

    return StateInfo(
        path=spath,
        max_computed=int(max_computed),
        ulam_count=int(ulam_count),
        track_pairs=track_pairs,
    )


def _peek_track_pairs(path: str, ulam_count: int, version: int) -> bool:
    """
    Layout after the 24-byte header:
      V2: ulams[](8*count)  sc_cap(8)  sc[](1*sc_cap)  track_pairs(4)
      V3: ulams[](8*count)  k_cap(8)   k[](4*k_cap)    track_pairs(4)
      V4: a_off[](4*count)  k_cap(8)   k[](4*k_cap)    track_pairs(4)
    """
    try:
        with open(path, "rb") as fh:
            fh.seek(24)
            seq_bytes = ulam_count * 4 if version == _VERSION else ulam_count * 8
            fh.seek(seq_bytes, 1)
            cap_raw = fh.read(8)
            if len(cap_raw) < 8:
                return False
            (cap,) = struct.unpack_from("<Q", cap_raw)
            skip = cap if version == _VERSION2 else cap * 4
            fh.seek(skip, 1)
            tp_raw = fh.read(4)
            if len(tp_raw) < 4:
                return False
            (tp,) = struct.unpack_from("<i", tp_raw)
            return bool(tp)
    except OSError:
        return False
