import ctypes
import subprocess
import sys
from pathlib import Path

_lib: ctypes.CDLL | None = None
_PKG_DIR = Path(__file__).parent
_ROOT_DIR = _PKG_DIR.parent
_SRC_FILE = _ROOT_DIR / "c_src" / "ulam.c"


def _lib_filename() -> str:
    if sys.platform == "win32":
        return "ulam.dll"
    if sys.platform == "darwin":
        return "libulam.dylib"
    return "libulam.so"


def _lib_path() -> Path:
    return _PKG_DIR / _lib_filename()


def build(verbose: bool = True) -> Path:
    out = _lib_path()
    flags = ["-O3", "-std=c11"]

    if sys.platform == "win32":
        cmd = ["gcc", *flags, "-shared", "-o", str(out), str(_SRC_FILE)]
    elif sys.platform == "darwin":
        cmd = ["gcc", *flags, "-shared", "-fPIC", "-o", str(out), str(_SRC_FILE)]
    else:
        cmd = ["gcc", *flags, "-shared", "-fPIC", "-o", str(out), str(_SRC_FILE), "-lm"]

    if verbose:
        print(f"Building ulam C library: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile ulam C library.\nCommand: {' '.join(cmd)}\nstderr:\n{result.stderr}"
        )
    if verbose:
        print(f"  → {out}")
    return out


def _setup_signatures(lib: ctypes.CDLL) -> None:
    u64 = ctypes.c_uint64
    u64p = ctypes.POINTER(ctypes.c_uint64)
    vp = ctypes.c_void_p
    c_int = ctypes.c_int
    c_char_p = ctypes.c_char_p

    lib.ulam_init.argtypes = [c_int]
    lib.ulam_init.restype = vp

    lib.ulam_free.argtypes = [vp]
    lib.ulam_free.restype = None

    lib.ulam_compute_up_to.argtypes = [vp, u64]
    lib.ulam_compute_up_to.restype = c_int

    lib.ulam_is_ulam.argtypes = [vp, u64]
    lib.ulam_is_ulam.restype = c_int

    lib.ulam_save.argtypes = [vp, c_char_p]
    lib.ulam_save.restype = c_int

    lib.ulam_load.argtypes = [c_char_p]
    lib.ulam_load.restype = vp

    lib.ulam_max_computed.argtypes = [vp]
    lib.ulam_max_computed.restype = u64

    lib.ulam_count.argtypes = [vp]
    lib.ulam_count.restype = u64

    lib.ulam_get.argtypes = [vp, u64]
    lib.ulam_get.restype = u64

    lib.ulam_get_ptr.argtypes = [vp]
    lib.ulam_get_ptr.restype = u64p

    lib.ulam_find_adjacent_pairs.argtypes = [vp, ctypes.POINTER(u64)]
    lib.ulam_find_adjacent_pairs.restype = u64p

    lib.ulam_free_ptr.argtypes = [vp]
    lib.ulam_free_ptr.restype = None

    lib.ulam_open_adj_stream.argtypes = [vp, c_char_p]
    lib.ulam_open_adj_stream.restype = c_int

    lib.ulam_close_adj_stream.argtypes = [vp]
    lib.ulam_close_adj_stream.restype = None

    lib.ulam_track_pairs.argtypes = [vp]
    lib.ulam_track_pairs.restype = c_int

    lib.ulam_pairs_count.argtypes = [vp]
    lib.ulam_pairs_count.restype = u64

    lib.ulam_get_pair.argtypes = [vp, u64, u64p, u64p]
    lib.ulam_get_pair.restype = None

    lib.ulam_pairs_ptr.argtypes = [vp]
    lib.ulam_pairs_ptr.restype = u64p


def get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib

    p = _lib_path()
    if not p.exists():
        build(verbose=True)

    _lib = ctypes.CDLL(str(p))
    _setup_signatures(_lib)
    return _lib
