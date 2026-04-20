import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent
SRC  = ROOT / "c_src" / "ulam.c"
OUT  = ROOT / "ulam_number_generator"


def lib_name() -> str:
    if sys.platform == "win32":
        return "ulam.dll"
    if sys.platform == "darwin":
        return "libulam.dylib"
    return "libulam.so"


def build(verbose: bool = True) -> bool:
    out  = OUT / lib_name()
    base = ["gcc", "-O3", "-std=c11"]

    if sys.platform == "win32":
        cmd = [*base, "-shared", "-o", str(out), str(SRC)]
    elif sys.platform == "darwin":
        cmd = [*base, "-shared", "-fPIC", "-o", str(out), str(SRC)]
    else:
        cmd = [*base, "-shared", "-fPIC", "-o", str(out), str(SRC), "-lm"]

    if verbose:
        print(f"$ {' '.join(cmd)}")

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("Build FAILED:", file=sys.stderr)
        print(r.stderr, file=sys.stderr)
        return False

    if verbose:
        print(f"Built: {out}")
    return True


if __name__ == "__main__":
    sys.exit(0 if build() else 1)
