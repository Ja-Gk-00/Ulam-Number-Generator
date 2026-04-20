# Ulam Number Generator

Fast Ulam-number generator. The heavy lifting is done by a C shared library
;Python provides a clean high-level API anda full CLI.

A number `k` is **Ulam** if it can be written as the sum of two distinct Ulam
numbers in exactly one way. The sequence starts `1, 2, 3, 4, 6, 8, 11, 13 …`
([OEIS A002858](https://oeis.org/A002858)).

---
## Requirements

| Tool | Version |
|------|---------|
| Python | >= 3.14 |
| gcc | any recent version |
| uv | >= 0.4 (package manager) |

---

## Installation

```bash
# clone and enter the project
git clone <repo-url>
cd Ulam-Number-Generator

# create venv + install dev dependencies with uv
uv sync --group dev

# compile the C shared library (one-time)
uv run python build.py
# or: python -m ulam build
```

---

## Python API

### `UlamComputer`

Stateful generator.  The same object can be reused for multiple queries and
extended incrementally.

```python
from ulam_number_generator import UlamComputer

comp = UlamComputer(track_pairs=False)   # track_pairs=True to record consecutive pairs
comp.compute_up_to(1_000_000)            # idempotent; extends previous computation

comp.is_ulam(47)          # True  | False  | None (not yet evaluated)
47 in comp                # same, raises ValueError if not yet evaluated
comp.count                # number of Ulam numbers found
comp.max_computed         # highest candidate evaluated
comp.get_ulam(0)          # 1  (0-based index)
comp.get_all_ulams()      # list[int]
list(comp)                # iterate in order

comp.get_pairs()          # list[tuple[int,int]]  - requires track_pairs=True

comp.save("state.bin")
comp2 = UlamComputer.from_file("state.bin")
```

### `calculate_all_ulams`

Function that handles loading, computing, saving, and signal-safe
auto-save on `SIGINT` / `SIGTERM`.

```python
from ulam_number_generator import calculate_all_ulams

ulams: list[int] = calculate_all_ulams(
    n,                          # compute all Ulam numbers <= n
    save_file="ck.bin",         # save state here on completion or interruption
    load_file="ck.bin",         # resume from this file if it exists
    track_pairs=True,           # record consecutive pairs
    pairs_file="pairs.csv",     # write pairs as CSV (requires track_pairs)
    verbose=True,               # print progress
)
```

Pass the **same path** to `save_file` and `load_file` for checkpoint/resume:

```python
# first run - computes and saves
calculate_all_ulams(10_000_000, save_file="ck.bin", load_file="ck.bin")

# interrupted and restarted - picks up exactly where it left off
calculate_all_ulams(10_000_000, save_file="ck.bin", load_file="ck.bin")
```

### `check_if_ulam`

```python
from ulam_number_generator import check_if_ulam

check_if_ulam(47)                          # computes up to 47, returns True
check_if_ulam(47, state_file="ck.bin")    # O(1) if 47 <= max_computed in file
```

---
## CLI

Run as `python -m ulam_number_generator <command> [options]`.

### `compute` - generate Ulam numbers up to N

```
python -m ulam compute N [options]
```

| Option | Description |
|--------|-------------|
| `--save FILE` | Save full state to FILE on completion **and** on Ctrl-C / SIGTERM |
| `--load FILE` | Resume from FILE if it exists; start fresh otherwise |
| `--output FILE` | Write resulting Ulam numbers (one per line) to FILE |
| `--pairs` | Track and store every consecutive pair `(U_k, U_{k+1})` |
| `--pairs-out FILE` | Write consecutive pairs as CSV to FILE (requires `--pairs`) |
| `--quiet` / `-q` | Suppress all progress output |

```bash
# compute, save checkpoint, track pairs
python -m ulam_number_generator compute 1000000 --save ck.bin --pairs --pairs-out pairs.csv

# resume a previous run
python -m ulam_number_generator compute 5000000 --load ck.bin --save ck.bin
```

### `check` - test a single number

```
python -m ulam check N [--state FILE]
```

Exits `0` if N is Ulam, `1` otherwise.

| Option | Description |
|--------|-------------|
| `--state FILE` | Use saved state for O(1) lookup when N ≤ max\_computed |

```bash
python -m ulam_number_generator check 47 --state ck.bin
# 47 IS a Ulam number.
```

### `list` - print Ulam numbers from a saved state

```
python -m ulam_number_generator list --state FILE [--limit K] [-q]
```

| Option | Description |
|--------|-------------|
| `--state FILE` | *(required)* State file to read from |
| `--limit K` | Print at most K numbers (0 = all; default 50) |

### `pairs` - print consecutive Ulam pairs from a saved state

```
python -m ulam_number_generator pairs --state FILE [--limit K] [-q]
```

Requires the state to have been computed with `--pairs` (or `track_pairs=True`).

| Option | Description |
|--------|-------------|
| `--state FILE` | *(required)* State file to read from |
| `--limit K` | Print at most K pairs (0 = all; default 50) |

### `info` - summarise a saved state file

```
python -m ulam_number_generator info --state FILE
```

Prints `max_computed`, Ulam count, whether pairs are tracked, and the first/last Ulam number.

### `build` - (re-)compile the C library

```
python -m ulam_number_generator build
# equivalent to: python build.py
```

### `adj-pairs` - print adjacent-integer pairs from a saved state

```
python -m ulam_number_generator adj-pairs --state FILE [--limit K] [-q]
```

Lists pairs `(n, n+1)` where both integers are Ulam numbers (e.g. `1 2`, `2 3`, `3 4`, `47 48`).

| Option | Description |
|--------|-------------|
| `--state FILE` | *(required)* State file to read from |
| `--limit K` | Print at most K pairs (0 = all; default 50) |

### `status` - quick progress peek (no full load)

```
python -m ulam_number_generator status FILE
```

Reads only the file header - instant even for multi-GB state files.

```
max_computed : 1,000,000,000
ulam_count   : 74,007,081
track_pairs  : False
```

---

## State files (`.bin`)

When you pass `--save` (or call `comp.save()`), the full computation state is
written to a binary file.  Loading it later with `--load` (or
`UlamComputer.from_file()`) resumes *exactly* where you left off - no
recomputation needed.

### What is stored

A state file contains three things:

| Contents | Purpose |
|----------|---------|
| **`max_computed`** | The upper bound `N` that has been fully evaluated. Every integer from 1 to `N` has been tested. |
| **Ulam sequence** | Every Ulam number found, stored as compressed int32 offsets from a linear trend (4 bytes each instead of 8 - halves memory for large runs). |
| **Membership bitset** | A packed bitset that answers "is `k` a Ulam number?" in O(1) for any `k <= max_computed`. Needed to resume computation correctly. |
| **Consecutive pairs** *(optional)* | Every pair `(U_k, U_{k+1})` - i.e. two Ulam numbers that are adjacent in the sequence. Only present if the run was started with `--pairs`. |

### What you can do with a loaded state

```python
comp = UlamComputer.from_file("state.bin")

# How far has this been computed?
comp.max_computed          # e.g. 1_000_000_000
comp.count                 # e.g. 74_007_081  (Ulam numbers found)

# O(1) membership test for any k <= max_computed
comp.is_ulam(47)           # True
comp.is_ulam(50)           # False
comp.is_ulam(2_000_000_000)  # None - not yet computed

# Look up Ulam numbers by index (0-based)
comp.get_ulam(0)           # 1
comp.get_ulam(1)           # 2
comp.get_ulam(99_999)      # 1_351_223  (the 100,000th Ulam number)

# Get all Ulam numbers as a list
comp.get_all_ulams()       # [1, 2, 3, 4, 6, 8, 11, …]

# Adjacent-integer pairs (n, n+1) - always available
comp.get_adjacent_pairs()  # [(1,2), (2,3), (3,4), (47,48), …]

# Sequence-consecutive pairs (U_k, U_{k+1}) - only if saved with --pairs
comp.get_pairs()           # [(1,2), (2,3), (3,4), (4,6), (6,8), …]

# Extend the computation and re-save
comp.compute_up_to(2_000_000_000)
comp.save("state.bin")
```

### Quick inspection without loading

```bash
# Instant header-only peek - works on multi-GB files
python -m ulam_number_generator status state.bin

# Full summary (loads the state, prints first/last Ulam, pair count)
python -m ulam_number_generator info --state state.bin
```

### File size

| `max_computed` | Ulam count | Approx. file size |
|---|---|---|
| 1,000,000 | ~74,000 | ~0.9 MB |
| 100,000,000 | ~7,400,000 | ~85 MB |
| 1,000,000,000 | ~74,000,000 | ~850 MB |
| 10,000,000,000 | ~740,000,000 | ~8.5 GB |

*(With `--pairs` enabled the file grows by an additional 16 bytes per Ulam number.)*

### Compatibility

Files from older runs are accepted automatically:

| Version | Format | Status |
|---------|--------|--------|
| V2 | `uint8 sc[]` (original sieve) | Loaded and converted |
| V3 | `uint64 ulams[]` + bitset | Loaded and converted |
| V4 | `int32 a_off[]` + bitset *(current)* | Native |

---

## Development

```bash
# lint
uv run ruff check .

# auto-fix lint issues
uv run ruff check --fix .

# format
uv run ruff format .

# tests
uv run pytest
```

---

## How it works

The algorithm is based on Philip Gibbs, *"An efficient method for computing
Ulam numbers"* ([arXiv:1508.0085](https://arxiv.org/abs/1508.0085)).