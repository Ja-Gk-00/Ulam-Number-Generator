import os
import tempfile

import pytest

from ulam_number_generator import (
    StateInfo,
    UlamComputer,
    calculate_all_ulams,
    check_if_ulam,
    peek_state,
)

# OEIS A002858 - first 28 Ulam numbers
KNOWN_ULAMS = [
    1, 2, 3, 4, 6, 8, 11, 13, 16, 18,
    26, 28, 36, 38, 47, 48, 53, 57, 62, 69,
    72, 77, 82, 87, 97, 99, 102, 106,
]


class TestUlamComputer:

    def test_seeds(self):
        comp = UlamComputer()
        assert comp.max_computed == 2
        assert comp.count == 2
        assert comp.get_ulam(0) == 1
        assert comp.get_ulam(1) == 2

    def test_first_28(self):
        comp = UlamComputer()
        comp.compute_up_to(KNOWN_ULAMS[-1])
        result = comp.get_all_ulams()
        assert result == KNOWN_ULAMS

    def test_is_ulam_known(self):
        comp = UlamComputer()
        comp.compute_up_to(110)
        for u in KNOWN_ULAMS:
            assert comp.is_ulam(u) is True, f"Expected {u} to be Ulam"

    def test_is_ulam_false(self):
        comp = UlamComputer()
        comp.compute_up_to(30)
        not_ulam = {5, 7, 9, 10, 12, 14, 15, 17, 19, 20}
        for v in not_ulam:
            assert comp.is_ulam(v) is False, f"Expected {v} NOT to be Ulam"

    def test_is_ulam_not_computed(self):
        comp = UlamComputer()
        assert comp.is_ulam(5) is None

    def test_contains_operator(self):
        comp = UlamComputer()
        comp.compute_up_to(50)
        assert 47 in comp
        assert 10 not in comp

    def test_iter(self):
        comp = UlamComputer()
        comp.compute_up_to(20)
        assert list(comp) == [1, 2, 3, 4, 6, 8, 11, 13, 16, 18]

    def test_incremental_extension(self):
        comp = UlamComputer()
        comp.compute_up_to(20)
        first = comp.get_all_ulams()
        comp.compute_up_to(110)
        extended = comp.get_all_ulams()
        assert extended[:len(first)] == first
        assert extended[:len(KNOWN_ULAMS)] == KNOWN_ULAMS

    def test_idempotent(self):
        comp = UlamComputer()
        comp.compute_up_to(100)
        u1 = comp.get_all_ulams()
        comp.compute_up_to(50)
        u2 = comp.get_all_ulams()
        assert u1 == u2

    def test_repr(self):
        comp = UlamComputer()
        assert "UlamComputer" in repr(comp)

    def test_large(self):
        comp = UlamComputer()
        comp.compute_up_to(10_000)
        assert comp.count == 827

class TestOEISRegression:
    """
    Verify the Gibbs algorithm against the reference values published by
    Philip Gibbs (arXiv:1508.0085v2) and independently confirmed on OEIS.

    a(n) denotes the n-th Ulam number (1-indexed).
    Spot-check values:
        a(827)       =    9 985   (last Ulam <= 10 000)
        a(3 869)     =   49 989   (last Ulam <= 50 000)
        a(7 584)     =   99 933   (last Ulam <= 100 000)
        a(100 000)   =1 351 223   (Gibbs paper, Table 1)
        a(1 000 000) =13 509 072  (Gibbs paper, Table 1)
    """

    def test_count_at_10k(self):
        comp = UlamComputer()
        comp.compute_up_to(10_000)
        assert comp.count == 827
        assert comp.get_ulam(826) == 9_985

    def test_count_at_50k(self):
        comp = UlamComputer()
        comp.compute_up_to(50_000)
        assert comp.count == 3_869
        assert comp.get_ulam(3_868) == 49_989

    def test_count_at_100k(self):
        comp = UlamComputer()
        comp.compute_up_to(100_000)
        assert comp.count == 7_584
        assert comp.get_ulam(7_583) == 99_933

    def test_gibbs_100k_th_ulam(self):
        comp = UlamComputer()
        comp.compute_up_to(1_351_223)
        assert comp.count == 100_000
        assert comp.get_ulam(99_999) == 1_351_223

    def test_gibbs_1m_th_ulam(self):
        import pytest
        pytest.importorskip("time")
        comp = UlamComputer()
        comp.compute_up_to(13_509_072)
        assert comp.count == 1_000_000
        assert comp.get_ulam(999_999) == 13_509_072

class TestPairs:

    def test_pairs_basic(self):
        comp = UlamComputer(track_pairs=True)
        comp.compute_up_to(20)
        pairs = comp.get_pairs()
        ulams = comp.get_all_ulams()
        expected = [(ulams[i], ulams[i + 1]) for i in range(len(ulams) - 1)]
        assert pairs[0] == (1, 2)
        assert pairs == expected

    def test_pairs_not_tracked_raises(self):
        comp = UlamComputer(track_pairs=False)
        comp.compute_up_to(20)
        with pytest.raises(ValueError, match="track_pairs"):
            comp.get_pairs()

    def test_pairs_count_equals_ulams_minus_one(self):
        comp = UlamComputer(track_pairs=True)
        comp.compute_up_to(100)
        assert len(comp.get_pairs()) == comp.count - 1


class TestAdjacentPairs:

    def test_known_adjacent_pairs(self):
        comp = UlamComputer()
        comp.compute_up_to(110)
        adj = comp.get_adjacent_pairs()
        assert (1, 2) in adj
        assert (2, 3) in adj
        assert (3, 4) in adj
        assert (47, 48) in adj

    def test_no_non_adjacent(self):
        comp = UlamComputer()
        comp.compute_up_to(50)
        adj = comp.get_adjacent_pairs()
        for a, b in adj:
            assert b == a + 1, f"Expected b == a+1, got ({a}, {b})"

    def test_both_members_are_ulam(self):
        comp = UlamComputer()
        comp.compute_up_to(200)
        adj = comp.get_adjacent_pairs()
        ulam_set = set(comp.get_all_ulams())
        for a, b in adj:
            assert a in ulam_set
            assert b in ulam_set

    def test_single_ulam_has_no_pairs(self):
        comp = UlamComputer()
        comp.compute_up_to(2)
        adj = comp.get_adjacent_pairs()
        assert adj == [(1, 2)]

    def test_adjacent_pairs_file(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as f:
            path = f.name
        try:
            calculate_all_ulams(110, adjacent_pairs_file=path)
            with open(path) as fh:
                lines = fh.read().splitlines()
            assert lines[0] == "n,n+1"
            data = {tuple(int(x) for x in ln.split(",")) for ln in lines[1:]}
            assert (1, 2) in data
            assert (47, 48) in data
        finally:
            os.unlink(path)

class TestPersistence:

    def test_save_load_roundtrip(self):
        comp = UlamComputer()
        comp.compute_up_to(200)
        original = comp.get_all_ulams()

        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            loaded = UlamComputer.from_file(path)
            assert loaded.get_all_ulams() == original
            assert loaded.max_computed == comp.max_computed
        finally:
            os.unlink(path)

    def test_save_load_with_pairs(self):
        comp = UlamComputer(track_pairs=True)
        comp.compute_up_to(100)
        pairs_orig = comp.get_pairs()

        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            loaded = UlamComputer.from_file(path)
            assert loaded.track_pairs is True
            assert loaded.get_pairs() == pairs_orig
        finally:
            os.unlink(path)

    def test_resume_computation(self):
        comp = UlamComputer()
        comp.compute_up_to(100)
        partial = comp.get_all_ulams()

        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            resumed = UlamComputer.from_file(path)
            resumed.compute_up_to(200)
            full = resumed.get_all_ulams()
            assert full[:len(partial)] == partial
        finally:
            os.unlink(path)

    def test_load_bad_file_raises(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            f.write(b"\x00" * 16)
            path = f.name
        try:
            with pytest.raises(ValueError, match="Failed to load"):
                UlamComputer.from_file(path)
        finally:
            os.unlink(path)

    def test_load_missing_file_raises(self):
        with pytest.raises(ValueError, match="Failed to load"):
            UlamComputer.from_file("/nonexistent/path/ulam.bin")

class TestCalculateAllUlams:

    def test_basic(self):
        result = calculate_all_ulams(20)
        assert result == [1, 2, 3, 4, 6, 8, 11, 13, 16, 18]

    def test_save_and_load(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            r1 = calculate_all_ulams(150, save_file=path)
            r2 = calculate_all_ulams(150, load_file=path)
            assert r1 == r2
        finally:
            os.unlink(path)

    def test_checkpoint_resume(self):
        """Same file used for both save and load — checkpoint pattern."""
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            r1 = calculate_all_ulams(100, save_file=path, load_file=path)
            r2 = calculate_all_ulams(200, save_file=path, load_file=path)
            assert r2[:len(r1)] == r1
        finally:
            os.unlink(path)

    def test_pairs_file(self):
        with (
            tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as sf,
            tempfile.NamedTemporaryFile(delete=False, suffix=".csv",
                                        mode="w") as pf,
        ):
            spath, ppath = sf.name, pf.name
        try:
            calculate_all_ulams(
                50,
                save_file=spath,
                track_pairs=True,
                pairs_file=ppath,
            )
            with open(ppath) as fh:
                lines = fh.read().splitlines()
            assert lines[0] == "a,b"
            assert len(lines) > 2
        finally:
            os.unlink(spath)
            os.unlink(ppath)

class TestCheckIfUlam:

    @pytest.mark.parametrize("n", KNOWN_ULAMS)
    def test_known_ulams(self, n):
        assert check_if_ulam(n) is True

    @pytest.mark.parametrize("n", [5, 7, 9, 10, 12, 14, 15, 17])
    def test_known_non_ulams(self, n):
        assert check_if_ulam(n) is False

    def test_zero(self):
        assert check_if_ulam(0) is False

    def test_negative(self):
        assert check_if_ulam(-5) is False

    def test_with_state_file_in_scope(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            calculate_all_ulams(200, save_file=path)
            assert check_if_ulam(47, state_file=path) is True
            assert check_if_ulam(10, state_file=path) is False
        finally:
            os.unlink(path)

    def test_with_state_file_out_of_scope(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            calculate_all_ulams(20, save_file=path)
            assert check_if_ulam(47, state_file=path) is True
        finally:
            os.unlink(path)

    def test_with_missing_state_file(self):
        assert check_if_ulam(47, state_file="/nonexistent/path.bin") is True

class TestPeekState:

    def test_basic_fields(self):
        comp = UlamComputer()
        comp.compute_up_to(500)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            info = peek_state(path)
            assert isinstance(info, StateInfo)
            assert info.max_computed == 500
            assert info.ulam_count == comp.count
            assert info.track_pairs is False
            assert info.path == path
        finally:
            os.unlink(path)

    def test_track_pairs_true(self):
        comp = UlamComputer(track_pairs=True)
        comp.compute_up_to(100)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            info = peek_state(path)
            assert info.track_pairs is True
        finally:
            os.unlink(path)

    def test_repr(self):
        comp = UlamComputer()
        comp.compute_up_to(50)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            path = f.name
        try:
            comp.save(path)
            assert "max_computed" in repr(peek_state(path))
        finally:
            os.unlink(path)

    def test_missing_file_raises(self):
        with pytest.raises(ValueError, match="Cannot open"):
            peek_state("/nonexistent/file.bin")

    def test_bad_magic_raises(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            f.write(b"\x00" * 24)
            path = f.name
        try:
            with pytest.raises(ValueError, match="bad magic"):
                peek_state(path)
        finally:
            os.unlink(path)

    def test_too_small_raises(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as f:
            f.write(b"\x00" * 8)
            path = f.name
        try:
            with pytest.raises(ValueError, match="too small"):
                peek_state(path)
        finally:
            os.unlink(path)
