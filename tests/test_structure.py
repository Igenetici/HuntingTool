import pytest

from huntingtool.structure import classify_adjacent, fragment_equivalent, parse_read, revcomp
from tests.simulate import ADJ, amplicon


@pytest.mark.parametrize("cfg", list(ADJ))
@pytest.mark.parametrize("strand", ["+", "-"])
def test_parse_each_configuration(cfg, strand):
    s = amplicon(19, cfg, 9)
    if strand == "-":
        s = revcomp(s)
    st = parse_read("GATTACA" + s)
    assert st.status == "complete"
    assert st.strand == strand
    assert st.config == cfg
    assert st.cag == 19
    assert st.ccg == 9


def test_gain_positions():
    for k in range(1, 9):
        adj = "CAA" + "CAG" * (k - 1) + "CAACAGCCGCCA"
        assert classify_adjacent(adj) == ("CAACAG-dup" if k == 2 else f"CAA-gain G{k}")
    assert classify_adjacent("CAACAGCCGCCT") == "other"


def test_expected_amplicon_size():
    # canonical (CCG)7: 107 + 3n bp (143 bp for 12 CAG, 356 bp for 83 CAG)
    assert len(amplicon(12)) == 143
    assert len(amplicon(83)) == 356


def test_fragment_equivalent():
    # CAA-loss: pure CAG = fragment CAG + 2; CAACAG-dup: pure CAG = fragment CAG - 2
    assert fragment_equivalent(19, ADJ["canonical"]) == 19
    assert fragment_equivalent(23, ADJ["CAA-loss"]) == 21
    assert fragment_equivalent(15, ADJ["CAACAG-dup"]) == 17
    assert fragment_equivalent(16, ADJ["CAA-gain G1"]) == 17
    assert fragment_equivalent(16, ADJ["CCA-loss"]) is None


def test_no_repeat():
    assert parse_read("ACGT" * 50) is None


def test_truncated_read():
    s = amplicon(20)[:-30]
    assert parse_read(s).status in ("incomplete_ccg", "no_ccg")
