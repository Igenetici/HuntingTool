"""Caller behaviour on read-count patterns observed in real amplicon data."""
from collections import Counter

from huntingtool.caller import call_alleles

ADJ = {"canonical": "CAACAGCCGCCA"}


def _stutter(cag, ccg, n, cfg="canonical", tail=0.5):
    """Mode plus a left PCR-stutter tail (n-1 ~ tail*n, n-2 ~ tail^2*n ...)."""
    h = Counter()
    for d in range(0, 6):
        h[(cfg, ccg, cag - d)] += int(n * tail ** d)
    h[(cfg, ccg, cag + 1)] += int(n * 0.1)
    return h


def test_chimeric_ccg_does_not_create_third_allele():
    # short allele CCG10, long allele CCG7; recombinant reads carry 18 CAG with CCG7
    h = _stutter(18, 10, 1000) + _stutter(18, 7, 900) + _stutter(40, 7, 3000, tail=0.7)
    al, notes = call_alleles(h, ADJ)
    assert [(a.cag, a.ccg) for a in al] == [(18, 10), (40, 7)]


def test_broad_tail_does_not_mask_short_allele():
    # strongly under-represented short allele next to a long allele with a broad tail
    h = _stutter(40, 7, 6000, tail=0.8)
    for c in range(20, 35):
        h[("canonical", 7, c)] += 60
    h += _stutter(17, 7, 900)   # ~6% of reads, as in the least represented real short allele
    al, _ = call_alleles(h, ADJ)
    assert [a.cag for a in al] == [17, 40]


def test_under_amplified_long_allele_called():
    h = _stutter(22, 7, 3000)
    for c in range(30, 76):
        h[("canonical", 7, c)] += 40
    h += _stutter(82, 7, 300, tail=0.85)
    al, _ = call_alleles(h, ADJ)
    assert [a.cag for a in al] == [22, 82]


def test_unclassified_configuration_reported():
    h = _stutter(12, 7, 1000)
    other = Counter({("CAGCAGCAACAGCAACAG", 4, 15): 900, ("CAGCAGCAACAGCAACAG", 4, 14): 350,
                     ("CAGCAGCAACAGCCGCCA", 7, 13): 15})
    al, _ = call_alleles(h, ADJ, other_hap=other)
    assert [(a.cag, a.config) for a in al] == [(12, "canonical"), (15, "unclassified [CAGCAGCAACAGCAACAG]")]
    assert al[1].fe_cag is None


def test_sporadic_unresolved_sequences_ignored():
    h = _stutter(17, 7, 1000) + _stutter(21, 7, 1000)
    other = Counter({("CAGCAGCAACAGCCACCA", 7, 21): 60, ("AGCAGCAACAGCCGCCAC", 7, 17): 50})
    al, _ = call_alleles(h, ADJ, other_hap=other)
    assert [(a.cag, a.config) for a in al] == [(17, "canonical"), (21, "canonical")]


def test_adjacent_alleles_shift_between_ccg_classes_noted():
    h = _stutter(16, 10, 800) + _stutter(17, 9, 1000)
    al, notes = call_alleles(h, ADJ)
    assert any("CAG mode differs between CCG classes" in n or len(al) == 2 for n in notes)


def test_dup_conversion_product_not_called():
    # CAACAG-dup 15 + canonical 20; A>G errors in the first CAA of the dup allele give canonical 17 reads
    h = _stutter(15, 8, 1800, cfg="CAACAG-dup", tail=0.3) + _stutter(17, 8, 210, tail=0.6) + _stutter(20, 7, 180)
    al, notes = call_alleles(h, {"canonical": "CAACAGCCGCCA", "CAACAG-dup": "CAACAGCAACAGCCGCCA"})
    assert [(a.cag, a.config) for a in al] == [(15, "CAACAG-dup"), (20, "canonical")]


def _expanded(n_true, n, ccg=7, minus1=0.97, plus1=0.3):
    """Expanded-allele profile observed in real data: n-1 almost as high as n, sharp drop at n+1."""
    h = Counter()
    for d, f in ((-4, 0.2), (-3, 0.38), (-2, 0.66), (-1, minus1), (0, 1.0), (1, plus1), (2, 0.07)):
        h[("canonical", ccg, n_true + d)] += int(n * f)
    return h


def test_stutter_aware_call_when_n_minus_1_is_the_mode():
    # as in real 17/40 samples: 39 CAG slightly more frequent than 40 CAG
    h = _stutter(17, 7, 3000) + _expanded(40, 1200, minus1=1.04)
    al, _ = call_alleles(h, ADJ)
    assert [a.cag for a in al] == [17, 40]
    assert al[1].mode_cag == 39 and any("stutter-aware" in f for f in al[1].flags)
    al, _ = call_alleles(h, ADJ, stutter_ratio=0)
    assert [a.cag for a in al] == [17, 39]


def test_stutter_aware_call_keeps_mode_when_it_is_the_longest_abundant_count():
    h = _stutter(17, 7, 3000) + _expanded(40, 1200, minus1=0.95)
    al, _ = call_alleles(h, ADJ)
    assert [a.cag for a in al] == [17, 40]
    assert al[1].mode_cag == 40 and not any("stutter-aware" in f for f in al[1].flags)


def test_stutter_aware_call_not_applied_below_36():
    h = _stutter(17, 7, 3000) + Counter({("canonical", 7, 26): 1000, ("canonical", 7, 27): 900,
                                         ("canonical", 7, 25): 500, ("canonical", 7, 28): 100})
    al, _ = call_alleles(h, ADJ)
    assert [a.cag for a in al] == [17, 26]
