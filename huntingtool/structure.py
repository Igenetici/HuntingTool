"""Read-level parsing of the HTT exon 1 repeat region.

Each read is oriented, anchored on the 5' flank immediately upstream of the CAG
tract and parsed into three components:

    5'-flank | (CAG)n | repeat-adjacent sequence | (CCG)m | CCT...

The repeat-adjacent sequence (the bases between the end of the uninterrupted CAG
tract and the start of the CCG repeat) defines the configuration:

    canonical            (CAG)n CAACAGCCGCCA (CCG)m
    CAA-loss             (CAG)n CCGCCA (CCG)m
    CCA-loss             (CAG)n CAACAG (CCG)m
    CAACAGCCGCCA-loss    (CAG)n (CCG)m          (also termed CAA/CCA-loss)
    CAACAG-duplication   (CAG)n CAACAGCAACAGCCGCCA (CCG)m   (= CAA-gain G2)
    CAA-gain Gk          (CAG)n CAA(CAG)k-1 CAACAGCCGCCA (CCG)m, k = 1..8

Nomenclature follows the GeM-HD Consortium (Nat Genet 2025).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

# 12 nt immediately upstream of the CAG tract in HTT exon 1
# (...GAG TCC CTC AAG TCC TTC | CAG CAG ...)
FLANK5 = "CTCAAGTCCTTC"
FLANK5_SHORT = "AGTCCTTC"
CANONICAL_ADJ = "CAACAGCCGCCA"

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_GAIN = re.compile(r"^CAA((?:CAG){0,7})CAACAGCCGCCA$")

FIXED = {
    CANONICAL_ADJ: "canonical",
    "CCGCCA": "CAA-loss",
    "CAACAG": "CCA-loss",
    "": "CAACAGCCGCCA-loss",
}

# display order and colours used across reports
CONFIG_ORDER = ["canonical", "CAA-loss", "CCA-loss", "CAACAGCCGCCA-loss", "CAACAG-dup"] + [
    f"CAA-gain G{k}" for k in (1, 3, 4, 5, 6, 7, 8)
] + ["other"]


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def classify_adjacent(adj: str) -> str:
    """Name the repeat-adjacent configuration from the adjacent sequence."""
    if adj in FIXED:
        return FIXED[adj]
    m = _GAIN.match(adj)
    if m:
        k = len(m.group(1)) // 3 + 1
        return "CAACAG-dup" if k == 2 else f"CAA-gain G{k}"
    return "other"


def fragment_equivalent(cag: int, adj: str) -> Optional[int]:
    """CAG count expected from fragment analysis with a primer ending in the
    CAA-CAG-CCG-CCA region, i.e. the number of codons between the 5' flank and the
    CCGCCA motif minus the four codons of the canonical CAACAGCCGCCA sequence.
    Defined only when the adjacent sequence ends with CCGCCA."""
    if not adj.endswith("CCGCCA"):
        return None
    return cag + len(adj) // 3 - 4


def _find_flank(seq: str, limit: int = 400) -> int:
    """Return the index just after the 5' flank, or -1."""
    i = seq.find(FLANK5, 0, limit)
    if i >= 0:
        return i + len(FLANK5)
    i = seq.find(FLANK5_SHORT, 0, limit)
    if i >= 0:
        return i + len(FLANK5_SHORT)
    # one mismatch tolerated on the 12-mer immediately preceding the first CAG triplet run
    p = seq.find("CAGCAGCAG", 0, limit + 30)
    L = len(FLANK5)
    if p >= L:
        mism = sum(1 for a, b in zip(seq[p - L:p], FLANK5) if a != b)
        if mism <= 1:
            return p
    return -1


# known repeat-adjacent sequences, longest first (suffix matching on the CAG tract region)
_KNOWN_ADJ = sorted(
    ["CAA" + "CAG" * (k - 1) + CANONICAL_ADJ for k in range(1, 9)] + [CANONICAL_ADJ, "CAACAG", "CCGCCA"],
    key=len, reverse=True)
_CCT_ANCHORS = (("CCTCCTCAG", 0), ("CCTCAGCTTCC", -3), ("CAGCTTCCTCAG", -6))


@dataclass
class ReadStructure:
    status: str              # 'complete', 'incomplete_ccg', 'no_flank', 'no_ccg'
    strand: str              # '+' or '-'
    cag: int = 0             # uninterrupted CAG length (codons between flank and adjacent sequence)
    adj: str = ""            # repeat-adjacent sequence
    config: str = ""         # configuration name
    ccg: Optional[int] = None
    pure: bool = False       # CAG tract basecalled without any error


def orient(seq: str) -> tuple[Optional[str], str]:
    """Return the read in HTT sense orientation, or (None, '') if no repeat is seen."""
    f = seq.count("CAGCAGCAG")
    r = seq.count("CTGCTGCTG")
    if f == 0 and r == 0:
        return None, ""
    if f >= r:
        return seq, "+"
    return revcomp(seq), "-"


def split_tract(region: str) -> tuple[str, str]:
    """Split the region between the 5' flank and the CCG repeat into
    (CAG tract, repeat-adjacent sequence) by exact suffix matching."""
    for adj in _KNOWN_ADJ:
        if region.endswith(adj) and region[: -len(adj)].endswith("CAGCAG"):
            return region[: -len(adj)], adj
    if region.endswith("CAGCAGCAG"):
        return region, ""                      # CAACAGCCGCCA-loss: (CAG)n directly followed by (CCG)n
    return region[:-12], None                  # unresolved junction (sequencing error or novel variant)


def parse_read(seq: str) -> Optional[ReadStructure]:
    """Parse one read; None if it carries no CAG/CTG repeat at all.

    The CAG-tract length is measured between two anchors (the 5' flank and the start of
    the CCG repeat) and divided by three, which makes it tolerant to isolated
    basecalling errors inside the tract; the configuration is read from the exact
    sequence at the 3' end of the tract."""
    s, strand = orient(seq)
    if s is None:
        return None
    start = _find_flank(s)
    if start < 0:
        return ReadStructure("no_flank", strand)
    c = s.find("CCGCCG", start)
    if c < 0:
        return ReadStructure("no_ccg", strand)
    tract, adj = split_tract(s[start:c])
    cag = round(len(tract) / 3)
    pure = tract == "CAG" * cag
    if adj is None:
        cfg, adj = "other", s[max(start, c - 18):c]
    else:
        cfg = classify_adjacent(adj)
    ccg = None
    for anchor, off in _CCT_ANCHORS:
        e = s.find(anchor, c, c + 120)
        if e >= 0:
            ccg = round((e + off - c) / 3)
            break
    return ReadStructure("complete" if ccg is not None else "incomplete_ccg", strand, cag, adj, cfg, ccg, pure)


def longest_run(seq: str, motif: str = "CAG") -> int:
    """Longest uninterrupted run of `motif` on either strand (generic motif mode)."""
    best = 0
    for s in (seq, revcomp(seq)):
        for m in re.finditer(f"(?:{re.escape(motif)})+", s):
            n = (m.end() - m.start()) // len(motif)
            if n > best:
                best = n
    return best
