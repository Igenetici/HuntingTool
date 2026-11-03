"""Synthetic HTT exon-1 amplicon reads for testing (no patient data)."""
from __future__ import annotations

import gzip
import random

from huntingtool.structure import revcomp

UP = "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC"   # 51 nt, forward primer to flank
DOWN = "CCTCCTCAGCTTCCTCAGCCGCC"                              # CCT repeat to end of reverse primer
ADJ = {
    "canonical": "CAACAGCCGCCA",
    "CAA-loss": "CCGCCA",
    "CCA-loss": "CAACAG",
    "CAACAGCCGCCA-loss": "",
    "CAACAG-dup": "CAACAGCAACAGCCGCCA",
    "CAA-gain G1": "CAACAACAGCCGCCA",
}


def amplicon(cag: int, config: str = "canonical", ccg: int = 7) -> str:
    return UP + "CAG" * cag + ADJ[config] + "CCG" * ccg + DOWN


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = []
    for b in seq:
        r = rng.random()
        if r < rate / 3:                      # substitution
            out.append(rng.choice("ACGT".replace(b, "")))
        elif r < 2 * rate / 3:                # deletion
            continue
        elif r < rate:                        # insertion
            out.append(b + rng.choice("ACGT"))
        else:
            out.append(b)
    return "".join(out)


def reads(alleles: list[tuple[int, str, int, float]], n: int, err: float = 0.01, seed: int = 1):
    """alleles: (cag, config, ccg, weight). Yields (seq, qual)."""
    rng = random.Random(seed)
    tot = sum(a[3] for a in alleles)
    for _ in range(n):
        x = rng.random() * tot
        for cag, cfg, ccg, w in alleles:
            x -= w
            if x <= 0:
                break
        s = mutate(amplicon(cag, cfg, ccg), err, rng)
        s = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 15))) + s
        if rng.random() < 0.5:
            s = revcomp(s)
        q = chr(33 + 15) * len(s)
        yield s, q


def write_fastq(path: str, alleles, n: int, err: float = 0.01, seed: int = 1):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        for i, (s, q) in enumerate(reads(alleles, n, err, seed)):
            fh.write(f"@read{i}\n{s}\n+\n{q}\n")
