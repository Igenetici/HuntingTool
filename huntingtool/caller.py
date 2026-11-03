"""Allele calling and phasing from read-level structures.

1. Allele discovery: modes of the uninterrupted-CAG distribution of all reads whose
   repeat-adjacent configuration was resolved.
2. Phasing: within the reads assigned to each CAG mode, the configuration and the CCG
   length are taken by majority; if two configurations or two CCG lengths are both well
   supported, the mode is split into two alleles (e.g. alleles with similar CAG length
   but different CCG length).
Because CAG, configuration and CCG are measured on the same read, each called allele
carries its own configuration and CCG length.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Allele:
    cag: int
    config: str
    ccg: Optional[int]
    support: int
    fraction: float = 0.0
    fe_cag: Optional[int] = None           # fragment-analysis-equivalent CAG count
    config_purity: float = 0.0             # % of resolved reads of this allele carrying its configuration
    neighbour_ratio: float = 0.0           # max(count[mode-1], count[mode+1]) / count[mode]
    flags: list[str] = field(default_factory=list)

    @property
    def label(self) -> str:
        c = "" if self.config == "canonical" else f" {self.config}"
        g = f" CCG{self.ccg}" if self.ccg else ""
        return f"{self.cag}{c}{g}"


def _smooth(dist: Counter) -> Counter:
    return Counter({x: dist[x - 1] + 2 * dist[x] + dist[x + 1] for x in range(min(dist) - 1, max(dist) + 2)})


def _peaks(dist: Counter, rel: float = 0.15, valley: float = 0.75, smooth: bool = False) -> list[int]:
    """Local maxima of a count distribution, separated from each other by a real valley."""
    dist = Counter({k: v for k, v in dist.items() if v > 0})
    if not dist:
        return []
    d = _smooth(dist) if smooth else dist
    xs = range(min(d) - 1, max(d) + 2)
    top = max(d.values())
    cand = [x for x in xs if d[x] > 0 and d[x] >= d[x - 1] and d[x] > d[x + 1] and d[x] >= rel * top]
    cand.sort(key=lambda x: -d[x])
    kept: list[int] = []
    for x in cand:
        if all(min(d[y] for y in range(min(x, k), max(x, k) + 1)) <= valley * min(d[x], d[k]) for k in kept):
            kept.append(x)
    if smooth:  # move each peak to the raw mode in its neighbourhood
        kept = [max((y for y in (x - 1, x, x + 1)), key=lambda y: dist[y]) for x in kept]
    return sorted(set(kept))


def _partition(dist: Counter, peaks: list[int]) -> dict[int, int]:
    sup = {p: 0 for p in peaks}
    for x, n in dist.items():
        sup[min(peaks, key=lambda q: (abs(q - x), -dist[q]))] += n
    return sup


def _nearest(x, peaks):
    return min(peaks, key=lambda q: abs(q - x))


def call_alleles(hap_counts: Counter, adj_of: dict, min_fraction: float = 0.10,
                 min_reads: int = 10, neighbour_flag: float = 0.6,
                 split_ratio: float = 0.25) -> tuple[list[Allele], list[str]]:
    """hap_counts: Counter[(config, ccg, cag)] of complete reads.
    adj_of: config -> representative adjacent sequence (for fragment equivalents).
    Returns (alleles sorted by CAG, sample-level notes)."""
    from .structure import fragment_equivalent

    notes: list[str] = []
    res = Counter({k: v for k, v in hap_counts.items() if k[0] != "other"})
    N = sum(res.values())
    if N == 0:
        return [], ["no read with a resolved repeat-adjacent configuration"]
    thr = max(min_reads, min_fraction * N)
    cag_dist: Counter = Counter()
    for (cfg, ccg, cag), n in res.items():
        cag_dist[cag] += n
    # raw modes resolve close non-expanded alleles; smoothing stabilises the broader
    # distributions of expanded alleles
    peaks = sorted(set([x for x in _peaks(cag_dist) if x <= 35] +
                       [x for x in _peaks(cag_dist, smooth=True) if x > 35]))

    # assign reads to CAG modes
    groups: dict[int, Counter] = defaultdict(Counter)
    for (cfg, ccg, cag), n in res.items():
        groups[_nearest(cag, peaks)][(cfg, ccg, cag)] += n

    cands: list[tuple[Allele, Counter]] = []
    for pk, g in groups.items():
        # split by configuration when two configurations are both well supported
        by_cfg = Counter()
        for (cfg, ccg, cag), n in g.items():
            by_cfg[cfg] += n
        top = by_cfg.most_common()
        cfgs = [top[0][0]] + [c for c, n in top[1:2] if n >= split_ratio * top[0][1] and n >= thr]
        for cfg in cfgs:
            sub = Counter({k: v for k, v in g.items() if k[0] == cfg}) if len(cfgs) > 1 else g
            ccg_dist = Counter()
            for (c, ccg, cag), n in sub.items():
                if c == cfg and ccg is not None:
                    ccg_dist[ccg] += n
            cpk = _peaks(ccg_dist, rel=split_ratio, valley=0.6) or [None]
            for cp in cpk:
                ss = Counter({k: v for k, v in sub.items()
                              if cp is None or k[1] is None or _nearest(k[1], cpk) == cp})
                if len(cpk) > 1 and sum(ss.values()) < thr:
                    continue
                cd = Counter()
                for (c, ccg, cag), n in ss.items():
                    cd[cag] += n
                mode = max(cd, key=lambda x: (cd[x], -abs(x - pk)))
                a = Allele(mode, cfg, cp, sum(ss.values()))
                a.neighbour_ratio = max(cd[mode - 1], cd[mode + 1]) / cd[mode]
                same = sum(n for (c, ccg, cag), n in ss.items() if c == cfg)
                a.config_purity = 100.0 * same / a.support
                cands.append((a, ss))

    called = [a for a, _ in cands if a.support >= thr]
    called.sort(key=lambda a: -a.support)
    if len(called) > 2:
        notes.append("more than two allele candidates above threshold: " +
                     ", ".join(a.label + f" ({a.support})" for a in called))
        called = called[:2]
    tot = sum(a.support for a in called) or 1
    for a in called:
        a.fraction = a.support / tot
        a.fe_cag = fragment_equivalent(a.cag, adj_of.get(a.config, ""))
        if a.neighbour_ratio >= neighbour_flag:
            a.flags.append(f"neighbour bin >= {int(neighbour_flag * 100)}% of mode: possible second allele "
                           "differing by one CAG repeat")
        if a.config_purity < 70:
            a.flags.append(f"configuration supported by {a.config_purity:.0f}% of the allele's resolved reads")
    if len(called) == 1:
        a = called[0]
        notes.append("single allele detected (homozygous or unresolved)")
        if a.neighbour_ratio >= neighbour_flag:
            other = a.cag - 1 if cag_dist[a.cag - 1] >= cag_dist[a.cag + 1] else a.cag + 1
            notes.append(f"possible genotype {min(a.cag, other)}/{max(a.cag, other)} "
                         "(two alleles differing by one CAG repeat, not resolved by CCG or configuration)")
    if len(called) == 2:
        mn = min(a.support for a in called) / max(a.support for a in called)
        if mn < 0.33:
            notes.append(f"marked allelic imbalance (minor/major read ratio {mn:.2f})")
    called.sort(key=lambda a: (a.cag, a.support))
    return called, notes


def canonical_percentages(cfg_counts: Counter) -> dict[str, float]:
    """v1-compatible metric: canonical/(canonical + X) * 100 for each non-canonical X."""
    can = cfg_counts.get("canonical", 0)
    out = {}
    for cfg, n in cfg_counts.items():
        if cfg in ("canonical", "other"):
            continue
        out[cfg] = 100.0 * can / (can + n) if (can + n) else float("nan")
    return out
