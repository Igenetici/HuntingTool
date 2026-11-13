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
    neighbour_ratio: float = 0.0           # max(count[cag-1], count[cag+1]) / count[cag]
    mode_cag: Optional[int] = None         # most frequent CAG count of the allele's reads (before stutter-aware calling)
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


def _caa_shift(adj: str) -> Optional[int]:
    """CAG gain of the canonical read produced by a CAA>CAG error in the first CAA codon of a
    CAA-gain adjacent sequence (CAACAG-duplication: 2; CAA-gain Gk: k); None otherwise."""
    from .structure import CANONICAL_ADJ
    if not adj.startswith("CAA") or adj == CANONICAL_ADJ:
        return None
    conv = "CAG" + adj[3:]
    if conv.endswith(CANONICAL_ADJ) and len(conv) > len(CANONICAL_ADJ):
        return (len(conv) - len(CANONICAL_ADJ)) // 3
    return None


def _window(dist: Counter, x: int, w: int) -> int:
    return sum(dist[y] for y in range(x - w, x + w + 1))


def _cfg_peaks(dist: Counter) -> list[int]:
    # raw modes resolve close non-expanded alleles; smoothing stabilises the broader
    # distributions of expanded alleles
    # a low relative height is accepted here: candidate modes are filtered afterwards on
    # their read support, which keeps a long allele that is under-amplified relative to a
    # short one (e.g. >80 CAG) from being missed
    return sorted(set([x for x in _peaks(dist, rel=0.02) if x <= 35] +
                      [x for x in _peaks(dist, rel=0.02, smooth=True) if x > 35]))


def _ccg_call(sub: Counter) -> tuple[Optional[int], float, Counter]:
    """Majority CCG length of a set of reads and the share of reads supporting it."""
    cd = Counter()
    for (c, ccg, cag), n in sub.items():
        if ccg is not None:
            cd[ccg] += n
    if not cd:
        return None, 0.0, cd
    best, nb = cd.most_common(1)[0]
    return best, nb / sum(cd.values()), cd


def call_alleles(hap_counts: Counter, adj_of: dict, min_fraction: float = 0.05,
                 min_reads: int = 10, neighbour_flag: float = 0.6,
                 split_ratio: float = 0.25, config_error_ratio: float = 0.2,
                 other_hap: Optional[Counter] = None, unclassified_fraction: float = 0.2,
                 stutter_ratio: float = 0.8, stutter_min_cag: int = 36
                 ) -> tuple[list[Allele], list[str]]:
    """hap_counts: Counter[(config, ccg, cag)] of complete reads.
    adj_of: config -> representative adjacent sequence (for fragment equivalents).
    other_hap: optional Counter[(adjacent_sequence, ccg, cag)] of complete reads whose
    repeat-adjacent sequence matched no known configuration. A sequence carried by at least
    `unclassified_fraction` of the resolved reads (far above the basecalling-error
    background) is treated as an allele with an unclassified configuration.
    Returns (alleles sorted by CAG, sample-level notes).

    1. CAG modes are searched separately within each repeat-adjacent configuration, so
       alleles with different configurations are never merged, and each mode is scored by
       the reads within +/-1 CAG (+/-2 above 35 CAG) of it. Broad tails (PCR stutter and
       truncated products) therefore do not inflate the support of a neighbouring allele.
    2. Up to two modes are retained; a minor mode of a non-canonical configuration lying
       within 3 CAG of a much larger mode is treated as a configuration miscall.
    3. The CCG length of each allele is the majority CCG length of its reads. A single CAG
       mode is split into two alleles only when its reads show two well-separated CCG
       lengths (two alleles of equal or adjacent CAG length distinguished by CCG); when two
       CAG modes are present no CCG split is attempted, because recombinant (chimeric)
       PCR products carry the CAG tract of one allele and the CCG tract of the other.
    4. Stutter-aware sizing of expanded alleles (mode >= `stutter_min_cag`): PCR stutter
       mainly removes repeats, and in expanded alleles the n-1 product is almost as abundant
       as the full-length product, so that the most frequent count falls on n-1 in about
       half of the samples. The allele is therefore assigned the longest CAG count, within
       2 repeats of the mode, whose read count is at least `stutter_ratio` of that of the
       mode. The mode is kept in `mode_cag`, and the allele is flagged when the two differ.
       `stutter_ratio=0` disables this step (v2.1.0 behaviour).
    """
    from .structure import fragment_equivalent

    notes: list[str] = []
    res = Counter({k: v for k, v in hap_counts.items() if k[0] != "other"})
    N = sum(res.values())
    if N == 0:
        return [], ["no read with a resolved repeat-adjacent configuration"]
    thr = max(min_reads, min_fraction * N)
    if other_hap:
        tot_adj = Counter()
        for (adj, ccg, cag), n in other_hap.items():
            tot_adj[adj] += n
        for adj, n in tot_adj.items():
            if n >= max(min_reads, unclassified_fraction * N):
                for (a2, ccg, cag), k in other_hap.items():
                    if a2 == adj:
                        res[(f"unclassified [{adj}]", ccg, cag)] += k
    by_cfg: dict[str, Counter] = defaultdict(Counter)
    for (cfg, ccg, cag), n in res.items():
        by_cfg[cfg][cag] += n

    cands = []  # (support, cfg, mode)
    for cfg, dist in by_cfg.items():
        if sum(dist.values()) < min_reads:
            continue
        for pk in _cfg_peaks(dist):
            cands.append((_window(dist, pk, 1 if pk <= 35 else 2), cfg, pk))
    cands.sort(key=lambda c: (-c[0], c[2]))

    accepted: list[tuple[int, str, int]] = []
    others: list[tuple[int, str, int]] = []
    for s, cfg, pk in cands:
        if s < thr:
            continue
        miscall = any(cfg != "canonical" and acfg != cfg and abs(apk - pk) <= 3 and s < config_error_ratio * asup
                      for asup, acfg, apk in accepted)
        # a single CAA>CAG basecalling error in the first CAA codon of a CAA-gain allele
        # (e.g. CAACAG-duplication) yields a canonical read k repeats longer
        miscall = miscall or any(cfg == "canonical" and _caa_shift(adj_of.get(acfg, "")) == pk - apk
                                 and s < config_error_ratio * asup for asup, acfg, apk in accepted)
        if miscall or len(accepted) == 2:
            others.append((s, cfg, pk))
        else:
            accepted.append((s, cfg, pk))
    if others:
        notes.append("additional minor modes not called: " +
                     ", ".join(f"{pk}{'' if cfg == 'canonical' else ' ' + cfg} ({s})" for s, cfg, pk in others))

    def reads_of(cfg, pk):
        w = 1 if pk <= 35 else 2
        return Counter({k: v for k, v in res.items() if k[0] == cfg and abs(k[2] - pk) <= w})

    called: list[Allele] = []
    for s, cfg, pk in accepted:
        sub = reads_of(cfg, pk)
        ccg, share, cd = _ccg_call(sub)
        split = None
        if len(accepted) == 1:
            cpk = _peaks(cd, rel=split_ratio, valley=0.6)
            cpk = [c for c in cpk if c is not None]
            if len(cpk) >= 2:
                cpk = sorted(cpk, key=lambda c: -cd[c])[:2]
                if abs(cpk[0] - cpk[1]) >= 2 and cd[cpk[1]] >= thr:
                    split = cpk
        groups = []
        if split:
            # widen the CAG window so that two alleles one CAG apart are both captured
            wide = Counter({k: v for k, v in res.items() if k[0] == cfg and abs(k[2] - pk) <= 2})
            for cp in split:
                groups.append((cp, Counter({k: v for k, v in wide.items()
                                            if k[1] is not None and _nearest(k[1], split) == cp})))
        else:
            groups.append((ccg, sub))
        for cp, g in groups:
            cdist = Counter()
            for (c, gg, cag), n in g.items():
                cdist[cag] += n
            mode = max(cdist, key=lambda x: (cdist[x], -abs(x - pk)))
            full = by_cfg[cfg]
            raw_mode = mode
            if stutter_ratio and mode >= stutter_min_cag:
                top = max(full[k] for k in range(mode - 2, mode + 3))
                mode = max(k for k in range(mode - 2, mode + 3) if full[k] >= stutter_ratio * top)
            a = Allele(mode, cfg, cp, sum(g.values()))
            a.mode_cag = raw_mode
            if mode != raw_mode:
                a.flags.append("stutter-aware call: longest CAG count with at least "
                               f"{int(round(stutter_ratio * 100))}% of the reads of the most frequent count "
                               f"({raw_mode} CAG: {full[raw_mode]} reads; {mode} CAG: {full[mode]} reads); "
                               "repeat count uncertain by 1 repeat")
            a.neighbour_ratio = max(full[mode - 1], full[mode + 1]) / full[mode] if full[mode] else 0.0
            _, sh, _ = _ccg_call(g)
            a.config_purity = 100.0 * sum(n for (c, gg, cag), n in res.items() if abs(cag - mode) <= (1 if mode <= 35 else 2)
                                          and c == cfg) / max(1, sum(n for (c, gg, cag), n in res.items()
                                                                     if abs(cag - mode) <= (1 if mode <= 35 else 2)))
            if not split and sh < 0.6:
                a.flags.append(f"CCG length supported by {sh * 100:.0f}% of the allele's reads")
            called.append(a)

    tot = sum(a.support for a in called) or 1
    for a in called:
        a.fraction = a.support / tot
        if a.config.startswith("unclassified"):
            a.fe_cag = None
            a.flags.append("repeat-adjacent sequence not matching any known configuration "
                           "(sequence shown); confirm by an orthogonal method")
        else:
            a.fe_cag = fragment_equivalent(a.cag, adj_of.get(a.config, ""))
        if a.cag <= 35 and a.neighbour_ratio >= neighbour_flag:
            a.flags.append(f"neighbour bin >= {int(neighbour_flag * 100)}% of mode: possible second allele "
                           "differing by one CAG repeat")
    if len(called) == 1:
        a = called[0]
        notes.append("single allele detected (homozygous or unresolved)")
        if a.cag <= 35 and a.neighbour_ratio >= neighbour_flag:
            full = by_cfg[a.config]
            other = a.cag - 1 if full[a.cag - 1] >= full[a.cag + 1] else a.cag + 1
            notes.append(f"possible genotype {min(a.cag, other)}/{max(a.cag, other)} "
                         "(two alleles differing by one CAG repeat, not resolved by CCG or configuration)")
        if a.cag <= 35:
            # CAG mode within each well-represented CCG length: a one-repeat shift between
            # CCG classes supports two alleles of adjacent CAG length
            near = Counter({k: v for k, v in res.items() if k[0] == a.config and abs(k[2] - a.cag) <= 2
                            and k[1] is not None})
            tot_near = sum(near.values())
            per: dict[int, Counter] = defaultdict(Counter)
            for (c, ccg, cag), n in near.items():
                per[ccg][cag] += n
            modes = {g: max(cd, key=cd.get) for g, cd in per.items()
                     if sum(cd.values()) >= max(min_reads, 0.15 * tot_near)}
            if len(set(modes.values())) > 1:
                notes.append("CAG mode differs between CCG classes (" +
                             ", ".join(f"CCG{g}: {m}" for g, m in sorted(modes.items())) +
                             "): supports two alleles of adjacent CAG length")
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
