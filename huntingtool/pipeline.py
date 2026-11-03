"""Per-sample analysis and multi-sample orchestration."""
from __future__ import annotations

import csv
import gzip
import json
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Optional

from . import __version__
from .caller import Allele, call_alleles, canonical_percentages, _peaks, _partition
from .fastq import Sample, iter_fastq, mean_qscore
from .structure import longest_run, parse_read


@dataclass
class Options:
    mode: str = "htt"            # 'htt' (structure-aware) or 'motif' (generic longest-run counting)
    motif: str = "CAG"
    min_length: int = 100
    max_length: int = 0          # 0 = no limit
    min_q: float = 0.0           # minimum mean read Q (0 = no filter)
    min_repeats: int = 5         # motif mode only
    min_fraction: float = 0.10   # minimum fraction of informative reads to call an allele
    min_reads: int = 10
    write_reads: bool = False
    html: bool = True


def analyze_sample(sample: Sample, outdir: str, opt: Options) -> dict:
    t0 = time.time()
    out = Path(outdir) / sample.name
    out.mkdir(parents=True, exist_ok=True)
    n = n_pass = n_rep = 0
    status = Counter()
    strand = Counter()
    hap = Counter()            # (config, ccg, cag) of complete reads
    adj_of: dict[str, str] = {}
    other_adj = Counter()
    len_all = Counter()        # read length bins (10 bp) of all reads
    len_rep = Counter()        # read length bins of repeat-containing reads
    q_all = Counter()          # mean Q (integer) of all reads
    motif_counts = Counter()
    cag_by_cfg_incomplete = Counter()

    rw = None
    if opt.write_reads:
        rw = gzip.open(out / "reads.tsv.gz", "wt")
        rw.write("read_index\tlength\tmean_q\tstrand\tstatus\tcag\tconfiguration\tadjacent\tccg\n")

    for f in sample.files:
        for seq, qual in iter_fastq(f):
            n += 1
            L = len(seq)
            len_all[min(L, 2000) // 10] += 1
            q = mean_qscore(qual)
            q_all[int(q)] += 1
            if L < opt.min_length or (opt.max_length and L > opt.max_length) or q < opt.min_q:
                continue
            n_pass += 1
            if opt.mode == "motif":
                r = longest_run(seq, opt.motif)
                if r >= opt.min_repeats:
                    n_rep += 1
                    len_rep[min(L, 2000) // 10] += 1
                    motif_counts[r] += 1
                continue
            st = parse_read(seq)
            if st is None:
                continue
            n_rep += 1
            len_rep[min(L, 2000) // 10] += 1
            status[st.status] += 1
            strand[st.strand] += 1
            if st.status == "complete":
                hap[(st.config, st.ccg, st.cag)] += 1
                adj_of.setdefault(st.config, st.adj)
                if st.config == "other":
                    other_adj[st.adj] += 1
            elif st.status == "incomplete_ccg":
                cag_by_cfg_incomplete[(st.config, st.cag)] += 1
            if rw:
                rw.write(f"{n}\t{L}\t{q:.1f}\t{st.strand}\t{st.status}\t{st.cag}\t{st.config}\t{st.adj}\t"
                         f"{'' if st.ccg is None else st.ccg}\n")
    if rw:
        rw.close()

    res: dict = dict(sample=sample.name, files=[str(f) for f in sample.files], version=__version__,
                     options=asdict(opt), n_reads=n, n_pass_filters=n_pass, n_repeat_reads=n_rep)
    if opt.mode == "motif":
        pk = _peaks(motif_counts)
        sup = _partition(motif_counts, pk) if pk else {}
        tot = sum(motif_counts.values()) or 1
        al = [Allele(p, f"motif {opt.motif}", None, sup[p]) for p in pk if sup[p] >= max(opt.min_reads, opt.min_fraction * tot)]
        al = sorted(sorted(al, key=lambda a: -a.support)[:2], key=lambda a: a.cag)
        s = sum(a.support for a in al) or 1
        for a in al:
            a.fraction = a.support / s
        res.update(alleles=[asdict(a) for a in al], notes=[], motif_counts=dict(motif_counts),
                   genotype="/".join(str(a.cag) for a in al) or "NA")
    else:
        alleles, notes = call_alleles(hap, adj_of, opt.min_fraction, opt.min_reads)
        cfg_counts = Counter()
        for (cfg, ccg, cag), k in hap.items():
            cfg_counts[cfg] += k
        n_complete = sum(hap.values())
        if n_rep and n_complete / n_rep < 0.5:
            notes.append(f"only {100*n_complete/n_rep:.0f}% of repeat reads fully parsed "
                         "(consider higher-accuracy basecalling)")
        res.update(
            alleles=[asdict(a) for a in alleles], notes=notes,
            genotype="/".join(str(a.cag) for a in alleles) or "NA",
            n_complete=n_complete, status=dict(status), strand=dict(strand),
            config_counts=dict(cfg_counts), canonical_pct=canonical_percentages(cfg_counts),
            other_adjacent=other_adj.most_common(10),
            hap_counts=[[c, g, a, k] for (c, g, a), k in sorted(hap.items(), key=lambda x: (str(x[0][0]), x[0][1] or 0, x[0][2]))],
        )
        with open(out / "haplotype_counts.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["configuration", "ccg", "cag", "reads"])
            for c, g, a, k in res["hap_counts"]:
                w.writerow([c, g, a, k])
    res.update(read_length_hist=sorted(len_all.items()), repeat_read_length_hist=sorted(len_rep.items()),
               mean_q_hist=sorted(q_all.items()), runtime_s=round(time.time() - t0, 2))
    with open(out / "result.json", "w") as fh:
        json.dump(res, fh, indent=1)
    if opt.html:
        from .report import sample_figures
        res["figures"] = sample_figures(res, out)
    return res


def run(samples: list[Sample], outdir: str, opt: Options, threads: int = 1,
        progress: Optional[Callable[[str], None]] = None) -> list[dict]:
    log = progress or (lambda m: None)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    if threads <= 1 or len(samples) == 1:
        for i, s in enumerate(samples, 1):
            log(f"[{i}/{len(samples)}] {s.name}: {len(s.files)} file(s)")
            r = analyze_sample(s, outdir, opt)
            log(f"    genotype {r['genotype']}  ({r['n_repeat_reads']} repeat reads, {r['runtime_s']} s)")
            results.append(r)
    else:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            fut = {ex.submit(analyze_sample, s, outdir, opt): s for s in samples}
            for i, f in enumerate(as_completed(fut), 1):
                r = f.result()
                log(f"[{i}/{len(samples)}] {r['sample']}: genotype {r['genotype']} "
                    f"({r['n_repeat_reads']} repeat reads, {r['runtime_s']} s)")
                results.append(r)
    results.sort(key=lambda r: r["sample"])
    write_summary(results, outdir, opt)
    if opt.html:
        from .report import write_html
        p = write_html(results, outdir, opt)
        log(f"Report: {p}")
    return results


def write_summary(results: list[dict], outdir: str, opt: Options) -> Path:
    p = Path(outdir) / "summary.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "genotype", "allele1_cag", "allele1_configuration", "allele1_ccg", "allele1_reads",
                    "allele1_fraction", "allele1_fragment_equivalent_cag", "allele2_cag", "allele2_configuration",
                    "allele2_ccg", "allele2_reads", "allele2_fraction", "allele2_fragment_equivalent_cag",
                    "reads_total", "reads_pass_filters", "repeat_reads", "complete_reads", "notes", "flags"])
        for r in results:
            al = r.get("alleles", [])
            row = [r["sample"], r["genotype"]]
            for i in range(2):
                if i < len(al):
                    a = al[i]
                    row += [a["cag"], a["config"], a["ccg"] or "", a["support"], f"{a['fraction']:.3f}",
                            "" if a["fe_cag"] is None else a["fe_cag"]]
                else:
                    row += ["", "", "", "", "", ""]
            flags = "; ".join(f"{a['cag']}: {fl}" for a in al for fl in a.get("flags", []))
            row += [r["n_reads"], r["n_pass_filters"], r["n_repeat_reads"], r.get("n_complete", ""),
                    "; ".join(r.get("notes", [])), flags]
            w.writerow(row)
    return p
