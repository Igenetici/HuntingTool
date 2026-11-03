"""Figures (PNG) and a self-contained HTML report."""
from __future__ import annotations

import base64
import html
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import __version__  # noqa: E402
from .structure import CONFIG_ORDER  # noqa: E402

COLORS = {
    "canonical": "#3B6FB6", "CAA-loss": "#C8424F", "CCA-loss": "#E39B2D", "CAACAGCCGCCA-loss": "#7B3C8C",
    "CAACAG-dup": "#21917F", "other": "#A3AAB1",
}
_GAIN_COL = ["#6E9A3A", "#9C6B3E", "#4F8FA8", "#B05C8A", "#7D7D3B", "#5A6FA0", "#A0703A"]


def color(cfg: str) -> str:
    if cfg in COLORS:
        return COLORS[cfg]
    if cfg.startswith("CAA-gain G"):
        return _GAIN_COL[(int(cfg[-1]) - 1) % len(_GAIN_COL)]
    return "#A3AAB1"


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#E6E8EB", lw=0.6)
    ax.set_axisbelow(True)


def amplicon_size(cag: int, adj_len: int, ccg: int) -> int:
    """Expected amplicon size for the reference primer pair (51 nt 5' of the CAG tract,
    23 nt from the CCT repeat to the end of the reverse primer)."""
    return 51 + 3 * cag + adj_len + 3 * ccg + 23


def sample_figures(res: dict, out: Path) -> dict:
    figs = {}
    if res["options"]["mode"] == "motif":
        mc = Counter({int(k): v for k, v in res["motif_counts"].items()})
        fig, ax = plt.subplots(figsize=(8, 3.6), dpi=130)
        if mc:
            xs = sorted(mc)
            ax.bar(xs, [mc[x] for x in xs], color="#3B6FB6", width=0.85)
        for a in res["alleles"]:
            ax.axvline(a["cag"], ls="--", color="#333", lw=0.8)
        ax.set_xlabel(f"Longest uninterrupted {res['options']['motif']} run (repeats)")
        ax.set_ylabel("Reads")
        ax.set_title(res["sample"], loc="left", fontsize=10)
        _style(ax)
        fig.tight_layout()
        p = out / "repeat_distribution.png"
        fig.savefig(p)
        plt.close(fig)
        figs["main"] = str(p)
        return figs

    hap = res["hap_counts"]
    by_cfg: dict[str, Counter] = defaultdict(Counter)
    grid: Counter = Counter()
    for cfg, ccg, cag, k in hap:
        by_cfg[cfg][cag] += k
        if cfg != "other":
            grid[(cag, ccg)] += k
    if grid:
        floor = max(3, 0.002 * sum(grid.values()))
        grid = Counter({k: v for k, v in grid.items() if v >= floor})
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), dpi=130)
    # A: CAG distribution stacked by configuration
    ax = axs[0, 0]
    allx = sorted({x for d in by_cfg.values() for x in d}) or [0]
    xs = list(range(min(allx), max(allx) + 1))
    bottom = [0] * len(xs)
    cfgs = [c for c in CONFIG_ORDER if c in by_cfg] + [c for c in by_cfg if c not in CONFIG_ORDER]
    for c in cfgs:
        ys = [by_cfg[c][x] for x in xs]
        ax.bar(xs, ys, bottom=bottom, color=color(c), width=0.85, label=f"{c} ({sum(ys)})")
        bottom = [b + y for b, y in zip(bottom, ys)]
    ax.axvspan(35.5, 39.5, color="#F4C7C9", alpha=0.45, lw=0)
    ax.axvspan(39.5, max(xs[-1], 40) + 0.5, color="#E8999E", alpha=0.35, lw=0)
    for a in res["alleles"]:
        ax.axvline(a["cag"], ls="--", color="#333", lw=0.8)
        ax.text(a["cag"], max(bottom) * 1.02 if bottom else 1, str(a["cag"]), ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Uninterrupted CAG repeats (complete reads; grey = adjacent sequence not resolved)")
    ax.set_ylabel("Reads")
    ax.set_title("A  CAG distribution by repeat-adjacent configuration", loc="left", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)
    # B: CAG x CCG
    ax = axs[0, 1]
    if grid:
        cags = [k[0] for k in grid]
        ccgs = [k[1] for k in grid]
        vals = list(grid.values())
        mx = max(vals)
        ax.scatter(cags, ccgs, s=[8 + 220 * v / mx for v in vals], c="#3B6FB6", alpha=0.6, edgecolors="white", lw=0.5)
        for a in res["alleles"]:
            if a["ccg"]:
                ax.scatter([a["cag"]], [a["ccg"]], s=120, facecolors="none", edgecolors="#C8424F", lw=1.4)
        ax.set_ylim(min(ccgs) - 1, max(ccgs) + 1)
    ax.set_xlabel("Uninterrupted CAG repeats")
    ax.set_ylabel("CCG repeats")
    ax.set_title("B  Resolved haplotypes (marker area = reads; red = called alleles)", loc="left", fontsize=9)
    _style(ax)
    # C: read lengths
    ax = axs[1, 0]
    la = dict(res["read_length_hist"])
    lr = dict(res["repeat_read_length_hist"])
    if la:
        xs2 = range(0, min(max(la) + 2, 101))
        ax.bar([x * 10 for x in xs2], [la.get(x, 0) for x in xs2], width=9, color="#C9CED4", label="all reads")
        ax.bar([x * 10 for x in xs2], [lr.get(x, 0) for x in xs2], width=9, color="#3B6FB6", label="repeat reads")
    for a in res["alleles"]:
        if a["ccg"]:
            adj = 12 if a["config"] == "canonical" else None
            if adj is not None:
                ax.axvline(amplicon_size(a["cag"], adj, a["ccg"]), ls="--", color="#333", lw=0.8)
    ax.set_xlabel("Read length (bp; dashed = expected amplicon size of called canonical alleles)")
    ax.set_ylabel("Reads")
    ax.set_title("C  Read-length distribution", loc="left", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)
    # D: parsing status and configuration counts
    ax = axs[1, 1]
    st = res.get("status", {})
    labels = ["complete", "incomplete_ccg", "no_ccg", "no_flank"]
    vals = [st.get(k, 0) for k in labels]
    ax.barh(labels[::-1], vals[::-1], color="#3B6FB6")
    for i, v in enumerate(vals[::-1]):
        ax.text(v, i, f" {v}", va="center", fontsize=7)
    ax.set_xlabel("Repeat-containing reads")
    ax.set_title("D  Read parsing status", loc="left", fontsize=9)
    _style(ax)
    fig.suptitle(f"{res['sample']}  –  genotype {res['genotype']}", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    p = out / "report.png"
    fig.savefig(p)
    plt.close(fig)
    figs["main"] = str(p)
    return figs


def _img(path: str) -> str:
    try:
        b = base64.b64encode(Path(path).read_bytes()).decode()
        return f'<img src="data:image/png;base64,{b}" alt="figure">'
    except OSError:
        return ""


CSS = """
:root{--bg:#fff;--fg:#1d2329;--mut:#5d6873;--line:#e3e6ea;--acc:#3B6FB6;--warn:#fff4e0}
@media (prefers-color-scheme:dark){:root{--bg:#15191d;--fg:#e8ecef;--mut:#9aa6b1;--line:#2c333a;--acc:#7aa7e6;--warn:#3a3020}}
body{background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0;padding:24px;max-width:1200px}
h1{font-size:20px;margin:0 0 4px} h2{font-size:16px;margin:28px 0 8px;border-top:1px solid var(--line);padding-top:18px}
.mut{color:var(--mut)} table{border-collapse:collapse;width:100%;margin:8px 0;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}
th{color:var(--mut);font-weight:600} img{max-width:100%;background:#fff;border-radius:6px}
.note{background:var(--warn);padding:6px 10px;border-radius:6px;margin:4px 0;font-size:13px}
a{color:var(--acc)} code{font-size:12px}
"""


def write_html(results: list[dict], outdir: str, opt) -> Path:
    o = Path(outdir)
    rows = []
    for r in results:
        al = r.get("alleles", [])
        desc = "; ".join(f"{a['cag']} {a['config']}" + (f" CCG{a['ccg']}" if a.get('ccg') else "") +
                         f" ({a['support']} reads, {100*a['fraction']:.0f}%)" for a in al)
        warn = " ⚠" if r.get("notes") or any(a.get("flags") for a in al) else ""
        rows.append(f"<tr><td><a href='#{html.escape(r['sample'])}'>{html.escape(r['sample'])}</a></td>"
                    f"<td><b>{html.escape(r['genotype'])}</b>{warn}</td><td>{html.escape(desc)}</td>"
                    f"<td>{r['n_reads']}</td><td>{r['n_repeat_reads']}</td><td>{r.get('n_complete','')}</td></tr>")
    secs = []
    for r in results:
        al = r.get("alleles", [])
        at = "".join(
            f"<tr><td>{a['cag']}</td><td>{html.escape(a['config'])}</td><td>{a.get('ccg') or ''}</td>"
            f"<td>{a['support']} ({100*a['fraction']:.0f}%)</td><td>{'' if a.get('fe_cag') is None else a['fe_cag']}</td>"
            f"<td>{a.get('config_purity',0):.1f}%</td><td>{a.get('neighbour_ratio',0):.2f}</td></tr>" for a in al)
        notes = "".join(f"<div class='note'>{html.escape(n)}</div>" for n in r.get("notes", []))
        notes += "".join(f"<div class='note'>Allele {a['cag']}: {html.escape(f)}</div>" for a in al for f in a.get("flags", []))
        cp = r.get("canonical_pct") or {}
        cpt = ", ".join(f"{k}: {v:.1f}%" for k, v in sorted(cp.items())) or "–"
        oth = ", ".join(f"<code>{html.escape(s) or '∅'}</code> ({n})" for s, n in r.get("other_adjacent", [])[:5]) or "–"
        fig = _img(r.get("figures", {}).get("main", ""))
        secs.append(f"""<h2 id="{html.escape(r['sample'])}">{html.escape(r['sample'])} — {html.escape(r['genotype'])}</h2>
<div class='mut'>{len(r['files'])} file(s); reads {r['n_reads']}, passing filters {r['n_pass_filters']},
repeat-containing {r['n_repeat_reads']}, fully parsed {r.get('n_complete','–')}; strand {r.get('strand',{})}; {r['runtime_s']} s</div>
{notes}
<table><tr><th>CAG (uninterrupted)</th><th>Configuration</th><th>CCG</th><th>Reads</th><th>Fragment-equivalent CAG</th>
<th>Configuration purity</th><th>Neighbour ratio</th></tr>{at}</table>
<div class='mut'>Canonical percentage canonical/(canonical+X) [v1 metric]: {cpt}<br>Unclassified adjacent sequences: {oth}</div>
{fig}""")
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HuntingTool report</title><style>{CSS}</style></head><body>
<h1>HuntingTool v{__version__} report</h1>
<div class='mut'>{datetime.now():%Y-%m-%d %H:%M} · mode {opt.mode} · min length {opt.min_length} bp · min mean Q {opt.min_q} ·
minimum allele fraction {opt.min_fraction:.2f}. For research use; results must be interpreted by qualified personnel.</div>
<h2>Summary</h2><table><tr><th>Sample</th><th>Genotype (CAG)</th><th>Alleles</th><th>Reads</th><th>Repeat reads</th><th>Parsed</th></tr>
{''.join(rows)}</table>{''.join(secs)}</body></html>"""
    p = o / "report.html"
    p.write_text(page, encoding="utf-8")
    return p
