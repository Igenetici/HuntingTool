"""Draws docs/workflow.png (HuntingTool v2 analysis workflow)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BLUE, INK, MUT, LINE = "#3B6FB6", "#1D2329", "#5D6873", "#C9CED4"
steps = [
    ("Input", "FASTQ / FASTQ.GZ files or MinKNOW fastq_pass folder\n(one sample per barcode)"),
    ("Read filters", "minimum read length (default 100 bp); optional minimum mean read Q"),
    ("Orientation", "reads oriented on the HTT sense strand (CAG vs CTG triplets):\nboth strands are used"),
    ("Anchoring", "5′ flank CTCAAGTCCTTC  ·  start of the CCG repeat  ·  CCT anchor"),
    ("Read structure", "(CAG)n  |  repeat-adjacent sequence  |  (CCG)m\nCAG length = anchor-to-anchor distance / 3 (tolerant to isolated errors)"),
    ("Configuration", "canonical · CAA-loss · CCA-loss · CAACAGCCGCCA-loss ·\nCAACAG-duplication (G2) · CAA-gain G1, G3–G8 · unresolved"),
    ("Allele discovery", "modes of the CAG distribution of reads with a resolved configuration\n(minimum 10% of reads per allele)"),
    ("Phasing", "configuration and CCG length of each allele by majority of its reads;\nsplit when two configurations or CCG lengths are both supported"),
    ("Output", "summary.tsv · result.json · haplotype counts · HTML report with\nCAG, CAG×CCG, read-length and parsing-status plots"),
]
fig, ax = plt.subplots(figsize=(7.2, 9.6), dpi=200)
ax.set_xlim(0, 10); ax.set_ylim(0, len(steps) * 1.1 + 0.2); ax.axis("off")
for i, (t, d) in enumerate(steps):
    y = len(steps) * 1.1 - i * 1.1 - 0.75
    ax.add_patch(FancyBboxPatch((0.2, y), 9.6, 0.82, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="white", ec=BLUE if i in (4, 5, 7) else LINE, lw=1.4))
    ax.text(0.45, y + 0.41, t, fontsize=9.5, fontweight="bold", color=INK, va="center")
    ax.text(2.75, y + 0.41, d, fontsize=7.6, color=MUT, va="center", linespacing=1.35)
    if i < len(steps) - 1:
        ax.annotate("", xy=(5, y - 0.26), xytext=(5, y + 0.0), arrowprops=dict(arrowstyle="-|>", color=MUT, lw=1))
fig.tight_layout()
fig.savefig("docs/workflow.png")
fig.savefig("docs/workflow.pdf")
