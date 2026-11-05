"""Command-line interface.

    huntingtool                      -> opens the graphical interface
    huntingtool gui                  -> opens the graphical interface
    huntingtool run -i IN [IN ...] -o OUTDIR [options]
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="huntingtool",
                                description="Sequence-resolved HTT CAG-repeat genotyping from amplicon long reads.")
    p.add_argument("--version", action="version", version=f"HuntingTool {__version__}")
    sub = p.add_subparsers(dest="cmd")
    sub.add_parser("gui", help="open the graphical interface")
    r = sub.add_parser("run", help="analyse FASTQ files or folders",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    r.add_argument("-i", "--input", nargs="+", required=True,
                   help="FASTQ(.gz) files and/or folders (a MinKNOW fastq_pass folder gives one sample per barcode)")
    r.add_argument("-o", "--outdir", required=True, help="output folder")
    r.add_argument("--mode", choices=["htt", "motif"], default="htt",
                   help="'htt': full structure parsing of HTT exon 1; 'motif': generic longest-run counting (other loci)")
    r.add_argument("--motif", default="CAG", help="repeat motif for --mode motif")
    r.add_argument("--min-length", type=int, default=100, help="minimum read length (bp)")
    r.add_argument("--max-length", type=int, default=0, help="maximum read length (bp, 0 = no limit)")
    r.add_argument("--min-q", type=float, default=0.0, help="minimum mean read Phred quality (0 = no filter)")
    r.add_argument("--min-repeats", type=int, default=5, help="minimum repeat count (motif mode)")
    r.add_argument("--min-allele-fraction", type=float, default=0.05,
                   help="minimum fraction of fully parsed reads supporting an allele")
    r.add_argument("--min-reads", type=int, default=10, help="minimum reads supporting an allele")
    r.add_argument("-t", "--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="samples analysed in parallel")
    r.add_argument("--merge", action="store_true", help="merge all inputs into a single sample")
    r.add_argument("--name", help="sample name when --merge is used")
    r.add_argument("--write-reads", action="store_true", help="write per-read results (reads.tsv.gz)")
    r.add_argument("--no-html", action="store_true", help="skip figures and HTML report")
    return p


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] == "gui":
        from .gui import main as gui_main
        return gui_main()
    args = build_parser().parse_args(argv)
    if args.cmd != "run":
        build_parser().print_help()
        return 1
    from .fastq import discover_samples
    from .pipeline import Options, run
    samples = discover_samples(args.input, merge=args.merge, name=args.name)
    if not samples:
        print("No FASTQ input found.", file=sys.stderr)
        return 1
    opt = Options(mode=args.mode, motif=args.motif.upper(), min_length=args.min_length,
                  max_length=args.max_length, min_q=args.min_q, min_repeats=args.min_repeats,
                  min_fraction=args.min_allele_fraction, min_reads=args.min_reads,
                  write_reads=args.write_reads, html=not args.no_html)
    print(f"HuntingTool {__version__}: {len(samples)} sample(s) -> {args.outdir}")
    run(samples, args.outdir, opt, threads=args.threads, progress=print)
    print(f"Summary: {os.path.join(args.outdir, 'summary.tsv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
