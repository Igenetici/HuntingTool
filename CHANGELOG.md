# Changelog

## 2.0.0 - 2026-11

Complete rewrite.

- Read-level parsing of the whole HTT exon 1 repeat region: uninterrupted CAG, repeat-adjacent sequence, CCG repeat.
- New configurations: CCA-loss and CAACAGCCGCCA-loss (CAA/CCA-loss); CAA-gain G2 reported as CAACAG-duplication; unclassified junctions listed.
- Nomenclature aligned with GeM-HD Consortium (2025): "WT" → canonical, "LOSS" → CAA-loss, "LOI/DOI" no longer used.
- Both read orientations analysed (v1 used forward-strand reads only).
- CAG length measured between anchors, tolerant to isolated basecalling errors.
- Automatic allele calling and phasing (configuration and CCG per allele); flags for unresolved alleles one repeat apart and for marked allelic imbalance.
- Optional mean read-quality filter; read-length and read-quality distributions reported.
- Command-line interface (`huntingtool run`) and redesigned graphical interface (`huntingtool` / `huntingtool-gui`), double-click launchers for Windows, macOS and Linux.
- MinKNOW `fastq_pass` folders supported (one sample per barcode); multi-sample parallel processing.
- Outputs: `summary.tsv`, per-sample JSON and haplotype counts, optional per-read table, self-contained HTML report.
- Generic motif mode for other repeat loci.
- Installable package (`pip install`), test suite with simulated reads, continuous integration.
- v1.0.0 kept as `legacy/huntingtool_v1.py`.

## 1.0.0 - 2026-06-12

- Initial public version.
