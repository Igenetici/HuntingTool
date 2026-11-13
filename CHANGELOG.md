# Changelog

## 2.2.0 - 2026-11

Stutter-aware sizing of expanded alleles.

- In expanded alleles the n−1 PCR-stutter product is almost as abundant as the full-length product, and reverse-strand reads (not analysed by v1.0.0) are on average slightly shorter than forward-strand reads; in v2.1.0 the most frequent CAG count therefore fell one repeat below the fragment-analysis size in about half of the expanded alleles. For alleles of 36 CAG or more, the allele is now assigned the longest CAG count, within 2 repeats of the most frequent count, that is supported by at least 80% of the reads of the most frequent count.
- The most frequent count is reported as `mode_cag`; alleles whose call differs from it are flagged ("stutter-aware call ... repeat count uncertain by 1 repeat").
- Alleles of 35 CAG or fewer are called as in v2.1.0. `stutter_ratio=0` in `call_alleles` restores the v2.1.0 behaviour.
- The 80% threshold was set on the study datasets (results unchanged between 70% and 80%); it depends on the PCR protocol and should be verified when the protocol changes.

## 2.1.0 - 2026-11

Allele caller revised after re-analysis of the study datasets (R10.4.1, 63 barcodes).

- CAG modes are searched separately within each repeat-adjacent configuration and scored on the reads within ±1 CAG (±2 above 35 CAG), so the broad stutter/truncation tail of a long allele no longer inflates or masks a neighbouring short allele.
- No CCG-based split when two CAG modes are present: recombinant (chimeric) PCR products carrying the CAG tract of one allele and the CCG tract of the other no longer produce spurious third alleles. A single CAG mode is still split when two CCG lengths ≥2 repeats apart are both supported.
- Default minimum allele fraction lowered from 10% to 5% of resolved reads; low-height modes are retained as candidates, so under-amplified long alleles (e.g. >80 CAG) are called.
- Minor non-canonical modes within 3 CAG of a ≥5-fold larger mode are treated as configuration miscalls; likewise a minor canonical mode located k repeats above a CAA-gain allele (CAACAG-duplication: k = 2), which is produced by a single CAA>CAG basecalling error in its first CAA codon.
- Repeat-adjacent sequences that match no known configuration but are carried by ≥20% of resolved reads are reported as alleles with an "unclassified" configuration (sequence shown).
- New note when the CAG mode differs between CCG classes of a single called allele (evidence for two alleles one repeat apart).
- Neighbour-bin flag restricted to alleles ≤35 CAG; CCG call flagged when supported by <60% of the allele's reads.
- New unit tests reproducing the read-count patterns observed in real data.

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
