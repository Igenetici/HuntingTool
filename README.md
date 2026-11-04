# HuntingTool

**HuntingTool** genotypes the *HTT* exon 1 repeat region from amplicon long reads (Oxford Nanopore). For each read it measures the uninterrupted CAG repeat, identifies the repeat-adjacent sequence configuration and measures the CCG repeat. Because all three are read from the same molecule, each allele is reported with its own configuration and CCG length (phasing).

It runs from a **graphical interface** or from the **command line**, on Linux, macOS and Windows.

> Research use. Results must be interpreted by qualified personnel and, when clinically relevant, confirmed with validated diagnostic methods.

## What's new in version 2.0

| | v1.0 | v2.0 |
|---|---|---|
| Strands used | forward reads only | both strands (reads are oriented automatically) |
| CAG sizing | longest pure CAG run anywhere in the read | anchored between the 5′ flank and the CCG repeat; tolerant to isolated basecalling errors |
| Configurations | canonical, CAA-loss, eight CAA-gain positions | canonical, **CAA-loss, CCA-loss, CAACAGCCGCCA-loss**, CAACAG-duplication, CAA-gain G1 and G3–G8; novel sequences are listed |
| CCG repeat | not analysed | measured on every read |
| Genotype | modal peaks read from the plots | automatic allele calling with phasing of configuration and CCG length |
| Close alleles | not resolved | separated by CCG length or configuration; flagged when they differ by one CAG repeat |
| Quality filter | none in the tool | optional minimum mean read Q (e.g. Q20) |
| Input | files selected one by one | files, folders, or a MinKNOW `fastq_pass` folder (one sample per barcode) |
| Speed | single process, GUI-bound | streaming, parallel samples (about 80,000 reads/s per core on simulated data) |
| Interface | GUI only | GUI **and** command line |
| Output | TSV + PNG | `summary.tsv`, per-sample JSON and haplotype counts, optional per-read table, self-contained HTML report |

Version 1.0 is kept in `legacy/huntingtool_v1.py` for reproducibility.

## Installation

You need Python 3.9 or later.

**Option A – install the package (recommended)**

```bash
pip install "git+https://github.com/Igenetici/HuntingTool.git"
```

This installs two commands: `huntingtool` (command line, or the window if run without arguments) and `huntingtool-gui` (the window).

**Option B – without installing**

Download the repository (green *Code* button → *Download ZIP*), unzip it, and install the only dependency:

```bash
pip install matplotlib
```

Tkinter (the window toolkit) is included with Python on Windows and macOS. On Ubuntu/Debian install it with `sudo apt install python3-tk`.

## Starting HuntingTool

**Graphical interface**

* Installed: run `huntingtool` (or `huntingtool-gui`).
* Not installed: double-click the launcher in the `launchers/` folder:
  * Windows: `HuntingTool.bat`
  * macOS: `HuntingTool.command` (the first time, right-click → *Open*)
  * Linux: `HuntingTool.sh`; to add HuntingTool to the application menu run `bash launchers/install_linux_shortcut.sh` once.
* Or from a terminal: `python huntingtool.py`

In the window: add FASTQ files or folders, choose the output folder, and click **Run analysis**. When the analysis has finished, click **Open report**.

**Command line**

```bash
# one sample per barcode from a MinKNOW run
huntingtool run -i /path/to/run/fastq_pass -o results/

# individual files, Q20 filter, per-read table, 4 samples in parallel
huntingtool run -i sample1.fastq.gz sample2.fastq.gz -o results/ --min-q 20 --write-reads -t 4

# several files of the same sample
huntingtool run -i chunk1.fastq.gz chunk2.fastq.gz --merge --name patient01 -o results/

# other loci: generic longest-run counting of any motif on both strands
huntingtool run -i reads.fastq.gz -o results/ --mode motif --motif CAG
```

Without installing, replace `huntingtool` with `python huntingtool.py`. All options are listed by `huntingtool run --help`.

| Option | Default | Meaning |
|---|---|---|
| `--min-length` | 100 | minimum read length (bp) |
| `--max-length` | 0 | maximum read length (0 = no limit) |
| `--min-q` | 0 | minimum mean read Phred quality (0 = no filter) |
| `--min-allele-fraction` | 0.05 | minimum fraction of resolved reads supporting an allele |
| `--min-reads` | 10 | minimum reads supporting an allele |
| `-t / --threads` | CPUs − 1 | samples analysed in parallel |
| `--write-reads` | off | write `reads.tsv.gz` (one row per read) |
| `--mode` | `htt` | `htt` = structure-aware; `motif` = generic repeat counting |

## Output

```
results/
├── report.html               summary table and one section per sample (open in a browser)
├── summary.tsv               one row per sample: genotype, configuration and CCG per allele, read support
└── <sample>/
    ├── result.json           all metrics
    ├── haplotype_counts.tsv  reads per (configuration, CCG, CAG)
    ├── report.png            CAG distribution, CAG×CCG haplotypes, read lengths, parsing status
    └── reads.tsv.gz          (optional) per-read structure
```

Main fields per allele:

* **CAG**: uninterrupted CAG repeats, i.e. codons between the 5′ flank and the repeat-adjacent sequence.
* **Configuration**: the repeat-adjacent configuration (see below).
* **CCG**: CCG repeats.
* **Fragment-equivalent CAG**: the value expected from PCR fragment analysis with a primer ending in the CAA-CAG-CCG-CCA region. For example, CAA-loss gives +2 by sequencing and CAACAG-duplication gives −2.
* **Configuration purity**: the percentage of the allele's reads that carry its configuration.
* **Neighbour ratio**: reads at CAG ± 1 divided by reads at the mode. A high value (≥0.6) is flagged as a possible second allele differing by one repeat.

## Configurations

Nomenclature follows the GeM-HD Consortium (*Nat Genet* 2025).

| Name | Structure |
|---|---|
| canonical | (CAG)n CAACAGCCGCCA (CCG)m |
| CAA-loss | (CAG)n CCGCCA (CCG)m |
| CCA-loss | (CAG)n CAACAG (CCG)m |
| CAACAGCCGCCA-loss (CAA/CCA-loss) | (CAG)n (CCG)m |
| CAACAG-dup (= CAA-gain G2) | (CAG)n CAACAGCAACAGCCGCCA (CCG)m |
| CAA-gain Gk (k = 1, 3–8) | (CAG)n CAA(CAG)k−1 CAACAGCCGCCA (CCG)m |

## How it works

1. Each read is oriented on the *HTT* sense strand and anchored on the 12-nt flank immediately upstream of the CAG tract (CTCAAGTCCTTC; one mismatch tolerated), on the start of the CCG repeat, and on the downstream CCT sequence.
2. The sequence between the flank and the CCG repeat is split into the CAG tract and the repeat-adjacent sequence. The adjacent sequence is identified by exact matching of the known configurations, preceded by at least two CAG codons. The CAG length is the length of the tract divided by three and rounded, which tolerates isolated basecalling errors inside long tracts. Reads whose junction cannot be resolved exactly are counted but are not used to assign configurations.
3. Allele discovery is done separately within each repeat-adjacent configuration, so alleles with different configurations are never merged. Candidate alleles are the modes of the CAG distribution (raw counts up to 35 repeats; smoothed counts above 35, where distributions are broader). Each mode is scored by the reads within ±1 CAG (±2 above 35 CAG), so the broad left tail produced by PCR stutter and truncated products of a long allele does not inflate the support of a shorter allele.
4. An allele is called when its score reaches 5% of resolved reads (and at least 10 reads); at most two alleles are reported and additional modes are listed in the notes. A minor mode of a non-canonical configuration lying within 3 CAG of a mode at least five times larger is treated as a configuration miscall (a single basecalling error in the CAA or CCA codon) and is not called; the same applies to a minor canonical mode located k repeats above a CAA-gain allele (k = 2 for CAACAG-duplication), the product of a CAA>CAG error in its first CAA codon. Repeat-adjacent sequences matching no known configuration but carried by ≥20% of resolved reads are reported as an allele with an *unclassified* configuration.
5. The CCG length of each allele is the majority CCG length of its reads (flagged when supported by <60% of them). When a single CAG mode is found, it is split into two alleles if its reads show two CCG lengths at least two repeats apart, each with ≥25% of the reads (alleles of equal or adjacent CAG length distinguished by CCG). When two CAG modes are present no CCG split is attempted, because recombinant (chimeric) PCR products carry the CAG tract of one allele and the CCG tract of the other.
6. For alleles ≤35 CAG, a neighbouring bin reaching 60% of the mode is flagged as a possible second allele differing by one repeat.

On simulated amplicon reads (random indels and substitutions, both strands, unbalanced alleles), the exact genotype, configuration and CCG were recovered in 100% of heterozygous cases at 2% and 5% per-base error. At 8% error, most discrepancies were ±1 repeat on expanded alleles. Run `pytest` to reproduce the tests.

**Limitations**: the amplicon does not allow reliable genotyping of the CCT repeat (the reverse primer overlaps it). Alleles with the same CCG length and configuration that differ by a single CAG repeat may not be resolved; they are flagged. Thresholds are user-adjustable and should be verified in each laboratory.

## Citation

If you use HuntingTool, please cite the software (see `CITATION.cff` and the Zenodo DOI of the release) and the accompanying article.

## License

MIT – see `LICENSE`.
