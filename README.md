# HuntingTool

**HuntingTool** is a Python/Tkinter tool for exploratory analysis of **HTT CAG repeat expansions** and **CAA interruption patterns** from FASTQ or FASTQ.GZ sequencing reads.

Repository description for GitHub:

> Python/Tkinter tool for HTT CAG repeat sizing and CAA interruption analysis from FASTQ/FASTQ.GZ reads.

## Main features

- Reads one or more `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` files.
- Counts the longest contiguous repeat run for a selected motif, with `CAG` as default.
- Filters reads by minimum read length and minimum number of repeats.
- Generates a TSV file with repeat counts per read.
- Generates a histogram of repeat-size distribution.
- Evaluates predefined sequence patterns for CAA interruption loss near the CAG repeat.
- Evaluates predefined positional CAA gain patterns across eight positions.
- Produces summary plots and text reports for CAA interruption analysis.

## Intended use

This tool was developed for in-house/research analysis of sequencing reads in the context of HTT repeat expansion assessment. It is not intended to be used as a standalone diagnostic device. Results should be interpreted by qualified personnel and, when clinically relevant, confirmed using validated diagnostic methods.

## Requirements

- Python 3.9 or newer is recommended.
- `matplotlib`
- `tkinter` for the graphical file-selection interface. On Windows and macOS it is usually included with Python. On some Linux distributions it may need to be installed separately.

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR-USERNAME/HuntingTool.git
cd HuntingTool
```

Create and activate a virtual environment.

On Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

On macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the Python dependency:

```bash
pip install -r requirements.txt
```

On Linux, if `tkinter` is missing, install it with your system package manager. For example, on Ubuntu/Debian:

```bash
sudo apt install python3-tk
```

## Usage

Run:

```bash
python huntingtool.py
```

The graphical interface will ask you to:

1. Select one or more FASTQ/FASTQ.GZ files.
2. Enter the repeat motif, usually `CAG`.
3. Enter the minimum read length, default `500`.
4. Select where to save the TSV output file.

The additional output files are automatically saved next to the TSV file.

## Output files

If the selected TSV output is named `sample.tsv`, HuntingTool may generate:

- `sample.tsv`: repeat counts per filtered read.
- `sample.png`: histogram of repeat counts.
- `sample_CAA_interruption.txt`: CAA interruption loss/gain report.
- `sample_CAA_summary.png`: summary plot for LOSS and GAIN positions.
- `sample_CAA_distributions.png`: CAG repeat distributions for LOSS and GAIN-positive reads.

## Notes on interpretation

The current version uses predefined sequence patterns for CAA interruption loss and eight positional CAA gain configurations. The main CAG count is based on the maximum contiguous run of the selected motif observed in each read. Thresholds and interpretation labels are implemented in the script and should be reviewed before use in a new analytical setting.

## Suggested GitHub topics

`huntington-disease`, `htt`, `cag-repeat`, `caa-interruption`, `short-tandem-repeats`, `fastq`, `ngs`, `bioinformatics`, `python`

## License

This project is released under the MIT License. See `LICENSE` for details.
