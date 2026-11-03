"""FASTQ input: fast (optionally parallel) decompression, read iteration, sample discovery."""
from __future__ import annotations

import gzip
import io
import math
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

# per-character error probability for Phred+33 qualities
_ERR = [10 ** (-(i - 33) / 10) if i >= 33 else 1.0 for i in range(256)]


def mean_qscore(qual: str) -> float:
    """Mean read quality as reported by ONT tools: -10*log10(mean error probability)."""
    if not qual:
        return 0.0
    e = sum(map(_ERR.__getitem__, qual.encode("ascii", "replace"))) / len(qual)
    return -10 * math.log10(e) if e > 0 else 60.0


def is_fastq(p: Path) -> bool:
    return p.is_file() and p.name.lower().endswith(FASTQ_EXT)


def _open(path: Path) -> io.TextIOBase:
    if not str(path).lower().endswith(".gz"):
        return open(path, "rt", buffering=1 << 20)
    # external decompressors are 2-4x faster than the gzip module
    for tool in ("pigz", "igzip", "gzip"):
        exe = shutil.which(tool)
        if exe:
            proc = subprocess.Popen([exe, "-dc", str(path)], stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, bufsize=1 << 20)
            return io.TextIOWrapper(proc.stdout, encoding="ascii", errors="replace")
    try:  # python-isal, if installed
        from isal import igzip  # type: ignore
        return io.TextIOWrapper(igzip.open(path, "rb"), encoding="ascii", errors="replace")
    except ImportError:
        return gzip.open(path, "rt", encoding="ascii", errors="replace")


def iter_fastq(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (sequence, quality) for every record of a FASTQ(.gz) file."""
    with _open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            if not h.startswith("@"):
                continue
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            yield seq.upper(), qual


@dataclass
class Sample:
    name: str
    files: list[Path] = field(default_factory=list)


def discover_samples(inputs: list[str], merge: bool = False, name: str | None = None) -> list[Sample]:
    """Turn files/folders into samples.

    * a FASTQ file                      -> one sample (named after the file)
    * a folder with barcodeXX/ subfolders (MinKNOW output, e.g. fastq_pass/)
                                        -> one sample per subfolder
    * a folder containing FASTQ files   -> one sample per folder (all chunks merged)
    * --merge                           -> everything into a single sample
    """
    samples: list[Sample] = []
    for raw in inputs:
        p = Path(raw).expanduser()
        if is_fastq(p):
            nm = p.name
            for ext in FASTQ_EXT:
                if nm.lower().endswith(ext):
                    nm = nm[: -len(ext)]
                    break
            samples.append(Sample(nm, [p]))
        elif p.is_dir():
            subs = sorted(d for d in p.iterdir() if d.is_dir() and any(is_fastq(f) for f in d.iterdir()))
            direct = sorted(f for f in p.iterdir() if is_fastq(f))
            for d in subs:
                if d.name.lower() == "unclassified":
                    continue
                samples.append(Sample(d.name, sorted(f for f in d.iterdir() if is_fastq(f))))
            if direct:
                samples.append(Sample(p.name, direct))
        else:
            raise FileNotFoundError(f"Not a FASTQ file or folder: {raw}")
    if merge and samples:
        allf = [f for s in samples for f in s.files]
        samples = [Sample(name or samples[0].name, allf)]
    # make names unique
    seen: dict[str, int] = {}
    for s in samples:
        if s.name in seen:
            seen[s.name] += 1
            s.name = f"{s.name}_{seen[s.name]}"
        else:
            seen[s.name] = 0
    return samples


def file_size(paths: list[Path]) -> int:
    return sum(os.path.getsize(p) for p in paths)
