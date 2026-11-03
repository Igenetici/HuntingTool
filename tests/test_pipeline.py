import json

import pytest

from huntingtool.fastq import discover_samples
from huntingtool.pipeline import Options, run
from tests.simulate import write_fastq

CASES = {
    # name: (alleles, expected sorted [(cag, config, ccg)])
    "het_canonical": ([(17, "canonical", 7, 1), (42, "canonical", 7, 1)], [(17, "canonical", 7), (42, "canonical", 7)]),
    "caa_loss": ([(19, "canonical", 7, 1), (27, "CAA-loss", 7, 1)], [(19, "canonical", 7), (27, "CAA-loss", 7)]),
    "dup": ([(15, "CAACAG-dup", 7, 1), (18, "canonical", 7, 1)], [(15, "CAACAG-dup", 7), (18, "canonical", 7)]),
    "one_apart_ccg": ([(17, "canonical", 7, 1), (18, "canonical", 10, 1)], [(17, "canonical", 7), (18, "canonical", 10)]),
    "caacagccgcca_loss": ([(18, "canonical", 7, 1), (44, "CAACAGCCGCCA-loss", 7, 1)],
                          [(18, "canonical", 7), (44, "CAACAGCCGCCA-loss", 7)]),
    "cca_loss": ([(17, "canonical", 7, 1), (41, "CCA-loss", 8, 1)], [(17, "canonical", 7), (41, "CCA-loss", 8)]),
    "imbalanced_expanded": ([(22, "canonical", 7, 1), (85, "canonical", 7, 0.4)], [(22, "canonical", 7), (85, "canonical", 7)]),
    "homozygous": ([(17, "canonical", 7, 1)], [(17, "canonical", 7)]),
}


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("err", [0.005, 0.03])
def test_genotypes(tmp_path, name, err):
    alleles, expected = CASES[name]
    fq = tmp_path / f"{name}.fastq.gz"
    write_fastq(str(fq), alleles, 2000, err=err, seed=7)
    res = run(discover_samples([str(fq)]), str(tmp_path / "out"), Options(html=False))
    got = [(a["cag"], a["config"], a["ccg"]) for a in res[0]["alleles"]]
    assert got == expected, (got, res[0]["notes"])


def test_barcode_folder_and_html(tmp_path):
    d = tmp_path / "fastq_pass"
    for b in ("barcode01", "barcode02"):
        (d / b).mkdir(parents=True)
        for c in range(2):
            write_fastq(str(d / b / f"chunk{c}.fastq.gz"), [(17, "canonical", 7, 1), (40, "canonical", 7, 1)], 300,
                        seed=c + 3)
    samples = discover_samples([str(d)])
    assert [s.name for s in samples] == ["barcode01", "barcode02"]
    assert all(len(s.files) == 2 for s in samples)
    res = run(samples, str(tmp_path / "out"), Options(), threads=2)
    assert (tmp_path / "out" / "report.html").exists()
    assert (tmp_path / "out" / "summary.tsv").exists()
    assert json.loads((tmp_path / "out" / "barcode01" / "result.json").read_text())["genotype"] == "17/40"
    assert res[1]["genotype"] == "17/40"


def test_motif_mode(tmp_path):
    fq = tmp_path / "m.fastq"
    write_fastq(str(fq), [(20, "canonical", 7, 1), (35, "canonical", 7, 1)], 800, err=0.005)
    res = run(discover_samples([str(fq)]), str(tmp_path / "o"), Options(mode="motif", html=False))
    # longest pure CAG run, both strands
    assert [a["cag"] for a in res[0]["alleles"]] == [20, 35]
