import pytest
from pathlib import Path
from typer.testing import CliRunner

from aegis.cli.split import app, classify_feature

runner = CliRunner()


@pytest.fixture
def populus_test_files(tmp_path):
    """Create synthetic Populus tomentosa genome FASTA and GFF3 files matching user's case."""
    fasta_content = (
        ">CM031969.1 Populus tomentosa isolate GM15 chromosome 1D, whole genome shotgun sequence\n"
        "ATGCATGCATGCATGC\n"
        ">CM031970.1 Populus tomentosa isolate GM15 chromosome 2A, whole genome shotgun sequence\n"
        "GCATGCATGCATGCAT\n"
        ">CM031971.1 Populus tomentosa isolate GM15 chromosome 2D, whole genome shotgun sequence\n"
        "TTGCAATTCGATCGAT\n"
        ">CM031972.1 Populus tomentosa isolate GM15 chromosome 3A, whole genome shotgun sequence\n"
        "AATTGGCCAATTGGCC\n"
        ">scaffold_999 Populus tomentosa isolate GM15 unplaced scaffold\n"
        "CCCCCCCCCCCCCCCC\n"
    )
    fasta_file = tmp_path / "populus.fasta"
    fasta_file.write_text(fasta_content)

    gff_content = (
        "##gff-version 3\n"
        "CM031969.1\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene1D;Name=gene1D\n"
        "CM031969.1\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna1D;Parent=gene1D\n"
        "CM031970.1\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene2A;Name=gene2A\n"
        "CM031970.1\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna2A;Parent=gene2A\n"
        "CM031971.1\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene2D;Name=gene2D\n"
        "CM031971.1\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna2D;Parent=gene2D\n"
        "CM031972.1\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene3A;Name=gene3A\n"
        "CM031972.1\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna3A;Parent=gene3A\n"
        "scaffold_999\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene_unplaced;Name=gene_unplaced\n"
        "scaffold_999\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna_unplaced;Parent=gene_unplaced\n"
    )
    gff_file = tmp_path / "populus.gff3"
    gff_file.write_text(gff_content)

    return fasta_file, gff_file


# ---------------------------------------------------------------------------
# Unit tests for classify_feature logic
# ---------------------------------------------------------------------------

def test_classify_feature_smart():
    """Verify classify_feature logic across diverse genomic naming patterns."""
    tag, warn = classify_feature(
        "CM031970.1",
        "Populus tomentosa isolate GM15 chromosome 2A, whole genome shotgun sequence",
        split_tags=["A", "D"],
    )
    assert tag == "A"
    assert warn is None

    tag, warn = classify_feature(
        "CM031969.1",
        "Populus tomentosa isolate GM15 chromosome 1D, whole genome shotgun sequence",
        split_tags=["A", "D"],
    )
    assert tag == "D"
    assert warn is None

    tag, warn = classify_feature(
        "scaffold_999",
        "Populus tomentosa isolate GM15 unplaced scaffold",
        split_tags=["A", "D"],
    )
    assert tag is None

    tag, _ = classify_feature("chr1A", "chr1A sequence", split_tags=["A", "B", "D"])
    assert tag == "A"


def test_classify_feature_case_sensitive():
    """Verify case sensitivity: lowercase 'a' in 'tomentosa' does not match tag 'A'."""
    tag, _ = classify_feature(
        "CM031969.1",
        "Populus tomentosa isolate GM15 chromosome 1D",
        split_tags=["A", "B"],
        ignore_case=False,
    )
    assert tag is None


def test_classify_feature_sweet_potato_cultivar_prefix():
    """Verify cultivar prefix contig pattern matching."""
    tags = ["A", "B", "C", "D", "E", "F", "u"]

    tag_a, warn_a = classify_feature("BrgdChr01A", "BrgdChr01A", split_tags=tags)
    assert tag_a == "A"
    assert warn_a is None

    tag_b, warn_b = classify_feature("BrgdChr01B", "BrgdChr01B", split_tags=tags)
    assert tag_b == "B"
    assert warn_b is None


# ---------------------------------------------------------------------------
# CLI Smoke tests
# ---------------------------------------------------------------------------

def test_split_coupled_genome_and_annotation_smoke(populus_test_files, tmp_path):
    """Smoke test: ensure split CLI runs on coupled FASTA + GFF3."""
    fasta_file, gff_file = populus_test_files
    out_dir = tmp_path / "split_out"

    args = [
        str(fasta_file),
        str(gff_file),
        "--split-by", "A,D",
        "-d", str(out_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"
    assert (out_dir / "populus_splitA.fasta").exists()
    assert (out_dir / "populus_splitA.gff3").exists()
    assert (out_dir / "populus_splitD.fasta").exists()
    assert (out_dir / "populus_splitD.gff3").exists()


def test_split_with_regex_smoke(populus_test_files, tmp_path):
    """Smoke test: ensure splitting with a regex capture group works."""
    fasta_file, gff_file = populus_test_files
    out_dir = tmp_path / "regex_out"

    args = [
        str(fasta_file),
        str(gff_file),
        "--regex", r"chromosome \d+([AD])",
        "-d", str(out_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"
    assert (out_dir / "populus_splitA.fasta").exists()
    assert (out_dir / "populus_splitD.fasta").exists()

