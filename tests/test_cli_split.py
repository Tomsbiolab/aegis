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



def test_split_by_chromosome_sends_listed_organelles_to_other_smoke(tmp_path):
    gff_file = tmp_path / "o.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chrM\ttest\tgene\t1\t9\t.\t+\t.\tID=g2\n"
        "chr9\ttest\tgene\t1\t9\t.\t+\t.\tID=g3\n"
    )
    out_dir = tmp_path / "split_out"
    result = runner.invoke(app, [str(gff_file), "--by-chromosome", "--mitochondria-chroms", "chr9", "-d", str(out_dir), "-q"])
    assert result.exit_code == 0, f"Error: {result.stdout}"
    assert sorted(f.name for f in out_dir.iterdir()) == ["o_split_other.gff3", "o_splitchr1.gff3"]


def test_split_paired_output_file_extensions(tmp_path):
    """Paired -o in split keeps distinct tags and proper GFF/FASTA extensions, including dotted prefixes."""
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1A\n" + "ACGT" * 25 + "\n>chr1D\n" + "TGCA" * 25 + "\n")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1A\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1A\tt\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1D\tt\tgene\t1\t30\t.\t+\t.\tID=g2\n"
        "chr1D\tt\tmRNA\t1\t30\t.\t+\t.\tID=t2;Parent=g2\n"
    )

    out1 = tmp_path / "out_dot"
    res1 = runner.invoke(app, [str(gff), str(fa), "-s", "A,D", "-d", str(out1), "-o", "release_v1.2", "-q"])
    assert res1.exit_code == 0, res1.output
    assert sorted(f.name for f in out1.iterdir()) == [
        "release_v1.2_splitA.fasta",
        "release_v1.2_splitA.gff3",
        "release_v1.2_splitD.fasta",
        "release_v1.2_splitD.gff3",
    ]

    out2 = tmp_path / "out_fa"
    res2 = runner.invoke(app, [str(gff), str(fa), "-s", "A,D", "-d", str(out2), "-o", "out.fasta", "-q"])
    assert res2.exit_code == 0, res2.output
    assert sorted(f.name for f in out2.iterdir()) == [
        "out_splitA.fasta",
        "out_splitA.gff3",
        "out_splitD.fasta",
        "out_splitD.gff3",
    ]
    assert (out2 / "out_splitA.fasta").read_text().startswith(">chr1A")
    assert (out2 / "out_splitD.fasta").read_text().startswith(">chr1D")
    assert (out2 / "out_splitA.gff3").read_text().startswith("##gff-version 3")


def test_split_write_empty_other_only_for_other_and_creates_fasta(tmp_path):
    """--write-empty-other writes empty _split_other FASTA/GFF3 and skips unmatched named tags."""
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1A\n" + "ACGT" * 25 + "\n")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1A\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1A\tt\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
    )
    out = tmp_path / "out"
    res = runner.invoke(app, [str(fa), str(gff), "-s", "A,Z", "--write-empty-other", "-d", str(out), "-q"])
    assert res.exit_code == 0, res.output
    assert sorted(f.name for f in out.iterdir()) == [
        "a_splitA.gff3",
        "a_split_other.gff3",
        "g_splitA.fasta",
        "g_split_other.fasta",
    ]
    assert (out / "g_split_other.fasta").read_text() == ""
