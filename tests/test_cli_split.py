import pytest
import os
from pathlib import Path
from typer.testing import CliRunner
from Bio import SeqIO

from aegis.cli.split import app, classify_feature
from aegis.cli.tidy_genome import app as tidy_genome_app
from aegis import Genome, Annotation

runner = CliRunner()


@pytest.fixture
def populus_test_files(tmp_path):
    """
    Create synthetic Populus tomentosa genome FASTA and GFF3 files matching user's case.
    """
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


def test_classify_feature_smart():
    """Verify classify_feature logic across diverse genomic naming patterns."""
    # Description with "chromosome 2A"
    tag, warn = classify_feature(
        "CM031970.1",
        "Populus tomentosa isolate GM15 chromosome 2A, whole genome shotgun sequence",
        split_tags=["A", "D"],
    )
    assert tag == "A"
    assert warn is None

    # Description with "chromosome 1D"
    tag, warn = classify_feature(
        "CM031969.1",
        "Populus tomentosa isolate GM15 chromosome 1D, whole genome shotgun sequence",
        split_tags=["A", "D"],
    )
    assert tag == "D"
    assert warn is None

    # Unplaced scaffold
    tag, warn = classify_feature(
        "scaffold_999",
        "Populus tomentosa isolate GM15 unplaced scaffold",
        split_tags=["A", "D"],
    )
    assert tag is None

    # ID directly contains haplotype / chromosome tag
    tag, _ = classify_feature("chr1A", "chr1A sequence", split_tags=["A", "B", "D"])
    assert tag == "A"

    tag, _ = classify_feature("chr1B", "chr1B sequence", split_tags=["A", "B", "D"])
    assert tag == "B"

    tag, _ = classify_feature("Chr01_hap1", "Haplotype 1", split_tags=["hap1", "hap2"])
    assert tag == "hap1"


def test_classify_feature_case_sensitive():
    """Verify case sensitivity: lowercase 'a' in 'tomentosa' does not match tag 'A'."""
    tag, _ = classify_feature(
        "CM031969.1",
        "Populus tomentosa isolate GM15 chromosome 1D",
        split_tags=["A", "B"],
        ignore_case=False,
    )
    # Since only 'A' and 'B' are tags, and description has 1D (no A or B), should be None
    assert tag is None


def test_split_coupled_genome_and_annotation(populus_test_files, tmp_path):
    """Test splitting coupled genome FASTA and GFF3 files based on FASTA descriptions."""
    fasta_file, gff_file = populus_test_files
    out_dir = tmp_path / "split_out"

    args = [
        str(fasta_file),
        str(gff_file),
        "--split-by", "A,D",
        "-d", str(out_dir),
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    # Check generated files
    assert (out_dir / "populus_splitA.fasta").exists()
    assert (out_dir / "populus_splitA.gff3").exists()
    assert (out_dir / "populus_splitD.fasta").exists()
    assert (out_dir / "populus_splitD.gff3").exists()
    assert (out_dir / "populus_split_other.fasta").exists()
    assert (out_dir / "populus_split_other.gff3").exists()

    # Validate split A genome contents
    g_a = Genome("pop_a", str(out_dir / "populus_splitA.fasta"), quiet=True)
    assert set(g_a.scaffolds.keys()) == {"CM031970.1", "CM031972.1"}

    # Validate split D genome contents
    g_d = Genome("pop_d", str(out_dir / "populus_splitD.fasta"), quiet=True)
    assert set(g_d.scaffolds.keys()) == {"CM031969.1", "CM031971.1"}

    # Validate split other genome contents
    g_other = Genome("pop_other", str(out_dir / "populus_split_other.fasta"), quiet=True)
    assert set(g_other.scaffolds.keys()) == {"scaffold_999"}

    # Validate split A annotation contents
    a_a = Annotation(str(out_dir / "populus_splitA.gff3"), quiet=True)
    assert set(a_a.chrs.keys()) == {"CM031970.1", "CM031972.1"}
    assert set(a_a.all_gene_ids.keys()) == {"gene2A", "gene3A"}

    # Validate split D annotation contents
    a_d = Annotation(str(out_dir / "populus_splitD.gff3"), quiet=True)
    assert set(a_d.chrs.keys()) == {"CM031969.1", "CM031971.1"}
    assert set(a_d.all_gene_ids.keys()) == {"gene1D", "gene2D"}

    # Validate split other annotation contents
    a_other = Annotation(str(out_dir / "populus_split_other.gff3"), quiet=True)
    assert set(a_other.chrs.keys()) == {"scaffold_999"}
    assert set(a_other.all_gene_ids.keys()) == {"gene_unplaced"}


def test_split_genome_only_with_keep_description(populus_test_files, tmp_path):
    """Test splitting only a genome file, verifying --keep-description preserves header."""
    fasta_file, _ = populus_test_files
    out_dir = tmp_path / "genome_split_out"

    args = [
        "-g", str(fasta_file),
        "-s", "A,D",
        "-d", str(out_dir),
        "--keep-description",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_fasta_a = out_dir / "populus_splitA.fasta"
    assert out_fasta_a.exists()

    # Read back records to verify full description is preserved
    records = list(SeqIO.parse(str(out_fasta_a), "fasta"))
    assert len(records) == 2
    assert "chromosome 2A" in records[0].description


def test_split_annotation_only(populus_test_files, tmp_path):
    """Test splitting only an annotation file."""
    # Create GFF with explicit chromosome names like chr1A, chr1B
    gff_content = (
        "##gff-version 3\n"
        "chr1A\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene1A\n"
        "chr1A\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna1A;Parent=gene1A\n"
        "chr1B\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene1B\n"
        "chr1B\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna1B;Parent=gene1B\n"
        "scaffold_un\tGenBank\tgene\t1\t10\t.\t+\t.\tID=gene_un\n"
        "scaffold_un\tGenBank\tmRNA\t1\t10\t.\t+\t.\tID=rna_un;Parent=gene_un\n"
    )
    gff_file = tmp_path / "polyploid.gff3"
    gff_file.write_text(gff_content)

    out_dir = tmp_path / "annot_only_out"
    args = [
        "-a", str(gff_file),
        "--split-by", "A,B",
        "-d", str(out_dir),
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    assert (out_dir / "polyploid_splitA.gff3").exists()
    assert (out_dir / "polyploid_splitB.gff3").exists()
    assert (out_dir / "polyploid_split_other.gff3").exists()

    a_a = Annotation(str(out_dir / "polyploid_splitA.gff3"), quiet=True)
    assert list(a_a.chrs.keys()) == ["chr1A"]
    assert list(a_a.all_gene_ids.keys()) == ["gene1A"]


def test_split_positional_reverse_order(populus_test_files, tmp_path):
    """Test positional arguments when passed as (gff, fasta) instead of (fasta, gff)."""
    fasta_file, gff_file = populus_test_files
    out_dir = tmp_path / "rev_out"

    args = [
        str(gff_file),
        str(fasta_file),
        "--split-by", "A,D",
        "-d", str(out_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    assert (out_dir / "populus_splitA.fasta").exists()
    assert (out_dir / "populus_splitA.gff3").exists()


def test_split_with_regex(populus_test_files, tmp_path):
    """Test splitting with a regular expression capture group."""
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
    assert (out_dir / "populus_split_other.fasta").exists()


def test_split_with_split_map(populus_test_files, tmp_path):
    """Test splitting with a TSV map file."""
    fasta_file, gff_file = populus_test_files
    map_file = tmp_path / "mapping.tsv"
    map_file.write_text(
        "CM031969.1\thapD\n"
        "CM031970.1\thapA\n"
        "CM031971.1\thapD\n"
        "CM031972.1\thapA\n"
    )
    out_dir = tmp_path / "map_out"

    args = [
        str(fasta_file),
        str(gff_file),
        "--split-map", str(map_file),
        "-d", str(out_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    assert (out_dir / "populus_splithapA.fasta").exists()
    assert (out_dir / "populus_splithapD.fasta").exists()
    assert (out_dir / "populus_split_other.fasta").exists()


def test_tidy_genome_keep_description(populus_test_files, tmp_path):
    """Test that tidy-genome CLI also supports and honors --keep-description."""
    fasta_file, _ = populus_test_files
    out_dir = tmp_path / "tidy_out"

    # 1. Default (no keep-description): header contains scaffold ID
    args_default = [
        str(fasta_file),
        "-d", str(out_dir / "default"),
    ]
    res1 = runner.invoke(tidy_genome_app, args_default)
    assert res1.exit_code == 0
    records1 = list(SeqIO.parse(str(out_dir / "default/populus_tidy.fasta"), "fasta"))
    assert records1[0].description == "CM031969.1"

    # 2. With --keep-description: header retains full description
    args_keep = [
        str(fasta_file),
        "-d", str(out_dir / "keep"),
        "--keep-description",
    ]
    res2 = runner.invoke(tidy_genome_app, args_keep)
    assert res2.exit_code == 0
    records2 = list(SeqIO.parse(str(out_dir / "keep/populus_tidy.fasta"), "fasta"))
    assert "chromosome 1D" in records2[0].description


def test_split_with_punctuation_and_jaawwd(tmp_path):
    """
    Test splitting with 'A,' and 'D,' punctuation patterns to avoid ambiguity
    with accessions like JAAWWD010000446.1 that contain both 'A' and 'D'.
    """
    fasta_content = (
        ">CM031969.1 Populus tomentosa isolate GM15 chromosome 1D, whole genome shotgun sequence\n"
        "ATGCATGCATGCATGC\n"
        ">CM031970.1 Populus tomentosa isolate GM15 chromosome 2A, whole genome shotgun sequence\n"
        "GCATGCATGCATGCAT\n"
        ">JAAWWD010000446.1 Populus tomentosa isolate GM15 unplaced scaffold JAAWWD010000446.1\n"
        "TTTTTTTTTTTTTTTT\n"
    )
    fasta_file = tmp_path / "genome.fasta"
    fasta_file.write_text(fasta_content)

    # 1. Test with multiple -s options: -s "A," -s "D,"
    out1 = tmp_path / "out1"
    args1 = [str(fasta_file), "-s", "A,", "-s", "D,", "-d", str(out1), "-q"]
    res1 = runner.invoke(app, args1)
    assert res1.exit_code == 0, f"Error: {res1.stdout}"
    assert (out1 / "genome_splitA.fasta").exists()
    assert (out1 / "genome_splitD.fasta").exists()
    assert (out1 / "genome_split_other.fasta").exists()

    g_a1 = Genome("g_a1", str(out1 / "genome_splitA.fasta"), quiet=True)
    assert list(g_a1.scaffolds.keys()) == ["CM031970.1"]

    g_d1 = Genome("g_d1", str(out1 / "genome_splitD.fasta"), quiet=True)
    assert list(g_d1.scaffolds.keys()) == ["CM031969.1"]

    g_other1 = Genome("g_other1", str(out1 / "genome_split_other.fasta"), quiet=True)
    assert list(g_other1.scaffolds.keys()) == ["JAAWWD010000446.1"]

    # 2. Test with semicolon delimiter: -s "A,; D,"
    out2 = tmp_path / "out2"
    args2 = [str(fasta_file), "-s", "A,; D,", "-d", str(out2), "-q"]
    res2 = runner.invoke(app, args2)
    assert res2.exit_code == 0
    assert (out2 / "genome_splitA.fasta").exists()
    assert (out2 / "genome_splitD.fasta").exists()
    assert (out2 / "genome_split_other.fasta").exists()

    # 3. Test with explicit label mapping: -s "A,:A; D,:D"
    out3 = tmp_path / "out3"
    args3 = [str(fasta_file), "-s", "A,:A; D,:D", "-d", str(out3), "-q"]
    res3 = runner.invoke(app, args3)
    assert res3.exit_code == 0
    assert (out3 / "genome_splitA.fasta").exists()
    assert (out3 / "genome_splitD.fasta").exists()
    assert (out3 / "genome_split_other.fasta").exists()

    # 4. Test with double quotes and single quotes around A,;D,
    out4 = tmp_path / "out4"
    args4 = [str(fasta_file), "--split-by", '"A,;D,"', "-d", str(out4), "-q"]
    res4 = runner.invoke(app, args4)
    assert res4.exit_code == 0
    assert (out4 / "genome_splitA.fasta").exists()
    assert (out4 / "genome_splitD.fasta").exists()

    out5 = tmp_path / "out5"
    args5 = [str(fasta_file), "--split-by", "'A,;D,'", "-d", str(out5), "-q"]
    res5 = runner.invoke(app, args5)
    assert res5.exit_code == 0
    assert (out5 / "genome_splitA.fasta").exists()
    assert (out5 / "genome_splitD.fasta").exists()

    # 5. Test with unquoted double-comma A,,D, and slash A,/D,
    out6 = tmp_path / "out6"
    args6 = [str(fasta_file), "--split-by", "A,,D,", "-d", str(out6), "-q"]
    res6 = runner.invoke(app, args6)
    assert res6.exit_code == 0
    assert (out6 / "genome_splitA.fasta").exists()
    assert (out6 / "genome_splitD.fasta").exists()

    out7 = tmp_path / "out7"
    args7 = [str(fasta_file), "--split-by", "A,/D,", "-d", str(out7), "-q"]
    res7 = runner.invoke(app, args7)
    assert res7.exit_code == 0
    assert (out7 / "genome_splitA.fasta").exists()
    assert (out7 / "genome_splitD.fasta").exists()

    # 6. Test that standard -s "A,D" also does NOT falsely match JAAWWD
    out8 = tmp_path / "out8"
    args8 = [str(fasta_file), "-s", "A,D", "-d", str(out8), "-q"]
    res8 = runner.invoke(app, args8)
    assert res8.exit_code == 0
    g_other8 = Genome("g_other8", str(out8 / "genome_split_other.fasta"), quiet=True)
    assert "JAAWWD010000446.1" in g_other8.scaffolds

