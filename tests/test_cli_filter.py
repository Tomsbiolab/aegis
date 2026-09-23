import pytest
from typer.testing import CliRunner

from aegis.cli.filter import app as filter_app
from aegis.cli.tidy import app as tidy_app
from aegis.annotation import Annotation

runner = CliRunner()


def test_filter_coding_only(rich_gff3_file, tmp_path):
    output_dir = tmp_path / "filter_out"
    output_file = "coding_only.gff3"

    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "--coding-only",
        "-q",
    ]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_gff = output_dir / output_file
    assert out_gff.exists()

    annot = Annotation(str(out_gff), quiet=True)
    # geneR2 is pure coding -> kept
    assert "geneR2" in annot.all_gene_ids
    # geneR3 is pure non-coding -> removed
    assert "geneR3" not in annot.all_gene_ids
    # geneR1 was mixed -> kept, but only has mRNA transcript now
    assert "geneR1" in annot.all_gene_ids
    gene = annot.chrs["chr1"]["geneR1"]
    for t in gene.transcripts.values():
        assert t.coding is True
        assert t.feature == "mRNA"


def test_filter_non_coding_only(rich_gff3_file, tmp_path):
    output_dir = tmp_path / "filter_out"
    output_file = "non_coding_only.gff3"

    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "--non-coding-only",
        "-q",
    ]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_gff = output_dir / output_file
    assert out_gff.exists()

    annot = Annotation(str(out_gff), quiet=True)
    # geneR2 is pure coding -> removed
    assert "geneR2" not in annot.all_gene_ids
    # geneR3 is pure non-coding -> kept
    assert "geneR3" in annot.all_gene_ids
    # geneR1 was mixed -> kept, but only has lnc_RNA transcript now
    assert "geneR1" in annot.all_gene_ids
    gene = annot.chrs["chr1"]["geneR1"]
    for t in gene.transcripts.values():
        assert t.coding is False
        assert t.feature == "lnc_RNA"


def test_filter_rna_classes(rich_gff3_file, tmp_path):
    output_dir = tmp_path / "filter_out"
    output_file = "lnc_only.gff3"

    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "-r", "lnc_RNA",
        "-q",
    ]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_gff = output_dir / output_file
    assert out_gff.exists()

    annot = Annotation(str(out_gff), quiet=True)
    # geneR2 (mRNA only) should have been removed entirely (no empty gene left)
    assert "geneR2" not in annot.all_gene_ids
    assert "geneR3" in annot.all_gene_ids
    assert "geneR1" in annot.all_gene_ids


def test_filter_pseudogenes(pseudogene_gff3_file, rich_gff3_file, tmp_path):
    combined = tmp_path / "combined.gff3"
    with open(combined, "w") as out:
        with open(rich_gff3_file) as f1:
            out.write(f1.read())
        with open(pseudogene_gff3_file) as f2:
            for line in f2:
                if not line.startswith("#"):
                    out.write(line)

    output_dir = tmp_path / "filter_out"

    # 1. Test --skip-pseudogenes
    args = [str(combined), "-d", str(output_dir), "-o", "no_pseudo.gff3", "--skip-pseudogenes", "-q"]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0
    annot = Annotation(str(output_dir / "no_pseudo.gff3"), quiet=True)
    assert "gene_ps1" not in annot.all_gene_ids
    assert "geneR1" in annot.all_gene_ids

    # 2. Test --pseudogenes-only
    args = [str(combined), "-d", str(output_dir), "-o", "only_pseudo.gff3", "--pseudogenes-only", "-q"]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0
    annot = Annotation(str(output_dir / "only_pseudo.gff3"), quiet=True)
    assert "gene_ps1" in annot.all_gene_ids
    assert "geneR1" not in annot.all_gene_ids
    assert "geneR2" not in annot.all_gene_ids

def test_tidy_removes_empty_genes_with_features_flag(rich_gff3_file, tmp_path):
    output_dir = tmp_path / "tidy_out"
    output_file = "tidy_lnc.gff3"

    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "-f", "lnc_RNA",
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code == 0

    out_gff = output_dir / output_file
    annot = Annotation(str(out_gff), quiet=True)
    # geneR2 had only mRNA -> must be removed, not left as an empty gene
    assert "geneR2" not in annot.all_gene_ids
    assert "geneR3" in annot.all_gene_ids
    assert "geneR1" in annot.all_gene_ids
