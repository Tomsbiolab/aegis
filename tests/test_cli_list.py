import pytest
from typer.testing import CliRunner
from aegis.cli.list import app

runner = CliRunner(env={"COLUMNS": "120"})


def test_list_genes_smoke(test_data_dir, tmp_path):
    """Ensure aegis list genes runs and outputs gene table."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    out_file = tmp_path / "genes.tsv"

    result = runner.invoke(app, ["genes", str(a1), "-o", str(out_file), "-d", str(tmp_path), "-q"])
    assert result.exit_code == 0, f"list genes failed: {result.stdout}"
    assert out_file.exists()
    content = out_file.read_text()
    assert "gene_id" in content


def test_list_genes_biotype_filter(test_data_dir, tmp_path):
    """Ensure aegis list genes accepts -b / --biotypes."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    out_file = tmp_path / "genes_mrna.tsv"

    result = runner.invoke(app, ["genes", str(a1), "-b", "mRNA", "-o", str(out_file), "-d", str(tmp_path), "-q"])
    assert result.exit_code == 0, f"list genes with -b failed: {result.stdout}"
    assert out_file.exists()


def test_list_transcripts_with_gene_id(test_data_dir, tmp_path):
    """Ensure aegis list transcripts supports --gene-id."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    out_file = tmp_path / "transcripts.tsv"

    result = runner.invoke(app, ["transcripts", str(a1), "--gene-id", "-o", str(out_file), "-d", str(tmp_path), "-q"])
    assert result.exit_code == 0, f"list transcripts failed: {result.stdout}"
    assert out_file.exists()
    content = out_file.read_text()
    assert "gene_id" in content
    assert "transcript_id" in content


def test_list_dynamic_template_expansion(test_data_dir, tmp_path):
    """Ensure {annotation-name} in output filename is replaced dynamically."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"

    result = runner.invoke(app, [
        "genes", str(a1),
        "-an", "myannot",
        "-o", "custom_{annotation-name}_list.tsv",
        "-d", str(tmp_path),
        "-q",
    ])
    assert result.exit_code == 0
    expected_file = tmp_path / "custom_myannot_list.tsv"
    assert expected_file.exists()
    assert not (tmp_path / "custom_{annotation-name}_list.tsv").exists()
