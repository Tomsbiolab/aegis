import pytest
from typer.testing import CliRunner

from aegis.cli.subset import app

runner = CliRunner()


def test_cli_subset_gene_cap_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure subset CLI runs with --gene-cap without errors."""
    gff3_path = test_data_dir / "input/annotation/arabidopsis_araport11.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "subset_cap.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--no-chr-cap",
        "--gene-cap", "50",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert (output_dir / output_annot).exists()


def test_cli_subset_chr_cap_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure subset CLI runs with --chr-cap without errors."""
    gff3_path = test_data_dir / "input/annotation/for_merge_2.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "subset_chr.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--chr-cap", "2",
        "--no-min-genes",
        "--no-gene-cap",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert (output_dir / output_annot).exists()
