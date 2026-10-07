import pytest
from typer.testing import CliRunner

from aegis.cli.filter import app as filter_app

runner = CliRunner()


def test_filter_smoke(rich_gff3_file, tmp_path):
    """Smoke test: ensure filter CLI produces an output file without errors."""
    output_dir = tmp_path / "filter_out"
    output_file = "filtered.gff3"

    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "--coding-only",
        "-q",
    ]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0, f"Filter CLI failed: {result.stdout}"
    assert (output_dir / output_file).exists()


def test_filter_conflicting_options(rich_gff3_file, tmp_path):
    """Ensure mutually exclusive CLI options trigger validation errors."""
    output_dir = tmp_path / "filter_out"

    # Conflicting coding flags
    res1 = runner.invoke(filter_app, [str(rich_gff3_file), "-d", str(output_dir), "--coding-only", "--non-coding-only"])
    assert res1.exit_code != 0

    # Conflicting pseudogene flags
    res2 = runner.invoke(filter_app, [str(rich_gff3_file), "-d", str(output_dir), "--skip-pseudogenes", "--pseudogenes-only"])
    assert res2.exit_code != 0

    # Conflicting transposable flags
    res3 = runner.invoke(filter_app, [str(rich_gff3_file), "-d", str(output_dir), "--skip-te", "--te-only"])
    assert res3.exit_code != 0


def test_filter_invalid_inputs(rich_gff3_file, tmp_path):
    """Ensure invalid parameter values are caught."""
    output_dir = tmp_path / "filter_out"

    # Invalid RNA class
    res_rna = runner.invoke(filter_app, [str(rich_gff3_file), "-d", str(output_dir), "-r", "invalid_class"])
    assert res_rna.exit_code != 0

    # Non-positive min CDS size
    res_cds = runner.invoke(filter_app, [str(rich_gff3_file), "-d", str(output_dir), "--min-cds-size", "0"])
    assert res_cds.exit_code != 0
