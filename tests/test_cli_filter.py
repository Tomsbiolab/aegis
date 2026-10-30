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


def test_filter_chromosomes(rich_gff3_file, tmp_path):
    """Ensure filtering by chromosome creates an output file with matching records."""
    output_dir = tmp_path / "filter_out_chr"
    output_file = "filtered_chr.gff3"
    args = [
        str(rich_gff3_file),
        "-d", str(output_dir),
        "-o", output_file,
        "-c", "chr1",
        "-q",
    ]
    result = runner.invoke(filter_app, args)
    assert result.exit_code == 0, f"Filter by chromosome failed: {result.stdout}"
    assert (output_dir / output_file).exists()


def test_filter_keeps_models_as_read_smoke(exons_to_collapse_gff3_file, tmp_path):
    """Smoke test: filter passes models through without collapsing overlapping exons or CDS rows."""
    output_dir = tmp_path / "filter_out"
    result = runner.invoke(filter_app, [exons_to_collapse_gff3_file, "-d", str(output_dir), "-o", "out.gff3", "--coding-only", "-q"])
    assert result.exit_code == 0, f"Filter CLI failed: {result.stdout}"
    rows = [line.split("\t") for line in (output_dir / "out.gff3").read_text().splitlines() if line and not line.startswith("#")]
    exons = [(int(r[3]), int(r[4])) for r in rows if r[2] == "exon" and "Parent=mRNA_ov1" in r[8]]
    assert exons == [(1000, 2500), (2000, 3500), (3500, 4000), (5000, 6000)]
