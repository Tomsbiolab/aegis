import pytest
from typer.testing import CliRunner

from aegis.cli.tidy import app as tidy_app

runner = CliRunner()


def test_tidy_rework_cds_requires_genome(tmp_path):
    """Ensure options that depend on genome FASTA fail when --genome-file is missing."""
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )

    # 1. --rework-all-cds without genome
    res1 = runner.invoke(tidy_app, [str(gff_file), "--rework-all-cds"])
    assert res1.exit_code != 0
    assert "genome" in res1.output.lower()

    # 2. --infer-missing-cds without genome
    res2 = runner.invoke(tidy_app, [str(gff_file), "--infer-missing-cds"])
    assert res2.exit_code != 0
    assert "genome" in res2.output.lower()


def test_tidy_cli_smoke(tmp_path):
    """Smoke test: ensure tidy CLI runs with genome and --rework-all-cds without errors."""
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    fa_file = tmp_path / "test.fasta"
    fa_file.write_text(">chr1\nATGAAAGGGAAAGGGAAAGGGAAATGATAA\n")

    output_dir = tmp_path / "tidy_out"
    output_file = "reworked.gff3"

    args = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--rework-all-cds",
        "-d", str(output_dir),
        "-o", output_file,
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code == 0, f"Tidy CLI failed: {result.stdout}"
    assert (output_dir / output_file).exists()


def test_tidy_cli_validation_errors(tmp_path):
    """Ensure invalid parameters are caught by CLI validation."""
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    output_dir = tmp_path / "tidy_out"

    res_shifts = runner.invoke(
        tidy_app,
        [str(gff_file), "--adjust-internal-shifts", "invalid_shift", "-d", str(output_dir), "-q"],
    )
    assert res_shifts.exit_code != 0
