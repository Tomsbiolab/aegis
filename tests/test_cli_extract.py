import pytest
from typer.testing import CliRunner

from aegis.cli.extract import app

runner = CliRunner()


@pytest.mark.parametrize(
    "feature_type, extra_args",
    [
        ("gene", []),
        ("protein", ["-gc", "1", "--auto-organelle-codes"]),
        ("CDS", ["--feature-id", "CDS"]),
    ],
)
def test_extract_cli_smoke(test_data_dir, tmp_path, feature_type, extra_args):
    """Smoke test: ensure extract CLI generates output for main feature types."""
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "extract_out"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-f", feature_type,
        "-d", str(output_dir),
        "-q",
    ] + extra_args

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Extract CLI failed for {feature_type}: {result.stdout}"
    assert any(output_dir.iterdir()), f"No output files created in {output_dir}"


def test_extract_cli_validation_errors(test_data_dir, tmp_path):
    """Ensure invalid parameters are rejected with non-zero exit codes."""
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "extract_out"

    # Nonexistent mito chromosome
    res_mito = runner.invoke(
        app,
        [str(gff3_path), str(fasta_path), "-f", "protein", "-d", str(output_dir), "--mitochondria-chroms", "bad_chr", "-q"],
    )
    assert res_mito.exit_code != 0

    # Invalid adjust_internal_shifts mode
    res_shifts = runner.invoke(
        app,
        [str(gff3_path), str(fasta_path), "-f", "protein", "-d", str(output_dir), "--adjust-internal-shifts", "bad_mode", "-q"],
    )
    assert res_shifts.exit_code != 0


def test_extract_cli_infer_missing_cdss_smoke(tmp_path):
    """Smoke test: ensure --infer-missing-cds runs without errors on annotation lacking CDSs."""
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nATGGCCGTTTAAAAGGGCCC\n")
    gff = tmp_path / "no_cds.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t20\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    out_dir = tmp_path / "out"

    res = runner.invoke(app, [str(gff), str(fa), "--infer-missing-cds", "-d", str(out_dir), "-q"])
    assert res.exit_code == 0
