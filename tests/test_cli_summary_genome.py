from typer.testing import CliRunner
from aegis.cli.summary_genome import app

runner = CliRunner()


def test_cli_summary_genome_smoke(test_data_dir):
    """
    Smoke test: ensure summary-genome CLI runs without error on single and multiple genomes.
    No brittle formatting checks so the CLI presentation can evolve freely.
    """
    g1 = test_data_dir / "input/fasta/minimal.fasta"
    g2 = test_data_dir / "input/fasta/extract_test.fasta"

    # Single genome
    result_single = runner.invoke(app, [str(g1)])
    assert result_single.exit_code == 0, f"Single genome run failed: {result_single.stdout}"

    # Multiple genomes comparison
    result_multi = runner.invoke(app, [str(g1), str(g2)])
    assert result_multi.exit_code == 0, f"Multi-genome run failed: {result_multi.stdout}"


def test_cli_summary_genome_export(test_data_dir, tmp_path):
    """Ensure summary-genome export option successfully creates a non-empty output file."""
    g1 = test_data_dir / "input/fasta/minimal.fasta"
    out_file = tmp_path / "summary.tsv"

    result = runner.invoke(app, [str(g1), "-o", str(out_file), "-q"])
    assert result.exit_code == 0, f"Export run failed: {result.stdout}"
    assert out_file.exists(), "Export file was not created"
    assert out_file.stat().st_size > 0, "Export file is empty"
