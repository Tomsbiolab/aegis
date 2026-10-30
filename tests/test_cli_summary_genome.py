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


def test_cli_summary_genome_synonym_pairing(tmp_path):
    """
    Ensure renamed chromosomes (e.g. chr1 vs PN1) pair correctly
    in standard mode (reporting synonym) and reference mode (reporting '= ref').
    """
    g1 = tmp_path / "g1.fa"
    g2 = tmp_path / "g2.fa"

    g1.write_text(">chr1\nACGTACGTACGTACGT\n>chr2\nAACCGGTTAACCGGTT\n")
    # PN1 has identical sequence to chr1
    g2.write_text(">PN1\nACGTACGTACGTACGT\n>chr2\nAACCGGTTAACCGGTT\n")

    # Standard mode: synonym PN1 appears in output
    res_std = runner.invoke(app, [str(g1), str(g2)])
    assert res_std.exit_code == 0, f"Standard run failed: {res_std.stdout}"
    assert "PN1" in res_std.stdout, "Expected synonym PN1 to appear in output"

    # Reference mode (-r): '= ref' appears and PN1 is noted
    res_ref = runner.invoke(app, [str(g1), str(g2), "-r"])
    assert res_ref.exit_code == 0, f"Reference run failed: {res_ref.stdout}"
    assert "= ref" in res_ref.stdout, "Expected '= ref' in reference mode"
    assert "PN1" in res_ref.stdout, "Expected synonym PN1 in reference mode"

    # Diff-only mode with reference: runs cleanly
    res_diff = runner.invoke(app, [str(g1), str(g2), "-r", "--diff-only"])
    assert res_diff.exit_code == 0, f"Diff-only run failed: {res_diff.stdout}"


def test_cli_summary_genome_ref_and_no_seq_options(tmp_path):
    """Ensure --ref-genome, --no-seq, and sort options run without error."""
    g1 = tmp_path / "g1.fa"
    g2 = tmp_path / "g2.fa"

    g1.write_text(">chr1\nACGTACGT\n")
    g2.write_text(">chr1\nACGTACGT\n")

    # Test --ref-genome
    res_ref_name = runner.invoke(app, [str(g1), str(g2), "--ref-genome", "g1"])
    assert res_ref_name.exit_code == 0, f"Ref-genome run failed: {res_ref_name.stdout}"

    # Test --no-seq
    res_no_seq = runner.invoke(app, [str(g1), str(g2), "-r", "--no-seq"])
    assert res_no_seq.exit_code == 0, f"No-seq run failed: {res_no_seq.stdout}"


def test_cli_summary_genome_soft_masked(tmp_path):
    """Ensure Soft-masked (%) row appears when soft masking is present,
    and soft-masking does not impair chromosome pairing or flag sequence diffs."""
    g1 = tmp_path / "g1.fa"
    g2 = tmp_path / "g2.fa"
    g1.write_text(">chr1\nACGTACGTACGT\n")
    g2.write_text(">PN1\nacgtacgtacgt\n")  # soft-masked + renamed

    # Standard mode: pairs chr1 and PN1 via hash matching
    res_std = runner.invoke(app, [str(g1), str(g2)])
    assert res_std.exit_code == 0, f"Run failed: {res_std.stdout}"
    assert "PN1" in res_std.stdout
    assert "Soft-masked (%)" in res_std.stdout

    # Reference mode: recognized as = ref (PN1) with zero sequence diff
    res_ref = runner.invoke(app, [str(g1), str(g2), "-r"])
    assert res_ref.exit_code == 0, f"Ref run failed: {res_ref.stdout}"
    assert "= ref (PN1)" in res_ref.stdout
    assert "*" not in res_ref.stdout  # no sequence difference marker


