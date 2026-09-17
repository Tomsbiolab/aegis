from typer.testing import CliRunner
from aegis.cli.summary import app

runner = CliRunner()


def test_cli_summary_single_annot_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure summary CLI runs without error on single GFF without genome."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"

    result = runner.invoke(app, [str(a1), "-d", str(tmp_path)])
    assert result.exit_code == 0, f"Single annotation run failed: {result.stdout}"


def test_cli_summary_with_genome_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure summary CLI runs without error with matching GFF and FASTA genome."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    g1 = test_data_dir / "input/fasta/minimal.fasta"

    # Positional GFF + FASTA (backwards compatibility test)
    res_pos = runner.invoke(app, [str(a1), str(g1), "-d", str(tmp_path)])
    assert res_pos.exit_code == 0, f"Positional run failed: {res_pos.stdout}"

    # With -g / --genome option
    res_opt = runner.invoke(app, [str(a1), "-g", str(g1), "-d", str(tmp_path)])
    assert res_opt.exit_code == 0, f"Option run failed: {res_opt.stdout}"


def test_cli_summary_multi_annot_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure summary CLI runs without error on multiple annotations and reference mode."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    a2 = test_data_dir / "input/annotation/extract_test.gff3"

    # Standard multi-annotation
    res_multi = runner.invoke(app, [str(a1), str(a2), "-a", "Ann1,Ann2", "-d", str(tmp_path)])
    assert res_multi.exit_code == 0, f"Multi-annotation run failed: {res_multi.stdout}"

    # Reference mode (-r)
    res_ref = runner.invoke(app, [str(a1), str(a2), "-a", "Ann1,Ann2", "-r", "-d", str(tmp_path)])
    assert res_ref.exit_code == 0, f"Reference run failed: {res_ref.stdout}"


def test_cli_summary_export(test_data_dir, tmp_path):
    """Ensure summary export option creates a non-empty output file."""
    a1 = test_data_dir / "input/annotation/minimal.gff3"
    out_file = tmp_path / "annot_summary.tsv"

    result = runner.invoke(app, [str(a1), "-o", str(out_file), "-d", str(tmp_path), "-q"])
    assert result.exit_code == 0, f"Export run failed: {result.stdout}"
    assert out_file.exists(), "Export file was not created"
    assert out_file.stat().st_size > 0, "Export file is empty"

def test_cli_summary_fatal_mismatch_halt(tmp_path):
    """
    Ensure total mismatch (zero matching contigs between annotation and genome)
    halts execution with exit code 1 and emits an error message.
    """
    gff_content = (
        "##gff-version 3\n"
        "chrA\taegis\tgene\t10\t100\t.\t+\t.\tID=geneA\n"
        "chrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=txA;Parent=geneA\n"
    )
    gff_file = tmp_path / "mismatch.gff3"
    gff_file.write_text(gff_content)

    fa_content = ">chrB\n" + "A" * 500 + "\n"
    fa_file = tmp_path / "mismatch.fasta"
    fa_file.write_text(fa_content)

    result = runner.invoke(app, [str(gff_file), "-g", str(fa_file)])
    assert result.exit_code != 0, "Expected non-zero exit code on fatal mismatch"
    output = (result.stderr or "") + (result.stdout or "")
    assert "match" in output.lower() or "error" in output.lower()
