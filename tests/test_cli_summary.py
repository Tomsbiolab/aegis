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

    result = runner.invoke(app, [str(gff_file), "-g", str(fa_file), "-d", str(tmp_path)])
    assert result.exit_code != 0, "Expected non-zero exit code on fatal mismatch"
    output = (result.stderr or "") + (result.stdout or "")
    assert "match" in output.lower() or "error" in output.lower()
    assert "Hint:" in output


def test_cli_summary_multi_annot_fatal_mismatch_hint(tmp_path):
    """
    Ensure fatal mismatch in multi-annotation mode provides multi-annotation specific hint.
    """
    gff1_content = "##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=geneA\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=txA;Parent=geneA\n"
    gff1 = tmp_path / "a1.gff3"
    gff1.write_text(gff1_content)

    gff2_content = "##gff-version 3\nchrB\taegis\tgene\t10\t100\t.\t+\t.\tID=geneB\nchrB\taegis\tmRNA\t10\t100\t.\t+\t.\tID=txB;Parent=geneB\n"
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text(gff2_content)

    fa_content = ">chrA\n" + "A" * 500 + "\n"
    fa_file = tmp_path / "genome.fasta"
    fa_file.write_text(fa_content)

    result = runner.invoke(app, [str(gff1), str(gff2), "-g", str(fa_file), "-d", str(tmp_path)])
    assert result.exit_code != 0
    output = (result.stderr or "") + (result.stdout or "")
    assert "applies the provided genome assembly to ALL input annotations" in output
    assert "--summary-only" in output


def test_cli_summary_disjoint_contigs_notice(tmp_path):
    """
    Ensure that when multiple annotations share zero contigs without a genome,
    an informative notice banner is emitted.
    """
    gff1_content = "##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=geneA\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=txA;Parent=geneA\n"
    gff1 = tmp_path / "a1.gff3"
    gff1.write_text(gff1_content)

    gff2_content = "##gff-version 3\nchrB\taegis\tgene\t10\t100\t.\t+\t.\tID=geneB\nchrB\taegis\tmRNA\t10\t100\t.\t+\t.\tID=txB;Parent=geneB\n"
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text(gff2_content)

    result = runner.invoke(app, [str(gff1), str(gff2), "-d", str(tmp_path)])
    assert result.exit_code == 0, f"Run failed: {result.stdout}"
    assert "NOTICE: Input annotations share no common contig names" in result.stdout
    assert "--summary-only" in result.stdout


def test_cli_summary_help_text():
    """Ensure --help outputs usage scenarios and guidance."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "USAGE SCENARIOS" in result.stdout
    assert "SAME genome assembly" in result.stdout
    assert "DIFFERENT genomes" in result.stdout
    assert "SYNONYMOUS assemblies" in result.stdout


def test_cli_summary_multi_genome_synonyms(tmp_path):
    """
    Ensure 1-to-1 multi-genome comparison matches synonymous contig names via sequence.
    """
    seq = "ATGC" * 100
    fa1 = tmp_path / "g1.fa"
    fa1.write_text(f">chrA\n{seq}\n")
    fa2 = tmp_path / "g2.fa"
    fa2.write_text(f">1\n{seq}\n")

    gff1 = tmp_path / "a1.gff3"
    gff1.write_text("##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=gA\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tA;Parent=gA\n")
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text("##gff-version 3\n1\taegis\tgene\t10\t100\t.\t+\t.\tID=gB\n1\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tB;Parent=gB\n")

    res = runner.invoke(app, [str(gff1), str(gff2), "-g", f"{fa1},{fa2}", "-d", str(tmp_path)])
    assert res.exit_code == 0, f"Run failed: {res.stdout}"
    assert "chrA (1)" in res.stdout


def test_cli_summary_multi_genome_asymmetric_contig(tmp_path):
    """
    Ensure multi-genome comparison pairs common contigs and retains unique contig (chr00).
    """
    seq = "ATGC" * 100
    fa1 = tmp_path / "g1.fa"
    fa1.write_text(f">chrA\n{seq}\n>chr00\n{seq}\n")
    fa2 = tmp_path / "g2.fa"
    fa2.write_text(f">1\n{seq}\n")

    gff1 = tmp_path / "a1.gff3"
    gff1.write_text(
        "##gff-version 3\n"
        "chrA\taegis\tgene\t10\t100\t.\t+\t.\tID=gA1\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tA1;Parent=gA1\n"
        "chr00\taegis\tgene\t10\t100\t.\t+\t.\tID=gA2\nchr00\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tA2;Parent=gA2\n"
    )
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text("##gff-version 3\n1\taegis\tgene\t10\t100\t.\t+\t.\tID=gB\n1\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tB;Parent=gB\n")

    res = runner.invoke(app, [str(gff1), str(gff2), "-g", f"{fa1},{fa2}", "-d", str(tmp_path)])
    assert res.exit_code == 0, f"Run failed: {res.stdout}"
    assert "chrA (1)" in res.stdout
    assert "chr00" in res.stdout


def test_cli_summary_multi_genome_mismatched_count(tmp_path):
    """
    Ensure an error is raised when genome count does not match annotation count.
    """
    gff1 = tmp_path / "a1.gff3"
    gff1.write_text("##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=gA\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tA;Parent=gA\n")
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text("##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=gB\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tB;Parent=gB\n")

    fa1 = tmp_path / "g1.fa"
    fa1.write_text(">chrA\nAAAA\n")
    fa2 = tmp_path / "g2.fa"
    fa2.write_text(">chrA\nAAAA\n")
    fa3 = tmp_path / "g3.fa"
    fa3.write_text(">chrA\nAAAA\n")

    res = runner.invoke(app, [str(gff1), str(gff2), "-g", f"{fa1},{fa2},{fa3}", "-d", str(tmp_path)])
    assert res.exit_code != 0
    assert "Number of genome files (3) must match the number of annotation files (2)" in ((res.stderr or "") + (res.stdout or ""))


def test_cli_summary_multi_genome_completely_different_species(tmp_path):
    """
    Ensure that passing two completely different species in multi-genome mode reports summary stats only.
    """
    fa1 = tmp_path / "g1.fa"
    fa1.write_text(">chrA\n" + "A" * 400 + "\n")
    fa2 = tmp_path / "g2.fa"
    fa2.write_text(">chrB\n" + "C" * 800 + "\n")

    gff1 = tmp_path / "a1.gff3"
    gff1.write_text("##gff-version 3\nchrA\taegis\tgene\t10\t100\t.\t+\t.\tID=gA\nchrA\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tA;Parent=gA\n")
    gff2 = tmp_path / "a2.gff3"
    gff2.write_text("##gff-version 3\nchrB\taegis\tgene\t10\t100\t.\t+\t.\tID=gB\nchrB\taegis\tmRNA\t10\t100\t.\t+\t.\tID=tB;Parent=gB\n")

    res = runner.invoke(app, [str(gff1), str(gff2), "-g", f"{fa1},{fa2}", "-d", str(tmp_path)])
    assert res.exit_code == 0, f"Run failed: {res.stdout}"
    assert "MULTI-GENOME ASSEMBLY WARNING" in res.stdout
    assert "Annotation Summary Statistics" in res.stdout

