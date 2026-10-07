import pytest
from typer.testing import CliRunner
from aegis.cli import app

runner = CliRunner(env={"COLUMNS": "120"})

COMMANDS = [
    "extract",
    "filter",
    "merge",
    "motifs",
    "motif-search",
    "orthology",
    "overlap",
    "prune",
    "reformat",
    "rename",
    "split",
    "subset",
    "summary",
    "summary-genome",
    "symbols",
    "tidy",
    "tidy-genome",
    "list",
]

@pytest.mark.parametrize("cmd", COMMANDS)
def test_all_cli_help(cmd):
    """Ensure every command runs --help successfully without crashing and outputs usage."""
    result = runner.invoke(app, [cmd, "--help"])
    assert result.exit_code == 0, f"Command '{cmd} --help' failed: {result.stdout}"
    assert "Usage:" in result.stdout or "Options" in result.stdout


def test_list_subcommands_help():
    """Ensure list subcommands genes and transcripts output help."""
    res_genes = runner.invoke(app, ["list", "genes", "--help"])
    assert res_genes.exit_code == 0
    assert "--annotation" in res_genes.stdout
    assert "--coding-only" in res_genes.stdout

    res_transcripts = runner.invoke(app, ["list", "transcripts", "--help"])
    assert res_transcripts.exit_code == 0
    assert "--main" in res_transcripts.stdout
    assert "--only-main" in res_transcripts.stdout


def test_motif_search_options():
    """Test motif-search help displays -ml / --motif-length and IO options."""
    res = runner.invoke(app, ["motifs", "--help"])
    assert res.exit_code == 0
    assert "--motif-length" in res.stdout
    assert "-ml" in res.stdout
    assert "-a" in res.stdout
    assert "-g" in res.stdout


def test_prune_keep_option():
    """Test prune help displays --keep / --whitelist."""
    res = runner.invoke(app, ["prune", "--help"])
    assert res.exit_code == 0
    assert "--keep" in res.stdout
    assert "--whitelist" in res.stdout


def test_merge_overlap_options():
    """Test merge help displays -go, -eo, -co options."""
    res = runner.invoke(app, ["merge", "--help"])
    assert res.exit_code == 0
    assert "--max-gene-overlap" in res.stdout
    assert "-go" in res.stdout
    assert "-eo" in res.stdout
    assert "-co" in res.stdout


def test_orthology_blast_panel():
    """Test orthology help displays BLAST options in BLAST panel."""
    res = runner.invoke(app, ["orthology", "--help"])
    assert res.exit_code == 0
    assert "BLASTp Options" in res.stdout
    assert "--skip-all-blasts" in res.stdout
    assert "--skip-RBHs" in res.stdout


def test_list_options_suite():
    """Test list genes and transcripts options."""
    res_genes = runner.invoke(app, ["list", "genes", "--help"])
    assert res_genes.exit_code == 0
    assert "--biotypes" in res_genes.stdout
    assert "--rna-classes" in res_genes.stdout

    res_transcripts = runner.invoke(app, ["list", "transcripts", "--help"])
    assert res_transcripts.exit_code == 0
    assert "--biotypes" in res_transcripts.stdout
    assert "--gene-id" in res_transcripts.stdout


def test_filter_chromosomes_help():
    """Test filter help displays -c / --chromosomes."""
    res = runner.invoke(app, ["filter", "--help"])
    assert res.exit_code == 0
    assert "--chromosomes" in res.stdout
    assert "-c" in res.stdout


def test_summary_all_help():
    """Test summary help displays -A / --all."""
    res = runner.invoke(app, ["summary", "--help"])
    assert res.exit_code == 0
    assert "--all" in res.stdout
    assert "-A" in res.stdout


def test_motif_promoter_panel():
    """Test motifs help displays -m / --motif in promoter panel."""
    res = runner.invoke(app, ["motifs", "--help"])
    assert res.exit_code == 0
    assert "Promoter & Motif Options" in res.stdout
    assert "--motif" in res.stdout
    assert "--motif-tag" in res.stdout
