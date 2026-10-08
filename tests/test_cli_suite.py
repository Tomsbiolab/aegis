import pytest
from typer.testing import CliRunner

from aegis.cli import app

runner = CliRunner()


def test_suite_panel_orders_all_tools():
    """Verify that in all tools, Input / Output Options is penultimate and Execution & Debugging is last."""
    commands_to_check = [
        "extract", "filter", "merge", "motifs", "orthology", "overlap",
        "prune", "reformat", "rename", "split", "subset", "summary",
        "summary-genome", "symbols", "tidy", "tidy-genome",
    ]

    for cmd in commands_to_check:
        res = runner.invoke(app, [cmd, "--help"])
        assert res.exit_code == 0, f"Command {cmd} --help failed: {res.stdout}"
        
        p_io = res.stdout.find("Input / Output Options")
        p_exec = res.stdout.find("Execution & Debugging")
        
        assert p_io != -1, f"Command {cmd} missing 'Input / Output Options' panel"
        assert p_exec != -1, f"Command {cmd} missing 'Execution & Debugging' panel"
        assert p_io < p_exec, f"In {cmd}, IO_PANEL must precede EXEC_PANEL"

        # If Reference FASTA Options is present, it must precede IO_PANEL
        p_fasta = res.stdout.find("Reference FASTA Options")
        if p_fasta != -1:
            assert p_fasta < p_io, f"In {cmd}, Reference FASTA Options must precede IO_PANEL"

        # If Genetic Codes is present, it must precede IO_PANEL
        p_gc = res.stdout.find("Genetic Codes")
        if p_gc != -1:
            assert p_gc < p_io, f"In {cmd}, Genetic Codes must precede IO_PANEL"


def test_suite_exec_panel_quiet_before_verbose():
    """Verify that in all commands, --quiet precedes --verbose in help output."""
    commands_to_check = [
        ["extract", "--help"],
        ["filter", "--help"],
        ["merge", "--help"],
        ["motifs", "--help"],
        ["orthology", "--help"],
        ["overlap", "--help"],
        ["prune", "--help"],
        ["reformat", "--help"],
        ["rename", "--help"],
        ["split", "--help"],
        ["subset", "--help"],
        ["summary", "--help"],
        ["summary-genome", "--help"],
        ["symbols", "--help"],
        ["tidy", "--help"],
        ["tidy-genome", "--help"],
        ["list", "genes", "--help"],
        ["list", "transcripts", "--help"],
    ]

    for cmd_args in commands_to_check:
        res = runner.invoke(app, cmd_args)
        assert res.exit_code == 0
        p_quiet = res.stdout.find("-q, --quiet")
        if p_quiet == -1:
            p_quiet = res.stdout.find("--quiet")
        p_verbose = res.stdout.find("-v, --verbose")
        if p_verbose == -1:
            p_verbose = res.stdout.find("--verbose ")
        assert -1 < p_quiet < p_verbose, f"In {' '.join(cmd_args)}, --quiet must precede --verbose"


def test_extract_mode_default_and_help():
    """Verify extract mode default is ['all'] and doesn't contain 'main'."""
    res = runner.invoke(app, ["extract", "--help"])
    assert res.exit_code == 0
    assert "[default: all]" in res.stdout or "all" in res.stdout
    assert "['all', 'main']" not in res.stdout


def test_orthology_panel_title():
    """Verify orthology does not duplicate 'Output' panel names."""
    res = runner.invoke(app, ["orthology", "--help"])
    assert res.exit_code == 0
    assert "Result Filtering & Table Formatting" in res.stdout
    assert "Output & Filtering Options" not in res.stdout


def test_split_write_empty_other_panel():
    """Verify --write-empty-other is in Split Criteria, not Execution & Debugging."""
    res = runner.invoke(app, ["split", "--help"])
    assert res.exit_code == 0
    p_split = res.stdout.find("Split Criteria")
    p_exec = res.stdout.find("Execution & Debugging")
    p_opt = res.stdout.find("--write-empty-other")
    assert -1 < p_split < p_opt < p_exec


def test_separator_aliases():
    """Verify --sep and --separator are accepted across symbols, rename, and list."""
    res_sym = runner.invoke(app, ["symbols", "--help"])
    assert res_sym.exit_code == 0
    assert "--sep" in res_sym.stdout
    assert "--separator" in res_sym.stdout

    res_ren = runner.invoke(app, ["rename", "--help"])
    assert res_ren.exit_code == 0
    assert "--sep" in res_ren.stdout
    assert "--separator" in res_ren.stdout

    res_lg = runner.invoke(app, ["list", "genes", "--help"])
    assert res_lg.exit_code == 0
    assert "--sep" in res_lg.stdout
    assert "--separator" in res_lg.stdout
    assert "with or without extension" in res_lg.stdout


def test_motifs_positional_signature():
    """Verify motifs only takes 2 positional files, and rejects missing -m/--motif gracefully."""
    res = runner.invoke(app, ["motifs", "tests/test_data/input/annotation/minimal.gff3", "tests/test_data/input/fasta/minimal.fasta"])
    assert res.exit_code != 0
    assert "Missing required gene list file" in res.output or "Missing required motif" in res.output
