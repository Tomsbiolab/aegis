import pytest
from typer.testing import CliRunner

from aegis.cli.extract import app

runner = CliRunner()


@pytest.mark.parametrize(
    "feature_type, extra_args",
    [
        ("gene", []),
        ("protein", ["-gc", "1", "--auto-organelle-detection"]),
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

    # Invalid internal_stops mode
    res_stops = runner.invoke(
        app,
        [str(gff3_path), str(fasta_path), "-f", "protein", "-d", str(output_dir), "--internal-stops", "drop", "-q"],
    )
    assert res_stops.exit_code != 0


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


def test_extract_warns_about_the_organelle_genetic_code_smoke(tmp_path):
    """Smoke test: proteins of an organelle contig come with a warning naming the genetic code used."""
    gff_file = tmp_path / "m.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chrM\tt\tgene\t1\t9\t.\t+\t.\tID=g2\n"
        "chrM\tt\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\n"
        "chrM\tt\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=t2\n"
    )
    fa_file = tmp_path / "m.fa"
    fa_file.write_text(">chrM\nATGTGATAA\n")
    output_dir = tmp_path / "extract_out"
    result = runner.invoke(app, [str(gff_file), str(fa_file), "-f", "protein", "-d", str(output_dir), "-q", "-tax", "vertebrate"])
    assert result.exit_code == 0, result.output
    assert "NCBI genetic code 2 (Vertebrate Mitochondrial)" in result.output


@pytest.mark.parametrize("internal_stops, keep_stop, expected", [
    ("keep", False, {"M*K", "MK"}),
    ("x", False, {"MXK", "MK"}),
    ("x", True, {"MXK*", "MXK", "MK*"}),
    ("skip", False, {"MK"}),
])
def test_extract_cli_internal_stops(tmp_path, internal_stops, keep_stop, expected):
    """Internal stops are exported as *, as X or not at all; only an existing terminal stop is kept with --keep-stop."""
    # t1: complete with an internal stop; t2: 3' partial with an internal stop; t3: complete
    seq = "ATGTAAAAATAA" + "ATGTAAAAA" + "ATGAAATAA"
    fa = tmp_path / "test.fa"
    fa.write_text(f">chr1\n{seq}\n")
    lines = ["##gff-version 3", f"chr1\ttest\tgene\t1\t{len(seq)}\t.\t+\t.\tID=g1"]
    for n, (start, end) in enumerate(((1, 12), (13, 21), (22, 30)), 1):
        lines += [
            f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID=t{n};Parent=g1",
            f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tID=e{n};Parent=t{n}",
            f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID=c{n};Parent=t{n}",
        ]
    gff = tmp_path / "test.gff3"
    gff.write_text("\n".join(lines) + "\n")
    output_dir = tmp_path / "out"

    args = [str(gff), str(fa), "-f", "protein", "-m", "all", "-d", str(output_dir), "--internal-stops", internal_stops, "-q"]
    if keep_stop:
        args.append("--keep-stop")
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.stdout
    (fasta_out,) = output_dir.iterdir()
    seqs = {line for line in fasta_out.read_text().splitlines() if line and not line.startswith(">")}
    assert seqs == expected
