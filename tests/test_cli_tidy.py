import pytest
from typer.testing import CliRunner

from aegis.cli.tidy import app as tidy_app
from aegis.annotation import Annotation

runner = CliRunner()


def test_tidy_rework_cds_requires_genome(tmp_path):
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )

    # 1. --rework-all-CDSs without --genome-file should fail
    args = [str(gff_file), "--rework-all-CDSs"]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code != 0
    assert "A genome FASTA file must be provided" in result.output

    # 2. --infer-missing-CDSs without --genome-file should also fail
    args = [str(gff_file), "--infer-missing-CDSs"]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code != 0
    assert "A genome FASTA file must be provided" in result.output


def test_tidy_rework_all_cds_with_genome(tmp_path):
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
        "--rework-all-CDSs",
        "-d", str(output_dir),
        "-o", output_file,
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_gff = output_dir / output_file
    annot = Annotation(str(out_gff), quiet=True)
    t = annot.chrs["chr1"]["g1"].transcripts["t1"]
    assert t.coding is True
    assert len(t.CDSs) >= 1


def test_tidy_rework_cds_fallback_to_trim(tmp_path):
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )

    # 30 bp with no start/stop codon
    fa_file = tmp_path / "test.fasta"
    fa_file.write_text(">chr1\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    output_dir = tmp_path / "tidy_out"

    # Without --fallback-to-trim -> remains non-coding
    args_no_fallback = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--rework-all-CDSs",
        "-d", str(output_dir),
        "-o", "no_fallback.gff3",
        "-q",
    ]
    result = runner.invoke(tidy_app, args_no_fallback)
    assert result.exit_code == 0
    annot_no = Annotation(str(output_dir / "no_fallback.gff3"), quiet=True)
    assert annot_no.chrs["chr1"]["g1"].transcripts["t1"].coding is False

    # With --fallback-to-trim -> generates CDS
    args_fallback = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--rework-all-CDSs",
        "--fallback-to-trim",
        "-d", str(output_dir),
        "-o", "fallback.gff3",
        "-q",
    ]
    result = runner.invoke(tidy_app, args_fallback)
    assert result.exit_code == 0
    annot_fb = Annotation(str(output_dir / "fallback.gff3"), quiet=True)
    assert annot_fb.chrs["chr1"]["g1"].transcripts["t1"].coding is True
