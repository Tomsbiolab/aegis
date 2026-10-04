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


def test_tidy_cli_translation_options(tmp_path):
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

    args = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--rework-all-CDSs",
        "-gc", "1",
        "--auto-organelle-codes",
        "--mito-code", "2",
        "--plastid-code", "11",
        "--adjust-internal-shifts", "intra_exon",
        "-d", str(output_dir),
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"


def test_tidy_cli_unknown_mito_contig(tmp_path):
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

    args = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--rework-all-CDSs",
        "--mitochondria-chroms", "nonexistent_mito",
        "-d", str(output_dir),
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code != 0
    assert "Specified mitochondrial chromosome 'nonexistent_mito' was not found" in str(result.exception or result.stdout)


def test_tidy_cli_invalid_adjust_shifts(tmp_path):
    gff_file = tmp_path / "test.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    output_dir = tmp_path / "tidy_out"

    args = [
        str(gff_file),
        "--adjust-internal-shifts", "invalid_shift",
        "-d", str(output_dir),
        "-q",
    ]
    result = runner.invoke(tidy_app, args)
    assert result.exit_code != 0
    assert "Invalid adjust_internal_shifts" in result.output


def test_tidy_cli_min_codon_len(tmp_path):
    """Ensure --min-codon-len suppresses micro-ORFs shorter than the threshold."""
    gff_file = tmp_path / "micro.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    fa_file = tmp_path / "micro.fasta"
    fa_file.write_text(">chr1\nATG" + "AAA" * 8 + "TAA\n")

    out_dir = tmp_path / "out"

    # With --min-codon-len 15, the 10-codon ORF should be suppressed (no CDS)
    args = [
        str(gff_file),
        "--genome-file", str(fa_file),
        "--infer-missing-CDSs",
        "--min-codon-len", "15",
        "-d", str(out_dir),
        "-o", "filtered.gff3",
        "-q",
    ]
    res = runner.invoke(tidy_app, args)
    assert res.exit_code == 0
    annot = Annotation(str(out_dir / "filtered.gff3"), quiet=True)
    t = annot.chrs["chr1"]["g1"].transcripts["t1"]
    assert len(t.CDSs) == 0


def test_tidy_cli_recalculate_phases(tmp_path):
    """Ensure --recalculate-phases fixes intron phase mismatches while preserving 5' phase."""
    gff_file = tmp_path / "phase.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1000\t2099\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1000\t2099\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1000\t1099\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\ttest\texon\t2000\t2099\t.\t+\t.\tID=e2;Parent=t1\n"
        "chr1\ttest\tCDS\t1000\t1099\t.\t+\t1\tID=c1;Parent=t1\n"
        "chr1\ttest\tCDS\t2000\t2099\t.\t+\t2\tID=c1;Parent=t1\n"
    )
    out_dir = tmp_path / "out"

    args = [
        str(gff_file),
        "--recalculate-phases",
        "-d", str(out_dir),
        "-o", "recalc.gff3",
        "-q",
    ]
    res = runner.invoke(tidy_app, args)
    assert res.exit_code == 0
    annot = Annotation(str(out_dir / "recalc.gff3"), quiet=True)
    cds = list(annot.chrs["chr1"]["g1"].transcripts["t1"].CDSs.values())[0]
    segs = cds.CDS_segments
    assert segs[0].phase == 1  # 5' partial phase preserved!
    assert segs[1].phase == 0  # recalculated to match leftover: (100-1)%3 = 0
    assert len(annot.warnings["phase_mismatch_across_intron"]) == 0


def test_tidy_cli_reset_phases_zero(tmp_path):
    """Ensure --reset-phases-zero resets initial phase to 0 and recalculates downstream phases."""
    gff_file = tmp_path / "phase.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1000\t2099\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1000\t2099\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1000\t1099\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\ttest\texon\t2000\t2099\t.\t+\t.\tID=e2;Parent=t1\n"
        "chr1\ttest\tCDS\t1000\t1099\t.\t+\t1\tID=c1;Parent=t1\n"
        "chr1\ttest\tCDS\t2000\t2099\t.\t+\t2\tID=c1;Parent=t1\n"
    )
    out_dir = tmp_path / "out"

    args = [
        str(gff_file),
        "--reset-phases-zero",
        "-d", str(out_dir),
        "-o", "reset.gff3",
        "-q",
    ]
    res = runner.invoke(tidy_app, args)
    assert res.exit_code == 0
    annot = Annotation(str(out_dir / "reset.gff3"), quiet=True)
    cds = list(annot.chrs["chr1"]["g1"].transcripts["t1"].CDSs.values())[0]
    segs = cds.CDS_segments
    assert segs[0].phase == 0  # reset to 0!
    assert segs[1].phase == 2  # (100-0)%3 = 1 -> leftover 1 -> next phase = 2
    assert len(annot.warnings["phase_mismatch_across_intron"]) == 0
