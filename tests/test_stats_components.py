from aegis.annotation import Annotation
from aegis.genome import Genome


def test_calculate_transcript_masking_requires_hard_masked_genome(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 25 + "\n", encoding="utf-8")
    hard_fa = tmp_path / "g_hard.fa"
    hard_fa.write_text(">chr1\n" + "N" * 50 + "ACGT" * 12 + "AC\n", encoding="utf-8")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1\t40\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tmRNA\t1\t40\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tsrc\texon\t1\t40\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tsrc\tCDS\t1\t36\t.\t+\t0\tID=c1;Parent=t1\n",
        encoding="utf-8",
    )

    genome = Genome("g", str(fa), quiet=True)
    # Without hard_masked_genome, calculate_transcript_masking must not crash
    annot_no_hard = Annotation(str(gff), genome=genome, quiet=True)
    annot_no_hard.stats.calculate_transcript_masking()
    assert annot_no_hard.chrs["chr1"]["g1"].quality.masked_fraction == 0

    # With hard_masked_genome, masking fraction is computed
    hard_genome = Genome("g_hard", str(hard_fa), quiet=True)
    annot_hard = Annotation(str(gff), genome=genome, hard_masked_genome=hard_genome, quiet=True)
    annot_hard.stats.calculate_transcript_masking()
    assert annot_hard.chrs["chr1"]["g1"].quality.masked_fraction == 1.0


def test_contig_stats_ignores_out_of_bounds_intervals_and_update_counts_transcriptless_genes(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 25 + "\n>chr2\n" + "ACGT" * 25 + "\n", encoding="utf-8")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # chr1: coding gene (10..50 = 41 bp), pseudogene (60..80 = 21 bp), out-of-bounds gene (150..200)
        "chr1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tmRNA\t10\t50\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tsrc\texon\t10\t50\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tsrc\tCDS\t10\t45\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tsrc\tpseudogene\t60\t80\t.\t+\t.\tID=pg1\n"
        "chr1\tsrc\tgene\t150\t200\t.\t+\t.\tID=g_oob\n"
        # chr2: ONLY an out-of-bounds gene -> coverage_pct must be 0.0, not negative
        "chr2\tsrc\tgene\t150\t200\t.\t+\t.\tID=g_oob2\n",
        encoding="utf-8",
    )

    genome = Genome("g", str(fa), quiet=True)
    annot = Annotation(str(gff), genome=genome, quiet=True)

    contig_stats = {row["contig"]: row for row in annot.stats.get_contig_stats()}
    assert contig_stats["chr1"]["coverage_pct"] == 62.0
    assert contig_stats["chr2"]["coverage_pct"] == 0.0

    annot.stats.update(quiet=True)
    assert annot.stats.data["coding_genes"] == 1
    assert annot.stats.data["noncoding_genes"] == 3

