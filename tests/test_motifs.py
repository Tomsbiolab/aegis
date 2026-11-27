import matplotlib

matplotlib.use("Agg")

from aegis.annotation import Annotation
from aegis.genome import Genome


def _make_motif_annot(tmp_path, seq, gff_rows, promoter_size=200):
    fa = tmp_path / "genome.fa"
    fa.write_text(f">chr1\n{seq}\n", encoding="utf-8")
    gff = tmp_path / "annot.gff3"
    gff.write_text("##gff-version 3\n" + "\n".join(gff_rows) + "\n", encoding="utf-8")
    genome = Genome("g", str(fa), quiet=True)
    annot = Annotation(str(gff), genome=genome, quiet=True)
    annot.generate_promoters(promoter_size=promoter_size)
    return genome, annot


def test_find_and_plot_handles_short_first_promoter_and_small_promoter_size(tmp_path):
    seq = "C" * 199 + "TATAAA" + "C" * 295
    gff_rows = [
        # First gene starts at base 10 -> promoter is 1..9 (9 bp < 30)
        "chr1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1",
        "chr1\tsrc\tmRNA\t10\t50\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\tsrc\texon\t10\t50\t.\t+\t.\tID=e1;Parent=t1",
        "chr1\tsrc\tCDS\t10\t50\t.\t+\t0\tID=c1;Parent=t1",
        # Second gene starts at base 300 -> promoter is 100..299 (200 bp)
        "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g2",
        "chr1\tsrc\tmRNA\t300\t400\t.\t+\t.\tID=t2;Parent=g2",
        "chr1\tsrc\texon\t300\t400\t.\t+\t.\tID=e2;Parent=t2",
        "chr1\tsrc\tCDS\t300\t400\t.\t+\t0\tID=c2;Parent=t2",
    ]
    _, annot = _make_motif_annot(tmp_path / "case1" if False else tmp_path, seq, gff_rows, promoter_size=200)
    out1 = tmp_path / "out1"
    res = annot.motifs.find_and_plot(
        query_genes=["g2"],
        motif="TATAAA",
        motif_length=6,
        glistname="q",
        tf_motif_tag="TATA",
        output_dir=str(out1),
        quiet=True,
    )
    assert res[0] == 1  # interest_count

    # Also test when promoter_size itself is < 30 bp
    annot.generate_promoters(promoter_size=20)
    out2 = tmp_path / "out2"
    annot.motifs.find_and_plot(
        query_genes=["g2"],
        motif="TATAAA",
        motif_length=6,
        glistname="q",
        tf_motif_tag="TATA",
        output_dir=str(out2),
        quiet=True,
    )
