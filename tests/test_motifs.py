import matplotlib
import pandas as pd

matplotlib.use("Agg")

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.utils.genefunctions import reverse_complement


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


def test_find_and_plot_genomic_coordinates_both_strands_and_orientations(tmp_path):
    motif = "ACGTTT"
    rc_motif = reverse_complement(motif)  # AAACGT
    seq = list("C" * 600)
    # + strand gene at 250..300 -> promoter 50..249
    seq[100:106] = list(motif)      # 1-based 101..106 -> same
    seq[150:156] = list(rc_motif)   # 1-based 151..156 -> different
    # - strand gene at 320..350 -> promoter 351..550
    seq[400:406] = list(rc_motif)   # 1-based 401..406 -> same on '-' strand
    seq[450:456] = list(motif)      # 1-based 451..456 -> different on '-' strand
    seq_str = "".join(seq)

    gff_rows = []
    for gid, s, e, st in [
        ("g_plus", 250, 300, "+"),
        ("g_minus", 320, 350, "-"),
        ("g_bg1", 1, 20, "+"),
        ("g_bg2", 560, 590, "+"),
    ]:
        gff_rows.extend([
            f"chr1\tsrc\tgene\t{s}\t{e}\t.\t{st}\t.\tID={gid}",
            f"chr1\tsrc\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID=t_{gid};Parent={gid}",
            f"chr1\tsrc\texon\t{s}\t{e}\t.\t{st}\t.\tID=e_{gid};Parent=t_{gid}",
            f"chr1\tsrc\tCDS\t{s}\t{e}\t.\t{st}\t0\tID=c_{gid};Parent=t_{gid}",
        ])

    _, annot = _make_motif_annot(tmp_path, seq_str, gff_rows, promoter_size=200)
    out = tmp_path / "out"
    annot.motifs.find_and_plot(
        query_genes=["g_plus", "g_minus"],
        motif=motif,
        motif_length=len(motif),
        glistname="q",
        tf_motif_tag="M",
        output_dir=str(out),
        quiet=True,
    )

    query_csv = next(p for p in out.iterdir() if p.name.endswith("_q.csv"))
    df = pd.read_csv(query_csv, sep="\t")
    assert len(df) == 4

    by_key = {}
    for _, row in df.iterrows():
        gs, ge = int(row["genomic_start"]), int(row["genomic_end"])
        assert gs <= ge
        genome_slice = seq_str[gs - 1 : ge]
        if row["gene_strand"] == "-":
            gene_oriented = reverse_complement(genome_slice)
        else:
            gene_oriented = genome_slice
        expected_seq = motif if row["orientation with respect to gene"] == "same" else rc_motif
        assert gene_oriented == expected_seq
        by_key[(row["gene_id"], row["orientation with respect to gene"])] = (gs, ge)

    assert by_key[("g_plus", "same")] == (101, 106)
    assert by_key[("g_plus", "different")] == (151, 156)
    assert by_key[("g_minus", "same")] == (401, 406)
    assert by_key[("g_minus", "different")] == (451, 456)

