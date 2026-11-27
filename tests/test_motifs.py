import matplotlib
import pandas as pd
from typer.testing import CliRunner

matplotlib.use("Agg")

from aegis.annotation import Annotation
from aegis.cli.motif_search import app as motif_app
from aegis.genome import Genome
from aegis.utils.genefunctions import reverse_complement

runner = CliRunner()


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


def test_find_and_plot_output_paths_never_overwrite_and_handle_tsv_and_backlist(tmp_path):
    seq = "C" * 100 + "TATAAA" + "C" * 394
    gff_rows = []
    for i, (s, e) in enumerate([(200, 230), (260, 290), (320, 350), (380, 410)], 1):
        gff_rows.extend([
            f"chr1\tsrc\tgene\t{s}\t{e}\t.\t+\t.\tID=g{i}",
            f"chr1\tsrc\tmRNA\t{s}\t{e}\t.\t+\t.\tID=t{i};Parent=g{i}",
            f"chr1\tsrc\texon\t{s}\t{e}\t.\t+\t.\tID=e{i};Parent=t{i}",
            f"chr1\tsrc\tCDS\t{s}\t{e}\t.\t+\t0\tID=c{i};Parent=t{i}",
        ])
    _, annot = _make_motif_annot(tmp_path, seq, gff_rows, promoter_size=150)

    # 1. Explicit .tsv filename: PDF saves as .pdf and query table is preserved
    out_tsv = tmp_path / "out_tsv"
    res = annot.motifs.find_and_plot(
        ["g1"], "TATAAA", 6, "q", "TATA", filename="results.tsv", output_dir=str(out_tsv), quiet=True
    )
    assert sorted(p.name for p in out_tsv.iterdir()) == [
        "results.pdf",
        "results.tsv",
        "results_random.pdf",
        "results_random.tsv",
        "results_whole_genome.pdf",
    ]
    assert res[-1] == out_tsv / "results.tsv"
    df_tsv = pd.read_csv(out_tsv / "results.tsv", sep="\t")
    assert df_tsv["gene_id"].tolist() == ["g1"]

    # 2. Extensionless filename: distinct files, query not overwritten by random/background
    out_extless = tmp_path / "out_extless"
    annot.motifs.find_and_plot(
        ["g1"], "TATAAA", 6, "q", "TATA", filename="myrun", output_dir=str(out_extless), quiet=True
    )
    assert sorted(p.name for p in out_extless.iterdir()) == [
        "myrun.csv",
        "myrun.pdf",
        "myrun_random.csv",
        "myrun_random.pdf",
        "myrun_whole_genome.pdf",
    ]
    df_extless = pd.read_csv(out_extless / "myrun.csv", sep="\t")
    assert df_extless["gene_id"].tolist() == ["g1"]

    # 3. Custom backlist: query & random outputs are STILL written alongside background PDF
    out_bg = tmp_path / "out_bg"
    annot.motifs.find_and_plot(
        ["g1"], "TATAAA", 6, "q", "TATA", backlist=["g2", "g3"], backlistname="bg", output_dir=str(out_bg), quiet=True
    )
    bg_files = sorted(p.name for p in out_bg.iterdir())
    assert len(bg_files) == 5
    assert any(name.endswith("_TATA_q.csv") for name in bg_files)
    assert any(name.endswith("_TATA_q.pdf") for name in bg_files)
    assert any(name.endswith("_TATA_q_random.csv") for name in bg_files)
    assert any(name.endswith("_TATA_q_random.pdf") for name in bg_files)
    assert any(name.endswith("_TATA_bg_as_background.pdf") for name in bg_files)

    # 4. CLI -o with .tsv works end-to-end
    genelist = tmp_path / "genes.tsv"
    genelist.write_text("g1\tquery\n", encoding="utf-8")
    out_cli = tmp_path / "out_cli"
    cli_res = runner.invoke(
        motif_app,
        [
            str(tmp_path / "annot.gff3"),
            str(tmp_path / "genome.fa"),
            "-l", str(genelist),
            "-m", "TATAAA",
            "-ps", "150",
            "-d", str(out_cli),
            "-o", "cli_out.tsv",
            "-q",
        ],
    )
    assert cli_res.exit_code == 0, cli_res.stdout
    assert (out_cli / "cli_out.tsv").exists()
    assert (out_cli / "cli_out.pdf").exists()
    assert (out_cli / "cli_out_whole_genome.pdf").exists()


