import networkx as nx

from aegis.annotation import Annotation
from aegis.hits import BlastHit
from aegis.misc_features import Protein


def test_add_better_ab_initio_models_picks_strongest_candidate(tmp_path):
    gff = tmp_path / "alt.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tstringtie\tgene\t100\t1000\t.\t+\t.\tID=g_keep\n"
        "chr1\tstringtie\tmRNA\t100\t1000\t.\t+\t.\tID=t_keep;Parent=g_keep\n"
        "chr1\tstringtie\texon\t100\t500\t.\t+\t.\tID=e1;Parent=t_keep\n"
        "chr1\tstringtie\texon\t600\t1000\t.\t+\t.\tID=e2;Parent=t_keep\n"
        "chr1\tstringtie\tCDS\t600\t800\t.\t+\t0\tID=c_keep;Parent=t_keep\n"
        "chr1\taugustus\tgene\t550\t900\t.\t+\t.\tID=alt_strong\n"
        "chr1\taugustus\tmRNA\t550\t900\t.\t+\t.\tID=t_strong;Parent=alt_strong\n"
        "chr1\taugustus\texon\t550\t900\t.\t+\t.\tID=es;Parent=t_strong\n"
        "chr1\taugustus\tCDS\t550\t900\t.\t+\t0\tID=cs;Parent=t_strong\n"
        "chr1\tgenemark\tgene\t560\t880\t.\t+\t.\tID=alt_weak\n"
        "chr1\tgenemark\tmRNA\t560\t880\t.\t+\t.\tID=t_weak;Parent=alt_weak\n"
        "chr1\tgenemark\texon\t560\t880\t.\t+\t.\tID=ew;Parent=t_weak\n"
        "chr1\tgenemark\tCDS\t560\t880\t.\t+\t0\tID=cw;Parent=t_weak\n",
        encoding="utf-8",
    )

    annot = Annotation(str(gff), quiet=True)
    annot.overlaps.detect(quiet=True)

    for gid, score, is_removed in [("g_keep", 100.0, False), ("alt_strong", 500.0, True), ("alt_weak", 200.0, True)]:
        g = annot.chrs["chr1"][gid]
        g.quality.remove = is_removed
        t = next(iter(g.transcripts.values()))
        c = next(iter(t.CDSs.values()))
        p = Protein(
            prot_id=f"{c.id}.prot",
            sequence="MAAA*",
            nuc_seq="ATGGCTGCTGCTTAA",
            chrom="chr1",
            start=c.start,
            end=c.end,
            readthrough=False,
        )
        p.blast_hits.append(BlastHit("db", score, 1e-10))
        c.protein = p

    annot.redundancy.add_better_ab_initio_models_as_alternative_transcripts(["db"], quiet=True)

    alt_transcripts = [t for t in annot.chrs["chr1"]["g_keep"].transcripts.values() if t.start != 100]
    assert [(t.start, t.end) for t in alt_transcripts] == [(550, 900)]


def test_remove_alternative_iterates_chromosomes_for_graphs_and_components(tmp_path):
    gff = tmp_path / "net.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tgene\t150\t300\t.\t+\t.\tID=g2\n"
        "chr1\tsrc\tgene\t250\t400\t.\t+\t.\tID=g3\n",
        encoding="utf-8",
    )
    annot = Annotation(str(gff), quiet=True)
    annot.overlaps.detect(quiet=True)

    # 1. With nx.Graph objects in networks
    graph = nx.Graph()
    graph.add_edges_from([("g1", "g2"), ("g2", "g3")])
    annot.overlaps.networks = {"chr1": [graph]}
    annot.redundancy.remove_alternative()
    assert set(graph.nodes()) == {"g1", "g3"}

    # 2. With connected-component sets from overlaps.as_networks()
    annot.overlaps.as_networks(self_mode=True)
    annot.redundancy.remove_alternative()
    assert annot.overlaps.networks["chr1"] == [{"g1", "g3"}]
