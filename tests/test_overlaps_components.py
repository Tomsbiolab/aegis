from aegis.annotation import Annotation


def test_overlaps_clear_respects_keep_self_and_keep_other(test_data_dir):
    gff = str(test_data_dir / "input/annotation/self_overlapping_genes.gff3")
    a1 = Annotation(gff, name="a1", quiet=True)
    a2 = Annotation(gff, name="a2", quiet=True)

    a1.overlaps.detect(quiet=True)
    a1.overlaps.detect(a2, clear=False, quiet=True)
    assert len(a1.overlaps.self_genes) > 0
    assert len(a1.overlaps.other_genes) > 0

    # keep_self=True, keep_other=True -> clears nothing
    a1.overlaps.clear(keep_self=True, keep_other=True)
    assert len(a1.overlaps.self_genes) > 0
    assert len(a1.overlaps.other_genes) > 0
    assert any(g.overlaps["self"] for genes in a1.chrs.values() for g in genes.values())
    assert any(g.overlaps["other"] for genes in a1.chrs.values() for g in genes.values())

    # keep_self=True, keep_other=False -> clears only other
    a1.overlaps.clear(keep_self=True, keep_other=False)
    assert len(a1.overlaps.self_genes) > 0
    assert len(a1.overlaps.other_genes) == 0
    assert a1.overlapped_annotations == set()
    assert any(g.overlaps["self"] for genes in a1.chrs.values() for g in genes.values())
    assert all(g.overlaps["other"] == [] for genes in a1.chrs.values() for g in genes.values())

    # keep_self=False, keep_other=True -> clears only self
    a2.overlaps.clear()
    a1.overlaps.detect(a2, clear=False, quiet=True)
    a1.overlaps.clear(keep_self=False, keep_other=True)
    assert len(a1.overlaps.self_genes) == 0
    assert len(a1.overlaps.other_genes) > 0
    assert all(g.overlaps["self"] == [] for genes in a1.chrs.values() for g in genes.values())
    assert any(g.overlaps["other"] for genes in a1.chrs.values() for g in genes.values())


def test_overlaps_and_redundancy_guard_none_percents_on_exonless_and_noncoding(tmp_path):
    gff = tmp_path / "pseudo_noncoding.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tstringtie\tpseudogene\t100\t500\t.\t+\t.\tID=pg1\n"
        "chr1\taugustus\tgene\t100\t500\t.\t+\t.\tID=nc1\n"
        "chr1\taugustus\tlnc_RNA\t100\t500\t.\t+\t.\tID=nc1.1;Parent=nc1\n"
        "chr1\taugustus\texon\t100\t500\t.\t+\t.\tID=nc1.e1;Parent=nc1.1\n",
        encoding="utf-8",
    )
    annot = Annotation(str(gff), quiet=True)
    annot.overlaps.detect(quiet=True)
    # Confirm OverlapHit actually has None percents
    pg1_hits = annot.chrs["chr1"]["pg1"].overlaps["self"]
    assert len(pg1_hits) == 1
    assert pg1_hits[0].exon_query_percent is None
    assert pg1_hits[0].CDS_query_percent is None

    # None of these should raise TypeError
    annot.overlaps.add_qualitative_info(quiet=True)
    annot.redundancy.mark_transcriptomic_supported_genes(quiet=True)
    annot.redundancy.mark_abinitio_supported_genes(quiet=True)
    annot.redundancy.find_best_gene_model_nested_overlaps(["stringtie", "augustus"])
    annot.redundancy.remove_exon_overlaps(["stringtie", "augustus"])
    annot.redundancy.remove_UTRs_from_exon_overlaps()

