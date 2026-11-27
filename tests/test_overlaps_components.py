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
