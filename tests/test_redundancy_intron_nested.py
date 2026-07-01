from aegis.annotation_components.redundancy import AnnotationRedundancy
from aegis.feature import Feature
from aegis.gene import Gene
from aegis.hits import OverlapHit
from aegis.subfeatures import CDS, Exon
from aegis.transcript import Transcript


def _score11_hit(gene_id: str) -> OverlapHit:
    return OverlapHit(
        ID=gene_id,
        origin="test",
        orientation=True,
        gene_query_percent=100,
        gene_target_percent=100,
        exons_in_both=True,
        exon_query_percent=100,
        exon_target_percent=100,
        CDSs_in_both=True,
        CDS_query_percent=100,
        CDS_target_percent=100,
        protein_query_percent=100,
        protein_target_percent=100,
        target_synteny_conserved=False,
        target_copy=False,
    )


def _build_gene(gene_id: str, cds_len: int, exon_len: int, remove: bool = False) -> Gene:
    gene = Gene(False, False, gene_id, "chr1", "test", "gene", "+", 1000, 1000 + exon_len, ".")
    gene.quality.intron_nested = True
    gene.quality.UTR_intron_nested = False
    gene.quality.remove = remove

    transcript_id = f"{gene_id}_t001"
    transcript = Transcript(transcript_id, "chr1", "test", "mRNA", "+", 1000, 1000 + exon_len, ".", [gene_id])
    transcript.exons = [
        Exon("e1", "chr1", "test", "exon", "+", 1000, 1000 + exon_len, ".", [transcript_id])
    ]
    cds_segment = Feature(
        "s1", "chr1", "test", "CDS", "+", 1000, 1000 + cds_len - 1, ".", [transcript_id]
    )
    cds = CDS(
        [cds_segment],
        f"{gene_id}_cds",
        "chr1",
        "test",
        "CDS",
        "+",
        cds_segment.start,
        cds_segment.end,
        ".",
        [transcript_id],
    )
    cds.main = True
    transcript.CDSs = {cds.id: cds}
    transcript.main = True
    gene.transcripts = {transcript.id: transcript}
    gene.update()
    return gene


class _DummyAnnotation:
    def __init__(self, genes):
        self.chrs = {"chr1": {gene.id: gene for gene in genes}}
        self.all_gene_ids = {gene.id: "chr1" for gene in genes}


def _link_identical_cds_pair(gene_a: Gene, gene_b: Gene) -> None:
    gene_a.overlaps["self"] = [_score11_hit(gene_b.id)]
    gene_b.overlaps["self"] = [_score11_hit(gene_a.id)]


def test_remove_fully_intron_nested_keeps_best_same_cds_model():
    bloated = _build_gene("bloated_g", cds_len=600, exon_len=2000, remove=False)
    compact = _build_gene("compact_g", cds_len=600, exon_len=600, remove=True)
    _link_identical_cds_pair(bloated, compact)

    redundancy = AnnotationRedundancy(_DummyAnnotation([bloated, compact]))
    redundancy.remove_fully_intron_nested_genes()

    assert bloated.quality.remove is True
    assert compact.quality.remove is False
    assert compact.quality.rescue is False


def test_remove_fully_intron_nested_keeps_only_one_qualifying_model():
    bloated = _build_gene("bloated_g", cds_len=600, exon_len=2000, remove=False)
    compact_a = _build_gene("compact_a", cds_len=600, exon_len=700, remove=True)
    compact_b = _build_gene("compact_b", cds_len=600, exon_len=600, remove=True)
    _link_identical_cds_pair(bloated, compact_a)
    _link_identical_cds_pair(bloated, compact_b)
    _link_identical_cds_pair(compact_a, compact_b)

    redundancy = AnnotationRedundancy(_DummyAnnotation([bloated, compact_a, compact_b]))
    redundancy.remove_fully_intron_nested_genes()

    assert bloated.quality.remove is True
    assert compact_a.quality.remove is True
    assert compact_b.quality.remove is False


def test_remove_fully_intron_nested_removes_cluster_when_none_qualify():
    weak_a = _build_gene("weak_a", cds_len=600, exon_len=2000, remove=False)
    weak_b = _build_gene("weak_b", cds_len=600, exon_len=1800, remove=True)
    _link_identical_cds_pair(weak_a, weak_b)

    redundancy = AnnotationRedundancy(_DummyAnnotation([weak_a, weak_b]))
    redundancy.remove_fully_intron_nested_genes()

    assert weak_a.quality.remove is True
    assert weak_b.quality.remove is True
