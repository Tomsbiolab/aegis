from aegis.annotation_components.redundancy import AnnotationRedundancy
from aegis.feature import Feature
from aegis.gene import Gene
from aegis.hits import OverlapHit
from aegis.subfeatures import CDS, Exon
from aegis.transcript import Transcript


def _score10_hit(gene_id: str) -> OverlapHit:
    hit = OverlapHit(
        ID=gene_id,
        origin="test",
        orientation=True,
        gene_query_percent=100,
        gene_target_percent=50,
        exons_in_both=True,
        exon_query_percent=100,
        exon_target_percent=50,
        CDSs_in_both=True,
        CDS_query_percent=100,
        CDS_target_percent=50,
        protein_query_percent=100,
        protein_target_percent=50,
        target_synteny_conserved=False,
        target_copy=False,
    )
    assert hit.score == 10
    return hit


def _build_gene(gene_id: str, cds_len: int, exon_len: int) -> Gene:
    gene = Gene(False, False, gene_id, "chr1", "test", "gene", "+", 1000, 1000 + exon_len, ".")
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


def test_compete_within_alt_transcript_groups_keeps_highest_reliable_score():
    winner = _build_gene("keeper", cds_len=600, exon_len=600)
    winner.quality.potential_alt_transcript_group = 42
    winner.quality.reliable_score = 13

    loser = _build_gene("low_score_isoform", cds_len=600, exon_len=5000)
    loser.quality.potential_alt_transcript_group = 42
    loser.quality.reliable_score = 2

    redundancy = AnnotationRedundancy(_DummyAnnotation([winner, loser]))
    redundancy.compete_within_alt_transcript_groups(source_priority=[])

    assert not winner.quality.remove
    assert loser.quality.remove
    assert not loser.quality.rescue


def test_print_changes_handles_missing_overlap_between_alt_transcript_group_members(capsys):
    winner = _build_gene("keeper", cds_len=600, exon_len=600)
    winner.quality.potential_alt_transcript_group = 42
    winner.quality.reliable_score = 13

    loser = _build_gene("low_score_isoform", cds_len=600, exon_len=5000)
    loser.quality.potential_alt_transcript_group = 42
    loser.quality.reliable_score = 2

    redundancy = AnnotationRedundancy(_DummyAnnotation([winner, loser]))
    snapshot = redundancy._take_snapshot()
    redundancy.compete_within_alt_transcript_groups(source_priority=[])
    redundancy._print_changes("compete_within_alt_transcript_groups", snapshot)

    captured = capsys.readouterr()
    assert "overlap=n/a" in captured.out
    assert redundancy._pairwise_decisions[0]["overlap"]["score"] is None


def test_find_best_gene_model_skips_same_alt_transcript_group_after_group_pick():
    keeper = _build_gene("keeper", cds_len=600, exon_len=600)
    keeper.quality.potential_alt_transcript_group = 7
    keeper.quality.reliable_score = 13

    isoform = _build_gene("isoform", cds_len=600, exon_len=5000)
    isoform.quality.potential_alt_transcript_group = 7
    isoform.quality.reliable_score = 2
    isoform.quality.remove = True

    keeper.overlaps["self"] = [_score10_hit(isoform.id)]
    isoform.overlaps["self"] = [_score10_hit(keeper.id)]

    redundancy = AnnotationRedundancy(_DummyAnnotation([keeper, isoform]))
    redundancy.find_best_gene_model(source_priority=[], just_with_reliables=True)

    assert not keeper.quality.remove
    assert isoform.quality.remove


def test_compete_within_alt_transcript_groups_removes_lower_score_despite_partial_cds_overlap():
    strong = _build_gene("STRG_winner", cds_len=900, exon_len=900)
    strong.quality.potential_alt_transcript_group = 11929
    strong.quality.reliable_score = 13

    weak = _build_gene("chr_isoform", cds_len=900, exon_len=9000)
    weak.quality.potential_alt_transcript_group = 11929
    weak.quality.reliable_score = 2

    strong.overlaps["self"] = [_score10_hit(weak.id)]
    weak.overlaps["self"] = [_score10_hit(strong.id)]

    redundancy = AnnotationRedundancy(_DummyAnnotation([strong, weak]))
    redundancy.compete_within_alt_transcript_groups(source_priority=[])

    assert not strong.quality.remove
    assert weak.quality.remove
