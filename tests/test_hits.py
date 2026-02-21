"""
Tests for aegis.hits — OverlapHit and BlastHit classes.
"""

import pytest

from aegis.hits import OverlapHit, BlastHit


# ============================================================
# BlastHit
# ============================================================

class TestBlastHit:
    def test_init(self):
        bh = BlastHit(source="blast_prot", score=250.0, evalue=1e-50)
        assert bh.source == "blast_prot"
        assert bh.score == 250.0
        assert bh.evalue == 1e-50


# ============================================================
# OverlapHit
# ============================================================

class TestOverlapHit:
    def test_high_cds_overlap_same_orientation(self):
        """CDSs in both, 100% overlap, same orientation → score 11."""
        hit = OverlapHit(
            ID="g1", origin="annot1", orientation=True,
            gene_query_percent=100, gene_target_percent=100,
            exons_in_both=True, exon_query_percent=100, exon_target_percent=100,
            CDSs_in_both=True, CDS_query_percent=100, CDS_target_percent=100,
            protein_query_percent=100, protein_target_percent=100,
            target_synteny_conserved=True, target_copy=False
        )
        assert hit.score == 11

    def test_low_cds_overlap(self):
        """CDSs in both, 5% overlap → score 5."""
        hit = OverlapHit(
            ID="g2", origin="annot1", orientation=True,
            gene_query_percent=50, gene_target_percent=50,
            exons_in_both=True, exon_query_percent=20, exon_target_percent=20,
            CDSs_in_both=True, CDS_query_percent=5, CDS_target_percent=5,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 5

    def test_no_cds_no_exon_gene_only(self):
        """No exons, no CDSs, gene overlap 80% → score 7."""
        hit = OverlapHit(
            ID="g3", origin="annot1", orientation=True,
            gene_query_percent=80, gene_target_percent=80,
            exons_in_both=False, exon_query_percent=None, exon_target_percent=None,
            CDSs_in_both=False, CDS_query_percent=None, CDS_target_percent=None,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 7

    def test_exon_only_overlap(self):
        """Exons in both but no CDSs, 95% exon overlap → score 9."""
        hit = OverlapHit(
            ID="g4", origin="annot1", orientation=True,
            gene_query_percent=95, gene_target_percent=95,
            exons_in_both=True, exon_query_percent=95, exon_target_percent=95,
            CDSs_in_both=False, CDS_query_percent=None, CDS_target_percent=None,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 9

    def test_antisense_orientation(self):
        """Antisense overlap uses antiscore, not score."""
        hit = OverlapHit(
            ID="g5", origin="annot1", orientation=False,
            gene_query_percent=100, gene_target_percent=100,
            exons_in_both=False, exon_query_percent=None, exon_target_percent=None,
            CDSs_in_both=False, CDS_query_percent=None, CDS_target_percent=None,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 1  # default, not updated for antisense
        assert hit.antiscore == 9

    def test_zero_cds_overlap(self):
        """CDSs in both but 0% → score 3."""
        hit = OverlapHit(
            ID="g6", origin="annot1", orientation=True,
            gene_query_percent=50, gene_target_percent=50,
            exons_in_both=True, exon_query_percent=20, exon_target_percent=20,
            CDSs_in_both=True, CDS_query_percent=0, CDS_target_percent=0,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 3

    def test_zero_exon_overlap(self):
        """Exons in both but 0% → score 2."""
        hit = OverlapHit(
            ID="g7", origin="annot1", orientation=True,
            gene_query_percent=50, gene_target_percent=50,
            exons_in_both=True, exon_query_percent=0, exon_target_percent=0,
            CDSs_in_both=False, CDS_query_percent=None, CDS_target_percent=None,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=False, target_copy=False
        )
        assert hit.score == 2

    def test_metadata_fields(self):
        hit = OverlapHit(
            ID="g8", origin="annot1", orientation=True,
            gene_query_percent=80, gene_target_percent=90,
            exons_in_both=False, exon_query_percent=None, exon_target_percent=None,
            CDSs_in_both=False, CDS_query_percent=None, CDS_target_percent=None,
            protein_query_percent=None, protein_target_percent=None,
            target_synteny_conserved=True, target_copy=True
        )
        assert hit.id == "g8"
        assert hit.origin == "annot1"
        assert hit.extra_copy is True
        assert hit.target_synteny_conserved is True
        assert hit.min_gene_percent == 80.0
