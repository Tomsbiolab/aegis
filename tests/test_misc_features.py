"""
Tests for aegis.misc_features — Protein and Promoter classes.
"""

import pytest

from ..misc_features import Protein, Promoter
from ..feature import Feature


# ============================================================
# Protein
# ============================================================

class TestProtein:
    def test_standard_protein(self):
        # ATG AAA TAA  ->  M K *
        nuc = "ATGAAATAA"
        p = Protein("prot1", nuc, "chr1")
        assert p.id == "prot1"
        assert p.ch == "chr1"
        assert "M" in p.seq
        assert p.start == "present"
        assert p.end_stop is True
        assert p.size > 0

    def test_partial_protein_no_start(self):
        # GGG AAA TAA  -> no M at start
        nuc = "GGGAAATAA"
        p = Protein("prot2", nuc, "chr1")
        assert p.partial is True
        assert p.start in ("absent", "late")

    def test_truncated_protein(self):
        # ATG TAA AAA TAA -> early stop at position 2
        nuc = "ATGTAAAAATAA"
        p = Protein("prot3", nuc, "chr1")
        assert p.truncated is True
        assert p.early_stop is True

    def test_protein_copy(self):
        nuc = "ATGAAATAA"
        p = Protein("prot1", nuc, "chr1")
        p2 = p.copy()
        p2.id = "changed"
        assert p.id == "prot1"

    def test_blast_hits_initially_empty(self):
        nuc = "ATGAAATAA"
        p = Protein("prot1", nuc, "chr1")
        assert p.blast_hits == []


# ============================================================
# Promoter
# ============================================================

class TestPromoter:
    def test_init(self):
        pr = Promoter(
            promoter_type="standard",
            feature_id="prom1",
            ch="chr1",
            source="aegis",
            feature="promoter",
            strand="+",
            start=2000,
            end=4000,
            score=".",
            phase=".",
            attributes="ID=prom1"
        )
        assert pr.type == "standard"
        assert pr.id == "prom1"
        assert pr.size == 2001
        assert isinstance(pr, Feature)

    def test_promoter_atg_type(self):
        pr = Promoter(
            promoter_type="ATG",
            feature_id="prom2",
            ch="chr1",
            source="aegis",
            feature="promoter",
            strand="-",
            start=10000,
            end=12000,
            score=".",
            phase=".",
            attributes="ID=prom2"
        )
        assert pr.type == "ATG"
        assert pr.strand == "-"
