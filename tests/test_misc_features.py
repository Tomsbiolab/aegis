"""
Tests for aegis.misc_features — Protein and Promoter classes.
"""

import pytest

from aegis.misc_features import Protein, Promoter
from aegis.feature import Feature


# ============================================================
# Protein
# ============================================================

# class TestProtein:
#     def test_standard_protein(self):
#         # ATG AAA TAA  ->  M K *
#         nuc = "ATGAAATAA"
#         p = Protein("prot1", nuc, "chr1")
#         assert p.id == "prot1"
#         assert p.ch == "chr1"
#         assert "M" in p.seq
#         assert p.start == "present"
#         assert p.end_stop is True
#         assert p.size > 0

#     def test_partial_protein_no_start(self):
#         # GGG AAA TAA  -> no M at start
#         nuc = "GGGAAATAA"
#         p = Protein("prot2", nuc, "chr1")
#         assert p.partial is True
#         assert p.start in ("absent", "late")

#     def test_truncated_protein(self):
#         # ATG TAA AAA TAA -> early stop at position 2
#         nuc = "ATGTAAAAATAA"
#         p = Protein("prot3", nuc, "chr1")
#         assert p.truncated is True
#         assert p.early_stop is True

#     def test_protein_copy(self):
#         nuc = "ATGAAATAA"
#         p = Protein("prot1", nuc, "chr1")
#         p2 = p.copy()
#         p2.id = "changed"
#         assert p.id == "prot1"

#     def test_blast_hits_initially_empty(self):
#         nuc = "ATGAAATAA"
#         p = Protein("prot1", nuc, "chr1")
#         assert p.blast_hits == []


# ============================================================
# Promoter
# ============================================================

class TestProteinTags:
    @staticmethod
    def _protein(seq, nuc_seq, **kwargs):
        return Protein("p1", seq, "chr1", 1, len(nuc_seq), "end", nuc_seq=nuc_seq, **kwargs)

    def test_complete_protein_only_has_mode_and_table(self):
        assert self._protein("MK*", "ATGAAATAA").tags == {"readthrough": "end", "table": "1"}

    @pytest.mark.parametrize("seq, nuc_seq, trimmed_5p, expected", [
        ("MK", "ATGAAA", 0, "3'"),
        ("MK*", "ATGAAATAA", 1, "5'"),  # 5' bases skipped by phase
        ("GK", "GGGAAA", 0, "5',3'"),
    ])
    def test_partial(self, seq, nuc_seq, trimmed_5p, expected):
        assert self._protein(seq, nuc_seq, trimmed_5p=trimmed_5p).tags["partial"] == expected

    def test_alt_start_shows_codon(self):
        # CTG is an NCBI initiator of table 1 but not a curated start codon
        assert self._protein("LK*", "CTGAAATAA").tags["alt_start"] == "CTG"

    def test_internal_stops_and_ambiguous_residues(self):
        tags = self._protein("MK*KX*", "ATGAAATAAAAANNNTAA").tags
        assert (tags["internal_stops"], tags["ambiguous_residues"]) == ("1", "1")

    def test_frameshifts_are_signed(self):
        assert self._protein("MK*", "ATGAAATAA", frameshifts=(-1, 2)).tags["frameshifts"] == "-1,+2"

    def test_custom_table(self):
        assert self._protein("MK*", "ATGAAATAA", table={"ATG": "M"}).tags["table"] == "custom"

    def test_key_order(self):
        p = self._protein("GK*KX", "GGGAAATAAAAANNN", frameshifts=(-1,), exception="ribosomal slippage")
        assert list(p.tags) == ["partial", "internal_stops", "ambiguous_residues", "frameshifts", "exception", "readthrough", "table"]


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
            score="."
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
            score="."
        )
        assert pr.type == "ATG"
        assert pr.strand == "-"
