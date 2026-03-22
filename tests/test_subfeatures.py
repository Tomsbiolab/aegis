"""
Tests for aegis.subfeatures — CDS, Exon, UTR, Intron classes.
"""

import pytest

from aegis.subfeatures import Exon, UTR, Intron
from aegis.feature import Feature


# ============================================================
# CDS
# ============================================================

class TestCDS:
    def test_init(self, make_CDS):
        cds = make_CDS()
        assert cds.id == "cds1"
        assert len(cds.CDS_segments) == 2
        assert cds.start == 1200
        assert cds.end == 4500

    def test_update_size(self, make_CDS):
        cds = make_CDS()
        # Size should be sum of segment sizes
        expected = (2000 - 1200 + 1) + (4500 - 3000 + 1)
        assert cds.size == expected

    def test_update_phase(self, make_CDS, make_CDS_segment):
        segments = [
            make_CDS_segment("s1", 100, 400),
            make_CDS_segment("s2", 500, 700),
        ]
        cds = make_CDS(segments=segments)
        cds.phase = None
        cds.update_phase()
        # Phase of second segment should be computed
        assert cds.CDS_segments[1].phase is not None

    def test_equal_segments_same(self, make_CDS):
        cds1 = make_CDS()
        cds2 = make_CDS()
        assert cds1.equal_segments(cds2) is True

    def test_equal_segments_different(self, make_CDS, make_CDS_segment):
        cds1 = make_CDS()
        cds2 = make_CDS(segments=[make_CDS_segment("s1", 100, 200)])
        assert cds1.equal_segments(cds2) is False

    def test_clear_utrs(self, make_CDS):
        cds = make_CDS()
        # Set up UTR state first so clear_UTRs has something to delete
        cds.UTRs = ["utr1", "utr2"] # type: ignore
        cds.full_UTR_exons = 2
        cds.clear_UTRs()
        assert cds.UTRs == []
        assert cds.full_UTR_exons == 0


# ============================================================
# Exon
# ============================================================

class TestExon:
    def test_inherits_from_feature(self):
        e = Exon(
            feature_id="exon1",
            ch="chr1",
            source="aegis",
            feature="exon",
            strand="+",
            start=1000,
            end=2000,
            score=".",
            parents=["mRNA1"]
        )
        assert e.id == "exon1"
        assert e.size == 1001
        assert isinstance(e, Feature)


# ============================================================
# UTR
# ============================================================

class TestUTR:
    def test_default_prime(self):
        u = UTR(
            feature_id="utr1",
            ch="chr1",
            source="aegis",
            feature="three_prime_UTR",
            strand="+",
            start=4501,
            end=5000,
            score=".",
            parents=["mRNA1"]
        )
        assert u.prime == "3'"

    def test_inherits_from_feature(self):
        u = UTR(
            feature_id="utr1",
            ch="chr1",
            source="aegis",
            feature="five_prime_UTR",
            strand="+",
            start=1000,
            end=1199,
            score=".",
            parents=["mRNA1"]
        )
        assert u.size == 200
        assert isinstance(u, Feature)


# ============================================================
# Intron
# ============================================================

class TestIntron:
    def test_init_defaults(self):
        i = Intron(
            feature_id="intron1",
            ch="chr1",
            source="aegis",
            feature="intron",
            strand="+",
            start=2001,
            end=2999,
            score=".",
            parents=["mRNA1"]
        )
        assert i.id == "intron1"
        assert i.intra_coding is False

    def test_size(self):
        i = Intron(
            feature_id="intron1",
            ch="chr1",
            source="aegis",
            feature="intron",
            strand="+",
            start=2001,
            end=2999,
            score="."
        )
        assert i.size == 999

    def test_inherits_from_feature(self):
        i = Intron(
            feature_id="intron1",
            ch="chr1",
            source="aegis",
            feature="intron",
            strand="+",
            start=2001,
            end=2999,
            score="."
        )
        assert isinstance(i, Feature)
