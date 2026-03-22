"""
Tests for aegis.subfeatures — CDS, Exon, UTR, Intron classes.
"""

import pytest

from aegis.subfeatures import CDS, Exon, UTR, Intron
from aegis.feature import Feature


# ============================================================
# Helpers
# ============================================================

def make_feature_segment(feature_id="seg1", start=100, end=300):
    """Create a Feature segment to use as a CDS segment."""
    return Feature(
        feature_id=feature_id,
        ch="chr1",
        source="aegis",
        feature="CDS",
        strand="+",
        start=start,
        end=end,
        score=".",
        parents=["mRNA1"]
    )


def make_cds(segments=None, feature_id="cds1"):
    if segments is None:
        segments = [
            make_feature_segment("seg1", 1200, 2000),
            make_feature_segment("seg2", 3000, 4500),
        ]
    first = segments[0]
    last = segments[-1]
    return CDS(
        CDS_segments=segments,
        feature_id=feature_id,
        ch="chr1",
        source="aegis",
        feature="CDS",
        strand="+",
        start=first.start,
        end=last.end,
        score=".",
        parents=["mRNA1"]
    )


# ============================================================
# CDS
# ============================================================

class TestCDS:
    def test_init(self):
        cds = make_cds()
        assert cds.id == "cds1"
        assert len(cds.CDS_segments) == 2
        assert cds.start == 1200
        assert cds.end == 4500

    def test_update_size(self):
        cds = make_cds()
        # Size should be sum of segment sizes
        expected = (2000 - 1200 + 1) + (4500 - 3000 + 1)
        assert cds.size == expected

    def test_update_phase(self):
        segments = [
            make_feature_segment("s1", 100, 400),
            make_feature_segment("s2", 500, 700),
        ]
        cds = make_cds(segments=segments)
        cds.phase = None
        cds.update_phase()
        # Phase of second segment should be computed
        assert cds.CDS_segments[1].phase is not None

    def test_equal_segments_same(self):
        cds1 = make_cds()
        cds2 = make_cds()
        assert cds1.equal_segments(cds2) is True

    def test_equal_segments_different(self):
        cds1 = make_cds()
        cds2 = make_cds(segments=[make_feature_segment("s1", 100, 200)])
        assert cds1.equal_segments(cds2) is False

    def test_clear_utrs(self):
        cds = make_cds()
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
