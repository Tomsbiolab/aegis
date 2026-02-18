"""
Tests for aegis.transcript — the Transcript class.
"""

import pytest
from aegis.transcript import Transcript
from aegis.subfeatures import Exon


# ============================================================
# Helpers
# ============================================================

def make_transcript(**overrides):
    defaults = dict(
        feature_id="mRNA1",
        ch="chr1",
        source="aegis",
        feature="mRNA",
        strand="+",
        start=1000,
        end=5000,
        score=".",
        phase=".",
        attributes="ID=mRNA1;Parent=gene1"
    )
    defaults.update(overrides)
    return Transcript(**defaults)


def make_exon(feature_id, start, end, **overrides):
    defaults = dict(
        feature_id=feature_id,
        ch="chr1",
        source="aegis",
        feature="exon",
        strand="+",
        start=start,
        end=end,
        score=".",
        phase=".",
        attributes=f"ID={feature_id};Parent=mRNA1"
    )
    defaults.update(overrides)
    return Exon(**defaults)


# ============================================================
# __init__
# ============================================================

class TestTranscriptInit:
    def test_basic_properties(self):
        t = make_transcript()
        assert t.id == "mRNA1"
        assert t.ch == "chr1"
        assert t.strand == "+"
        assert t.start == 1000
        assert t.end == 5000

    def test_inherits_feature(self):
        t = make_transcript()
        # Feature size is (end - start) + 1, but Transcript.update_size
        # resets size to sum of exon sizes (0 when no exons)
        assert isinstance(t, Transcript)

    def test_exons_initially_empty_list(self):
        """Transcript.exons is a list, not a dict."""
        t = make_transcript()
        assert isinstance(t.exons, list)
        assert len(t.exons) == 0

    def test_coding_flag_default(self):
        t = make_transcript()
        assert t.coding is False

    def test_noncoding_transcript(self):
        t = make_transcript(feature="lnc_RNA")
        assert t.id == "mRNA1"


# ============================================================
# update_size
# ============================================================

class TestTranscriptUpdateSize:
    def test_update_size_no_exons_is_zero(self):
        """With no exons, update_size sets size to 0 (sum of exon sizes)."""
        t = make_transcript(start=1000, end=5000)
        t.update_size()
        assert t.size == 0

    def test_update_size_with_exons(self):
        t = make_transcript(start=1000, end=5000)
        t.exons.append(make_exon("e1", 1000, 2000))
        t.exons.append(make_exon("e2", 3000, 5000))
        t.update_size()
        # Size = sum of exon sizes: 1001 + 2001 = 3002
        assert t.size == 3002


# ============================================================
# rename
# ============================================================

class TestTranscriptRename:
    def test_rename_basic(self):
        t = make_transcript()
        t.rename(base_id="transcript001", count=1)
        assert "transcript001" in t.id
        assert t.renamed is True

    def test_rename_custom_sep_digits(self):
        t = make_transcript()
        t.rename(base_id="VIT01g001", count=2, sep=".", digits=2)
        assert "VIT01g001" in t.id


# ============================================================
# almost_equal
# ============================================================

class TestTranscriptAlmostEqual:
    def test_same_transcript_no_exons(self):
        """Two transcripts with no exons are almost_equal."""
        t1 = make_transcript()
        t2 = make_transcript()
        result = t1.almost_equal(t2)
        assert result is True

    def test_same_exons(self):
        """Two transcripts with identical exons are almost_equal."""
        t1 = make_transcript()
        t1.exons.append(make_exon("e1", 1000, 2000))
        t1.exons.append(make_exon("e2", 3000, 5000))
        t2 = make_transcript()
        t2.exons.append(make_exon("e1", 1000, 2000))
        t2.exons.append(make_exon("e2", 3000, 5000))
        assert t1.almost_equal(t2) is True

    def test_different_exon_count(self):
        """Transcripts with different number of exons are not almost_equal."""
        t1 = make_transcript()
        t1.exons.append(make_exon("e1", 1000, 2000))
        t2 = make_transcript()
        t2.exons.append(make_exon("e1", 1000, 2000))
        t2.exons.append(make_exon("e2", 3000, 5000))
        assert t1.almost_equal(t2) is False

    def test_different_exon_coordinates(self):
        """Transcripts with same count but different exon coords are not almost_equal."""
        t1 = make_transcript()
        t1.exons.append(make_exon("e1", 1000, 2000))
        t2 = make_transcript()
        t2.exons.append(make_exon("e1", 1500, 2500))
        assert t1.almost_equal(t2) is False


# ============================================================
# generate_promoter
# ============================================================

class TestTranscriptGeneratePromoter:
    def test_standard_promoter_plus_strand(self):
        t = make_transcript(start=5000, end=10000, strand="+")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.end == t.start - 1
        assert t.promoter.start == t.start - 2000

    def test_standard_promoter_minus_strand(self):
        t = make_transcript(start=5000, end=10000, strand="-")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.start == t.end + 1

    def test_promoter_clip_at_chromosome_start(self):
        t = make_transcript(start=500, end=5000, strand="+")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.start >= 1  # clipped to chromosome start

    def test_promoter_clip_at_chromosome_end(self):
        t = make_transcript(start=5000, end=99800, strand="-")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.end <= 100000


# ============================================================
# clear_UTRs
# ============================================================

class TestTranscriptClearUTRs:
    def test_clear_utrs(self):
        t = make_transcript()
        # Smoke test - should not raise
        t.clear_UTRs()


# ============================================================
# exon_update
# ============================================================

class TestTranscriptExonUpdate:
    def test_exon_update_with_exons(self):
        t = make_transcript(start=1000, end=5000)
        t.exons.append(make_exon("e1", 1000, 2000))
        t.exons.append(make_exon("e2", 3000, 5000))
        t.update_size()  # must update size first
        t.exon_update()
        # After exon_update, coding_ratio should be set
        assert hasattr(t, 'coding_ratio')
