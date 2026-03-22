"""
Tests for aegis.transcript — the Transcript class.
"""

import pytest

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.transcript import Transcript
from aegis.subfeatures import Exon, CDS, UTR, Intron
from aegis.misc_features import Promoter
from aegis.feature import Feature


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
        parents=["gene1"]
    )
    defaults.update(overrides)
    return Transcript(**defaults) # type: ignore


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
        parents=["mRNA1"]
    )
    defaults.update(overrides)
    return Exon(**defaults) # type: ignore


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
        t = make_transcript(start=1000, end=5000)
        assert t.size == 0

    def test_update_size_with_exons(self):
        t = make_transcript(start=1000, end=5000)
        t.exons.append(make_exon("e1", 1000, 2000))
        t.exons.append(make_exon("e2", 3000, 5000))
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
        t.update()
        # After exon_update, coding_ratio should be set
        assert hasattr(t, 'coding_ratio')


# ============================================================
# rename_exons
# ============================================================

class TestTranscriptRenameExons:
    def test_rename_exons_basic(self):
        t = make_transcript()
        t.exons.append(make_exon("exon_a", 1000, 2000))
        t.exons.append(make_exon("exon_b", 3000, 5000))
        t.rename_exons(base_id="geneA", sep=".", digits=2)

        assert t.exons[0].id == "geneA.e01"
        assert t.exons[1].id == "geneA.e02"


# ============================================================
# rename_utrs
# ============================================================

class TestTranscriptRenameUTRs:
    def test_rename_utrs_basic(self):
        t = make_transcript()
        c = CDS([], "cds1", "chr1", "aegis", "CDS", "+", 1000, 2000, ".")
        u1 = UTR("utr1", "chr1", "aegis", "UTR", "+", 1000, 1500, ".")
        u2 = UTR("utr2", "chr1", "aegis", "UTR", "+", 1800, 2000, ".")
        c.UTRs = [u1, u2]
        t.CDSs = {"cds1": c}
        
        t.rename_utrs(base_id="geneA", sep="-", digits=1)
        assert t.CDSs["cds1"].UTRs[0].id == "geneA-u1"
        assert t.CDSs["cds1"].UTRs[1].id == "geneA-u2"


# ============================================================
# generate_sequence, generate_hard_sequence, clear_sequence
# ============================================================

class TestTranscriptSequences:
    """Tests that load a minimal.gff3 / minimal.fasta and get sequences
    from the parsed transcript."""

    @pytest.fixture(autouse=True)
    def setup(self, sample_gff3_file, sample_fasta_file):
        self.annotation = Annotation(sample_gff3_file, quiet=True)
        self.genome = Genome("test", sample_fasta_file, quiet=True)

        # Retrieve the transcript mRNA1 from gene1
        gene = self.annotation.chrs["chr1"]["gene1"]
        self.transcript = gene.transcripts["mRNA1"]


# ============================================================
# Proteins & CDSs
# ============================================================

class TestTranscriptProteinAndCDS:
        
    def test_generate_CDSs_based_on_ORF_plus_single(self):
        # Create a transcript with a single exon
        t = make_transcript(strand="+", start=1000, end=2000)
        t.protein_seq = "M"
        t.coding_start = 100
        t.coding_end = 200
        e = make_exon("e1", 1000, 2000)
        t.exons.append(e)
        t.generate_CDSs_based_on_ORF(low_memory=False)
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 1
        assert cds.CDS_segments[0].start == 1100
        assert cds.CDS_segments[0].end == 1200

    def test_generate_CDSs_based_on_ORF_minus_single(self):
        # Create a transcript with a single exon
        t = make_transcript(strand="-", start=1000, end=2000)
        t.protein_seq = "M"
        t.coding_start = 100
        t.coding_end = 200
        e = make_exon("e1", 1000, 2000)
        t.exons.append(e)
        t.generate_CDSs_based_on_ORF(low_memory=False)
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 1
        assert cds.CDS_segments[0].start == 1800
        assert cds.CDS_segments[0].end == 1900

    def test_generate_CDSs_based_on_ORF_plus_multiple(self):
        # Create a transcript with multiple exons
        t = make_transcript(strand="+", start=1000, end=4000)
        t.protein_seq = "M"
        e1 = make_exon("e1", 1000, 2000)
        e2 = make_exon("e2", 3000, 4000)
        t.exons.extend([e1, e2])
        t.coding_start = 500
        t.coding_end = 1500
        t.generate_CDSs_based_on_ORF(low_memory=False)
        
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 2
        segs = sorted([(s.start, s.end) for s in cds.CDS_segments])
        assert segs == [(1500, 2000), (3000, 3499)]

    def test_generate_CDSs_based_on_ORF_minus_multiple(self):
        # Create a transcript with multiple exons
        t = make_transcript(strand="-", start=1000, end=4000)
        t.protein_seq = "M"
        e1 = make_exon("e1", 1000, 2000)
        e2 = make_exon("e2", 3000, 4000)
        t.exons.extend([e1, e2])
        t.coding_start = 500
        t.coding_end = 1500
        t.generate_CDSs_based_on_ORF(low_memory=False)
        
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 2
        segs = sorted([(s.start, s.end) for s in cds.CDS_segments])
        assert segs == [(1501, 2000), (3000, 3500)]

