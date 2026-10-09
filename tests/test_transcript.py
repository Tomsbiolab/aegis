"""
Tests for aegis.transcript — the Transcript class.
"""

import pytest

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.transcript import Transcript
from aegis.feature import Feature
from aegis.subfeatures import CDS, UTR

# ============================================================
# __init__
# ============================================================

class TestTranscriptInit:
    def test_basic_properties(self, make_transcript):
        t = make_transcript()
        assert t.id == "mRNA1"
        assert t.ch == "chr1"
        assert t.strand == "+"
        assert t.start == 1000
        assert t.end == 5000

    def test_inherits_feature(self, make_transcript):
        t = make_transcript()
        assert isinstance(t, Transcript)

    def test_exons_initially_empty_list(self, make_transcript):
        """Transcript.exons is a list, not a dict."""
        t = make_transcript()
        assert isinstance(t.exons, list)
        assert len(t.exons) == 0

    def test_coding_flag_default(self, make_transcript):
        t = make_transcript()
        assert t.coding is False

    def test_noncoding_transcript(self, make_transcript):
        t = make_transcript(feature="lnc_RNA")
        assert t.id == "mRNA1"


# ============================================================
# update_size
# ============================================================

class TestTranscriptUpdateSize:
    def test_update_size_no_exons_is_zero(self, make_transcript):
        t = make_transcript()
        assert t.size == 0

    def test_update_size_with_exons(self, make_transcript, make_exon):
        t = make_transcript()
        t.exons.append(make_exon("e1", 1000, 2000))
        t.exons.append(make_exon("e2", 3000, 5000))
        # Size = sum of exon sizes: 1001 + 2001 = 3002
        assert t.size == 3002


# ============================================================
# rename
# ============================================================

class TestTranscriptRename:
    def test_rename_basic(self, make_transcript):
        t = make_transcript()
        t.rename(base_id="transcript001", count=1)
        assert "transcript001" in t.id
        assert t.renamed is True

    def test_rename_custom_sep_digits(self, make_transcript):
        t = make_transcript()
        t.rename(base_id="VIT01g001", count=2, sep=".", digits=2)
        assert "VIT01g001" in t.id


# ============================================================
# almost_equal
# ============================================================

class TestTranscriptAlmostEqual:
    def test_same_transcript_no_exons(self, make_transcript):
        """Two transcripts with no exons are almost_equal."""
        t1 = make_transcript()
        t2 = make_transcript()
        result = t1.almost_equal(t2)
        assert result is True

    def test_same_exons(self, make_transcript, make_exon):
        """Two transcripts with identical exons are almost_equal."""
        t1 = make_transcript()
        t1.exons.append(make_exon("e1", 1000, 2000))
        t1.exons.append(make_exon("e2", 3000, 5000))
        t2 = make_transcript()
        t2.exons.append(make_exon("e1", 1000, 2000))
        t2.exons.append(make_exon("e2", 3000, 5000))
        assert t1.almost_equal(t2) is True

    def test_different_exon_count(self, make_transcript, make_exon):
        """Transcripts with different number of exons are not almost_equal."""
        t1 = make_transcript()
        t1.exons.append(make_exon("e1", 1000, 2000))
        t2 = make_transcript()
        t2.exons.append(make_exon("e1", 1000, 2000))
        t2.exons.append(make_exon("e2", 3000, 5000))
        assert t1.almost_equal(t2) is False

    def test_different_exon_coordinates(self, make_transcript, make_exon):
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
    def test_standard_promoter_plus_strand(self, make_transcript):
        t = make_transcript(start=5000, end=10000, strand="+")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.end == t.start - 1
        assert t.promoter.start == t.start - 2000

    def test_standard_promoter_minus_strand(self, make_transcript):
        t = make_transcript(start=5000, end=10000, strand="-")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.start == t.end + 1

    def test_promoter_clip_at_chromosome_start(self, make_transcript):
        t = make_transcript(start=500, end=5000, strand="+")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.start >= 1  # clipped to chromosome start

    def test_promoter_clip_at_chromosome_end(self, make_transcript):
        t = make_transcript(start=5000, end=99800, strand="-")
        t.generate_promoter(promoter_size=2000, ch_size=100000)
        assert t.promoter is not None
        assert t.promoter.end <= 100000


# ============================================================
# clear_UTRs
# ============================================================

class TestTranscriptClearUTRs:
    def test_clear_utrs(self, make_transcript):
        t = make_transcript()
        # Smoke test - should not raise
        t.clear_UTRs()
        assert t.temp_UTRs is None


# ============================================================
# exon_update
# ============================================================

class TestTranscriptExonUpdate:
    def test_exon_update_with_exons(self, make_transcript, make_exon):
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
    def test_rename_exons_basic(self, make_transcript, make_exon):
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
    def test_rename_utrs_basic(self, make_transcript, make_exon):
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
        self.genome = Genome("test", sample_fasta_file, quiet=True)
        self.annotation = Annotation(sample_gff3_file, genome=self.genome, quiet=True)

        # Retrieve the transcript mRNA1 from gene1
        gene = self.annotation.chrs["chr1"]["gene1"]
        self.transcript = gene.transcripts["mRNA1"]

    def test_transcript_sequence_access(self):
        assert len(self.transcript.seq) > 0


# ============================================================
# generate_best_protein
# ============================================================

class MockScaffold:
    def __init__(self, seq: str):
        self.seq = seq


class MockGenome:
    def __init__(self, seq_dict: dict[str, str]):
        self.name = "mock_genome"
        self.scaffolds = {k: MockScaffold(v) for k, v in seq_dict.items()}


class TestTranscriptGenerateBestProtein:
    """Tests for Transcript.generate_best_protein()."""

    @pytest.fixture
    def setup_mock_genome(self):
        saved = Feature._ACTIVE_GENOME
        def _activate(seq: str, chrom: str = "chr1"):
            Feature._ACTIVE_GENOME = MockGenome({chrom: seq})
        yield _activate
        Feature._ACTIVE_GENOME = saved

    def test_generate_best_protein_success(self, setup_mock_genome, make_transcript, make_exon):
        # 6bp 5' UTR (CCCGGG) + ATG AAA TAA (9bp CDS) + 6bp 3' UTR (CCCGGG) = 21bp total
        seq = "CCCGGGATGAAATAACCCGGG"
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=21, strand="+")
        e = make_exon(feature_id="e1", start=1, end=21, strand="+")
        t.exons = [e]

        t.generate_best_protein(mode="orf")

        assert t.coding is True
        assert "t1_CDS1" in t.CDSs
        cds = t.CDSs["t1_CDS1"]
        assert cds.protein is not None
        assert cds.protein.seq == "MK*"
        assert cds.protein.start == 7
        assert cds.protein.end == 15
        assert cds.start == 7
        assert cds.end == 15
        assert t.CDS_size == 9
        assert t.coding_ratio == round(9 / 21, 2)
        assert len(cds.UTRs) == 2
        assert cds.UTRs[0].prime == "5'"
        assert cds.UTRs[1].prime == "3'"

    def test_generate_best_protein_failure_strict_orf(self, setup_mock_genome, make_transcript, make_exon):
        # Noncoding sequence without start or stop codon
        seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=30, strand="+")
        e = make_exon(feature_id="e1", start=1, end=30, strand="+")
        t.exons = [e]

        t.generate_best_protein(mode="orf")

        assert t.coding is False
        assert len(t.CDSs) == 0
        assert t.CDS_size == 0
        assert t.coding_ratio == 0
        assert e.coding is False

    def test_generate_best_protein_orf_or_end_fallback(self, setup_mock_genome, make_transcript, make_exon):
        # Sequence with no canonical start/stop, but mode="orf_or_end" falls back to 3' trim
        seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"  # 30bp, multiple of 3
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=30, strand="+")
        e = make_exon(feature_id="e1", start=1, end=30, strand="+")
        t.exons = [e]

        t.generate_best_protein(mode="orf_or_end")

        assert t.coding is True
        assert "t1_CDS1" in t.CDSs
        cds = t.CDSs["t1_CDS1"]
        assert cds.protein is not None
        assert cds.protein.seq == "G" * 10
        assert t.CDS_size == 30
        assert t.coding_ratio == 1.0

    def test_generate_best_protein_clears_old_cdss(self, setup_mock_genome, make_transcript, make_exon):
        seq = "CCCGGGATGAAATAACCCGGG"
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=21, strand="+")
        e = make_exon(feature_id="e1", start=1, end=21, strand="+")
        t.exons = [e]

        # Add pre-existing old CDS
        old_cds = CDS([], "old_CDS", "chr1", "aegis", "CDS", "+", 1, 21, ".")
        t.CDSs = {"old_CDS": old_cds}

        t.generate_best_protein(mode="orf")

        assert "old_CDS" not in t.CDSs
        assert "t1_CDS1" in t.CDSs

    def test_generate_best_protein_multi_exon_splicing(self, setup_mock_genome, make_transcript, make_exon):
        # Exon 1: 1..10  (CCCGGGATGA: 5' UTR 1..6, ATG 7..9, 1st base of AAG at 10)
        # Intron: 11..20 (NNNNNNNNNN)
        # Exon 2: 21..30 (AATAACCCGG: rest of AAG at 21..22, TAA 23..25, 3' UTR 26..30)
        # Spliced: CCCGGG ATG AAA TAA CCCGG (20bp)
        seq = "CCCGGGATGA" + ("N" * 10) + "AATAACCCGG"
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=30, strand="+")
        e1 = make_exon(feature_id="e1", start=1, end=10, strand="+")
        e2 = make_exon(feature_id="e2", start=21, end=30, strand="+")
        t.exons = [e1, e2]

        t.generate_best_protein(mode="orf")

        assert t.coding is True
        cds = t.CDSs["t1_CDS1"]
        assert cds.protein is not None
        assert cds.protein.seq == "MK*"
        assert cds.protein.start == 7
        assert cds.protein.end == 25
        assert len(cds.CDS_segments) == 2
        assert cds.CDS_segments[0].start == 7
        assert cds.CDS_segments[0].end == 10
        assert cds.CDS_segments[1].start == 21
        assert cds.CDS_segments[1].end == 25

    def test_generate_best_protein_strand_resolution(self, setup_mock_genome, make_transcript, make_exon):
        # On minus strand: reverse complement of ATG AAA TAA is TTA TTT CAT
        # Let genomic be: CCC TTA TTT CAT GGG (15bp)
        # Revcomp: CCC ATG AAA TAA GGG
        seq = "CCCTTATTTCATGGG"
        setup_mock_genome(seq)

        t = make_transcript(feature_id="t1", start=1, end=15, strand=".")
        e = make_exon(feature_id="e1", start=1, end=15, strand=".")
        t.exons = [e]

        t.generate_best_protein(mode="orf", always_resolve_strand=True)

        assert t.strand == "-"
        assert e.strand == "-"
        cds = t.CDSs["t1_CDS1"]
        assert cds.protein is not None
        assert cds.protein.seq == "MK*"


# ============================================================
# Exon sorting & boundary sync
# ============================================================

class TestTranscriptExonSorting:
    def test_exons_sorted_upon_update(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=3000, end=4000)
        e2 = make_exon(feature_id="e2", start=1000, end=2000)
        t.exons = [e1, e2]

        t.update()

        assert [e.start for e in t.exons] == [1000, 3000]
        assert t.start == 1000
        assert t.end == 4000

    def test_rename_exons_orders_before_naming_plus_strand(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=5000, end=6000)
        e2 = make_exon(feature_id="e2", start=1000, end=2000)
        t.exons = [e1, e2]

        t.rename_exons(base_id="T1", digits=3)

        assert t.exons[0].start == 1000
        assert t.exons[0].id == "T1_e001"
        assert t.exons[1].start == 5000
        assert t.exons[1].id == "T1_e002"

    def test_rename_exons_orders_before_naming_minus_strand(self, make_transcript, make_exon):
        t = make_transcript(strand="-")
        e1 = make_exon(feature_id="e1", start=1000, end=2000, strand="-")
        e2 = make_exon(feature_id="e2", start=5000, end=6000, strand="-")
        # Appended in ascending order, but minus strand should assign e001 to highest coordinate (5' end)
        t.exons = [e1, e2]

        t.rename_exons(base_id="T1", digits=3)

        assert t.exons[0].start == 1000
        assert t.exons[0].id == "T1_e002"
        assert t.exons[1].start == 5000
        assert t.exons[1].id == "T1_e001"

    def test_generate_introns_with_unordered_exons(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=3000, end=4000)
        e2 = make_exon(feature_id="e2", start=1000, end=2000)
        t.exons = [e1, e2]

        t.generate_introns()

        assert len(t.introns) == 1
        assert t.introns[0].start == 2001
        assert t.introns[0].end == 2999


# ============================================================
# collapse_exons
# ============================================================

class TestCollapseExons:
    def test_collapse_contiguous_exons(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=1000, end=2000)
        e2 = make_exon(feature_id="e2", start=2001, end=3000)
        t.exons = [e1, e2]

        t.collapse_exons()

        assert t.collapsed_exons is True
        assert len(t.exons) == 1
        assert t.exons[0].start == 1000
        assert t.exons[0].end == 3000
        assert t.start == 1000
        assert t.end == 3000

    def test_collapse_overlapping_exons(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=1000, end=2000)
        e2 = make_exon(feature_id="e2", start=1500, end=2500)
        t.exons = [e1, e2]

        t.collapse_exons()

        assert t.collapsed_exons is True
        assert len(t.exons) == 1
        assert t.exons[0].start == 1000
        assert t.exons[0].end == 2500

    def test_collapse_contained_exons(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=1000, end=3000)
        e2 = make_exon(feature_id="e2", start=1200, end=1800)
        t.exons = [e1, e2]

        t.collapse_exons()

        assert t.collapsed_exons is True
        assert len(t.exons) == 1
        assert t.exons[0].start == 1000
        assert t.exons[0].end == 3000

    def test_collapse_mixed_exons_maintains_order(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        # Mixed and unordered
        e_middle = make_exon(feature_id="em", start=4000, end=5000)
        e_first1 = make_exon(feature_id="ef1", start=1000, end=2000)
        e_first2 = make_exon(feature_id="ef2", start=2001, end=2500)
        e_last = make_exon(feature_id="el", start=7000, end=8000)
        t.exons = [e_middle, e_first2, e_last, e_first1]

        t.collapse_exons()

        assert t.collapsed_exons is True
        assert len(t.exons) == 3
        assert [(e.start, e.end) for e in t.exons] == [(1000, 2500), (4000, 5000), (7000, 8000)]
        assert t.start == 1000
        assert t.end == 8000

    def test_collapse_no_overlaps_unchanged(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        e1 = make_exon(feature_id="e1", start=1000, end=2000)
        e2 = make_exon(feature_id="e2", start=3000, end=4000)
        t.exons = [e1, e2]

        t.collapse_exons()

        assert t.collapsed_exons is False
        assert len(t.exons) == 2


# ============================================================
# collapse_CDS_segments
# ============================================================

class TestCollapseCDSSegments:
    def test_collapse_contiguous_cds_plus_strand(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # 3 segments on + strand:
        # seg1: 1000..1099 (100 bp)
        # seg2: 1100..1199 (100 bp, contiguous with seg1 -> merged to 1000..1199 = 200 bp)
        # seg3: 1500..1599 (100 bp, separate)
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1099)
        seg2 = make_CDS_segment("s2", strand="+", start=1100, end=1199)
        seg3 = make_CDS_segment("s3", strand="+", start=1500, end=1599)
        cds = make_CDS(segments=[seg1, seg2, seg3], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is True
        assert len(cds.CDS_segments) == 2
        # Segments must be ordered ascending
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1199
        assert cds.CDS_segments[1].start == 1500
        assert cds.CDS_segments[1].end == 1599
        assert cds.start == 1000
        assert cds.end == 1599

        # Phase check on + strand:
        # 1st segment: phase = 0, leftover = (200 - 0) % 3 = 2
        # 2nd segment: phase = 3 - 2 = 1, leftover = (100 - 1) % 3 = 0
        assert cds.CDS_segments[0].phase == 0
        assert cds.CDS_segments[1].phase == 1

        # Frame check on + strand:
        # 1st segment: (1000 + 0) % 3 = 1 -> frame = 1
        # 2nd segment: (1500 + 1) % 3 = 1 -> frame = 1
        assert cds.CDS_segments[0].frame == 1
        assert cds.CDS_segments[1].frame == 1

    def test_collapse_contiguous_cds_minus_strand(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="-")
        # 3 segments on - strand:
        # seg1: 1000..1099 (100 bp)
        # seg2: 1100..1199 (100 bp, contiguous with seg1 -> merged to 1000..1199 = 200 bp)
        # seg3: 2000..2099 (100 bp, 5' end where translation starts!)
        seg1 = make_CDS_segment("s1", strand="-", start=1000, end=1099)
        seg2 = make_CDS_segment("s2", strand="-", start=1100, end=1199)
        seg3 = make_CDS_segment("s3", strand="-", start=2000, end=2099)
        cds = make_CDS(segments=[seg1, seg2, seg3], strand="-", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is True
        assert len(cds.CDS_segments) == 2
        # Segments must be ordered ascending in genomic coordinates
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1199
        assert cds.CDS_segments[1].start == 2000
        assert cds.CDS_segments[1].end == 2099
        assert cds.start == 1000
        assert cds.end == 2099

        # Phase check on - strand (5' to 3' is reversed: seg[1] then seg[0]):
        # 5' segment is seg[1] (2000..2099, size 100):
        # phase = 0, leftover = (100 - 0) % 3 = 1
        # 3' segment is seg[0] (1000..1199, size 200):
        # phase = 3 - 1 = 2
        assert cds.CDS_segments[1].phase == 0
        assert cds.CDS_segments[0].phase == 2

        # Frame check on - strand:
        # seg[1]: (2000 + 0) % 3 = 2 -> 7 - 2 = 5
        # seg[0]: (1000 + 2) % 3 = 0 -> 3 -> 7 - 3 = 4
        assert cds.CDS_segments[1].frame == 5
        assert cds.CDS_segments[0].frame == 4

    def test_collapse_overlapping_cds(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1200)
        seg2 = make_CDS_segment("s2", strand="+", start=1150, end=1300)
        cds = make_CDS(segments=[seg1, seg2], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is True
        assert len(cds.CDS_segments) == 1
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1300
        assert cds.CDS_segments[0].phase == 0
        assert cds.start == 1000
        assert cds.end == 1300

    def test_collapse_cds_preserves_initial_phase_plus(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # seg1 and seg2 contiguous: 1000..1099 (size 100, phase 1 -> 0 leftover) and 1100..1199 (phase 0)
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1099, phase=1)
        seg2 = make_CDS_segment("s2", strand="+", start=1100, end=1199, phase=0)
        # seg3 separated: 2000..2099
        seg3 = make_CDS_segment("s3", strand="+", start=2000, end=2099, phase=0)
        cds = make_CDS(segments=[seg1, seg2, seg3], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is True
        assert len(cds.CDS_segments) == 2
        # Merged segment (1000..1199, size 200) retains initial phase 1
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1199
        assert cds.CDS_segments[0].phase == 1
        # Leftover from merged: (200 - 1) % 3 = 199 % 3 = 1
        # seg3 phase should be 3 - 1 = 2
        assert cds.CDS_segments[1].start == 2000
        assert cds.CDS_segments[1].end == 2099
        assert cds.CDS_segments[1].phase == 2

    def test_collapse_cds_preserves_initial_phase_minus(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="-")
        # seg1: 1000..1099 (3' segment on minus strand)
        seg1 = make_CDS_segment("s1", strand="-", start=1000, end=1099, phase=2)
        # seg2 & seg3 contiguous: 2000..2099 and 2100..2199 (5' is seg3 at 2100..2199, phase 2)
        seg2 = make_CDS_segment("s2", strand="-", start=2000, end=2099, phase=1)
        seg3 = make_CDS_segment("s3", strand="-", start=2100, end=2199, phase=2)
        cds = make_CDS(segments=[seg1, seg2, seg3], strand="-", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is True
        assert len(cds.CDS_segments) == 2
        # seg[1] is 2000..2199 (5' segment on minus strand, merged from s2 and s3)
        assert cds.CDS_segments[1].start == 2000
        assert cds.CDS_segments[1].end == 2199
        assert cds.CDS_segments[1].phase == 2
        # Leftover from merged 2000..2199 (size 200, phase 2): (200 - 2) % 3 = 198 % 3 = 0
        # Downstream 3' segment seg[0] (1000..1099) gets phase 0
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1099
        assert cds.CDS_segments[0].phase == 0

    def test_collapse_cds_does_not_merge_across_internal_phase_shift(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # seg1: 1000..1099 (size 100), phase 0 -> leftover = 1
        # seg2: 1100..1199 (size 100), contiguous, but phase 0 (frameshift! 1 + 0 = 1 != 0 mod 3)
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1099, phase=0)
        seg2 = make_CDS_segment("s2", strand="+", start=1100, end=1199, phase=0)
        cds = make_CDS(segments=[seg1, seg2], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        # Should NOT collapse because of the internal phase shift between s1 and s2
        assert len(cds.CDS_segments) == 2
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1099
        assert cds.CDS_segments[0].phase == 0
        assert cds.CDS_segments[1].start == 1100
        assert cds.CDS_segments[1].end == 1199
        assert cds.CDS_segments[1].phase == 0

    def test_collapse_cds_preserves_1bp_overlap_frameshift(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # seg1: 1000..1099 (size 100), phase 0 -> leftover = 1
        # seg2: 1099..1199 (1-bp overlap at 1099), phase 1 -> (1 - 1 + 1) % 3 = 1 != 0 (incompatible -1 frameshift)
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1099, phase=0)
        seg2 = make_CDS_segment("s2", strand="+", start=1099, end=1199, phase=1)
        cds = make_CDS(segments=[seg1, seg2], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        # Should NOT collapse because of incompatible phase in 1bp overlap
        assert len(cds.CDS_segments) == 2
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1099
        assert cds.CDS_segments[1].start == 1099
        assert cds.CDS_segments[1].end == 1199

    def test_collapse_cds_merges_1bp_overlap_compatible_phase(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # seg1: 1000..1099 (size 100), phase 0 -> leftover = 1
        # seg2: 1099..1199 (1-bp overlap at 1099), phase 0 -> (1 - 1 + 0) % 3 = 0 (compatible phase)
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1099, phase=0)
        seg2 = make_CDS_segment("s2", strand="+", start=1099, end=1199, phase=0)
        cds = make_CDS(segments=[seg1, seg2], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        # Should collapse because phase is compatible
        assert len(cds.CDS_segments) == 1
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1199


# ============================================================
# Boundary & overlap edge cases
# ============================================================

class TestCoordinateEdgeCases:
    def test_generate_CDSs_detects_1bp_boundary_overlap(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # Two distinct CDS segments sharing a single boundary nucleotide (1500)
        seg1 = make_CDS_segment("cds_a", strand="+", start=1000, end=1500)
        seg2 = make_CDS_segment("cds_b", strand="+", start=1500, end=2000)
        t.temp_CDSs = [seg1, seg2]

        t.generate_CDSs(quiet=True, consider_polycistronic=True)

        # Because they overlap at 1500, it should recognize multiple overlapping CDSs
        assert len(t.CDSs) == 2 or t.polycistronic in ("yes", "maybe")

    def test_generate_UTRs_skips_internal_exons(self, make_transcript, make_exon, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        e1 = make_exon("e1", start=1000, end=1200, strand="+")
        e2 = make_exon("e2", start=1300, end=1400, strand="+")  # entirely internal to CDS
        e3 = make_exon("e3", start=1500, end=1800, strand="+")
        t.exons = [e1, e2, e3]

        # CDS starts in e1 at 1100, spans all of e2 (1300..1400), and ends in e3 at 1600
        cs1 = make_CDS_segment("cs1", strand="+", start=1100, end=1200)
        cs2 = make_CDS_segment("cs2", strand="+", start=1300, end=1400)
        cs3 = make_CDS_segment("cs3", strand="+", start=1500, end=1600)
        cds = make_CDS(segments=[cs1, cs2, cs3], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.generate_UTRs()

        utr_coords = [(u.start, u.end, u.prime) for u in cds.UTRs]
        # Only 5' UTR from e1 (1000..1099) and 3' UTR from e3 (1601..1800)
        assert len(utr_coords) == 2
        assert utr_coords[0] == (1000, 1099, "5'")
        assert utr_coords[1] == (1601, 1800, "3'")

    def test_assign_UTRs_plus_and_minus_strand(self, make_transcript, make_CDS, make_CDS_segment):
        # 1. Plus strand
        t_plus = make_transcript(feature_id="t_plus", strand="+")
        cs_plus = make_CDS_segment("cs_plus", strand="+", start=1200, end=1500)
        cds_plus = make_CDS(segments=[cs_plus], strand="+", feature_id="t_plus_CDS1")
        t_plus.CDSs = {"t_plus_CDS1": cds_plus}

        # GFF UTRs without explicit prime:
        u1 = UTR("u1", "chr1", "test", "UTR", "+", 1000, 1199, ".")
        u2 = UTR("u2", "chr1", "test", "UTR", "+", 1501, 1700, ".")
        t_plus.temp_UTRs = [u1, u2]
        t_plus.assign_UTRs()

        assert len(cds_plus.UTRs) == 2
        assert cds_plus.UTRs[0].prime == "5'"
        assert cds_plus.UTRs[1].prime == "3'"

        # 2. Minus strand
        t_minus = make_transcript(feature_id="t_minus", strand="-")
        cs_minus = make_CDS_segment("cs_minus", strand="-", start=1200, end=1500)
        cds_minus = make_CDS(segments=[cs_minus], strand="-", feature_id="t_minus_CDS1")
        t_minus.CDSs = {"t_minus_CDS1": cds_minus}

        # For minus strand, 5' UTR is at higher genomic coordinates (1501..1700)
        # 3' UTR is at lower coordinates (1000..1199)
        u_minus_3p = UTR("u_minus_3p", "chr1", "test", "UTR", "-", 1000, 1199, ".")
        u_minus_5p = UTR("u_minus_5p", "chr1", "test", "UTR", "-", 1501, 1700, ".")
        t_minus.temp_UTRs = [u_minus_3p, u_minus_5p]
        t_minus.assign_UTRs()

        assert len(cds_minus.UTRs) == 2
        # UTRs are sorted by start coordinate: [1000..1199, 1501..1700]
        # On minus strand: 1000..1199 has end <= cds.end -> prime is 3'
        # 1501..1700 has end > cds.end -> prime is 5'
        assert cds_minus.UTRs[0].prime == "3'"
        assert cds_minus.UTRs[1].prime == "5'"



