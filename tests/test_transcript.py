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

