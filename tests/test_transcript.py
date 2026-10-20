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

    def test_overlapping_exons_count_once(self, make_transcript, make_exon):
        t = make_transcript()
        t.exons.append(make_exon("e1", 1000, 2000))
        t.exons.append(make_exon("e2", 1501, 2500))
        assert t.size == 1501
        assert t.genomic_span == t.end - t.start + 1


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
        # Transcript coordinates are left to Annotation coordinate polishing
        assert t.start == 1000
        assert t.end == 5000

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

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_generate_introns_nested_exon(self, make_transcript, make_exon, strand):
        t = make_transcript(strand=strand)
        t.exons = [
            make_exon("e1", start=1000, end=5000, strand=strand),
            make_exon("e2", start=2000, end=3000, strand=strand),  # nested in e1
            make_exon("e3", start=6000, end=7000, strand=strand),
        ]

        t.generate_introns()

        # Only the gap after the furthest exon end (5000) is an intron
        assert [(i.start, i.end) for i in t.introns] == [(5001, 5999)]

    def test_generate_introns_overlapping_and_adjacent_exons(self, make_transcript, make_exon):
        t = make_transcript(strand="+")
        t.exons = [
            make_exon("e1", start=1000, end=2000),
            make_exon("e2", start=1500, end=2500),  # overlaps e1
            make_exon("e3", start=2501, end=3000),  # adjacent to e2
            make_exon("e4", start=4000, end=5000),
        ]

        t.generate_introns()

        assert [(i.start, i.end) for i in t.introns] == [(3001, 3999)]
        assert all(i.start <= i.end for i in t.introns)

    def test_generate_introns_ids_are_consecutive(self, make_transcript, make_exon):
        t = make_transcript(feature_id="T1", strand="+")
        t.exons = [
            make_exon("e1", start=1000, end=2000),
            make_exon("e2", start=1500, end=1800),  # nested, produces no intron
            make_exon("e3", start=3000, end=4000),
            make_exon("e4", start=5000, end=6000),
        ]

        t.generate_introns()

        assert [i.id for i in t.introns] == ["T1_intron_1", "T1_intron_2"]
        assert [(i.start, i.end) for i in t.introns] == [(2001, 2999), (4001, 4999)]

    def test_generate_introns_minus_strand_numbered_5_to_3(self, make_transcript, make_exon):
        t = make_transcript(feature_id="T1", strand="-")
        t.exons = [
            make_exon("e1", start=1000, end=2000, strand="-"),
            make_exon("e2", start=3000, end=4000, strand="-"),
            make_exon("e3", start=5000, end=6000, strand="-"),
        ]

        t.generate_introns()

        # Kept in coordinate order, numbered from the 5' end (rightmost on the minus strand)
        assert [(i.start, i.end) for i in t.introns] == [(2001, 2999), (4001, 4999)]
        assert [i.id for i in t.introns] == ["T1_intron_2", "T1_intron_1"]
        assert all(i.original_id == i.id for i in t.introns)

    def test_generate_introns_numbering_matches_exons(self, make_transcript, make_exon):
        for strand in ("+", "-"):
            t = make_transcript(feature_id="T1", strand=strand)
            t.exons = [
                make_exon("e1", start=1000, end=2000, strand=strand),
                make_exon("e2", start=3000, end=4000, strand=strand),
                make_exon("e3", start=5000, end=6000, strand=strand),
            ]
            t.rename_exons("T1")
            t.generate_introns()
            exon_number = {e.start: int(e.id.split("_e")[-1]) for e in t.exons}
            for i in t.introns:
                upstream_exon_start = max(s for s in exon_number if s < i.start) if strand == "+" else min(s for s in exon_number if s > i.end)
                # Intron n lies between exon n and exon n+1 in transcription order
                assert int(i.id.split("_intron_")[-1]) == exon_number[upstream_exon_start], strand


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
        # Transcript coordinates are left to Annotation coordinate polishing
        assert t.start == 1000
        assert t.end == 5000

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
        # Transcript coordinates are left to Annotation coordinate polishing
        assert t.start == 1000
        assert t.end == 5000

    @pytest.mark.parametrize("strand, expected_id", [("+", "e1"), ("-", "e2")])
    def test_merged_exon_keeps_5prime_id_and_attributes(self, make_transcript, make_exon, strand, expected_id):
        t = make_transcript(strand=strand)
        e1 = make_exon(feature_id="e1", start=1000, end=2000, strand=strand, attributes={"Note": "first"})
        e2 = make_exon(feature_id="e2", start=1500, end=2500, strand=strand, attributes={"Note": "second"})
        e3 = make_exon(feature_id="e3", start=4000, end=5000, strand=strand)
        t.exons = [e1, e2, e3]

        t.collapse_exons()

        assert [(e.id, e.start, e.end) for e in t.exons] == [(expected_id, 1000, 2500), ("e3", 4000, 5000)]
        assert t.exons[0].misc_attributes == [f"Note={'first' if strand == '+' else 'second'}"]
        assert t.exons[1] is e3

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

    def test_collapse_cds_keeps_segment_attributes(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        attributes = {"transl_except": "(pos:1000..1002%2Caa:Sec)", "Note": "x"}
        seg1 = make_CDS_segment("s1", strand="+", start=1000, end=1098, phase=0, attributes=attributes)
        seg2 = make_CDS_segment("s2", strand="+", start=1099, end=1200, phase=0, attributes=attributes)
        cds = make_CDS(segments=[seg1, seg2], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert len(cds.CDS_segments) == 1
        assert cds.CDS_segments[0].misc_attributes == ["transl_except=(pos:1000..1002%2Caa:Sec)", "Note=x"]
        assert cds.transl_except == ((1000, 1002, "U"),)

    def test_collapse_cds_keeps_annotated_phases_plus(self, make_transcript, make_CDS, make_CDS_segment):
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
        # seg3 keeps its annotated phase 0, although the merged leftover (1) implies 2: collapsing does not
        # recalculate phases (see recalculate_phases), the mismatch across the intron is kept as annotated
        assert cds.CDS_segments[1].start == 2000
        assert cds.CDS_segments[1].end == 2099
        assert cds.CDS_segments[1].phase == 0

    def test_collapse_cds_keeps_annotated_phases_minus(self, make_transcript, make_CDS, make_CDS_segment):
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
        # The downstream 3' segment keeps its annotated phase 2 (the merged leftover implies 0)
        assert cds.CDS_segments[0].start == 1000
        assert cds.CDS_segments[0].end == 1099
        assert cds.CDS_segments[0].phase == 2

    def test_collapse_cds_keeps_phases_of_other_segments(self, make_transcript, make_CDS, make_CDS_segment):
        # A|B: contiguous phase mismatch (leftover 1 + phase 0), kept apart; C|D: one coding stretch split in two rows
        t = make_transcript(feature_id="t1", strand="+")
        segments = [make_CDS_segment(f"s{i}", strand="+", start=a, end=b, phase=p) for i, (a, b, p) in enumerate([(1, 10, 0), (11, 20, 0), (31, 40, 2), (41, 50, 1)])]
        cds = make_CDS(segments=segments, strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert [(s.start, s.end, s.phase) for s in cds.CDS_segments] == [(1, 10, 0), (11, 20, 0), (31, 50, 2)]

    # NCBI programmed frameshifts: the phases continue the frame by reading the overlapping bases twice,
    # so (upstream leftover + downstream phase) % 3 == 0 and merging would delete the shared bases
    @pytest.mark.parametrize("strand, segments", [
        ("+", [(1000, 1099, 0), (1099, 1199, 2)]),  # -1: leftover 1, phase 2
        ("+", [(1000, 1098, 0), (1098, 1199, 0)]),  # -1: leftover 0, phase 0 (PEG10)
        ("+", [(1000, 1099, 0), (1098, 1199, 2)]),  # -2: leftover 1, phase 2
        ("-", [(1000, 1101, 2), (1101, 1200, 0)]),  # -1 on minus: 5' segment 1101..1200 leftover 1, phase 2
        ("-", [(1000, 1101, 0), (1101, 1199, 0)]),  # -1 on minus: 5' segment 1101..1199 leftover 0, phase 0
    ])
    def test_collapse_cds_keeps_programmed_frameshift_overlaps(self, make_transcript, make_CDS, make_CDS_segment, strand, segments):
        t = make_transcript(feature_id="t1", strand=strand)
        segs = [make_CDS_segment(f"s{i}", strand=strand, start=a, end=b, phase=p) for i, (a, b, p) in enumerate(segments)]
        cds = make_CDS(segments=segs, strand=strand, feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert t.collapsed_CDS_segments is False
        assert [(s.start, s.end, s.phase) for s in cds.CDS_segments] == segments

    def test_collapse_cds_short_overlap_without_phases_kept(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        segs = [make_CDS_segment(f"s{i}", strand="+", start=a, end=b) for i, (a, b) in enumerate([(1000, 1099), (1099, 1199), (1200, 1300)])]
        cds = make_CDS(segments=segs, strand="+", feature_id="t1_CDS1")
        for s in cds.CDS_segments:
            s.phase = None
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        # The 1 bp overlap may be a programmed frameshift; the adjacent row merges
        assert [(s.start, s.end) for s in cds.CDS_segments] == [(1000, 1099), (1099, 1300)]

    def test_collapse_cds_minus_merged_row_takes_phase_of_its_5prime_end(self, make_transcript, make_CDS, make_CDS_segment):
        # Nested rows on the minus strand: the merged row starts (5') at the outer row's end
        t = make_transcript(feature_id="t1", strand="-")
        segs = [make_CDS_segment("outer", strand="-", start=1000, end=1299, phase=2), make_CDS_segment("nested", strand="-", start=1100, end=1200, phase=0)]
        cds = make_CDS(segments=segs, strand="-", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert [(s.start, s.end, s.phase) for s in cds.CDS_segments] == [(1000, 1299, 2)]

    def test_collapse_cds_merges_duplicated_and_nested_rows(self, make_transcript, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        segs = [make_CDS_segment(f"s{i}", strand="+", start=a, end=b, phase=0) for i, (a, b) in enumerate([(1000, 1299), (1000, 1299), (1100, 1200)])]
        cds = make_CDS(segments=segs, strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.collapse_CDS_segments()

        assert [(s.start, s.end, s.phase) for s in cds.CDS_segments] == [(1000, 1299, 0)]

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

class TestDetectPolycistronic:
    def test_single_segment_is_none(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        seg = make_CDS_segment("cds1", strand="+", start=1000, end=1500)
        t.temp_CDSs = [seg]
        assert t._detect_polycistronic() == []
        t.temp_CDSs = []
        assert t._detect_polycistronic() == []

    def test_auto_generated_ids_without_overlap_is_none(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        seg1 = make_CDS_segment("t1_CDS_1000_1200", strand="+", start=1000, end=1200)
        seg2 = make_CDS_segment("t1_CDS_1300_1500", strand="+", start=1300, end=1500)
        t.temp_CDSs = [seg1, seg2]
        assert t._detect_polycistronic() == []

    def test_shared_identical_id_without_overlap_is_none(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        seg1 = make_CDS_segment("cds1", strand="+", start=1000, end=1200)
        seg2 = make_CDS_segment("cds1", strand="+", start=1300, end=1500)
        t.temp_CDSs = [seg1, seg2]
        assert t._detect_polycistronic() == []

    def test_multi_segment_grouped_ids_is_multiple_cdss(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        s1 = make_CDS_segment("cds_a", strand="+", start=1000, end=1200)
        s2 = make_CDS_segment("cds_a", strand="+", start=1300, end=1500)
        s3 = make_CDS_segment("cds_b", strand="+", start=1600, end=1800)
        s4 = make_CDS_segment("cds_b", strand="+", start=1900, end=2100)
        t.temp_CDSs = [s1, s2, s3, s4]
        assert t._detect_polycistronic() == ["possible_polycistronic_transcript"]

    def test_overlapping_distinct_ids_is_multiple_cdss(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        s1 = make_CDS_segment("cds_a", strand="+", start=1000, end=1500)
        s2 = make_CDS_segment("cds_b", strand="+", start=1500, end=2000)
        t.temp_CDSs = [s1, s2]
        assert set(t._detect_polycistronic()) == {"overlapping_CDS_segments", "possible_polycistronic_transcript"}

    def test_overlapping_same_id_is_possible_polycistronic(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        s1 = make_CDS_segment("cds_a", strand="+", start=1000, end=1500)
        s2 = make_CDS_segment("cds_a", strand="+", start=1500, end=2000)
        t.temp_CDSs = [s1, s2]
        assert t._detect_polycistronic() == ["overlapping_CDS_segments"]

    def test_distinct_custom_ids_non_overlapping_is_possible_polycistronic(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        s1 = make_CDS_segment("cds_seg1", strand="+", start=1000, end=1200)
        s2 = make_CDS_segment("cds_seg2", strand="+", start=1300, end=1500)
        t.temp_CDSs = [s1, s2]
        assert t._detect_polycistronic() == ["independent_CDS_segment_ids"]

    def test_canonical_gff3_polycistronic_eden_example(self, make_transcript, make_CDS_segment):
        """
        Canonical GFF3 specification example (mRNA00003 from canonical_gene_example.gff3):
        Two logical CDSs (cds00003: 3301-3902, 5000-5500, 7000-7600 and
        cds00004: 3391-3902, 5000-5500, 7000-7600) encoded under the same mRNA.
        """
        t = make_transcript(feature_id="mRNA00003", strand="+")
        s1 = make_CDS_segment("cds00003", strand="+", start=3301, end=3902)
        s2 = make_CDS_segment("cds00003", strand="+", start=5000, end=5500)
        s3 = make_CDS_segment("cds00003", strand="+", start=7000, end=7600)
        s4 = make_CDS_segment("cds00004", strand="+", start=3391, end=3902)
        s5 = make_CDS_segment("cds00004", strand="+", start=5000, end=5500)
        s6 = make_CDS_segment("cds00004", strand="+", start=7000, end=7600)
        t.temp_CDSs = [s1, s2, s3, s4, s5, s6]
        assert set(t._detect_polycistronic()) == {
            "overlapping_CDS_segments",
            "possible_polycistronic_transcript",
        }


class TestCoordinateEdgeCases:
    def test_generate_CDSs_detects_1bp_boundary_overlap(self, make_transcript, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        # Two distinct CDS segments sharing a single boundary nucleotide (1500)
        seg1 = make_CDS_segment("cds_a", strand="+", start=1000, end=1500)
        seg2 = make_CDS_segment("cds_b", strand="+", start=1500, end=2000)
        t.temp_CDSs = [seg1, seg2]

        assert set(t._detect_polycistronic()) == {"overlapping_CDS_segments", "possible_polycistronic_transcript"}
        t.generate_CDSs(quiet=True)

        assert len(t.CDSs) == 1

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

    @pytest.mark.parametrize("strand,primes", [("+", ("5'", "3'")), ("-", ("3'", "5'"))])
    def test_generate_UTRs_nested_CDS_segments(self, make_transcript, make_exon, make_CDS, make_CDS_segment, strand, primes):
        t = make_transcript(feature_id="t1", strand=strand)
        t.exons = [make_exon("e1", start=1000, end=5000, strand=strand)]
        # Nested segment sorts last but the CDS really ends at 4000
        cs1 = make_CDS_segment("cs1", strand=strand, start=1200, end=4000)
        cs2 = make_CDS_segment("cs2", strand=strand, start=1500, end=1800)
        cds = make_CDS(segments=[cs1, cs2], strand=strand, feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}

        t.generate_UTRs()

        assert [(u.start, u.end, u.prime) for u in cds.UTRs] == [(1000, 1199, primes[0]), (4001, 5000, primes[1])]
        for u in cds.UTRs:
            for cs in cds.CDS_segments:
                assert u.end < cs.start or u.start > cs.end

    def test_generate_UTRs_uses_current_segments(self, make_transcript, make_exon, make_CDS, make_CDS_segment):
        t = make_transcript(feature_id="t1", strand="+")
        t.exons = [make_exon("e1", start=1000, end=5000, strand="+")]
        cs1 = make_CDS_segment("cs1", strand="+", start=1200, end=2000)
        cds = make_CDS(segments=[cs1], strand="+", feature_id="t1_CDS1")
        t.CDSs = {"t1_CDS1": cds}
        # Segment edited after construction without calling CDS.update()
        cs1.end = 4500

        t.generate_UTRs()

        assert [(u.start, u.end) for u in cds.UTRs] == [(1000, 1199), (4501, 5000)]
        assert cds.end == 4500

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



# ============================================================
# UTRs read from the GFF
# ============================================================

class TestReadUTRs:
    """
    Transcript as loaded from a GFF (rows staged in temp_CDSs/temp_UTRs): CDS 101..166 within exon 1..400
    on +, CDS 235..300 on -, so the derived UTRs reach the exon ends while the read ones stop short of them.
    """

    @staticmethod
    def _transcript(make_transcript, make_exon, make_CDS_segment, strand="+", exons=((1, 400),), cds=None, utrs=None, keep=False):
        t = make_transcript(feature_id="t1", strand=strand, start=1, end=400)
        t.exons = [make_exon(f"e{n}", start, end, strand=strand, parents=["t1"]) for n, (start, end) in enumerate(exons, 1)]
        if cds is None:
            cds = [(101, 166)] if strand == "+" else [(235, 300)]
        if utrs is None:
            utrs = [(50, 100), (167, 300)] if strand == "+" else [(101, 234), (301, 350)]
        t.temp_CDSs = [make_CDS_segment("c1", strand=strand, start=start, end=end, parents=["t1"], phase=0) for start, end in cds]
        t.temp_UTRs = [UTR(f"u{n}", "chr1", "aegis", "UTR", strand, start, end, ".", ["t1"]) for n, (start, end) in enumerate(utrs, 1)]
        if keep:
            t.keep_read_UTRs()
        return t

    @staticmethod
    def _utrs(t):
        return [(u.id, u.feature, u.start, u.end) for c in t.CDSs.values() for u in c.UTRs]

    def test_derived_by_default(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment)
        t.update()
        assert self._utrs(t) == [("c1_u1", "five_prime_UTR", 1, 100), ("c1_u2", "three_prime_UTR", 167, 400)]
        assert t.read_UTRs_CDS is None

    def test_derived_UTR_keeps_the_id_of_a_read_one_with_its_coordinates(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, utrs=[(1, 100), (167, 300)])
        t.update()
        assert self._utrs(t) == [("u1", "five_prime_UTR", 1, 100), ("c1_u1", "three_prime_UTR", 167, 400)]

    @pytest.mark.parametrize("strand, expected", [
        ("+", [("u1", "five_prime_UTR", 50, 100), ("u2", "three_prime_UTR", 167, 300)]),
        ("-", [("u1", "three_prime_UTR", 101, 234), ("u2", "five_prime_UTR", 301, 350)]),
    ])
    def test_kept_and_survive_updates(self, make_transcript, make_exon, make_CDS_segment, strand, expected):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, strand=strand, keep=True)
        t.update()
        assert self._utrs(t) == expected
        t.update()
        t.generate_CDSs()
        assert self._utrs(t) == expected
        assert t.temp_UTRs is None

    def test_kept_when_CDS_segments_are_collapsed_first(self, make_transcript, make_exon, make_CDS_segment):
        # Gene.collapse_subfeatures builds the CDS before Transcript.update
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, cds=[(101, 130), (131, 166)], keep=True)
        t.collapse_CDS_segments()
        assert t.collapsed_CDS_segments
        t.update()
        assert self._utrs(t) == [("u1", "five_prime_UTR", 50, 100), ("u2", "three_prime_UTR", 167, 300)]

    def test_renamed_UTRs_keep_their_new_ids(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, keep=True)
        t.update()
        t.rename_utrs(base_id="g1")
        t.update()
        assert [u_id for u_id, *_ in self._utrs(t)] == ["g1_u001", "g1_u002"]

    def test_each_CDS_gets_its_own_list(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, keep=True)
        read = t.temp_UTRs
        t.update()
        assert next(iter(t.CDSs.values())).UTRs is not read

    def test_kept_beyond_the_outer_exons(self, make_transcript, make_exon, make_CDS_segment):
        # Coordinate polishing extends the outer exons to them, as for CDS rows
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, exons=[(60, 300)], utrs=[(50, 100), (167, 350)], keep=True)
        t.update()
        assert self._utrs(t) == [("u1", "five_prime_UTR", 50, 100), ("u2", "three_prime_UTR", 167, 350)]
        assert t.read_UTRs_CDS == (101, 166)

    @pytest.mark.parametrize("exons, utrs, derived", [
        # Overlapping the CDS
        ([(1, 400)], [(50, 110), (167, 300)], [("c1_u1", 1, 100), ("c1_u2", 167, 400)]),
        # Spanning an intron; the derived 5' UTR has the coordinates of the read one, so it keeps its id
        ([(1, 200), (250, 400)], [(1, 100), (167, 260)], [("u1", 1, 100), ("c1_u1", 167, 200), ("c1_u2", 250, 400)]),
        # Within an intron
        ([(1, 200), (250, 400)], [(50, 100), (210, 240)], [("c1_u1", 1, 100), ("c1_u2", 167, 200), ("c1_u3", 250, 400)]),
    ])
    def test_inconsistent_ones_are_derived(self, make_transcript, make_exon, make_CDS_segment, exons, utrs, derived):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, exons=exons, utrs=utrs, keep=True)
        t.update()
        assert t.read_UTRs_CDS is None
        assert [(u_id, start, end) for u_id, _, start, end in self._utrs(t)] == derived

    def test_UTR_on_the_other_strand_is_inconsistent(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, keep=True)
        t.temp_UTRs[1].strand = "-"
        t.update()
        assert t.read_UTRs_CDS is None

    def test_derived_once_the_CDS_changes(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, keep=True)
        t.update()
        c = next(iter(t.CDSs.values()))
        c.CDS_segments[0].end = 157
        c.update_span()
        t.update()
        assert t.read_UTRs_CDS is None
        assert [(u.start, u.end) for u in c.UTRs] == [(1, 100), (158, 400)]

    def test_cleared_UTRs_are_derived(self, make_transcript, make_exon, make_CDS_segment):
        # As when a CDS is reworked (generate_best_protein)
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, keep=True)
        t.update()
        t.clear_UTRs()
        t.update()
        assert t.read_UTRs_CDS is None
        assert [(u.start, u.end) for c in t.CDSs.values() for u in c.UTRs] == [(1, 100), (167, 400)]

    def test_not_kept_without_CDS_rows(self, make_transcript, make_exon, make_CDS_segment):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, cds=[], keep=True)
        assert t.read_UTRs_CDS is None

    @pytest.mark.parametrize("keep", [False, True])
    @pytest.mark.parametrize("collapse_first", [False, True])
    def test_build_the_exons_of_a_transcript_without_exons(self, make_transcript, make_exon, make_CDS_segment, keep, collapse_first):
        t = self._transcript(make_transcript, make_exon, make_CDS_segment, exons=(), cds=[(101, 130), (131, 166)], keep=keep)
        if collapse_first:
            t.collapse_CDS_segments()
        t.update()
        assert [(e.start, e.end) for e in t.exons] == [(50, 300)]
        assert self._utrs(t) == [("u1", "five_prime_UTR", 50, 100), ("u2", "three_prime_UTR", 167, 300)]
        t.update()
        assert self._utrs(t) == [("u1", "five_prime_UTR", 50, 100), ("u2", "three_prime_UTR", 167, 300)]
