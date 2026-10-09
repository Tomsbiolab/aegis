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

    def test_update_phase_minus_strand_descending_input(self, make_CDS, make_CDS_segment):
        # Exact coordinates corresponding to Vitvi000006 on minus strand
        # Provided in descending (5' to 3') order
        seg_5p = make_CDS_segment("cds1", start=106282, end=106588, strand="-")  # length = 307, (307-0)%3 = 1 leftover
        seg_mid = make_CDS_segment("cds2", start=106094, end=106179, strand="-") # length = 86, phase must be 3-1 = 2
        seg_3p = make_CDS_segment("cds3", start=105672, end=105956, strand="-")  # length = 285, phase must be 0

        # Pass in descending biological order
        cds = make_CDS(segments=[seg_5p, seg_mid, seg_3p], strand="-")
        
        # Segments must be sorted in ascending genomic coordinate order
        assert cds.CDS_segments[0].start == 105672
        assert cds.CDS_segments[1].start == 106094
        assert cds.CDS_segments[2].start == 106282

        # In ascending order, segment [2] is 5' start, [1] is mid, [0] is 3' end
        assert cds.CDS_segments[2].phase == 0  # 5' start codon segment must have phase 0
        assert cds.CDS_segments[1].phase == 2  # (3 - 1) = 2
        assert cds.CDS_segments[0].phase == 0  # (86 - 2) % 3 = 0 leftover, so phase = 0

    def test_cds_init_sorts_unordered_segments(self, make_CDS, make_CDS_segment):
        s1 = make_CDS_segment("s1", start=5000, end=6000)
        s2 = make_CDS_segment("s2", start=1000, end=2000)
        s3 = make_CDS_segment("s3", start=3000, end=4000)
        cds = make_CDS(segments=[s1, s2, s3])
        assert [s.start for s in cds.CDS_segments] == [1000, 3000, 5000]
        assert cds.start == 1000
        assert cds.end == 6000

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

    def test_cds_phase_preservation_when_valid(self, make_CDS, make_CDS_segment):
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="+", phase=1)
        seg2 = make_CDS_segment("seg2", start=2000, end=2199, strand="+", phase=2)
        cds = make_CDS(segments=[seg1, seg2], strand="+")
        assert cds.phase == 1
        assert cds.CDS_segments[0].phase == 1
        assert cds.CDS_segments[1].phase == 2
        cds.update_phase(override=False)
        assert cds.CDS_segments[0].phase == 1
        assert cds.CDS_segments[1].phase == 2

    def test_cds_phase_recalculation_with_override(self, make_CDS, make_CDS_segment):
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="+", phase=1)  # size 100
        seg2 = make_CDS_segment("seg2", start=2000, end=2199, strand="+", phase=2)  # size 200
        cds = make_CDS(segments=[seg1, seg2], strand="+")
        cds.update_phase(override=True)
        # seg1 size 100, phase 1: leftover = (100 - 1) % 3 = 0.
        # seg2 phase should be recalculated to 0:
        assert cds.CDS_segments[0].phase == 1
        assert cds.CDS_segments[1].phase == 0

    def test_cds_phase_recalculation_with_full_override(self, make_CDS, make_CDS_segment):
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="+", phase=1)  # size 100
        seg2 = make_CDS_segment("seg2", start=2000, end=2199, strand="+", phase=2)  # size 200
        cds = make_CDS(segments=[seg1, seg2], strand="+")
        cds.update_phase(full_override=True)
        # With full_override=True, seg1 (5') phase is reset to 0
        # seg1 size 100, phase 0 -> leftover = (100 - 0) % 3 = 1 -> seg2 phase = 3 - 1 = 2
        assert cds.CDS_segments[0].phase == 0
        assert cds.CDS_segments[1].phase == 2
        assert cds.phase == 0

    def test_cds_phase_minus_strand_preservation(self, make_CDS, make_CDS_segment):
        # seg1: 1000..1099 (3' segment on minus strand), phase 2
        # seg2: 2000..2099 (5' segment on minus strand), phase 1
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="-", phase=2)
        seg2 = make_CDS_segment("seg2", start=2000, end=2099, strand="-", phase=1)
        cds = make_CDS(segments=[seg1, seg2], strand="-")
        assert cds.phase == 1
        assert cds.CDS_segments[1].phase == 1
        assert cds.CDS_segments[0].phase == 2

    def test_cds_relative_coding_intervals_plus_strand(self, make_CDS, make_CDS_segment, monkeypatch):
        from aegis.misc_features import Protein
        # Two segments on + strand: 1000..1099 (size 100), 2000..2099 (size 100) -> total size 200
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="+")
        seg2 = make_CDS_segment("seg2", start=2000, end=2099, strand="+")
        cds = make_CDS(segments=[seg1, seg2], strand="+")
        # Suppose protein spans 1010..1099 and 2000..2050
        cds.protein = Protein(
            prot_id="p1", sequence="M" * 47, chrom="chr1",
            start=1010, end=2050, readthrough="end",
            nuc_seq="ATG" * 47,
            segments=((1010, 1099), (2000, 2050))
        )
        intervals = cds.relative_coding_intervals
        assert intervals == ((10, 99), (100, 150))
        assert cds.relative_coding_start == 10
        assert cds.relative_coding_end == 150

        # Verify slice reconstruction matches nuc_seq
        class DummyScaffold:
            seq = "N" * 1000 + "A" * 10 + "ATG" * 30 + "N" * 900 + "ATG" * 17 + "C" * 49
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
            name = "dummy"
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        reconstructed = "".join(cds.seq[s:e+1] for s, e in intervals)
        assert len(reconstructed) == len(cds.protein.nuc_seq)

    def test_cds_relative_coding_intervals_minus_strand(self, make_CDS, make_CDS_segment, monkeypatch):
        from aegis.misc_features import Protein
        # Two segments on - strand: 1000..1099 (3'), 2000..2099 (5') -> total size 200
        seg1 = make_CDS_segment("seg1", start=1000, end=1099, strand="-")
        seg2 = make_CDS_segment("seg2", start=2000, end=2099, strand="-")
        cds = make_CDS(segments=[seg1, seg2], strand="-")
        # In transcription direction (5' to 3'), seg2 is first (offset 0..99), seg1 is second (offset 100..199)
        # Protein spans 2000..2090 (5' piece, length 91) and 1020..1099 (3' piece, length 80)
        cds.protein = Protein(
            prot_id="p1", sequence="M" * 57, chrom="chr1",
            start=1020, end=2090, readthrough="end",
            nuc_seq="ATG" * 57,
            segments=((1020, 1099), (2000, 2090))
        )
        intervals = cds.relative_coding_intervals
        # In seg2 (offset 0): 2000..2090 -> s = 0 + (2099 - 2090) = 9, e = 0 + (2099 - 2000) = 99
        # In seg1 (offset 100): 1020..1099 -> s = 100 + (1099 - 1099) = 100, e = 100 + (1099 - 1020) = 179
        assert intervals == ((9, 99), (100, 179))
        assert cds.relative_coding_start == 9
        assert cds.relative_coding_end == 179

        class DummyScaffold:
            seq = "N" * 1000 + "T" * 80 + "G" * 20 + "N" * 900 + "A" * 91 + "C" * 9
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
            name = "dummy"
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        reconstructed = "".join(cds.seq[s:e+1] for s, e in intervals)
        assert len(reconstructed) == len(cds.protein.nuc_seq)



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
        assert u.prime == "5'"
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

    def test_splice_site_case_insensitivity(self, monkeypatch):
        i = Intron(
            feature_id="intron1",
            ch="chr1",
            source="aegis",
            feature="intron",
            strand="+",
            start=1,
            end=10,
            score="."
        )
        class DummyScaffold:
            seq = "gtaaaaccag"
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
        
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        assert i.splice_site_donor == "GT"
        assert i.splice_site_acceptor == "AG"
        assert i.boundary == "GT-AG"
        assert i.canonical is True
