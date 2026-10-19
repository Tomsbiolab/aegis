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

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_span_nested_segments(self, make_CDS, make_CDS_segment, strand):
        # Sorted by (start, end), the nested segment comes last but does not hold the latest end
        segments = [
            make_CDS_segment("s1", strand=strand, start=1200, end=4000),
            make_CDS_segment("s2", strand=strand, start=1500, end=1800),
        ]
        cds = make_CDS(segments=segments, strand=strand)
        assert cds.CDS_segments[-1].end == 1800
        assert (cds.start, cds.end) == (1200, 4000)

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_span_frameshift_overlap(self, make_CDS, make_CDS_segment, strand):
        segments = [
            make_CDS_segment("s1", strand=strand, start=100, end=400),
            make_CDS_segment("s2", strand=strand, start=399, end=700),
        ]
        cds = make_CDS(segments=segments, strand=strand)
        assert (cds.start, cds.end) == (100, 700)

    def test_span_unsorted_segments(self, make_CDS, make_CDS_segment):
        segments = [
            make_CDS_segment("s2", start=3000, end=4500),
            make_CDS_segment("s1", start=1200, end=2000),
        ]
        cds = make_CDS(segments=segments)
        assert [s.start for s in cds.CDS_segments] == [1200, 3000]
        assert (cds.start, cds.end) == (1200, 4500)

    def test_update_recomputes_span(self, make_CDS):
        cds = make_CDS()
        cds.CDS_segments[0].start = 1100
        cds.CDS_segments[0].end = 5000  # now contains the second segment
        cds.update()
        assert (cds.start, cds.end) == (1100, 5000)

    def test_update_span_only_touches_coordinates(self, make_CDS, make_CDS_segment):
        segments = [
            make_CDS_segment("s1", start=1200, end=4000, phase=0),
            make_CDS_segment("s2", start=1500, end=1800, phase=2),
        ]
        cds = make_CDS(segments=segments)
        cds.start, cds.end = 1, 2
        cds.update_span()
        assert (cds.start, cds.end) == (1200, 4000)
        assert [s.phase for s in cds.CDS_segments] == [0, 2]

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

    def test_translation_attributes_read_from_segments(self, make_CDS, make_CDS_segment):
        attributes = {"transl_table": "2", "exception": "ribosomal slippage", "Note": "x"}
        segments = [
            make_CDS_segment("s1", start=100, end=300, attributes={**attributes, "transl_except": "(pos:100..102%2Caa:Sec)"}),
            make_CDS_segment("s2", start=400, end=600, attributes={**attributes, "transl_except": "(pos:400..402%2Caa:Sec),(pos:100..102%2Caa:Sec)"}),
        ]
        cds = make_CDS(segments=segments)
        assert cds.transl_table == 2
        assert cds.exception == "ribosomal slippage"
        # Entries repeated across rows are merged
        assert cds.transl_except == ((100, 102, "U"), (400, 402, "U"))

    def test_translation_attributes_decoded(self, make_CDS, make_CDS_segment):
        seg = make_CDS_segment("s1", attributes={"exception": "alternative start codon%2C annotated by transcript or proteomic data"})
        assert make_CDS(segments=[seg]).exception == "alternative start codon, annotated by transcript or proteomic data"

    def test_translation_attributes_absent(self, make_CDS):
        cds = make_CDS()
        assert (cds.transl_table, cds.transl_except, cds.exception) == (None, (), None)

    def test_corrected_CDS_segments_keep_attributes(self, make_CDS, make_CDS_segment, monkeypatch):
        class DummyScaffold:
            seq = "GG" + "ATGAAATAA" + "GG"
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
            name = "dummy"
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        seg = make_CDS_segment("s1", start=1, end=13, phase=0, attributes={"transl_table": "1", "Note": "x"})
        cds = make_CDS(segments=[seg])
        cds.generate_protein(mode="orf", correct_CDS=True)
        assert (cds.CDS_segments[0].start, cds.CDS_segments[0].end) == (3, 11)
        assert cds.CDS_segments[0].misc_attributes == ["transl_table=1", "Note=x"]

    @staticmethod
    def _shift_protein(make_CDS, make_CDS_segment, monkeypatch, strand, segments, adjust_internal_shifts, mode="end", seq="ACGGCA" * 10, attributes=None, enforce_start_codon=False):
        # Default sequence has no T, so no stop codons in any frame
        class DummyScaffold:
            pass
        DummyScaffold.seq = seq
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
            name = "dummy"
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        segs = [make_CDS_segment(f"s{i}", strand=strand, start=a, end=b, phase=p, attributes=attributes) for i, (a, b, p) in enumerate(segments)]
        cds = make_CDS(segments=segs, strand=strand)
        cds.generate_protein(mode=mode, enforce_start_codon=enforce_start_codon, adjust_internal_shifts=adjust_internal_shifts)
        return cds.protein, seq

    # Phase mismatch (leftover 1 + phase 0) both at a contiguous junction (1..10 | 11..20) and across an intron (11..20 | 31..42)
    SHIFT_SEGMENTS_PLUS = [(1, 10, 0), (11, 20, 0), (31, 42, 0)]

    def test_intra_exon_shift_adjusts_only_contiguous_junction(self, make_CDS, make_CDS_segment, monkeypatch):
        p, seq = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", self.SHIFT_SEGMENTS_PLUS, "intra_exon")
        # Base 10 dropped at the contiguous junction; base 20 opens the first codon across the intron
        assert p.nuc_seq == seq[0:9] + seq[10:20] + seq[30:41]
        assert p.segments == ((1, 9), (11, 20), (31, 41))
        assert p.frameshifts == (1,)

    def test_all_shift_adjusts_every_junction(self, make_CDS, make_CDS_segment, monkeypatch):
        p, seq = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", self.SHIFT_SEGMENTS_PLUS, "all")
        assert p.nuc_seq == seq[0:9] + seq[10:19] + seq[30:42]
        assert p.segments == ((1, 9), (11, 19), (31, 42))
        assert p.frameshifts == (1, 1)

    def test_intra_exon_shift_adjusts_only_contiguous_junction_minus(self, make_CDS, make_CDS_segment, monkeypatch):
        # Mirror of the plus case: 5' segment 33..42, contiguous 23..32, then intron down to 1..12
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "-", [(1, 12, 0), (23, 32, 0), (33, 42, 0)], "intra_exon")
        assert p.segments == ((2, 12), (23, 32), (34, 42))
        assert p.frameshifts == (1,)

    @pytest.mark.parametrize("strand, segments", [
        # PEG10-like -1 frameshift: 1-bp overlap with phases consistent with reading the shared base twice
        ("+", [(1, 10, 0), (10, 21, 2)]),
        ("-", [(10, 21, 2), (21, 30, 0)]),
    ])
    def test_annotated_overlap_is_negative_shift(self, make_CDS, make_CDS_segment, monkeypatch, strand, segments):
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, strand, segments, "intra_exon")
        assert len(p.nuc_seq) == 21  # 10 + 12 - 1 trailing base: the shared base is translated twice
        assert p.frameshifts == (-1,)

    def test_annotated_gap_is_positive_shift(self, make_CDS, make_CDS_segment, monkeypatch):
        # OAZ1-like +1 frameshift: base 11 skipped, phases consistent with the skip
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", [(1, 10, 0), (12, 22, 2)], "intra_exon")
        assert p.frameshifts == (1,)

    def test_shift_outside_ORF_not_reported(self, make_CDS, make_CDS_segment, monkeypatch):
        # The +1 gap at base 11 lies upstream of the ORF (ATG at 31), so the ORF protein has no shift
        seq = "ACGGCA" * 5 + "ATGAAACCCGGGTAA" + "ACGGCA" * 2 + "ACG"
        segments = [(1, 10, 0), (12, 60, 2)]
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", segments, "intra_exon", seq=seq)
        assert p.frameshifts == (1,)
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", segments, "intra_exon", mode="orf", seq=seq, enforce_start_codon=True)
        assert p.start == 31
        assert p.frameshifts == ()

    def test_exception_reaches_protein(self, make_CDS, make_CDS_segment, monkeypatch):
        p, _ = self._shift_protein(make_CDS, make_CDS_segment, monkeypatch, "+", [(1, 10, 0), (12, 22, 2)], "intra_exon", attributes={"exception": "ribosomal slippage"})
        assert p.exception == "ribosomal slippage"
        assert p.tags["exception"] == "ribosomal slippage"



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
