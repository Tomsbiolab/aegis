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

    def test_size_counts_overlapping_rows_once(self, make_CDS, make_CDS_segment):
        # A -1 ribosomal slippage: the shared base is read twice but covered once
        cds = make_CDS(segments=[make_CDS_segment("s1", start=1, end=30, phase=0), make_CDS_segment("s2", start=30, end=60, phase=0)])
        assert cds.size == 60
        assert cds.genomic_span == 60

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

    @pytest.mark.parametrize("cds_id, expected", [
        ("AT1G17000:CDS:1", 1),
        ("t1_CDS001", 1),
        ("cds-NP_065184.2", None),  # accession version, not a CDS number
    ])
    def test_id_number_follows_type_name(self, make_CDS, cds_id, expected):
        assert make_CDS(feature_id=cds_id).id_number == expected

    def test_segments_share_the_CDS_id_when_renamed(self, make_CDS, make_CDS_segment):
        segments = [make_CDS_segment("AT1G17000:CDS:1", start=100, end=300), make_CDS_segment("AT1G17000:CDS:2", start=400, end=600)]
        cds = make_CDS(segments=segments, feature_id="AT1G17000:CDS:1")
        cds.rename(base_id="t1", base_gene_id="AT1G17000", count=1, keep_numbering=True)
        assert [cds.id] + [cs.id for cs in cds.CDS_segments] == ["t1_CDS001"] * 3

    def test_unique_segment_ids_are_padded(self, make_CDS, make_CDS_segment):
        segments = [make_CDS_segment("c1", start=100, end=300), make_CDS_segment("c1", start=400, end=600)]
        cds = make_CDS(segments=segments, feature_id="c1")
        cds.main = True
        cds.rename(base_id="t1", base_gene_id="g1", count=1, cds_segment_ids=True)
        assert [cs.id for cs in cds.CDS_segments] == ["t1_CDS001_1", "t1_CDS001_2"]

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

    def test_unsupported_transl_table_ignored(self, make_CDS, make_CDS_segment):
        # Tables 27-31 have context-dependent stops and are not supported
        seg = make_CDS_segment("s1", attributes={"transl_table": "27"})
        assert make_CDS(segments=[seg]).transl_table is None

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
    def _transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, seq, strand, transl_except):
        class DummyScaffold:
            pass
        DummyScaffold.seq = seq
        class DummyGenome:
            scaffolds = {"chr1": DummyScaffold()}
            name = "dummy"
        monkeypatch.setattr(Feature, "_ACTIVE_GENOME", DummyGenome())
        seg = make_CDS_segment("s1", strand=strand, start=1, end=len(seq), phase=0, attributes={"transl_except": transl_except})
        cds = make_CDS(segments=[seg], strand=strand)
        cds.generate_protein(mode="end")
        return cds.protein

    @pytest.mark.parametrize("seq, strand, transl_except", [
        ("ATGTGAAAATAA", "+", "(pos:4..6%2Caa:Sec)"),
        ("TTATTTTCACAT", "-", "(pos:complement(7..9)%2Caa:Sec)"),  # reverse complement of the plus case
    ])
    def test_transl_except_recodes_codon(self, make_CDS, make_CDS_segment, monkeypatch, seq, strand, transl_except):
        p = self._transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, seq, strand, transl_except)
        assert p.seq == "MUK*"
        assert p.internal_stops == 0
        assert p.tags["transl_except"] == "(pos:4..6,aa:Sec)"

    def test_transl_except_readthrough(self, make_CDS, make_CDS_segment, monkeypatch):
        # VEGFA-like: a stop codon read as Ser
        p = self._transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, "ATGTGAAAATAA", "+", "(pos:4..6%2Caa:Ser)")
        assert p.seq == "MSK*"

    @pytest.mark.parametrize("seq, strand, transl_except", [
        ("ATGAAATA", "+", "(pos:7..8%2Caa:TERM)"),
        ("TATTTCAT", "-", "(pos:complement(1..2)%2Caa:TERM)"),
    ])
    def test_transl_except_completes_stop(self, make_CDS, make_CDS_segment, monkeypatch, seq, strand, transl_except):
        # ND1-like: the genome ends in TA, the stop codon is finished by polyadenylation
        p = self._transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, seq, strand, transl_except)
        assert p.seq == "MK*"
        # The polyA tail's A completes the stop codon in nuc_seq; the genomic segments end at the genome's TA
        assert p.nuc_seq == "ATGAAATAA"
        assert p.segments == ((1, 8),)
        assert (p.trimmed_3p, p.partial_3prime) == (0, False)
        assert p.tags["transl_except"] == "(pos:7..8,aa:TERM)"

    def test_transl_except_off_frame_ignored(self, make_CDS, make_CDS_segment, monkeypatch):
        p = self._transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, "ATGTGAAAATAA", "+", "(pos:5..7%2Caa:Sec)")
        assert p.seq == "M*K*"
        assert "transl_except" not in p.tags

    def test_transl_except_start_overrides_initiator_methionine(self, make_CDS, make_CDS_segment, monkeypatch):
        # PTEN-L-like: CTG start read as Leu instead of the initiator Met
        p = self._transl_except_protein(make_CDS, make_CDS_segment, monkeypatch, "CTGAAATAA", "+", "(pos:1..3%2Caa:Leu)")
        assert p.seq == "LK*"

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



class TestCDSMissingPhases:
    # Segments of 10, 11 and 30 bp in transcript order
    @staticmethod
    def _phases(make_CDS, make_CDS_segment, phases, strand="+"):
        coords = [(1, 10), (101, 111), (201, 230)] if strand == "+" else [(221, 230), (111, 121), (1, 30)]
        segments = [make_CDS_segment(f"s{i}", start=start, end=end, strand=strand, phase=phase) for i, ((start, end), phase) in enumerate(zip(coords, phases))]
        cds = make_CDS(segments=sorted(segments, key=lambda s: s.start), strand=strand)
        ordered = cds.CDS_segments if strand == "+" else list(reversed(cds.CDS_segments))
        return [cs.phase for cs in ordered], cds.phase

    @pytest.mark.parametrize("strand", ["+", "-"])
    @pytest.mark.parametrize("annotated, expected", [
        ([None, None, None], ([0, 2, 0], 0)),  # nothing annotated: from phase 0
        ([0, 2, 0], ([0, 2, 0], 0)),           # all annotated: kept
        ([0, None, 1], ([0, 2, 1], 0)),        # read on from the upstream segment; the annotated 1 is kept
        ([None, 1, 1], ([2, 1, 1], 2)),        # derived back: the first segment starts mid-codon
        ([None, 1, None], ([2, 1, 2], 2)),     # both ways from the one annotated phase
    ])
    def test_only_missing_phases_are_derived(self, make_CDS, make_CDS_segment, strand, annotated, expected):
        assert self._phases(make_CDS, make_CDS_segment, annotated, strand) == expected

    def test_kept_phases_report_their_inconsistency(self, make_CDS, make_CDS_segment):
        # The annotated phase 1 of the third segment disagrees with reading on from the second (phase 2 expected)
        segments = [make_CDS_segment("s1", start=1, end=10, phase=0), make_CDS_segment("s2", start=101, end=111, phase=None), make_CDS_segment("s3", start=201, end=230, phase=1)]
        assert make_CDS(segments=segments).inconsistent_junctions() == [(2, False)]


class TestCDSInconsistentJunctions:
    # The first segment in transcript order is 10 bp with phase 0, leaving 1 base over: the next phase must be 2
    @staticmethod
    def _cds(make_CDS, make_CDS_segment, strand, second, second_phase):
        # second: (start, end) of the downstream segment in transcript order
        first = (1, 10) if strand == "+" else (101, 110)
        segments = [make_CDS_segment("s1", start=first[0], end=first[1], strand=strand, phase=0), make_CDS_segment("s2", start=second[0], end=second[1], strand=strand, phase=second_phase)]
        return make_CDS(segments=sorted(segments, key=lambda s: s.start), strand=strand)

    @pytest.mark.parametrize("strand, second", [("+", (11, 30)), ("-", (81, 100))])
    def test_consistent_phases_give_no_junction(self, make_CDS, make_CDS_segment, strand, second):
        cds = self._cds(make_CDS, make_CDS_segment, strand, second, 2)
        assert cds.inconsistent_junctions() == []
        assert cds.has_internal_shift is False

    @pytest.mark.parametrize("strand, second, within_exon", [
        ("+", (11, 30), True),    # abutting
        ("+", (12, 30), True),    # 1 bp gap
        ("+", (13, 30), False),   # 2 bp gap
        ("+", (9, 30), True),     # overlap
        ("+", (500, 530), False), # intron
        ("-", (81, 100), True),
        ("-", (81, 99), True),
        ("-", (81, 98), False),
        ("-", (81, 102), True),
        ("-", (1, 30), False),
    ])
    def test_inconsistent_junction_within_exon_or_across_intron(self, make_CDS, make_CDS_segment, strand, second, within_exon):
        cds = self._cds(make_CDS, make_CDS_segment, strand, second, 0)
        assert cds.inconsistent_junctions() == [(1, within_exon)]
        assert cds.has_internal_shift is True

    def test_junctions_are_indexed_in_transcript_order(self, make_CDS, make_CDS_segment):
        # On the - strand the downstream segment of the inconsistent junction is the lowest one
        segments = [
            make_CDS_segment("s3", start=1, end=20, strand="-", phase=0),
            make_CDS_segment("s2", start=41, end=60, strand="-", phase=2),
            make_CDS_segment("s1", start=101, end=110, strand="-", phase=0),
        ]
        cds = make_CDS(segments=segments, strand="-")
        # s1 leaves 1 base, completed by s2 (phase 2); s2 leaves (20 - 2) % 3 = 0, but s3 has phase 0: consistent
        assert cds.inconsistent_junctions() == []
        cds.CDS_segments[0].phase = 1
        assert cds.inconsistent_junctions() == [(2, False)]

    def test_junctions_next_to_a_segment_without_phase_are_not_judged(self, make_CDS, make_CDS_segment):
        segments = [
            make_CDS_segment("s1", start=1, end=10, phase=0),
            make_CDS_segment("s2", start=20, end=40, phase=2),
            make_CDS_segment("s3", start=50, end=70, phase=0),
        ]
        cds = make_CDS(segments=segments)
        assert cds.inconsistent_junctions() == [(2, False)]
        # Unknown phases are derived when the CDS is built; set afterwards, the junctions next to it are skipped
        cds.CDS_segments[1].phase = None
        assert cds.inconsistent_junctions() == []


# ============================================================
# Exon
# ============================================================

class TestSubfeatureIdNumbers:
    @pytest.mark.parametrize("cls, feature_id, expected", [
        (Exon, "AT4G00005:exon:1", 1),
        (Exon, "g1_e003", 3),
        (Exon, "my_exon", None),
        (Exon, "exon-NM_000384.3-26", None),  # NCBI exon ids are numbered per transcript after the accession
        (UTR, "AT1G16460:five_prime_UTR:7", 7),
        (UTR, "t1_CDS1_u2", 2),
        (UTR, "g1_u001", 1),
    ])
    def test_id_number_follows_type_name(self, cls, feature_id, expected):
        assert cls(feature_id, "chr1", "test", "exon", "+", 1, 10, ".", ["t1"]).id_number == expected


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
