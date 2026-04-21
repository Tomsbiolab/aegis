"""
Tests for aegis.genefunctions — pure utility functions.
"""

import pytest
from unittest.mock import patch, PropertyMock

from pathlib import Path
from aegis.utils.genefunctions import (
    reverse_complement,
    find_ORFs,
    choose_orf,
    translate,
    trim_surplus
)

from aegis.feature import Feature
from aegis.utils.gtf_gff import parse_gff_parts, parse_gff_attributes
from aegis.utils.misc import count_occurrences, find_all_occurrences


# ============================================================
# parse_gff_attributes
# ============================================================

class TestParseGffAttributes:
    def test_standard_attributes(self):
        attrs = parse_gff_attributes("ID=gene1;Name=MyGene;Note=some note;Symbol=1,2;Alias=A1,A45")
        assert attrs["ID"] == "gene1"
        assert attrs["Name"] == ["MyGene"]
        assert attrs["Alias"] == ["A1", "A45"]
        assert attrs["Symbol"] == ["1", "2"]
        assert attrs["Note"] == "some note"

    def test_parent_attribute_splits_into_list(self):
        attrs = parse_gff_attributes("ID=exon1;Name=E123,E2;Parent=mRNA1,mRNA2")
        assert attrs["Name"] == ["E123", "E2"]
        assert attrs["Parent"] == ["mRNA1", "mRNA2"]

    def test_empty_string(self):
        assert parse_gff_attributes("") == {}

    def test_dot_string(self):
        assert parse_gff_attributes(".") == {}

    def test_single_attribute(self):
        attrs = parse_gff_attributes("ID=feat1")
        assert attrs["ID"] == "feat1"

    def test_derives_from_treated_as_parent(self):
        attrs = parse_gff_attributes("ID=cds1;Derives_from=mRNA1")
        assert attrs["Parent"] == ["mRNA1"]


# ============================================================
# parse_gff_parts
# ============================================================

class TestParseGffLine:
    LINES_FILE = Path(__file__).resolve().parent / "test_data" / "input" / "annotation" / "genefunctions_lines.gff3"

    @classmethod
    def _read_line(cls, index):
        """Read a specific line (0-indexed) from genefunctions_lines.gff3."""
        with open(cls.LINES_FILE) as f:
            lines = f.read().splitlines()
        return lines[index]

    def test_basic_parsing(self):
        line = self._read_line(0)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.ch == "chr1"
        assert entry.source == "aegis"
        assert entry.feature == "gene"
        assert entry.start == 1000
        assert entry.end == 5000
        assert entry.strand == "+"
        assert entry.id == "gene1"
        assert entry.pseudogene is False
        assert entry.transposable is False
        assert entry.decreasing_coordinates is False

    def test_pseudogene_detected(self):
        line = self._read_line(1)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.pseudogene is True

    def test_transposable_by_feature(self):
        line = self._read_line(2)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.transposable is True

    def test_transposable_by_attribute(self):
        line = self._read_line(3)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.transposable is True

    def test_decreasing_coordinates_swapped(self):
        line = self._read_line(4)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.decreasing_coordinates is True
        assert entry.start == 1000
        assert entry.end == 5000

    def test_multi_parent(self):
        line = self._read_line(5)
        line = line.strip().split("\t")
        entry = parse_gff_parts(line)
        assert entry.parents == ["mRNA1", "mRNA2"]



# ============================================================
# count_occurrences
# ============================================================

class TestCountOccurrences:
    def test_basic_count(self):
        assert count_occurrences("AATTGGCC", "A") == 2
        assert count_occurrences("AATTGGCC", "G") == 2

    def test_missing_character(self):
        assert count_occurrences("AATTGGCC", "X") == 0

    def test_empty_string(self):
        assert count_occurrences("", "A") == 0


# ============================================================
# find_all_occurrences
# ============================================================

class TestFindAllOccurrences:
    def test_multiple_matches(self):
        matches = find_all_occurrences("AT", "ATATGCAT")
        # AT at pos 0, 2, 6
        assert len(matches) >= 2
        assert all(m[2] == "AT" for m in matches)

    def test_no_matches(self):
        matches = find_all_occurrences("ZZZ", "ATGCATGC")
        assert matches == []


# ============================================================
# reverse_complement
# ============================================================

class TestReverseComplement:
    def test_simple_sequence(self):
        # ATGC -> reverse complement = GCAT
        assert reverse_complement("ATGC") == "GCAT"

    def test_palindrome(self):
        # ATAT -> reverse complement = ATAT
        assert reverse_complement("ATAT") == "ATAT"

    def test_longer_sequence(self):
        seq = "ATGCGATCGA"
        rc = reverse_complement(seq)
        assert reverse_complement(rc) == seq  # double RC == original


# ============================================================
# find_ORFs / longest_ORF
# ============================================================

class TestFindORFs:

    def test_orfs(self):
        # ATG + 6 codons + TAA  = 24 nt
        seq = "AGATATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAGA"
        orfs = find_ORFs(seq, must_have_stop=True, readthrough_stop=False)
        assert len(orfs) == 1
        # the ORF should start with ATG and end with a stop codon
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "TAA"
        assert orfs[0][1] == 4
        assert orfs[0][2] == 27
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAA"

        orfs = find_ORFs(seq, must_have_stop=True, readthrough_stop=True)
        assert len(orfs) == 2
        # the ORF should start with ATG and end with a stop codon
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "TAA"
        assert orfs[0][1] == 4
        assert orfs[0][2] == 27
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAA"

        orf_seq = orfs[1][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "TAG"
        assert orfs[1][1] == 4
        assert orfs[1][2] == 36
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAAAAAGATTAG"

        orfs = find_ORFs(seq, must_have_stop=False, readthrough_stop=False)
        assert len(orfs) == 7
        
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "AAA"
        assert orfs[0][1] == 4
        assert orfs[0][2] == 9
        assert orf_seq == "ATGAAA"

        orf_seq = orfs[-1][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "TAA"
        assert orfs[-1][1] == 4
        assert orfs[-1][2] == 27
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAA"

        orfs = find_ORFs(seq, must_have_stop=False, readthrough_stop=True)
        assert len(orfs) == 11
        
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "AAA"
        assert orfs[0][1] == 4
        assert orfs[0][2] == 9
        assert orf_seq == "ATGAAA"

        orf_seq = orfs[-1][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "AAG"
        assert orfs[-1][1] == 4
        assert orfs[-1][2] == 39
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAG"

        orfs = find_ORFs(seq, must_have_stop=False, readthrough_stop=False, min_codon_len=1)
        assert len(orfs) == 8
        
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "ATG"
        assert orfs[0][1] == 4
        assert orfs[0][2] == 6
        assert orf_seq == "ATG"

        seq = "AGATATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAGAA"

        orfs = find_ORFs(seq, must_have_stop=False, readthrough_stop=True)
        orf_seq = orfs[-1][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "AAG"
        assert orfs[-1][1] == 4
        assert orfs[-1][2] == 39
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAG"

        seq = "AGATATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAGAAA"

        orfs = find_ORFs(seq, must_have_stop=False, readthrough_stop=True)
        orf_seq = orfs[-1][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] == "AAA"
        assert orfs[-1][1] == 4
        assert orfs[-1][2] == 42
        assert orf_seq == "ATGAAACCCGGGTTGATTAACTAAAAAGATTAGAAGAAA"

    def test_no_start_codon(self):
        seq = "AAACCCGGGTTTTTTTAA"
        orfs = find_ORFs(seq)
        assert orfs == []

    def test_no_stop_codon_without_must_have_stop(self):
        seq = "ATGAAACCCGGG"
        orfs = find_ORFs(seq, must_have_stop=False)
        assert len(orfs) >= 1


class TestLongestORF:
    def test_single_orf(self):
        orfs = [("ATGAAATAA", 0, 8)]
        result = choose_orf(orfs, mode="longest")
        assert result[0] == "ATGAAATAA"

    def test_multiple_orfs(self):
        orfs = [
            ("ATGTAA", 0, 5),
            ("ATGAAAAAATAA", 10, 21),
            ("ATGCCCTAA", 30, 38),
        ]
        result = choose_orf(orfs, mode="longest")
        assert result[0] == "ATGAAAAAATAA"

    def test_empty_list(self):
        result = choose_orf([], mode="longest")
        assert result[0] == ""


# ============================================================
# trim_surplus
# ============================================================

class TestTrimSurplus:
    def test_divisible_by_3(self):
        seq = "ATGAAATAA"  # 9 nt, already divisible
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is False
        assert out == seq

    def test_surplus_trimmed(self):
        seq = "ATGAAATAAG"  # 10 nt, surplus = 1
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
    
    def test_surplus_other(self):

        seq = "AGGAGATGTAAGATGAT"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "AGGAGATGTAAGATG"

        seq = "AAATGTAA"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "ATGTAA"

        seq = "AATGTAA"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "ATGTAA"

        seq = "AAATGTAAAA"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "ATGTAA"

        seq = "AAATGTAAAAA"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "ATGTAA"

        seq = "AAAATGTAAAA"
        out, surplus, _, _ = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0
        assert out == "ATGTAA"


# ============================================================
# translate
# ============================================================

class TestTranslate:
    """Test the translate() utility function."""

    def test_simple_orf(self):
        # ATG AAA TAA = M K *
        assert translate("ATGAAATAA") == "MK*"

    def test_no_start_codon(self):
        # GGG AAA TAA = G K *
        prot = translate("GGGAAATAA")
        assert prot[0] != "M"
        assert prot == "GK*"

    def test_late_start_in_frame(self):
        # GGG ATG TAA = G M *
        prot = translate("GGGATGTAA")
        assert prot[0] != "M"
        assert "M" in prot[1:] # ATG is in-frame but not first codon

    def test_atg_not_in_frame(self):
        # GGA TGA TAA = G * * (the ATG straddles codons 1 and 2)
        prot = translate("GGATGATAA")
        assert "M" not in prot # no methionine at any codon boundary

    def test_ambiguous_codons_produce_gap(self):
        # N in sequence signals ambiguity
        # ATG NNA TAA = M - * (NN + A can't resolve to one AA)
        prot = translate("ATGNNATAA")
        assert prot[0] == "M"
        assert "-" in prot # ambiguous codon rendered as gap
        assert prot[-1] == "*"

    def test_longer_sequence(self):
        # ATG AAA CCC GGG TTT AA = M K P G F
        prot = translate("ATGAAACCCGGGTTTAA")
        # 17 nt is NOT a multiple of 3, but translate just stops when bytes run out
        assert len(prot) == 5


# ============================================================
# translate + trim_surplus
# ============================================================

class TestTranslatePipeline:
    """
    Test trim_surplus -> translate.
    """

    def test_orf_or_end_extracts_orf(self):
        # GGG ATG AAA TAA (12 nt, multiple of 3)
        seq = "GGGATGAAATAA"
        trimmed, surplus, cs, ce = trim_surplus(seq, mode="orf_or_end")
        assert trimmed == "ATGAAATAA"
        assert surplus is False
        prot = translate(trimmed)
        assert prot == "MK*"

    def test_orf_or_end_with_surplus(self):
        # 13 nt = surplus flag should be True
        seq = "GGGGATGAAATAA"
        trimmed, surplus, cs, ce = trim_surplus(seq, mode="orf_or_end")
        assert surplus is True
        prot = translate(trimmed)
        assert prot.startswith("M")
        assert prot.endswith("*")

    def test_end_mode_trims_from_3prime(self):
        # 10 nt = end-trim removes 1 nt
        seq = "ATGAAATAAG"
        trimmed, surplus, _, _ = trim_surplus(seq, mode="end")
        assert len(trimmed) % 3 == 0
        assert surplus is True
        prot = translate(trimmed)
        assert prot == "MK*"

    def test_no_atg_in_sequence(self):
        # No ATG = orf mode yields no ORF, fallback to end-trimming
        seq = "GGGCCCAAA"
        trimmed, surplus, _, _ = trim_surplus(seq, mode="orf_or_end")
        prot = translate(trimmed)
        assert "M" not in prot
        assert prot == "GPK"

    def test_internal_stop(self):
        # ATG AAA TAA CCC GGG (15 nt)
        seq = "ATGAAATAACCCGGG"
        trimmed, _, _, _ = trim_surplus(seq, mode="orf_or_end")
        assert trimmed == "ATGAAATAA"
        prot = translate(trimmed)
        assert prot == "MK*"

    def test_internal_stop_with_large_trim_fallback(self):
        # "GGG ATG TAA CCC GGG" = 15 nt, trim difference > 6 nt so fallback to end-trimming
        seq = "GGGATGTAACCCGGG"
        trimmed, _, _, _ = trim_surplus(seq, mode="orf_or_end")
        assert len(trimmed) == 15
        prot = translate(trimmed)
        assert "*" in prot[:-1]

    def test_start_mode_trims_from_5prime(self):
        # "AATGTAA" = start-trim removes 1 nt from 5' to get mult-of-3
        seq = "AATGTAA"
        trimmed, surplus, _, _ = trim_surplus(seq, mode="start")
        assert len(trimmed) % 3 == 0
        assert surplus is True
        assert trimmed.startswith("ATG")


# ============================================================
# CDS.generate_protein
# ============================================================

class TestGenerateProtein:
    """
    Test CDS.generate_protein() and resulting Protein properties.

    We use unittest.mock to patch the CDS.seq property so it returns
    a known nucleotide string without needing a genome loaded.
    """

    @staticmethod
    def _make_cds(make_CDS_segment, make_CDS, seq_len: int, strand: str = "+"):
        """Build a CDS spanning *seq_len* nucleotides."""
        seg = make_CDS_segment("seg1", start=1000, end=1000 + seq_len - 1, strand=strand)
        return make_CDS(segments=[seg], strand=strand)

    def test_standard_protein(self, make_CDS_segment, make_CDS):
        # ATG AAA TAA  ->  M K *
        seq = "ATGAAATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.seq == "MK*"
        assert p.ATG_start is True
        assert p.end_stop is True
        assert p.early_stop is False
        assert p.gaps is False
        assert p.partial is False
        assert p.truncated is False

    def test_no_start_codon(self, make_CDS_segment, make_CDS):
        # GGG AAA TAA = G K *
        seq = "GGGAAATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.ATG_start is False
        assert p.end_stop is True
        assert p.partial is True # no start

    def test_late_start_in_frame(self, make_CDS_segment, make_CDS):
        # GGG ATG TAA = G M *
        seq = "GGGATGTAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.ATG_start is False # first codon is not M
        assert p.ATG_late is True # M exists later in frame
        assert p.partial is True # no start

    def test_atg_not_in_frame(self, make_CDS_segment, make_CDS):
        # GGA TGA TAA = G * * (ATG straddles codon boundary)
        seq = "GGATGATAA" 
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.ATG_start is False
        assert p.ATG_late is False # no M anywhere in protein
        assert p.partial is True

    def test_ambiguous_codons(self, make_CDS_segment, make_CDS):
        # ATG NNA TAA = M - *
        seq = "ATGNNATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.gaps is True
        assert p.partial is True # gaps = partial

    def test_orf_extraction(self, make_CDS_segment, make_CDS):
        # GGG ATG AAA TAA (12 nt, ORF trim < 6)
        seq = "GGGATGAAATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="orf_or_end")
        p = cds.protein
        assert p.seq == "MK*"
        assert p.ATG_start is True
        assert p.end_stop is True
        assert p.early_stop is False

    def test_early_stop(self, make_CDS_segment, make_CDS):
        # ATG TAA AAA TAA = M * K *
        seq = "ATGTAAAAATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="end")
        p = cds.protein
        assert p.early_stop is True
        assert p.end_stop is True
        assert p.truncated is True

    def test_surplus_flagged(self, make_CDS_segment, make_CDS):
        seq = "GGGGATGAAATAA"
        cds = self._make_cds(make_CDS_segment, make_CDS, len(seq))
        with patch.object(type(cds), "seq", new_callable=PropertyMock, return_value=seq):
            cds.generate_protein(mode="orf_or_end")
        p = cds.protein
        assert p.nucleotide_surplus is True
        assert p.seq.startswith("M")
        assert p.seq.endswith("*")


# ============================================================
# overlap
# ============================================================

class TestOverlap:
    """Test the overlap function using simple features."""

    def test_overlapping_features(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t300\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t200\t400\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 101  # 200..300 inclusive

    def test_overlapping_features_displaced(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t500\t6000\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t2000\t8000\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 4001

    def test_non_overlapping_features(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t200\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t300\t400\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is False
        assert bp == 0

    def test_adjacent_features(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t200\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t201\t300\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is False
        assert bp == 0

    def test_small_overlap(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t200\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t200\t300\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 1

    def test_contained_feature(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t500\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t200\t300\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 101

    def test_identical_features(self, create_test_feature):
        f1 = create_test_feature("mock\tmock\tmock\t100\t200\t.\tmock\t.\tID=f1")
        f2 = create_test_feature("mock\tmock\tmock\t100\t200\t.\tmock\t.\tID=f2")
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 101

