"""
Tests for aegis.genefunctions — pure utility functions.
"""

import pytest

from pathlib import Path
from aegis.utils.genefunctions import (
    reverse_complement,
    find_ORFs,
    longest_ORF,
    trim_surplus,
    translate
)

from aegis.feature import Feature
from aegis.utils.gtf_gff import parse_gff_line, parse_gff_attributes
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
# parse_gff_line
# ============================================================

class TestParseGffLine:
    LINES_FILE = Path(__file__).resolve().parent / "test_data" / "genefunctions_lines.gff3"

    @classmethod
    def _read_line(cls, index):
        """Read a specific line (0-indexed) from genefunctions_lines.gff3."""
        with open(cls.LINES_FILE) as f:
            lines = f.read().splitlines()
        return lines[index]

    def test_basic_parsing(self):
        entry = parse_gff_line(self._read_line(0))  # gene line
        assert entry["ch"] == "chr1"
        assert entry["source"] == "aegis"
        assert entry["feature"] == "gene"
        assert entry["start"] == 1000
        assert entry["end"] == 5000
        assert entry["strand"] == "+"
        assert entry["id"] == "gene1"
        assert entry["pseudogene"] is False
        assert entry["transposable"] is False
        assert entry["decreasing_coordinates"] is False

    def test_pseudogene_detected(self):
        entry = parse_gff_line(self._read_line(1))  # pseudogene line
        assert entry["pseudogene"] is True

    def test_transposable_by_feature(self):
        entry = parse_gff_line(self._read_line(2))  # transposable_element_gene line
        assert entry["transposable"] is True

    def test_transposable_by_attribute(self):
        entry = parse_gff_line(self._read_line(3))  # gene with transposable=True attribute
        assert entry["transposable"] is True

    def test_decreasing_coordinates_swapped(self):
        entry = parse_gff_line(self._read_line(4))  # gene with start > end
        assert entry["decreasing_coordinates"] is True
        assert entry["start"] == 1000
        assert entry["end"] == 5000

    def test_multi_parent(self):
        entry = parse_gff_line(self._read_line(5))  # exon with two parents
        assert entry["parents"] == ["mRNA1", "mRNA2"]



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
    def test_simple_orf(self):
        # ATG + 6 codons + TAA  = 24 nt
        seq = "ATGAAACCCGGGTTTTTTAACTAA"
        orfs = find_ORFs(seq)
        assert len(orfs) >= 1
        # the ORF should start with ATG and end with a stop codon
        orf_seq = orfs[0][0]
        assert orf_seq.startswith("ATG")
        assert orf_seq[-3:] in ("TAA", "TAG", "TGA")

    def test_no_start_codon(self):
        seq = "AAACCCGGGTTTTTTTAA"
        orfs = find_ORFs(seq)
        assert orfs == []

    def test_no_stop_codon_with_must_have_stop(self):
        seq = "ATGAAACCCGGG"  # no stop codon
        orfs = find_ORFs(seq, must_have_stop=True)
        assert orfs == []

    def test_no_stop_codon_without_must_have_stop(self):
        seq = "ATGAAACCCGGG"
        orfs = find_ORFs(seq, must_have_stop=False)
        assert len(orfs) >= 1


class TestLongestORF:
    def test_single_orf(self):
        orfs = [("ATGAAATAA", 0, 8)]
        result = longest_ORF(orfs)
        assert result[0] == "ATGAAATAA"

    def test_multiple_orfs(self):
        orfs = [
            ("ATGTAA", 0, 5),
            ("ATGAAAAAATAA", 10, 21),
            ("ATGCCCTAA", 30, 38),
        ]
        result = longest_ORF(orfs)
        assert result[0] == "ATGAAAAAATAA"

    def test_empty_list(self):
        result = longest_ORF([])
        assert result[0] == ""


# ============================================================
# trim_surplus
# ============================================================

class TestTrimSurplus:
    def test_divisible_by_3(self):
        seq = "ATGAAATAA"  # 9 nt, already divisible
        out, surplus = trim_surplus(seq)
        assert surplus is False
        assert out == seq

    def test_surplus_trimmed(self):
        seq = "ATGAAATAAG"  # 10 nt, surplus = 1
        out, surplus = trim_surplus(seq)
        assert surplus is True
        assert len(out) % 3 == 0


# ============================================================
# translate
# ============================================================

class TestTranslate:
    def test_translate_both_readthrough(self):
        # ATG AAA TAA -> M K *
        seq = "ATGAAATAA"
        start, end_stop, early_stop, surplus, gaps, protein, cs, ce = translate(seq, readthrough="both")
        assert start == "present"
        assert end_stop is True
        assert early_stop is False
        assert "M" in protein
        assert "*" in protein

    def test_translate_no_start(self):
        # GGG AAA TAA -> no M
        seq = "GGGAAATAA"
        start, end_stop, early_stop, surplus, gaps, protein, cs, ce = translate(seq, readthrough="both")
        assert start in ("absent", "late")

    def test_translate_none_readthrough_with_orfs(self):
        # has an ORF inside
        seq = "GGGGATGAAATAAGGG"
        # "none" mode finds the longest ORF
        result = translate(seq, readthrough="none")
        # returns a tuple
        assert isinstance(result, tuple)

    def test_translate_ambiguous_codons(self):
        # N in sequence signals ambiguity
        seq = "ATGNNATAA"
        start, end_stop, early_stop, surplus, gaps, protein, cs, ce = translate(seq, readthrough="both")
        assert gaps is True


# ============================================================
# overlap
# ============================================================

class TestOverlap:
    """Test the overlap function using simple mock objects."""

    class MockFeature(Feature):
        def __init__(self, start, end):
            super().__init__(feature_id="mock", ch="mock", source="mock", feature="mock", strand="mock", start=start, end=end,
                             score=".", phase=".", parents=[], attributes={})

    def test_overlapping_features(self):
        f1 = self.MockFeature(100, 300)
        f2 = self.MockFeature(200, 400)
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 101  # 200..300 inclusive

    def test_non_overlapping_features(self):
        f1 = self.MockFeature(100, 200)
        f2 = self.MockFeature(300, 400)
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is False

    def test_adjacent_features(self):
        f1 = self.MockFeature(100, 200)
        f2 = self.MockFeature(201, 300)
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is False

    def test_contained_feature(self):
        f1 = self.MockFeature(100, 500)
        f2 = self.MockFeature(200, 300)
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True

    def test_identical_features(self):
        f1 = self.MockFeature(100, 200)
        f2 = self.MockFeature(100, 200)
        is_overlapping, bp = f1.overlap(f2)
        assert is_overlapping is True
        assert bp == 101
