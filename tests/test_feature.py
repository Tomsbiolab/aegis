"""
Tests for aegis.feature — the Feature base class.
"""

import pytest

from ..feature import Feature


# ============================================================
# Helpers
# ============================================================

def make_feature(**overrides):
    """Factory for Feature instances with sensible defaults."""
    defaults = dict(
        feature_id="feat001",
        ch="chr1",
        source="aegis",
        feature="gene",
        strand="+",
        start=100,
        end=500,
        score=".",
        phase=".",
        attributes="ID=feat001;Name=TestFeature;Alias=TF1,TF2;Symbol=TFS"
    )
    defaults.update(overrides)
    return Feature(**defaults)


# ============================================================
# __init__ and attribute parsing
# ============================================================

class TestFeatureInit:
    def test_basic_properties(self):
        f = make_feature()
        assert f.id == "feat001"
        assert f.original_id == "feat001"
        assert f.ch == "chr1"
        assert f.source == "aegis"
        assert f.feature == "gene"
        assert f.strand == "+"
        assert f.start == 100
        assert f.end == 500

    def test_size_calculation(self):
        f = make_feature(start=100, end=500)
        assert f.size == 401

    def test_names_parsed(self):
        f = make_feature()
        assert "TestFeature" in f.names

    def test_aliases_parsed(self):
        f = make_feature()
        assert "TF1" in f.aliases
        assert "TF2" in f.aliases

    def test_symbols_parsed(self):
        f = make_feature()
        assert "TFS" in f.symbols

    def test_parent_parsed(self):
        f = make_feature(attributes="ID=exon1;Parent=mRNA1")
        assert "mRNA1" in f.parents

    def test_multiple_parents(self):
        f = make_feature(attributes="ID=exon1;Parent=mRNA1,mRNA2")
        assert "mRNA1" in f.parents
        assert "mRNA2" in f.parents

    def test_id_number_extraction(self):
        f = make_feature(feature_id="gene123")
        assert f.id_number == 123
        assert f.original_id_number == 123

    def test_id_number_no_digits(self):
        f = make_feature(feature_id="gene_abc")
        assert f.id_number is None

    def test_dict_attributes(self):
        attrs = {"id": "feat1", "name": "Test", "parent": ["mRNA1"]}
        f = make_feature(attributes=attrs)
        # feature_id arg is used for self.id, dict is converted to attributes list
        assert f.id == "feat001"
        assert "mRNA1" in f.parents

    def test_list_attributes(self):
        attrs = ["ID=feat1", "Name=Test"]
        f = make_feature(attributes=attrs)
        # feature_id arg is used for self.id, list is stored as self.attributes
        assert f.id == "feat001"
        assert "Test" in f.names

    def test_misc_attributes_collected(self):
        f = make_feature(attributes="ID=feat1;Dbxref=GeneID:12345;custom=value")
        assert any("Dbxref" in a for a in f.misc_attributes)
        assert any("custom" in a for a in f.misc_attributes)


# ============================================================
# Methods
# ============================================================

class TestFeatureMethods:
    def test_update_size(self):
        f = make_feature(start=100, end=500)
        f.end = 600
        f.update_size()
        assert f.size == 501

    def test_print_gff_format(self):
        f = make_feature(attributes="ID=feat001")
        gff_line = f.print_gff()
        assert gff_line.startswith("chr1\taegis\tgene\t100\t500")
        assert "ID=feat001" in gff_line
        assert gff_line.endswith("\n")

    def test_copy_is_independent(self):
        f = make_feature()
        f2 = f.copy()
        f2.id = "changed"
        assert f.id == "feat001"

    def test_str(self):
        f = make_feature()
        assert str(f) == "feat001"

    def test_calculate_gc_content(self):
        f = make_feature(start=1, end=10)
        f.seq = "GGCCAATTGG"  # 6 GC out of 10
        f.calculate_gc_content()
        assert f.gc_content == 0.6

    def test_calculate_gc_content_empty_seq(self):
        f = make_feature()
        f.seq = ""
        f.calculate_gc_content()
        assert f.gc_content == 0

    def test_clear_sequence(self):
        f = make_feature()
        f.seq = "ATGC"
        f.hard_seq = "ATGC"
        f.clear_sequence()
        assert f.seq == ""
        assert f.hard_seq == ""

    def test_clear_sequence_just_hard(self):
        f = make_feature()
        f.seq = "ATGC"
        f.hard_seq = "NNNN"
        f.clear_sequence(just_hard=True)
        assert f.seq == "ATGC"  # preserved
        assert f.hard_seq == ""


# ============================================================
# Comparison operators
# ============================================================

class TestFeatureComparisons:
    def test_equal_features(self):
        f1 = make_feature()
        f2 = make_feature()
        assert f1 == f2

    def test_not_equal_different_id(self):
        f1 = make_feature(feature_id="a")
        f2 = make_feature(feature_id="b")
        assert not (f1 == f2)

    def test_lt_by_start(self):
        f1 = make_feature(start=100, end=500)
        f2 = make_feature(start=200, end=500)
        assert f1 < f2

    def test_lt_same_start_different_end(self):
        f1 = make_feature(start=100, end=400)
        f2 = make_feature(start=100, end=500)
        assert f1 < f2

    def test_le(self):
        f1 = make_feature(start=100, end=500)
        f2 = make_feature(start=100, end=500)
        assert f1 <= f2

    def test_equal_coordinates(self):
        f1 = make_feature(start=100, end=500, ch="chr1")
        f2 = make_feature(start=100, end=500, ch="chr1")
        assert f1.equal_coordinates(f2)

    def test_not_equal_coordinates_different_chr(self):
        f1 = make_feature(ch="chr1")
        f2 = make_feature(ch="chr2")
        assert not f1.equal_coordinates(f2)

    def test_equal_sequence_location(self):
        f1 = make_feature(start=100, end=500, ch="chr1", strand="+")
        f2 = make_feature(start=100, end=500, ch="chr1", strand="+")
        assert f1.equal_sequence(f2)

    def test_not_equal_sequence_different_strand(self):
        f1 = make_feature(strand="+")
        f2 = make_feature(strand="-")
        assert not f1.equal_sequence(f2)

    def test_longer(self):
        f1 = make_feature()
        f1.seq = "ATGCATGC"
        f2 = make_feature()
        f2.seq = "ATGC"
        assert f1.longer(f2) is True
        assert f2.longer(f1) is False
