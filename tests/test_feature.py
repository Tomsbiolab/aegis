"""
Tests for aegis.feature — the Feature base class.
"""

import pytest

from aegis.feature import Feature


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
        attributes={"Name":["TestFeature"], "Alias": ["TF1", "TF2"], "Symbol": ["TFS"]}
    )
    defaults.update(overrides)
    return Feature(**defaults) # type: ignore


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
        assert "TestFeature" in f.names # type: ignore

    def test_aliases_parsed(self):
        f = make_feature()
        assert "TF1" in f.aliases # type: ignore
        assert "TF2" in f.aliases # type: ignore

    def test_symbols_parsed(self):
        f = make_feature()
        assert "TFS" in f.symbols # type: ignore

    def test_id_number_extraction(self):
        f = make_feature(feature_id="gene123")
        assert f.id_number == 123
        assert f.original_id_number == 123

    def test_id_number_no_digits(self):
        f = make_feature(feature_id="gene_abc")
        assert f.id_number is None

    def test_dict_and_list_attributes(self):
        f = make_feature(feature_id="feat01", attributes={"Name": ["Test"]}, parents=["mRNA1"])
        # feature_id arg is used for self.id, dict is converted to attributes list
        assert f.id == "feat01"
        assert "mRNA1" in f.parents # type: ignore
        assert "Test" in f.names # type: ignore

    def test_misc_attributes_collected(self):
        f = make_feature(feature_id="feat1", attributes={"Dbxref": "GeneID:12345", "custom": "value"})
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
        f = make_feature(feature_id="feat001")
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


# ============================================================
# Comparison operators
# ============================================================

class TestFeatureComparisons:
    def test_equal_features(self):
        f1 = make_feature()
        f2 = make_feature()
        assert f1.identical(f2)

    def test_not_equal_different_id(self):
        f1 = make_feature(feature_id="a")
        f2 = make_feature(feature_id="b")
        assert not (f1.identical(f2))

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
