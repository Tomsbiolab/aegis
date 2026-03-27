"""
Tests for aegis.feature — the Feature base class.
"""

import pytest

# ============================================================
# __init__ and attribute parsing
# ============================================================

class TestFeatureInit:
    def test_basic_properties(self, make_feature):
        f = make_feature()
        assert f.id == "feat001"
        assert f.original_id == "feat001"
        assert f.ch == "chr1"
        assert f.source == "aegis"
        assert f.feature == "gene"
        assert f.strand == "+"
        assert f.start == 100
        assert f.end == 500

    def test_size_calculation(self, make_feature):
        f = make_feature(start=100, end=500)
        assert f.size == 401

    def test_names_parsed(self, make_feature):
        f = make_feature()
        assert "TestFeature" in f.names

    def test_aliases_parsed(self, make_feature):
        f = make_feature()
        assert "TF1" in f.aliases
        assert "TF2" in f.aliases

    def test_symbols_parsed(self, make_feature):
        f = make_feature()
        assert "TFS" in f.symbols

    def test_id_number_extraction(self, make_feature):
        f = make_feature(feature_id="gene123")
        assert f.id_number == 123
        assert f.original_id_number == 123

    def test_id_number_no_digits(self, make_feature):
        f = make_feature(feature_id="gene_abc")
        assert f.id_number is None

    def test_dict_and_list_attributes(self, make_feature):
        f = make_feature(feature_id="feat01", attributes={"Name": ["Test"]}, parents=["mRNA1"])
        # feature_id arg is used for self.id, dict is converted to attributes list
        assert f.id == "feat01"
        assert "mRNA1" in f.parents # type: ignore
        assert "Test" in f.names # type: ignore

    def test_misc_attributes_collected(self, make_feature):
        f = make_feature(feature_id="feat1", attributes={"Dbxref": "GeneID:12345", "custom": "value"})
        assert any("Dbxref" in a for a in f.misc_attributes)
        assert any("custom" in a for a in f.misc_attributes)


# ============================================================
# Methods
# ============================================================

class TestFeatureMethods:
    def test_update_size(self, make_feature):
        f = make_feature(start=100, end=500)
        f.end = 600
        assert f.size == 501

    def test_print_gff_format(self, make_feature):
        f = make_feature(feature_id="feat001")
        gff_line = f.print_gff()
        assert gff_line.startswith("chr1\taegis\tgene\t100\t500")
        assert "ID=feat001" in gff_line
        assert gff_line.endswith("\n")

    def test_copy_is_independent(self, make_feature):
        f = make_feature()
        f2 = f.copy()
        f2.id = "changed"
        assert f.id == "feat001"

    def test_str(self, make_feature):
        f = make_feature()
        assert str(f) == "feat001"


# ============================================================
# Comparison operators
# ============================================================

class TestFeatureComparisons:
    def test_equal_features(self, make_feature):
        f1 = make_feature()
        f2 = make_feature()
        assert f1.identical(f2)

    def test_not_equal_different_id(self, make_feature):
        f1 = make_feature(feature_id="a")
        f2 = make_feature(feature_id="b")
        assert not (f1.identical(f2))

    def test_lt_by_start(self, make_feature):
        f1 = make_feature(start=100, end=500)
        f2 = make_feature(start=200, end=500)
        assert f1 < f2

    def test_lt_same_start_different_end(self, make_feature):
        f1 = make_feature(start=100, end=400)
        f2 = make_feature(start=100, end=500)
        assert f1 < f2

    def test_le(self, make_feature):
        f1 = make_feature(start=100, end=500)
        f2 = make_feature(start=100, end=500)
        assert f1 <= f2

    def test_equal_coordinates(self, make_feature):
        f1 = make_feature(start=100, end=500, ch="chr1")
        f2 = make_feature(start=100, end=500, ch="chr1")
        assert f1.equal_coordinates(f2)

    def test_not_equal_coordinates_different_chr(self, make_feature):
        f1 = make_feature(ch="chr1")
        f2 = make_feature(ch="chr2")
        assert not f1.equal_coordinates(f2)

    def test_equal_sequence_location(self, make_feature):
        f1 = make_feature(start=100, end=500, ch="chr1", strand="+")
        f2 = make_feature(start=100, end=500, ch="chr1", strand="+")
        assert f1.equal_sequence(f2)

    def test_not_equal_sequence_different_strand(self, make_feature):
        f1 = make_feature(strand="+")
        f2 = make_feature(strand="-")
        assert not f1.equal_sequence(f2)
