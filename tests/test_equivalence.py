"""
Tests for aegis.equivalence — utility functions.
"""

import math
import pytest

from aegis.utils.evalue import parse_evalue, round_evalue


# ============================================================
# round_evalue
# ============================================================

class TestRoundEvalue:
    def test_small_evalue(self):
        result = round_evalue("1e-50")
        assert result == "1.00e-50"

    def test_zero(self):
        result = round_evalue("0")
        assert result == "0.00e+00"

    def test_large_value(self):
        result = round_evalue("10")
        assert result == "1.00e+01"

    def test_decimal(self):
        result = round_evalue("0.005")
        assert result == "5.00e-03"


# ============================================================
# parse_evalue
# ============================================================

class TestParseEvalue:
    def test_normal_evalue(self):
        result = parse_evalue("1e-50")
        assert result == pytest.approx(1e-50)

    def test_greater_than_prefix(self):
        result = parse_evalue(">10")
        assert result == pytest.approx(10.0)

    def test_na_string(self):
        result = parse_evalue("NA")
        assert math.isnan(result)

    def test_nan_string(self):
        result = parse_evalue("nan")
        assert math.isnan(result)

    def test_empty_string(self):
        result = parse_evalue("")
        assert math.isnan(result)

    def test_whitespace_stripped(self):
        result = parse_evalue("  1e-10  ")
        assert result == pytest.approx(1e-10)



# ============================================================
# Confidence tiers of filtered equivalences
# ============================================================

from aegis.equivalence import Equivalence, Simple_gene

BLAST_SCORE = "evalue=1e-50/1e-40, coverage=90.0/80.0, identity=70.0/60.0"

def make_equivalence(target, type_, score=BLAST_SCORE, target_annotation="annot2"):
    evalue = "1e-50" if type_ in ("rbbh", "rbh") or "blast" in type_ else None
    return Equivalence(target, type_, target_annotation, "species2", score, evalue)

def filtered(gene, target="B", **kwargs):
    gene.filter_equivalences(**kwargs)
    hits = [e for e in gene.filtered_equivalences if e.id == target]
    return (hits[0].score, hits[0].reliability) if hits else None

def tier(*equivalences, **kwargs):
    gene = Simple_gene("A", "chr1")
    gene.equivalences = list(equivalences)
    return filtered(gene, **kwargs)


class TestConfidenceTiers:
    def test_high_needs_overlap_of_7(self):
        assert tier(make_equivalence("B", "rec_liftoff_aegis", "7/7")) == ("high_confidence", "top_reliable")
        assert tier(make_equivalence("B", "fwd_lifton_aegis", "7")) == ("high_confidence", "top_reliable")
        assert tier(make_equivalence("B", "aegis_overlap", "9")) == ("high_confidence", "vvtop_reliable")
        assert tier(make_equivalence("B", "rec_liftoff_aegis", "7/7"), make_equivalence("B", "rbbh"), make_equivalence("B", "mcscan_anchors", "100")) == ("high_confidence", "vtop_reliable")

    def test_overlap_of_6_alone_is_lower_for_liftoff_and_lifton(self):
        assert tier(make_equivalence("B", "fwd_liftoff_aegis", "6"))[0] == "lower_confidence"
        assert tier(make_equivalence("B", "fwd_lifton_aegis", "6"))[0] == "lower_confidence"
        assert tier(make_equivalence("B", "aegis_overlap", "6"))[0] == "lower_confidence"

    def test_medium(self):
        assert tier(make_equivalence("B", "rbbh")) == ("medium_confidence", "vvvv_reliable")
        assert tier(make_equivalence("B", "fwd_liftoff_aegis", "6"), make_equivalence("B", "rbh")) == ("medium_confidence", "vvv_reliable")
        assert tier(make_equivalence("B", "rbh"), make_equivalence("B", "mcscan_anchors", "100")) == ("medium_confidence", "vv_reliable")

    def test_lower(self):
        assert tier(make_equivalence("B", "rbh")) == ("lower_confidence", "v_reliable")
        assert tier(make_equivalence("B", "mcscan_anchors", "100")) == ("lower_confidence", "v_reliable")
        assert tier(make_equivalence("B", "fwd_liftoff_aegis", "5"), make_equivalence("B", "orthofinder", "")) == ("lower_confidence", "reliable")


class TestSameAssemblyTiers:
    def test_rbbh_on_different_chromosome_is_lowered(self):
        targets = {"annot2": {"B": Simple_gene("B", "chr2")}}
        assert tier(make_equivalence("B", "rbbh"), same_assembly_targets=targets) == ("lower_confidence", "v_reliable")

    def test_rbbh_on_same_chromosome_stays_medium(self):
        targets = {"annot2": {"B": Simple_gene("B", "chr1")}}
        assert tier(make_equivalence("B", "rbbh"), same_assembly_targets=targets) == ("medium_confidence", "vvvv_reliable")

    def test_rbbh_conflicting_with_an_overlap_is_lowered(self):
        targets = {"annot2": {"B": Simple_gene("B", "chr1"), "C": Simple_gene("C", "chr1")}}
        gene = Simple_gene("A", "chr1")
        gene.equivalences = [make_equivalence("B", "rbbh"), make_equivalence("C", "aegis_overlap", "8")]
        assert filtered(gene, same_assembly_targets=targets) == ("lower_confidence", "v_reliable")
        assert filtered(gene, target="C", same_assembly_targets=targets)[0] == "high_confidence"

    def test_overlap_support_keeps_medium(self):
        targets = {"annot2": {"B": Simple_gene("B", "chr1"), "C": Simple_gene("C", "chr1")}}
        gene = Simple_gene("A", "chr1")
        gene.equivalences = [make_equivalence("B", "rbbh"), make_equivalence("B", "aegis_overlap", "6"), make_equivalence("C", "aegis_overlap", "8")]
        assert filtered(gene, same_assembly_targets=targets) == ("medium_confidence", "vvvv_reliable")

    def test_different_assemblies_are_not_compared(self):
        # Chromosome names of different assemblies are not comparable
        assert tier(make_equivalence("B", "rbbh"), same_assembly_targets={}) == ("medium_confidence", "vvvv_reliable")
        targets = {"annot3": {"B": Simple_gene("B", "chr2")}}
        assert tier(make_equivalence("B", "rbbh"), same_assembly_targets=targets) == ("medium_confidence", "vvvv_reliable")


class TestBlastScoreParsing:
    def test_coverage_and_identity_are_read_from_their_own_fields(self):
        eq = make_equivalence("B", "rbbh")
        assert (eq.coverage, eq.identity) == (90.0, 70.0)
        assert (eq.coverage_rev, eq.identity_rev) == (80.0, 60.0)

    def test_summary_type_parsing(self):
        eq = Equivalence("B", f"rbh [{BLAST_SCORE}]", "annot2", "species2", "medium_confidence")
        assert (eq.coverage, eq.identity) == (90.0, 70.0)
        assert (eq.coverage_rev, eq.identity_rev) == (80.0, 60.0)

    def test_thresholds_apply_to_the_right_value(self):
        low_identity = make_equivalence("B", "rbh", "evalue=1e-50/1e-40, coverage=90.0/80.0, identity=20.0/20.0")
        assert tier(low_identity, identity_threshold=30) is None
        low_identity = make_equivalence("B", "rbh", "evalue=1e-50/1e-40, coverage=90.0/80.0, identity=20.0/20.0")
        assert tier(low_identity, coverage_threshold=30) == ("lower_confidence", "v_reliable")
