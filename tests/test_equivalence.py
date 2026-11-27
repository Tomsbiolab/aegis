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


class TestReciprocalOverlapEquivalences:
    COLUMNS = ["gene_id_A", "gene_id_B", "gene_id_A_synteny_conserved", "gene_id_B_synteny_conserved", "same_strand", "min_gene_percent", "min_exon_percent", "min_CDS_percent", "gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A_copy", "gene_id_B_copy"]

    def write(self, path, rows):
        lines = ["\t".join(self.COLUMNS)] + ["\t".join(row) for row in rows]
        path.write_text("\n".join(lines) + "\n")

    def test_synteny_of_both_lifted_directions(self, tmp_path):
        from types import SimpleNamespace
        from aegis.equivalence import Simple_annotation

        # Lifted genes are gene_id_A in both files: Q genes lifted onto T, and T genes lifted onto Q
        self.write(tmp_path / "liftoff__Q__to__T_overlaps.tsv", [
            ["q1", "t1", "True", "False", "True", "100.0", "100.0", "100.0", "liftoff_Q", "T", "9", "False", "False"],
        ])
        self.write(tmp_path / "liftoff__T__to__Q_overlaps.tsv", [
            ["t1", "q1", "True", "False", "True", "90.0", "90.0", "90.0", "liftoff_T", "Q", "8", "False", "False"],
        ])
        annot = Simple_annotation("Q", SimpleNamespace(all_gene_ids={"q1": "chr1"}, target=False), "species1")
        annot.add_reciprocal_overlap_equivalences(str(tmp_path), "Q", "T", "species2", quiet=True, synteny_present=True)

        eq = annot.genes["q1"].equivalences[0]
        assert eq.type == "rec_liftoff_aegis"
        assert eq.score == "9 (synteny)/8 (synteny)"
        assert (eq.aegis_synteny, eq.aegis_synteny_rev) == (True, True)

    def test_native_gene_id_B_extra_copy_not_stripped(self, tmp_path):
        from types import SimpleNamespace
        from aegis.equivalence import Simple_annotation

        # Native genes in Q and T happen to have extra_copy_number>=1 (gene_id_B_copy == "True")
        # Only lifted gene_id_A should strip its synthetic _<copy> suffix
        self.write(tmp_path / "liftoff__Q__to__T_overlaps.tsv", [
            ["q_1_1", "t_1", "False", "False", "True", "100.0", "100.0", "100.0", "liftoff_Q", "T", "9", "True", "True"],
        ])
        self.write(tmp_path / "liftoff__T__to__Q_overlaps.tsv", [
            ["t_1_1", "q_1", "False", "False", "True", "95.0", "95.0", "95.0", "liftoff_T", "Q", "8", "True", "True"],
            ["t_2", "q_2", "False", "False", "True", "90.0", "90.0", "90.0", "liftoff_T", "Q", "7", "False", "True"],
        ])
        annot = Simple_annotation("Q", SimpleNamespace(all_gene_ids={"q_1": "chr1", "q_2": "chr1"}, target=False), "species1")
        annot.add_reciprocal_overlap_equivalences(str(tmp_path), "Q", "T", "species2", quiet=True, synteny_present=False)

        eq1 = annot.genes["q_1"].equivalences[0]
        assert eq1.id == "t_1"
        assert eq1.type == "rec_liftoff_aegis"
        assert eq1.score == "9 (copies)/8 (copies)"

        eq2 = annot.genes["q_2"].equivalences[0]
        assert eq2.id == "t_2"
        assert eq2.type == "rev_liftoff_aegis"
        assert eq2.score == "7"


class TestMergeScoreStrings:
    def test_identical_scores_are_kept(self):
        from aegis.cli.orthology import merge_score_strings
        assert merge_score_strings(["rec_liftoff_aegis [9 (synteny)/8]", "rec_liftoff_aegis [9 (synteny)/8]"]) == "rec_liftoff_aegis [9 (synteny)/8]"

    def test_each_direction_keeps_its_own_modifiers(self):
        from aegis.cli.orthology import merge_score_strings
        merged = merge_score_strings(["rec_liftoff_aegis [9 (synteny)/7]", "rec_liftoff_aegis [8/8 (copies) (multiple 2)]"])
        assert merged == "rec_liftoff_aegis [9 (synteny)/8 (copies) (multiple 2)]"

    def test_single_direction_scores(self):
        from aegis.cli.orthology import merge_score_strings
        assert merge_score_strings(["fwd_lifton_aegis [7 (synteny and copies)]", "fwd_lifton_aegis [9]"]) == "fwd_lifton_aegis [9 (synteny and copies)]"
        assert merge_score_strings(["aegis_overlap [6]", "aegis_overlap [8]"]) == "aegis_overlap [8]"


class TestEmptyToolOutputs:
    def test_read_headerless_tsv(self, tmp_path):
        from aegis.utils.misc import read_headerless_tsv
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        df = read_headerless_tsv(empty, ["a", "b"])
        assert df.empty and list(df.columns) == ["a", "b"]
        comments = tmp_path / "comments.anchors"
        comments.write_text("###\n###\n")
        assert read_headerless_tsv(comments, ["0", "1", "2"], comment="#").empty

    def test_blast_and_mcscan_files_without_hits(self, tmp_path):
        from types import SimpleNamespace
        from aegis.equivalence import Simple_annotation
        header = "\tquery_x\tsubject_x\tidentity_x\tcoverage_x\tqlength_x\tslength_x\talength_x\tbitscore_x\tE-value_x\tnorm_bitscore_x\tqcov_x\tscov_x\tquery_y\tsubject_y\tidentity_y\tcoverage_y\tqlength_y\tslength_y\talength_y\tbitscore_y\tE-value_y\tnorm_bitscore_y\tqcov_y\tscov_y\n"
        (tmp_path / "rbbh_Q__to__T.txt").write_text(header)
        (tmp_path / "rbh_Q__to__T.txt").write_text(header)
        (tmp_path / "single_Q__to__T.txt").write_text("")
        (tmp_path / "single_T__to__Q.txt").write_text("")
        (tmp_path / "Q.T.anchors").write_text("###\n")
        annot = Simple_annotation("Q", SimpleNamespace(all_gene_ids={"q1": "chr1"}, target=False), "species1")
        annot.add_blast_equivalences(str(tmp_path), "Q", "T", "species2", quiet=True)
        annot.add_mcscan_equivalences(str(tmp_path / "Q.T.anchors"), "0", "T", "species2")
        assert annot.genes["q1"].equivalences == []


class TestRBBHCoverage:
    def test_rbbh_covering_both_proteins_is_medium(self):
        assert tier(make_equivalence("B", "rbbh"), coverage_threshold=30) == ("medium_confidence", "vvvv_reliable")

    def test_partial_rbbh_is_kept_as_lower(self):
        partial = make_equivalence("B", "rbbh", "evalue=1e-50/1e-40, coverage=90.0/20.0, identity=70.0/60.0")
        assert tier(partial, coverage_threshold=30) == ("lower_confidence", "v_reliable")
        partial = make_equivalence("B", "rbbh", "evalue=1e-50/1e-40, coverage=20.0/90.0, identity=70.0/60.0")
        assert tier(partial, coverage_threshold=30) == ("lower_confidence", "v_reliable")

    def test_identity_does_not_affect_rbbh(self):
        divergent = make_equivalence("B", "rbbh", "evalue=1e-50/1e-40, coverage=90.0/80.0, identity=20.0/20.0")
        assert tier(divergent, coverage_threshold=30, identity_threshold=30) == ("medium_confidence", "vvvv_reliable")

    def test_partial_rbbh_does_not_change_high(self):
        partial = make_equivalence("B", "rbbh", "evalue=1e-50/1e-40, coverage=90.0/20.0, identity=70.0/60.0")
        assert tier(partial, make_equivalence("B", "rec_liftoff_aegis", "7/7"), coverage_threshold=30)[0] == "high_confidence"


class TestEquivalenceRanking:
    def test_na_evalue_does_not_break_tiebreakers(self):
        e_low = Equivalence("g_low", "fwd_liftoff_aegis [7]", "annot2", "sp", "high_confidence", "NA", "top_reliable")
        e_high = Equivalence("g_high", "fwd_liftoff_aegis [8]", "annot2", "sp", "high_confidence", "NA", "top_reliable")
        e_syn = Equivalence("g_syn", "fwd_liftoff_aegis [8 (synteny)]", "annot2", "sp", "high_confidence", "NA", "top_reliable")
        e_finite = Equivalence("g_fin", "fwd_liftoff_aegis [7], rbbh [evalue=1e-10/1e-10, coverage=90/90, identity=90/90]", "annot2", "sp", "high_confidence", "1e-10", "top_reliable")

        assert [e.id for e in sorted([e_low, e_high, e_syn, e_finite])] == ["g_fin", "g_syn", "g_high", "g_low"]
        assert e_syn == Equivalence("g_syn2", "fwd_liftoff_aegis [8 (synteny)]", "annot2", "sp", "high_confidence", "NA", "top_reliable")

    def test_missing_aegis_score_sorts_after_valid_scores(self):
        e_no_ov = Equivalence("g_no_ov", "rbbh [evalue=1e-50/1e-50, coverage=90/90, identity=90/90]", "annot2", "sp", "medium_confidence", "1e-50", "vvvv_reliable")
        e_with_ov = Equivalence("g_with_ov", "fwd_liftoff_aegis [9], rbbh [evalue=1e-50/1e-50, coverage=90/90, identity=90/90]", "annot2", "sp", "medium_confidence", "1e-50", "vvvv_reliable")
        e_rec_ov = Equivalence("g_rec_ov", "rec_liftoff_aegis [9/8], rbbh [evalue=1e-50/1e-50, coverage=90/90, identity=90/90]", "annot2", "sp", "medium_confidence", "1e-50", "vvvv_reliable")

        assert [e.id for e in sorted([e_no_ov, e_with_ov, e_rec_ov])] == ["g_rec_ov", "g_with_ov", "g_no_ov"]


class TestBlastEquivalenceKeys:
    HEADER = "\tquery_x\tsubject_x\tidentity_x\tcoverage_x\tqlength_x\tslength_x\talength_x\tbitscore_x\tE-value_x\tnorm_bitscore_x\tqcov_x\tscov_x\tquery_y\tsubject_y\tidentity_y\tcoverage_y\tqlength_y\tslength_y\talength_y\tbitscore_y\tE-value_y\tnorm_bitscore_y\tqcov_y\tscov_y\n"

    def test_directed_tuple_keys_avoid_shared_and_hyphen_collisions(self, tmp_path):
        from types import SimpleNamespace
        from aegis.equivalence import Simple_annotation

        (tmp_path / "rbbh_Q__to__T.txt").write_text(
            self.HEADER
            + "0\tg1\tg2\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\tg2\tg1\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\n"
            + "1\ta-b\tc\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\tc\ta-b\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\n"
        )
        (tmp_path / "rbh_Q__to__T.txt").write_text(
            self.HEADER
            + "0\tg1\tg2\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\tg2\tg1\t90\t90\t100\t100\t90\t200\t1e-50\t2.0\t0.9\t0.9\n"
            + "1\tg2\tg1\t80\t80\t100\t100\t80\t150\t1e-30\t1.5\t0.8\t0.8\tg1\tg2\t80\t80\t100\t100\t80\t150\t1e-30\t1.5\t0.8\t0.8\n"
            + "2\ta\tb-c\t80\t80\t100\t100\t80\t150\t1e-30\t1.5\t0.8\t0.8\tb-c\ta\t80\t80\t100\t100\t80\t150\t1e-30\t1.5\t0.8\t0.8\n"
        )
        (tmp_path / "single_Q__to__T.txt").write_text("g1\tg2\t90\t90\t100\t100\t90\t200\t1e-50\n")
        (tmp_path / "single_T__to__Q.txt").write_text("g3\tg2\t70\t70\t100\t100\t70\t120\t1e-20\n")

        annot = Simple_annotation("Q", SimpleNamespace(all_gene_ids={"g1": "c1", "g2": "c1", "a-b": "c1", "a": "c1"}, target=False), "sp1")
        annot.add_blast_equivalences(str(tmp_path), "Q", "T", "sp2", quiet=True)

        assert [(e.id, e.type) for e in annot.genes["g1"].equivalences] == [("g2", "rbbh")]
        assert [(e.id, e.type) for e in annot.genes["g2"].equivalences] == [("g1", "rbh"), ("g3", "rev_blastp")]
        assert [(e.id, e.type) for e in annot.genes["a-b"].equivalences] == [("c", "rbbh")]
        assert [(e.id, e.type) for e in annot.genes["a"].equivalences] == [("b-c", "rbh")]


class TestPairwiseOrthologySyntenyStripTags:
    def test_synteny_conserved_when_strip_gene_tags_enabled(self, tmp_path):
        from aegis.equivalence import pairwise_orthology
        from aegis.genome import Genome

        fa_path = tmp_path / "g.fa"
        fa_path.write_text(">chr1\n" + "A" * 5000 + "\n")
        g1_path = tmp_path / "g1.fa"
        g2_path = tmp_path / "g2.fa"
        g1_path.write_text(fa_path.read_text())
        g2_path.write_text(fa_path.read_text())

        gff_a1 = tmp_path / "a1.gff3"
        gff_a1.write_text(
            "##gff-version 3\n"
            "chr1\taegis\tgene\t100\t500\t.\t+\t.\tID=gene-G1\n"
            "chr1\taegis\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=gene-G1\n"
            "chr1\taegis\texon\t100\t500\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr1\taegis\tgene\t1000\t1500\t.\t+\t.\tID=gene-G2\n"
            "chr1\taegis\tmRNA\t1000\t1500\t.\t+\t.\tID=t2;Parent=gene-G2\n"
            "chr1\taegis\texon\t1000\t1500\t.\t+\t.\tID=e2;Parent=t2\n"
        )
        gff_a2 = tmp_path / "a2.gff3"
        gff_a2.write_text(gff_a1.read_text().replace("gene-G", "gene-H"))

        work = tmp_path / "work"
        liftoff_pair = work / "liftoff" / "A1__to__A2"
        liftoff_pair.mkdir(parents=True)
        # Pre-seed liftoff output with stripped IDs (as Liftoff produces from stripped GFF)
        (liftoff_pair / "liftoff__A1__to__A2.gff").write_text(gff_a1.read_text().replace("gene-", ""))

        g1 = Genome(name="G1", genome_file_path=str(g1_path), quiet=True)
        g2 = Genome(name="G2", genome_file_path=str(g2_path), quiet=True)

        pairwise_orthology(
            annot1_name="A1", annot2_name="A2",
            annot1_file=str(gff_a1), annot2_file=str(gff_a2),
            genome1=g1, genome2=g2,
            working_directory=work, num_threads=1, types="",
            synteny=True, strip_gene_tags=True,
            skip_liftoff=False, skip_lifton=True, skip_mcscan=True,
            skip_blasts=True, skip_orthofinder=True, quiet=True,
        )

        out_tsv = work / "liftoff" / "liftoff__A1__to__A2_overlaps.tsv"
        text = out_tsv.read_text()
        assert "G1\tH1\tTrue\t" in text
        assert "G2\tH2\tTrue\t" in text
