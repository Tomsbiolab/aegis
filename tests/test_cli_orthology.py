"""Tests for aegis orthology CLI helpers and edge cases."""

import pandas as pd

from aegis.cli.orthology import get_tiered_cardinality


class TestGetTieredCardinality:
    def test_basic_cardinalities_and_name_order_symmetry(self):
        # 1:1 (x1-y1), 1:N (x2-y2, x2-y3), N:1 (x3-y4, x4-y4), N:N path (x5-y5, x5-y6, x6-y6)
        rows = [
            ("x1", "y1"),
            ("x2", "y2"), ("x2", "y3"),
            ("x3", "y4"), ("x4", "y4"),
            ("x5", "y5"), ("x5", "y6"), ("x6", "y6"),
        ]
        expected = ["1:1", "1:N", "1:N", "N:1", "N:1", "N:N", "N:N", "N:N"]

        for name_a, name_b in (("Alpha", "Beta"), ("Zeta", "Beta")):
            df = pd.DataFrame({
                "annotation_A": [name_a] * len(rows),
                "annotation_B": [name_b] * len(rows),
                "gene_id_A": [a for a, _ in rows],
                "gene_id_B": [b for _, b in rows],
                "summary_score": ["high_confidence"] * len(rows),
            })
            assert list(get_tiered_cardinality(df, ["high_confidence"])) == expected

    def test_reversed_direction_flips_1N_and_N1_symmetrically(self):
        rows = [
            ("x1", "y1"),
            ("x2", "y2"), ("x2", "y3"),
            ("x3", "y4"), ("x4", "y4"),
            ("x5", "y5"), ("x5", "y6"), ("x6", "y6"),
        ]
        df_fwd = pd.DataFrame({
            "annotation_A": ["Alpha"] * len(rows),
            "annotation_B": ["Beta"] * len(rows),
            "gene_id_A": [a for a, _ in rows],
            "gene_id_B": [b for _, b in rows],
            "summary_score": ["high_confidence"] * len(rows),
        })
        df_rev = pd.DataFrame({
            "annotation_A": ["Beta"] * len(rows),
            "annotation_B": ["Alpha"] * len(rows),
            "gene_id_A": [b for _, b in rows],
            "gene_id_B": [a for a, _ in rows],
            "summary_score": ["high_confidence"] * len(rows),
        })
        both = pd.concat([df_fwd, df_rev], ignore_index=True)
        labels = list(get_tiered_cardinality(both, ["high_confidence"]))
        fwd_labels = labels[:len(rows)]
        rev_labels = labels[len(rows):]
        flip = {"1:1": "1:1", "1:N": "N:1", "N:1": "1:N", "N:N": "N:N"}
        assert [flip[x] for x in fwd_labels] == rev_labels


class TestZeroHitOrthologySchema:
    def _setup_non_overlapping_same_genome(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_text(">chr1\n" + ("ATGGCTAGCGATCGATCGAT" * 50) + "TAA" + ("A" * 500) + "\n")
        a1 = tmp_path / "a1.gff3"
        a1.write_text(
            "##gff-version 3\n"
            "chr1\taegis\tgene\t1\t303\t.\t+\t.\tID=g1\n"
            "chr1\taegis\tmRNA\t1\t303\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\taegis\texon\t1\t303\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr1\taegis\tCDS\t1\t303\t.\t+\t0\tID=c1;Parent=t1\n"
        )
        a2 = tmp_path / "a2.gff3"
        a2.write_text(
            "##gff-version 3\n"
            "chr1\taegis\tgene\t601\t903\t.\t+\t.\tID=h1\n"
            "chr1\taegis\tmRNA\t601\t903\t.\t+\t.\tID=u1;Parent=h1\n"
            "chr1\taegis\texon\t601\t903\t.\t+\t.\tID=f1;Parent=u1\n"
            "chr1\taegis\tCDS\t601\t903\t.\t+\t0\tID=d1;Parent=u1\n"
        )
        return a1, a2, fa

    def test_zero_hit_schema_consistency(self, tmp_path):
        from typer.testing import CliRunner
        from aegis.cli import app

        runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})
        a1, a2, fa = self._setup_non_overlapping_same_genome(tmp_path)
        common = [
            "orthology",
            "-a", f"{a1},{a2}",
            "-g", f"{fa},{fa}",
            "-an", "A1,A2",
            "-gn", "G1,G1",
            "--skip-lifton", "--skip-mcscan", "--skip-all-blasts", "--skip-orthofinder",
            "-q",
        ]

        # 1. 0-hit without -na: must include cardinality and drop species_A/species_B
        d1 = tmp_path / "d1"
        res = runner.invoke(app, [*common, "-d", str(d1), "-o", "out.tsv"])
        assert res.exit_code == 0, res.output
        df1 = pd.read_csv(d1 / "out.tsv", sep="\t")
        assert list(df1.columns) == ["gene_id_A", "gene_id_B", "summary_score", "cardinality", "annotation_A", "annotation_B", "score"]
        assert df1.empty

        # 2. 0-hit with -na: must include cardinality="NA" and drop species_A/species_B
        d2 = tmp_path / "d2"
        res = runner.invoke(app, [*common, "-na", "-d", str(d2), "-o", "out.tsv"])
        assert res.exit_code == 0, res.output
        df2 = pd.read_csv(d2 / "out.tsv", sep="\t", keep_default_na=False)
        assert list(df2.columns) == ["gene_id_A", "gene_id_B", "summary_score", "cardinality", "annotation_A", "annotation_B", "score"]
        assert list(df2["cardinality"]) == ["NA", "NA"]

        # 3. 0-hit with -na --tiered-cardinality: must include tiered columns
        d3 = tmp_path / "d3"
        res = runner.invoke(app, [*common, "-na", "--tiered-cardinality", "-d", str(d3), "-o", "out.tsv"])
        assert res.exit_code == 0, res.output
        df3 = pd.read_csv(d3 / "out.tsv", sep="\t", keep_default_na=False)
        assert list(df3.columns) == [
            "gene_id_A", "gene_id_B", "summary_score", "annotation_A", "annotation_B",
            "cardinality_strict", "cardinality_moderate", "cardinality_relaxed", "score",
        ]
        assert list(df3["cardinality_strict"]) == ["NA", "NA"]
        assert list(df3["cardinality_moderate"]) == ["NA", "NA"]
        assert list(df3["cardinality_relaxed"]) == ["NA", "NA"]

        # 4. 0-hit with --group-names: must keep species_A/species_B
        d4 = tmp_path / "d4"
        res = runner.invoke(app, [*common, "--group-names", "sp1,sp2", "-d", str(d4), "-o", "out.tsv"])
        assert res.exit_code == 0, res.output
        df4 = pd.read_csv(d4 / "out.tsv", sep="\t")
        assert "species_A" in df4.columns and "species_B" in df4.columns
