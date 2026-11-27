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


class TestMcscanGtfBedPreparation:
    def test_gtf_coding_transcripts_standardised_to_mRNA_for_mcscan(self, tmp_path, monkeypatch):
        from pathlib import Path
        from typer.testing import CliRunner
        from aegis.cli import app
        import aegis.cli.orthology as orthology_mod

        fa1 = tmp_path / "g1.fa"
        fa2 = tmp_path / "g2.fa"
        seq = ">chr1\n" + ("ATG" + "GCT" * 50 + "TAA") + ("A" * 200) + "\n"
        fa1.write_text(seq)
        fa2.write_text(seq)

        gtf_template = (
            'chr1\ttest\tgene\t1\t156\t.\t+\t.\tgene_id "{g}";\n'
            'chr1\ttest\ttranscript\t1\t156\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n'
            'chr1\ttest\texon\t1\t156\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n'
            'chr1\ttest\tCDS\t1\t156\t.\t+\t0\tgene_id "{g}"; transcript_id "{g}.1";\n'
        )
        a1 = tmp_path / "a1.gtf"
        a2 = tmp_path / "a2.gtf"
        a1.write_text(gtf_template.format(g="G1"))
        a2.write_text(gtf_template.format(g="H1"))

        mcscan_gffs_seen = []

        def fake_run_command(cwd, cmd):
            if "jcvi.formats.gff" in cmd and "bed" in cmd:
                gff_arg = Path(cmd[cmd.index("--primary_only") + 1])
                mcscan_gffs_seen.append(gff_arg.read_text())
                out_bed = Path(cmd[cmd.index("-o") + 1])
                out_bed.write_text("")
            elif "jcvi.formats.fasta" in cmd:
                Path(cmd[-1]).write_text("")

        monkeypatch.setattr(orthology_mod, "run_command", fake_run_command)
        monkeypatch.setattr(orthology_mod, "pairwise_orthology", lambda **kwargs: None)

        runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})
        out_dir = tmp_path / "out"
        res = runner.invoke(app, [
            "orthology",
            "-a", f"{a1},{a2}",
            "-g", f"{fa1},{fa2}",
            "-an", "A1,A2",
            "-gn", "G1,G2",
            "-d", str(out_dir),
            "--skip-liftoff", "--skip-lifton", "--skip-all-blasts", "--skip-orthofinder",
            "-q",
        ])
        assert res.exit_code == 0, res.output
        assert len(mcscan_gffs_seen) == 2
        for gff_text in mcscan_gffs_seen:
            assert "\tmRNA\t" in gff_text
            assert "\ttranscript\t" not in gff_text


class TestOrthologyNameDisambiguation:
    def _write_annot(self, path, gid, start=1, end=156):
        path.write_text(
            "##gff-version 3\n"
            f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}\n"
            f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.1;Parent={gid}\n"
            f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tID={gid}.1.e1;Parent={gid}.1\n"
            f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID={gid}.1.c1;Parent={gid}.1\n"
        )

    def test_duplicate_annotation_and_distinct_genome_basenames_disambiguated(self, tmp_path, monkeypatch):
        from typer.testing import CliRunner
        from aegis.cli import app
        import aegis.cli.orthology as orthology_mod

        (tmp_path / "sp1").mkdir()
        (tmp_path / "sp2").mkdir()
        seq = ">chr1\n" + ("ATG" + "GCT" * 50 + "TAA") + ("A" * 200) + "\n"
        fa1 = tmp_path / "sp1" / "genome.fa"
        fa2 = tmp_path / "sp2" / "genome.fa"
        fa1.write_text(seq)
        fa2.write_text(seq)

        a1 = tmp_path / "sp1" / "annot.gff3"
        a2 = tmp_path / "sp2" / "annot.gff3"
        self._write_annot(a1, "G1")
        self._write_annot(a2, "H1")

        seen_calls = []

        def fake_pairwise(**kwargs):
            seen_calls.append((
                kwargs["annot1_name"],
                kwargs["annot2_name"],
                kwargs["genome1"].name,
                kwargs["genome2"].name,
                kwargs["pairwise_orthofinder"],
            ))

        monkeypatch.setattr(orthology_mod, "pairwise_orthology", fake_pairwise)

        runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})
        out_dir = tmp_path / "out"
        res = runner.invoke(app, [
            "orthology",
            "-a", f"{a1},{a2}",
            "-g", f"{fa1},{fa2}",
            "-d", str(out_dir),
            "--skip-lifton", "--skip-mcscan", "--skip-all-blasts", "--skip-orthofinder",
            "-na",
            "-q",
        ])
        assert res.exit_code == 0, res.output
        assert seen_calls == [
            ("annot", "annot_2", "genome", "genome_2", False),
            ("annot_2", "annot", "genome_2", "genome", False),
        ]
        df = pd.read_csv(out_dir / "equivalences_with_NAs.tsv", sep="\t", keep_default_na=False)
        assert set(df["annotation_A"]) == {"annot", "annot_2"}

    def test_shared_and_distinct_genomes_with_same_basename(self, tmp_path, monkeypatch):
        """Two annotations sharing sp1/genome.fa keep 'genome', while sp2/genome.fa is disambiguated to 'genome_2'."""
        from typer.testing import CliRunner
        from aegis.cli import app
        import aegis.cli.orthology as orthology_mod

        (tmp_path / "sp1").mkdir()
        (tmp_path / "sp2").mkdir()
        seq = ">chr1\n" + ("ATG" + "GCT" * 50 + "TAA") + ("A" * 200) + "\n"
        fa1 = tmp_path / "sp1" / "genome.fa"
        fa2 = tmp_path / "sp2" / "genome.fa"
        fa1.write_text(seq)
        fa2.write_text(seq)

        a1 = tmp_path / "sp1" / "annot.gff3"
        a2 = tmp_path / "sp1" / "annot_alt.gff3"
        a3 = tmp_path / "sp2" / "annot.gff3"
        self._write_annot(a1, "G1")
        self._write_annot(a2, "G2")
        self._write_annot(a3, "H1")

        seen_genome_names = {}

        def fake_pairwise(**kwargs):
            seen_genome_names[kwargs["annot1_name"]] = kwargs["genome1"].name
            seen_genome_names[kwargs["annot2_name"]] = kwargs["genome2"].name
            assert kwargs["pairwise_orthofinder"] is True

        monkeypatch.setattr(orthology_mod, "pairwise_orthology", fake_pairwise)

        runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})
        out_dir = tmp_path / "out_mix"
        res = runner.invoke(app, [
            "orthology",
            "-a", f"{a1},{a2},{a3}",
            "-g", f"{fa1},{fa1},{fa2}",
            "-d", str(out_dir),
            "--skip-liftoff", "--skip-lifton", "--skip-mcscan", "--skip-all-blasts", "--skip-orthofinder",
        ])
        # All methods skipped -> should fail with No analysis methods selected, after name resolution.
        assert res.exit_code != 0

        res = runner.invoke(app, [
            "orthology",
            "-a", f"{a1},{a2},{a3}",
            "-g", f"{fa1},{fa1},{fa2}",
            "-d", str(out_dir),
            "--skip-lifton", "--skip-mcscan", "--skip-all-blasts", "--skip-orthofinder",
            "-q",
        ])
        assert res.exit_code == 0, res.output
        assert seen_genome_names == {
            "annot": "genome",
            "annot_alt": "genome",
            "annot_2": "genome_2",
        }

    def test_conflicting_explicit_tags_for_same_genome_file_rejected(self, tmp_path):
        from typer.testing import CliRunner
        from aegis.cli import app

        fa = tmp_path / "genome.fa"
        fa.write_text(">chr1\nATG" + "GCT" * 50 + "TAA\n")
        a1 = tmp_path / "a1.gff3"
        a2 = tmp_path / "a2.gff3"
        self._write_annot(a1, "G1")
        self._write_annot(a2, "G2")

        runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})
        res = runner.invoke(app, [
            "orthology",
            "-a", f"{a1},{a2}",
            "-g", f"{fa},{fa}",
            "-gn", "G1,G2",
            "-d", str(tmp_path / "out"),
            "--skip-lifton", "--skip-mcscan", "--skip-all-blasts", "--skip-orthofinder",
            "-q",
        ])
        assert res.exit_code != 0
        normalized = " ".join(res.output.replace("│", " ").split())
        assert "multiple different tags" in normalized
        assert "single consistent tag" in normalized
