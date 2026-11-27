"""Tests for aegis overlap CLI and export_group_equivalences."""

import pandas as pd
from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.cli import app
from aegis.genome import Genome
from aegis.utils.genefunctions import export_group_equivalences

runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})


def _write_simple_gff(path, genes):
    lines = ["##gff-version 3"]
    for gid, start, end in genes:
        lines.append(f"chr1\taegis\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}")
        lines.append(f"chr1\taegis\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.1;Parent={gid}")
        lines.append(f"chr1\taegis\texon\t{start}\t{end}\t.\t+\t.\tID={gid}.1.e1;Parent={gid}.1")
    path.write_text("\n".join(lines) + "\n")


class TestExportGroupEquivalences:
    def test_three_annotations_with_non_first_reference(self, tmp_path):
        for name in ("a1", "a2", "a3"):
            _write_simple_gff(tmp_path / f"{name}.gff3", [(f"{name}_g1", 100, 500)])

        annots = [Annotation(str(tmp_path / f"{name}.gff3"), name=name, quiet=True) for name in ("a1", "a2", "a3")]
        annots[1].target = True  # reference is a2 (index 1 != 0)

        out_dir = tmp_path / "out"
        export_group_equivalences(annots, out_dir, quiet=True)

        csv_files = list((out_dir / "overlaps").glob("*.csv"))
        assert len(csv_files) == 1
        df = pd.read_csv(csv_files[0], sep="\t")
        assert set(df["gene_id_A_origin"]) == {"a2"}
        assert set(df["gene_id_B_origin"]) == {"a1", "a3"}

    def test_genome_name_used_in_filename_and_equality(self, tmp_path):
        import pytest

        fa = tmp_path / "g.fa"
        fa.write_text(">chr1\n" + "A" * 1000 + "\n")
        _write_simple_gff(tmp_path / "a1.gff3", [("a1_g1", 100, 500)])
        _write_simple_gff(tmp_path / "a2.gff3", [("a2_g1", 100, 500)])

        g_shared = Genome("MyAssembly", str(fa), quiet=True)
        a1 = Annotation(str(tmp_path / "a1.gff3"), name="a1", genome=g_shared, quiet=True)
        a2 = Annotation(str(tmp_path / "a2.gff3"), name="a2", genome=g_shared, quiet=True)
        export_group_equivalences([a1, a2], tmp_path / "out_shared", quiet=True)

        out_files = [f.name for f in (tmp_path / "out_shared" / "overlaps").iterdir()]
        assert out_files == ["a1_a2_on_MyAssembly_overlaps_t6.csv"]

        # Two separate Genome instances with the same name must not raise ValueError
        g1 = Genome("MyAssembly", str(fa), quiet=True)
        g2 = Genome("MyAssembly", str(fa), quiet=True)
        a1b = Annotation(str(tmp_path / "a1.gff3"), name="a1", genome=g1, quiet=True)
        a2b = Annotation(str(tmp_path / "a2.gff3"), name="a2", genome=g2, quiet=True)
        export_group_equivalences([a1b, a2b], tmp_path / "out_same_name", quiet=True)
        assert [f.name for f in (tmp_path / "out_same_name" / "overlaps").iterdir()] == [
            "a1_a2_on_MyAssembly_overlaps_t6.csv"
        ]

        # Distinct genome names must still raise ValueError
        g_other = Genome("OtherAssembly", str(fa), quiet=True)
        a2c = Annotation(str(tmp_path / "a2.gff3"), name="a2", genome=g_other, quiet=True)
        with pytest.raises(ValueError, match="same genome version"):
            export_group_equivalences([a1b, a2c], tmp_path / "out_diff", quiet=True)


class TestOverlapCLI:
    def test_original_annotation_files_reports_synteny_conserved(self, tmp_path):
        _write_simple_gff(tmp_path / "orig1.gff3", [("g1", 100, 500), ("g2", 1000, 1500)])
        _write_simple_gff(tmp_path / "lifted1.gff3", [("g1", 100, 500), ("g2", 1000, 1500)])
        _write_simple_gff(tmp_path / "target2.gff3", [("h1", 100, 500), ("h2", 1000, 1500)])

        out_dir = tmp_path / "out"
        res = runner.invoke(app, [
            "overlap",
            str(tmp_path / "lifted1.gff3"),
            str(tmp_path / "target2.gff3"),
            "--original-annotation-files",
            f"{tmp_path / 'orig1.gff3'},NA",
            "-d", str(out_dir),
            "-q",
        ])
        assert res.exit_code == 0, res.output
        csv_files = list((out_dir / "overlaps").glob("*.csv"))
        assert len(csv_files) == 1
        df = pd.read_csv(csv_files[0], sep="\t")
        assert list(df["gene_id_A_synteny_conserved"]) == [True, True]
        assert list(df["gene_id_B_synteny_conserved"]) == [False, False]
