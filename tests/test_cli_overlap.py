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
