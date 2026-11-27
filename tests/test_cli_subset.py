import pytest
from typer.testing import CliRunner

from aegis.cli.subset import app

runner = CliRunner()


def test_cli_subset_gene_cap_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure subset CLI runs with --gene-cap without errors."""
    gff3_path = test_data_dir / "input/annotation/arabidopsis_araport11.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "subset_cap.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--no-chr-cap",
        "--gene-cap", "50",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert (output_dir / output_annot).exists()


def test_cli_subset_chr_cap_smoke(test_data_dir, tmp_path):
    """Smoke test: ensure subset CLI runs with --chr-cap without errors."""
    gff3_path = test_data_dir / "input/annotation/for_merge_2.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "subset_chr.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--chr-cap", "2",
        "--no-min-genes",
        "--no-gene-cap",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert (output_dir / output_annot).exists()


def test_cli_subset_paired_output_file_extensions(tmp_path):
    """Paired -o keeps proper GFF/FASTA extensions and preserves dotted prefixes like release_v1.2."""
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 25 + "\n")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tt\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
    )

    out1 = tmp_path / "out_dot"
    res1 = runner.invoke(app, [str(gff), str(fa), "-c", "chr1", "--no-gene-cap", "-d", str(out1), "-o", "release_v1.2", "-q"])
    assert res1.exit_code == 0, res1.output
    assert sorted(f.name for f in out1.iterdir()) == ["release_v1.2_subset.fasta", "release_v1.2_subset.gff3"]
    assert (out1 / "release_v1.2_subset.fasta").read_text().startswith(">chr1")
    assert (out1 / "release_v1.2_subset.gff3").read_text().startswith("##gff-version 3")

    out2 = tmp_path / "out_gff"
    res2 = runner.invoke(app, [str(gff), str(fa), "-c", "chr1", "--no-gene-cap", "-d", str(out2), "-o", "out.gff3", "-q"])
    assert res2.exit_code == 0, res2.output
    assert sorted(f.name for f in out2.iterdir()) == ["out.gff3", "out_subset.fasta"]
    assert (out2 / "out_subset.fasta").read_text().startswith(">chr1")
    assert (out2 / "out.gff3").read_text().startswith("##gff-version 3")
