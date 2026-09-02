import pytest
from typer.testing import CliRunner

from aegis.cli.subset import app
from aegis import Annotation

runner = CliRunner()


def test_cli_subset_no_gene_cap(test_data_dir, tmp_path):
    """Test aegis subset CLI with --no-gene-cap flag."""
    gff3_path = test_data_dir / "input/annotation/arabidopsis_araport11.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "araport11_no_cap.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--no-gene-cap",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"

    out_file = output_dir / output_annot
    assert out_file.exists()

    annot = Annotation(str(out_file), quiet=True)
    assert len(annot.all_gene_ids) == 3000


def test_cli_subset_gene_cap_zero(test_data_dir, tmp_path):
    """Test aegis subset CLI with --gene-cap 0."""
    gff3_path = test_data_dir / "input/annotation/arabidopsis_araport11.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "araport11_zero_cap.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--gene-cap", "0",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"

    out_file = output_dir / output_annot
    assert out_file.exists()

    annot = Annotation(str(out_file), quiet=True)
    assert len(annot.all_gene_ids) == 3000


def test_cli_subset_gene_cap_enforced(test_data_dir, tmp_path):
    """Test aegis subset CLI with an explicit numeric --gene-cap."""
    gff3_path = test_data_dir / "input/annotation/arabidopsis_araport11.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "araport11_200.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--gene-cap", "200",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"

    out_file = output_dir / output_annot
    assert out_file.exists()

    annot = Annotation(str(out_file), quiet=True)
    assert len(annot.all_gene_ids) == 200


def test_cli_subset_no_chr_cap(test_data_dir, tmp_path):
    """Test aegis subset CLI with --no-chr-cap and --no-gene-cap."""
    gff3_path = test_data_dir / "input/annotation/for_merge_2.gff3"
    output_dir = tmp_path / "subset_out"
    output_annot = "for_merge_all_chrs.gff3"

    args = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", output_annot,
        "--no-chr-cap",
        "--no-gene-cap",
        "-q",
    ]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Command failed: {result.stdout}"

    out_file = output_dir / output_annot
    assert out_file.exists()

    annot = Annotation(str(out_file), quiet=True)
    # for_merge_2 has 6 chromosomes: chr1, chr11, chr12, chr19, chr2, chr7
    assert len(annot.chrs) == 6


def test_cli_subset_chr_cap_and_seed(test_data_dir, tmp_path):
    """Test aegis subset CLI with --chr-cap, --no-min-genes, --no-gene-cap, and --seed."""
    gff3_path = test_data_dir / "input/annotation/for_merge_2.gff3"
    output_dir = tmp_path / "subset_out"

    args_1 = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", "run1.gff3",
        "--chr-cap", "2",
        "--no-min-genes",
        "--no-gene-cap",
        "--seed", "99",
        "-q",
    ]
    result_1 = runner.invoke(app, args_1)
    assert result_1.exit_code == 0, f"Command failed: {result_1.stdout}"

    args_2 = [
        str(gff3_path),
        "-d", str(output_dir),
        "-oa", "run2.gff3",
        "--chr-cap", "2",
        "--no-min-genes",
        "--no-gene-cap",
        "--seed", "99",
        "-q",
    ]
    result_2 = runner.invoke(app, args_2)
    assert result_2.exit_code == 0, f"Command failed: {result_2.stdout}"

    annot_1 = Annotation(str(output_dir / "run1.gff3"), quiet=True)
    annot_2 = Annotation(str(output_dir / "run2.gff3"), quiet=True)

    assert len(annot_1.chrs) == 2
    assert set(annot_1.chrs.keys()) == set(annot_2.chrs.keys())
    assert set(annot_1.all_gene_ids) == set(annot_2.all_gene_ids)
