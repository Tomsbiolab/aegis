from typer.testing import CliRunner

from aegis.cli.tidy_genome import app

runner = CliRunner()


def test_tidy_genome_paired_output_file_extensions(tmp_path):
    """Paired -o in tidy-genome preserves dotted prefixes and never swaps FASTA/GFF extensions."""
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 25 + "\n")
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tt\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
    )

    out1 = tmp_path / "out_dot"
    res1 = runner.invoke(app, [str(fa), str(gff), "-d", str(out1), "-o", "release_v1.2", "-q"])
    assert res1.exit_code == 0, res1.output
    assert sorted(f.name for f in out1.iterdir()) == [
        "release_v1.2_matching_genome_tidy.gff3",
        "release_v1.2_tidy.fasta",
    ]
    assert (out1 / "release_v1.2_tidy.fasta").read_text().startswith(">chr1")
    assert (out1 / "release_v1.2_matching_genome_tidy.gff3").read_text().startswith("##gff-version 3")

    out2 = tmp_path / "out_gff"
    res2 = runner.invoke(app, [str(fa), str(gff), "-d", str(out2), "-o", "out.gff3", "-q"])
    assert res2.exit_code == 0, res2.output
    assert sorted(f.name for f in out2.iterdir()) == ["out.gff3", "out_tidy.fasta"]
    assert (out2 / "out_tidy.fasta").read_text().startswith(">chr1")
    assert (out2 / "out.gff3").read_text().startswith("##gff-version 3")
