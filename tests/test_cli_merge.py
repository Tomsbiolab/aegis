from typer.testing import CliRunner

from aegis.cli.merge import app

runner = CliRunner()


def test_merge_quiet_and_no_skip_renaming(tmp_path):
    """aegis merge rejects removed --skip-renaming and stays silent with -q."""
    res_help = runner.invoke(app, ["--help"])
    assert res_help.exit_code == 0
    assert "--skip-renaming" not in res_help.stdout

    g1 = tmp_path / "g1.gff3"
    g1.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    g2 = tmp_path / "g2.gff3"
    g2.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t40\t70\t.\t+\t.\tID=g2\n"
        "chr1\ttest\tmRNA\t40\t70\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\ttest\texon\t40\t70\t.\t+\t.\tID=e2;Parent=t2\n"
    )

    out_dir = tmp_path / "out"
    res = runner.invoke(app, [str(g1), str(g2), "-d", str(out_dir), "-o", "merged.gff3", "-q"])
    assert res.exit_code == 0, res.output
    assert res.output == ""
    merged = (out_dir / "merged.gff3").read_text()
    assert "ID=g1" in merged and "ID=g2" in merged

    res_bad = runner.invoke(app, [str(g1), str(g2), "--skip-renaming", "-d", str(out_dir), "-q"])
    assert res_bad.exit_code != 0
