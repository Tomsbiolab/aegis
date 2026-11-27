from typer.testing import CliRunner

from aegis.cli.prune import app

runner = CliRunner()


def test_prune_keep_gene_and_transcript(tmp_path):
    """--keep retains only specified gene or transcript IDs without crashing."""
    gff = tmp_path / "test.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\ttest\tgene\t200\t300\t.\t+\t.\tID=g2\n"
        "chr1\ttest\tmRNA\t200\t300\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\ttest\texon\t200\t300\t.\t+\t.\tID=e2;Parent=t2\n"
        "chr1\ttest\tmRNA\t210\t290\t.\t+\t.\tID=t3;Parent=g2\n"
        "chr1\ttest\texon\t210\t290\t.\t+\t.\tID=e3;Parent=t3\n"
    )

    out_gene = tmp_path / "keep_gene"
    res_gene = runner.invoke(app, [str(gff), "g1", "--keep", "-d", str(out_gene), "-o", "out.gff3", "-q"])
    assert res_gene.exit_code == 0, res_gene.output
    text_gene = (out_gene / "out.gff3").read_text()
    assert "ID=g1" in text_gene
    assert "ID=g2" not in text_gene

    out_tx = tmp_path / "keep_tx"
    res_tx = runner.invoke(app, [str(gff), "t2", "-f", "transcript", "--keep", "-d", str(out_tx), "-o", "out.gff3", "-q"])
    assert res_tx.exit_code == 0, res_tx.output
    text_tx = (out_tx / "out.gff3").read_text()
    assert "ID=t2" in text_tx
    assert "ID=g2" in text_tx
    assert "ID=t3" not in text_tx
    assert "ID=g1" not in text_tx
