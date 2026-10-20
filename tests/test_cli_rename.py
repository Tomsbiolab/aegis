from typer.testing import CliRunner

from aegis.cli.rename import app

runner = CliRunner()


def test_rename_prefix_labels_organelle_genes_smoke(tmp_path):
    """Smoke test: genes on detected or listed organelle contigs are renamed as chrMT/chrCP."""
    gff_file = tmp_path / "o.gff3"
    gff_file.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chrM\ttest\tgene\t1\t9\t.\t+\t.\tID=g2\n"
        "contig_9\ttest\tgene\t1\t9\t.\t+\t.\tID=g3\n"
        "ChrC\ttest\tgene\t1\t9\t.\t+\t.\tID=g4\n"
    )
    out_dir = tmp_path / "rename_out"
    result = runner.invoke(app, [str(gff_file), "-p", "TST", "--spacer", "10", "--mitochondria-chroms", "contig_9", "-d", str(out_dir), "-o", "out.gff3", "-q"])
    assert result.exit_code == 0, f"Error: {result.stdout}"
    ids = {line.split("\t")[0]: line.split("\t")[8].strip() for line in (out_dir / "out.gff3").read_text().splitlines() if line and not line.startswith("#")}
    assert ids["chr1"] == "ID=TSTchr1g00010"
    assert ids["ChrC"] == "ID=TSTchrCPg00010"
    assert sorted([ids["chrM"], ids["contig_9"]]) == ["ID=TSTchrMTg00010", "ID=TSTchrMTg00020"]
