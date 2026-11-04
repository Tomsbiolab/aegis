"""
Feature standardisation (Annotation(standardise_features=True), aegis tidy --standard-features): genes become
'gene', transcripts 'mRNA' if they have a CDS and 'transcript' otherwise, whatever their type, as tools read
'mRNA' as protein-coding and may only look for the most common types.
"""

from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.cli.tidy import app as tidy_app

MIXED = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1\n"
    "chr1\tsrc\ttranscript\t100\t900\t.\t+\t.\tID=g1.t1;Parent=g1\n"
    "chr1\tsrc\texon\t100\t900\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1\n"
    "chr1\tsrc\tCDS\t200\t799\t.\t+\t0\tID=g1.t1.c;Parent=g1.t1\n"
    "chr1\tsrc\tgene\t1000\t1900\t.\t+\t.\tID=g2\n"
    "chr1\tsrc\ttranscript\t1000\t1900\t.\t+\t.\tID=g2.t1;Parent=g2\n"
    "chr1\tsrc\texon\t1000\t1900\t.\t+\t.\tID=g2.t1.e1;Parent=g2.t1\n"
    "chr1\tsrc\tpseudogene\t2000\t2900\t.\t+\t.\tID=g3\n"
    "chr1\tsrc\tpseudogenic_transcript\t2000\t2900\t.\t+\t.\tID=g3.t1;Parent=g3\n"
    "chr1\tsrc\tpseudogenic_exon\t2000\t2900\t.\t+\t.\tID=g3.t1.e1;Parent=g3.t1\n"
    "chr1\tsrc\ttransposable_element_gene\t3000\t3900\t.\t+\t.\tID=g4\n"
    "chr1\tsrc\tmRNA_TE_gene\t3000\t3900\t.\t+\t.\tID=g4.t1;Parent=g4\n"
    "chr1\tsrc\texon\t3000\t3900\t.\t+\t.\tID=g4.t1.e1;Parent=g4.t1\n"
    "chr1\tsrc\tgene\t4000\t4900\t.\t+\t.\tID=g5\n"
    "chr1\tsrc\tlnc_RNA\t4000\t4900\t.\t+\t.\tID=g5.t1;Parent=g5\n"
    "chr1\tsrc\texon\t4000\t4900\t.\t+\t.\tID=g5.t1.e1;Parent=g5.t1\n"
    "chr1\tsrc\tgene\t5000\t5900\t.\t+\t.\tID=g6\n"
    "chr1\tsrc\tprimary_transcript\t5000\t5900\t.\t+\t.\tID=g6.t1;Parent=g6\n"
    "chr1\tsrc\texon\t5000\t5900\t.\t+\t.\tID=g6.t1.e1;Parent=g6.t1\n"
    "chr1\tsrc\tmiRNA\t5100\t5121\t.\t+\t.\tID=g6.m;Parent=g6.t1\n"
)


def _write(tmp_path, text=MIXED):
    path = tmp_path / "in.gff3"
    path.write_text(text, encoding="utf-8")
    return str(path)


def _types(path):
    """{ID: (feature type, column 9)} of an exported GFF3."""
    out = {}
    for line in open(path, encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 9:
            ID = parts[8].split(";")[0].removeprefix("ID=")
            out.setdefault(ID, (parts[2], parts[8]))
    return out


def test_transcripts_are_mrna_only_with_a_cds(tmp_path):
    annot = Annotation(_write(tmp_path), quiet=True, standardise_features=True)
    transcripts = {t.id: t.feature for genes in annot.chrs.values() for g in genes.values() for t in g.transcripts.values()}
    assert transcripts == {"g1.t1": "mRNA", "g2.t1": "transcript", "g3.t1": "transcript",
                           "g4.t1": "transcript", "g5.t1": "transcript", "g6.t1": "transcript"}
    genes = {g.id: g.feature for genes in annot.chrs.values() for g in genes.values()}
    assert set(genes.values()) == {"gene"}
    assert "standardised_features" in annot.tags


def test_tidy_keeps_gene_status_and_subfeature_types(tmp_path):
    res = CliRunner().invoke(tidy_app, [_write(tmp_path), "-d", str(tmp_path), "-o", "out.gff3", "--standard-features", "-q"])
    assert res.exit_code == 0, res.output
    types = _types(tmp_path / "out.gff3")
    assert types["g3"] == ("gene", "ID=g3;Pseudogene=True")
    assert types["g4"] == ("gene", "ID=g4;Transposable=True")
    assert types["g3.t1.e1"][0] == "exon"
    assert types["g6.m"][0] == "miRNA"
    assert types["g1.t1"][0] == "mRNA" and types["g5.t1"][0] == "transcript"


def test_tidy_filters_biotypes_before_standardising(tmp_path):
    res = CliRunner().invoke(tidy_app, [_write(tmp_path), "-d", str(tmp_path), "-o", "out.gff3", "--standard-features", "-b", "lnc_RNA", "-q"])
    assert res.exit_code == 0, res.output
    types = _types(tmp_path / "out.gff3")
    assert {"g5", "g5.t1"} <= set(types) and types["g5.t1"][0] == "transcript"
    assert not {"g1", "g2", "g3", "g4", "g6"} & set(types)


def test_without_the_option_types_are_kept(tmp_path):
    res = CliRunner().invoke(tidy_app, [_write(tmp_path), "-d", str(tmp_path), "-o", "out.gff3", "-q"])
    assert res.exit_code == 0, res.output
    types = _types(tmp_path / "out.gff3")
    assert (types["g3"][0], types["g3.t1"][0], types["g5.t1"][0], types["g2.t1"][0]) == ("pseudogene", "pseudogenic_transcript", "lnc_RNA", "transcript")
