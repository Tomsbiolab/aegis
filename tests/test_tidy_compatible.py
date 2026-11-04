"""
aegis tidy --compatible: a preset of the options that make an annotation readable by most downstream tools
without losing names, attributes or gene models with exons. Options given explicitly win over the preset.
"""

from typer.testing import CliRunner

from aegis.cli.tidy import app as tidy_app

INPUT = (
    "##gff-version 3\n"
    # coding gene; the second CDS row has an inconsistent phase (51 bp first row: recalculated it would be 0)
    "chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1;Symbol=ABC1;Note=kept\n"
    "chr1\tsrc\ttranscript\t100\t900\t.\t+\t.\tID=g1.t1;Parent=g1\n"
    "chr1\tsrc\texon\t100\t150\t.\t+\t.\tID=g1.e1;Parent=g1.t1\n"
    "chr1\tsrc\texon\t301\t900\t.\t+\t.\tID=g1.e2;Parent=g1.t1\n"
    "chr1\tsrc\tCDS\t100\t150\t.\t+\t0\tID=g1.c;Parent=g1.t1\n"
    "chr1\tsrc\tCDS\t301\t602\t.\t+\t2\tID=g1.c;Parent=g1.t1\n"
    "chr1\tsrc\tprotein\t100\t602\t.\t+\t.\tID=g1.p;Derives_from=g1.t1\n"
    # non-coding gene
    "chr1\tsrc\tgene\t1000\t1900\t.\t+\t.\tID=g2;Name=LNC1\n"
    "chr1\tsrc\tlnc_RNA\t1000\t1900\t.\t+\t.\tID=g2.t1;Parent=g2\n"
    "chr1\tsrc\texon\t1000\t1900\t.\t+\t.\tID=g2.e1;Parent=g2.t1\n"
    # gene without transcripts
    "chr1\tsrc\tgene\t2000\t2900\t.\t+\t.\tID=g3\n"
)


def _run(tmp_path, *options):
    path = tmp_path / "in.gff3"
    path.write_text(INPUT, encoding="utf-8")
    res = CliRunner().invoke(tidy_app, [str(path), "-d", str(tmp_path), "-o", "out.gff3", *options])
    assert res.exit_code == 0, res.output
    rows = {}
    for line in open(tmp_path / "out.gff3", encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 9:
            rows.setdefault(parts[8].split(";")[0].removeprefix("ID="), []).append(parts)
    return rows, res.output


def _cds_phases(rows):
    return [r[7] for rs in rows.values() for r in rs if r[2] == "CDS"]


def test_compatible_preset(tmp_path):
    rows, output = _run(tmp_path, "--compatible")
    assert "--compatible enabled: --clean-features, --standard-features, --standardise-attribute-keys, --remove-genes-with-no-transcripts" in output
    # features: standard types, protein row dropped, empty models removed
    assert rows["g1.t1"][0][2] == "mRNA" and rows["g2.t1"][0][2] == "transcript"
    assert "g1.p" not in rows
    assert "g3" not in rows
    # attributes: key standardised, names and other attributes kept
    assert rows["g1"][0][8] == "ID=g1;symbol=ABC1;Note=kept"
    assert rows["g2"][0][8] == "ID=g2;Name=LNC1"
    # annotated phases are kept (only --recalculate-phases overrides them)
    assert _cds_phases(rows) == ["0", "2"]


def test_without_preset_nothing_changes(tmp_path):
    rows, _ = _run(tmp_path, "-q")
    assert rows["g1.t1"][0][2] == "transcript" and rows["g2.t1"][0][2] == "lnc_RNA"
    assert "g1.p" in rows and "g3" in rows
    assert "Symbol=ABC1" in rows["g1"][0][8]
    assert _cds_phases(rows) == ["0", "2"]


def test_explicit_options_win(tmp_path):
    rows, output = _run(tmp_path, "--compatible", "--no-clean-features", "--no-remove-genes-with-no-transcripts", "--no-standardise-attribute-keys")
    assert "--compatible enabled: --standard-features" in output
    assert "g1.p" in rows and "g3" in rows
    assert "Symbol=ABC1" in rows["g1"][0][8]
    assert rows["g1.t1"][0][2] == "mRNA"


def test_flag_and_its_negation_conflict(tmp_path):
    path = tmp_path / "in.gff3"
    path.write_text(INPUT, encoding="utf-8")
    res = CliRunner().invoke(tidy_app, [str(path), "-d", str(tmp_path), "--clean-features", "--no-clean-features"])
    assert res.exit_code != 0


def test_phases_recalculated_only_when_asked(tmp_path):
    rows, _ = _run(tmp_path, "--compatible", "--recalculate-phases", "-q")
    assert _cds_phases(rows) == ["0", "0"]


def test_remove_transcripts_with_no_exons_flag_is_gone(tmp_path):
    path = tmp_path / "in.gff3"
    path.write_text(INPUT, encoding="utf-8")
    res = CliRunner().invoke(tidy_app, [str(path), "-d", str(tmp_path), "--remove-transcripts-with-no-exons"])
    assert res.exit_code != 0
