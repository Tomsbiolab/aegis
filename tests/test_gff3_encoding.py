"""
GFF3 column 9 percent-encoding: values are decoded on read, stored decoded and encoded on export.
"""

import warnings

import pytest
from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.cli.tidy import app as tidy_app
from aegis.utils.gtf_gff import gff3_escape, gff3_unescape, format_gff3_values, parse_gff_attributes, format_gff3_attributes

TRICKY = ["PIP1;3", "a=b", "x,y", "50%", "ASKβ"]


def _gff(tmp_path, gene_attributes, mrna_attributes="ID=t1;Parent=g1", name="in.gff3"):
    path = tmp_path / name
    path.write_text(
        "##gff-version 3\n"
        f"chr1\tsrc\tgene\t1\t99\t.\t+\t.\t{gene_attributes}\n"
        f"chr1\tsrc\tmRNA\t1\t99\t.\t+\t.\t{mrna_attributes}\n"
        "chr1\tsrc\texon\t1\t99\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tsrc\tCDS\t1\t99\t.\t+\t0\tID=c1;Parent=t1\n",
        encoding="utf-8",
    )
    return str(path)


def _export(annot, tmp_path, name="out.gff3", **kwargs):
    annot.export.gff(output_dir=str(tmp_path), filename=name, quiet=True, **kwargs)
    return str(tmp_path / name)


def _line(path, feature):
    for line in open(path, encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 9 and parts[2] == feature:
            return parts[8]
    raise AssertionError(f"No {feature} line in {path}")


class TestEscape:
    def test_reserved_characters(self):
        assert gff3_escape("%;=&,") == "%25%3B%3D%26%2C"

    def test_control_characters(self):
        assert gff3_escape("a\tb\nc\rd\x00e\x1ff\x7f") == "a%09b%0Ac%0Dd%00e%1Ff%7F"

    def test_percent_is_encoded_once(self):
        assert gff3_escape("50%3B") == "50%253B"

    def test_left_as_is(self):
        value = "a b/c'd+e(f)β:g|h"
        assert gff3_escape(value) == value

    def test_unescape_leaves_plus_and_lone_percent(self):
        assert gff3_unescape("a+b") == "a+b"
        assert gff3_unescape("50%") == "50%"
        assert gff3_unescape("100%zz") == "100%zz"

    @pytest.mark.parametrize("value", TRICKY + ["\t", "%25", "a%3Bb", "+ %"])
    def test_round_trip(self, value):
        assert gff3_unescape(gff3_escape(value)) == value

    def test_multi_values_encoded_before_joining(self):
        assert format_gff3_values(["x,y", "z"]) == "x%2Cy,z"


class TestParse:
    def test_mapped_values_decoded_after_split(self):
        attrs = parse_gff_attributes("ID=g1;Name=PIP1%3B3,a%3Db;Symbol=x%2Cy;Alias=50%25,ASK%CE%B2")
        assert attrs["Name"] == ["PIP1;3", "a=b"]
        assert attrs["Symbol"] == ["x,y"]
        assert attrs["Alias"] == ["50%", "ASKβ"]

    def test_ids_parents_and_misc_kept_as_read(self):
        attrs = parse_gff_attributes("ID=g%2C1;Parent=p%3B1,p2;Note=keep%2Cme", gene=True)
        assert attrs["ID"] == "g%2C1"
        assert attrs["Parent"] == ["p%3B1", "p2"]
        assert attrs["Note"] == "keep%2Cme"

    def test_lowercase_keys_decoded(self):
        attrs = parse_gff_attributes("ID=g1;symbol=a%3Bb;gene_synonym=c%2Cd")
        assert attrs["Symbol"] == ["a;b"]
        assert attrs["Alias"] == ["c,d"]


class TestExportRoundTrip:
    def test_values_set_in_memory_round_trip(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1"), quiet=True)
        g = annot.chrs["chr1"]["g1"]
        g.names, g.symbols, g.aliases = TRICKY[:], TRICKY[:], TRICKY[:]

        out = _export(annot, tmp_path)
        assert _line(out, "gene") == (
            "ID=g1;Name=PIP1%3B3,a%3Db,x%2Cy,50%25,ASKβ;Alias=PIP1%3B3,a%3Db,x%2Cy,50%25,ASKβ;"
            "symbol=PIP1%3B3,a%3Db,x%2Cy,50%25,ASKβ"
        )

        g2 = Annotation(out, quiet=True).chrs["chr1"]["g1"]
        assert g2.names == TRICKY
        assert g2.symbols == TRICKY
        assert g2.aliases == TRICKY

    def test_file_read_export_read(self, tmp_path):
        attributes = "ID=g1;Name=PIP1%3B3,ASKβ;Symbol=a%3Db;Alias=x%2Cy,50%25;Note=keep%2Cme"
        annot = Annotation(_gff(tmp_path, attributes), quiet=True)
        g = annot.chrs["chr1"]["g1"]
        assert g.names == ["PIP1;3", "ASKβ"]
        assert g.symbols == ["a=b"]
        assert g.aliases == ["x,y", "50%"]

        out = _export(annot, tmp_path, clean_attributes=False)
        line = _line(out, "gene")
        assert "Name=PIP1%3B3,ASKβ" in line
        assert "Symbol=a%3Db" in line
        assert "Alias=x%2Cy,50%25" in line
        assert "Note=keep%2Cme" in line  # misc attributes are not re-encoded

        g2 = Annotation(out, quiet=True).chrs["chr1"]["g1"]
        assert (g2.names, g2.symbols, g2.aliases) == (g.names, g.symbols, g.aliases)

    def test_description_from_symbols_encoded(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1;Symbol=a%3Bb"), quiet=True)
        out = _export(annot, tmp_path, symbols=False, symbols_as_description=True)
        assert _line(out, "gene") == "ID=g1;Description=a%3Bb"

    def test_transcript_names_encoded(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1", mrna_attributes="ID=t1;Parent=g1;Name=t%3B1"), quiet=True)
        assert annot.chrs["chr1"]["g1"].transcripts["t1"].names == ["t;1"]
        assert "Name=t%3B1" in _line(_export(annot, tmp_path), "mRNA")


class TestNameExport:
    def test_identity_attributes_written_by_default(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1;Name=N1;Symbol=S1;Alias=A1"), quiet=True)
        assert _line(_export(annot, tmp_path), "gene") == "ID=g1;Name=N1;Alias=A1;Symbol=S1"

    def test_attributes_can_be_left_out(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1;Name=N1;Symbol=S1;Alias=A1"), quiet=True)
        out = _export(annot, tmp_path, names=False, symbols=False, aliases=False)
        assert _line(out, "gene") == "ID=g1"

    def test_copy_keeps_attributes(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=g1;Name=PIP1%3B3;Symbol=S1;Alias=A1"), quiet=True)
        copied = annot.copy()
        g, gc = annot.chrs["chr1"]["g1"], copied.chrs["chr1"]["g1"]
        assert (gc.names, gc.symbols, gc.aliases) == (["PIP1;3"], ["S1"], ["A1"])
        assert gc.names is not g.names
        assert _line(_export(copied, tmp_path), "gene") == "ID=g1;Name=PIP1%3B3;Alias=A1;Symbol=S1"

    def test_tidy_remove_names(self, tmp_path):
        gff = _gff(tmp_path, "ID=g1;Name=N1;Symbol=S1")
        out_dir = tmp_path / "out"
        res = CliRunner().invoke(tidy_app, [gff, "-d", str(out_dir), "-o", "tidy.gff3", "--remove-names", "-q"])
        assert res.exit_code == 0, res.output
        line = _line(str(out_dir / "tidy.gff3"), "gene")
        assert "Name=" not in line
        assert "Symbol=S1" in line


class TestRenameMovesStaleNames:
    def test_name_equal_to_old_id_becomes_alias(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=AT4G00005;Name=AT4G00005", mrna_attributes="ID=t1;Parent=AT4G00005;Name=t1"), quiet=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True)
        g = next(iter(annot.chrs["chr1"].values()))
        assert g.id != "AT4G00005"
        assert g.names is None
        assert g.aliases == ["AT4G00005"]
        t = next(iter(g.transcripts.values()))
        assert t.names is None and t.aliases == ["t1"]

        line = _line(_export(annot, tmp_path), "gene")
        assert "Name=" not in line
        assert "Alias=AT4G00005" in line

    def test_real_names_kept(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=AT4G00005;Name=NAC001,AT4G00005;Alias=X"), quiet=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True)
        g = next(iter(annot.chrs["chr1"].values()))
        assert g.names == ["NAC001"]
        assert g.aliases == ["X", "AT4G00005"]

    def test_unrenamed_names_untouched(self, tmp_path):
        annot = Annotation(_gff(tmp_path, "ID=AT4G00005;Name=AT4G00005"), quiet=True)
        g = annot.chrs["chr1"]["AT4G00005"]
        assert g.names == ["AT4G00005"]
        assert g.aliases is None


def test_gtf_conversion_encodes_values():
    result = format_gff3_attributes({"gene_id": "g1", "gene_name": "PIP1;3", "note": "a=b,c"}, "gene")
    assert result == "ID=g1;Symbol=PIP1%3B3;note=a%3Db%2Cc"
