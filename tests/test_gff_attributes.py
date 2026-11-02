"""
Attributes are kept as read by default (only tidy changes them on request), and gene names, aliases
and symbols can be set from a mapping or a table (aegis symbols).
"""

import warnings

import pytest
from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.cli.filter import app as filter_app
from aegis.cli.prune import app as prune_app
from aegis.cli.rename import app as rename_app
from aegis.cli.reformat import app as reformat_app
from aegis.cli.symbols import app as symbols_app
from aegis.cli.tidy import app as tidy_app
from aegis.utils.gtf_gff import parse_gff_attributes, EXTRA_ATTRIBUTES, SOURCE_KEYS

runner = CliRunner()

RICH = "ID=g1;symbol=S1;gene_synonym=A1;Alias=A2;Dbxref=x:1;Dbxref=y:2;Note=keep%2Cme"
RICH_OUT = "ID=g1;gene_synonym=A1;symbol=S1;Dbxref=x:1;Note=keep%2Cme;Alias=A2;Dbxref=y:2"


def _gene(gene_id, start, gene_attributes, feature="gene"):
    end = start + 98
    return (
        f"chr1\tsrc\t{feature}\t{start}\t{end}\t.\t+\t.\t{gene_attributes}\n"
        f"chr1\tsrc\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}.t1;Parent={gene_id}\n"
        f"chr1\tsrc\texon\t{start}\t{end}\t.\t+\t.\tID={gene_id}.e1;Parent={gene_id}.t1\n"
        f"chr1\tsrc\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene_id}.c1;Parent={gene_id}.t1\n"
    )


def _write(tmp_path, *genes, name="in.gff3"):
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + "".join(genes), encoding="utf-8")
    return str(path)


def _export(annot, tmp_path, name="out.gff3", **kwargs):
    annot.export.gff(output_dir=str(tmp_path), filename=name, quiet=True, **kwargs)
    return str(tmp_path / name)


def _attributes(path, gene_id):
    for line in open(path, encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 9 and parts[8].split(";")[0] == f"ID={gene_id}":
            return parts[8]
    raise AssertionError(f"No {gene_id} line in {path}")


class TestParsePreservesEverything:
    def test_key_spelling_and_repeats(self):
        attrs = parse_gff_attributes(RICH, gene=True)
        assert attrs["Alias"] == ["A1"]
        assert attrs["Symbol"] == ["S1"]
        assert attrs[SOURCE_KEYS] == {"Alias": "gene_synonym"}  # 'symbol' is the default spelling
        assert attrs[EXTRA_ATTRIBUTES] == ["Alias=A2", "Dbxref=y:2"]

    def test_parent_wins_over_id_like_fallback(self):
        attrs = parse_gff_attributes("gene_id=g;Parent=t1;ID=e1")
        assert attrs["Parent"] == ["t1"]
        assert attrs[EXTRA_ATTRIBUTES] == ["gene_id=g"]

    def test_id_wins_over_transcript_id(self):
        attrs = parse_gff_attributes("transcript_id=T;ID=t1;Parent=g1", transcript=True)
        assert attrs["ID"] == "t1"
        assert attrs[EXTRA_ATTRIBUTES] == ["transcript_id=T"]

    def test_fallbacks_still_used_when_needed(self):
        attrs = parse_gff_attributes("transcript_id=T;gene_id=G", transcript=True)
        assert attrs["ID"] == "T"
        assert attrs["Parent"] == ["G"]
        assert EXTRA_ATTRIBUTES not in attrs


class TestExportAsRead:
    def test_default_export_keeps_attributes(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, RICH, feature="pseudogene")), quiet=True)
        out = _export(annot, tmp_path)
        assert _attributes(out, "g1") == RICH_OUT  # no Pseudogene=True: the feature type says so

        again = _export(Annotation(out, quiet=True), tmp_path, name="again.gff3")
        assert _attributes(again, "g1") == RICH_OUT

    def test_capitalised_symbol_kept(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Symbol=S1;NAME=N1")), quiet=True)
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;NAME=N1;Symbol=S1"

    def test_transcript_id_kept(self, tmp_path):
        path = _write(tmp_path, _gene("g1", 1, "ID=g1").replace("ID=g1.t1;Parent=g1", "ID=g1.t1;Parent=g1;transcript_id=T1"))
        assert _attributes(_export(Annotation(path, quiet=True), tmp_path), "g1.t1") == "ID=g1.t1;Parent=g1;transcript_id=T1"

    def test_pseudogene_attribute_kept_as_read(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;pseudogene=True")), quiet=True)
        assert annot.chrs["chr1"]["g1"].pseudogene is True
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;pseudogene=True"

    def test_pseudogene_written_when_otherwise_lost(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;pseudogene=True")), quiet=True)
        assert _attributes(_export(annot, tmp_path, clean_attributes=True), "g1") == "ID=g1;Pseudogene=True"

    def test_clean_attributes_still_available(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, RICH)), quiet=True)
        assert _attributes(_export(annot, tmp_path, clean_attributes=True), "g1") == "ID=g1;gene_synonym=A1;symbol=S1"

    def test_standardise_attribute_keys(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Symbol=S1;gene_synonym=A1;NAME=N1")), quiet=True)
        annot.standardise_attribute_keys()
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;Name=N1;Alias=A1;symbol=S1"

    def test_rename_can_keep_old_names(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("AT4G00005", 1, "ID=AT4G00005;Name=AT4G00005")), quiet=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True, old_id_names="keep")
        g = next(iter(annot.chrs["chr1"].values()))
        assert g.names == ["AT4G00005"] and g.aliases is None


class TestCLIsKeepAttributes:
    """Commands other than tidy only change what they are about."""

    def _two_genes(self, tmp_path):
        return _write(tmp_path, _gene("g1", 1, RICH), _gene("g2", 1001, "ID=g2;Name=N2"))

    @pytest.mark.parametrize("app, extra", [
        (filter_app, ["--coding-only"]),
        (prune_app, ["-i", "g2"]),
        (reformat_app, ["--to", "gff3"]),
        (rename_app, []),
    ])
    def test_untouched_gene_attributes(self, tmp_path, app, extra):
        out_dir = tmp_path / "out"
        res = runner.invoke(app, [self._two_genes(tmp_path), "-d", str(out_dir), "-o", "out.gff3", "-q"] + extra)
        assert res.exit_code == 0, res.output
        assert _attributes(str(out_dir / "out.gff3"), "g1") == RICH_OUT

    def test_rename_moves_old_id_name_unless_kept(self, tmp_path):
        path = _write(tmp_path, _gene("AT4G00005", 1, "ID=AT4G00005;Name=AT4G00005;Note=n"))
        for flag, expected in (([], "Alias=AT4G00005;Note=n"), (["--keep-old-names"], "Name=AT4G00005;Note=n")):
            out_dir = tmp_path / f"out{len(flag)}"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                res = runner.invoke(rename_app, [path, "-p", "TST", "-d", str(out_dir), "-o", "out.gff3", "-q"] + flag)
            assert res.exit_code == 0, res.output
            gene_line = next(l for l in open(out_dir / "out.gff3") if "\tgene\t" in l)
            assert gene_line.rstrip("\n").split("\t")[8].split(";", 1)[1] == expected

    def test_tidy_standardise_attribute_keys(self, tmp_path):
        out_dir = tmp_path / "out"
        path = _write(tmp_path, _gene("g1", 1, "ID=g1;Symbol=S1;gene_synonym=A1"))
        res = runner.invoke(tidy_app, [path, "-d", str(out_dir), "-o", "out.gff3", "-q", "--standardise-attribute-keys"])
        assert res.exit_code == 0, res.output
        assert _attributes(str(out_dir / "out.gff3"), "g1") == "ID=g1;Alias=A1;symbol=S1"


class TestSetGeneAttribute:
    def _annot(self, tmp_path):
        return Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Name=g1;Alias=A1"), _gene("g2", 1001, "ID=g2;Name=Old2")), quiet=True)

    def test_default_replaces_name_of_mapped_genes_only(self, tmp_path):
        annot = self._annot(tmp_path)
        missing = annot.set_gene_attribute({"g1": "PIP1;3"})
        assert missing == []
        assert annot.chrs["chr1"]["g1"].names == ["PIP1;3"]
        assert annot.chrs["chr1"]["g2"].names == ["Old2"]
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;Name=PIP1%3B3;Alias=A1"

    def test_string_is_one_value_and_lists_are_kept_in_order(self, tmp_path):
        annot = self._annot(tmp_path)
        annot.set_gene_attribute({"g1": "x,y", "g2": ["b", "a", "b", "", None]})
        assert annot.chrs["chr1"]["g1"].names == ["x,y"]
        assert annot.chrs["chr1"]["g2"].names == ["b", "a"]

    def test_several_targets_append_and_clear(self, tmp_path):
        annot = self._annot(tmp_path)
        annot.set_gene_attribute({"g1": "S"}, targets=("alias", "SYMBOL"), mode="append", clear_unmapped=True)
        g1, g2 = annot.chrs["chr1"]["g1"], annot.chrs["chr1"]["g2"]
        assert g1.aliases == ["A1", "S"] and g1.symbols == ["S"] and g1.names == ["g1"]
        assert g2.aliases is None and g2.symbols is None and g2.names == ["Old2"]

    def test_missing_ids(self, tmp_path):
        annot = self._annot(tmp_path)
        with pytest.warns(UserWarning, match="1 gene ids not found"):
            assert annot.set_gene_attribute({"nope": "x", "g1": "y"}) == ["nope"]
        assert annot.set_gene_attribute({"nope": "x"}, on_missing="ignore") == ["nope"]
        with pytest.raises(KeyError):
            annot.set_gene_attribute({"nope": "x"}, on_missing="raise")

    def test_original_ids(self, tmp_path):
        annot = self._annot(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True)
        assert annot.set_gene_attribute({"g2": "N"}, id_type="original", on_missing="raise") == []
        assert [g.names for g in annot.chrs["chr1"].values()] == [None, ["N"]]

    def test_bad_target(self, tmp_path):
        with pytest.raises(ValueError):
            self._annot(tmp_path).set_gene_attribute({"g1": "x"}, targets="Description")

    def test_new_values_use_the_file_spelling(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Symbol=S1"), _gene("g2", 1001, "ID=g2")), quiet=True)
        annot.set_gene_attribute({"g2": "S2"}, targets="Symbol")
        assert _attributes(_export(annot, tmp_path), "g2") == "ID=g2;Symbol=S2"

    def test_new_symbols_default_to_lowercase_key(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1")), quiet=True)
        annot.set_gene_attribute({"g1": "S"}, targets="Symbol")
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;symbol=S"

    def test_copy_keeps_set_values(self, tmp_path):
        annot = self._annot(tmp_path)
        annot.set_gene_attribute({"g1": "N"})
        assert annot.copy().chrs["chr1"]["g1"].names == ["N"]


class TestSymbolsCLI:
    def _run(self, tmp_path, table, *flags):
        gff = _write(tmp_path, _gene("g1", 1, "ID=g1;Name=g1;Note=n"), _gene("g2", 1001, "ID=g2;Name=Old2"))
        tsv = tmp_path / "symbols.tsv"
        tsv.write_text(table, encoding="utf-8")
        out_dir = tmp_path / "out"
        res = runner.invoke(symbols_app, [gff, str(tsv), "-d", str(out_dir), "-o", "out.gff3", "-q", *flags])
        return res, out_dir / "out.gff3"

    def test_default_to_name(self, tmp_path):
        res, out = self._run(tmp_path, "g1\tPIP1;3\ng1\tASKβ\nmissing\tX\n")
        assert res.exit_code == 0, res.output
        assert _attributes(str(out), "g1") == "ID=g1;Name=PIP1%3B3,ASKβ;Note=n"
        assert _attributes(str(out), "g2") == "ID=g2;Name=Old2"

    def test_several_targets_and_clear(self, tmp_path):
        res, out = self._run(tmp_path, "g1\tS1\n", "--to", "Name,Symbol", "--clear")
        assert res.exit_code == 0, res.output
        assert _attributes(str(out), "g1") == "ID=g1;Name=S1;symbol=S1;Note=n"
        assert _attributes(str(out), "g2") == "ID=g2"

    def test_append(self, tmp_path):
        res, out = self._run(tmp_path, "g2\tNew2\n", "--append")
        assert res.exit_code == 0, res.output
        assert _attributes(str(out), "g2") == "ID=g2;Name=Old2,New2"

    def test_on_missing_raise(self, tmp_path):
        res, _ = self._run(tmp_path, "missing\tX\n", "--on-missing", "raise")
        assert res.exit_code != 0

    def test_bad_target(self, tmp_path):
        res, _ = self._run(tmp_path, "g1\tX\n", "--to", "Description")
        assert res.exit_code != 0
