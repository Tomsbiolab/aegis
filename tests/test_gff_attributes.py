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
from aegis.utils.gtf_gff import parse_gff_attributes, extra_attributes, SOURCE_KEYS

runner = CliRunner()

RICH = "ID=g1;symbol=S1;gene_synonym=A1;Alias=A2;Dbxref=x:1;Dbxref=y:2;Note=keep%2Cme"
RICH_OUT = RICH  # written back exactly as read, repeated keys included


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
        assert extra_attributes(attrs) == ["Alias=A2", "Dbxref=y:2"]

    def test_parent_wins_over_id_like_fallback(self):
        attrs = parse_gff_attributes("gene_id=g;Parent=t1;ID=e1")
        assert attrs["Parent"] == ["t1"]
        assert extra_attributes(attrs) == ["gene_id=g"]

    def test_id_wins_over_transcript_id(self):
        attrs = parse_gff_attributes("transcript_id=T;ID=t1;Parent=g1", transcript=True)
        assert attrs["ID"] == "t1"
        assert extra_attributes(attrs) == ["transcript_id=T"]

    def test_fallbacks_still_used_when_needed(self):
        attrs = parse_gff_attributes("transcript_id=T;gene_id=G", transcript=True)
        assert attrs["ID"] == "T"
        assert attrs["Parent"] == ["G"]
        assert extra_attributes(attrs) == []


class TestExportAsRead:
    def test_default_export_keeps_attributes(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, RICH, feature="pseudogene")), quiet=True)
        out = _export(annot, tmp_path)
        assert _attributes(out, "g1") == RICH_OUT  # no Pseudogene=True: the feature type says so

        again = _export(Annotation(out, quiet=True), tmp_path, name="again.gff3")
        assert _attributes(again, "g1") == RICH_OUT

    def test_capitalised_symbol_kept(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Symbol=S1;NAME=N1")), quiet=True)
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;Symbol=S1;NAME=N1"

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
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;symbol=S1;Alias=A1;Name=N1"  # keys only
        annot.reorder_attributes()
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


class TestColumnOrder:
    """Name, Alias and symbol are written where they were in column 9, not first."""

    def test_interleaved_order_kept(self, tmp_path):
        attributes = "ID=g1;Dbxref=x:1;Name=N1;gbkey=Gene;Alias=A1;symbol=S1;Note=n"
        annot = Annotation(_write(tmp_path, _gene("g1", 1, attributes)), quiet=True)
        assert _attributes(_export(annot, tmp_path), "g1") == attributes

    def test_default_layout_stores_no_order(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Name=N1;symbol=S1;Note=n")), quiet=True)
        assert annot.chrs["chr1"]["g1"].attribute_order is None

    def test_attributes_the_input_lacked_go_first(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Note=n;Name=N1")), quiet=True)
        annot.set_gene_attribute({"g1": "S1"}, targets="Symbol")
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;symbol=S1;Note=n;Name=N1"

    def test_clean_and_reorder_use_default_order(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Note=n;symbol=S1;Name=N1")), quiet=True)
        assert _attributes(_export(annot, tmp_path, clean_attributes=True), "g1") == "ID=g1;Name=N1;symbol=S1"
        annot.reorder_attributes()
        assert _attributes(_export(annot, tmp_path), "g1") == "ID=g1;Name=N1;symbol=S1;Note=n"

    def test_tidy_reorders_unless_asked_not_to(self, tmp_path):
        path = _write(tmp_path, _gene("g1", 1, "ID=g1;Note=n;symbol=S1;Name=N1"))
        for flag, expected in (([], "ID=g1;Name=N1;symbol=S1;Note=n"), (["--keep-attribute-order"], "ID=g1;Note=n;symbol=S1;Name=N1")):
            out_dir = tmp_path / f"out{len(flag)}"
            res = runner.invoke(tidy_app, [path, "-d", str(out_dir), "-o", "out.gff3", "-q"] + flag)
            assert res.exit_code == 0, res.output
            assert _attributes(str(out_dir / "out.gff3"), "g1") == expected

    def test_rename_puts_alias_where_name_was(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("AT1", 1, "ID=AT1;Dbxref=x;Name=AT1;Note=n")), quiet=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True)
        g = next(iter(annot.chrs["chr1"].values()))
        assert _attributes(_export(annot, tmp_path), g.id) == f"ID={g.id};Dbxref=x;Alias=AT1;Note=n"

    def test_description_goes_where_symbol_was(self, tmp_path):
        annot = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Note=n;symbol=S1")), quiet=True)
        assert _attributes(_export(annot, tmp_path, symbols=False, symbols_as_description=True), "g1") == "ID=g1;Note=n;description=S1"

    def test_order_is_not_part_of_equality(self, tmp_path):
        a = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Note=n;Name=N1"), name="a.gff3"), quiet=True).chrs["chr1"]["g1"]
        b = Annotation(_write(tmp_path, _gene("g1", 1, "ID=g1;Name=N1;Note=n"), name="b.gff3"), quiet=True).chrs["chr1"]["g1"]
        assert a._attributes == b._attributes


class TestSymbolsOrNames:
    """--has-symbol and --gene-symbols also see symbols placed in Name, but not a Name that is the ID."""

    def _annot(self, tmp_path):
        return Annotation(_write(tmp_path,
            _gene("AT1G01010", 1, "ID=AT1G01010;Name=AT1G01010"),
            _gene("AT1G01020", 1001, "ID=AT1G01020;Name=ARV1"),
            _gene("AT1G01030", 2001, "ID=AT1G01030;Name=AT1G01030;symbol=NGA3"),
        ), quiet=True)

    def test_symbols_or_names(self, tmp_path):
        genes = self._annot(tmp_path).chrs["chr1"]
        assert [g.symbols_or_names() for g in genes.values()] == [[], ["ARV1"], ["NGA3"]]

    def test_remove_genes_without_symbols(self, tmp_path):
        annot = self._annot(tmp_path)
        annot.remove_genes_without_symbols(quiet=True)
        assert sorted(annot.all_gene_ids) == ["AT1G01020", "AT1G01030"]

    def test_old_id_kept_as_name_is_not_a_symbol(self, tmp_path):
        annot = self._annot(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            annot.rename_ids(prefix="TST", quiet=True, old_id_names="keep")
        assert [g.symbols_or_names() for g in annot.chrs["chr1"].values()] == [[], ["ARV1"], ["NGA3"]]

    def test_gene_list_symbol_column_is_never_missing(self, tmp_path):
        annot = self._annot(tmp_path)
        annot.export.gene_list(output_dir=str(tmp_path), filename="genes.txt", gene_symbols=True, quiet=True)
        rows = [l.rstrip("\n").split("\t") for l in open(tmp_path / "genes.txt", encoding="utf-8")]
        assert rows[1:] == [["AT1G01010", ""], ["AT1G01020", "ARV1"], ["AT1G01030", "NGA3"]]


class TestTranscripts:
    """Transcripts follow their gene in set_gene_attribute / aegis symbols."""

    def _gff(self, tmp_path):
        g1 = _gene("g1", 1, "ID=g1;Name=g1").replace("ID=g1.t1;Parent=g1", "ID=g1.t1;Parent=g1;Name=g1.1") + \
            "chr1\tsrc\tmRNA\t1\t99\t.\t+\t.\tID=g1.t2;Parent=g1;Name=g1.2\n"
        g2 = _gene("g2", 1001, "ID=g2;Name=Old2").replace("ID=g2.t1;Parent=g2", "ID=g2.t1;Parent=g2;Name=Old2.1")
        return _write(tmp_path, g1, g2)

    def _names(self, annot):
        return {t.id: t.names for g in annot.chrs["chr1"].values() for t in g.transcripts.values()}

    def test_replace_empties_transcripts_of_mapped_genes(self, tmp_path):
        annot = Annotation(self._gff(tmp_path), quiet=True)
        annot.set_gene_attribute({"g1": "PIP1;3"})
        assert self._names(annot) == {"g1.t1": None, "g1.t2": None, "g2.t1": ["Old2.1"]}
        assert _attributes(_export(annot, tmp_path), "g1.t1") == "ID=g1.t1;Parent=g1"

    def test_append_leaves_transcripts(self, tmp_path):
        annot = Annotation(self._gff(tmp_path), quiet=True)
        annot.set_gene_attribute({"g1": "X"}, mode="append")
        assert self._names(annot) == {"g1.t1": ["g1.1"], "g1.t2": ["g1.2"], "g2.t1": ["Old2.1"]}

    def test_clear_unmapped_empties_their_transcripts(self, tmp_path):
        annot = Annotation(self._gff(tmp_path), quiet=True)
        annot.set_gene_attribute({"g1": "X"}, clear_unmapped=True)
        assert self._names(annot)["g2.t1"] is None

    def test_transcript_names_take_gene_values_or_their_own(self, tmp_path):
        annot = Annotation(self._gff(tmp_path), quiet=True)
        assert annot.set_gene_attribute({"g1": "PIP1;3", "g1.t2": "PIP1;3b", "g2.t1": "T"}, transcript_names=True) == []
        assert self._names(annot) == {"g1.t1": ["PIP1;3"], "g1.t2": ["PIP1;3b"], "g2.t1": ["T"]}
        assert annot.chrs["chr1"]["g2"].names == ["Old2"]  # a transcript id does not set its gene
        assert "Name=PIP1%3B3" in _attributes(_export(annot, tmp_path), "g1.t1")

    def test_transcript_ids_need_the_flag(self, tmp_path):
        annot = Annotation(self._gff(tmp_path), quiet=True)
        with pytest.warns(UserWarning, match="1 of them are transcript ids"):
            assert annot.set_gene_attribute({"g1.t1": "X"}) == ["g1.t1"]
        assert self._names(annot)["g1.t1"] == ["g1.1"]

    def test_cli(self, tmp_path):
        tsv = tmp_path / "symbols.tsv"
        tsv.write_text("g1\tNEW\ng1.t2\tNEWb\n", encoding="utf-8")
        for flag, t1, t2 in (([], "ID=g1.t1;Parent=g1", "ID=g1.t2;Parent=g1"), (["--transcript-names"], "ID=g1.t1;Parent=g1;Name=NEW", "ID=g1.t2;Parent=g1;Name=NEWb")):
            out_dir = tmp_path / f"out{len(flag)}"
            res = runner.invoke(symbols_app, [self._gff(tmp_path), str(tsv), "-d", str(out_dir), "-o", "out.gff3", "-q"] + flag)
            assert res.exit_code == 0, res.output
            out = str(out_dir / "out.gff3")
            assert (_attributes(out, "g1.t1"), _attributes(out, "g1.t2"), _attributes(out, "g2.t1")) == (t1, t2, "ID=g2.t1;Parent=g2;Name=Old2.1")


class TestRenameCorrespondences:
    def test_only_with_flag(self, tmp_path):
        path = _write(tmp_path, _gene("g1", 1, "ID=g1"))
        for flag, expected in (([], False), (["--gene-id-correspondences"], True)):
            out_dir = tmp_path / f"out{len(flag)}"
            with warnings.catch_warnings():
                warnings.simplefilter("error", UserWarning)  # no 'Correspondences ... requested' warning
                res = runner.invoke(rename_app, [path, "-p", "TST", "-d", str(out_dir), "-o", "out.gff3", "-q"] + flag)
            assert res.exit_code == 0, res.output
            assert (out_dir / "out_correspondences.tsv").exists() is expected

    def test_no_warning_when_gene_ids_do_not_change(self, tmp_path):
        path = _write(tmp_path, _gene("g1", 1, "ID=g1"))
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            res = runner.invoke(rename_app, [path, "-f", "transcript", "-d", str(tmp_path / "out"), "-o", "out.gff3", "-q"])
        assert res.exit_code == 0, res.output
