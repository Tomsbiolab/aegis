import pytest
import typer
from aegis.annotation import Annotation
from aegis.cli.utils import split_callback, apply_organelle_chroms, warn_organelle_genetic_codes


def test_split_callback_string_comma():
    assert split_callback("gene,CDS,protein") == ["gene", "CDS", "protein"]
    assert split_callback(" gene , CDS , protein ") == ["gene", "CDS", "protein"]


def test_split_callback_single_string():
    assert split_callback("gene") == ["gene"]
    assert split_callback(" gene ") == ["gene"]


def test_split_callback_empty_and_none():
    assert split_callback("") == []
    assert split_callback(None) == []
    assert split_callback([]) == []


def test_split_callback_list_inputs():
    assert split_callback(["gene", "CDS"]) == ["gene", "CDS"]
    assert split_callback(["gene,CDS", "protein"]) == ["gene", "CDS", "protein"]
    assert split_callback([" gene , CDS ", " protein "]) == ["gene", "CDS", "protein"]


def _annotation(tmp_path, name, chroms):
    gff_file = tmp_path / f"{name}.gff3"
    gff_file.write_text("##gff-version 3\n" + "".join(f"{c}\ttest\tgene\t1\t9\t.\t+\t.\tID={name}_{c}\n" for c in chroms))
    return Annotation(str(gff_file), quiet=True)


def test_apply_organelle_chroms_gives_each_annotation_its_listed_chromosomes(tmp_path):
    a1 = _annotation(tmp_path, "a1", ["chr1", "mito_a"])
    a2 = _annotation(tmp_path, "a2", ["chr1", "plastid_b"])
    apply_organelle_chroms([a1, a2], ["mito_a"], ["plastid_b"])
    assert (a1.mitochondria_chroms, a1.chloroplast_chroms) == (("mito_a",), ())
    assert (a2.mitochondria_chroms, a2.chloroplast_chroms) == ((), ("plastid_b",))
    assert a1.chromosome_compartment("mito_a") == "mitochondria"
    assert a2.chromosome_compartment("plastid_b") == "chloroplast"


def test_apply_organelle_chroms_rejects_unknown_chromosomes(tmp_path):
    a1 = _annotation(tmp_path, "a1", ["chr1"])
    with pytest.raises(typer.BadParameter, match="typo_chr"):
        apply_organelle_chroms([a1], ["typo_chr"], [])


def test_warn_organelle_genetic_codes_names_the_tables(tmp_path, capsys):
    annot = _annotation(tmp_path, "a1", ["chr1", "chrM", "ChrC"])
    warn_organelle_genetic_codes(annot)
    err = capsys.readouterr().err
    assert err.startswith("Warning: proteins on mitochondrial contigs (chrM) are translated with NCBI genetic code 1 (Standard)")
    assert "proteins on plastid contigs (ChrC) are translated with NCBI genetic code 11" in err
    assert "taxonomy 'plant'" in err
    warn_organelle_genetic_codes(annot, show_name=True)
    assert capsys.readouterr().err.startswith("Warning: a1: proteins on mitochondrial contigs")


def test_warn_organelle_genetic_codes_follows_the_taxonomy_and_lists(tmp_path, capsys):
    gff_file = tmp_path / "v.gff3"
    gff_file.write_text("##gff-version 3\ncontig_9\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n")
    warn_organelle_genetic_codes(Annotation(str(gff_file), quiet=True, taxonomy="vertebrate", mitochondria_chroms=["contig_9"]))
    assert "proteins on mitochondrial contigs (contig_9) are translated with NCBI genetic code 2 (Vertebrate Mitochondrial)" in capsys.readouterr().err


def test_warn_organelle_genetic_codes_is_silent_without_organelles(tmp_path, capsys):
    warn_organelle_genetic_codes(_annotation(tmp_path, "a1", ["chr1"]))
    assert capsys.readouterr().err == ""
