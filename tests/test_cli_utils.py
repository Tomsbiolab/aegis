import pytest
import typer
from aegis.annotation import Annotation
from aegis.cli.utils import (
    ANNOT_EXTENSIONS,
    FASTA_EXTENSIONS,
    apply_organelle_chroms,
    derive_bio_name,
    resolve_bio_name,
    resolve_unique_names,
    split_callback,
    warn_organelle_genetic_codes,
)


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


def test_derive_and_resolve_bio_name():
    assert derive_bio_name("/a/b/sample.gff3.gz", ANNOT_EXTENSIONS) == "sample"
    assert derive_bio_name("/a/b/release_v1.2.gff3", ANNOT_EXTENSIONS) == "release_v1.2"
    assert derive_bio_name("/a/b/genome.FASTA.GZ", FASTA_EXTENSIONS) == "genome"
    assert derive_bio_name("/a/b/assembly.fna.gz", FASTA_EXTENSIONS) == "assembly"
    assert resolve_bio_name("/a/b/sample.gff3.gz", "{annotation-file}", extensions=ANNOT_EXTENSIONS) == "sample"
    assert resolve_bio_name("/a/b/sample.gff3.gz", "custom_tag", extensions=ANNOT_EXTENSIONS) == "custom_tag"


def test_resolve_unique_names_disambiguation_rules(tmp_path, capsys):
    assert resolve_unique_names(
        ["v1/annot.gff3", "v2/annot.gff3", "v3/annot.gff3"],
        extensions=ANNOT_EXTENSIONS,
        kind="annotation",
        quiet=True,
    ) == ["annot", "annot_2", "annot_3"]

    assert resolve_unique_names(
        ["other.fa", "v1/genome.fa", "v2/genome.fa"],
        extensions=FASTA_EXTENSIONS,
        kind="genome",
        quiet=True,
    ) == ["other", "genome", "genome_2"]

    # Preserves an entry that is genuinely named genome_2 instead of mangling it
    assert resolve_unique_names(
        ["v1/genome.fa", "v2/genome.fa", "genome_2.fa"],
        extensions=FASTA_EXTENSIONS,
        kind="genome",
        quiet=True,
    ) == ["genome", "genome_3", "genome_2"]

    assert resolve_unique_names(
        ["genome_2.fa", "v1/genome.fa", "v2/genome.fa"],
        extensions=FASTA_EXTENSIONS,
        kind="genome",
        quiet=True,
    ) == ["genome_2", "genome", "genome_3"]

    # Same stem across different/gzipped extensions
    assert resolve_unique_names(
        ["annot.gff3", "annot.gtf", "annot.gff3.gz"],
        extensions=ANNOT_EXTENSIONS,
        kind="annotation",
        quiet=True,
    ) == ["annot", "annot_2", "annot_3"]

    # Partial explicit name colliding with derived name
    assert resolve_unique_names(
        ["a.gff3", "b.gff3"],
        provided_names=["b"],
        extensions=ANNOT_EXTENSIONS,
        kind="annotation",
        quiet=True,
    ) == ["b", "b_2"]

    # Emits warning to stderr when quiet=False
    names = resolve_unique_names(
        ["v1/annot.gff3", "v2/annot.gff3"],
        extensions=ANNOT_EXTENSIONS,
        kind="annotation",
        quiet=False,
    )
    assert names == ["annot", "annot_2"]
    err = capsys.readouterr().err
    assert "Duplicate annotation name 'annot' disambiguated to 'annot_2'" in err


def test_resolve_unique_names_validation():
    with pytest.raises(typer.BadParameter, match="does not match"):
        resolve_unique_names(["a.gff3", "b.gff3"], provided_names=["only_one"], strict_count=True, quiet=True)

    with pytest.raises(typer.BadParameter, match="Too many"):
        resolve_unique_names(["a.gff3"], provided_names=["x", "y"], strict_count=False, quiet=True)


def test_resolve_unique_names_share_by_file(tmp_path):
    (tmp_path / "v1").mkdir()
    (tmp_path / "v2").mkdir()
    g1 = tmp_path / "v1" / "genome.fa"
    g2 = tmp_path / "v2" / "genome.fa"
    g1.write_text(">chr1\nACGT\n")
    g2.write_text(">chr1\nTGCA\n")

    # Identical resolved file shares its name; distinct file with same basename is disambiguated
    assert resolve_unique_names(
        [str(g1), str(g1), str(g2)],
        extensions=FASTA_EXTENSIONS,
        kind="genome",
        share_by_file=True,
        quiet=True,
    ) == ["genome", "genome", "genome_2"]

    # Explicit identical names for distinct files are disambiguated
    assert resolve_unique_names(
        [str(g1), str(g2)],
        provided_names=["G", "G"],
        extensions=FASTA_EXTENSIONS,
        kind="genome",
        share_by_file=True,
        quiet=True,
    ) == ["G", "G_2"]

    # Conflicting explicit names for the same file are rejected
    with pytest.raises(typer.BadParameter, match="multiple different tags"):
        resolve_unique_names(
            [str(g1), str(g1)],
            provided_names=["G1", "G2"],
            extensions=FASTA_EXTENSIONS,
            kind="genome",
            share_by_file=True,
            quiet=True,
        )
