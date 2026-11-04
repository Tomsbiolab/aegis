"""
Load-time handling and warnings for cases from AGAT's gff_syntax suite and Araport11:
- AGAT 15: a file of match/match_part lines only reads no gene model, which is now said;
- AGAT 6: mRNAs linked by locus_tag join the childless gene with their span that lacks the locus_tag;
- AGAT 36: a transcript trimmed to its exons (exon 1 missing) gets a message pointing at the missing exon;
- Araport11 AT4G15260.1: UTR rows replaced by derived UTRs are counted when loading.
"""

import pytest

from aegis.annotation import Annotation

HEADER = "##gff-version 3\n"


def write(tmp_path, rows, name="in.gff3"):
    f = tmp_path / name
    f.write_text(HEADER + "\n".join("\t".join(r.split("|")) for r in rows) + "\n")
    return str(f)


def genes(annot):
    return {g.id: sorted(g.transcripts) for gs in annot.chrs.values() for g in gs.values()}


class TestNoGeneModels:

    def test_match_only_file_warns(self, tmp_path, capsys):
        f = write(tmp_path, [
            "s1|maker|match|100|900|.|+|.|ID=m1;Parent=G1",
            "s1|maker|match_part|100|300|.|+|.|ID=m1p1;Parent=m1;Target=m1 1 201",
            "s1|maker|match_part|500|900|.|+|.|ID=m1p2;Parent=m1;Target=m1 202 602",
        ])
        annot = Annotation(f)
        out = capsys.readouterr().out
        assert genes(annot) == {} and len(annot.atypical_features) == 3
        assert "no gene models were read" in out and "only its 3 lines of other feature types (match_part: 2, match: 1) are kept" in out
        assert "match/match_part lines are alignments" in out

        annot = Annotation(f, skip_atypical_features=True)
        out = capsys.readouterr().out
        assert annot.atypical_features == []
        assert "3 lines of feature types that are not genes, transcripts or their subfeatures were not read (match_part: 2, match: 1)" in out
        assert "nothing but its header will be exported" in out and "match/match_part lines are alignments" in out

    def test_skipped_types_counted_only_on_chosen_chromosomes(self, tmp_path, capsys):
        f = write(tmp_path, [
            "s1|src|gene|100|900|.|+|.|ID=g1",
            "s1|src|mRNA|100|900|.|+|.|ID=t1;Parent=g1",
            "s1|src|exon|100|900|.|+|.|ID=e1;Parent=t1",
            "s1|src|repeat_region|1000|1200|.|+|.|ID=r1",
            "s2|src|repeat_region|1000|1200|.|+|.|ID=r2",
        ])
        Annotation(f, chosen_chromosomes=("s1",), skip_atypical_features=True)
        out = capsys.readouterr().out
        assert "1 lines of feature types that are not genes, transcripts or their subfeatures were not read (repeat_region: 1)" in out
        assert "no gene models were read" not in out


class TestChildlessGeneTakesLinkedFeatures:

    ROWS = [
        "s789|maker|gene|558184|564780|.|+|.|ID=CLUHARG00000003852;Name=PF11_0240",
        "s789|maker|mRNA|558184|564780|.|+|.|ID=T6146;locus_tag=AAAAA",
        "s789|maker|exon|558184|560123|.|+|.|ID=T6146:exon:1;locus_tag=AAAAA",
        "s789|maker|exon|564372|564780|.|+|.|ID=T6146:exon:2;locus_tag=AAAAA",
        "s789|maker|mRNA|558184|564780|.|+|.|ID=T6147;locus_tag=AAAAA",
        "s789|maker|exon|558184|560123|.|+|.|ID=T6147:exon:1;locus_tag=AAAAA",
        "s789|maker|exon|562057|564780|.|+|.|ID=T6147:exon:2;locus_tag=AAAAA",
    ]

    def test_agat_6_mrnas_join_the_gene_line(self, tmp_path):
        annot = Annotation(write(tmp_path, self.ROWS), quiet=True)
        assert genes(annot) == {"CLUHARG00000003852": ["T6146", "T6147"]}
        assert annot.warnings["linked_to_childless_gene_by_span"] == {"CLUHARG00000003852"}
        assert annot.chrs["s789"]["CLUHARG00000003852"].names == ["PF11_0240"]

    def test_other_span_or_strand_still_creates_a_gene(self, tmp_path):
        rows = [self.ROWS[0].replace("558184|564780|.|+", "558000|564780|.|+")] + self.ROWS[1:]
        assert set(genes(Annotation(write(tmp_path, rows), quiet=True))) == {"CLUHARG00000003852", "AAAAA"}
        rows = [self.ROWS[0].replace("|+|", "|-|")] + self.ROWS[1:]
        assert set(genes(Annotation(write(tmp_path, rows, name="minus.gff3"), quiet=True))) == {"CLUHARG00000003852", "AAAAA"}

    def test_gene_with_its_own_locus_tag_or_children_is_not_taken(self, tmp_path):
        rows = [self.ROWS[0] + ";locus_tag=BBBBB"] + self.ROWS[1:]
        assert set(genes(Annotation(write(tmp_path, rows), quiet=True))) == {"CLUHARG00000003852", "AAAAA"}
        rows = self.ROWS + [
            "s789|maker|mRNA|558184|564780|.|+|.|ID=T1;Parent=CLUHARG00000003852",
            "s789|maker|exon|558184|564780|.|+|.|ID=T1:exon:1;Parent=T1",
        ]
        assert genes(Annotation(write(tmp_path, rows, name="children.gff3"), quiet=True)) == {"CLUHARG00000003852": ["T1"], "AAAAA": ["T6146", "T6147"]}

    def test_two_childless_genes_with_the_span_are_ambiguous(self, tmp_path):
        rows = [self.ROWS[0], self.ROWS[0].replace("CLUHARG00000003852", "OTHER")] + self.ROWS[1:]
        annot = Annotation(write(tmp_path, rows), quiet=True)
        assert set(genes(annot)) == {"CLUHARG00000003852", "OTHER", "AAAAA"}
        assert annot.warnings["linked_to_childless_gene_by_span"] == set()


class TestTrimmedTranscripts:

    def test_transcript_with_missing_first_exon_points_at_it(self, tmp_path, capsys):
        f = write(tmp_path, [
            "c1|HAVANA|gene|540236|543139|.|+|.|ID=g1",
            "c1|HAVANA|mRNA|540236|543139|.|+|.|ID=t1;Parent=g1",
            "c1|HAVANA|exon|542087|543139|.|+|.|ID=t1:exon:2;Parent=t1",
        ])
        annot = Annotation(f)
        t = annot.chrs["c1"]["g1"].transcripts["t1"]
        assert (t.start, t.end) == (542087, 543139)
        assert "1 transcripts extended past their first or last exon and were trimmed to them; an exon line may be missing from the file" in capsys.readouterr().out


class TestReadUTRsReplaced:

    def test_utr_row_past_its_exon_is_counted(self, tmp_path, capsys):
        f = write(tmp_path, [
            "Chr4|A|gene|8713689|8715380|.|+|.|ID=AT4G15260",
            "Chr4|A|mRNA|8713689|8715380|.|+|.|ID=AT4G15260.1;Parent=AT4G15260",
            "Chr4|A|exon|8713689|8715339|.|+|.|ID=e1;Parent=AT4G15260.1",
            "Chr4|A|CDS|8714065|8715144|.|+|0|ID=c1;Parent=AT4G15260.1",
            "Chr4|A|five_prime_UTR|8713689|8714064|.|+|.|ID=u5;Parent=AT4G15260.1",
            "Chr4|A|three_prime_UTR|8715145|8715380|.|+|.|ID=u3;Parent=AT4G15260.1",
        ])
        annot = Annotation(f)
        assert annot.warnings["read_utrs_differ_from_derived"] == {"AT4G15260.1"}
        assert "1 transcripts have UTR rows that differ from the UTRs derived from their exons and CDS" in capsys.readouterr().out


# ------------------------------------------------------------------
# Atypical features (e.g. Araport11 protein, transposable_element) are kept by aegis tidy
# ------------------------------------------------------------------

from typer.testing import CliRunner
from aegis.cli.tidy import app as tidy_app

ATYPICAL_ROWS = [
    "Chr1|Araport11|gene|3631|5899|.|+|.|ID=AT1G01010;Name=AT1G01010",
    "Chr1|Araport11|mRNA|3631|5899|.|+|.|ID=AT1G01010.1;Parent=AT1G01010",
    "Chr1|Araport11|protein|3760|5630|.|+|.|ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1",
    "Chr1|Araport11|exon|3631|5899|.|+|.|ID=AT1G01010:exon:1;Parent=AT1G01010.1",
    "Chr1|Araport11|transposable_element|11897|11976|.|+|.|ID=AT1TE00010;Name=AT1TE00010",
    "Chr1|Araport11|transposon_fragment|11897|11976|.|+|.|ID=AT1TE00010:fragment:1;Parent=AT1TE00010",
]


def _rows(text, feature):
    return [line.split("\t")[8] for line in text.splitlines() if line.count("\t") == 8 and line.split("\t")[2] == feature]


class TestAtypicalFeaturesKept:

    def test_tidy_keeps_them_by_default(self, tmp_path):
        f = write(tmp_path, ATYPICAL_ROWS)
        res = CliRunner().invoke(tidy_app, [f, "-d", str(tmp_path / "out"), "-o", "out.gff3"])
        assert res.exit_code == 0, res.output
        text = (tmp_path / "out" / "out.gff3").read_text()
        assert _rows(text, "protein") == ["ID=AT1G01010.1-Protein;Derives_from=AT1G01010.1;Name=AT1G01010.1"]
        assert _rows(text, "transposon_fragment") == ["ID=AT1TE00010:fragment:1;Parent=AT1TE00010"]
        assert len(_rows(text, "transposable_element")) == 1

    def test_clean_features_skips_them(self, tmp_path):
        f = write(tmp_path, ATYPICAL_ROWS)
        res = CliRunner().invoke(tidy_app, [f, "--clean-features", "-d", str(tmp_path / "out"), "-o", "out.gff3"])
        assert res.exit_code == 0, res.output
        text = (tmp_path / "out" / "out.gff3").read_text()
        assert _rows(text, "protein") == [] and _rows(text, "transposable_element") == []
        assert "not genes, transcripts or their subfeatures were not read" in res.output

    def test_derives_from_survives_standardised_keys(self, tmp_path):
        annot = Annotation(write(tmp_path, ATYPICAL_ROWS), quiet=True, skip_atypical_features=False)
        annot.standardise_attribute_keys()
        protein = next(f for f in annot.atypical_features if f.feature == "protein")
        assert "Derives_from=AT1G01010.1" in protein.print_gff()

    def test_renaming_relinks_their_parents(self, tmp_path):
        annot = Annotation(write(tmp_path, ATYPICAL_ROWS), quiet=True, skip_atypical_features=False)
        annot.rename_ids(prefix="TEST", quiet=True)
        (new_transcript,) = [t.id for gs in annot.chrs.values() for g in gs.values() for t in g.transcripts.values()]
        assert new_transcript != "AT1G01010.1"
        protein = next(f for f in annot.atypical_features if f.feature == "protein")
        assert protein.parents == [new_transcript]
        fragment = next(f for f in annot.atypical_features if f.feature == "transposon_fragment")
        assert fragment.parents == ["AT1TE00010"]

    def test_written_next_to_the_features_they_point_to(self, tmp_path):
        annot = Annotation(write(tmp_path, ATYPICAL_ROWS + [
            "Chr1|Araport11|gene|20000|21000|.|+|.|ID=AT1G01030",
            "Chr1|Araport11|mRNA|20000|21000|.|+|.|ID=AT1G01030.1;Parent=AT1G01030",
            "Chr1|Araport11|exon|20000|21000|.|+|.|ID=AT1G01030:exon:1;Parent=AT1G01030.1",
        ]), quiet=True, skip_atypical_features=False)
        annot.export.gff(output_dir=str(tmp_path), filename="out.gff3", quiet=True)
        lines = [l for l in (tmp_path / "out.gff3").read_text().splitlines() if l == "###" or not l.startswith("#")]
        features = [l.split("\t")[2] if "\t" in l else l for l in lines]
        # each '###' closes a group: nothing after it refers to a feature before it
        assert features == ["gene", "mRNA", "exon", "protein", "###", "gene", "mRNA", "exon", "###", "transposable_element", "transposon_fragment"]


# ------------------------------------------------------------------
# Removals take their atypical features with them; gene spans follow removed transcripts
# ------------------------------------------------------------------

from aegis.cli import app as aegis_app

ISOFORM_ROWS = [
    "Chr1|A|gene|100|900|.|+|.|ID=g1",
    "Chr1|A|mRNA|100|900|.|+|.|ID=t1;Parent=g1",
    "Chr1|A|exon|100|900|.|+|.|ID=t1:e1;Parent=t1",
    "Chr1|A|CDS|130|870|.|+|0|ID=t1:c1;Parent=t1",
    "Chr1|A|protein|130|870|.|+|.|ID=t1-Protein;Derives_from=t1",
    "Chr1|A|mRNA|200|700|.|+|.|ID=t2;Parent=g1",
    "Chr1|A|exon|200|700|.|+|.|ID=t2:e1;Parent=t2",
    "Chr1|A|CDS|230|670|.|+|0|ID=t2:c1;Parent=t2",
    "Chr1|A|protein|230|670|.|+|.|ID=t2-Protein;Derives_from=t2",
    "Chr1|A|signal_peptide|230|280|.|+|.|ID=t2-sp;Parent=t2-Protein",
    "Chr2|A|transposable_element|100|500|.|+|.|ID=TE1",
    "Chr2|A|transposon_fragment|100|500|.|+|.|ID=TE1:f1;Parent=TE1",
]


def atypical_ids(annot):
    return sorted(f.id for f in annot.atypical_features)


class TestRemovals:

    def test_removed_transcript_takes_its_atypical_features(self, tmp_path, capsys):
        annot = Annotation(write(tmp_path, ISOFORM_ROWS), quiet=True)
        annot.remove_transcripts({"t2"})
        assert atypical_ids(annot) == ["TE1", "TE1:f1", "t1-Protein"]
        assert "2 atypical features (e.g. protein rows) were removed with the genes or transcripts they belong to" in capsys.readouterr().out

    def test_removed_gene_takes_them_too(self, tmp_path):
        annot = Annotation(write(tmp_path, ISOFORM_ROWS), quiet=True)
        annot.remove_genes({"g1"}, quiet=True)
        assert atypical_ids(annot) == ["TE1", "TE1:f1"]

    def test_feature_with_another_parent_keeps_it(self, tmp_path):
        rows = ISOFORM_ROWS[:4] + ["Chr1|A|protein|130|870|.|+|.|ID=shared;Derives_from=t1,t2"] + ISOFORM_ROWS[5:8]
        annot = Annotation(write(tmp_path, rows), quiet=True)
        annot.remove_transcripts({"t2"}, quiet=True)
        (shared,) = [f for f in annot.atypical_features if f.id == "shared"]
        assert shared.parents == ["t1"]

    def test_removed_chromosome_takes_its_atypical_features(self, tmp_path):
        annot = Annotation(write(tmp_path, ISOFORM_ROWS), quiet=True)
        annot.subset(chosen_features={"Chr1"}, gene_cap=None, no_gene_cap=True, min_genes=0, quiet=True)
        assert atypical_ids(annot) == ["t1-Protein", "t2-Protein", "t2-sp"]

    def test_gene_trimmed_to_remaining_transcripts(self, tmp_path):
        annot = Annotation(write(tmp_path, ISOFORM_ROWS), quiet=True)
        annot.remove_transcripts({"t1"}, quiet=True)
        g = annot.chrs["Chr1"]["g1"]
        assert (g.start, g.end) == (200, 700)

    def test_gene_extension_beyond_transcripts_is_kept(self, tmp_path):
        rows = [ISOFORM_ROWS[0].replace("|100|900|", "|50|900|")] + ISOFORM_ROWS[1:]
        annot = Annotation(write(tmp_path, rows), quiet=True)
        annot.remove_transcripts({"t1"}, quiet=True)
        g = annot.chrs["Chr1"]["g1"]
        assert (g.start, g.end) == (50, 700)  # the gene started before its transcripts: that side is kept

    def test_main_only_export_leaves_out_atypical_features_of_other_isoforms(self, tmp_path):
        annot = Annotation(write(tmp_path, ISOFORM_ROWS), quiet=True)
        main = [t.id for t in annot.chrs["Chr1"]["g1"].transcripts.values() if t.main]
        assert main == ["t1"]
        annot.export.gff(output_dir=str(tmp_path), filename="main.gff3", main_only=True, quiet=True)
        written = [l.split("\t")[8].split(";")[0] for l in (tmp_path / "main.gff3").read_text().splitlines() if l.count("\t") == 8]
        assert "ID=t1-Protein" in written and "ID=t2-Protein" not in written and "ID=t2-sp" not in written
        assert "ID=TE1:f1" in written


class TestCleanFeaturesAcrossCommands:

    # (reformat writes GTF from a GFF, and GTF holds gene models only)
    @pytest.mark.parametrize("command", [["rename"], ["prune", "-i", "t2", "-f", "transcript"]])
    def test_kept_by_default_and_removed_with_clean_features(self, tmp_path, command):
        f = write(tmp_path, ISOFORM_ROWS)
        for clean in (False, True):
            out = tmp_path / ("clean" if clean else "kept")
            res = CliRunner().invoke(aegis_app, command[:1] + [f] + command[1:] + (["--clean-features"] if clean else []) + ["-d", str(out), "-o", "out.gff3"])
            assert res.exit_code == 0, res.output
            text = (out / "out.gff3").read_text()
            assert ("\ttransposable_element\t" in text) is not clean
            if command[0] == "prune" and not clean:
                assert "ID=t2-Protein" not in text and "ID=t1-Protein" in text


class TestMergeCarriesAtypicalFeatures:

    BASE = [
        "Chr1|A|gene|100|900|.|+|.|ID=g1",
        "Chr1|A|mRNA|100|900|.|+|.|ID=t1;Parent=g1",
        "Chr1|A|exon|100|900|.|+|.|ID=t1:e1;Parent=t1",
        "Chr1|A|CDS|130|870|.|+|0|ID=t1:c1;Parent=t1",
        "Chr1|A|protein|130|870|.|+|.|ID=base-Protein;Derives_from=t1",
        "Chr2|A|transposable_element|100|500|.|+|.|ID=TE1;Note=base",
    ]
    # same gene and transcript ids as the base, elsewhere: they clash until renaming
    OTHER = [
        "Chr1|B|gene|5000|5900|.|-|.|ID=g1",
        "Chr1|B|mRNA|5000|5900|.|-|.|ID=t1;Parent=g1",
        "Chr1|B|exon|5000|5900|.|-|.|ID=t1:e1;Parent=t1",
        "Chr1|B|CDS|5030|5870|.|-|0|ID=t1:c1;Parent=t1",
        "Chr1|B|protein|5030|5870|.|-|.|ID=other-Protein;Derives_from=t1",
        "Chr2|B|transposable_element|100|500|.|+|.|ID=TE1;Note=other",
        "Chr2|B|transposable_element|900|1500|.|+|.|ID=TE2",
        "Chr2|B|transposon_fragment|900|1500|.|+|.|ID=TE2:f1;Parent=TE2",
    ]

    def test_atypical_features_follow_their_genes(self, tmp_path):
        base = Annotation(write(tmp_path, self.BASE, name="base.gff3"), quiet=True)
        other = Annotation(write(tmp_path, self.OTHER, name="other.gff3"), quiet=True)
        base.merge(other, quiet=True)
        transcripts = {t.id: (t.start, t.strand) for gs in base.chrs.values() for g in gs.values() for t in g.transcripts.values()}
        assert len(transcripts) == 2
        by_id = {f.id: f for f in base.atypical_features}
        assert sorted(by_id) == ["TE1", "TE2", "TE2:f1", "base-Protein", "other-Protein"]
        # each protein points at its own transcript, under its new id
        assert transcripts[by_id["base-Protein"].parents[0]] == (100, "+")
        assert transcripts[by_id["other-Protein"].parents[0]] == (5000, "-")
        assert "Note=base" in by_id["TE1"].print_gff()  # self has priority on a shared id
        assert by_id["TE2:f1"].parents == ["TE2"]


class TestRemovedParentReferences:

    def test_one_summary_line_instead_of_one_per_reference(self, tmp_path, capsys):
        f = write(tmp_path, [
            "Chr1|A|gene|100|900|.|+|.|ID=g1",
            "Chr1|A|mRNA|100|900|.|+|.|ID=t1;Parent=g1",
            "Chr1|A|mRNA|100|900|.|+|.|ID=t2;Parent=g1",
            "Chr1|A|exon|100|300|.|+|.|ID=e1;Parent=t1,t2",
            "Chr1|A|exon|500|900|.|+|.|ID=e2;Parent=t1,t2",
            "Chr1|A|exon|600|900|.|+|.|ID=e3;Parent=t2",
        ])
        annot = Annotation(f, quiet=True)
        capsys.readouterr()
        annot.remove_transcripts({"t2"})
        out = capsys.readouterr().out
        assert "Transcript level id" not in out
        assert "parent references to transcripts no longer in the annotation were removed from exon and CDS rows" in out
        assert all(e.parents == ["t1"] for e in annot.chrs["Chr1"]["g1"].transcripts["t1"].exons)


class TestRenameMatureMiRNAs:

    def test_mature_mirna_follows_its_renamed_transcript(self, tmp_path):
        annot = Annotation(write(tmp_path, [
            "Chr1|A|miRNA_primary_transcript_gene|100|400|.|+|.|ID=g1",
            "Chr1|A|gene|100|400|.|+|.|ID=g1",
            "Chr1|A|miRNA_primary_transcript|100|400|.|+|.|ID=t1;Parent=g1",
            "Chr1|A|exon|100|400|.|+|.|ID=t1:e1;Parent=t1",
            "Chr1|A|miRNA|200|221|.|+|.|ID=miR1;Parent=t1",
        ][1:]), quiet=True)
        annot.rename_ids(prefix="TEST", quiet=True)
        (t,) = [t for gs in annot.chrs.values() for g in gs.values() for t in g.transcripts.values()]
        assert t.id != "t1" and [m.parents for m in t.miRNAs] == [[t.id]]


class TestMergedGenesKeepAtypicalFeatures:

    def test_atypical_features_follow_merged_genes_through_chains(self, tmp_path):
        rows = [
            "Chr1|A|gene|100|900|.|+|.|ID=gA",
            "Chr1|A|mRNA|100|900|.|+|.|ID=tA;Parent=gA",
            "Chr1|A|exon|100|900|.|+|.|ID=tA:e1;Parent=tA",
            "Chr1|A|region_of_interest|100|200|.|+|.|ID=roi1;Parent=gA",
            "Chr1|A|region_of_interest|300|400|.|+|.|ID=roi2;Parent=gB,gC",
            "Chr1|A|protein|100|900|.|+|.|ID=pA;Derives_from=tA",
        ]
        annot = Annotation(write(tmp_path, rows), quiet=True)
        # gA merged into gB, which was merged into gC in turn
        annot._repoint_atypical_parents({"gA": "gB", "gB": "gC"})
        by_id = {f.id: f.parents for f in annot.atypical_features}
        assert by_id == {"roi1": ["gC"], "roi2": ["gC"], "pA": ["tA"]}


class TestReformatRemoveNames:

    ROWS = [
        "Chr1|A|gene|100|900|.|+|.|ID=g1;Name=GENE1",
        "Chr1|A|mRNA|100|900|.|+|.|ID=t1;Parent=g1;Name=GENE1.1",
        "Chr1|A|exon|100|900|.|+|.|ID=t1:e1;Parent=t1",
    ]

    @pytest.mark.parametrize("ext", ["gtf", "gff3"])
    def test_names_written_unless_removed(self, tmp_path, ext):
        f = write(tmp_path, self.ROWS)
        for remove in (False, True):
            out = tmp_path / ("removed" if remove else "kept")
            res = CliRunner().invoke(aegis_app, ["reformat", f, "-d", str(out), "-o", f"out.{ext}"] + (["--remove-names"] if remove else []) + ([] if ext == "gtf" else ["--to", "gff"]))
            assert res.exit_code == 0, res.output
            text = next(out.glob("out.*")).read_text()
            name = 'gene_name "GENE1"' if ext == "gtf" else "Name=GENE1"
            assert (name in text) is not remove
            assert ("transcript_name" in text or "Name=GENE1.1" in text) is not remove
