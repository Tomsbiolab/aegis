"""
Loading GFF3/GTF files whose features are not (or not fully) linked by Parent attributes:
crash fixes for CDS rows without Parent, ##FASTA sections, AUGUSTUS bare IDs, repeated lines,
and linking by a shared attribute (locus_tag by default, see Annotation(link_by=...)).
Several inputs are reduced versions of AGAT's gff_syntax test cases.
"""

import pytest

from aegis.annotation import Annotation
from aegis.utils.gtf_gff import detect_file_format

HEADER = "##gff-version 3\n"


def write(tmp_path, name, rows, header=HEADER):
    f = tmp_path / name
    f.write_text(header + "\n".join("\t".join(r.split("|")) for r in rows) + "\n")
    return str(f)


def models(annot):
    """{gene_id: {transcript_id: (exon chain, CDS chain)}}"""
    out = {}
    for genes in annot.chrs.values():
        for g in genes.values():
            out[g.id] = {}
            for t in g.transcripts.values():
                exons = tuple(sorted((e.start, e.end) for e in t.exons))
                cds = tuple(sorted((s.start, s.end) for c in t.CDSs.values() for s in c.CDS_segments))
                out[g.id][t.id] = (exons, cds)
    return out


# ------------------------------------------------------------------
# Step 0: crashes on the default tidy path
# ------------------------------------------------------------------

class TestCrashFixes:

    def test_cds_without_parent_with_parent_reference_cleanup(self, tmp_path):
        f = write(tmp_path, "cds_no_parent.gff3", [
            "chr1|src|CDS|100|200|.|+|0|ID=cds1",
            "chr1|src|CDS|300|400|.|+|1|ID=cds1",
        ])
        annot = Annotation(f, quiet=True, remove_missing_transcript_parent_references=True)
        t = annot.chrs["chr1"]["cds1_gene"].transcripts["cds1_transcript"]
        for c in t.CDSs.values():
            for s in c.CDS_segments:
                assert s.parents == ["cds1_transcript"]

    def test_exon_parent_in_another_gene_with_parent_reference_cleanup(self, tmp_path):
        f = write(tmp_path, "cross_gene_exon.gff3", [
            "chr1|src|gene|100|500|.|+|.|ID=g1",
            "chr1|src|mRNA|100|500|.|+|.|ID=t1;Parent=g1",
            "chr1|src|gene|100|500|.|+|.|ID=g2",
            "chr1|src|mRNA|100|500|.|+|.|ID=t2;Parent=g2",
            "chr1|src|exon|100|500|.|+|.|ID=e1;Parent=t1,t2",
        ])
        annot = Annotation(f, quiet=True, remove_missing_transcript_parent_references=True)
        assert set(annot.all_transcript_ids) == {"t1", "t2"}


# ------------------------------------------------------------------
# Step 1: parsing gaps
# ------------------------------------------------------------------

class TestParsing:

    def test_fasta_section_ends_the_annotation(self, tmp_path):
        f = write(tmp_path, "with_fasta.gff3", [
            "chr1|src|gene|100|500|.|+|.|ID=g1",
            "chr1|src|mRNA|100|500|.|+|.|ID=t1;Parent=g1",
            "chr1|src|exon|100|500|.|+|.|ID=e1;Parent=t1",
            "##FASTA",
            ">chr1",
            "ACGTACGTACGT",
        ])
        annot = Annotation(f, quiet=True)
        assert list(annot.all_gene_ids) == ["g1"]
        assert not any(h.startswith(">") or h == "ACGTACGTACGT" for h in annot.gff_header)

    def test_augustus_bare_gene_id_gff3(self, tmp_path):
        f = write(tmp_path, "augustus.gff3", [
            "4|AUGUSTUS|gene|100|900|0.01|+|.|g2",
            "4|AUGUSTUS|transcript|100|900|0.01|+|.|g2.t1;Parent=g2",
            "4|AUGUSTUS|exon|100|300|.|+|.|Parent=g2.t1",
            "4|AUGUSTUS|exon|500|900|.|+|.|Parent=g2.t1",
        ], header="")
        assert detect_file_format(f, "utf-8") == "gff3"
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"g2": {"g2.t1": (((100, 300), (500, 900)), ())}}

    def test_augustus_bare_transcript_id_gtf(self, tmp_path):
        f = write(tmp_path, "tsebra.gtf", [
            "s1|AUGUSTUS|gene|100|900|.|+|.|g_1",
            "s1|AUGUSTUS|transcript|100|900|0.7|+|.|anno1.g5.t1",
            's1|AUGUSTUS|CDS|100|300|.|+|0|transcript_id "anno1.g5.t1"; gene_id "g_1";',
            's1|AUGUSTUS|exon|100|300|.|+|.|transcript_id "anno1.g5.t1"; gene_id "g_1";',
            's1|AUGUSTUS|CDS|500|900|.|+|2|transcript_id "anno1.g5.t1"; gene_id "g_1";',
            's1|AUGUSTUS|exon|500|900|.|+|.|transcript_id "anno1.g5.t1"; gene_id "g_1";',
        ], header="")
        annot = Annotation(f, quiet=True)
        assert list(models(annot)) == ["g_1"]
        assert list(models(annot)["g_1"]) == ["anno1.g5.t1"]

    def test_gtf_subfeature_with_only_gene_id_joins_the_single_transcript(self, tmp_path):
        f = write(tmp_path, "jgi.gtf", [
            's5|JGI|exon|100|300|.|+|.|gene_id "e_gw1"; transcript_id 3359',
            's5|JGI|CDS|100|300|.|+|0|gene_id "e_gw1"; protein_id 3358; exonNumber 1',
            's5|JGI|exon|500|900|.|+|.|gene_id "e_gw1"; transcript_id 3359',
            's5|JGI|CDS|500|900|.|+|2|gene_id "e_gw1"; protein_id 3358; exonNumber 2',
        ], header="")
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"e_gw1": {"3359": (((100, 300), (500, 900)), ((100, 300), (500, 900)))}}

    def test_processed_transcript_is_a_transcript(self, tmp_path):
        f = write(tmp_path, "ensembl.gff3", [
            "1|ens|gene|100|500|.|+|.|ID=gene:G1",
            "1|ens|processed_transcript|100|500|.|+|.|ID=transcript:T1;Parent=gene:G1",
            "1|ens|exon|100|500|.|+|.|Parent=transcript:T1",
        ])
        annot = Annotation(f, quiet=True)
        assert list(models(annot)["gene:G1"]) == ["transcript:T1"]

    def test_exact_repeat_of_a_gene_block_is_skipped(self, tmp_path):
        block = [
            "1|ens|gene|100|500|.|+|.|ID=gene:G1",
            "1|ens|mRNA|100|500|.|+|.|ID=transcript:T1;Parent=gene:G1",
            "1|ens|exon|100|500|.|+|.|Parent=transcript:T1",
        ]
        f = write(tmp_path, "repeated.gff3", block + block)
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"gene:G1": {"transcript:T1": (((100, 500),), ())}}
        assert annot.warnings["repeated_feature_line"] == {"gene:G1", "transcript:T1"}

    def test_repeated_gene_id_with_different_line_is_still_kept(self, tmp_path):
        f = write(tmp_path, "repeated_id.gff3", [
            "1|src|gene|100|500|.|+|.|ID=G1",
            "1|src|gene|900|1500|.|+|.|ID=G1",
        ])
        annot = Annotation(f, quiet=True)
        assert len(annot.all_gene_ids) == 2

    def test_ensembl_pseudogene_line_with_parent_is_a_transcript(self, tmp_path):
        f = write(tmp_path, "ensembl_pseudo.gff3", [
            "1|havana|pseudogene|100|500|.|+|.|ID=gene:G1",
            "1|havana|pseudogene|100|500|.|+|.|ID=transcript:T1;Parent=gene:G1",
            "1|havana|exon|100|500|.|+|.|Parent=transcript:T1",
        ])
        annot = Annotation(f, quiet=True)
        assert list(models(annot)) == ["gene:G1"]
        assert list(models(annot)["gene:G1"]) == ["transcript:T1"]

    def test_records_without_id_or_parent_do_not_collide(self, tmp_path):
        f = write(tmp_path, "mitos.gff3", [
            "mt|mitfi|tRNA|15|82|.|+|.|Name=trnF",
            "mt|mitfi|rRNA|82|1063|.|+|.|Name=rrnS",
            "mt|mitfi|tRNA|1063|1132|.|+|.|Name=trnV",
        ])
        annot = Annotation(f, quiet=True)
        assert len(annot.all_transcript_ids) == 3

    def test_orphan_exon_gets_its_own_gene(self, tmp_path):
        f = write(tmp_path, "orphan_exon.gff3", [
            "chr1|src|exon|100|200|.|+|.|ID=lonely_exon",
        ])
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"lonely_exon_gene": {"lonely_exon_transcript": (((100, 200),), ())}}


# ------------------------------------------------------------------
# Step 2: linking by a shared attribute
# ------------------------------------------------------------------

MAKER_LOCUS = [
    "s789|maker|gene|1000|9000|.|+|.|ID=G1;locus_tag=AAA",
    "s789|maker|mRNA|1000|9000|.|+|.|ID=T1;locus_tag=AAA",
    "s789|maker|exon|1000|2000|.|+|.|ID=T1:exon:1;locus_tag=AAA",
    "s789|maker|exon|5000|9000|.|+|.|ID=T1:exon:2;locus_tag=AAA",
    "s789|maker|CDS|1500|2000|.|+|0|ID=T1:cds;locus_tag=AAA",
    "s789|maker|CDS|5000|6000|.|+|1|ID=T1:cds;locus_tag=AAA",
]


class TestLinkByAttribute:

    def test_locus_tag_links_transcript_and_subfeatures(self, tmp_path):
        f = write(tmp_path, "maker.gff3", MAKER_LOCUS)
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"G1": {"T1": (((1000, 2000), (5000, 9000)), ((1500, 2000), (5000, 6000)))}}
        assert {"T1", "T1:exon:1", "T1:exon:2", "T1:cds"} <= annot.warnings["linked_by_shared_attribute"]

    def test_spread_lines_are_linked(self, tmp_path):
        other = ["s625|maker|gene|1|500|.|+|.|ID=G2;locus_tag=BBB",
                 "s625|maker|mRNA|1|500|.|+|.|ID=T2;locus_tag=BBB",
                 "s625|maker|exon|1|500|.|+|.|ID=T2:exon:1;locus_tag=BBB"]
        f = write(tmp_path, "spread.gff3", MAKER_LOCUS[4:] + other + MAKER_LOCUS[:4])
        annot = Annotation(f, quiet=True)
        assert models(annot)["G1"] == {"T1": (((1000, 2000), (5000, 9000)), ((1500, 2000), (5000, 6000)))}
        assert models(annot)["G2"] == {"T2": (((1, 500),), ())}

    def test_linking_can_be_turned_off(self, tmp_path):
        f = write(tmp_path, "maker.gff3", MAKER_LOCUS)
        annot = Annotation(f, quiet=True, link_by=())
        assert "T1" not in models(annot).get("G1", {})
        assert not annot.warnings["linked_by_shared_attribute"]

    def test_prokka_cds_gets_a_gene_named_after_its_locus_tag(self, tmp_path):
        f = write(tmp_path, "prokka.gff3", [
            "contig1|Prodigal|CDS|476|670|.|-|0|ID=Tob1_00001;locus_tag=Tob1_00001;product=hypothetical protein",
            "contig1|Prodigal|CDS|900|1200|.|+|0|ID=Tob1_00002;locus_tag=Tob1_00002",
        ])
        annot = Annotation(f, quiet=True)
        m = models(annot)
        assert set(m) == {"Tob1_00001", "Tob1_00002"}
        assert list(m["Tob1_00001"].values()) == [(((476, 670),), ((476, 670),))]

    def test_refseq_gene_without_id_and_without_transcripts(self, tmp_path):
        f = write(tmp_path, "refseq.gff3", [
            "NC_1|RefSeq|gene|100|900|.|-|.|locus_tag=AT1G01440",
            "NC_1|RefSeq|exon|700|900|.|-|.|gbkey=mRNA;locus_tag=AT1G01440",
            "NC_1|RefSeq|exon|100|300|.|-|.|gbkey=mRNA;locus_tag=AT1G01440",
            "NC_1|RefSeq|CDS|700|800|.|-|0|locus_tag=AT1G01440",
            "NC_1|RefSeq|CDS|200|300|.|-|1|locus_tag=AT1G01440",
            "NC_1|RefSeq|gene|2000|2900|.|+|.|locus_tag=AT1G01355",
            "NC_1|RefSeq|exon|2000|2900|.|+|.|gbkey=mRNA;locus_tag=AT1G01355",
            "NC_1|RefSeq|CDS|2000|2900|.|+|0|locus_tag=AT1G01355",
        ])
        annot = Annotation(f, quiet=True)
        m = models(annot)
        assert set(m) == {"AT1G01440", "AT1G01355"}
        assert list(m["AT1G01440"].values()) == [(((100, 300), (700, 900)), ((200, 300), (700, 800)))]
        assert list(m["AT1G01355"].values()) == [(((2000, 2900),), ((2000, 2900),))]

    def test_isoforms_without_transcript_lines_are_split(self, tmp_path):
        f = write(tmp_path, "isoforms.gff3", [
            "s789|maker|exon|1000|2000|.|+|.|ID=a:exon:1;locus_tag=AAA",
            "s789|maker|exon|5000|9000|.|+|.|ID=a:exon:2;locus_tag=AAA",
            "s789|maker|exon|1000|2000|.|+|.|ID=b:exon:1;locus_tag=AAA",
            "s789|maker|exon|3000|3500|.|+|.|ID=b:exon:2;locus_tag=AAA",
            "s789|maker|exon|5000|9000|.|+|.|ID=b:exon:3;locus_tag=AAA",
        ])
        annot = Annotation(f, quiet=True)
        m = models(annot)
        assert list(m) == ["AAA"]
        assert sorted(chain for chain, _ in m["AAA"].values()) == [
            ((1000, 2000), (3000, 3500), (5000, 9000)),
            ((1000, 2000), (5000, 9000)),
        ]

    def test_subfeature_picks_the_transcript_named_in_its_id(self, tmp_path):
        f = write(tmp_path, "two_mrnas.gff3", [
            "s1|maker|gene|1000|9000|.|+|.|ID=G1;locus_tag=AAA",
            "s1|maker|mRNA|1000|9000|.|+|.|ID=T1;Parent=G1;locus_tag=AAA",
            "s1|maker|mRNA|1000|9000|.|+|.|ID=T2;Parent=G1;locus_tag=AAA",
            "s1|maker|exon|1000|9000|.|+|.|ID=T1:exon:1;locus_tag=AAA",
            "s1|maker|exon|1000|2000|.|+|.|ID=T2:exon:1;locus_tag=AAA",
            "s1|maker|exon|5000|9000|.|+|.|ID=T2:exon:2;locus_tag=AAA",
        ])
        annot = Annotation(f, quiet=True)
        assert models(annot)["G1"] == {"T1": (((1000, 9000),), ()), "T2": (((1000, 2000), (5000, 9000)), ())}

    def test_value_shared_by_several_genes_is_not_used(self, tmp_path):
        f = write(tmp_path, "shared_locus.gff3", [
            "pHB9|Manual|gene|100|900|.|+|.|ID=EBNA2;locus_tag=EBNAs",
            "pHB9|Manual|mRNA|100|900|.|+|.|ID=t2;Parent=EBNA2",
            "pHB9|Manual|exon|100|900|.|+|.|ID=e2;Parent=t2",
            "pHB9|Manual|gene|100|2000|.|+|.|ID=EBNA1;locus_tag=EBNAs",
            "pHB9|Manual|mRNA|100|2000|.|+|.|ID=t1;Parent=EBNA1",
            "pHB9|Manual|exon|100|2000|.|+|.|ID=e1;Parent=t1",
            "pHB9|Manual|exon|1500|1600|.|+|.|ID=Qp-exon;Name=Qp;locus_tag=EBNAs",
        ])
        annot = Annotation(f, quiet=True)
        m = models(annot)
        assert m["EBNA2"] == {"t2": (((100, 900),), ())}
        assert m["EBNA1"] == {"t1": (((100, 2000),), ())}
        assert "Qp-exon" in annot.warnings["shared_attribute_not_linked"]
        # the exon is kept on its own, as any exon without Parent
        assert m["Qp-exon_gene"] == {"Qp-exon_transcript": (((1500, 1600),), ())}

    def test_mixed_strands_are_not_linked(self, tmp_path):
        f = write(tmp_path, "strands.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1;locus_tag=AAA",
            "s1|src|mRNA|100|900|.|-|.|ID=T1;locus_tag=AAA",
        ])
        annot = Annotation(f, quiet=True)
        assert "T1" in annot.warnings["shared_attribute_not_linked"]
        assert "T1" not in models(annot)["G1"]

    def test_link_by_name(self, tmp_path):
        f = write(tmp_path, "by_name.gff3", [
            "mt|mitos|gene|100|900|.|+|.|Name=nad1",
            "mt|mitos|CDS|100|900|.|+|0|Name=nad1",
        ])
        annot = Annotation(f, quiet=True, link_by=("locus_tag", "Name"))
        m = models(annot)
        assert len(m) == 1
        assert list(next(iter(m.values())).values()) == [(((100, 900),), ((100, 900),))]

    def test_link_by_accepts_a_comma_separated_string(self, tmp_path):
        f = write(tmp_path, "maker.gff3", MAKER_LOCUS)
        annot = Annotation(f, quiet=True, link_by="locus_tag, Name")
        assert "T1" in models(annot)["G1"]

    def test_row_without_parent_takes_the_parent_of_rows_sharing_its_id(self, tmp_path):
        f = write(tmp_path, "typo.gff3", [
            "s1|maker|gene|1000|9000|.|+|.|ID=G1",
            "s1|maker|mRNA|1000|9000|.|+|.|ID=T1;Parent=G1",
            "s1|maker|exon|1000|2000|.|+|.|ID=T1:exon:1;Parent=T1",
            "s1|maker|exon|5000|9000|.|+|.|ID=T1:exon:2;Parent=T1",
            "s1|maker|CDS|1500|2000|.|+|0|ID=T1:cds;Parent=T1",
            "s1|maker|CDS|5000|6000|.|+|1|ID=T1:cds;lParent=T1",
        ])
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"G1": {"T1": (((1000, 2000), (5000, 9000)), ((1500, 2000), (5000, 6000)))}}

    def test_files_with_parent_attributes_are_not_changed(self, tmp_path):
        f = write(tmp_path, "parented.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1;locus_tag=AAA",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1;locus_tag=AAA",
            "s1|src|exon|100|900|.|+|.|ID=E1;Parent=T1;locus_tag=AAA",
            "s1|src|gene|2000|2900|.|+|.|ID=G2;locus_tag=AAA",
        ])
        annot = Annotation(f, quiet=True)
        assert models(annot) == {"G1": {"T1": (((100, 900),), ())}, "G2": {}}
        assert not annot.warnings["linked_by_shared_attribute"]
        assert not annot.warnings["shared_attribute_not_linked"]


# ------------------------------------------------------------------
# Step 3: linking by file order (opt-in last resort)
# ------------------------------------------------------------------

NO_PARENTS = [
    "s625|maker|gene|1000|9000|.|+|.|ID=G1",
    "s625|maker|mRNA|1000|9000|.|+|.|ID=T1",
    "s625|maker|exon|1000|2000|.|+|.|ID=T1:exon:1",
    "s625|maker|exon|5000|9000|.|+|.|ID=T1:exon:2",
    "s625|maker|CDS|1500|2000|.|+|0|ID=T1:cds",
    "s625|maker|CDS|5000|6000|.|+|1|ID=T1:cds",
    "s789|maker|gene|100|900|.|-|.|ID=G2",
    "s789|maker|mRNA|100|900|.|-|.|ID=T2",
    "s789|maker|exon|100|900|.|-|.|ID=T2:exon:1",
]


class TestLinkByFileOrder:

    def test_off_by_default_with_a_hint(self, tmp_path, capsys):
        f = write(tmp_path, "no_parents.gff3", NO_PARENTS)
        annot = Annotation(f)
        assert "T1" not in models(annot).get("G1", {})
        assert {"T1", "T1:exon:1", "T1:cds", "T2"} <= annot.warnings["linkable_by_file_order"]
        assert not annot.warnings["linked_by_file_order"]
        assert "--link-by-file-order" in capsys.readouterr().out

    def test_opt_in_links_sorted_records(self, tmp_path):
        f = write(tmp_path, "no_parents.gff3", NO_PARENTS)
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert models(annot) == {
            "G1": {"T1": (((1000, 2000), (5000, 9000)), ((1500, 2000), (5000, 6000)))},
            "G2": {"T2": (((100, 900),), ())},
        }
        assert {"T1", "T1:exon:1", "T1:exon:2", "T1:cds", "T2", "T2:exon:1"} == annot.warnings["linked_by_file_order"]
        assert not annot.warnings["linkable_by_file_order"]

    def test_no_hint_when_everything_has_a_parent(self, tmp_path):
        f = write(tmp_path, "parented.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1",
            "s1|src|exon|100|900|.|+|.|ID=E1;Parent=T1",
        ])
        annot = Annotation(f, quiet=True)
        assert not annot.warnings["linkable_by_file_order"]

    def test_feature_outside_the_preceding_transcript_is_not_linked(self, tmp_path):
        f = write(tmp_path, "outside.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1",
            "s1|src|exon|100|900|.|+|.|ID=E1;Parent=T1",
            "s1|src|exon|5000|6000|.|+|.|ID=far_exon",
        ])
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert models(annot)["G1"] == {"T1": (((100, 900),), ())}
        assert "far_exon" not in annot.warnings["linked_by_file_order"]

    def test_feature_on_the_other_strand_is_not_linked(self, tmp_path):
        f = write(tmp_path, "strand.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1",
            "s1|src|exon|100|900|.|+|.|ID=E1;Parent=T1",
            "s1|src|exon|200|300|.|-|.|ID=antisense_exon",
        ])
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert models(annot)["G1"] == {"T1": (((100, 900),), ())}

    def test_overlapping_subfeature_of_the_same_type_ends_the_chain(self, tmp_path):
        # a second isoform written without its mRNA line: file order cannot place it
        f = write(tmp_path, "isoform.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1",
            "s1|src|exon|100|300|.|+|.|ID=E1;Parent=T1",
            "s1|src|exon|600|900|.|+|.|ID=E2;Parent=T1",
            "s1|src|exon|100|300|.|+|.|ID=iso2_e1",
            "s1|src|exon|400|500|.|+|.|ID=iso2_e2",
        ])
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert models(annot)["G1"] == {"T1": (((100, 300), (600, 900)), ())}
        assert not annot.warnings["linked_by_file_order"] & {"iso2_e1", "iso2_e2"}

    def test_shared_attribute_comes_before_file_order(self, tmp_path):
        # the exon follows G2's lines but its locus_tag names G1
        f = write(tmp_path, "hierarchy.gff3", [
            "s1|src|gene|100|9000|.|+|.|ID=G1;locus_tag=AAA",
            "s1|src|mRNA|100|9000|.|+|.|ID=T1;Parent=G1",
            "s1|src|gene|200|800|.|+|.|ID=G2",
            "s1|src|mRNA|200|800|.|+|.|ID=T2;Parent=G2",
            "s1|src|exon|300|400|.|+|.|ID=tagged_exon;locus_tag=AAA",
        ])
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert models(annot)["G1"] == {"T1": (((300, 400),), ())}
        assert "tagged_exon" in annot.warnings["linked_by_shared_attribute"]
        assert "tagged_exon" not in annot.warnings["linked_by_file_order"]

    def test_ambiguous_shared_attribute_is_not_linked_by_file_order(self, tmp_path):
        f = write(tmp_path, "ambiguous.gff3", [
            "s1|src|gene|100|900|.|+|.|ID=G1;locus_tag=SHARED",
            "s1|src|mRNA|100|900|.|+|.|ID=T1;Parent=G1",
            "s1|src|gene|100|900|.|+|.|ID=G2;locus_tag=SHARED",
            "s1|src|mRNA|100|900|.|+|.|ID=T2;Parent=G2",
            "s1|src|exon|300|400|.|+|.|ID=which_exon;locus_tag=SHARED",
        ])
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        assert "which_exon" in annot.warnings["shared_attribute_not_linked"]
        assert "which_exon" not in annot.warnings["linked_by_file_order"]

    def test_gtf_lines_naming_nothing_keep_their_file_order(self, tmp_path):
        # exonerate: gene lines with gene_id, then cds/exon lines with no gene_id or transcript_id
        f = write(tmp_path, "exonerate.gtf", [
            "c1|exonerate|gene|1000|2000|1150|-|.|gene_id 1 ; sequence XP_1",
            "c1|exonerate|cds|1600|2000|.|-|.|",
            "c1|exonerate|exon|1600|2000|.|-|.|insertions 0 ; deletions 0",
            "c1|exonerate|cds|1000|1500|.|-|.|",
            "c1|exonerate|exon|1000|1500|.|-|.|insertions 0 ; deletions 0",
            "c1|exonerate|gene|5000|6000|1131|-|.|gene_id 2 ; sequence XP_1",
            "c1|exonerate|cds|5000|6000|.|-|.|",
            "c1|exonerate|exon|5000|6000|.|-|.|insertions 0 ; deletions 0",
        ], header="")
        annot = Annotation(f, quiet=True, link_by_file_order=True)
        m = models(annot)
        assert list(m["1"].values()) == [(((1000, 1500), (1600, 2000)), ((1000, 1500), (1600, 2000)))]
        assert list(m["2"].values()) == [(((5000, 6000),), ((5000, 6000),))]


class TestColumnEdgeCases:

    def test_empty_attribute_column_is_kept(self, tmp_path):
        f = tmp_path / "empty_col9.gff3"
        f.write_text(HEADER + "s1\tsrc\tgene\t100\t900\t.\t+\t.\tID=G1\n"
                     "s1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=T1;Parent=G1\n"
                     "s1\tsrc\tCDS\t100\t300\t.\t+\t0\t\n")
        annot = Annotation(str(f), quiet=True, link_by_file_order=True)
        assert models(annot)["G1"]["T1"][1] == ((100, 300),)

    def test_line_without_seqid_is_skipped(self, tmp_path):
        f = tmp_path / "no_seqid.gff3"
        f.write_text(HEADER + "\tsrc\tgene\t100\t900\t.\t+\t.\tID=G0\n"
                     "s1\tsrc\tgene\t100\t900\t.\t+\t.\tID=G1\n")
        annot = Annotation(str(f), quiet=True)
        assert list(annot.all_gene_ids) == ["G1"]
