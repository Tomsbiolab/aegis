"""
Malformed GFF3 lines from Araport11 (Araport11_GFF3_genes_transposons.current.gff):
- attribute columns broken by a line break (AT2G12430.1 mRNA and AT2G12430:exon:3, lines 266915/266917),
  joined back instead of being written as header lines;
- text after an unescaped ';' in Note (e.g. AT1G12013 'Note=SNOR111; snoRNA; snoRNA'), kept as '%3B';
- other lines with fewer than 9 columns, skipped with a warning.
"""

from aegis.annotation import Annotation
from aegis.utils.gtf_gff import parse_gff_attributes, FRAGMENTS

DESCRIPTION = "hAT-like transposase family (hobo/Ac/Tam3), has a 4.6e-33 P-value blast match to GB:CAA29005 ORFa "
WRAPPED = "of Maize Ac (hAT-element) (Zea mays)(source:TAIR10)"

ARAPORT_WRAPPED = (
    "##gff-version 3\n"
    "#chromosome assembly: TAIR10.1\n"
    "Chr2\tAraport11\ttransposable_element_gene\t5031199\t5033508\t.\t-\t.\tID=AT2G12430;Name=AT2G12430\n"
    f"Chr2\tAraport11\tmRNA\t5031199\t5033508\t.\t-\t.\tID=AT2G12430.1;Parent=AT2G12430;computational_description={DESCRIPTION}\n"
    f"{WRAPPED};Note=transposable element gene;symbol=None;curator_summary=None\n"
    f"Chr2\tAraport11\texon\t5031199\t5031296\t.\t-\t.\tID=AT2G12430:exon:3;Parent=AT2G12430.1;computational_description={DESCRIPTION}\n"
    f"{WRAPPED}\n"
    "Chr2\tAraport11\texon\t5031381\t5033508\t.\t-\t.\tID=AT2G12430:exon:2;Parent=AT2G12430.1\n"
)

ARAPORT_SEMICOLONS = (
    "##gff-version 3\n"
    "Chr1\tAraport11\tgene\t4058287\t4058416\t69\t+\t.\tID=AT1G12013;Note=SNOR111; snoRNA; snoRNA;symbol=SNOR111;Name=AT1G12013\n"
    "Chr1\tAraport11\tsnoRNA\t4058287\t4058416\t69\t+\t.\tID=AT1G12013.1;Parent=AT1G12013;Note=SNOR111; snoRNA\n"
    "Chr1\tAraport11\texon\t4058287\t4058416\t69\t+\t.\tID=AT1G12013:exon:1;Parent=AT1G12013.1\n"
)


def _write(tmp_path, text, name="in.gff3"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def _export(annot, tmp_path, name="out.gff3"):
    annot.export.gff(output_dir=str(tmp_path), filename=name, quiet=True)
    return (tmp_path / name).read_text(encoding="utf-8")


def _column9(text, feature_id):
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) == 9 and f"ID={feature_id};" in parts[8] + ";":
            return parts[8]
    raise AssertionError(f"No {feature_id} line")


class TestWrappedLines:

    def test_continuation_is_joined_to_its_feature(self, tmp_path):
        annot = Annotation(_write(tmp_path, ARAPORT_WRAPPED), quiet=True)
        assert annot.gff_header == ["##gff-version 3", "#chromosome assembly: TAIR10.1"]
        assert annot.warnings["wrapped_line_joined"] == {"AT2G12430.1", "AT2G12430:exon:3"}
        assert annot.warnings["malformed_line_skipped"] == set()

        out = _export(annot, tmp_path)
        assert all(line.startswith("#") or line.count("\t") == 8 for line in out.splitlines() if line)
        mrna = _column9(out, "AT2G12430.1")
        # the trailing space before the line break is kept: the text reads as in the unbroken exon:2 line of Araport11
        assert "ORFa of Maize Ac (hAT-element)" in mrna
        assert "Note=transposable element gene" in mrna and "symbol=None" in mrna and "curator_summary=None" in mrna

    def test_break_without_spaces_reads_as_a_space(self, tmp_path):
        text = (
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=g1;Note=first\n"
            "half\n"
            "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tsrc\texon\t100\t500\t.\t+\t.\tID=e1;Parent=t1\n"
        )
        out = _export(Annotation(_write(tmp_path, text), quiet=True), tmp_path)
        assert "Note=first half" in _column9(out, "g1")

    def test_several_continuation_lines(self, tmp_path):
        rows = list(Annotation._unwrapped_gff_lines(["a\tb\tc\td\te\tf\tg\th\tNote=x \n", "y \n", "z\n", "#c\n"]))
        assert rows[0] == (1, "a\tb\tc\td\te\tf\tg\th\tNote=x y z", [2, 3])
        assert rows[1] == (4, "#c\n", [])

    def test_short_lines_that_continue_nothing_are_skipped(self, tmp_path):
        text = (
            "##gff-version 3\n"
            "stray text before any feature\n"
            "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\tsrc\tmRNA\t100\t500\n"
            "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tsrc\texon\t100\t500\t.\t+\t.\tID=e1;Parent=t1\n"
        )
        annot = Annotation(_write(tmp_path, text), quiet=True)
        assert annot.gff_header == ["##gff-version 3"]
        assert annot.warnings["malformed_line_skipped"] == {"line 2", "line 4"}
        assert annot.warnings["wrapped_line_joined"] == set()
        assert set(annot.all_transcript_ids) == {"t1"}

    def test_header_export_only_writes_comment_lines(self, tmp_path):
        annot = Annotation(_write(tmp_path, ARAPORT_WRAPPED), quiet=True)
        annot.gff_header.append("not a comment")
        out = _export(annot, tmp_path)
        assert "not a comment" not in out


class TestTextAfterSemicolon:

    def test_parser_appends_pieces_to_the_previous_value(self):
        attrs = parse_gff_attributes("ID=AT1G12013;Note=SNOR111; snoRNA; snoRNA;symbol=SNOR111", gene=True)
        assert attrs["Note"] == "SNOR111%3B snoRNA%3B snoRNA"
        assert attrs["ID"] == "AT1G12013"
        assert attrs[FRAGMENTS] == (2, 0)

    def test_pieces_after_id_or_parent_are_not_appended(self):
        attrs = parse_gff_attributes("ID=e1;stray;Parent=t1")
        assert attrs["ID"] == "e1" and attrs["Parent"] == ["t1"]
        assert attrs[FRAGMENTS] == (0, 1)

    def test_augustus_bare_ids_still_read(self):
        attrs = parse_gff_attributes("g2.t1;Parent=g2", transcript=True)
        assert attrs["ID"] == "g2.t1" and attrs["Parent"] == ["g2"]
        assert FRAGMENTS not in attrs

    def test_repeated_key_is_continued_too(self):
        attrs = parse_gff_attributes("ID=g1;Note=a;Note=b; c")
        assert [v for k, v in attrs.items() if k.startswith("\x00extra")] == ["Note=b%3B c"]

    def test_araport_note_survives_tidy_and_reload(self, tmp_path):
        annot = Annotation(_write(tmp_path, ARAPORT_SEMICOLONS), quiet=True)
        assert annot.warnings["attribute_text_after_semicolon_joined"] == {"AT1G12013", "AT1G12013.1"}
        out = _export(annot, tmp_path)
        assert "Note=SNOR111%3B snoRNA%3B snoRNA" in _column9(out, "AT1G12013")
        assert "Note=SNOR111%3B snoRNA" in _column9(out, "AT1G12013.1")

        again = Annotation(_write(tmp_path, out, name="again.gff3"), quiet=True)
        assert again.warnings["attribute_text_after_semicolon_joined"] == set()
        assert "Note=SNOR111%3B snoRNA%3B snoRNA" in _column9(_export(again, tmp_path, name="again_out.gff3"), "AT1G12013")


ARAPORT_SYMBOL = (
    "##gff-version 3\n"
    "Chr1\tAraport11\tgene\t29333366\t29337588\t.\t-\t.\tID=AT1G78000;Name=AT1G78000;symbol=SULTR1;2;full_name=sulfate transporter 1;2;Alias=SEL1;x\n"
    "Chr1\tAraport11\tmRNA\t29333366\t29337588\t.\t-\t.\tID=AT1G78000.1;Parent=AT1G78000;Name=AT1G78000.1;RNA binding\n"
    "Chr1\tAraport11\texon\t29333366\t29337588\t.\t-\t.\tID=AT1G78000:exon:1;Parent=AT1G78000.1\n"
)


class TestSymbolAndAliasWithSemicolon:

    def test_parser_appends_to_the_last_symbol_or_alias_value(self):
        attrs = parse_gff_attributes("ID=g1;symbol=A1,SULTR1;2;2;Alias=X%2CY;1", gene=True)
        assert attrs["Symbol"] == ["A1", "SULTR1;2;2"]
        assert attrs["Alias"] == ["X,Y;1"]
        assert attrs[FRAGMENTS] == (3, 0)

    def test_repeated_symbol_key_continues_the_value_just_read(self):
        attrs = parse_gff_attributes("ID=g1;symbol=A;symbol=B;1", gene=True)
        assert attrs["Symbol"] == ["A", "B;1"]

    def test_text_after_name_is_still_dropped(self):
        attrs = parse_gff_attributes("ID=g1;Name=AT1G01760;RNA binding;adenosine deaminases", gene=True)
        assert attrs["Name"] == ["AT1G01760"]
        assert attrs[FRAGMENTS] == (0, 2)

    def test_leftover_of_a_lost_note_is_not_a_name(self):
        # Araport11 AT1G20860: '8' is what is left of a Note, not part of the Name
        attrs = parse_gff_attributes("ID=AT1G20860;Name=AT1G20860;8;8;8, a member of the Pht1 family;locus=1", gene=True)
        assert attrs["Name"] == ["AT1G20860"]
        assert attrs[FRAGMENTS] == (0, 3)

    def test_name_cut_at_semicolon_is_joined_when_it_continues_like_the_id(self):
        # Araport11 CDS/UTR rows name the symbol: 'Name=CYCA2;3:CDS:1' is 'CYCA2;3' + ':CDS:1'
        attrs = parse_gff_attributes("ID=AT1G15570:CDS:1;Parent=AT1G15570.1;Name=CYCA2;3:CDS:1;Note=CYCLIN A2;3")
        assert attrs["Name"] == ["CYCA2;3:CDS:1"]
        assert attrs["Note"] == "CYCLIN A2%3B3"
        assert attrs[FRAGMENTS] == (2, 0)
        # CDS shared by isoforms: the ID carries an extra '.2'
        attrs = parse_gff_attributes("ID=AT1G16780:CDS:14.2;Parent=AT1G16780.2;Name=VHP2;2:CDS:14")
        assert attrs["Name"] == ["VHP2;2:CDS:14"]
        attrs = parse_gff_attributes("Name=PIP1;3:five_prime_UTR:2;ID=AT1G01030:five_prime_UTR:2;Parent=t1")
        assert attrs["Name"] == ["PIP1;3:five_prime_UTR:2"]

    def test_name_is_not_joined_without_a_matching_id(self):
        for column9 in ("ID=g1;Name=CYCA2;3:CDS:1", "ID=g1:CDS:1;Name=CYCA2;3:1", "ID=g1:CDS:1;Name=CYCA2;3:CDS:1 x"):
            attrs = parse_gff_attributes(column9)
            assert attrs["Name"] == ["CYCA2"], column9
            assert attrs[FRAGMENTS] == (0, 1), column9

    def test_bare_id_after_name_still_read_as_id(self):
        attrs = parse_gff_attributes("Name=x;g2", gene=True)
        assert attrs["ID"] == "g2" and attrs["Name"] == ["x"]
        assert FRAGMENTS not in attrs

    def test_araport_cds_name_survives_export_and_reload(self, tmp_path):
        text = (
            "##gff-version 3\n"
            "Chr1\tAraport11\tgene\t100\t400\t.\t+\t.\tID=AT1G15570;Name=AT1G15570\n"
            "Chr1\tAraport11\tmRNA\t100\t400\t.\t+\t.\tID=AT1G15570.1;Parent=AT1G15570;Name=AT1G15570.1\n"
            "Chr1\tAraport11\texon\t100\t400\t.\t+\t.\tID=AT1G15570:exon:1;Parent=AT1G15570.1\n"
            "Chr1\tAraport11\tCDS\t100\t399\t.\t+\t0\tID=AT1G15570:CDS:1;Parent=AT1G15570.1;Name=CYCA2;3:CDS:1\n"
        )
        annot = Annotation(_write(tmp_path, text), quiet=True)
        assert annot.warnings["attribute_text_without_key_dropped"] == set()
        out = _export(annot, tmp_path)
        assert "Name=CYCA2%3B3:CDS:1" in out
        again = Annotation(_write(tmp_path, out, name="again.gff3"), quiet=True)
        assert again.warnings["attribute_text_after_semicolon_joined"] == set()
        assert _export(again, tmp_path, name="again_out.gff3") == out

    def test_symbol_survives_tidy_reload_and_gtf(self, tmp_path):
        annot = Annotation(_write(tmp_path, ARAPORT_SYMBOL), quiet=True)
        gene = annot.chrs["Chr1"]["AT1G78000"]
        assert gene.symbols == ["SULTR1;2"] and gene.aliases == ["SEL1;x"]
        assert annot.warnings["attribute_text_after_semicolon_joined"] == {"AT1G78000"}
        assert annot.warnings["attribute_text_without_key_dropped"] == {"AT1G78000.1"}

        out = _export(annot, tmp_path)
        column9 = _column9(out, "AT1G78000")
        assert "symbol=SULTR1%3B2" in column9 and "Alias=SEL1%3Bx" in column9
        assert "full_name=sulfate transporter 1%3B2" in column9

        again = Annotation(_write(tmp_path, out, name="again.gff3"), quiet=True)
        assert again.chrs["Chr1"]["AT1G78000"].symbols == ["SULTR1;2"]

        annot.export.gtf(output_dir=str(tmp_path), filename="out.gtf", quiet=True)
        gtf = (tmp_path / "out.gtf").read_text(encoding="utf-8")
        # ';' is written as '%3B' in GTF too: readers that split column 9 on ';' (e.g. gffread) would cut 'SULTR1;2'
        assert 'gene_name "SULTR1%3B2"' in gtf and ";2" not in gtf
        from_gtf = Annotation(str(tmp_path / "out.gtf"), quiet=True)
        assert from_gtf.chrs["Chr1"]["AT1G78000"].symbol_values() == ["SULTR1;2"]

    def test_gtf_id_with_semicolon_does_not_break_the_line(self, tmp_path):
        gtf = tmp_path / "in.gtf"
        gtf.write_text(
            'Chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id "g%3B1"; transcript_id "t1"; gene_name "SULTR1%3B2";\n'
            'Chr1\tsrc\texon\t100\t500\t.\t+\t.\tgene_id "g%3B1"; transcript_id "t1";\n'
        )
        annot = Annotation(str(gtf), quiet=True)
        (gene,) = [g for gs in annot.chrs.values() for g in gs.values()]
        assert gene.id == "g%3B1" and sorted(gene.transcripts) == ["t1"]
