from aegis import *

a = Annotation("file_path.gff3")

extra_3_exon = 500

for chrom, genes in a.chrs.items():
    for g in genes.items():
        for t in g.transcripts.items():
            if t.strand == "+":
                start = t.exons[-1].end + 1
                end = t.exons[-1].end + extra_3_exon
            elif t.strand == "-":
                start = t.exons[0].start - extra_3_exon
                end = t.exons[0].start - 1
            attributes = f"ID=temp_id;Parent={t.id}"

            t.exons.append(Exon("temp_id", chrom, t.source, "exon", t.strand, start, end, t.score, t.phase, attributes))

# crucial to generate UTRs themselves
a.update()

a.rename_ids(features=["exon", "UTR"])

a.export.gff("file_path_out.gff3", UTRs=True)


