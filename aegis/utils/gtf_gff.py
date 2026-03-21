import sys
import re

from ..conf import default_features, default_subfeatures
from collections import OrderedDict
from dataclasses import dataclass
from .misc import open_file

ALT_ID_gene = {"gene_id", "geneid"}
ALT_ID_transcript = {"transcript_id", "transcriptid"}
ID_OR_PARENT_KEYS = ALT_ID_gene.union(ALT_ID_transcript)
PARENT_KEYS = {"parent", "parents", "derives_from"}
ALIAS_KEYS = {"alias", "gene_synonym"}

INT_ID = sys.intern("ID")
INT_PARENT = sys.intern("Parent")
INT_ALIAS = sys.intern("Alias")
INT_NAME = sys.intern("Name")
INT_SYMBOL = sys.intern("Symbol")

MAPPED_ATTRIBUTES = {
    "parent": INT_PARENT,
    "parents": INT_PARENT,
    "derives_from": INT_PARENT,
    "alias": INT_ALIAS,
    "gene_synonym": INT_ALIAS,
    "name": INT_NAME,
    "symbol": INT_SYMBOL,
}

@dataclass(slots=True)
class GffEntry:
    ch: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    attributes: dict[str, str]
    id: str
    parents: list[str]
    pseudogene: bool
    transposable: bool
    decreasing_coordinates: bool

def parse_gff_parts(parts):

    # Interning high-frequency strings
    source = sys.intern(parts[1])
    feature = sys.intern(parts[2])
    strand = sys.intern(parts[6])

    ch = sys.intern(parts[0].partition(":")[0]) if feature == "nucleotide_to_protein_match" else sys.intern(parts[0])

    pseudogene = "pseudo" in feature
    transposable = "transposable" in feature or "transposon" in feature

    start = int(parts[3])
    end = int(parts[4])

    decreasing_coordinates = start > end
    if decreasing_coordinates:
        actual_start, actual_end = end, start
    else:
        actual_start, actual_end = start, end

    gene = feature in default_features["gene"]
    transcript = feature in default_features["transcript"]

    attr_dict = parse_gff_attributes(parts[8], gene, transcript)
    
    if not transposable:
        if attr_dict.get("transposable") == "True" or attr_dict.get("transposon") == "True":
            transposable = True

    entry_id = attr_dict.pop(INT_ID, "")
    parents = attr_dict.pop(INT_PARENT, [])

    if not entry_id and parents:
        entry_id = f"{parents[0]}_{feature}_{actual_start}_{actual_end}"

    return GffEntry(
        ch=ch,
        source=source,
        feature=feature,
        start=actual_start,
        end=actual_end,
        score=parts[5],
        strand=strand,
        attributes=attr_dict,
        id=entry_id,
        parents=parents,
        pseudogene=pseudogene,
        transposable=transposable,
        decreasing_coordinates=decreasing_coordinates
    )

def parse_gff_attributes(attributes, gene:bool=False, transcript:bool=False):

    if not attributes or attributes == ".":
        return {}

    parsed = {}
    
    for p in attributes.split(";"):

        key, sep, val = p.partition("=")
        if not sep: # Means '=' was not found
            continue
        
        key = key.strip()
        val = val.strip()

        lower_key = key.lower()

        if lower_key == "id":
            if INT_ID not in parsed:
                parsed[INT_ID] = val
            continue
            
        # Promote gene_id/transcript_id ONLY if correct feature type and ID is missing
        if (gene and lower_key in ALT_ID_gene) or (transcript and lower_key in ALT_ID_transcript):
            if INT_ID not in parsed:
                parsed[INT_ID] = val
            continue

        # --- Priority 2: Mapped logic (Parents, Names, Symbols, etc.) ---
        if lower_key in MAPPED_ATTRIBUTES:
            target_key = MAPPED_ATTRIBUTES[lower_key]
            if target_key not in parsed:
                parsed[target_key] = [x.strip() for x in val.split(",") if x.strip()]
            continue
            
        # Handle the "ID_OR_PARENT" fallback logic (e.g. gene_id as Parent)
        if lower_key in ID_OR_PARENT_KEYS:
            if INT_PARENT not in parsed:
                parsed[INT_PARENT] = [x.strip() for x in val.split(",") if x.strip()]
            continue

        # Otherwise, it's a normal attribute
        parsed[sys.intern(key)] = val

    return parsed

def parse_gtf_attributes(attr_string):
    """
    Parses the GFF/GTF attribute column and returns a dictionary.
    Handles GTF-specific format where values are quoted.
    """
    attributes = {}

    for match in re.finditer(r'(\w+)\s+"([^"]+)"', attr_string):
        key, value = match.groups()
        attributes[key] = value
    return attributes

def format_gff3_attributes(attrs, feature_type):
    """
    Formats a dictionary of attributes into a GFF3-compliant string.
    """

    special_keys = {'gene_id', 'transcript_id', 'exon_id', 'exon_number'}
    
    gff3_attrs = []

    if feature_type in default_features["gene"]:
        if 'gene_id' in attrs:
            gff3_attrs.append(f"ID={attrs['gene_id']}")
    
    elif feature_type in default_features["transcript"]:
        if 'transcript_id' in attrs:
            gff3_attrs.append(f"ID={attrs['transcript_id']}")
        if 'gene_id' in attrs:
            gff3_attrs.append(f"Parent={attrs['gene_id']}")
            
    elif feature_type in default_subfeatures:
        parent_id = f"{attrs['transcript_id']}"
        gff3_attrs.append(f"Parent={parent_id}")

        feature_id = f"{attrs['transcript_id']}"

        if 'exon_number' in attrs:
            feature_id += f"_e{attrs['exon_number']}"
            gff3_attrs.append(f"ID={feature_id}")
        elif feature_type in default_features["CDS"]:
            gff3_attrs.append(f"ID={feature_id}_CDS1")

    if 'gene_name' in attrs:
        gff3_attrs.append(f"Symbol={attrs['gene_name']}")
    elif 'transcript_name' in attrs:
        gff3_attrs.append(f"Symbol={attrs['transcript_name']}")

    for key, value in attrs.items():
        if key not in special_keys and key not in ['gene_name', 'transcript_name']:
            gff3_attrs.append(f"{key}={value}")
            
    return ";".join(gff3_attrs)

def convert_gtf_to_gff3(gtf_file, gff3_file, encoding, quiet:bool=False):
    """
    Reads a GTF file, converts it to GFF3 format, and writes to an output file.

    Uses a two-pass approach:
      Pass 1 — Collect all lines and infer gene/transcript boundaries from
               subfeature attributes (gene_id, transcript_id) when no explicit
               gene or transcript rows exist.
      Pass 2 — Emit gene, transcript, and subfeature lines in GFF3 format.
    """

    comments = []
    seen_genes = set()
    seen_transcripts = set()
    subfeature_lines = []

    inferred_genes: dict[str, dict] = OrderedDict()
    inferred_transcripts: dict[str, dict] = OrderedDict()

    explicit_gene_lines = []
    explicit_transcript_lines = []

    # PASS 1 — read file
    with open_file(gtf_file, 'r', encoding=encoding) as infile:
        for line in infile:
            if line.startswith('#'):
                if not line.startswith('##'):
                    comments.append(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                if not quiet:
                    sys.stderr.write(f"Warning: Skipping malformed line: {line.strip()}\n")
                continue

            seqname, source, feature, start, end, score, strand, frame, attr_string = parts
            attributes = parse_gtf_attributes(attr_string)
            start_int = int(start)
            end_int = int(end)

            gene_id = attributes.get('gene_id', '')
            transcript_id = attributes.get('transcript_id', '')

            if feature in default_features["gene"]:
                if gene_id and gene_id not in seen_genes:
                    explicit_gene_lines.append((parts, attributes))
                    seen_genes.add(gene_id)
                continue

            if feature in default_features["transcript"]:
                if transcript_id and transcript_id not in seen_transcripts:
                    explicit_transcript_lines.append((parts, attributes))
                    seen_transcripts.add(transcript_id)
                continue

            subfeature_lines.append((parts, attributes, feature))

            if gene_id:
                if gene_id not in inferred_genes:
                    inferred_genes[gene_id] = {
                        'seqname': seqname, 'source': source,
                        'strand': strand, 'start': start_int,
                        'end': end_int, 'gene_name': attributes.get('gene_name', ''),
                        'gene_biotype': attributes.get('gene_biotype', ''),
                        'transcript_ids': OrderedDict()
                    }
                else:
                    g = inferred_genes[gene_id]
                    if start_int < g['start']:
                        g['start'] = start_int
                    if end_int > g['end']:
                        g['end'] = end_int

            if transcript_id:
                if transcript_id not in inferred_transcripts:
                    inferred_transcripts[transcript_id] = {
                        'gene_id': gene_id, 'seqname': seqname,
                        'source': source, 'strand': strand,
                        'start': start_int, 'end': end_int,
                    }
                else:
                    t = inferred_transcripts[transcript_id]
                    if start_int < t['start']:
                        t['start'] = start_int
                    if end_int > t['end']:
                        t['end'] = end_int

                if gene_id and gene_id in inferred_genes:
                    inferred_genes[gene_id]['transcript_ids'][transcript_id] = True

    # PASS 2 — write GFF3
    with open_file(gff3_file, 'w', encoding=encoding) as outfile:
        outfile.write("##gff-version 3\n")

        for c in comments:
            outfile.write(f"#{c}")

        for parts, attributes in explicit_gene_lines:
            gene_attrs = {'gene_id': attributes['gene_id']}
            if 'gene_name' in attributes:
                gene_attrs['gene_name'] = attributes['gene_name']
            if 'gene_biotype' in attributes:
                gene_attrs['gene_biotype'] = attributes['gene_biotype']
            gene_attr_str = format_gff3_attributes(gene_attrs, 'gene')
            outfile.write("\t".join([parts[0], parts[1], 'gene', parts[3], parts[4],
                                     parts[5], parts[6], parts[7], gene_attr_str]) + '\n')

        for gene_id, g in inferred_genes.items():
            if gene_id in seen_genes:
                continue
            gene_attrs = {'gene_id': gene_id}
            if g['gene_name']:
                gene_attrs['gene_name'] = g['gene_name']
            if g['gene_biotype']:
                gene_attrs['gene_biotype'] = g['gene_biotype']
            gene_attr_str = format_gff3_attributes(gene_attrs, 'gene')
            outfile.write("\t".join([g['seqname'], g['source'], 'gene',
                                     str(g['start']), str(g['end']), '.',
                                     g['strand'], '.', gene_attr_str]) + '\n')
            seen_genes.add(gene_id)

        for parts, attributes in explicit_transcript_lines:
            tx_attrs = {k: v for k, v in attributes.items()
                        if 'transcript' in k or 'gene' in k}
            tx_feature_type = 'transcript'
            if 'transcript_biotype' in attributes:
                if 'RNA' in attributes['transcript_biotype']:
                    tx_feature_type = attributes['transcript_biotype']
            tx_attr_str = format_gff3_attributes(tx_attrs, tx_feature_type)
            outfile.write("\t".join([parts[0], parts[1], tx_feature_type,
                                     parts[3], parts[4], parts[5], parts[6],
                                     parts[7], tx_attr_str]) + '\n')

        for t_id, t in inferred_transcripts.items():
            if t_id in seen_transcripts:
                continue
            tx_attrs = {'transcript_id': t_id}
            if t['gene_id']:
                tx_attrs['gene_id'] = t['gene_id']
            tx_attr_str = format_gff3_attributes(tx_attrs, 'mRNA')
            outfile.write("\t".join([t['seqname'], t['source'], 'mRNA',
                                     str(t['start']), str(t['end']), '.',
                                     t['strand'], '.', tx_attr_str]) + '\n')
            seen_transcripts.add(t_id)

        subfeature_counters: dict[tuple[str, str], int] = {}

        for parts, attributes, feature in subfeature_lines:
            transcript_id = attributes.get('transcript_id', '')
            gene_name = attributes.get('gene_name', '')

            gff3_attrs = []

            if transcript_id:
                key = (transcript_id, feature)
                subfeature_counters[key] = subfeature_counters.get(key, 0) + 1
                count = subfeature_counters[key]
                feature_id = f"{transcript_id}_{feature}{count}"
                gff3_attrs.append(f"ID={feature_id}")
                gff3_attrs.append(f"Parent={transcript_id}")
            else:
                gff3_attr_string = format_gff3_attributes(attributes, feature)
                outfile.write("\t".join([parts[0], parts[1], feature,
                                         parts[3], parts[4], parts[5], parts[6],
                                         parts[7], gff3_attr_string]) + '\n')
                continue

            if gene_name:
                gff3_attrs.append(f"Symbol={gene_name}")

            gff3_attr_string = ";".join(gff3_attrs)
            outfile.write("\t".join([parts[0], parts[1], feature,
                                     parts[3], parts[4], parts[5], parts[6],
                                     parts[7], gff3_attr_string]) + '\n')

    if not quiet:
        print(f"Successfully converted {gtf_file} to a gff format")

def detect_file_format(file_path, encoding, lines_to_check=20):
    """
    Detects if a file is likely GTF or GFF3 format.
    """
    try:
        with open_file(file_path, 'r', encoding=encoding) as f:

            i = 0

            for line in f:

                line = line.strip()
                if not line:
                    continue

                if line.startswith("##gff-version 3"):
                    return 'gff3'

                if line.startswith('#'):
                    continue

                i += 1

                if i >= lines_to_check:
                    break

                parts = line.split('\t')
                if len(parts) == 9:
                    attributes = parts[8]

                    if re.search(r'\w+=', attributes):
                        return 'gff3'

                    if re.search(r'\w+\s+"[^"]+";', attributes):
                        return 'gtf'
            
            # unknown format returns gff3
            return 'gff3'
            
    except Exception as e:
        sys.stderr.write(f"Error reading file for format detection: {e}\n")
        sys.exit(1)