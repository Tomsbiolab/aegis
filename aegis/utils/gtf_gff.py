import sys
import os
import re

from ..conf import default_features, default_subfeatures
from collections import OrderedDict
from dataclasses import dataclass
from urllib.parse import unquote
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

# Free-text values read percent-decoded and written percent-encoded (GFF3 column 9).
# ID and Parent are kept as opaque tokens, exactly as read.
DECODED_ATTRIBUTES = {INT_NAME, INT_SYMBOL, INT_ALIAS}

# Symbol and Alias values may hold an unescaped ';' (Araport11 'symbol=SULTR1;2'): the text after it is
# appended back to the value. Not Name, after which such text is a description that lost its key.
CONTINUED_ATTRIBUTES = {INT_SYMBOL, INT_ALIAS}

# Keys written for these attributes unless the input spelt them differently (e.g. 'gene_synonym', 'Symbol').
# Lowercase 'symbol': GFF3 reserves capitalised attribute names and Symbol is not one of them.
DEFAULT_WRITE_KEYS = {INT_NAME: "Name", INT_ALIAS: "Alias", INT_SYMBOL: "symbol", INT_PARENT: "Parent"}

# Entries of parse_gff_attributes() that are not attributes ('\x00' cannot start a GFF key)
SOURCE_KEYS = "\x00source_keys"    # {INT_NAME/INT_ALIAS/INT_SYMBOL: key as spelt in the input}
EXTRA_ATTRIBUTES = "\x00extra"     # prefix of numbered keys holding one raw 'key=value' (repeated or unused key) each, at its column-9 position
PARENT_KEY = "\x00parent_key"      # key the Parent attribute was read under when not 'Parent' (e.g. 'Derives_from'), see GffEntry.parent_key
FRAGMENTS = "\x00fragments"        # (joined, dropped): pieces of column 9 with no '=' (text after an unescaped ';'), see parse_gff_attributes

def extra_attributes(parsed: dict) -> list[str]:
    """Raw 'key=value' strings kept under EXTRA_ATTRIBUTES keys, in column-9 order."""
    return [v for k, v in parsed.items() if k.startswith(EXTRA_ATTRIBUTES)]

# GFF3 reserved characters plus control characters; one translate pass, so '%' needs no special ordering
_GFF3_ESCAPES = {i: f"%{i:02X}" for i in range(0x20)}
_GFF3_ESCAPES.update({ord(c): f"%{ord(c):02X}" for c in "%;=&,\x7f"})

def gff3_escape(value: str) -> str:
    """Percent-encodes one GFF3 attribute value."""
    return value.translate(_GFF3_ESCAPES)

def gff3_unescape(value: str) -> str:
    """Percent-decodes one GFF3 attribute value ('+' is left as is)."""
    return unquote(value) if "%" in value else value

_GTF_UNESCAPE = re.compile(r'\\(["\\])')

_GTF_SEMICOLON = re.compile("%3B", re.IGNORECASE)

def gtf_quote(value: str) -> str:
    """
    Quotes one GTF attribute value, backslash-escaping '"' and '\\' so the value cannot end early.
    ';' is written as '%3B' (e.g. the symbol 'SULTR1;2'), as GTF readers that split column 9 on ';'
    (e.g. gffread) would otherwise cut the value there; gtf_unquote reads it back.
    """
    if ";" in value:
        value = value.replace(";", "%3B")
    if '"' in value or "\\" in value:
        value = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{value}"'

def gtf_unquote(value: str) -> str:
    """Undoes gtf_quote escaping on the text between the quotes ('%3B' is read as ';')."""
    if "%" in value:
        value = _GTF_SEMICOLON.sub(";", value)
    return _GTF_UNESCAPE.sub(r"\1", value) if "\\" in value else value

def format_gff3_values(values) -> str:
    """Encodes each value and joins them as a multi-value GFF3 attribute."""
    return ",".join(gff3_escape(v) for v in values)

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
    phase: int | None = None
    line: int = 0
    joined_fragments: int = 0     # pieces with no '=' appended to the free-text value before them
    dropped_fragments: int = 0    # pieces with no '=' that could not be appended (nor used as a bare ID)
    parent_key: str = "Parent"    # spelling of the Parent attribute as read (e.g. 'Derives_from' in Araport11 protein rows)

def parse_gff_parts(parts) -> GffEntry:

    # Interning high-frequency strings
    source = sys.intern(parts[1])
    feature = sys.intern(parts[2])
    strand = sys.intern(parts[6])
    score = sys.intern(parts[5])

    raw_phase = parts[7]
    phase = int(raw_phase) if (raw_phase in ("0", "1", "2") and feature in default_features["CDS"]) else None

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

    if not pseudogene:
        if attr_dict.get("pseudogene") == "True" or attr_dict.get("pseudo") == "True":
            pseudogene = True

    entry_id = attr_dict.pop(INT_ID, "")
    parents = attr_dict.pop(INT_PARENT, [])
    joined_fragments, dropped_fragments = attr_dict.pop(FRAGMENTS, (0, 0))
    parent_key = attr_dict.pop(PARENT_KEY, INT_PARENT)

    if not entry_id and parents:
        entry_id = f"{parents[0]}_{feature}_{actual_start}_{actual_end}"

    return GffEntry(
        ch=ch,
        source=source,
        feature=feature,
        start=actual_start,
        end=actual_end,
        score=score,
        strand=strand,
        attributes=attr_dict,
        id=entry_id,
        parents=parents,
        pseudogene=pseudogene,
        transposable=transposable,
        decreasing_coordinates=decreasing_coordinates,
        phase=phase,
        joined_fragments=joined_fragments,
        dropped_fragments=dropped_fragments,
        parent_key=parent_key,
    )

_ID_SEPARATORS = re.compile(r"[:._|-]")

def _continues_id(piece: str, feature_id: str) -> bool:
    """
    Whether a piece with no '=' right after Name is the rest of that Name, cut at an unescaped ';':
    it has no whitespace, and its part from the first separator on holds a letter and appears in the
    line's ID as a whole token. E.g. '3:CDS:1' after 'Name=CYCA2' on 'ID=AT1G15570:CDS:1' (or on
    'ID=AT1G16780:CDS:14.2' for '2:CDS:14'), but not '8' after 'Name=AT1G20860' (left of a lost Note).
    """
    if not feature_id or any(c.isspace() for c in piece):
        return False
    m = _ID_SEPARATORS.search(piece)
    if m is None:
        return False
    tail = piece[m.start():]
    if not any(c.isalpha() for c in tail):
        return False
    return re.search(re.escape(tail) + r"(?=$|[:._|-])", feature_id) is not None

def parse_gff_attributes(attributes, gene:bool=False, transcript:bool=False):
    """
    Parses a GFF3 attribute column. Name, Symbol and Alias (under any of their MAPPED_ATTRIBUTES
    spellings) are split on ',' and percent-decoded, and repeats of one are merged into its first
    occurrence; their spelling is recorded under SOURCE_KEYS when it differs from DEFAULT_WRITE_KEYS
    so it can be written back as read. Nothing is dropped: repeated
    keys, and gene_id/transcript_id-like keys not used as ID or Parent, are kept raw under numbered
    EXTRA_ATTRIBUTES keys. Dict order is column-9 order, so attributes can be written back in it.

    A piece with no '=' after a free-text attribute is text that followed an unescaped ';' in its
    value (e.g. Araport11 'Note=SNOR111; snoRNA; snoRNA'), so it is appended back to that value as
    '%3B' + piece; after Symbol or Alias it is appended to their last value as ';' + piece (decoded
    values, encoded on export), e.g. 'symbol=SULTR1;2'. Text after Name is usually a description that
    lost its key (e.g. Araport11 'Name=AT1G01760;RNA binding'), so it is only appended back to the Name
    when it is the piece right after it and visibly continues it like the line's ID does (see
    _continues_id): Araport11 'ID=AT1G15570:CDS:1;Name=CYCA2;3:CDS:1' reads Name 'CYCA2;3:CDS:1'.
    A bare first piece on a gene or transcript line is still read as its ID (AUGUSTUS).
    Pieces that are none of these are dropped. Both are counted under FRAGMENTS as (joined, dropped).
    """

    if not attributes or attributes == ".":
        return {}

    parsed = {}
    source_keys = {}
    fallback_id = None
    fallback_parent = None
    bare_id = None
    # key of the last plain 'key=value' read (an attribute or an EXTRA_ATTRIBUTES slot), which a piece with no '=' continues
    free_text_key = None
    # (values, index) of the last Symbol/Alias value read, which a piece with no '=' continues (e.g. Araport11 'symbol=SULTR1;2')
    continued_value = None
    # (values, index) of the Name value just read, and the piece with no '=' right after it, if any
    name_value = None
    name_piece = None
    joined = dropped = 0

    def keep_extra(raw):
        slot = f"{EXTRA_ATTRIBUTES}{len(parsed)}"
        parsed[slot] = raw
        return slot

    for p in attributes.split(";"):

        key, sep, val = p.partition("=")
        if not sep: # Means '=' was not found
            token = p.strip()
            if not token:
                continue
            if free_text_key is not None:
                parsed[free_text_key] += "%3B" + p.rstrip()
                joined += 1
            elif continued_value is not None:
                values, i = continued_value
                values[i] += ";" + gff3_unescape(p.rstrip())  # decoded, like the value it continues
                joined += 1
            elif name_value is not None:
                # decided once the ID is known (see _continues_id below)
                name_piece = (*name_value, token)
                name_value = None
            # AUGUSTUS writes the ID of gene and transcript lines as a bare value ('g2', 'g2.t1;Parent=g2')
            elif (gene or transcript) and bare_id is None and not any(c.isspace() for c in token):
                bare_id = token
            else:
                dropped += 1
            continue

        free_text_key = continued_value = name_value = None
        key = key.strip()
        val = val.strip()

        lower_key = key.lower()

        if lower_key == "id":
            if INT_ID not in parsed:
                parsed[INT_ID] = val
            continue
            
        # gene_id/transcript_id only stand in for a missing ID on the right feature type
        if (gene and lower_key in ALT_ID_gene) or (transcript and lower_key in ALT_ID_transcript):
            if fallback_id is None:
                fallback_id = (key, val, keep_extra(None)) # slot filled in below if not used as ID
            else:
                keep_extra(f"{key}={val}")
            continue

        if lower_key in MAPPED_ATTRIBUTES:
            target_key = MAPPED_ATTRIBUTES[lower_key]
            values = [x.strip() for x in val.split(",") if x.strip()]
            if target_key in DECODED_ATTRIBUTES:
                values = [gff3_unescape(x) for x in values]
                if target_key in parsed:
                    # e.g. 'Name=a;Name=b' or 'gene_synonym=A1;Alias=A2': one attribute, so that
                    # replacing or removing it leaves no stale copy; first key spelling and position
                    merged = parsed[target_key]
                    merged.extend(v for v in values if v not in merged)
                    if values and target_key in CONTINUED_ATTRIBUTES:
                        continued_value = (merged, merged.index(values[-1]))
                    elif values and target_key == INT_NAME:
                        name_value = (merged, merged.index(values[-1]))
                    continue
                if key != DEFAULT_WRITE_KEYS[target_key]:
                    source_keys[target_key] = sys.intern(key)
                if values and target_key in CONTINUED_ATTRIBUTES:
                    continued_value = (values, len(values) - 1)
                elif values and target_key == INT_NAME:
                    name_value = (values, len(values) - 1)
            elif target_key in parsed:
                keep_extra(f"{key}={val}")
                continue
            elif target_key == INT_PARENT and key != INT_PARENT:
                parsed[PARENT_KEY] = sys.intern(key)
            parsed[target_key] = values
            continue
            
        # e.g. gene_id as Parent, used only when there is no Parent attribute
        if lower_key in ID_OR_PARENT_KEYS:
            if fallback_parent is None:
                fallback_parent = (key, val, keep_extra(None))
            elif not (gene or transcript) and lower_key in ALT_ID_transcript and fallback_parent[0].lower() in ALT_ID_gene:
                # a subfeature belongs to its transcript, so transcript_id wins over gene_id whatever their order
                old_key, old_val, old_slot = fallback_parent
                parsed[old_slot] = f"{old_key}={old_val}"
                fallback_parent = (key, val, keep_extra(None))
            else:
                keep_extra(f"{key}={val}")
            continue

        if key in parsed:
            free_text_key = keep_extra(f"{key}={val}")
        else:
            free_text_key = sys.intern(key)
            parsed[free_text_key] = val

    if fallback_id is not None:
        key, val, slot = fallback_id
        if INT_ID not in parsed:
            del parsed[slot]
            parsed[INT_ID] = val
        else:
            parsed[slot] = f"{key}={val}"

    if bare_id is not None:
        if INT_ID not in parsed:
            parsed[INT_ID] = bare_id
        else:
            dropped += 1

    if name_piece is not None:
        values, i, token = name_piece
        if _continues_id(token, parsed.get(INT_ID, "")):
            values[i] += ";" + gff3_unescape(token)
            joined += 1
        elif (gene or transcript) and INT_ID not in parsed and not any(c.isspace() for c in token):
            parsed[INT_ID] = token  # a bare ID, as before Name was continued
        else:
            dropped += 1

    if fallback_parent is not None:
        key, val, slot = fallback_parent
        if INT_PARENT not in parsed:
            del parsed[slot]
            parsed[INT_PARENT] = [x.strip() for x in val.split(",") if x.strip()]
        else:
            parsed[slot] = f"{key}={val}"

    if source_keys:
        parsed[SOURCE_KEYS] = source_keys

    if joined or dropped:
        parsed[FRAGMENTS] = (joined, dropped)

    return parsed

def parse_gtf_attributes(attr_string, default_key="id"):
    """
    Parses GTF attribute strings into a dictionary.
    Handles:
      1. Standard GTF: key "value"
      2. Unquoted GTF: key value
      3. Bare/Keyless values: "value" or value (stored under default_key)
    """
    attr_string = attr_string.strip()
    if not attr_string:
        return {}

    attributes = {}

    # Quoted values may hold backslash-escaped quotes (see gtf_quote)
    pattern = r'(?:(\w+)\s+"((?:[^"\\]|\\.)*)")|(?:(\w+)\s+([^\s;]+))|(?:"((?:[^"\\]|\\.)*)"|([^\s;]+))'

    for match in re.finditer(pattern, attr_string):
        k1, v1, k2, v2, v3_quoted, v3_unquoted = match.groups()

        if k1 is not None:
            attributes[k1] = gtf_unquote(v1)
        elif k2 is not None:
            attributes[k2] = v2
        else:
            bare_val = gtf_unquote(v3_quoted) if v3_quoted is not None else v3_unquoted
            attributes[default_key] = bare_val

    return attributes

def format_gff3_attributes(attrs, feature_type):
    """
    Formats a dictionary of attributes into a GFF3-compliant string.
    """

    special_keys = {'gene_id', 'transcript_id', 'exon_id', 'exon_number'}

    # IDs are written as read, but a ';' (e.g. a GTF '%3B', see gtf_quote) would end the GFF3 attribute
    if any(";" in attrs.get(k, "") for k in ('gene_id', 'transcript_id')):
        attrs = dict(attrs)
        for k in ('gene_id', 'transcript_id'):
            if k in attrs:
                attrs[k] = attrs[k].replace(";", "%3B")

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
        transcript_id = attrs.get('transcript_id')
        if transcript_id:
            gff3_attrs.append(f"Parent={transcript_id}")
            feature_id = f"{transcript_id}"

            if 'exon_number' in attrs:
                feature_id += f"_e{attrs['exon_number']}"
                gff3_attrs.append(f"ID={feature_id}")
            elif feature_type in default_features["CDS"]:
                gff3_attrs.append(f"ID={feature_id}_CDS1")

    # gene_name names the gene and transcript_name the transcript, so each becomes the Name (what
    # 'aegis symbols' sets and browsers label) of its own feature only; other lines repeat them
    if feature_type in default_features["gene"]:
        name = attrs.get('gene_name')
    elif feature_type in default_features["transcript"]:
        name = attrs.get('transcript_name')
    else:
        name = None
    # GTF values are unencoded (quoted), so everything but IDs is percent-encoded for GFF3
    if name:
        gff3_attrs.append(f"Name={gff3_escape(name)}")

    for key, value in attrs.items():
        if key not in special_keys and key not in ('gene_name', 'transcript_name'):
            gff3_attrs.append(f"{key}={gff3_escape(value)}")
            
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
    # Subfeatures naming neither gene nor transcript are written right after the gene/transcript line they
    # follow in the input, so that the file order is kept for Annotation(link_by_file_order=True)
    anchor = None
    unnamed_after: dict[tuple[str, int], list] = {}

    # PASS 1 — read file
    with open_file(gtf_file, 'r', encoding=encoding) as infile:
        for line in infile:
            if line.startswith('##FASTA') or line.startswith('>'):
                break
            if line.startswith('#'):
                if not line.startswith('##'):
                    comments.append(line)
                continue

            parts = line.rstrip('\r\n').strip(' ').split('\t')
            if len(parts) == 8:
                parts.append('')  # empty attribute column without its tab
            if len(parts) == 9 and not parts[0].strip():
                parts = []  # no seqid
            if len(parts) != 9:
                if not quiet:
                    sys.stderr.write(f"Warning: Skipping malformed line: {line.strip()}\n")
                continue

            seqname, source, feature, start, end, score, strand, frame, attr_string = parts

            if feature in default_features["gene"]:
                fallback_key = "gene_id"
            elif feature in default_features["transcript"]:
                fallback_key = "transcript_id"
            else:
                fallback_key = "id"

            attributes = parse_gtf_attributes(attr_string, default_key=fallback_key)

            if len(attributes) == 1 and fallback_key in attributes:
                bare_val = attributes[fallback_key]
                if feature in default_features["transcript"]:
                    # e.g. AUGUSTUS/TSEBRA 'anno1.g20064.t1': the transcript's own ID; its gene_id is taken
                    # from the subfeatures that name it (below), not made equal to the transcript ID
                    attributes = {'transcript_id': bare_val}
                else:
                    attributes['gene_id'] = bare_val
                    if feature not in default_features["gene"]:
                        attributes['transcript_id'] = bare_val
                    
            start_int = int(start)
            end_int = int(end)

            gene_id = attributes.get('gene_id', '')
            transcript_id = attributes.get('transcript_id', '')

            if feature in default_features["gene"]:
                if gene_id and gene_id not in seen_genes:
                    explicit_gene_lines.append((parts, attributes))
                    seen_genes.add(gene_id)
                    anchor = ('gene', len(explicit_gene_lines) - 1)
                continue

            if feature in default_features["transcript"]:
                if transcript_id and transcript_id not in seen_transcripts:
                    explicit_transcript_lines.append((parts, attributes))
                    seen_transcripts.add(transcript_id)
                    anchor = ('transcript', len(explicit_transcript_lines) - 1)
                continue

            if not gene_id and not transcript_id and anchor is not None:
                unnamed_after.setdefault(anchor, []).append((parts, attributes, feature))
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
                        'transcript_name': attributes.get('transcript_name', ''),
                    }
                else:
                    t = inferred_transcripts[transcript_id]
                    if start_int < t['start']:
                        t['start'] = start_int
                    if end_int > t['end']:
                        t['end'] = end_int

                if gene_id and gene_id in inferred_genes:
                    inferred_genes[gene_id]['transcript_ids'][transcript_id] = True

    # Subfeatures naming only their gene (e.g. JGI CDS lines with gene_id and protein_id but no transcript_id)
    # belong to that gene's transcript when it has exactly one
    gene_transcripts: dict[str, set] = {}
    for parts, attributes in explicit_transcript_lines:
        if attributes.get('gene_id'):
            gene_transcripts.setdefault(attributes['gene_id'], set()).add(attributes['transcript_id'])
    for t_id, t in inferred_transcripts.items():
        if t['gene_id']:
            gene_transcripts.setdefault(t['gene_id'], set()).add(t_id)
    for parts, attributes, feature in subfeature_lines:
        if attributes.get('transcript_id') or not attributes.get('gene_id'):
            continue
        candidates = gene_transcripts.get(attributes['gene_id'], set())
        if len(candidates) == 1:
            t_id = next(iter(candidates))
            attributes['transcript_id'] = t_id
            if t_id in inferred_transcripts:
                t = inferred_transcripts[t_id]
                t['start'] = min(t['start'], int(parts[3]))
                t['end'] = max(t['end'], int(parts[4]))

    # PASS 2 — write GFF3
    with open_file(gff3_file, 'w', encoding=encoding) as outfile:
        outfile.write("##gff-version 3\n")

        for c in comments:
            outfile.write(f"#{c}")

        def write_unnamed(key):
            for parts, attributes, feature in unnamed_after.get(key, ()):
                outfile.write("\t".join(parts[:8] + [format_gff3_attributes(attributes, feature)]) + '\n')

        for i, (parts, attributes) in enumerate(explicit_gene_lines):
            gene_attrs = {'gene_id': attributes['gene_id']}
            if 'gene_name' in attributes:
                gene_attrs['gene_name'] = attributes['gene_name']
            if 'gene_biotype' in attributes:
                gene_attrs['gene_biotype'] = attributes['gene_biotype']
            gene_attr_str = format_gff3_attributes(gene_attrs, 'gene')
            outfile.write("\t".join([parts[0], parts[1], 'gene', parts[3], parts[4],
                                     parts[5], parts[6], parts[7], gene_attr_str]) + '\n')
            write_unnamed(('gene', i))

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

        for i, (parts, attributes) in enumerate(explicit_transcript_lines):
            tx_attrs = {k: v for k, v in attributes.items()
                        if 'transcript' in k or 'gene' in k}
            if 'gene_id' not in tx_attrs:
                inferred_gene = inferred_transcripts.get(attributes['transcript_id'], {}).get('gene_id')
                if inferred_gene:
                    tx_attrs['gene_id'] = inferred_gene
            tx_feature_type = 'transcript'
            if 'transcript_biotype' in attributes:
                if 'RNA' in attributes['transcript_biotype']:
                    tx_feature_type = attributes['transcript_biotype']
            tx_attr_str = format_gff3_attributes(tx_attrs, tx_feature_type)
            outfile.write("\t".join([parts[0], parts[1], tx_feature_type,
                                     parts[3], parts[4], parts[5], parts[6],
                                     parts[7], tx_attr_str]) + '\n')
            write_unnamed(('transcript', i))

        for t_id, t in inferred_transcripts.items():
            if t_id in seen_transcripts:
                continue
            tx_attrs = {'transcript_id': t_id}
            if t['gene_id']:
                tx_attrs['gene_id'] = t['gene_id']
            if t['transcript_name']:
                tx_attrs['transcript_name'] = t['transcript_name']
            tx_attr_str = format_gff3_attributes(tx_attrs, 'mRNA')
            outfile.write("\t".join([t['seqname'], t['source'], 'mRNA',
                                     str(t['start']), str(t['end']), '.',
                                     t['strand'], '.', tx_attr_str]) + '\n')
            seen_transcripts.add(t_id)

        subfeature_counters: dict[tuple[str, str], int] = {}

        for parts, attributes, feature in subfeature_lines:
            transcript_id = attributes.get('transcript_id', '')

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

            gff3_attr_string = ";".join(gff3_attrs)
            outfile.write("\t".join([parts[0], parts[1], feature,
                                     parts[3], parts[4], parts[5], parts[6],
                                     parts[7], gff3_attr_string]) + '\n')

    if not quiet:
        print(f"Successfully converted {gtf_file} to a gff format")

def detect_file_format(file_path, encoding, lines_to_check=20):
    """
    Detects if a file is GTF or GFF3 format.
    
    Order of operations:
      1. Inspect content (headers & column 9 attributes).
      2. Fallback to file extension (.gtf vs .gff/.gff3), handling .gz/.bz2.
      3. Default to 'gff3' if still unknown.
    """
    try:
        with open_file(file_path, 'r', encoding=encoding) as f:
            i = 0
            bare_values = False
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("##FASTA") or line.startswith(">"):
                    break

                if line.startswith("##gff-version 3"):
                    return 'gff3'

                if line.startswith('#'):
                    continue

                i += 1
                parts = line.split('\t')

                if len(parts) == 9:
                    attributes = parts[8].strip()

                    if re.search(r'\b\w+=', attributes):
                        return 'gff3'

                    if re.search(r'\b\w+\s+(?:"[^"]+"|[^\s;]+)', attributes):
                        return 'gtf'

                    # A bare value (e.g. AUGUSTUS gene lines 'g2') fits both formats: let the next lines decide
                    if re.match(r'^[\w\.\-]+$', attributes):
                        bare_values = True

                if i >= lines_to_check:
                    break

            if bare_values:
                return 'gtf'

    except Exception as e:
        sys.stderr.write(f"Error reading file for format detection: {e}\n")
        sys.exit(1)

    clean_path = str(file_path).lower()

    for comp_ext in ['.gz', '.bz2', '.xz']:
        if clean_path.endswith(comp_ext):
            clean_path = clean_path[:-len(comp_ext)]

    ext = os.path.splitext(clean_path)[1]

    if ext == '.gtf':
        return 'gtf'
    elif ext in ['.gff', '.gff3']:
        return 'gff3'

    return 'gff3'