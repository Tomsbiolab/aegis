import sys

def parse_gff_line(line):
    parts = line.strip().split("\t")

    # Interning high-frequency strings

    source = sys.intern(parts[1])
    feature = sys.intern(parts[2])
    strand = sys.intern(parts[6])
    phase = sys.intern(parts[7])

    if feature == "nucleotide_to_protein_match":
        ch = sys.intern(parts[0].split(":")[0])
    
    else:
        ch = sys.intern(parts[0])

    if "pseudo" in feature:
        pseudogene = True
    else:
        pseudogene = False

    start = int(parts[3])
    end = int(parts[4])

    if start > end:
        decreasing_coordinates = True
        actual_start, actual_end = (int(parts[4]), int(parts[3]))
    else:
        decreasing_coordinates = False
        actual_start, actual_end = (start, end)

    if "transposable" in feature or "transposon" in feature:
        transposable = True
    else:
        transposable = False

    attr_dict = parse_gff_attributes(parts[8])
    
    if not transposable:
        if attr_dict.get("transposable") == "True" or attr_dict.get("transposon") == "True":
            transposable = True

    entry = {
        "ch": ch,
        "source": source,
        "feature": feature,
        "start": actual_start,
        "end": actual_end,
        "score": parts[5],
        "strand": strand,
        "phase": phase,
        "attributes": attr_dict,
        "id": attr_dict.get("id", ""),
        "parents": attr_dict.get("parent", []),
        "pseudogene": pseudogene,
        "transposable": transposable,
        "decreasing_coordinates": decreasing_coordinates
    }


    return entry

def parse_gff_attributes(attributes):

    if not attributes or attributes == ".":
        return {}

    parsed = {}
    
    for p in attributes.split(";"):
        p = p.strip()
        if not p: 
            continue

        if "=" in p:
            key, val = p.split("=", 1)

            key = key.strip().lower()
            val = val.strip()

            if key in ["parent", "parents", "derives_from"]:
                parent_key = sys.intern("parent")
                parsed[parent_key] = [x.strip() for x in val.split(",") if x.strip()]
            else:
                key = sys.intern(key)
                parsed[key] = val

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
    """
    with open(gtf_file, 'r', encoding=encoding) as infile, open(gff3_file, 'w', encoding=encoding) as outfile:

        outfile.write("##gff-version 3\n")

        seen_genes = set()
        seen_transcripts = set()

        for line in infile:
            if line.startswith('#'):
                if not line.startswith('##'):
                    outfile.write(f"#{line}")
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                if not quiet:
                    sys.stderr.write(f"Warning: Skipping malformed line: {line.strip()}\n")
                continue

            seqname, source, feature, start, end, score, strand, frame, attr_string = parts
            attributes = parse_gtf_attributes(attr_string)

            if 'gene_id' in attributes and attributes['gene_id'] not in seen_genes and feature in default_features["gene"]:
                gene_attrs = {'gene_id': attributes['gene_id']}
                if 'gene_name' in attributes:
                    gene_attrs['gene_name'] = attributes['gene_name']
                if 'gene_biotype' in attributes:
                    gene_attrs['gene_biotype'] = attributes['gene_biotype']
                gene_attr_str = format_gff3_attributes(gene_attrs, 'gene')
                gene_line = "\t".join([seqname, source, 'gene', start, end, score, strand, frame, gene_attr_str])
                outfile.write(gene_line + '\n')
                seen_genes.add(attributes['gene_id'])

            if 'transcript_id' in attributes and attributes['transcript_id'] not in seen_transcripts and feature in default_features["transcript"]:
                tx_attrs = {k: v for k, v in attributes.items() if 'transcript' in k or 'gene' in k}
                tx_feature_type = 'transcript'
                if 'transcript_biotype' in attributes:
                    if 'RNA' in attributes['transcript_biotype']:
                        tx_feature_type = attributes['transcript_biotype']

                tx_attr_str = format_gff3_attributes(tx_attrs, tx_feature_type)
                tx_line = "\t".join([seqname, source, tx_feature_type, start, end, score, strand, frame, tx_attr_str])
                outfile.write(tx_line + '\n')
                seen_transcripts.add(attributes['transcript_id'])

            if feature not in default_features["transcript"] and feature not in default_features["gene"]:
                continue

            gff3_attr_string = format_gff3_attributes(attributes, feature)

            gff3_line = "\t".join([seqname, source, feature, start, end, score, strand, frame, gff3_attr_string])
            outfile.write(gff3_line + '\n')

    if not quiet:
        print(f"Successfully converted {gtf_file} to a gff format")

def detect_file_format(file_path, encoding, lines_to_check=20):
    """
    Detects if a file is likely GTF or GFF3 format.
    """
    try:
        with open(file_path, 'r', encoding=encoding) as f:

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