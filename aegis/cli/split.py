import typer
import os
import re
import copy
import warnings

from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from ..utils.misc import open_file

app = typer.Typer(add_completion=False, no_args_is_help=True)


def parse_split_specs(values: list[str] | str | None) -> list[tuple[str, str]]:
    """
    Parse split specifications into a list of (pattern, label) tuples.

    Supports:
    - Multiple -s options: -s "A," -s "D,"
    - Semicolon delimited: -s "A,; D,"
    - Pipe delimited: -s "A, | D,"
    - Escaped commas: -s "A\\,, D\\,"
    - Double commas: -s "A,, D," or -s "A,,D,"
    - Explicit label mapping: -s "A,:A; D,:D" or -s "chromosome 2A,:A; chromosome 1D,:D"
    - Standard comma separated: -s "A, D" or -s "hap1, hap2"
    """
    if not values:
        return []
    if isinstance(values, str):
        values = [values]

    raw_items = []
    for val in values:
        val = val.strip()
        if not val:
            continue

        if ";" in val:
            items = val.split(";")
        elif "|" in val:
            items = val.split("|")
        elif "/" in val:
            items = val.split("/")
        elif "\\," in val:
            placeholder = "\x00ESCAPED_COMMA\x00"
            escaped = val.replace("\\,", placeholder)
            items = [item.replace(placeholder, ",") for item in escaped.split(",")]
        elif re.search(r",\s*,", val):
            parts = re.split(r",\s*,", val)
            items = [p + "," if not p.endswith(",") else p for p in parts]
        elif len(values) > 1:
            if val.endswith(",") and val.count(",") == 1:
                items = [val]
            elif "," in val:
                items = val.split(",")
            else:
                items = [val]
        else:
            if val.endswith(",") and val.count(",") == 1:
                items = [val]
            else:
                items = val.split(",")

        for item in items:
            item = item.strip().strip("'\"")
            if item:
                raw_items.append(item)

    specs = []
    for item in raw_items:
        if ":" in item:
            pat, lab = item.split(":", 1)
            pat, lab = pat.strip().strip("'\""), lab.strip().strip("'\"")
        elif "=" in item:
            pat, lab = item.split("=", 1)
            pat, lab = pat.strip().strip("'\""), lab.strip().strip("'\"")
        else:
            pat = item.strip().strip("'\"")
            # Auto-clean label for filenames (strip punctuation like commas, colons, semicolons)
            lab = re.sub(r'^[_\-\.:,;\s]+|[_\-\.:,;\s]+$', '', pat)
            if not lab:
                lab = pat
        specs.append((pat, lab))
    return specs


def detect_file_type(filepath: str) -> str:
    """Detect whether a file is FASTA or GFF/GTF annotation."""
    if not filepath or not os.path.exists(filepath):
        return "unknown"
    fasta_exts = (".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fasta.gz", ".fna.gz")
    annot_exts = (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")
    lower = filepath.lower()
    if lower.endswith(fasta_exts):
        return "fasta"
    if lower.endswith(annot_exts):
        return "annotation"
    try:
        with open_file(filepath, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("##") or "\t" in line:
                    return "annotation"
                break
    except Exception:
        pass
    return "unknown"


def classify_feature(
    name: str,
    description: str,
    split_specs: list[tuple[str, str]] | list[str] | None = None,
    regex_pattern: str | None = None,
    match_mode: str = "smart",
    ignore_case: bool = False,
    split_map: dict[str, str] | None = None,
    allowed_tags: set[str] | None = None,
    split_tags: list[tuple[str, str]] | list[str] | None = None,
) -> tuple[Optional[str], Optional[str]]:
    """
    Classify a genomic feature (scaffold/contig/chromosome) into a split tag label.

    Returns:
        tuple (assigned_tag, warning_or_reason)
    """
    if split_specs is None and split_tags is not None:
        split_specs = split_tags
    # 1. Direct dictionary map check
    if split_map:
        if name in split_map:
            return split_map[name], None

    flags = re.IGNORECASE if ignore_case else 0
    full_text = f"{name} {description}" if description and description != name else name

    # Normalize split_specs to list of (pattern, label) tuples
    normalized_specs: list[tuple[str, str]] = []
    if split_specs:
        for item in split_specs:
            if isinstance(item, tuple):
                normalized_specs.append(item)
            elif isinstance(item, str):
                lab = re.sub(r'^[_\-\.:,;\s]+|[_\-\.:,;\s]+$', '', item)
                normalized_specs.append((item, lab if lab else item))

    # 2. Custom Regex pattern
    if regex_pattern:
        try:
            match = re.search(regex_pattern, full_text, flags)
            if match:
                if match.groups():
                    captured = match.group(1).strip()
                    if allowed_tags is not None:
                        for tag in allowed_tags:
                            if (captured.lower() == tag.lower()) if ignore_case else (captured == tag):
                                return tag, None
                        return None, None
                    return captured, None
                else:
                    return normalized_specs[0][1] if normalized_specs else "matched", None
        except re.error as e:
            raise ValueError(f"Invalid regular expression pattern '{regex_pattern}': {e}")
        return None, None

    if not normalized_specs:
        return None, None

    matched_labels = []

    for pattern, label in normalized_specs:
        has_punctuation = any(c in pattern for c in ',;:.|_- ')

        if match_mode == "exact":
            p = rf"\b{re.escape(pattern)}\b" if not has_punctuation else re.escape(pattern)
            if re.search(p, full_text, flags):
                matched_labels.append(label)

        elif match_mode == "substring":
            target_text = full_text.lower() if ignore_case else full_text
            target_pat = pattern.lower() if ignore_case else pattern
            if target_pat in target_text:
                matched_labels.append(label)

        else:  # "smart" mode
            # If pattern has punctuation/spaces (e.g. 'A,', 'chromosome 2A,', '_A'):
            if has_punctuation:
                target_text = full_text.lower() if ignore_case else full_text
                target_pat = pattern.lower() if ignore_case else pattern
                if target_pat in target_text:
                    matched_labels.append(label)
            else:
                # Discrete pattern for alphanumeric tags (like 'A', 'D', 'hap1')
                # 1. Check feature name (ID)
                p_name = rf"(?:^|[_\-\.])(?:chr\d*|\d+|hap(?:lotype)?\d*|scaffold\d*|contig\d*)?{re.escape(pattern)}(?:[_\-\.]|$|\b)"
                if re.search(p_name, name, flags) or re.search(rf"\b\d*{re.escape(pattern)}\b", name, flags):
                    matched_labels.append(label)
                    continue

                # 2. Check description context (chromosome 2A, 1D, etc.)
                p_desc = rf"(?:chromosome\s*\d*|chr\s*\d*|scaffold\s*\d*|contig\s*\d*|hap(?:lotype)?\s*|[_\-\.]|\b\d+)\s*{re.escape(pattern)}(?:[,\s;:\._\-]|$|\b)"
                if re.search(p_desc, description, flags) or re.search(rf"\b\d*{re.escape(pattern)}\b", description, flags):
                    matched_labels.append(label)
                    continue

                # 3. Only if length > 1, allow substring fallback
                if len(pattern) > 1:
                    target_text = full_text.lower() if ignore_case else full_text
                    target_pat = pattern.lower() if ignore_case else pattern
                    if target_pat in target_text:
                        matched_labels.append(label)

    unique_matches = list(dict.fromkeys(matched_labels))
    if len(unique_matches) == 1:
        return unique_matches[0], None
    elif len(unique_matches) > 1:
        return None, f"Scaffold '{name}' matched multiple tags: {unique_matches}"

    return None, None


def resolve_split_filename(template: str, name: str, tag: str, default_ext: str) -> str:
    """Format output filename using template, name, tag and extension."""
    tag_str = tag
    if tag == "other":
        if "_split{tag}" in template:
            filename = template.replace("_split{tag}", "_split_other")
        elif "split{tag}" in template:
            filename = template.replace("split{tag}", "split_other")
        elif "{tag}" in template:
            filename = template.replace("{tag}", "_other")
        else:
            filename = f"{template}_other"
    else:
        filename = template.replace("{tag}", tag_str)

    filename = filename.replace("{annotation-name}", name).replace("{genome-name}", name)

    if not os.path.splitext(filename)[1]:
        filename = f"{filename}{default_ext}"
    return filename


@app.command()
def main(
    file1: Annotated[str, typer.Argument(
        help="Input annotation (GFF/GTF) or genome (FASTA) file."
    )] = "",
    file2: Annotated[str, typer.Argument(
        help="Optional second input file (genome FASTA or annotation GFF/GTF)."
    )] = "",
    annotation_file: Annotated[str, typer.Option(
        "-a", "--annotation", help="Path to the input annotation GFF/GTF file. Overrides positional argument if provided."
    )] = "",
    genome_file: Annotated[str, typer.Option(
        "-g", "--genome", help="Path to the input genome FASTA file. Overrides positional argument if provided."
    )] = "",
    split_by: Annotated[Optional[list[str]], typer.Option(
        "-s", "--split-by", help="Tag(s) or pattern(s) to split by (e.g. 'A,B', 'A,;D,', 'A,,D,', or multiple -s A, -s D,).",
    )] = None,
    regex: Annotated[str, typer.Option(
        "-r", "--regex", help="Optional regular expression pattern to extract the split tag from sequence ID or description (e.g. 'chromosome \\d+([A-Z])')."
    )] = "",
    match_mode: Annotated[str, typer.Option(
        "-m", "--match-mode", help="Matching mode for split-by tags: 'smart' (default; checks chromosome/haplotype context and word boundaries), 'substring' (simple substring match), or 'exact' (exact word match)."
    )] = "smart",
    ignore_case: Annotated[bool, typer.Option(
        "-i", "--ignore-case", help="Case-insensitive matching. By default, matching is case-sensitive."
    )] = False,
    split_map: Annotated[str, typer.Option(
        "--split-map", help="Path to a TSV file for explicit scaffold mapping. Format: 'scaffold_id<tab>split_group' per line."
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/split/",
    output_annot_file: Annotated[str, typer.Option(
        "-oa", "--output-annot-file", help="Template for output annotation filename. Use '{tag}' or '{annotation-name}'. [default: '{annotation-name}_split{tag}.gff3']"
    )] = "{annotation-name}_split{tag}.gff3",
    output_genome_file: Annotated[str, typer.Option(
        "-og", "--output-genome-file", help="Template for output genome filename. Use '{tag}' or '{genome-name}'. [default: '{genome-name}_split{tag}.fasta']"
    )] = "{genome-name}_split{tag}.fasta",
    keep_description: Annotated[bool, typer.Option(
        "--keep-description/--no-keep-description", help="Preserve full FASTA header descriptions in output genome files."
    )] = False,
    write_empty_other: Annotated[bool, typer.Option(
        "--write-empty-other", help="Write '_split_other' files even if no scaffolds are unassigned."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Split an annotation (GFF/GTF) and/or genome assembly (FASTA) into distinct files
    based on genomic feature names, descriptions, haplotypes, or regex patterns.
    """
    # 1. Resolve input files (positional vs options)
    pos_files = [f for f in (file1, file2) if f]

    if not annotation_file and not genome_file:
        if not pos_files:
            typer.secho("Error: At least one input file (annotation GFF/GTF or genome FASTA) must be provided.", fg=typer.colors.RED)
            raise typer.BadParameter("Missing input file(s). Provide an annotation or genome file.")

        for f in pos_files:
            ftype = detect_file_type(f)
            if ftype == "fasta":
                if not genome_file:
                    genome_file = f
                elif not annotation_file:
                    annotation_file = f
            elif ftype == "annotation":
                if not annotation_file:
                    annotation_file = f
                elif not genome_file:
                    genome_file = f
            else:
                if not annotation_file:
                    annotation_file = f
                elif not genome_file:
                    genome_file = f
    else:
        # Fill missing file from positional arguments if available
        for f in pos_files:
            if f == annotation_file or f == genome_file:
                continue
            ftype = detect_file_type(f)
            if ftype == "fasta" and not genome_file:
                genome_file = f
            elif ftype == "annotation" and not annotation_file:
                annotation_file = f
            elif not genome_file:
                genome_file = f
            elif not annotation_file:
                annotation_file = f

    # Validate that at least one file is present
    if not annotation_file and not genome_file:
        typer.secho("Error: Neither annotation nor genome file could be determined.", fg=typer.colors.RED)
        raise typer.BadParameter("No valid input files provided.")

    if annotation_file and not os.path.exists(annotation_file):
        raise FileNotFoundError(f"Annotation file not found: {annotation_file}")
    if genome_file and not os.path.exists(genome_file):
        raise FileNotFoundError(f"Genome file not found: {genome_file}")

    # 2. Validate split criteria
    split_specs = parse_split_specs(split_by)

    tsv_map = None
    if split_map:
        if not os.path.exists(split_map):
            raise FileNotFoundError(f"Split map file not found: {split_map}")
        with open(split_map, "r", encoding="utf-8") as f:
            tsv_map = {}
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    scaff_id = parts[0].strip()
                    grp = parts[1].strip()
                    tsv_map[scaff_id] = grp
                    if grp not in [lab for _, lab in split_specs]:
                        split_specs.append((grp, grp))

    if not split_specs and not regex and not tsv_map:
        typer.secho("Error: Please specify split criteria via --split-by, --regex, or --split-map.", fg=typer.colors.RED)
        raise typer.BadParameter("Missing split criteria.")

    # 3. Derive names
    if annotation_file:
        if annotation_name == "{annotation-file}":
            annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
            if annotation_name.endswith(".gff") or annotation_name.endswith(".gtf"):
                annotation_name = os.path.splitext(annotation_name)[0]
    else:
        annotation_name = "annotation"

    if genome_file:
        if genome_name == "{genome-file}":
            genome_name = os.path.splitext(os.path.basename(genome_file))[0]
            if genome_name.endswith(".fa") or genome_name.endswith(".fasta"):
                genome_name = os.path.splitext(genome_name)[0]
    else:
        genome_name = "genome"

    os.makedirs(output_dir, exist_ok=True)

    # 4. Load objects
    genome_obj = None
    annot_obj = None

    if genome_file:
        genome_obj = Genome(name=genome_name, genome_file_path=genome_file, quiet=quiet)

    if annotation_file:
        if genome_obj:
            annot_obj = Annotation(annot_file_path=annotation_file, name=annotation_name, genome=genome_obj, quiet=quiet)
        else:
            annot_obj = Annotation(annot_file_path=annotation_file, name=annotation_name, quiet=quiet)

    # 5. Classify features (scaffolds/chromosomes)
    # Collect all feature identifiers and their descriptions
    feature_info = {}  # ft_id -> description

    if genome_obj:
        for scaffold_id, scaffold in genome_obj.scaffolds.items():
            desc = getattr(scaffold, "description", scaffold_id)
            feature_info[scaffold_id] = desc

    if annot_obj:
        for chrom in annot_obj.chrs:
            if chrom not in feature_info:
                feature_info[chrom] = chrom

    split_labels = list(dict.fromkeys([lab for _, lab in split_specs]))
    tag_to_features: dict[str, set[str]] = {t: set() for t in split_labels}
    tag_to_features["other"] = set()

    target_allowed_tags = set(split_labels) if split_labels else None

    for ft_id, desc in feature_info.items():
        tag, warn = classify_feature(
            name=ft_id,
            description=desc,
            split_specs=split_specs,
            regex_pattern=regex,
            match_mode=match_mode,
            ignore_case=ignore_case,
            split_map=tsv_map,
            allowed_tags=target_allowed_tags,
        )
        if warn and not quiet:
            warnings.warn(warn, category=UserWarning)

        if tag:
            if tag not in tag_to_features:
                tag_to_features[tag] = set()
                if tag not in split_labels:
                    split_labels.append(tag)
            tag_to_features[tag].add(ft_id)
        else:
            tag_to_features["other"].add(ft_id)

    # 6. Execute partition export
    all_tags = list(split_labels)
    if tag_to_features["other"] or write_empty_other:
        all_tags.append("other")

    if not quiet:
        typer.echo(f"\n--- AEGIS Split Summary ---")
        if genome_file:
            typer.echo(f"Genome: {genome_file} ({len(genome_obj.scaffolds)} scaffolds)")
        if annotation_file:
            typer.echo(f"Annotation: {annotation_file} ({len(annot_obj.chrs)} chromosomes, {len(annot_obj.all_gene_ids)} genes)")
        typer.echo(f"Output directory: {output_dir}")

    for tag in all_tags:
        target_features = tag_to_features.get(tag, set())

        if not target_features and not write_empty_other:
            if not quiet:
                typer.echo(f"  Tag '{tag}': 0 features matched, skipping export.")
            continue

        num_scaffolds = len(target_features.intersection(genome_obj.scaffolds)) if genome_obj else len(target_features)
        num_genes = 0

        # Export genome subset
        if genome_obj:
            genome_fts = target_features.intersection(genome_obj.scaffolds)
            if genome_fts or write_empty_other:
                g_split = copy.copy(genome_obj)
                g_split.scaffolds = {k: v.copy() for k, v in genome_obj.scaffolds.items() if k in genome_fts}
                g_split.update()

                out_g_filename = resolve_split_filename(output_genome_file, genome_name, tag, ".fasta")
                g_split.export(output_dir=output_dir, filename=out_g_filename, keep_description=keep_description, quiet=True)

        # Export annotation subset
        if annot_obj:
            annot_fts = target_features.intersection(annot_obj.chrs)
            if annot_fts or write_empty_other:
                if annot_fts:
                    a_split = annot_obj.copy()
                    a_split.subset(chosen_features=annot_fts, gene_cap=None, no_gene_cap=True, min_genes=0, quiet=True)
                    num_genes = len(a_split.all_gene_ids)
                else:
                    a_split = annot_obj.copy()
                    a_split.chrs = {}
                    a_split.all_gene_ids = {}

                # Partition atypical and orphaned features if present
                if hasattr(a_split, "atypical_features"):
                    a_split.atypical_features = [
                        ft for ft in a_split.atypical_features
                        if getattr(ft, "chromosome", getattr(ft, "ch", None)) in annot_fts
                    ]
                if hasattr(a_split, "orphaned_features"):
                    a_split.orphaned_features = [
                        ft for ft in a_split.orphaned_features
                        if getattr(ft, "chromosome", getattr(ft, "ch", None)) in annot_fts
                    ]

                out_a_filename = resolve_split_filename(output_annot_file, annotation_name, tag, ".gff3")
                is_gtf = out_a_filename.lower().endswith(".gtf")
                if is_gtf:
                    a_split.export.gtf(output_dir=output_dir, filename=out_a_filename, subfolder=False, quiet=True)
                else:
                    a_split.export.gff(output_dir=output_dir, filename=out_a_filename, subfolder=False, skip_atypical_fts=False, quiet=True)

        if not quiet:
            label = "other (unassigned)" if tag == "other" else f"Tag '{tag}'"
            typer.echo(f"  - {label}: {num_scaffolds} scaffolds" + (f", {num_genes} genes" if annot_obj else ""))

    if not quiet:
        typer.secho("\nSplit complete!\n", fg=typer.colors.GREEN)


if __name__ == "__main__":
    app()
