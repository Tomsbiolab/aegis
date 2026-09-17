import os
import re
import typer
from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

from ..genome import Genome, Scaffold

app = typer.Typer(add_completion=False, no_args_is_help=True)


def format_human_readable(size: int | float, is_bp: bool = True) -> str:
    """Format a base-pair size or count into a human-readable string."""
    abs_size = abs(size)
    sign = "-" if size < 0 else ""
    if abs_size >= 1e9:
        return f"{sign}{abs_size / 1e9:.2f} Gb"
    if abs_size >= 1e6:
        return f"{sign}{abs_size / 1e6:.2f} Mb"
    if abs_size >= 1e3:
        return f"{sign}{abs_size / 1e3:.2f} Kb"
    suffix = " bp" if is_bp else ""
    return f"{sign}{int(abs_size)}{suffix}"


def parse_size_str(size_str: str) -> int:
    """Parse size string with optional K/M/G suffix (e.g. '135M', '1.2G', '500k')."""
    s = size_str.strip().upper()
    multipliers = {
        "KB": 1_000, "K": 1_000,
        "MB": 1_000_000, "M": 1_000_000,
        "GB": 1_000_000_000, "G": 1_000_000_000,
        "BP": 1
    }
    for suffix, mult in sorted(multipliers.items(), key=lambda x: -len(x[0])):
        if s.endswith(suffix):
            return int(float(s[:-len(suffix)].strip()) * mult)
    return int(float(s))


def format_number(val: int | float | None, human_readable: bool = False, is_terminal: bool = True, is_pct: bool = False) -> str:
    """Format a metric value for terminal or tabular output."""
    if val is None:
        return "-"
    if is_pct:
        return f"{val:.2f}%" if is_terminal else f"{val:.2f}"
    if human_readable:
        return format_human_readable(val)
    if is_terminal:
        if isinstance(val, float):
            return f"{val:,.2f}"
        return f"{val:,}"
    return str(val)


def format_diff(diff: int | float | None, human_readable: bool = False, is_terminal: bool = True, is_pct: bool = False) -> str:
    """Format a difference value with explicit sign."""
    if diff is None:
        return "-"
    if diff == 0:
        if is_pct:
            return "0.00%" if is_terminal else "0.00"
        if human_readable:
            return "0 bp"
        return "0"

    sign = "+" if diff > 0 else ""
    if is_pct:
        return f"{sign}{diff:.2f}%" if is_terminal else f"{sign}{diff:.2f}"
    if human_readable:
        return f"{sign}{format_human_readable(diff)}"
    if is_terminal:
        if isinstance(diff, float):
            return f"{sign}{diff:,.2f}"
        return f"{sign}{diff:,}"
    return f"{sign}{diff}"


def get_natural_sort_key(name: str, genomes: list[Genome]):
    """Natural sort key respecting chromosome > unknown > organelle > scaffold hierarchy."""
    # Check category across all genomes
    category = 4
    for g in genomes:
        scf = g.get_scaffold(name)
        if scf is not None:
            if scf.chromosome and not scf.organelle and not scf.unknown_chromosome:
                category = min(category, 1)
            elif scf.unknown_chromosome:
                category = min(category, 2)
            elif scf.organelle:
                category = min(category, 3)

    if category == 4 and name.lower().startswith("ch"):
        category = 1

    parts = [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', name)]
    return (category, parts)


def normalize_chr_name(name: str) -> str:
    """Normalize chromosome names for heuristic matching (e.g. 'chr01', 'Chr1', 'chromosome_1' -> '1')."""
    s = name.strip().lower()
    for prefix in ["chromosome_", "chromosome", "scaffold_", "scaffold", "contig_", "contig", "chr_", "chr"]:
        if s.startswith(prefix):
            s = s[len(prefix):]
            break
    s = s.lstrip("0")
    return s if s else name.strip().lower()


class PairedFeature:
    """Represents a unified chromosome or scaffold feature across compared genomes."""
    def __init__(self, primary_name: str):
        self.primary_name = primary_name
        self.scaffolds: dict[str, Scaffold] = {}
        self.synonyms: dict[str, str] = {}
        self.match_types: dict[str, str] = {}
        self.seq_diff: dict[str, bool] = {}


def pair_genome_features(
    genomes: list[Genome],
    ref_idx: int = 0,
    check_sequence: bool = True,
    include_all: bool = False,
    chromosomes_only: bool = False,
) -> list[PairedFeature]:
    """
    Pairs chromosomes and scaffolds across multiple genomes using hierarchical matching:
    1. Exact Name match
    2. Exact Sequence match (via MD5 / reverse-complement MD5)
    3. Unequivocal unique size match (identifiable unique sizes across candidate contigs)
    4. Normalized Name match (e.g. 'chr01' vs '1')
    5. Unique unpaired contigs
    """
    show_all = include_all and not chromosomes_only

    def is_candidate(scf: Scaffold) -> bool:
        if show_all:
            return True
        if chromosomes_only:
            return bool(scf.chromosome and not scf.organelle and not scf.unknown_chromosome)
        return bool(scf.chromosome or scf.unknown_chromosome or scf.organelle)

    ref_genome = genomes[ref_idx]
    ref_candidates = [scf for scf in ref_genome.scaffolds.values() if is_candidate(scf)]
    if not ref_candidates and not chromosomes_only:
        ref_candidates = list(ref_genome.scaffolds.values())

    features: list[PairedFeature] = []
    for ref_scf in ref_candidates:
        pf = PairedFeature(primary_name=ref_scf.name)
        pf.scaffolds[ref_genome.name] = ref_scf
        pf.synonyms[ref_genome.name] = ref_scf.name
        pf.match_types[ref_genome.name] = "anchor"
        features.append(pf)

    # Process each non-reference genome against features
    for g in genomes:
        if g == ref_genome:
            continue

        available_scaffolds = dict(g.scaffolds)
        unmatched_pfs = [pf for pf in features if g.name not in pf.scaffolds]

        # Stage 1: Exact Name Match
        for pf in list(unmatched_pfs):
            if pf.primary_name in available_scaffolds:
                scf = available_scaffolds.pop(pf.primary_name)
                pf.scaffolds[g.name] = scf
                pf.synonyms[g.name] = scf.name
                pf.match_types[g.name] = "name"
                unmatched_pfs.remove(pf)

        # Stage 2: Exact Sequence Match (if check_sequence)
        if check_sequence and unmatched_pfs and available_scaffolds:
            for pf in list(unmatched_pfs):
                anchor_scf = pf.scaffolds.get(ref_genome.name) or next(iter(pf.scaffolds.values()))
                h = anchor_scf.seq_hash
                rc_h = anchor_scf.rc_seq_hash

                exact_matches = [s for s in available_scaffolds.values() if s.seq_hash == h]
                if len(exact_matches) == 1:
                    matched = exact_matches[0]
                    available_scaffolds.pop(matched.name)
                    pf.scaffolds[g.name] = matched
                    pf.synonyms[g.name] = matched.name
                    pf.match_types[g.name] = "exact_seq"
                    pf.seq_diff[g.name] = False
                    unmatched_pfs.remove(pf)
                    continue

                rc_matches = [s for s in available_scaffolds.values() if s.seq_hash == rc_h]
                if len(rc_matches) == 1:
                    matched = rc_matches[0]
                    available_scaffolds.pop(matched.name)
                    pf.scaffolds[g.name] = matched
                    pf.synonyms[g.name] = matched.name
                    pf.match_types[g.name] = "rc_seq"
                    pf.seq_diff[g.name] = False
                    unmatched_pfs.remove(pf)
                    continue

        # Stage 3: Unequivocal Unique Size Match
        if unmatched_pfs and available_scaffolds:
            for pf in list(unmatched_pfs):
                anchor_scf = pf.scaffolds.get(ref_genome.name) or next(iter(pf.scaffolds.values()))
                target_size = anchor_scf.size

                ref_same_size = [p for p in unmatched_pfs if (p.scaffolds.get(ref_genome.name) or next(iter(p.scaffolds.values()))).size == target_size]
                k_same_size = [s for s in available_scaffolds.values() if s.size == target_size]

                if len(ref_same_size) == 1 and len(k_same_size) == 1:
                    matched = k_same_size[0]
                    available_scaffolds.pop(matched.name)
                    pf.scaffolds[g.name] = matched
                    pf.synonyms[g.name] = matched.name
                    pf.match_types[g.name] = "unique_size"
                    if check_sequence:
                        is_ident = (matched.seq_hash == anchor_scf.seq_hash or matched.seq_hash == anchor_scf.rc_seq_hash)
                        pf.seq_diff[g.name] = not is_ident
                    else:
                        pf.seq_diff[g.name] = False
                    unmatched_pfs.remove(pf)

        # Stage 4: Normalized Name Match
        if unmatched_pfs and available_scaffolds:
            for pf in list(unmatched_pfs):
                norm_ref = normalize_chr_name(pf.primary_name)
                norm_matches = [s for s in available_scaffolds.values() if normalize_chr_name(s.name) == norm_ref]
                if len(norm_matches) == 1:
                    matched = norm_matches[0]
                    available_scaffolds.pop(matched.name)
                    pf.scaffolds[g.name] = matched
                    pf.synonyms[g.name] = matched.name
                    pf.match_types[g.name] = "norm_name"
                    anchor_scf = pf.scaffolds.get(ref_genome.name) or next(iter(pf.scaffolds.values()))
                    if check_sequence:
                        pf.seq_diff[g.name] = (matched.size == anchor_scf.size and matched.seq_hash != anchor_scf.seq_hash and matched.seq_hash != anchor_scf.rc_seq_hash)
                    unmatched_pfs.remove(pf)

        # Stage 5: Remaining scaffolds qualifying on their own
        for scf_name, scf in available_scaffolds.items():
            if is_candidate(scf):
                pf = PairedFeature(primary_name=scf.name)
                pf.scaffolds[g.name] = scf
                pf.synonyms[g.name] = scf.name
                pf.match_types[g.name] = "unique"
                features.append(pf)

        # If any scaffolds in g matched a nuclear chromosome in the anchor, promote them
        promoted = False
        for pf in features:
            if g.name in pf.scaffolds:
                scf = pf.scaffolds[g.name]
                anchor = pf.scaffolds.get(ref_genome.name)
                if anchor and anchor.chromosome and not anchor.organelle and not anchor.unknown_chromosome:
                    if not scf.chromosome:
                        scf.chromosome = True
                        promoted = True
        if promoted:
            g.update()

    return features


def render_terminal_table(headers: list[str], rows: list[list[str]], summary_rows: list[list[str]]) -> str:
    """Renders a cleanly formatted ASCII table with aligned columns."""
    all_rows = rows + summary_rows
    num_cols = len(headers)

    col_widths = [0] * num_cols

    # Col 0 accounts for 2-space indentation in summary rows
    max_c0 = len(headers[0])
    for r in rows:
        if len(r) > 0:
            max_c0 = max(max_c0, len(r[0]))
    for r in summary_rows:
        if len(r) > 0:
            max_c0 = max(max_c0, len(r[0]) + 2)
    col_widths[0] = max_c0 + 2

    # Remaining columns
    for i in range(1, num_cols):
        max_ci = len(headers[i])
        for r in all_rows:
            val = r[i] if i < len(r) else ""
            max_ci = max(max_ci, len(val))
        col_widths[i] = max_ci + 2

    total_width = sum(col_widths) + (num_cols - 1)

    lines = []
    lines.append("=" * total_width)

    # Header
    header_line = headers[0].ljust(col_widths[0])
    for i in range(1, num_cols):
        header_line += " " + headers[i].rjust(col_widths[i])
    lines.append(header_line)
    lines.append("-" * total_width)

    # Feature rows
    for row in rows:
        line = row[0].ljust(col_widths[0])
        for i in range(1, num_cols):
            val = row[i] if i < len(row) else "-"
            line += " " + val.rjust(col_widths[i])
        lines.append(line)

    # Summary section
    if summary_rows:
        if rows:
            lines.append("-" * total_width)
        lines.append(" Assembly Summary:")
        for row in summary_rows:
            line = f"  {row[0]}".ljust(col_widths[0])
            for i in range(1, num_cols):
                val = row[i] if i < len(row) else "-"
                line += " " + val.rjust(col_widths[i])
            lines.append(line)

    lines.append("=" * total_width)
    return "\n".join(lines)


@app.command()
def main(
    genome_files: Annotated[List[str], typer.Argument(
        help="Path to one or more input genome FASTA file(s)."
    )],
    genome_names: Annotated[str, typer.Option(
        "-g", "--genome-names", help="Comma-separated genome names/tags (e.g. 'Col-0,Cvi-0'). Defaults to filenames."
    )] = "",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to output TSV/CSV file."
    )] = "",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Output directory if output file name is not explicitly specified."
    )] = "",
    include_all: Annotated[bool, typer.Option(
        "-a", "--all", help="Include all scaffolds and contigs in the table, not just chromosomes."
    )] = False,
    summary_only: Annotated[bool, typer.Option(
        "--summary-only", help="Report only assembly-level summary statistics without listing individual chromosomes."
    )] = False,
    contigs_only: Annotated[bool, typer.Option(
        "--contigs-only", "--scaffolds-only", help="Report only chromosome/scaffold statistics without the summary table."
    )] = False,
    chromosomes_only: Annotated[bool, typer.Option(
        "--chromosomes-only", help="Report only chromosomes in table and exclude unplaced scaffolds/contigs."
    )] = False,
    human_readable: Annotated[bool, typer.Option(
        "-H", "--human-readable", help="Display sizes in human-readable units (e.g., Kb, Mb, Gb)."
    )] = False,
    sort_by: Annotated[str, typer.Option(
        "--sort-by", help="Sort order for chromosomes: 'name' (natural sort, default), 'size' (descending size), or 'order' (FASTA order)."
    )] = "name",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Suppress terminal output (useful when exporting to file)."
    )] = False,
    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID')."
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)')."
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers."
    )] = False,
    genome_size: Annotated[str, typer.Option(
        "-s", "--genome-size", help="Estimated/expected genome size (e.g. '135M', '1.2G', or '135000000') for calculating NG50, LG50, and auNG."
    )] = "",
    reference: Annotated[bool, typer.Option(
        "-r", "--reference", help="Use first genome as reference (or specified via --ref-genome) and report relative differences (showing '= ref' for identical features/metrics)."
    )] = False,
    ref_genome: Annotated[str, typer.Option(
        "--ref-genome", help="Specify a particular genome name or 1-based index to use as reference (automatically activates reference mode)."
    )] = "",
    diff_only: Annotated[bool, typer.Option(
        "--diff-only", help="Report only features and summary statistics where genomes differ from reference (automatically activates reference mode; hides rows that are '= ref')."
    )] = False,
    no_seq: Annotated[bool, typer.Option(
        "--no-seq", help="Disable sequence-level hash matching (rely on name and unequivocal size matching only)."
    )] = False,
):
    """
    Summarise and compare chromosome sizes and assembly statistics for one or more genomes.
    """
    if not genome_files:
        typer.echo("Error: At least one genome FASTA file must be provided.", err=True)
        raise typer.Exit(code=1)

    if summary_only and contigs_only:
        typer.echo("Error: Cannot specify both --summary-only and --contigs-only.", err=True)
        raise typer.Exit(code=1)

    # Parse estimated genome size if provided
    exp_size = None
    if genome_size:
        try:
            exp_size = parse_size_str(genome_size)
        except Exception:
            typer.echo(f"Warning: Could not parse --genome-size '{genome_size}'. Ignoring.", err=True)

    # Parse names
    names = [item.strip() for item in genome_names.split(",")] if genome_names else []
    
    genomes: list[Genome] = []
    for idx, gfile in enumerate(genome_files):
        if idx < len(names) and names[idx]:
            gname = names[idx]
        else:
            gname = Path(gfile).stem
            if gname.endswith(".fasta") or gname.endswith(".fa"):
                gname = Path(gname).stem

        g = Genome(
            name=gname,
            genome_file_path=gfile,
            quiet=True,
            header_id_tag=header_id_tag if header_id_tag else None,
            header_id_regex=header_id_regex if header_id_regex else None,
            gwh=gwh,
        )
        genomes.append(g)

    is_two_genomes = len(genomes) == 2

    # Determine reference genome
    ref_idx = 0
    if ref_genome:
        if ref_genome.isdigit():
            idx = int(ref_genome) - 1
            if 0 <= idx < len(genomes):
                ref_idx = idx
            else:
                typer.echo(f"Warning: Reference index '{ref_genome}' out of range. Defaulting to 1st genome.", err=True)
        else:
            found = False
            for idx, g in enumerate(genomes):
                if g.name.lower() == ref_genome.lower():
                    ref_idx = idx
                    found = True
                    break
            if not found:
                typer.echo(f"Warning: Reference genome '{ref_genome}' not found in genomes. Defaulting to 1st genome.", err=True)

    if diff_only and len(genomes) < 2:
        typer.echo("Warning: --diff-only requires at least two genomes to compare. Ignoring.", err=True)
        diff_only = False

    is_ref_mode = bool(reference or ref_genome or diff_only)
    ref_name = genomes[ref_idx].name

    # Collect and pair features
    paired_features = pair_genome_features(
        genomes,
        ref_idx=ref_idx,
        check_sequence=not no_seq,
        include_all=include_all,
        chromosomes_only=chromosomes_only,
    )

    stats_list = [g.get_stats(estimated_genome_size=exp_size) for g in genomes]

    # Sort features
    if sort_by == "order":
        sorted_features = paired_features
    elif sort_by == "size":
        def max_feat_size(pf: PairedFeature):
            return max([s.size for s in pf.scaffolds.values()] or [0])
        sorted_features = sorted(paired_features, key=lambda pf: (-max_feat_size(pf), pf.primary_name))
    else:
        # Default: natural sort
        sorted_features = sorted(paired_features, key=lambda pf: get_natural_sort_key(pf.primary_name, genomes))

    # Prepare Headers
    unit_label = "" if human_readable else " (bp)"
    headers = ["Feature / Metric"]
    for idx, g in enumerate(genomes):
        if is_ref_mode and idx == ref_idx:
            headers.append(f"{g.name} [Ref]")
        else:
            headers.append(f"{g.name}{unit_label}")
    if is_two_genomes and not is_ref_mode:
        headers.append(f"Diff{unit_label}")

    # Build Feature Rows
    terminal_rows = []
    file_rows = []

    if not summary_only:
        for feat in sorted_features:
            ref_scf = feat.scaffolds.get(ref_name) if is_ref_mode else None
            all_same_as_ref = True

            term_row = [feat.primary_name]
            file_row = [feat.primary_name]

            if not is_ref_mode:
                # Standard Mode
                sizes = []
                for g in genomes:
                    scf = feat.scaffolds.get(g.name)
                    if scf is None:
                        sizes.append(None)
                        term_row.append("-")
                        file_row.append("-")
                    else:
                        sizes.append(scf.size)
                        synonym = feat.synonyms.get(g.name, scf.name)
                        seq_diff = feat.seq_diff.get(g.name, False)
                        term_val = format_number(scf.size, human_readable=human_readable, is_terminal=True)
                        file_val = format_number(scf.size, human_readable=human_readable, is_terminal=False)
                        if synonym != feat.primary_name:
                            suffix = f" ({synonym}, seq diff)" if seq_diff else f" ({synonym})"
                            term_row.append(f"{term_val}{suffix}")
                            file_row.append(f"{file_val}{suffix}")
                        else:
                            term_row.append(term_val)
                            file_row.append(file_val)

                if is_two_genomes:
                    s1, s2 = sizes[0], sizes[1]
                    diff = None
                    if s1 is not None and s2 is not None:
                        diff = s2 - s1
                    elif s1 is not None and s2 is None:
                        diff = -s1
                    elif s1 is None and s2 is not None:
                        diff = s2
                    term_row.append(format_diff(diff, human_readable=human_readable, is_terminal=True))
                    file_row.append(format_diff(diff, human_readable=human_readable, is_terminal=False))

            else:
                # Reference Mode
                for idx, g in enumerate(genomes):
                    scf = feat.scaffolds.get(g.name)
                    if idx == ref_idx:
                        if scf is None:
                            term_row.append("-")
                            file_row.append("-")
                        else:
                            term_row.append(format_number(scf.size, human_readable=human_readable, is_terminal=True))
                            file_row.append(format_number(scf.size, human_readable=human_readable, is_terminal=False))
                    else:
                        if scf is None:
                            all_same_as_ref = False
                            term_row.append("-")
                            file_row.append("-")
                        elif ref_scf is None:
                            all_same_as_ref = False
                            t_val = format_number(scf.size, human_readable=human_readable, is_terminal=True)
                            f_val = format_number(scf.size, human_readable=human_readable, is_terminal=False)
                            term_row.append(f"{t_val} (unique)")
                            file_row.append(f"{f_val} (unique)")
                        else:
                            diff = scf.size - ref_scf.size
                            synonym = feat.synonyms.get(g.name, scf.name)
                            is_synonym = (synonym != feat.primary_name)
                            seq_diff = feat.seq_diff.get(g.name, False)
                            match_type = feat.match_types.get(g.name, "")

                            if diff == 0:
                                if not is_synonym:
                                    term_row.append("= ref")
                                    file_row.append("= ref")
                                else:
                                    if match_type == "rc_seq":
                                        syn_str = f"= ref ({synonym}, revcomp)"
                                    elif seq_diff:
                                        syn_str = f"= ref ({synonym}, seq diff)"
                                    else:
                                        syn_str = f"= ref ({synonym})"
                                    term_row.append(syn_str)
                                    file_row.append(syn_str)
                            else:
                                all_same_as_ref = False
                                t_diff = format_diff(diff, human_readable=human_readable, is_terminal=True)
                                f_diff = format_diff(diff, human_readable=human_readable, is_terminal=False)
                                t_val = format_number(scf.size, human_readable=human_readable, is_terminal=True)
                                f_val = format_number(scf.size, human_readable=human_readable, is_terminal=False)
                                if is_synonym:
                                    term_row.append(f"{t_val} ({t_diff}) ({synonym})")
                                    file_row.append(f"{f_val} ({f_diff}) ({synonym})")
                                else:
                                    term_row.append(f"{t_val} ({t_diff})")
                                    file_row.append(f"{f_val} ({f_diff})")

            if is_ref_mode and diff_only and all_same_as_ref:
                continue

            terminal_rows.append(term_row)
            file_rows.append(file_row)

    # Build Summary Rows
    summary_metrics = [
        ("Total Chromosome Size", "chromosome_size", False),
        ("Total Scaffolds / Contigs", "num_scaffolds", False),
        ("Total Scaffold Size", "scaffold_size", False),
        ("Total Assembly Size", "total_size", False),
        ("N50", "n50", False),
        ("L50", "l50", False),
        ("auN", "auN", False),
        ("N90", "n90", False),
        ("L90", "l90", False),
    ]

    if exp_size is not None:
        summary_metrics.extend([
            ("Estimated Genome Size", "estimated_genome_size", False),
            ("NG50", "ng50", False),
            ("LG50", "lg50", False),
            ("auNG", "auNG", False),
        ])

    summary_metrics.extend([
        ("GC Content (%)", "gc_content", True),
        ("Gap Content (N%)", "gap_content", True),
    ])

    if any((s.get("soft_masked_pct") or 0) > 0 for s in stats_list):
        summary_metrics.append(("Soft-masked (%)", "soft_masked_pct", True))

    terminal_summary_rows = []
    file_summary_rows = []

    if not contigs_only:
        for label, key, is_pct in summary_metrics:
            term_row = [label]
            file_key = label.upper().replace(" ", "_").replace("_(%)", "").replace("_(N%)", "")
            file_row = [file_key]

            is_size = "size" in key.lower() or key.lower() in ("n50", "n90", "aun", "ng50", "aung")
            vals = [s.get(key) for s in stats_list]

            if not is_ref_mode:
                for val in vals:
                    term_row.append(format_number(val, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct))
                    file_row.append(format_number(val, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct))

                if is_two_genomes:
                    v1, v2 = vals[0], vals[1]
                    diff = (v2 - v1) if (v1 is not None and v2 is not None) else None
                    term_row.append(format_diff(diff, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct))
                    file_row.append(format_diff(diff, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct))

            else:
                ref_val = vals[ref_idx]
                all_metrics_same = True

                for idx, val in enumerate(vals):
                    if idx == ref_idx:
                        term_row.append(format_number(ref_val, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct))
                        file_row.append(format_number(ref_val, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct))
                    else:
                        if val == ref_val:
                            term_row.append("= ref")
                            file_row.append("= ref")
                        else:
                            all_metrics_same = False
                            diff = (val - ref_val) if (val is not None and ref_val is not None) else None
                            t_diff = format_diff(diff, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct)
                            f_diff = format_diff(diff, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct)
                            t_val = format_number(val, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct)
                            f_val = format_number(val, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct)
                            term_row.append(f"{t_val} ({t_diff})")
                            file_row.append(f"{f_val} ({f_diff})")

                if diff_only and all_metrics_same:
                    continue

            terminal_summary_rows.append(term_row)
            file_summary_rows.append(file_row)

    # Render Terminal
    if not quiet:
        output_str = render_terminal_table(headers, terminal_rows, terminal_summary_rows)
        typer.echo(output_str)

    # Export to File if requested
    export_path = None
    if output_file:
        export_path = Path(output_file)
    elif output_dir:
        export_path = Path(output_dir) / "genome_summary.tsv"

    if export_path:
        export_path.parent.mkdir(parents=True, exist_ok=True)
        delimiter = "," if export_path.suffix.lower() == ".csv" else "\t"

        file_headers = ["Feature"]
        for idx, g in enumerate(genomes):
            if is_ref_mode and idx == ref_idx:
                file_headers.append(f"{g.name} [Ref]")
            else:
                file_headers.append(g.name)
        if is_two_genomes and not is_ref_mode:
            file_headers.append("Diff")

        with open(export_path, "w", encoding="utf-8") as f_out:
            f_out.write(delimiter.join(file_headers) + "\n")
            for row in file_rows:
                f_out.write(delimiter.join(row) + "\n")
            for row in file_summary_rows:
                f_out.write(delimiter.join(row) + "\n")

        if not quiet:
            typer.echo(f"\nGenome summary exported to: {export_path}")


if __name__ == "__main__":
    app()
