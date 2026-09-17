import os
import re
import typer
from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

from ..genome import Genome

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

    parts = [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', name)]
    return (category, parts)


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
):
    """
    Summarise and compare chromosome sizes and assembly statistics for one or more genomes.
    """
    if not genome_files:
        typer.echo("Error: At least one genome FASTA file must be provided.", err=True)
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
    stats_list = [g.get_stats(estimated_genome_size=exp_size) for g in genomes]

    # Collect features
    feature_set = set()
    feature_sizes: dict[str, dict[str, int]] = {}

    for g in genomes:
        show_all = include_all and not chromosomes_only
        for name, scf in g.scaffolds.items():
            is_chr = scf.chromosome or scf.unknown_chromosome or scf.organelle
            if show_all or is_chr:
                feature_set.add(name)
                if name not in feature_sizes:
                    feature_sizes[name] = {}
                feature_sizes[name][g.name] = scf.size

    # Sort features
    if sort_by == "order":
        ordered = []
        for g in genomes:
            for name in g.scaffolds:
                if name in feature_set and name not in ordered:
                    ordered.append(name)
        sorted_features = ordered
    elif sort_by == "size":
        def avg_size(fname):
            sizes = [feature_sizes[fname].get(g.name, 0) for g in genomes]
            return max(sizes)
        sorted_features = sorted(feature_set, key=lambda f: (-avg_size(f), f))
    else:
        # Default: natural sort
        sorted_features = sorted(feature_set, key=lambda f: get_natural_sort_key(f, genomes))

    # Prepare Headers
    unit_label = "" if human_readable else " (bp)"
    headers = ["Feature / Metric"]
    for g in genomes:
        headers.append(f"{g.name}{unit_label}")
    if is_two_genomes:
        headers.append(f"Diff{unit_label}")

    # Build Feature Rows
    terminal_rows = []
    file_rows = []

    if not summary_only:
        for feat in sorted_features:
            term_row = [feat]
            file_row = [feat]
            sizes = []
            for g in genomes:
                sz = feature_sizes.get(feat, {}).get(g.name, None)
                sizes.append(sz)
                term_row.append(format_number(sz, human_readable=human_readable, is_terminal=True))
                file_row.append(format_number(sz, human_readable=human_readable, is_terminal=False))

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

    terminal_summary_rows = []
    file_summary_rows = []

    for label, key, is_pct in summary_metrics:
        term_row = [label]
        file_key = label.upper().replace(" ", "_").replace("_(%)", "").replace("_(N%)", "")
        file_row = [file_key]

        is_size = "size" in key.lower() or key.lower() in ("n50", "n90", "aun", "ng50", "aung")

        vals = [s.get(key) for s in stats_list]
        for val in vals:
            term_row.append(format_number(val, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct))
            file_row.append(format_number(val, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct))

        if is_two_genomes:
            v1, v2 = vals[0], vals[1]
            diff = (v2 - v1) if (v1 is not None and v2 is not None) else None
            term_row.append(format_diff(diff, human_readable=(human_readable and is_size), is_terminal=True, is_pct=is_pct))
            file_row.append(format_diff(diff, human_readable=(human_readable and is_size), is_terminal=False, is_pct=is_pct))

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

        file_headers = ["Feature"] + [g.name for g in genomes]
        if is_two_genomes:
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
