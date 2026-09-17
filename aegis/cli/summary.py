import os
import re
import typer
from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from .summary_genome import pair_genome_features, PairedFeature, normalize_chr_name


app = typer.Typer(add_completion=False, no_args_is_help=True)


def is_fasta_path(filepath: str) -> bool:
    """Check whether a filepath appears to be a FASTA genome file."""
    lower = filepath.lower()
    for ext in [".fasta", ".fa", ".fna", ".faa", ".fasta.gz", ".fa.gz", ".fna.gz", ".faa.gz"]:
        if lower.endswith(ext):
            return True
    if os.path.isfile(filepath) and not lower.endswith(".gz"):
        try:
            with open(filepath, "rt", encoding="utf-8", errors="ignore") as f:
                first_line = f.readline().strip()
                if first_line.startswith(">"):
                    return True
        except Exception:
            pass
    return False


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


def format_number(val: int | float | None, human_readable: bool = False, is_terminal: bool = True, is_pct: bool = False) -> str:
    """Format a metric value for terminal or tabular output."""
    if val is None:
        return "-"
    if is_pct:
        return f"{val:.2f}%" if is_terminal else f"{val:.2f}"
    if human_readable and isinstance(val, (int, float)):
        return format_human_readable(val)
    if is_terminal:
        if isinstance(val, float):
            return f"{val:,.2f}"
        if isinstance(val, int):
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
        if isinstance(diff, int):
            return f"{sign}{diff:,}"
    return f"{sign}{diff}"


def render_terminal_table(headers: list[str], rows: list[list[str]], section_title: str = "", summary_rows: list[list[str]] | None = None) -> str:
    """Renders a cleanly formatted ASCII table with aligned columns."""
    all_rows = list(rows)
    if summary_rows:
        all_rows += summary_rows
    num_cols = len(headers)
    col_widths = [0] * num_cols

    for i in range(num_cols):
        max_ci = len(headers[i])
        for r in all_rows:
            val = r[i] if i < len(r) else ""
            max_ci = max(max_ci, len(val))
        col_widths[i] = max_ci + 2

    total_width = sum(col_widths) + (num_cols - 1)

    lines = []
    lines.append("=" * total_width)
    if section_title:
        lines.append(f" {section_title}")
        lines.append("-" * total_width)

    # Header line
    header_line = headers[0].ljust(col_widths[0])
    for i in range(1, num_cols):
        header_line += " " + headers[i].rjust(col_widths[i])
    lines.append(header_line)
    lines.append("-" * total_width)

    # Rows
    for row in rows:
        line = row[0].ljust(col_widths[0])
        for i in range(1, num_cols):
            val = row[i] if i < len(row) else "-"
            line += " " + val.rjust(col_widths[i])
        lines.append(line)

    # Summary rows if present
    if summary_rows:
        if rows:
            lines.append("-" * total_width)
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
    files: Annotated[List[str], typer.Argument(
        help="Path to one or more annotation GFF/GTF file(s). (Optional: a single genome FASTA can be provided as the last argument, or explicitly via -g/--genome)."
    )],
    genome: Annotated[Optional[List[str]], typer.Option(
        "-g", "--genome", "--genome-file", help="Path to input genome FASTA file(s). Provide 1 file for shared assembly, or 1-to-1 matching annotations (comma-separated or repeated -g)."
    )] = None,
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", "--annotation-name", help="Comma-separated annotation names or tags (defaults to filenames)."
    )] = "",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag (comma-separated if multiple genomes)."
    )] = "{genome-file}",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to output summary table (TSV/CSV)."
    )] = "",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder for stats and reports."
    )] = "./aegis_output/stats/",
    reference: Annotated[bool, typer.Option(
        "-r", "--reference", help="Use first annotation as reference (or specified via --ref-annotation) and report relative differences."
    )] = False,
    ref_annotation: Annotated[str, typer.Option(
        "--ref-annotation", "--ref-annot", help="Specify an annotation name or 1-based index to use as reference."
    )] = "",
    diff_only: Annotated[bool, typer.Option(
        "--diff-only", help="Report only features and summary statistics where annotations differ from reference (automatically activates reference mode; hides rows that are '= ref')."
    )] = False,
    summary_only: Annotated[bool, typer.Option(
        "--summary-only", help="Report only overall summary statistics without listing individual contigs."
    )] = False,
    contigs_only: Annotated[bool, typer.Option(
        "--contigs-only", help="Report only contig-level statistics without the summary table."
    )] = False,
    chromosomes_only: Annotated[bool, typer.Option(
        "--chromosomes-only", help="Report only chromosomes in table and exclude unplaced scaffolds/contigs."
    )] = False,
    human_readable: Annotated[bool, typer.Option(
        "-H", "--human-readable", help="Display sizes in human-readable units (e.g., Kb, Mb, Gb)."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    plots: Annotated[bool, typer.Option(
        "--plots", help="Export distribution barplots and pie charts into output directory."
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
):
    """
    Outputs summary statistics and chromosome breakdown for one or more annotations,
    with optional assembly reconciliation checks against one or more genome FASTA files.

    USAGE SCENARIOS:
      1. Single annotation with genome reconciliation:
         aegis summary annot.gff -g genome.fasta

      2. Comparing multiple annotations on the SAME genome assembly (e.g. NCBI vs Ensembl):
         aegis summary ncbi.gff ensembl.gff -g genome.fasta

      3. Comparing multiple annotations on SYNONYMOUS assemblies (automatic contig matching):
         aegis summary ncbi.gff ensembl.gff -g ncbi.fa,ensembl.fa

      4. Comparing annotations from DIFFERENT genomes/species (macro statistics):
         aegis summary speciesA.gff speciesB.gff --summary-only
    """
    if not files:
        typer.echo("Error: At least one annotation GFF/GTF file must be provided.", err=True)
        raise typer.Exit(code=1)

    if summary_only and contigs_only:
        typer.echo("Error: Cannot specify both --summary-only and --contigs-only.", err=True)
        raise typer.Exit(code=1)

    # 1. Disambiguate positional arguments vs genome file
    raw_genome_files: list[str] = []
    if genome:
        g_list = [genome] if isinstance(genome, str) else genome
        for g_arg in g_list:
            for part in g_arg.split(","):
                part = part.strip()
                if part:
                    raw_genome_files.append(part)

    if raw_genome_files:
        annot_files = list(files)
        genome_files = raw_genome_files
    else:
        if len(files) == 1:
            annot_files = [files[0]]
            genome_files = []
        elif len(files) == 2:
            if is_fasta_path(files[1]):
                annot_files = [files[0]]
                genome_files = [files[1]]
            else:
                annot_files = list(files)
                genome_files = []
        else:
            if is_fasta_path(files[-1]):
                annot_files = list(files[:-1])
                genome_files = [files[-1]]
            else:
                annot_files = list(files)
                genome_files = []

    num_annots = len(annot_files)
    num_genomes = len(genome_files)

    if num_genomes > 1 and num_genomes != num_annots:
        typer.echo(
            f"Error: Number of genome files ({num_genomes}) must match the number of annotation files ({num_annots}) or be exactly 1.",
            err=True,
        )
        raise typer.Exit(code=1)

    is_multi_genome = (num_genomes > 1 and num_genomes == num_annots)
    os.makedirs(output_dir, exist_ok=True)

    # 2. Load genome assembly if specified
    gnames = [item.strip() for item in genome_name.split(",")] if genome_name and genome_name != "{genome-file}" else []
    genome_objs: list[Genome] = []

    for g_idx, g_file in enumerate(genome_files):
        if g_idx < len(gnames) and gnames[g_idx]:
            cur_gname = gnames[g_idx]
        else:
            g_stem = Path(g_file).stem
            if g_stem.endswith(".fa") or g_stem.endswith(".fasta"):
                g_stem = Path(g_stem).stem
            cur_gname = g_stem

        if any(g.name == cur_gname for g in genome_objs):
            cur_gname = f"{cur_gname}_{g_idx + 1}"

        g_obj = Genome(
            name=cur_gname,
            genome_file_path=g_file,
            quiet=True,
            header_id_tag=header_id_tag if header_id_tag != "" else None,
            header_id_regex=header_id_regex if header_id_regex != "" else None,
            gwh=gwh,
        )
        genome_objs.append(g_obj)

    genome_obj: Genome | None = genome_objs[0] if (num_genomes == 1) else None
    has_genomes = bool(genome_obj is not None or (is_multi_genome and genome_objs))

    # 3. Parse annotation names
    names = [item.strip() for item in annotation_names.split(",")] if annotation_names else []

    # 4. Load annotation objects
    annotations: list[Annotation] = []
    for idx, afile in enumerate(annot_files):
        if idx < len(names) and names[idx]:
            aname = names[idx]
        else:
            aname = Path(afile).stem
            if aname.endswith(".gff3") or aname.endswith(".gff") or aname.endswith(".gtf"):
                aname = Path(aname).stem

        if is_multi_genome:
            assigned_genome = genome_objs[idx]
        elif genome_objs:
            assigned_genome = genome_objs[0]
        else:
            assigned_genome = None

        def emit_mismatch_hint():
            if len(annot_files) > 1 and not is_multi_genome:
                typer.echo(
                    "\nHint: 'aegis summary' applies the provided genome assembly to ALL input annotations "
                    "(assuming they correspond to the same genome build).\n"
                    "• If your annotations use different assemblies or contig naming conventions, pass their respective genomes: '-g g1.fa,g2.fa'.\n"
                    "• If comparing annotations from different species without genomes, run without '-g/--genome'.\n"
                    "• To compare overall structural metrics across different assemblies without contig matching, pass '--summary-only'.",
                    err=True,
                )
            elif is_multi_genome:
                typer.echo(
                    f"\nHint: None of the contig identifiers in annotation '{aname}' match sequences in its assigned genome FASTA '{assigned_genome.name}'.\n"
                    "• Check if chromosome naming differs (e.g. 'chr1' vs '1', or accession IDs vs chromosome names).\n"
                    "• For GWH or tagged FASTA headers, try '--gwh' or '--header-id-tag'.",
                    err=True,
                )
            else:
                typer.echo(
                    "\nHint: None of the contig identifiers in this annotation match sequences in the genome FASTA.\n"
                    "• Check if chromosome naming differs (e.g. 'chr1' vs '1', or accession IDs vs chromosome names).\n"
                    "• For GWH or tagged FASTA headers, try '--gwh' or '--header-id-tag'.",
                    err=True,
                )

        try:
            annot = Annotation(
                name=aname,
                annot_file_path=afile,
                genome=assigned_genome,
                quiet=True
            )
        except ValueError as e:
            typer.echo(f"Error: {e}", err=True)
            if "match" in str(e).lower() or "scaffold" in str(e).lower() or "chromosome" in str(e).lower():
                emit_mismatch_hint()
            raise typer.Exit(code=1)

        # Check for fatal mismatch (zero matching contigs between annotation and genome)
        v = getattr(annot, "genome_validation", {})
        if v.get("fatal_mismatch"):
            typer.echo(f"Error: {v.get('summary_message')}", err=True)
            emit_mismatch_hint()
            raise typer.Exit(code=1)

        # Update stats and export basic/full stats CSV files
        annot.stats.update(output_dir=output_dir, export=True, quiet=True)
        annotations.append(annot)

    num_annots = len(annotations)
    is_multi = num_annots > 1
    is_two = num_annots == 2

    # 5. Determine reference annotation if reference mode requested
    ref_idx = 0
    if ref_annotation:
        if ref_annotation.isdigit():
            idx = int(ref_annotation) - 1
            if 0 <= idx < num_annots:
                ref_idx = idx
            else:
                typer.echo(f"Warning: Reference index '{ref_annotation}' out of range. Defaulting to 1st annotation.", err=True)
        else:
            found = False
            for idx, a in enumerate(annotations):
                if a.name.lower() == ref_annotation.lower():
                    ref_idx = idx
                    found = True
                    break
            if not found:
                typer.echo(f"Warning: Reference annotation '{ref_annotation}' not found. Defaulting to 1st annotation.", err=True)

    if diff_only and not is_multi:
        typer.echo("Warning: --diff-only requires at least two annotations to compare. Ignoring.", err=True)
        diff_only = False

    is_ref_mode = bool(reference or ref_annotation or diff_only) and is_multi

    # 6. Assembly Reconciliation & Genome Pairing
    paired_features: list[PairedFeature] = []
    multi_genome_warning = ""
    if is_multi_genome:
        paired_features = pair_genome_features(
            genome_objs,
            ref_idx=ref_idx,
            check_sequence=True,
            include_all=not chromosomes_only,
            chromosomes_only=chromosomes_only,
        )
        shared_count = sum(1 for pf in paired_features if len(pf.scaffolds) > 1)
        if shared_count == 0 and len(paired_features) > 0:
            multi_genome_warning = (
                "⚠️  MULTI-GENOME ASSEMBLY WARNING:\n"
                "   Zero contigs or sequences could be matched between the provided genome assemblies.\n"
                "   The annotations appear to be from completely different species or assemblies without synonymous contigs.\n"
                "   Reporting summary statistics only."
            )
            summary_only = True

    warning_banners = []
    if not is_multi_genome and genome_obj is not None:
        for a in annotations:
            v = getattr(a, "genome_validation", {})
            if v.get("has_warning"):
                out_of_bounds = v.get("out_of_bounds_features", [])
                missing_chroms = v.get("missing_chromosomes", [])
                msg_lines = [f"⚠️  ASSEMBLY RECONCILIATION WARNING for '{a.name}' on '{genome_obj.name}':"]
                if missing_chroms:
                    if len(missing_chroms) == 1:
                        msg_lines.append(f"   • 1 contig in annotation not found in genome FASTA: '{missing_chroms[0]}'.")
                    else:
                        sample = f" ({', '.join(missing_chroms[:3])}{'...' if len(missing_chroms) > 3 else ''})"
                        msg_lines.append(f"   • {len(missing_chroms):,} contigs in annotation not found in genome FASTA{sample}.")
                if out_of_bounds:
                    if len(out_of_bounds) == 1:
                        msg_lines.append(f"   • 1 feature has coordinates exceeding contig boundaries: '{out_of_bounds[0]}'.")
                    else:
                        sample = f" (e.g., {', '.join(out_of_bounds[:3])}{'...' if len(out_of_bounds) > 3 else ''})"
                        msg_lines.append(f"   • {len(out_of_bounds):,} features have coordinates exceeding contig boundaries{sample}.")
                msg_lines.append("   Possible assembly mismatch! Please verify that annotation and genome builds correspond.")
                warning_banners.append("\n".join(msg_lines))
    elif is_multi_genome:
        for idx, a in enumerate(annotations):
            g_obj = genome_objs[idx]
            v = getattr(a, "genome_validation", {})
            if v.get("has_warning"):
                out_of_bounds = v.get("out_of_bounds_features", [])
                missing_chroms = v.get("missing_chromosomes", [])
                msg_lines = [f"⚠️  ASSEMBLY RECONCILIATION WARNING for '{a.name}' on '{g_obj.name}':"]
                if missing_chroms:
                    if len(missing_chroms) == 1:
                        msg_lines.append(f"   • 1 contig in annotation not found in genome FASTA: '{missing_chroms[0]}'.")
                    else:
                        sample = f" ({', '.join(missing_chroms[:3])}{'...' if len(missing_chroms) > 3 else ''})"
                        msg_lines.append(f"   • {len(missing_chroms):,} contigs in annotation not found in genome FASTA{sample}.")
                if out_of_bounds:
                    if len(out_of_bounds) == 1:
                        msg_lines.append(f"   • 1 feature has coordinates exceeding contig boundaries: '{out_of_bounds[0]}'.")
                    else:
                        sample = f" (e.g., {', '.join(out_of_bounds[:3])}{'...' if len(out_of_bounds) > 3 else ''})"
                        msg_lines.append(f"   • {len(out_of_bounds):,} features have coordinates exceeding contig boundaries{sample}.")
                msg_lines.append("   Possible assembly mismatch! Please verify that annotation and genome builds correspond.")
                warning_banners.append("\n".join(msg_lines))

    # 7. Collect Contig Statistics
    annot_contig_stats = [a.stats.get_contig_stats() for a in annotations]
    all_contig_names = []
    contig_maps = []
    for cstats in annot_contig_stats:
        cmap = {item["contig"]: item for item in cstats}
        contig_maps.append(cmap)
        for cname in cmap:
            if cname not in all_contig_names:
                all_contig_names.append(cname)

    # Check for disjoint contigs when multiple annotations are compared without a genome
    disjoint_contig_banner = ""
    if is_multi and not is_multi_genome and genome_obj is None and not summary_only:
        contig_sets = [set(cmap.keys()) for cmap in contig_maps if cmap]
        if len(contig_sets) >= 2:
            common_contigs = set.intersection(*contig_sets)
            if not common_contigs:
                disjoint_contig_banner = (
                    "ℹ️  NOTICE: Input annotations share no common contig names "
                    "(likely from different genomes or differing chromosome naming conventions).\n"
                    "   Tip: Pass '--summary-only' to suppress the disjoint contig breakdown and view overall comparative metrics."
                )

    def contig_sort_key(name: str):
        nl = name.lower()
        if "mit" in nl or "mt" in nl or "pt" in nl or "chlor" in nl or "cp" in nl:
            cat = 3
        elif nl.startswith("chr") or any(nl.startswith(p) for p in ["ch", "scaffold", "contig"]) or name.isdigit():
            cat = 1
        else:
            cat = 2
        parts = [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', name)]
        return (cat, parts)

    if not is_multi_genome:
        all_contig_names.sort(key=contig_sort_key)

        if chromosomes_only:
            filtered_contigs = []
            for cname in all_contig_names:
                if genome_obj:
                    scf = genome_obj.get_scaffold(cname)
                    if scf:
                        if scf.chromosome and not scf.organelle and not scf.unknown_chromosome:
                            filtered_contigs.append(cname)
                        continue
                nl = cname.lower()
                is_organelle = any(p in nl for p in ["mit", "mt", "pt", "chlor", "cp"])
                if (nl.startswith("chr") or nl.startswith("chromosome") or cname.isdigit()) and not is_organelle:
                    filtered_contigs.append(cname)
            all_contig_names = filtered_contigs
    else:
        paired_features.sort(key=lambda pf: contig_sort_key(pf.primary_name))

    contig_headers = []
    contig_rows = []
    contig_file_rows = []

    if not summary_only:
        if not is_multi:
            # Single annotation table
            a = annotations[0]
            if genome_obj is not None:
                contig_headers = ["Contig", "Size", "Genes", "Coding", "Non-coding", "Transcripts", "Density (g/Mb)"]
            else:
                contig_headers = ["Contig", "Genes", "Coding", "Non-coding", "Transcripts"]

            for cname in all_contig_names:
                info = contig_maps[0].get(cname)
                if not info:
                    continue
                row = [cname]
                if genome_obj is not None:
                    row.append(format_number(info["size"], human_readable=human_readable))
                row.append(str(info["genes"]))
                row.append(str(info["coding_genes"]))
                row.append(str(info["noncoding_genes"]))
                row.append(str(info["transcripts"]))
                if genome_obj is not None:
                    dens_str = f"{info['density_genes_per_mb']:.2f}" if info["density_genes_per_mb"] is not None else "-"
                    row.append(dens_str)
                contig_rows.append(row)
                contig_file_rows.append(row)

        else:
            # Multi-annotation table
            if has_genomes:
                contig_headers = ["Contig", "Size"]
            else:
                contig_headers = ["Contig"]

            if is_ref_mode:
                for idx, a in enumerate(annotations):
                    if idx == ref_idx:
                        contig_headers.append(f"{a.name} (ref)")
                    else:
                        contig_headers.append(f"{a.name} vs ref")
            elif is_two:
                for a in annotations:
                    contig_headers.append(f"{a.name} Genes")
                contig_headers.append("Diff")
            else:
                for a in annotations:
                    contig_headers.append(f"{a.name} Genes")

            if is_multi_genome:
                for pf in paired_features:
                    syn_parts = []
                    for g_obj in genome_objs:
                        if g_obj.name in pf.synonyms:
                            syn_val = pf.synonyms[g_obj.name]
                            if syn_val != pf.primary_name and syn_val not in syn_parts:
                                syn_parts.append(syn_val)

                    if len(pf.scaffolds) == 1:
                        owner_g = next(iter(pf.scaffolds.keys()))
                        owner_annot_name = annotations[0].name
                        for a_idx, g_obj in enumerate(genome_objs):
                            if g_obj.name == owner_g:
                                owner_annot_name = annotations[a_idx].name
                                break
                        c_term = f"{pf.primary_name} [{owner_annot_name} only]"
                        c_file = f"{pf.primary_name} [{owner_annot_name}]"
                    elif syn_parts:
                        c_term = f"{pf.primary_name} ({', '.join(syn_parts)})"
                        c_file = f"{pf.primary_name} ({', '.join(syn_parts)})"
                    else:
                        c_term = pf.primary_name
                        c_file = pf.primary_name

                    term_row = [c_term]
                    file_row = [c_file]

                    ref_scf = pf.scaffolds.get(genome_objs[ref_idx].name) or next(iter(pf.scaffolds.values()))
                    scf_sz = ref_scf.size if ref_scf else None
                    term_row.append(format_number(scf_sz, human_readable=human_readable, is_terminal=True))
                    file_row.append(format_number(scf_sz, human_readable=human_readable, is_terminal=False))

                    gene_counts = []
                    for a_idx, g_obj in enumerate(genome_objs):
                        cname_in_g = pf.synonyms.get(g_obj.name)
                        cmap = contig_maps[a_idx]
                        info = cmap.get(cname_in_g) if cname_in_g else None
                        gene_counts.append(info["genes"] if info else 0)

                    if is_ref_mode:
                        ref_count = gene_counts[ref_idx]
                        all_same = True
                        for idx, cnt in enumerate(gene_counts):
                            if idx == ref_idx:
                                term_row.append(str(cnt))
                                file_row.append(str(cnt))
                            else:
                                diff = cnt - ref_count
                                if diff == 0:
                                    term_row.append("= ref")
                                    file_row.append("= ref")
                                else:
                                    all_same = False
                                    term_row.append(format_diff(diff))
                                    file_row.append(format_diff(diff, is_terminal=False))
                        if diff_only and all_same:
                            continue

                    elif is_two:
                        term_row.append(str(gene_counts[0]))
                        file_row.append(str(gene_counts[0]))
                        term_row.append(str(gene_counts[1]))
                        file_row.append(str(gene_counts[1]))
                        diff = gene_counts[1] - gene_counts[0]
                        term_row.append(format_diff(diff))
                        file_row.append(format_diff(diff, is_terminal=False))
                        if diff_only and diff == 0:
                            continue

                    else:
                        for cnt in gene_counts:
                            term_row.append(str(cnt))
                            file_row.append(str(cnt))

                    contig_rows.append(term_row)
                    contig_file_rows.append(file_row)

            else:
                for cname in all_contig_names:
                    term_row = [cname]
                    file_row = [cname]

                    if genome_obj is not None:
                        scf = genome_obj.scaffolds.get(cname)
                        scf_sz = scf.size if scf else None
                        term_row.append(format_number(scf_sz, human_readable=human_readable, is_terminal=True))
                        file_row.append(format_number(scf_sz, human_readable=human_readable, is_terminal=False))

                    gene_counts = []
                    for cmap in contig_maps:
                        info = cmap.get(cname)
                        gene_counts.append(info["genes"] if info else 0)

                    if is_ref_mode:
                        ref_count = gene_counts[ref_idx]
                        all_same = True
                        for idx, cnt in enumerate(gene_counts):
                            if idx == ref_idx:
                                term_row.append(str(cnt))
                                file_row.append(str(cnt))
                            else:
                                diff = cnt - ref_count
                                if diff == 0:
                                    term_row.append("= ref")
                                    file_row.append("= ref")
                                else:
                                    all_same = False
                                    term_row.append(format_diff(diff))
                                    file_row.append(format_diff(diff, is_terminal=False))
                        if diff_only and all_same:
                            continue

                    elif is_two:
                        term_row.append(str(gene_counts[0]))
                        file_row.append(str(gene_counts[0]))
                        term_row.append(str(gene_counts[1]))
                        file_row.append(str(gene_counts[1]))
                        diff = gene_counts[1] - gene_counts[0]
                        term_row.append(format_diff(diff))
                        file_row.append(format_diff(diff, is_terminal=False))
                        if diff_only and diff == 0:
                            continue

                    else:
                        for cnt in gene_counts:
                            term_row.append(str(cnt))
                            file_row.append(str(cnt))

                    contig_rows.append(term_row)
                    contig_file_rows.append(file_row)

    # 8. Collect Overall Annotation Summary Metrics
    summary_metrics = [
        ("Total Genes", "total_genes", False),
        ("Protein-coding Genes", "coding_genes", False),
        ("Non-coding Genes", "noncoding_genes", False),
        ("Total Transcripts", "total_transcripts", False),
        ("Mean Transcripts / Gene", "mean_transcripts", True),
        ("Mean Exons / Transcript", "mean_exons", True),
        ("Mean Gene Size (bp)", "mean_gene_size", False),
        ("Mean Transcript Size (bp)", "mean_transcript_size", False),
        ("Mean CDS Size (bp)", "mean_CDS_size", False),
        ("Mean Intron Size (bp)", "mean_intron_size", False),
        ("Total Gene Span", "total_length_gene", False),
        ("Total mRNA Length", "total_length_mRNA", False),
    ]
    if has_genomes:
        summary_metrics.extend([
            ("Out-of-bounds Features", "out_of_bounds", False),
            ("Contigs Missing in Genome", "missing_chroms", False),
            ("Unannotated Scaffolds", "unannotated_scaffolds", False),
        ])

    def get_annot_metric_val(annot: Annotation, metric_key: str):
        stats_data = annot.stats.data
        if metric_key == "total_genes":
            return sum(len(genes) for genes in annot.chrs.values())
        if metric_key == "coding_genes":
            val = stats_data.get("coding_genes", 0)
            return len(val) if isinstance(val, (list, set)) else val
        if metric_key == "noncoding_genes":
            val = stats_data.get("noncoding_genes", 0)
            return len(val) if isinstance(val, (list, set)) else val
        if metric_key == "total_transcripts":
            return sum(len(g.transcripts) for genes in annot.chrs.values() for g in genes.values())
        if metric_key == "out_of_bounds":
            return len(annot.warnings.get("feature_exceeds_scaffold_length", set()))
        if metric_key == "missing_chroms":
            return len(annot.warnings.get("chromosome_not_in_genome", set()))
        if metric_key == "unannotated_scaffolds":
            return len(annot.genome_validation.get("unannotated_scaffolds", [])) if hasattr(annot, "genome_validation") else 0
        return stats_data.get(metric_key, 0)

    summary_headers = ["Metric"]
    if is_ref_mode:
        for idx, a in enumerate(annotations):
            if idx == ref_idx:
                summary_headers.append(f"{a.name} (ref)")
            else:
                summary_headers.append(f"{a.name} vs ref")
    elif is_two:
        summary_headers.extend([annotations[0].name, annotations[1].name, "Diff"])
    else:
        for a in annotations:
            summary_headers.append(a.name)

    summary_table_rows = []
    summary_file_rows = []

    if not contigs_only:
        for label, mkey, is_float in summary_metrics:
            vals = [get_annot_metric_val(a, mkey) for a in annotations]
            term_row = [label]
            file_row = [label]

            is_size_metric = "size" in mkey.lower() or "length" in mkey.lower() or "span" in mkey.lower()

            if not is_multi:
                val = vals[0]
                if is_size_metric and human_readable:
                    term_row.append(format_human_readable(val))
                    file_row.append(str(val))
                elif is_float:
                    term_row.append(f"{val:.2f}")
                    file_row.append(f"{val:.2f}")
                else:
                    term_row.append(format_number(val, human_readable=False, is_terminal=True))
                    file_row.append(str(val))

            elif is_ref_mode:
                ref_val = vals[ref_idx]
                all_same = True
                for idx, val in enumerate(vals):
                    if idx == ref_idx:
                        if is_size_metric and human_readable:
                            term_row.append(format_human_readable(val))
                        elif is_float:
                            term_row.append(f"{val:.2f}")
                        else:
                            term_row.append(format_number(val, human_readable=False, is_terminal=True))
                        file_row.append(str(val))
                    else:
                        diff = val - ref_val
                        if diff == 0:
                            term_row.append("= ref")
                            file_row.append("= ref")
                        else:
                            all_same = False
                            if is_size_metric and human_readable:
                                term_row.append(format_diff(diff, human_readable=True))
                            elif is_float:
                                term_row.append(f"{diff:+.2f}")
                            else:
                                term_row.append(format_diff(diff))
                            file_row.append(f"{diff:+}")

                if diff_only and all_same:
                    continue

            elif is_two:
                v1, v2 = vals[0], vals[1]
                diff = v2 - v1
                for v in [v1, v2]:
                    if is_size_metric and human_readable:
                        term_row.append(format_human_readable(v))
                    elif is_float:
                        term_row.append(f"{v:.2f}")
                    else:
                        term_row.append(format_number(v, human_readable=False, is_terminal=True))
                    file_row.append(str(v))

                if is_size_metric and human_readable:
                    term_row.append(format_diff(diff, human_readable=True))
                elif is_float:
                    term_row.append(f"{diff:+.2f}")
                else:
                    term_row.append(format_diff(diff))
                file_row.append(f"{diff:+}")

                if diff_only and diff == 0:
                    continue

            else:
                for val in vals:
                    if is_size_metric and human_readable:
                        term_row.append(format_human_readable(val))
                    elif is_float:
                        term_row.append(f"{val:.2f}")
                    else:
                        term_row.append(format_number(val, human_readable=False, is_terminal=True))
                    file_row.append(str(val))

            summary_table_rows.append(term_row)
            summary_file_rows.append(file_row)

    # 9. Terminal Output Rendering
    if not quiet:
        output_blocks = []

        if multi_genome_warning:
            border = "=" * 80
            output_blocks.append(f"\n{border}\n{multi_genome_warning}\n{border}")

        if disjoint_contig_banner:
            border = "-" * 80
            output_blocks.append(f"\n{border}\n{disjoint_contig_banner}\n{border}")

        for banner in warning_banners:
            border = "=" * 80
            output_blocks.append(f"\n{border}\n{banner}\n{border}")

        if not summary_only and contig_rows:
            section_lbl = "Contig Breakdown (Genes & Density)" if has_genomes else "Contig Breakdown"
            output_blocks.append(render_terminal_table(contig_headers, contig_rows, section_title=section_lbl))

        if not contigs_only and summary_table_rows:
            section_lbl = "Annotation Summary Statistics"
            output_blocks.append(render_terminal_table(summary_headers, summary_table_rows, section_title=section_lbl))

        if output_blocks:
            typer.echo("\n".join(output_blocks))

    # 10. Export Tables to File if requested
    export_path: Path | None = None
    if output_file:
        export_path = Path(output_file)
    elif output_dir:
        export_path = Path(output_dir) / "annotation_summary.tsv"

    if export_path is not None:
        export_path.parent.mkdir(parents=True, exist_ok=True)
        with open(export_path, "w", encoding="utf-8") as f:
            if not summary_only and contig_file_rows:
                f.write("# Contig Breakdown\n")
                f.write("\t".join(contig_headers) + "\n")
                for r in contig_file_rows:
                    f.write("\t".join(r) + "\n")
                f.write("\n")

            if not contigs_only and summary_file_rows:
                f.write("# Summary Statistics\n")
                f.write("\t".join(summary_headers) + "\n")
                for r in summary_file_rows:
                    f.write("\t".join(r) + "\n")

        if not quiet:
            typer.echo(f"\nSummary table exported to: {export_path}")


if __name__ == "__main__":
    app()
