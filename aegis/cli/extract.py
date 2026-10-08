import typer
import os

from typing import List, Optional
from typing_extensions import Annotated

from ..genome import Genome
from ..annotation import Annotation
from .utils import (
    split_callback,
    detect_file_type,
    TaxonomyOption,
    GeneticCodeOption,
    AutoOrganelleCodesOption,
    MitoCodeOption,
    PlastidCodeOption,
    MitochondriaChromsOption,
    ChloroplastChromsOption,
    InitiatorMethionineOption,
    IO_PANEL,
    EXEC_PANEL,
    FASTA_HEADER_PANEL,
    CDS_PANEL,
    OUTPUT_HEADER_PANEL,
    EXTRACTION_PANEL,
)

FEATURES = ["gene", "transcript", "CDS", "protein", "promoter"]

VALID_IDS = ["gene", "transcript", "CDS", "feature"]

EXTRACTION_MODES = ["all", "main", "unique", "unique_per_gene"]

PROMOTER_TYPES = ["standard", "upstream_ATG", "standard_plus_up_to_ATG"]
from ..conf import RNA_CLASSES

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = "",
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file (or provide via -g/--genome)."
    )] = "",

    # 1. Feature Extraction Options
    features: Annotated[List[str], typer.Option(
        "-f", "--features", help=f"Feature type(s) to extract, as a comma-separated list. Available options: {', '.join(FEATURES)}.",
        callback=split_callback,
        rich_help_panel=EXTRACTION_PANEL,
    )] = ["gene"],
    mode: Annotated[List[str], typer.Option(
        "-m", "--mode", help=f"""Extraction mode(s), as a comma-separated list. Controls filtering of features.\n\n
        - 'all': Extract all features (e.g., all transcripts for a gene).\n
        - 'unique_per_gene': Keep one copy of each unique protein/CDS sequence per gene.\n
        - 'main': Extract only the main variant (e.g., the longest transcript).\n
        - 'unique': Keep only one copy of each unique protein/CDS sequence across the entire output.""",
        callback=split_callback,
        rich_help_panel=EXTRACTION_PANEL,
    )] = ["all"],
    rna_classes: Annotated[List[str], typer.Option(
        "-b", "-r", "--biotypes", "--rna-classes", help=f"Filter transcripts by biotype (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included.",
        callback=split_callback,
        rich_help_panel=EXTRACTION_PANEL,
    )] = [],
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help=f"Size of the promoter region in base pairs (bp). Used only if 'promoter' is a selected feature.",
        rich_help_panel=EXTRACTION_PANEL,
    )] = 2000,
    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help="""\
                Defines the reference point for extracting promoter regions. Used only if 'promoter' is selected.\n\n
                Options:\n
                - 'standard': Upstream of the transcript's start site (TSS).\n
                - 'upstream_ATG': Upstream of the main CDS's start codon (ATG). Falls back to 'standard' if no CDS is present.\n
                - 'standard_plus_up_to_ATG': The 'standard' promoter plus the 5' UTR (sequence between TSS and ATG). Falls back to 'standard' if no CDS.
                """,
        rich_help_panel=EXTRACTION_PANEL,
    )] = "standard",
    raw_cds: Annotated[bool, typer.Option(
        "--raw-cds", help="Export raw spliced genomic CDS sequences instead of in-frame protein-oriented coding sequences.",
        rich_help_panel=EXTRACTION_PANEL,
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons.",
        rich_help_panel=EXTRACTION_PANEL,
    )] = False,
    no_collapse_cds: Annotated[bool, typer.Option(
        "--no-collapse-cds", help="Do not merge overlapping/adjacent CDS segments.",
        rich_help_panel=EXTRACTION_PANEL,
    )] = False,
    strip_stop: Annotated[Optional[bool], typer.Option(
        "--strip-stop/--no-strip-stop", help="Global override: explicitly enable or disable stripping stop codons across both proteins and CDSs.",
        rich_help_panel=EXTRACTION_PANEL,
    )] = None,
    keep_stop: Annotated[bool, typer.Option(
        "--keep-stop", "--keep-protein-stop", help="Keep trailing stop codon (* in proteins) when exporting protein sequences (by default trailing stop codons are stripped).",
        rich_help_panel=EXTRACTION_PANEL,
    )] = False,
    strip_stop_cds: Annotated[bool, typer.Option(
        "--strip-stop-cds", "--strip-cds-stop", help="Strip trailing stop codon (terminal 3-nt stop codon) when exporting CDS sequences (by default CDS sequences retain the stop codon).",
        rich_help_panel=EXTRACTION_PANEL,
    )] = False,

    # 2. CDS Inference & Reworking
    infer_missing_cds: Annotated[bool, typer.Option(
        "--infer-missing-cds", help="Detects and creates CDSs where missing, without overriding existing CDS annotations.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    rework_all_cds: Annotated[bool, typer.Option(
        "--rework-all-cds", help="Recalculates ALL CDSs from the genome sequence, overriding existing ones. More aggressive than --infer-missing-cds.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    fallback_to_trim: Annotated[bool, typer.Option(
        "--fallback-to-trim", help="When recalculating CDSs, fallback to trimming unaligned ends if no high-ratio ORF is found.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    coding_ratio_threshold: Annotated[float, typer.Option(
        "--coding-ratio-threshold", help="Threshold ratio of coding sequence length to transcript length for rework CDS (default: 0.7).",
        rich_help_panel=CDS_PANEL,
    )] = 0.7,
    allow_internal_stops: Annotated[bool, typer.Option(
        "--allow-internal-stops/--no-allow-internal-stops", help="Allow internal stop codons in progressive rework fallback (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    allow_partial: Annotated[bool, typer.Option(
        "--allow-partial/--no-allow-partial", help="Allow partial ORFs without stop codon in progressive rework fallback (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    enforce_start_codon: Annotated[bool, typer.Option(
        "--enforce-start-codon/--no-enforce-start-codon", help="Require start codon (ATG) in initial rework passes (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    orf_choice_mode: Annotated[str, typer.Option(
        "--orf-choice-mode", help="ORF selection criteria: 'longest' or 'earliest' (default: 'longest').",
        rich_help_panel=CDS_PANEL,
    )] = "longest",
    min_codon_len: Annotated[int, typer.Option(
        "--min-codon-len", help="Minimum codon length required for predicted ORFs (e.g. 30 or 50 to suppress micro-ORFs, default: 2).",
        rich_help_panel=CDS_PANEL,
    )] = 2,
    adjust_internal_shifts: Annotated[str, typer.Option(
        "--adjust-internal-shifts", help="Frameshift / phase handling mode for multi-segment CDS translation: 'intra_exon' (default: adjust shifts only across contiguous/overlapping exon segments, preserving continuous splicing across introns), 'all' (adjust across all junctions), or 'none' (translate continuous spliced sequence).",
        rich_help_panel=CDS_PANEL,
    )] = "intra_exon",
    polish_coordinates: Annotated[bool, typer.Option(
        "--polish-coordinates/--skip-coordinate-polishing", help="Mutate feature coordinates when boundaries differ [default: disabled; preserves original coordinates].",
        rich_help_panel=CDS_PANEL,
    )] = False,
    recalculate_phases: Annotated[bool, typer.Option(
        "--recalculate-phases", help="Recalculate CDS segment phases based on segment lengths and splicing leftover, preserving 5' initial phase for partial CDSs.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    reset_phases_zero: Annotated[bool, typer.Option(
        "--reset-phases-zero", help="Reset initial CDS phase to 0 and recalculate all downstream segment phases.",
        rich_help_panel=CDS_PANEL,
    )] = False,

    # 3. Genetic Codes
    taxonomy: TaxonomyOption = "plant",
    genetic_code: GeneticCodeOption = 1,
    auto_organelle_codes: AutoOrganelleCodesOption = True,
    mito_code: MitoCodeOption = None,
    plastid_code: PlastidCodeOption = None,
    mitochondria_chroms: MitochondriaChromsOption = [],
    chloroplast_chroms: ChloroplastChromsOption = [],
    initiator_methionine: InitiatorMethionineOption = "canonical",

    # 4. Output Sequence Header Options
    feature_id: Annotated[str, typer.Option(
        "--feature-id", help=f"Specifies which feature ID to use in FASTA headers. E.g., use 'gene' to label all outputs (transcripts, proteins) with their parent gene ID. 'feature' uses the most specific ID available. Available: {', '.join(VALID_IDS)}.",
        rich_help_panel=OUTPUT_HEADER_PANEL,
    )] = "feature",
    detailed_headers: Annotated[bool, typer.Option(
        "-dh", "--detailed-headers", "--verbose-headers", help=f"Add extra details in fasta headers; scaffold/chromosome number, genome co-ordinates, and/or protein tags if applicable.",
        rich_help_panel=OUTPUT_HEADER_PANEL,
    )] = False,

    # 5. Reference FASTA Options
    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID').",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)').",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers.",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = False,

    # 6. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    genome_file_opt: Annotated[str, typer.Option(
        "-g", "--genome", "--genomes", "--genome-file", help="Path to input genome FASTA file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", help="A name or tag for the annotation version (e.g., 'Araport11'). [default: a name derived from the annotation filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", "--genome-names", help="A name or tag for the genome assembly (e.g., 'TAIR10'). [default: a name derived from the genome FASTA filename]",
        rich_help_panel=IO_PANEL,
    )] = "{genome-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the directory where output FASTA files will be saved.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/features/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", "--output-prefix", help="Optional output filename or prefix.",
        rich_help_panel=IO_PANEL,
    )] = "",

    # 7. Execution & Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Enable detailed console output.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Extract sequences from a genome based on an annotation file.

    This command supports multiple output formats and allows selecting
    specific features (e.g. gene, transcript, CDS, protein, promoter).
    
    Use the --mode and --feature-id flags to control sequence filtering
    and ID labeling. Promoter generation supports multiple strategies,
    including upstream of TSS or ATG.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    genome_in = genome_file_opt if genome_file_opt else genome_file

    if not annot_in or not genome_in:
        if not annot_in and not genome_in:
            raise typer.BadParameter("Both an annotation GFF/GTF file and a genome FASTA file must be provided.")
        elif not annot_in:
            raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
        else:
            raise typer.BadParameter("Missing required genome file. Provide as positional argument or via -g/--genome.")

    # Swap if accidentally provided in reverse order
    if detect_file_type(annot_in) == "fasta" and detect_file_type(genome_in) == "annotation":
        annot_in, genome_in = genome_in, annot_in

    annotation_file = annot_in
    genome_file = genome_in

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_cds
    skip_coordinate_polishing = not polish_coordinates

    for f_type in features:
        if f_type not in FEATURES:
            raise typer.BadParameter(f"Invalid feature type: {f_type}. Choose from: {FEATURES}")

    for m_type in mode:
        if m_type not in EXTRACTION_MODES:
            raise typer.BadParameter(f"Invalid mode: {m_type}. Choose from: {EXTRACTION_MODES}")

    if feature_id not in VALID_IDS:
        raise typer.BadParameter(f"Invalid feature ID: {feature_id}. Choose from: {VALID_IDS}")
    
    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    if promoter_type not in PROMOTER_TYPES:
        raise typer.BadParameter(f"Invalid promoter type: '{promoter_type}'. Choose from: {PROMOTER_TYPES}")

    for rna_class in rna_classes:
        if rna_class not in RNA_CLASSES:
            raise typer.BadParameter(f"Invalid rna class: {rna_class}. Choose from: {RNA_CLASSES}")

    if adjust_internal_shifts not in ("intra_exon", "all", "none"):
        raise typer.BadParameter(f"Invalid adjust_internal_shifts: '{adjust_internal_shifts}'. Choose from: 'intra_exon', 'all', 'none'.")

    if orf_choice_mode not in ("longest", "earliest"):
        raise typer.BadParameter(f"Invalid orf_choice_mode: '{orf_choice_mode}'. Choose from: 'longest', 'earliest'.")


    genome = Genome(
        name=genome_name,
        genome_file_path=genome_file,
        quiet=quiet,
        header_id_tag=header_id_tag if header_id_tag != "" else None,
        header_id_regex=header_id_regex if header_id_regex != "" else None,
        gwh=gwh,
    )

    annotation = Annotation(
        name=annotation_name,
        annot_file_path=annotation_file,
        genome=genome,
        rework_all_CDSs=rework_all_cds,
        work_out_missing_CDSs=infer_missing_cds,
        fallback_to_trim=fallback_to_trim,
        quiet=quiet,
        collapse_exons=collapse_exons,
        collapse_CDSs=collapse_CDSs,
        adjust_internal_shifts=adjust_internal_shifts,
        taxonomy=taxonomy,
        table=genetic_code,
        auto_organelle_codes=auto_organelle_codes,
        mito_table=mito_code,
        plastid_table=plastid_code,
        mitochondria_chroms=mitochondria_chroms,
        chloroplast_chroms=chloroplast_chroms,
        initiator_methionine=initiator_methionine,
        skip_coordinate_polishing=skip_coordinate_polishing,
        coding_ratio_threshold=coding_ratio_threshold,
        allow_internal_stops=allow_internal_stops,
        allow_partial=allow_partial,
        enforce_start_codon=enforce_start_codon,
        orf_choice_mode=orf_choice_mode,
        min_codon_len=min_codon_len,
        recalculate_phases=recalculate_phases,
        reset_phases_zero=reset_phases_zero,
    )

    has_any_cds = any(bool(t.CDSs) for genes in annotation.chrs.values() for g in genes.values() for t in g.transcripts.values())
    if not has_any_cds and not infer_missing_cds and not rework_all_cds:
        typer.secho(
            "Notice: The input annotation does not contain any annotated CDS features. No protein or CDS sequences were extracted.\n"
            "Tip: Pass '--infer-missing-cds' to automatically detect and predict CDSs across transcripts.",
            fg=typer.colors.YELLOW,
            err=True,
        )

    def resolve_export_filename(feat_name: str, mode_name: str) -> Optional[str]:
        if not output_file:
            return None
        fn = output_file
        if "{feature}" in fn:
            fn = fn.replace("{feature}", feat_name)
        if "{mode}" in fn:
            fn = fn.replace("{mode}", mode_name)
        if "{annotation-name}" in fn:
            fn = fn.replace("{annotation-name}", annotation_name)
        if "{genome-name}" in fn:
            fn = fn.replace("{genome-name}", genome_name)
        if fn == output_file and (len(features) > 1 or len(mode) > 1) and not ("{feature}" in output_file or "{mode}" in output_file):
            stem, ext = os.path.splitext(output_file)
            ext = ext if ext else ".fasta"
            fn = f"{stem}_{feat_name}_{mode_name}{ext}"
        elif fn == output_file and not fn.endswith(".fasta") and not fn.endswith(".fa"):
            fn = f"{fn}.fasta"
        return fn

    if "gene" in features:

        annotation.export.genes(output_dir=output_dir, verbose=detailed_headers, filename=resolve_export_filename("gene", "all"))

    if "transcript" in features:

        if "gene" in feature_id:
            used_id = "gene"
        else:
            used_id = "transcript"

        if "unique_per_gene" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, output_dir=output_dir, used_id=used_id, rna_classes=rna_classes, unique_transcripts_per_gene=True, filename=resolve_export_filename("transcript", "unique_per_gene")) #type: ignore
        elif "unique" in mode:
            annotation.export.unique_transcripts(output_dir=output_dir, quiet=quiet, rna_classes=rna_classes, filename=resolve_export_filename("transcript", "unique")) #type: ignore
        elif "all" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, output_dir=output_dir, used_id=used_id, rna_classes=rna_classes, filename=resolve_export_filename("transcript", "all")) #type: ignore
        else:
            annotation.export.transcripts(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, rna_classes=rna_classes, filename=resolve_export_filename("transcript", "main")) #type: ignore

    if strip_stop is not None:
        protein_strip = strip_stop
        cds_strip = strip_stop
    else:
        protein_strip = not keep_stop
        cds_strip = strip_stop_cds

    if "protein" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        elif "CDS" in feature_id:
            used_id = "CDS"
        else:
            used_id = "protein"

        if "unique_per_gene" in mode:
            annotation.export.proteins(only_main=False, output_dir=output_dir, verbose=detailed_headers, unique_proteins_per_gene=True, used_id=used_id, strip_stop=protein_strip, filename=resolve_export_filename("protein", "unique_per_gene"))
        elif "unique" in mode:
            annotation.export.unique_proteins(output_dir=output_dir, quiet=quiet, strip_stop=protein_strip, filename=resolve_export_filename("protein", "unique"))
        elif "all" in mode:
            annotation.export.proteins(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False, strip_stop=protein_strip, filename=resolve_export_filename("protein", "all"))
        else:
            annotation.export.proteins(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, strip_stop=protein_strip, filename=resolve_export_filename("protein", "main"))

    if "CDS" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "CDS"

        protein_oriented = not raw_cds

        if "unique_per_gene" in mode:
            annotation.export.CDSs(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, unique_CDSs_per_gene=True, protein_oriented=protein_oriented, strip_stop=cds_strip, filename=resolve_export_filename("CDS", "unique_per_gene"))
        elif "unique" in mode:
            annotation.export.unique_CDSs(output_dir=output_dir, quiet=quiet, protein_oriented=protein_oriented, strip_stop=cds_strip, filename=resolve_export_filename("CDS", "unique"))
        elif "all" in mode:
            annotation.export.CDSs(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False, protein_oriented=protein_oriented, strip_stop=cds_strip, filename=resolve_export_filename("CDS", "all"))
        else:
            annotation.export.CDSs(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, protein_oriented=protein_oriented, strip_stop=cds_strip, filename=resolve_export_filename("CDS", "main"))

    if "promoter" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "promoter"

        if "all" in mode or "unique_per_gene" in mode or "unique" in mode:
            annotation.export.promoters(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, promoter_type=promoter_type, promoter_size=promoter_size, quiet=quiet, filename=resolve_export_filename("promoter", "all"))
        else:
            annotation.export.promoters(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, promoter_type=promoter_type, promoter_size=promoter_size, quiet=quiet, filename=resolve_export_filename("promoter", "main"))

if __name__ == "__main__":
    app()