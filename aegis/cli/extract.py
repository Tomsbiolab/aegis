import typer
import os

from typing import List, Optional
from typing_extensions import Annotated

from ..genome import Genome
from ..annotation import Annotation
from ..utils.genefunctions import TAXONOMY_ORGANELLE_CODES, resolve_taxonomy_tables, NCBI_GENETIC_CODES
from .utils import split_callback

FEATURES = ["gene", "transcript", "CDS", "protein", "promoter"]

VALID_IDS = ["gene", "transcript", "CDS", "feature"]

EXTRACTION_MODES = ["all", "main", "unique", "unique_per_gene"]

PROMOTER_TYPES = ["standard", "upstream_ATG", "standard_plus_up_to_ATG"]

RNA_CLASSES = ["mRNA", "antisense_lncRNA", "antisense_RNA", 
                "miRNA_primary_transcript", "ncRNA", "lncRNA",
                "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
                "SRP_RNA", "RNase_MRP_RNA"]

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="A name or tag for the genome assembly (e.g., 'TAIR10'). [default: a name derived from the genome FASTA filename]"
    )] = "{genome-file}",
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="A name or tag for the annotation version (e.g., 'Araport11'). [default: a name derived from the annotation filename]"
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the directory where output FASTA files will be saved."
    )] = "./aegis_output/features/",
    features: Annotated[List[str], typer.Option(
        "-f", "--features", help=f"Feature type(s) to extract, as a comma-separated list. Available options: {', '.join(FEATURES)}.",
        callback=split_callback
    )] = ["gene"],

    mode: Annotated[List[str], typer.Option(
        "-m", "--mode", help=f"""Extraction mode(s), as a comma-separated list. Controls filtering of features.\n\n

        - 'all': Extract all features (e.g., all transcripts for a gene).\n
        - 'unique_per_gene': Keep one copy of each unique protein/CDS sequence per gene.\n
        - 'main': Extract only the main variant (e.g., the longest transcript).\n
        - 'unique': Keep only one copy of each unique protein/CDS sequence across the entire output.""",
        callback=split_callback
    )] = ["all", "main"],
    rna_classes: Annotated[List[str], typer.Option(
        "-r", "--rna-classes", help=f"Filter transcripts by biotype (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included.",
        callback=split_callback
    )] = [],
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help=f"Size of the promoter region in base pairs (bp). Used only if 'promoter' is a selected feature."
    )] = 2000,

    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help="""\
                Defines the reference point for extracting promoter regions. Used only if 'promoter' is selected.\n\n
                
                Options:\n
                - 'standard': Upstream of the transcript's start site (TSS).\n
                - 'upstream_ATG': Upstream of the main CDS's start codon (ATG). Falls back to 'standard' if no CDS is present.\n
                - 'standard_plus_up_to_ATG': The 'standard' promoter plus the 5' UTR (sequence between TSS and ATG). Falls back to 'standard' if no CDS.
                """
    )] = "standard",

    detailed_headers: Annotated[bool, typer.Option(
        "-dh", "--detailed-headers", help=f"Add extra details in fasta headers; scaffold/chromosome number, genome co-ordinates, and/or protein tags if applicable."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons."
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments."
    )] = False,
    feature_id: Annotated[str, typer.Option(
        "--feature-id", help=f"Specifies which feature ID to use in FASTA headers. E.g., use 'gene' to label all outputs (transcripts, proteins) with their parent gene ID. 'feature' uses the most specific ID available. Available: {', '.join(VALID_IDS)}."
    )] = "feature",
    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID')."
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)')."
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers."
    )] = False,
    raw_cds: Annotated[bool, typer.Option(
        "--raw-cds", help="Export raw spliced genomic CDS sequences instead of in-frame protein-oriented coding sequences."
    )] = False,
    adjust_internal_shifts: Annotated[str, typer.Option(
        "--adjust-internal-shifts", help="Frameshift / phase handling mode for multi-segment CDS translation: 'intra_exon' (default: adjust shifts only across contiguous/overlapping exon segments, preserving continuous splicing across introns), 'all' (adjust across all junctions), or 'none' (translate continuous spliced sequence)."
    )] = "intra_exon",
    taxonomy: Annotated[str, typer.Option(
        "-tax", "--taxonomy", help="Taxonomic group preset for organelle genetic codes: 'plant' (default: nuclear=1, mito=1, plastid=11), 'vertebrate' (nuclear=1, mito=2), 'invertebrate' (nuclear=1, mito=5), or 'yeast' (nuclear=1, mito=3). Specific codes can be individually customized with --genetic-code, --mito-code, or --plastid-code."
    )] = "plant",
    genetic_code: Annotated[int, typer.Option(
        "-gc", "--genetic-code", help="NCBI genetic code table number to use for nuclear translation (default: 1)."
    )] = 1,
    auto_organelle_codes: Annotated[bool, typer.Option(
        "--auto-organelle-codes/--no-auto-organelle-codes", help="Automatically use mitochondrial and plastid genetic codes for organelle contigs (default: True)."
    )] = True,
    mito_code: Annotated[Optional[int], typer.Option(
        "-mc", "--mito-code", help="NCBI genetic code table number for mitochondrial contigs (overrides --taxonomy default)."
    )] = None,
    plastid_code: Annotated[Optional[int], typer.Option(
        "-pc", "--plastid-code", help="NCBI genetic code table number for plastid/chloroplast contigs (overrides --taxonomy default)."
    )] = None,
    mitochondria_chroms: Annotated[List[str], typer.Option(
        "--mitochondria-chroms", help="Explicit list or comma-separated names of mitochondrial chromosomes/scaffolds to translate with --mito-code.",
        callback=split_callback
    )] = [],
    chloroplast_chroms: Annotated[List[str], typer.Option(
        "--chloroplast-chroms", help="Explicit list or comma-separated names of chloroplast/plastid chromosomes/scaffolds to translate with --plastid-code.",
        callback=split_callback
    )] = [],
    polish_coordinates: Annotated[bool, typer.Option(
        "--polish-coordinates/--skip-coordinate-polishing", help="Mutate feature coordinates when boundaries differ (default: False, preserves original coordinates)."
    )] = False,
    strip_stop: Annotated[Optional[bool], typer.Option(
        "--strip-stop/--no-strip-stop", help="Explicitly enable/disable stripping stop codons across both proteins and CDSs (overrides defaults)."
    )] = None,
    keep_stop: Annotated[bool, typer.Option(
        "--keep-stop", help="Keep trailing stop codon (* in proteins) when exporting protein sequences (by default trailing stop codons are stripped)."
    )] = False,
    strip_stop_cds: Annotated[bool, typer.Option(
        "--strip-stop-cds", help="Strip trailing stop codon (terminal 3-nt stop codon) when exporting CDS sequences (by default CDS sequences retain the stop codon)."
    )] = False,
    rework_all_CDSs: Annotated[bool, typer.Option(
        "--rework-all-CDSs", help="Recalculates ALL CDSs from the genome sequence, overriding existing ones. More aggressive than --infer-missing-CDSs."
    )] = False,
    infer_missing_CDSs: Annotated[bool, typer.Option(
        "--infer-missing-CDSs", help="Detects and creates CDSs where missing, without overriding existing CDS annotations."
    )] = False,
    fallback_to_trim: Annotated[bool, typer.Option(
        "--fallback-to-trim", help="When recalculating CDSs, fallback to trimming unaligned ends if no high-ratio ORF is found."
    )] = False,
    coding_ratio_threshold: Annotated[float, typer.Option(
        "--coding-ratio-threshold", help="Threshold ratio of coding sequence length to transcript length for rework CDS (default: 0.7)."
    )] = 0.7,
    allow_internal_stops: Annotated[bool, typer.Option(
        "--allow-internal-stops/--no-allow-internal-stops", help="Allow internal stop codons in progressive rework fallback (default: True)."
    )] = True,
    allow_partial: Annotated[bool, typer.Option(
        "--allow-partial/--no-allow-partial", help="Allow partial ORFs without stop codon in progressive rework fallback (default: True)."
    )] = True,
    enforce_start_codon: Annotated[bool, typer.Option(
        "--enforce-start-codon/--no-enforce-start-codon", help="Require start codon (ATG) in initial rework passes (default: True)."
    )] = True,
    orf_choice_mode: Annotated[str, typer.Option(
        "--orf-choice-mode", help="ORF selection criteria: 'longest' or 'earliest' (default: 'longest')."
    )] = "longest",
    min_codon_len: Annotated[int, typer.Option(
        "--min-codon-len", help="Minimum codon length required for predicted ORFs (e.g. 30 or 50 to suppress micro-ORFs, default: 2)."
    )] = 2,
    recalculate_phases: Annotated[bool, typer.Option(
        "--recalculate-phases", help="Recalculate CDS segment phases based on segment lengths and splicing leftover, preserving 5' initial phase for partial CDSs."
    )] = False,
    reset_phases_zero: Annotated[bool, typer.Option(
        "--reset-phases-zero", help="Reset initial CDS phase to 0 and recalculate all downstream segment phases."
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

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_CDSs
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

    mito_chroms = mitochondria_chroms if len(mitochondria_chroms) > 0 else None
    chloro_chroms = chloroplast_chroms if len(chloroplast_chroms) > 0 else None

    genome = Genome(
        name=genome_name,
        genome_file_path=genome_file,
        quiet=quiet,
        header_id_tag=header_id_tag if header_id_tag != "" else None,
        header_id_regex=header_id_regex if header_id_regex != "" else None,
        gwh=gwh,
    )
    res_gc, res_mito, res_plastid = resolve_taxonomy_tables(
        taxonomy=taxonomy,
        table=genetic_code,
        mito_table=mito_code,
        plastid_table=plastid_code,
    )

    annotation = Annotation(
        name=annotation_name,
        annot_file_path=annotation_file,
        genome=genome,
        rework_all_CDSs=rework_all_CDSs,
        work_out_missing_CDSs=infer_missing_CDSs,
        fallback_to_trim=fallback_to_trim,
        quiet=quiet,
        collapse_exons=collapse_exons,
        collapse_CDSs=collapse_CDSs,
        adjust_internal_shifts=adjust_internal_shifts,
        taxonomy=taxonomy,
        table=res_gc,
        auto_organelle_codes=auto_organelle_codes,
        mito_table=res_mito,
        plastid_table=res_plastid,
        mitochondria_chroms=mito_chroms,
        chloroplast_chroms=chloro_chroms,
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
    if not has_any_cds and not infer_missing_CDSs and not rework_all_CDSs:
        typer.secho(
            "Notice: The input annotation does not contain any annotated CDS features. No protein or CDS sequences were extracted.\n"
            "Tip: Pass '--infer-missing-CDSs' to automatically detect and predict CDSs across transcripts.",
            fg=typer.colors.YELLOW,
            err=True,
        )

    if "gene" in features:

        annotation.export.genes(output_dir=output_dir, verbose=detailed_headers)

    if "transcript" in features:

        if "gene" in feature_id:
            used_id = "gene"
        else:
            used_id = "transcript"

        if "unique_per_gene" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, output_dir=output_dir, used_id=used_id, rna_classes=rna_classes, unique_transcripts_per_gene=True) #type: ignore
        elif "unique" in mode:
            annotation.export.unique_transcripts(output_dir=output_dir, quiet=quiet, rna_classes=rna_classes) #type: ignore
        elif "all" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, output_dir=output_dir, used_id=used_id, rna_classes=rna_classes) #type: ignore
        else:
            annotation.export.transcripts(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, rna_classes=rna_classes) #type: ignore

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
            annotation.export.proteins(only_main=False, output_dir=output_dir, verbose=detailed_headers, unique_proteins_per_gene=True, used_id=used_id, strip_stop=protein_strip)
        elif "unique" in mode:
            annotation.export.unique_proteins(output_dir=output_dir, quiet=quiet, strip_stop=protein_strip)
        elif "all" in mode:
            annotation.export.proteins(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False, strip_stop=protein_strip)
        else:
            annotation.export.proteins(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, strip_stop=protein_strip)

    if "CDS" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "CDS"

        protein_oriented = not raw_cds

        if "unique_per_gene" in mode:
            annotation.export.CDSs(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, unique_CDSs_per_gene=True, protein_oriented=protein_oriented, strip_stop=cds_strip)
        elif "unique" in mode:
            annotation.export.unique_CDSs(output_dir=output_dir, quiet=quiet, protein_oriented=protein_oriented, strip_stop=cds_strip)
        elif "all" in mode:
            annotation.export.CDSs(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False, protein_oriented=protein_oriented, strip_stop=cds_strip)
        else:
            annotation.export.CDSs(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, protein_oriented=protein_oriented, strip_stop=cds_strip)

    if "promoter" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "promoter"

        if "all" in mode or "unique_per_gene" in mode or "unique" in mode:
            annotation.export.promoters(only_main=False, output_dir=output_dir, verbose=detailed_headers, used_id=used_id, promoter_type=promoter_type, promoter_size=promoter_size, quiet=quiet)
        else:
            annotation.export.promoters(output_dir=output_dir, verbose=detailed_headers, used_id=used_id, promoter_type=promoter_type, promoter_size=promoter_size, quiet=quiet)

if __name__ == "__main__":
    app()