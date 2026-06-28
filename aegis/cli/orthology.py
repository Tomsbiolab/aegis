import typer
import os
import pandas as pd
import numpy as np
import warnings
import shutil
import sys
import re

from pathlib import Path
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from ..feature import Feature
from ..equivalence import Simple_annotation, pairwise_orthology, run_command
from time import time

app = typer.Typer(add_completion=False, no_args_is_help=True)

def parse_score_column(score_str):
    if pd.isna(score_str) or score_str == "NA":
        return {
            "Liftoff": "NA", "LiftOn": "NA", "AEGIS_Overlap": "NA", 
            "MCscan": "NA", "BLAST": "NA", "OrthoFinder": "NA"
        }
    
    res = {"Liftoff": [], "LiftOn": [], "AEGIS_Overlap": [], "MCscan": [], "BLAST": [], "OrthoFinder": []}
    
    pattern = re.compile(r'([a-zA-Z0-9_]+)\s+\[(.*?)\]')
    matches = pattern.findall(str(score_str))
    
    for tool, val in matches:
        entry = f"{tool} [{val}]"
        if "liftoff" in tool:
            res["Liftoff"].append(entry)
        elif "lifton" in tool:
            res["LiftOn"].append(entry)
        elif tool == "aegis_overlap":
            res["AEGIS_Overlap"].append(entry)
        elif "mcscan" in tool:
            res["MCscan"].append(entry)
        elif "blast" in tool or "rbh" in tool or "rbbh" in tool:
            res["BLAST"].append(entry)
        elif "orthofinder" in tool:
            res["OrthoFinder"].append(entry)

    return {k: (", ".join(v) if v else "NA") for k, v in res.items()}

def get_tiered_cardinality(df, allowed_scores):
    if df.empty:
        return pd.Series(dtype=str)
    
    is_ordered = df['annotation_A'] <= df['annotation_B']
    
    base_df = pd.DataFrame({
        'a1': np.where(is_ordered.to_numpy(), df['annotation_A'].to_numpy(), df['annotation_B'].to_numpy()),
        'a2': np.where(is_ordered.to_numpy(), df['annotation_B'].to_numpy(), df['annotation_A'].to_numpy()),
        'g1': np.where(is_ordered.to_numpy(), df['gene_id_A'].to_numpy(), df['gene_id_B'].to_numpy()),
        'g2': np.where(is_ordered.to_numpy(), df['gene_id_B'].to_numpy(), df['gene_id_A'].to_numpy()),
        'score': df['summary_score']
    })
    
    sub_df = base_df[base_df['score'].isin(allowed_scores)].copy()
    sub_df = sub_df.drop_duplicates(subset=['a1', 'a2', 'g1', 'g2'])
    
    if sub_df.empty:
        return pd.Series(['NA'] * len(df), index=df.index)

    c1 = sub_df.groupby(['a1', 'a2', 'g1'])['g2'].count().reset_index(name='d1')
    c2 = sub_df.groupby(['a1', 'a2', 'g2'])['g1'].count().reset_index(name='d2')
    
    t1 = sub_df.merge(c2, on=['a1', 'a2', 'g2']).groupby(['a1', 'a2', 'g1'])['d2'].max().reset_index(name='max_d2')
    
    sub_df = sub_df.merge(c1, on=['a1', 'a2', 'g1'])
    sub_df = sub_df.merge(t1, on=['a1', 'a2', 'g1'])
    
    deg_A = sub_df['d1'].to_numpy()
    deg_B = sub_df['max_d2'].to_numpy()
    
    condlist = [
        (deg_A == 1) & (deg_B == 1),
        (deg_A > 1) & (deg_B == 1),
        (deg_A == 1) & (deg_B > 1),
        (deg_A > 1) & (deg_B > 1)
    ]
    
    sub_df['label_ordered'] = np.select(condlist, ['1:1', '1:N', 'N:1', 'N:N'], default='NA')
    sub_df['label_reversed'] = np.select(condlist, ['1:1', 'N:1', '1:N', 'N:N'], default='NA')
    
    sub_df = sub_df.drop_duplicates(subset=['a1', 'a2', 'g1', 'g2']).set_index(['a1', 'a2', 'g1', 'g2'])
    base_df = base_df.join(sub_df[['label_ordered', 'label_reversed']], on=['a1', 'a2', 'g1', 'g2'])
    
    final_labels = np.where(
        is_ordered.to_numpy(),
        base_df['label_ordered'].fillna('NA'),
        base_df['label_reversed'].fillna('NA')
    )
    
    return pd.Series(final_labels, index=df.index)

def flip_score_string(score_str):
    if pd.isna(score_str) or score_str == "NA":
        return score_str
    
    def replace_tool(match):
        tool = match.group(1)
        vals = match.group(2)
        
        if tool.startswith("fwd_"):
            tool = tool.replace("fwd_", "rev_", 1)
        elif tool.startswith("rev_"):
            tool = tool.replace("rev_", "fwd_", 1)
            
        if "rbh" in tool or "rbbh" in tool:
            new_vals_list = []
            for metric in vals.split(", "):

                if "=" in metric:
                    k, v = metric.split("=", 1)
                    if "/" in v:
                        parts = v.split("/")
                        if len(parts) == 2:
                            v = f"{parts[1]}/{parts[0]}"
                    new_vals_list.append(f"{k}={v}")
                else:
                    new_vals_list.append(metric)
            vals = ", ".join(new_vals_list)
        
        elif tool.startswith("rec_"):
            if "/" in vals:
                parts = vals.split("/")
                if len(parts) == 2:
                    vals = f"{parts[1]}/{parts[0]}"

        return f"{tool} [{vals}]"
        
    return re.sub(r'([a-zA-Z0-9_]+)\s+\[(.*?)\]', replace_tool, str(score_str))

def flip_masked_rows(df, mask):
    if not mask.any():
        return df
        
    df.loc[mask, ['gene_id_A', 'gene_id_B']] = df.loc[mask, ['gene_id_B', 'gene_id_A']].values
    df.loc[mask, ['annotation_A', 'annotation_B']] = df.loc[mask, ['annotation_B', 'annotation_A']].values
    
    if "species_A" in df.columns and "species_B" in df.columns:
        df.loc[mask, ['species_A', 'species_B']] = df.loc[mask, ['species_B', 'species_A']].values
        
    for col in ['cardinality', 'cardinality_strict', 'cardinality_moderate', 'cardinality_relaxed']:
        if col in df.columns:
            df.loc[mask, col] = df.loc[mask, col].replace({'1:N': 'tmp', 'N:1': '1:N'}).replace({'tmp': 'N:1'})
            
    if 'score' in df.columns:
        df.loc[mask, 'score'] = df.loc[mask, 'score'].apply(flip_score_string)
        
    return df

def merge_score_strings(series):
    """Combines the score strings of identical edges so no tool evidence is lost."""
    combined = set()
    for s in series:
        if pd.isna(s) or s == "NA":
            continue
        # Extract every "tool_name [values]" block and add it to a set
        for match in re.finditer(r'([a-zA-Z0-9_]+)\s+\[(.*?)\]', str(s)):
            combined.add(match.group(0))
    return ", ".join(sorted(list(combined))) if combined else "NA"

def best_summary_score(series):
    """When merging twins, keeps the highest confidence tier found."""
    vals = set(series.dropna().values)
    if "high_confidence" in vals: return "high_confidence"
    if "medium_confidence" in vals: return "medium_confidence"
    if "lower_confidence" in vals: return "lower_confidence"
    return "NA"

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_files: Annotated[list[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    
    # ==========================================
    # CORE INPUT / OUTPUT CONFIGURATION
    # ==========================================
    genome_files: Annotated[str, typer.Option(
        "-g", "--genome-files", 
        help="Genome assemblies corresponding to annotation files. Provide them in the same number and order, separated by commas. e.g. -g genomefile1,genomefile2,genomefile3,genomefile4",
        callback=split_callback,
        rich_help_panel="Core Input/Output Configuration"
    )],
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", 
        help="Annotation versions, names or tags otherwise they will just be the annotation file basename without the extension. Provide them in the same number and order as the corresponding annotation files, separated by commas. e.g. --annotation-names name1,name2,name3,name4",
        callback=split_callback,
        rich_help_panel="Core Input/Output Configuration"
    )] = "{annotation-filename(s)}",
    genome_names: Annotated[str, typer.Option(
        "--genome-names", 
        help="Genome versions, names or tags otherwise they will just be the genome file basename without the extension. Provide them in the same number and order as the corresponding genome files, separated by commas. e.g. --genome-names name1,name2,name3,name4",
        callback=split_callback,
        rich_help_panel="Core Input/Output Configuration"
    )] = "{genome-filename(s)}",
    group_names: Annotated[str, typer.Option(
        "-gn", "--group-names", 
        help="Optional grouping of input annotations, into species for example. Use NA as a placemarker for annotation files without a group label. e.g. --group-names group1,NA,group1,group2",
        callback=split_callback,
        rich_help_panel="Core Input/Output Configuration"
    )] = "",
    reference_annotation: Annotated[str, typer.Option(
        "--reference-annotation", 
        help="Select a single annotation, by providing its name/tag or filename, to use as a reference. Only matches to and from this annotation will be reported. Otherwise matches are reported between all annotations.",
        rich_help_panel="Core Input/Output Configuration"
    )] = "None",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", 
        help="Path to the output folder.",
        rich_help_panel="Core Input/Output Configuration"
    )] = "./aegis_output/orthologues/",
    output_filename: Annotated[str, typer.Option(
        "-o", "--output-file", 
        help="Output filename to be saved to output folder without extension. The '.tsv' extension will be added to the filename.",
        rich_help_panel="Core Input/Output Configuration"
    )] = "equivalences{other_tags}.tsv",

    # ==========================================
    # Orthology Tool Options
    # ==========================================
    skip_liftoff: Annotated[bool, typer.Option(
        "--skip-liftoff", 
        help="Skip Liftoff.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_lifton: Annotated[bool, typer.Option(
        "--skip-lifton", 
        help="Skip LiftOn.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_copies: Annotated[bool, typer.Option(
        "--skip-copies", 
        help="Liftoff and LiftOn are run in copies mode by default, flag to deactivate.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_mcscan: Annotated[bool, typer.Option(
        "--skip-mcscan", 
        help="Skip the JCVI toolkit synteny and collinearity analysis (MCscan). Useful when JCVI is causing compatibility issues.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_synteny: Annotated[bool, typer.Option(
        "--skip-synteny", 
        help="Skip conservation of synteny metrics whenever an annotation is lifted over to another genome.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_orthofinder: Annotated[bool, typer.Option(
        "--skip-orthofinder", 
        help="Skip the OrthoFinder analysis.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    pairwise_orthofinder: Annotated[bool, typer.Option(
        "--pairwise-orthofinder", 
        help="Execute OrthoFinder on independent annotation pairs. Overrides the default multi-annotation bulk analysis. Recommended for highly divergent taxa or targeted 1:1 orthologue mapping.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    lift_feature_types: Annotated[str, typer.Option(
        "--lift-feature-types", 
        help="All feature types within an annotation files are lifted over by default, however a more restrictive set can be used, separated by commas, such as 'gene,mRNA,exon,CDS,pseudogene,pseudogenic_exon,pseudogenic_transcript'.", 
        callback=split_callback,
        rich_help_panel="Orthology Tool Options"
    )] = "ALL",
    include_single_blasts: Annotated[bool, typer.Option(
        "-b", "--include-single-blasts", 
        help="Decide whether to report unidirectional (i.e. just fw or rv) blasts in the orthologue summary.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_rbhs: Annotated[bool, typer.Option(
        "--skip-RBHs", 
        help="Decide whether to skip RBHs which are not RBBHs, these are reported by default in the orthologue summary.",
        rich_help_panel="Orthology Tool Options"
    )] = False,
    skip_all_blasts: Annotated[bool, typer.Option(
        "--skip-all-blasts", 
        help="Skip all BLASTs.",
        rich_help_panel="Orthology Tool Options"
    )] = False,

    # ==========================================
    # BLAST Options
    # ==========================================
    identity: Annotated[float, typer.Option(
        "-i", "--identity", 
        help="Minimum identity threshold for BLAST hits.",
        rich_help_panel="BLAST Options"
    )] = 30.0,
    coverage: Annotated[float, typer.Option(
        "-c", "--coverage", 
        help="Minimum coverage threshold for BLAST hits.",
        rich_help_panel="BLAST Options"
    )] = 30.0,
    evalue: Annotated[float, typer.Option(
        "-e", "--evalue", 
        help="Maximum e-value threshold for BLAST hits.",
        rich_help_panel="BLAST Options"
    )] = 0.00001,

    # ==========================================
    # Output Options
    # ==========================================
    confidence: Annotated[str, typer.Option(
        "--confidence", 
        help="Filter the final output by confidence levels. Options: high, medium, lower. Separate by commas.",
        callback=split_callback,
        rich_help_panel="Output Options"
    )] = "high,medium,lower",
    include_NAs: Annotated[bool, typer.Option(
        "-na", "--include-NAs", 
        help="Append all genes that have no equivalences (or were filtered out) at the end of the output.",
        rich_help_panel="Output Options"
    )] = False,
    skip_cardinality: Annotated[bool, typer.Option(
        "-sc", "--skip-cardinality", 
        help="Skip the cardinality analysis (which marks gene pairs as 1:N, N:1, N:N, or 1:1) in the final output table (after the confidence level filter).",
        rich_help_panel="Output Options"
    )] = False,
    tiered_cardinality: Annotated[bool, typer.Option(
        "--tiered-cardinality", 
        help="Report three separate cardinality columns (strict: just looking at high-confidence orthologues, moderate: high- and medium-confidence orthologues, relaxed: high-, medium- and lower-confidence orthologues).",
        rich_help_panel="Output Options"
    )] = False,
    include_duplicates: Annotated[bool, typer.Option(
        "--include-duplicates", 
        help="Report equivalences from both from gene_id_A to gene_id_B as well as from gene_id_B to gene_id_A. These 'duplicate gene pairs' are not included by default.",
        rich_help_panel="Output Options"
    )] = False,
    split_scores: Annotated[bool, typer.Option(
        "--split-scores", 
        help="Split the aggregated 'score' column into individual columns for Liftoff, LiftOn, Overlap, MCscan, BLAST, and OrthoFinder.",
        rich_help_panel="Output Options"
    )] = False,

    # ==========================================
    # Execution/Debugging
    # ==========================================
    threads: Annotated[int, typer.Option(
        "-t", "--threads", 
        help="Number of threads.",
        rich_help_panel="Execution/Debugging"
    )] = 1,
    parallel_pairs: Annotated[bool, typer.Option(
        "--parallel-pairs", 
        help="Run independent pairwise comparisons concurrently. Strongly recommended for multiple genomes to optimize CPU scaling. The threads parameter will be distributed across the pairwise comparisons.",
        rich_help_panel="Execution/Debugging"
    )] = False,
    keep_intermediate: Annotated[bool, typer.Option(
        "-k", "--keep-intermediate", 
        help="Keep intermediate files, useful for identifying errors.",
        rich_help_panel="Execution/Debugging"
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", 
        help="Verbose logging, useful if encountering a problem or error.",
        rich_help_panel="Execution/Debugging"
    )] = False,
):
    """
    Provides a set of orthologues relationships leveraging external and internal tools such as Litfoff + AEGIS overlaps, LiftOn + AEGIS overlaps, MCscan, orthofinder and BLAST. Wherever relevant, tools are run reciprocally for an extra confidence mark in orthologous relationships.

    Script to perform gene correspondence analysis between multiple genomes.

    This script automates a bioinformatics pipeline that uses several tools
    (Liftoff, Lifton, DIAMOND, JCVI) to find orthologous genes between all possible
    pairs from a given list of genomes. For each pair, it performs:

    1.  Annotation Liftover: Transfers annotations between the two genomes based
        on sequence homology (using Liftoff and Lifton).
    2.  Reciprocal Best Hit (RBH): Uses DIAMOND (a fast alternative to BLAST) to
        find the best reciprocal homologs at the protein level.
    3.  Synteny and Collinearity: Uses the JCVI toolkit to identify conserved gene
        blocks (synteny) and find orthologs within that context.

    After all pairwise comparisons, it collects all unique proteomes and runs
    OrthoFinder to infer orthogroups across all species.

    The script is designed to be scalable, processing any number of genomes provided
    in the configuration section. All results for each pairwise comparison are stored
    in a separate, clearly named directory.

    All external tool commands are executed via Docker to ensure a reproducible
    environment.
    """

    start = time()

    print(f"Running command: {' '.join(sys.argv)}")

    valid_confidence_levels = {"high", "medium", "lower", "high_confidence", "medium_confidence", "lower_confidence", "highest"}
    for conf in confidence:
        if conf.strip().lower() not in valid_confidence_levels:
            raise typer.BadParameter(
                f"Invalid confidence level: '{conf}'. Allowed options are 'high', 'medium', 'lower', 'high_confidence', 'medium_confidence', or 'lower_confidence' (separated by commas)."
            )

    expected_extensions = {".gff", ".gtf", ".gff3", ".fasta", ".fa"}
    files_to_check = []
    
    if annotation_files and annotation_names == ["{annotation-filename(s)}"]:
        files_to_check.extend(annotation_files)
    if genome_files and genome_names == ["{genome-filename(s)}"]:
        files_to_check.extend(genome_files)
        
    invalid_extensions = []
    for f in files_to_check:
        _, ext = os.path.splitext(f)
        if ext.lower() not in expected_extensions:
            invalid_extensions.append(f)
            
    if invalid_extensions:
        print("\nWarning: Proper file extensions were not found for the following file(s):")
        for f in invalid_extensions:
            print(f"  - {os.path.basename(f)}")
        print("Expected extensions: .gff, .gtf, .gff3, .fasta, .fa\n"
              "Some naming or downstream processing functionality may be compromised if extensions are missing.\n")

    quiet=not(verbose)
    synteny=not(skip_synteny)
    
    if len(annotation_files) < 2:
        raise typer.BadParameter(f"At least 2 annotation-files must be provided.")
    
    if len(annotation_files) > 1 and annotation_files[-1].lower() in ("true", "false"):
        typer.echo(
            "⚠️  Detected extra value 'true' or 'false' at the end of positional arguments.\n"
            "👉 Did you mean to use the '--include_NAs' or '--simple' flags? Use them like this: '-n' or '-s' (no 'true' needed).",
            err=True,
        )
        raise typer.Exit(code=1)

    if annotation_names == ["{annotation-filename(s)}"]:
        annotation_names = [] # type: ignore
        for annotation_file in annotation_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0]) # type: ignore

    if genome_names == ["{genome-filename(s)}"]:
        genome_names = [] # type: ignore
        for genome_file in genome_files:
            genome_names.append(os.path.splitext(os.path.basename(genome_file))[0]) # type: ignore

    genome_files = [str(Path(g).resolve()) for g in genome_files] # type: ignore
    annotation_files = [str(Path(g).resolve()) for g in annotation_files] # type: ignore

    if len(annotation_files) != len(annotation_names):
        raise typer.BadParameter(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if len(genome_files) != len(genome_names):
        raise typer.BadParameter(f"The provided number of genome name(s)/tag(s) do not match the number of genome file(s).")

    name_to_files = {}
    file_to_names = {}
    for g_file, g_name in zip(genome_files, genome_names):
        name_to_files.setdefault(g_name, set()).add(g_file)
        file_to_names.setdefault(g_file, set()).add(g_name)

    for name, files in name_to_files.items():
        if len(files) > 1:
            base_files = [os.path.basename(f) for f in files]
            raise typer.BadParameter(
                f"The genome tag '{name}' is mapped to multiple distinct genome files: {base_files}. "
                "Each distinct assembly file must have a unique tag."
            )

    for g_file, names in file_to_names.items():
        if len(names) > 1:
            raise typer.BadParameter(
                f"The genome file '{os.path.basename(g_file)}' is mapped to multiple different tags: {list(names)}. "
                "The same assembly file must use a single consistent tag."
            )

    for annotation_name in annotation_names:
        if "__to__" in annotation_name:
            raise typer.BadParameter(f"The provided annotation name/tag '{annotation_name}' has an incompatible term: '__to__' as it is used internally for temporary file naming.")

    for genome_name in genome_names:
        if "__to__" in genome_name:
            raise typer.BadParameter(f"The provided genome name/tag '{genome_name}' has an incompatible term: '__to__' as it is used internally for temporary file naming.")

    if skip_mcscan and skip_liftoff and skip_all_blasts and skip_orthofinder and skip_lifton:
        raise typer.BadParameter("No analysis methods selected. Please select at least one analysis method. Run 'aegis orthology --help' for more information.")

    if len(annotation_names) != len(set(annotation_names)):
        raise typer.BadParameter("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annotation_files) != len(set(annotation_files)):
        raise typer.BadParameter("Avoid repeated annotation filename(s).")

    if len(genome_files) != len(annotation_files):
        raise typer.BadParameter(f"A single genome file must be provided for each annotation file. Genome files: {genome_files}. Annotation files: {annotation_files}")

    if len(genome_names) != len(set(genome_names)):
        repeated_genomes = {}
        for idx, genome in enumerate(genome_names):
            repeated_genomes.setdefault(genome, []).append(annotation_names[idx])
            
        repeated_msg = []
        for genome, a_names in repeated_genomes.items():
            if len(a_names) > 1:
                repeated_msg.append(f"  - Genome '{genome}' used by annotations: {', '.join(a_names)}")
                
        print("\nWarning: Repeated genome assemblies.\n"
              "\nBear in mind that annotations associated to the same genome will only be compared at the level of 'aegis overlap' and BLAST results. Also, for this reason, OrthoFinder will only be run in pairwise mode, regardless of the flags used.\n"
              "\nRepeated assignments found:\n" + "\n".join(repeated_msg) + "\n")
        
        pairwise_orthofinder = True
    
    if group_names != []:
        if len(annotation_files) != len(group_names):
            raise typer.BadParameter(f"The provided number of groups do not match the number of annotation file(s).")
        
    else:
        group_names = ["NA"] * len(annotation_files) # type: ignore

    if reference_annotation != "None":
        if reference_annotation not in annotation_files and reference_annotation not in annotation_names:
            raise typer.BadParameter(f"The provided reference-annotation = {reference_annotation} is not present neither in annotation-files ({annotation_files}) nor annotation-names ({annotation_names}).")

    skip_unidirectional_blasts = not (include_single_blasts)

    if skip_rbhs and not skip_unidirectional_blasts:
        raise typer.BadParameter(f"Do not include single blasts if rbhs are to be skipped as these provide higher support for orthology.")

    if parallel_pairs and threads == 1:
        print("\n⚠️  Warning: '--parallel-pairs' is enabled but '--threads' is set to 1.")
        print("   Under this configuration, tasks will still execute sequentially (one at a time)")
        print("   because only 1 worker can run. To run multiple pairwise comparisons concurrently,")
        print("   please set '--threads' (or '-t') to a value greater than 1 (e.g., -t 4).\n")

    genome_name_map = {path: name for path, name in zip(genome_files, genome_names)}    

    genomes:dict[str, Genome] = { g: Genome(name=genome_name_map[g], genome_file_path=g, quiet=quiet) for g in set(genome_files) }

    annotations:list[Annotation] = []

    for n, annotation_file in enumerate(annotation_files):

        annotations.append(Annotation(name=annotation_names[n], genome=genomes[genome_files[n]], annot_file_path=annotation_file, quiet=quiet))

        annotations[-1].rename_ids(strip_gene_tag=True, quiet=quiet)

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True

    output_dir_path = Path(output_dir).resolve()

    if output_dir_path.exists():
        if any(output_dir_path.iterdir()):
            warnings.warn(f"The folder '{output_dir_path}' already exists and is not empty. Conflict may arise with existing output and/or temp folder files.")

    output_dir_path.mkdir(parents=True, exist_ok=True)

    results_directory = output_dir_path / "temp"
    
    if results_directory.exists() and any(results_directory.iterdir()):
        print("\n" + "!" * 80)
        print("⚠️  NOTICE: Existing intermediate files detected in the temporary directory:")
        print(f"   {results_directory}")
        print("\n   The pipeline is configured to resume analysis by reusing matching files.")
        print("   At the user's discretion: if the previous run was aborted, cancelled,")
        print("   or crashed, some intermediate files may be incomplete or corrupted.")
        print("   This tool does not validate the internal completeness of pre-existing")
        print("   intermediates.")
        print("\n   👉 To guarantee a completely clean run, manually delete the 'temp/' folder")
        print("      before executing.")
        print("!" * 80 + "\n")

    protein_path = results_directory / "proteins"
    protein_path.mkdir(parents=True, exist_ok=True)

    CDS_path = results_directory / "CDSs"
    CDS_path.mkdir(parents=True, exist_ok=True)

    diamond_path = results_directory / "diamond"
    diamond_path.mkdir(parents=True, exist_ok=True)

    gff_path = results_directory / "gffs"
    gff_path.mkdir(parents=True, exist_ok=True)

    if not skip_lifton:
        lifton_path = results_directory / "lifton"
        lifton_path.mkdir(parents=True, exist_ok=True)

    liftoff_path = results_directory / "liftoff"
    liftoff_path.mkdir(parents=True, exist_ok=True)

    if not skip_mcscan:
        mcscan_path = results_directory / "mcscan"
        mcscan_path.mkdir(parents=True, exist_ok=True)

    liftless_overlaps_dir = results_directory / "litfless_overlaps"

    if len(genome_files) != len(set(genome_files)):
        liftless_overlaps_dir.mkdir(parents=True, exist_ok=True)

    if lift_feature_types == ["ALL"]:
        lift_feature_types = ["gene", "mRNA", "exon", "CDS", "pseudogene", "pseudogenic_exon", "pseudogenic_transcript"] # type: ignore
    
    lift_feature_types_file = results_directory / "chosen_liftover_features.txt"
    lift_feature_types_file = str(lift_feature_types_file)

    f_in = open(lift_feature_types_file, "w", encoding="utf-8")

    for ft in lift_feature_types:
        ft = ft.strip()
        f_in.write(f"{ft}\n")
    f_in.close()


    # Create gff, protein, CDS files, mcscan, and diamond databases in a non-redundant manner
    for n, a in enumerate(annotations):

        mcscan_name = a.name.replace(".", "_")

        gff_file1 = gff_path / f"{a.name}.gff3"
        gff_file2 = gff_path / f"{mcscan_name}.gff3"

        if not gff_file1.exists():
            a.export.gff(output_dir=str(gff_path), filename=f"{a.name}.gff3", subfolder=False, quiet=quiet)
        else:
            print(f"\n\tExisting GFF file for {a.name}. Skipping.")
        if not gff_file2.exists():
            a.export.gff(output_dir=str(gff_path), filename=f"{mcscan_name}.gff3", subfolder=False, quiet=quiet)
        else:
            print(f"\n\tExisting GFF mcscan file for {a.name}. Skipping.")
            
        if not skip_lifton:
            lifton_prep_file = gff_path / f"{a.name}__for__lifton.gff3"
            if not lifton_prep_file.exists():
                a_lifton = a.copy()
                a_lifton.CDS_to_CDS_segment_ids(override=True)
                a_lifton.export.gff(output_dir=str(gff_path), filename=f"{a_lifton.name}__for__lifton.gff3", subfolder=False, quiet=quiet)
                del a_lifton
            else:
                print(f"\n\tExisting lifton prep file for {a.name}. Skipping.")

        Feature._ACTIVE_GENOME = a.genome

        protein_fasta = protein_path / f"{a.name}_proteins_g_id_main.fasta"
        if not protein_fasta.exists():
            a.export.proteins(only_main=True, output_dir=str(protein_path), used_id="gene", verbose=False, filename=f"{a.name}_proteins_g_id_main.fasta")
        else:
            print(f"\n\tExisting protein fasta for {a.name}. Skipping.")
        a.clear_proteins()
        
        cds_fasta = CDS_path / f"{mcscan_name}_CDSs_g_id_main.fasta"
        if not cds_fasta.exists():
            a.export.CDSs(only_main=True, output_dir=str(CDS_path), used_id="gene", verbose=False, filename=f"{mcscan_name}_CDSs_g_id_main.fasta")
        else:
            print(f"\n\tExisting CDS fasta for {a.name}. Skipping.")


        if not skip_all_blasts:
            diamond_db_file = diamond_path / f"{a.name}_diamond_db.dmnd"
            if not diamond_db_file.exists():
                makedb_cmd = [
                    "diamond", "makedb", "-p", str(threads), "--in", str(protein_fasta), "--db", str(diamond_path / f"{a.name}_diamond_db")
                ]
                run_command(diamond_path, makedb_cmd)
            else:
                print(f"\n\tExisting DIAMOND database for {a.name}. Skipping.")

        if not skip_mcscan:
            cleaned_cds = mcscan_path / Path(f"{mcscan_name}.cds")
            bed_file = mcscan_path / Path(f"{mcscan_name}.bed")

            if not cleaned_cds.exists():
                jcvi_format_cmd_1 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta), str(cleaned_cds)]
                run_command(mcscan_path, jcvi_format_cmd_1)
            else:
                print(f"\n\tExisting cleaned CDS for {a.name}. Skipping.")

            if not bed_file.exists():
                gff_to_bed_cmd_1 = [
                    "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
                    "--key=Parent", "--primary_only", f"{gff_path}/{mcscan_name}.gff3", "-o", str(bed_file)
                ]
                run_command(mcscan_path, gff_to_bed_cmd_1)
            else:
                print(f"\n\tExisting BED file for {a.name}. Skipping.")


    tasks = []

    for n1, a1 in enumerate(annotations):

        for n2, a2 in enumerate(annotations):

            if n1 == n2:
                continue

            if reference_annotation != "None" and not a1.target and not a2.target:
                continue

            tasks.append((a1, a2, genomes[genome_files[n1]], genomes[genome_files[n2]]))

    if parallel_pairs and len(tasks) > 1:
        from concurrent.futures import ThreadPoolExecutor

        # Dynamically split total threads among concurrent workers
        num_workers = min(threads, len(tasks))
        threads_per_task = max(1, threads // num_workers)

        print(f"\nRunning {len(tasks)} pairwise comparisons in parallel using {num_workers} workers "
            f"({threads_per_task} thread(s) per worker)...")

        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            for ta1, ta2, tg1, tg2 in tasks:
                f = executor.submit(
                    pairwise_orthology,
                    annot1=ta1,
                    annot2=ta2,
                    genome1=tg1,
                    genome2=tg2,
                    working_directory=results_directory,
                    num_threads=threads_per_task,
                    copies=not(skip_copies),
                    synteny=synteny,
                    skip_liftoff=skip_liftoff,
                    skip_lifton=skip_lifton,
                    skip_mcscan=skip_mcscan,
                    types=lift_feature_types_file,
                    coverage=coverage,
                    evalue=evalue,
                    skip_blasts=skip_all_blasts,
                    pairwise_orthofinder=pairwise_orthofinder,
                    skip_orthofinder=skip_orthofinder,
                    quiet=True
                )
                futures.append(f)

            # Gather results to raise any exceptions that occurred
            for f in futures:
                f.result()
    else:
        # Standard sequential execution
        print(f"\nRunning pairwise comparisons sequentially ({threads} thread(s) per run)...")
        for ta1, ta2, tg1, tg2 in tasks:
            pairwise_orthology(
                annot1=ta1,
                annot2=ta2,
                genome1=tg1,
                genome2=tg2,
                working_directory=results_directory,
                num_threads=threads,
                copies=not(skip_copies),
                synteny=synteny,
                skip_liftoff=skip_liftoff,
                skip_lifton=skip_lifton,
                skip_mcscan=skip_mcscan,
                types=lift_feature_types_file,
                coverage=coverage,
                evalue=evalue,
                skip_blasts=skip_all_blasts,
                pairwise_orthofinder=pairwise_orthofinder,
                skip_orthofinder=skip_orthofinder,
                quiet=quiet
            )


    if not skip_all_blasts:
        # Obtaining RBHs and RBBHs from single blast results
        checked_pairs = []
        for n1, a1 in enumerate(annotations):
            for n2, a2 in enumerate(annotations):
                if n1 == n2:
                    continue
                pair = [n1, n2]
                pair.sort()
                if pair in checked_pairs:
                    continue
                checked_pairs.append(pair)

                rbh_out = diamond_path / f"rbh_{a1.name}__to__{a2.name}.txt"
                rbbh_out = diamond_path / f"rbbh_{a1.name}__to__{a2.name}.txt"


                if rbh_out.exists() and rbbh_out.exists():
                    if not quiet:
                        print(f"\n\tExisting RBH/RBBH tables found for {a1.name} and {a2.name}. Skipping generation.")
                    continue

                print(f"\nProcessing RBH and RBBHs for {a1.name} and {a2.name}")

                fwd_in = diamond_path / f"single_{a1.name}__to__{a2.name}.txt"
                rev_in = diamond_path / f"single_{a2.name}__to__{a1.name}.txt"
                fwd_best_in = diamond_path / f"single_best_{a1.name}__to__{a2.name}.txt"
                rev_best_in = diamond_path / f"single_best_{a2.name}__to__{a1.name}.txt"
                
                fwd_results = pd.read_csv(fwd_in, sep="\t", header=None)
                rev_results = pd.read_csv(rev_in, sep="\t", header=None)

                headers = ["query", "subject", "identity", "coverage", "qlength", "slength", "alength", "bitscore", "E-value"]
                fwd_results.columns = headers
                rev_results.columns = headers

                fwd_results = fwd_results[(fwd_results["identity"] >= identity)]
                rev_results = rev_results[(rev_results["identity"] >= identity)]

                # Create a new column in both dataframes: normalised bitscore
                fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
                rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

                # Create query and subject coverage columns in both dataframes
                fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
                rev_results['qcov'] = rev_results.alength/rev_results.qlength
                fwd_results['scov'] = fwd_results.alength/fwd_results.slength
                rev_results['scov'] = rev_results.alength/rev_results.slength

                # Clip maximum coverage values at 1.0
                fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
                rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
                fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
                rev_results['scov'] = rev_results['scov'].clip(upper=1)

                # Merge forward and reverse results
                rbh = pd.merge(fwd_results, rev_results, left_on=['subject', 'query'], right_on=['query', 'subject'], how='inner')

                rbh.to_csv(rbh_out, sep = '\t')

                del rbh

                fwd_results = pd.read_csv(fwd_best_in, sep="\t", header=None)
                rev_results = pd.read_csv(rev_best_in, sep="\t", header=None)

                headers = ["query", "subject", "identity", "coverage", "qlength", "slength", "alength", "bitscore", "E-value"]
                fwd_results.columns = headers
                rev_results.columns = headers

                # Create a new column in both dataframes: normalised bitscore
                fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
                rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

                # Create query and subject coverage columns in both dataframes
                fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
                rev_results['qcov'] = rev_results.alength/rev_results.qlength
                fwd_results['scov'] = fwd_results.alength/fwd_results.slength
                rev_results['scov'] = rev_results.alength/rev_results.slength

                # Clip maximum coverage values at 1.0
                fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
                rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
                fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
                rev_results['scov'] = rev_results['scov'].clip(upper=1)

                # Merge forward and reverse results
                rbbh = pd.merge(fwd_results, rev_results, left_on=['subject', 'query'], right_on=['query', 'subject'], how='inner')

                rbbh.to_csv(rbbh_out, sep = '\t')

                duplicates = rbbh[rbbh.duplicated(subset=['query_x', 'subject_x'], keep=False)]
                if not duplicates.empty:
                    print(f"\nWarning: Duplicate rows found based on ['query_x', 'subject_x']: for {a1.name} and {a2.name} RBBHs")
                    print(duplicates)

                del rbbh

    if not skip_orthofinder:
        if not pairwise_orthofinder:
            orthofinder_results_parent = protein_path / "orthofinder"
            existing_results = list(orthofinder_results_parent.glob("Results*")) if orthofinder_results_parent.exists() else []

            if existing_results:
                print(f"\n\tExisting global OrthoFinder results folder found in '{orthofinder_results_parent}'. Skipping OrthoFinder execution.")
            else:
                print(f"\nRunning OrthoFinder (this can take a very long time) between all annotations {annotation_names}")
                orthofinder_cmd = [
                    "orthofinder",
                    "-f", str(protein_path),
                    "-t", str(threads),
                    "-a", str(threads),
                    "-o", f"{str(protein_path)}/orthofinder/"
                ]
                run_command(results_directory, orthofinder_cmd)

    simple_annotations = []

    for n, a in enumerate(annotations):
        simple_annotations.append(Simple_annotation(a.name, a, group_names[n]))

    del annotations

    extra_tag = ""

    for n1, a1 in enumerate(simple_annotations):

        for n2, a2 in enumerate(simple_annotations):

            if n1 == n2:
                continue
            
            if reference_annotation != "None":
                if not a1.target and not a2.target:
                    continue

            if genome_files[n1] == genome_files[n2]:
                a1.add_liftless_overlap_equivalences(str(liftless_overlaps_dir), a1.name, a2.name, group_names[n2], quiet=quiet)
            else:
                if not skip_mcscan:
                    a1_mcscan_name = a1.name.replace(".", "_")
                    a2_mcscan_name = a2.name.replace(".", "_")

                    a1.add_mcscan_equivalences(f"{mcscan_path}/{a1_mcscan_name}.{a2_mcscan_name}.anchors", "0", a2.name, group_names[n2])
                    a1.add_mcscan_equivalences(f"{mcscan_path}/{a1_mcscan_name}.{a2_mcscan_name}.last.filtered", "0", a2.name, group_names[n2])

                if not skip_orthofinder:
                    if pairwise_orthofinder:
                        orthofile_pattern = f"orthofinder_{a1.name}__to__{a2.name}/orthofinder/Results*/Orthologues/Orthologues_{a1.name}_proteins_g_id_main/{a1.name}_proteins_g_id_main__v__{a2.name}_proteins_g_id_main.tsv"
                    else:
                        orthofile_pattern = f"orthofinder/Results*/Orthologues/Orthologues_{a1.name}_proteins_g_id_main/{a1.name}_proteins_g_id_main__v__{a2.name}_proteins_g_id_main.tsv"
                    matching_files = list(protein_path.glob(orthofile_pattern))

                    if not matching_files:
                        warnings.warn(f"No orthofinder file for {a1.name} vs {a2.name} found! Orthofinder results not added.", category=UserWarning)
                    elif len(matching_files) > 1:
                        warnings.warn(f"More than one orthofinder file for {a1.name} vs {a2.name} found! Orthofinder results not added.", category=UserWarning)
                    else:
                        ortho_file_path = matching_files[0]
                        a1.add_orthofinder_equivalences(str(ortho_file_path), a2.name, group_names[n2])

                if not skip_liftoff:
                    a1.add_reciprocal_overlap_equivalences(liftoff_path, a1.name, a2.name, group_names[n2], quiet=quiet)
                if not skip_lifton:
                    a1.add_reciprocal_overlap_equivalences(lifton_path, a1.name, a2.name, group_names[n2], liftoff=False, quiet=quiet)

            if not skip_all_blasts:
                a1.add_blast_equivalences(str(diamond_path), a1.name, a2.name, group_names[n2], skip_rbhs=skip_rbhs, skip_unidirectional_blasts=skip_unidirectional_blasts, quiet=quiet)

        output_file = f"{output_dir}{a1.name}_equivalences{extra_tag}.tsv"
        output_file_filtered_just_rbbhs_and_rbhs = f"{output_dir}{a1.name}_equivalences_just_rbbhs_and_rbhs{extra_tag}.tsv"
        output_file_filtered_just_rbbhs = f"{output_dir}{a1.name}_equivalences_just_rbbhs{extra_tag}.tsv"

        if skip_rbhs and skip_unidirectional_blasts and not skip_all_blasts:
            df = a1.export_summary_equivalences(output_file_filtered_just_rbbhs, filtered=True, simple_rbh_blasts=False, unidirectional_blasts=False, verbose=False, quiet=quiet, return_df=True, export_csv=False)

        elif skip_unidirectional_blasts and not skip_all_blasts:
            df = a1.export_summary_equivalences(output_file_filtered_just_rbbhs_and_rbhs, filtered=True, unidirectional_blasts=False, coverage_threshold=coverage, identity_threshold=identity, verbose=False, quiet=quiet, return_df=True, export_csv=False)

        else:
            df = a1.export_summary_equivalences(output_file, filtered=True, coverage_threshold=coverage, identity_threshold=identity, verbose=False, quiet=quiet, return_df=True, export_csv=False)

        if n1 == 0:
            final_df = df.copy()
        else:
            final_df = pd.concat([final_df, df], ignore_index=True)

    post_processing_start = time()

    if not final_df.empty:
        mask_standardize = (final_df['annotation_A'] > final_df['annotation_B']) & (final_df['annotation_B'] != "NA")
        final_df = flip_masked_rows(final_df, mask_standardize)

        groupby_cols = ['annotation_A', 'annotation_B', 'gene_id_A', 'gene_id_B']
        agg_dict = {col: 'first' for col in final_df.columns if col not in groupby_cols}
        
        if 'score' in agg_dict:
            agg_dict['score'] = merge_score_strings
        if 'summary_score' in agg_dict:
            agg_dict['summary_score'] = best_summary_score
        
        final_df = final_df.groupby(groupby_cols, as_index=False, dropna=False).agg(agg_dict)

        if include_duplicates:
            df_cloned = final_df.copy()
            mask_all = pd.Series([True] * len(df_cloned), index=df_cloned.index)
            df_cloned = flip_masked_rows(df_cloned, mask_all)
            
            final_df = pd.concat([final_df, df_cloned], ignore_index=True)
        else:
            if reference_annotation != "None":
                mask_ref = (final_df['annotation_B'] == reference_annotation) & (final_df['annotation_A'] != "NA")
                final_df = flip_masked_rows(final_df, mask_ref)

    clean_confidences = [c.strip().lower().replace('_confidence', '') for c in confidence]
    conf_set = set(clean_confidences)
    conf_order = {"high": 1, "medium": 2, "lower": 3}
    sorted_confs = sorted(list(conf_set), key=lambda x: conf_order.get(x, 4))

    standard_sets = [
        {"high"},
        {"high", "medium"},
        {"high", "medium", "lower"}
    ]

    valid_confidences = [f"{c}_confidence" for c in sorted_confs]

    if not final_df.empty and not skip_cardinality:
        if tiered_cardinality:
            if conf_set not in standard_sets:
                final_df['cardinality'] = get_tiered_cardinality(final_df, valid_confidences)

            final_df['cardinality_strict'] = get_tiered_cardinality(final_df, ['high_confidence'])
            final_df['cardinality_moderate'] = get_tiered_cardinality(final_df, ['high_confidence', 'medium_confidence'])
            final_df['cardinality_relaxed'] = get_tiered_cardinality(final_df, ['high_confidence', 'medium_confidence', 'lower_confidence'])
            
        else:
            final_df['cardinality'] = get_tiered_cardinality(final_df, valid_confidences)
        
    if conf_set != {"high", "medium", "lower"}:
        final_df = final_df[final_df["summary_score"].isin(valid_confidences)].copy()
        extra_tag += f"_confidence{'_'.join(sorted_confs)}"

    if include_NAs:
        extra_tag += "_with_NAs"
        matched_set = set()
        if not final_df.empty:
            matched_set.update(zip(final_df['annotation_A'], final_df['gene_id_A']))
            matched_set.update(zip(final_df['annotation_B'], final_df['gene_id_B']))
        
        na_rows = []

        for a in simple_annotations:
            for gene in a.genes.keys():
                if (a.name, gene) not in matched_set:
                    row_dict = {
                        "gene_id_A": gene,
                        "gene_id_B": "NA",
                        "score": "NA",
                        "summary_score": "NA",
                        "annotation_A": a.name,
                        "annotation_B": "NA",
                        "species_A": a.species,
                        "species_B": "NA"
                    }
                    if "cardinality" in final_df.columns:
                        row_dict["cardinality"] = "NA"
                    if "cardinality_strict" in final_df.columns:
                        row_dict["cardinality_strict"] = "NA"
                        row_dict["cardinality_moderate"] = "NA"
                        row_dict["cardinality_relaxed"] = "NA"
                    if "reliability" in final_df.columns:
                        row_dict["reliability"] = "NA"
                    na_rows.append(row_dict)
        
        if na_rows:
            na_df = pd.DataFrame(na_rows)
            final_df = pd.concat([final_df, na_df], ignore_index=True)

    if output_filename != "equivalences{other_tags}.tsv":
        final_output_file = f"{output_dir}{output_filename}"
        if not final_output_file.endswith(".tsv"):
            final_output_file += ".tsv"
    elif skip_rbhs and skip_unidirectional_blasts and not skip_all_blasts:
        final_output_file = f"{output_dir}equivalences_just_rbbhs{extra_tag}.tsv"
    elif skip_unidirectional_blasts and not skip_all_blasts:
        final_output_file = f"{output_dir}equivalences_just_rbbhs_and_rbhs{extra_tag}.tsv"
    else:
        final_output_file = f"{output_dir}equivalences{extra_tag}.tsv"

    if not final_df.empty:

        confidence_order = ["high_confidence", "medium_confidence", "lower_confidence", "NA"]
        final_df['summary_score'] = pd.Categorical(final_df['summary_score'], categories=confidence_order, ordered=True)
        
        species_order = sorted([s for s in final_df['annotation_B'].unique() if s != "NA"]) + ["NA"]
        final_df['annotation_B'] = pd.Categorical(final_df['annotation_B'], categories=species_order, ordered=True)

        sort_cols = ['annotation_A', 'annotation_B', 'gene_id_A', 'summary_score', 'gene_id_B']
        final_df = final_df.sort_values(by=sort_cols)
        
        final_df['summary_score'] = final_df['summary_score'].astype(str)
        final_df['annotation_B'] = final_df['annotation_B'].astype(str)

    if not final_df.empty and "species_A" in final_df.columns and "species_B" in final_df.columns:
        if (final_df["species_A"] == "NA").all() and (final_df["species_B"] == "NA").all():
            final_df = final_df.drop(columns=["species_A", "species_B"])

    if not final_df.empty and split_scores:

        if not quiet:
            print("\nSplitting aggregate scores into individual tool columns.")

        score_data = [parse_score_column(s) for s in final_df["score"]]
        split_df = pd.DataFrame(score_data, index=final_df.index)

        split_df = split_df.loc[:, (split_df != "NA").any(axis=0)]
        
        try:
            col_idx = list(final_df.columns).index("score")
        except ValueError:
            col_idx = 0

        final_df = final_df.drop(columns=["score"], errors="ignore")
        final_df = pd.concat([final_df, split_df], axis=1)

    cols = list(final_df.columns)
    card_cols = [c for c in cols if 'cardinality' in c]

    for col in card_cols:
        cols.remove(col)
    idx = cols.index('summary_score') + 1

    for i, col in enumerate(card_cols):
        cols.insert(idx + i, col)

    final_df = final_df[cols]

    print(f"\nPost processing took: {round((time() - post_processing_start) / 60, 2)} minutes.")

    final_df.to_csv(final_output_file, sep="\t", encoding="utf-8", index=False)

    if not keep_intermediate:
        if os.path.exists(str(results_directory)):
            shutil.rmtree(str(results_directory))

    end = time()

    print(f"\naegis orthology command ('{' '.join(sys.argv)}') complete. Total time: {round((end - start) / 60, 2)} minutes.")

if __name__ == "__main__":
    try:
        app()
    except Exception as e:
        import traceback
        typer.echo("aegis orthology crashed with an unexpected error:", err=True)
        typer.echo(str(e), err=True)
        typer.echo(traceback.format_exc(), err=True)
        raise
