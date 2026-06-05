import typer
import os
import pandas as pd
import numpy as np
import warnings
import shutil

from pathlib import Path
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from ..feature import Feature
from ..equivalence import Simple_annotation, pairwise_orthology, run_command

app = typer.Typer(add_completion=False, no_args_is_help=True)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_files: Annotated[list[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    genome_files: Annotated[str, typer.Option(
        "-g", "--genome-fastas", help="Genome assemblies corresponding to annotation files. Provide them in the same number and order, separated by commas. e.g. -g genomefile1,genomefile2,genomefile3,genomefile4",
        callback=split_callback
    )] = "",
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", help="Annotation versions, names or tags. Provide them in the same number and order as the corresponding annotation files, separated by commas. e.g. --annotation-names name1,name2,name3,name4",
        callback=split_callback
    )] = "{annotation-filename(s)}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    output_filename: Annotated[str, typer.Option(
        "-o", "--output-file", help="Output filename to be saved to output folder without extension. The '.tsv' extension will be added to the filename."
    )] = "equivalences{other_tags}.tsv",
    group_names: Annotated[str, typer.Option(
        "-gn", "--group-names", help="Optional grouping of input annotations, into species for example. Use NA as a placemarker for annotation files without a group label. e.g. --group-names group1,NA,group1,group2",
        callback=split_callback
    )] = "",
    skip_synteny: Annotated[bool, typer.Option(
        "--skip-synteny", help="Skip conservation of synteny metrics whenever an annotation is lifted over to another genome."
    )] = False,
    reference_annotation: Annotated[str, typer.Option(
        "--reference-annotation", help="Select a single annotation, by providing its name/tag or filename, to use as a reference. Only matches to and from this annotation will be reported. Otherwise matches are reported between all annotations."
    )] = "None",
    include_single_blasts: Annotated[bool, typer.Option(
        "-b", "--include-single-blasts", help="Decide whether to report unidirectional (i.e. just fw or rv) blasts in the orthologue summary."
    )] = False,
    threads: Annotated[int, typer.Option(
        "-t", "--threads", help="Number of threads."
    )] = 1,
    skip_rbhs: Annotated[bool, typer.Option(
        "--skip-RBHs", help="Decide whether to skip RBHs which are not RBBHs, these are reported by default in the orthologue summary."
    )] = False,
    lift_feature_types: Annotated[str, typer.Option(
        "--lift-feature-types", help="All feature types within an annotation files are lifted over by default, however a more restrictive set can be used, separated by commas, such as 'gene,mRNA,exon,CDS,pseudogene,pseudogenic_exon,pseudogenic_transcript'.", callback=split_callback
    )] = "ALL",
    skip_lifton: Annotated[bool, typer.Option(
        "--skip-lifton", help="Skip LiftOn."
    )] = False,
    skip_liftoff: Annotated[bool, typer.Option(
        "--skip-liftoff", help="Skip Liftoff."
    )] = False,
    skip_copies: Annotated[bool, typer.Option(
        "--skip-copies", help="Liftoff and LiftOn are run in copies mode my default, flag to deactivate."
    )] = False,
    skip_mcscan: Annotated[bool, typer.Option(
        "--skip-mcscan", help="Skip the JCVI toolkit synteny and collinearity analysis (MCScan). Useful when JCVI is causing compatibility issues."
    )] = False,
    skip_orthofinder: Annotated[bool, typer.Option(
        "--skip-orthofinder", help="Skip the OrthoFinder analysis."
    )] = False,
    pairwise_orthofinder: Annotated[bool, typer.Option(
        "--pairwise-orthofinder", help="Execute OrthoFinder on independent annotation pairs. Overrides the default multi-annotation bulk analysis. Recommended for highly divergent taxa or targeted 1:1 orthologue mapping."
    )] = False,
    skip_all_blasts: Annotated[bool, typer.Option(
        "--skip-all-blasts", help="Skip all BLASTs."
    )] = False,
    keep_intermediate: Annotated[bool, typer.Option(
        "-k", "--keep-intermediate", help="Keep intermediate files, useful for identifying errors."
    )] = False,
    include_duplicates: Annotated[bool, typer.Option(
        "--include-duplicates", help="Report equivalences from both from gene_id_A to gene_id_B as well as from gene_id_B to gene_id_A. These 'duplicate gene pairs' are not included by default."
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Verbose logging, useful if encountering a problem or error."
    )] = False,
    identity: Annotated[float, typer.Option(
        "-i", "--identity", help="Minimum identity threshold for BLAST hits."
    )] = 30.0,
    coverage: Annotated[float, typer.Option(
        "-c", "--coverage", help="Minimum coverage threshold for BLAST hits."
    )] = 30.0,
    evalue: Annotated[float, typer.Option(
        "-e", "--evalue", help="Maximum e-value threshold for BLAST hits."
    )] = 0.00001


):
    """
    Provides a set of orthologues relationships leveraging external and internal tools such as Litfoff + AEGIS overlaps, LiftOn + AEGIS overlaps, MCScan, orthofinder and BLAST. Wherever relevant, tools are run reciprocally for an extra confidence mark in orthologous relationships.

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

    if len(annotation_files) != len(annotation_names):
        raise typer.BadParameter(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    for annotation_name in annotation_names:
        if "__to__" in annotation_name:
            raise typer.BadParameter(f"The provided annotation name/tag '{annotation_name}' has an incompatible term: '__to__' as it is used internally for temporary file naming.")

    if skip_mcscan and skip_liftoff and skip_all_blasts and skip_orthofinder and skip_lifton:
        raise typer.BadParameter("No analysis methods selected. Please select at least one analysis method. Run 'aegis orthology --help' for more information.")

    if len(annotation_names) != len(set(annotation_names)):
        raise typer.BadParameter("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annotation_files) != len(set(annotation_files)):
        raise typer.BadParameter("Avoid repeated annotation filename(s).")

    if len(genome_files) != len(annotation_files):
        raise typer.BadParameter(f"A single genome file must be provided for each annotation file. Genome files: {genome_files}. Annotation files: {annotation_files}")
    
    if len(genome_files) != len(set(genome_files)):
        repeated_genomes = {}
        for idx, genome in enumerate(genome_files):
            repeated_genomes.setdefault(genome, []).append(annotation_names[idx])
            
        repeated_msg = []
        for genome, a_names in repeated_genomes.items():
            if len(a_names) > 1:
                repeated_msg.append(f"  - Genome '{genome}' used by annotations: {', '.join(a_names)}")
                
        print("\n⚠️  Warning: Repeated genome assemblies.\n"
              "Bear in mind that annotations associated to the same genome will only be compared at the level of 'aegis overlap' and BLAST results. Also, OrthoFinder will only be run in pairwise mode, except for the annotation pairs which share the same genome.\n"
              "Repeated assignments found:\n" + "\n".join(repeated_msg) + "\n")
        
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

    genomes:dict[str, Genome] = { g: Genome(name=g, genome_file_path=g, quiet=quiet) for g in set(genome_files) }

    annotations:list[Annotation] = []

    for n, annotation_file in enumerate(annotation_files):

        annotations.append(Annotation(name=annotation_names[n], genome=genomes[genome_files[n]], annot_file_path=annotation_file, quiet=quiet))

        annotations[-1].rename_ids(strip_gene_tag=True, quiet=quiet)

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True

    output_dir_path = Path(output_dir).resolve() / "orthologues"
    output_dir = str(output_dir_path) + "/"

    if output_dir_path.exists():
        warnings.warn(f"The folder '{output_dir}' already exists. Please be aware that conflict may arise with existing output and/or temp folder files.")

    output_dir_path.mkdir(parents=True, exist_ok=True)

    results_directory = Path(f"{output_dir}temp/")
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
        a.export.gff(output_dir=str(gff_path), filename=f"{a.name}.gff3", subfolder=False, quiet=quiet)
        a.export.gff(output_dir=str(gff_path), filename=f"{mcscan_name}.gff3", subfolder=False, quiet=quiet)

        if not skip_lifton:

            a_lifton = a.copy()
            a_lifton.CDS_to_CDS_segment_ids(override=True)
            a_lifton.export.gff(output_dir=str(gff_path), filename=f"{a_lifton.name}__for__lifton.gff3", subfolder=False, quiet=quiet)

            del a_lifton

        Feature._ACTIVE_GENOME = a.genome

        a.export.proteins(only_main=True, custom_path=str(protein_path), used_id="gene", verbose=False, custom_filename=f"{a.name}_proteins_g_id_main.fasta")
        a.clear_proteins()
        
        a.export.CDSs(only_main=True, custom_path=str(CDS_path), used_id="gene", verbose=False, custom_filename=f"{mcscan_name}_CDSs_g_id_main.fasta")
        protein_fasta = protein_path / f"{a.name}_proteins_g_id_main.fasta"


        if not skip_all_blasts:
            diamond_db_file = diamond_path / f"{a.name}_diamond_db"
            makedb_cmd = [
                "diamond", "makedb", "-p", str(threads), "--in", str(protein_fasta), "--db", str(diamond_db_file)
            ]
            run_command(diamond_path, makedb_cmd)

        if not skip_mcscan:
            cds_fasta = CDS_path / f"{mcscan_name}_CDSs_g_id_main.fasta"            

            cleaned_cds = mcscan_path / Path(f"{mcscan_name}.cds")

            jcvi_format_cmd_1 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta), str(cleaned_cds)]
            run_command(mcscan_path, jcvi_format_cmd_1)

            bed_file = mcscan_path / Path(f"{mcscan_name}.bed")
            gff_to_bed_cmd_1 = [
                "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
                "--key=Parent", "--primary_only", f"{gff_path}/{mcscan_name}.gff3", "-o", str(bed_file)
            ]
            run_command(mcscan_path, gff_to_bed_cmd_1)


    for n1, a1 in enumerate(annotations):

        for n2, a2 in enumerate(annotations):

            if n1 == n2:
                continue

            pairwise_orthology(annot1=a1, annot2=a2, genome1=genomes[genome_files[n1]], genome2=genomes[genome_files[n2]], working_directory=results_directory, num_threads=threads, copies=not(skip_copies), synteny=synteny, skip_liftoff=skip_liftoff, skip_lifton=skip_lifton, skip_mcscan=skip_mcscan, types=lift_feature_types_file, coverage=coverage, evalue=evalue, skip_blasts=skip_all_blasts, pairwise_orthofinder=pairwise_orthofinder, skip_orthofinder=skip_orthofinder,quiet=quiet)

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

                print(f"\nProcessing RBH and RBBHs for {a1.name} and {a2.name}")

                fwd_in = diamond_path / f"single_{a1.name}__to__{a2.name}.txt"
                rev_in = diamond_path / f"single_{a2.name}__to__{a1.name}.txt"
                fwd_best_in = diamond_path / f"single_best_{a1.name}__to__{a2.name}.txt"
                rev_best_in = diamond_path / f"single_best_{a2.name}__to__{a1.name}.txt"
                rbh_out = diamond_path / f"rbh_{a1.name}__to__{a2.name}.txt"
                rbbh_out = diamond_path / f"rbbh_{a1.name}__to__{a2.name}.txt"

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

                    a1.add_mcscan_equivalences(f"{mcscan_path}/{a1_mcscan_name}.{a2_mcscan_name}.anchors", "0", a2_mcscan_name, group_names[n2])
                    a1.add_mcscan_equivalences(f"{mcscan_path}/{a1_mcscan_name}.{a2_mcscan_name}.last.filtered", "0", a2_mcscan_name, group_names[n2])

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

    if output_filename != "equivalences{other_tags}.tsv":
        final_output_file = f"{output_dir}{output_filename}.tsv"
    elif skip_rbhs and skip_unidirectional_blasts and not skip_all_blasts:
        final_output_file = f"{output_dir}equivalences_just_rbbhs{extra_tag}.tsv"
    elif skip_unidirectional_blasts and not skip_all_blasts:
        final_output_file = f"{output_dir}equivalences_just_rbbhs_and_rbhs{extra_tag}.tsv"
    else:
        final_output_file = f"{output_dir}equivalences{extra_tag}.tsv"

    if not include_duplicates:

        final_df['gene_id_tuple'] = pd.DataFrame(np.sort(final_df[['gene_id_A', 'gene_id_B']], axis=1), index=final_df.index).agg(tuple, axis=1)
        final_df['annotation_tuple'] = pd.DataFrame(np.sort(final_df[['annotation_A', 'annotation_B']], axis=1), index=final_df.index).agg(tuple, axis=1)
        final_df['species_tuple'] = pd.DataFrame(np.sort(final_df[['species_A', 'species_B']], axis=1), index=final_df.index).agg(tuple, axis=1)

        subset_for_duplicates = ["gene_id_tuple", "annotation_tuple", "species_tuple"]

        final_df.drop_duplicates(subset=subset_for_duplicates, keep='first', inplace=True)

        final_df.drop(subset_for_duplicates, axis=1, inplace=True)

    final_df.to_csv(final_output_file, sep="\t", encoding="utf-8", index=False)

    if not keep_intermediate:
        if os.path.exists(str(results_directory)):
            shutil.rmtree(str(results_directory))

    print(f"aegis orthology run complete.")

if __name__ == "__main__":
    try:
        app()
    except Exception as e:
        import traceback
        typer.echo("aegis orthology crashed with an unexpected error:", err=True)
        typer.echo(str(e), err=True)
        typer.echo(traceback.format_exc(), err=True)
        raise
