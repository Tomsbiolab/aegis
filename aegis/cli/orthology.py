import typer
import os
import pandas as pd
import warnings
import shutil
from pathlib import Path
from typing import List
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.equivalence import Simple_annotation, pairwise_orthology, run_command

app = typer.Typer(add_completion=False)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    genome_files: Annotated[str, typer.Option(
        "-g", "--genome-fastas", help="Genome assemblies corresponding to annotation files. Provide them in the same number and order, separated by commas. e.g. genomefile1,genomefile2,genomefile3,genomefile4",
        callback=split_callback
    )],
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", help="Annotation versions, names or tags. Provide them in the same number and order as the corresponding annotation files, separated by commas. e.g. name1,name2,name3,name4",
        callback=split_callback
    )] = "{annotation-filename(s)}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    group_names: Annotated[str, typer.Option(
        "-gn", "--group-names", help="Optional grouping of input annotations, into species for example. Use NA as a placemarker for annotation files without a group label. e.g. '-g group1,NA,group1,group2'",
        callback=split_callback
    )] = "",
    original_annotation_files: Annotated[str, typer.Option(
        "-ot", "--original-annotation-files", help="Should some of the annotations be a result of a liftover or coordinate transfer, you can optionally provide a list of the original files before the transfer, separated by commas. If at least 2 annotation files are being compared, conservation of synteny will be calculated wherever possible based on gene order before/after transfer. These original annotation files must be in the same number and order as the corresponding annotation files. Use NA as a placemarker for annotation files without an original annotation file. e.g. '-t original_file_1,NA,original_file_3'",
        callback=split_callback
    )] = "",
    reference_annotation: Annotated[str, typer.Option(
        "-r", "--reference-annotation", help="Select a single annotation, by providing its name/tag or filename, to use as a reference. Only matches to and from this annotation will be reported. Otherwise matches are reported between all annotations."
    )] = "None",
    include_single_blasts: Annotated[bool, typer.Option(
        "-b", "--include_single_blasts", help="Decide whether to report unidirectional (i.e. just fw or rv) blasts in the orthologue summary."
    )] = False,
    threads: Annotated[int, typer.Option(
        "-t", "--threads", help="Number of threads."
    )] = 5,
    skip_rbhs: Annotated[bool, typer.Option(
        "-rb", "--skip_RBHs", help="Decide whether to skip RBHs which are not RBBHs, these are reported by default in the orthologue summary."
    )] = False,
    lift_feature_types: Annotated[str, typer.Option(
        "-lt", "--lift_feature_types", help="All feature types within an annotation files are lifted over by default, however a more restrictive set can be used, separated by commas, such as 'gene,mRNA,exon,CDS,pseudogene,pseudogenic_exon,pseudogenic_transcript'.", callback=split_callback
    )] = "ALL",
    skip_lifton: Annotated[bool, typer.Option(
        "-sl", "--skip_lifton", help="Skip LiftOn, use flag in case LiftOn is causing compatibility issues."
    )] = False,
    skip_copies: Annotated[bool, typer.Option(
        "-cl", "--skip_copies", help="Liftoff and Lifton are run in copies mode my default, flag to deactivate."
    )] = False,
    keep_intermediate: Annotated[bool, typer.Option(
        "-k", "--keep_intermediate", help="Keep intermediate files, useful for identifying errors."
    )] = False
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

    if len(annotation_files) < 2:
        raise typer.BadParameter(f"At least 2 annotation-files must be provided.")
    
    if len(annotation_files) > 1 and annotation_files[-1].lower() in ("true", "false"):
        typer.echo(
            "⚠️  Detected extra value 'true' or 'false' at the end of positional arguments.\n"
            "👉 Did you mean to use the '--include_NAs' or '--simple' flags? Use them like this: '-n' or '-s' (no 'true' needed).",
            err=True,
        )
        raise typer.Exit(code=1)

    if annotation_names == "{annotation-filename(s)}":
        annotation_names = []
        for annotation_file in annotation_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0])

    if len(annotation_files) != len(annotation_names):
        raise ValueError(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if len(annotation_names) != len(set(annotation_names)):
        raise ValueError("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annotation_files) != len(set(annotation_files)):
        raise ValueError("Avoid repeated annotation filename(s).")
    
    if len(genome_files) != len(set(genome_files)):
        raise ValueError("Avoid repeated genome assemblies. If looking to compare annotation versions associated to the same genome assembly, 'aegis-overlap' may be more appropriate.")

    if original_annotation_files:
        synteny = True
        if len(annotation_files) != len(original_annotation_files):
            raise ValueError(f"The provided number of original annotation files do not match the number of annotation file(s).")
        
    else:
        synteny = False
        original_annotation_files = ["NA"] * len(annotation_files)
    
    if group_names:
        if len(annotation_files) != len(group_names):
            raise ValueError(f"The provided number of groups do not match the number of annotation file(s).")
        
    else:
        group_names = ["NA"] * len(annotation_files)

    if reference_annotation != "None":
        if reference_annotation not in annotation_files and reference_annotation not in annotation_names:
            raise ValueError(f"The provided reference-annotation = {reference_annotation} is not present neither in annotation-files ({annotation_files}) nor annotation-names ({annotation_names}).")

    skip_unidirectional_blasts = not (include_single_blasts)

    if skip_rbhs and not skip_unidirectional_blasts:
        raise ValueError(f"Do not include single blasts if rbhs are to be skipped as these provide higher support for orthology.")

    genomes = [ Genome(name=f"{annotation_names[n]}_genome", genome_file_path=g) for n, g in enumerate(genome_files) ]

    annotations = []

    for n, annotation_file in enumerate(annotation_files):

        if original_annotation_files[n].lower() != "na":
            original_annotation = Annotation(name=f"{annotation_names[n]}_original", genome=genome_files[n], annot_file_path=original_annotation_files[n])
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, original_annotation=original_annotation))
        else:
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file))

        annotations[-1].rename_ids(strip_gene_tag=True, quiet=True)

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True

    output_folder = Path(output_folder).resolve() / "orthologues"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_folder = str(output_folder) + "/"

    results_directory = Path(f"{output_folder}temp/")
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

    mcscan_path = results_directory / "mcscan"
    mcscan_path.mkdir(parents=True, exist_ok=True)

    if lift_feature_types == ["ALL"]:
        lift_feature_types = ["gene", "mRNA", "exon", "CDS", "pseudogene", "pseudogenic_exon", "pseudogenic_transcript"]
    
    lift_feature_types_file = results_directory / "chosen_liftover_features.txt"
    lift_feature_types_file = str(lift_feature_types_file)

    f_in = open(lift_feature_types_file, "w", encoding="utf-8")

    for ft in lift_feature_types:
        ft = ft.strip()
        f_in.write(f"{ft}\n")
    f_in.close()

    identity = 30
    coverage = 30

    # Create gff, protein, CDS files, mcscan, and diamond databases in a non-redundant manner
    for n, a in enumerate(annotations):

        a.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False, quiet=True)
        a.export_gff(custom_path=str(gff_path), tag=f"{a.name}.gff3", subfolder=False, quiet=True)

        if not skip_lifton:

            a_lifton = a.copy()
            a_lifton.CDS_to_CDS_segment_ids(quiet=True)
            a_lifton.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False, quiet=True)
            a_lifton.export_gff(custom_path=str(gff_path), tag=f"{a_lifton.name}_for_lifton.gff3", subfolder=False, no_1bp_features=True, quiet=True)

            del a_lifton

        a.generate_sequences(genomes[n])
        a.export_proteins(only_main=True, custom_path=str(protein_path), used_id="gene", verbose=False)
        a.export_CDSs(only_main=True, custom_path=str(CDS_path), used_id="gene", verbose=False)

        protein_fasta = protein_path / f"{a.name}_proteins_g_id_main.fasta"

        diamond_db_file = diamond_path / f"{a.name}_diamond_db"
        makedb_cmd = [
            "diamond", "makedb", "-p", str(threads), "--in", str(protein_fasta), "--db", str(diamond_db_file)
        ]
        run_command(diamond_path, makedb_cmd)

        cds_fasta = CDS_path / f"{a.name}_CDSs_g_id_main.fasta"
        cleaned_cds = mcscan_path / f"{a.name}.cds"

        jcvi_format_cmd_1 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta), str(cleaned_cds)]
        run_command(mcscan_path, jcvi_format_cmd_1)

        bed_file = mcscan_path / f"{a.name}.bed"
        gff_to_bed_cmd_1 = [
            "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
            "--key=Parent", "--primary_only", f"{gff_path}/{a.name}.gff3", "-o", str(bed_file)
        ]
        run_command(mcscan_path, gff_to_bed_cmd_1)

    for n1, a1 in enumerate(annotations):

        for n2, a2 in enumerate(annotations):

            if n1 == n2:
                continue

            original_annotation = original_annotation_files[n1].lower()
            if original_annotation == "na":
                original_annotation = None
            else:
                original_annotation = Annotation(original_annotation_files[n1])
            
            pairwise_orthology(annot1=a1, annot2=a2, genome1=genomes[n1], genome2=genomes[n2], working_directory=results_directory, num_threads=threads, original_annot1=original_annotation, copies=not(skip_copies), synteny=synteny, skip_lifton=skip_lifton, types=lift_feature_types_file)

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

            fwd_in = diamond_path / f"single_{a1.name}_to_{a2.name}.txt"
            rev_in = diamond_path / f"single_{a2.name}_to_{a1.name}.txt"
            fwd_best_in = diamond_path / f"single_best_{a1.name}_to_{a2.name}.txt"
            rev_best_in = diamond_path / f"single_best_{a2.name}_to_{a1.name}.txt"
            rbh_out = diamond_path / f"rbh_{a1.name}_to_{a2.name}.txt"
            rbbh_out = diamond_path / f"rbbh_{a1.name}_to_{a2.name}.txt"

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

            a1.add_mcscan_equivalences(f"{mcscan_path}/{a1.name}.{a2.name}.anchors", "0", a2.name, group_names[n2])
            a1.add_mcscan_equivalences(f"{mcscan_path}/{a1.name}.{a2.name}.last.filtered", "0", a2.name, group_names[n2])

            orthofile_pattern = f"orthofinder/Results*/Orthologues/Orthologues_{a1.name}_proteins_g_id_main/{a1.name}_proteins_g_id_main__v__{a2.name}_proteins_g_id_main.tsv"
            matching_files = list(protein_path.glob(orthofile_pattern))

            if not matching_files:
                warnings.warn(f"No orthofinder file for {a1.name} vs {a2.name} found! Orthofinder results not added.", category=UserWarning)
            elif len(matching_files) > 1:
                warnings.warn(f"More than one orthofinder file for {a1.name} vs {a2.name} found! Orthofinder results not added.", category=UserWarning)
            else:
                ortho_file_path = matching_files[0]
                a1.add_orthofinder_equivalences(str(ortho_file_path), a2.name, group_names[n2])

            a1.add_reciprocal_overlap_equivalences(liftoff_path, a1.name, a2.name, group_names[n2], quiet=True)
            if not skip_lifton:
                a1.add_reciprocal_overlap_equivalences(lifton_path, a1.name, a2.name, group_names[n2], liftoff=False, quiet=True)

            a1.add_blast_equivalences(str(diamond_path), a1.name, a2.name, group_names[n2], skip_rbhs=skip_rbhs, skip_unidirectional_blasts=skip_unidirectional_blasts, quiet=True)

        output_file = f"{output_folder}{a1.name}_equivalences{extra_tag}.tsv"
        output_file_filtered_just_rbbhs_and_rbhs = f"{output_folder}{a1.name}_equivalences_filtered_just_rbbhs_and_rbhs{extra_tag}.tsv"
        output_file_filtered_just_rbbhs = f"{output_folder}{a1.name}_equivalences_filtered_just_rbbhs{extra_tag}.tsv"


        if skip_rbhs and skip_unidirectional_blasts:
            a1.export_summary_equivalences(output_file_filtered_just_rbbhs, filtered=True, simple_rbh_blasts=False, unidirectional_blasts=False, verbose=False, quiet=True)

        elif skip_unidirectional_blasts:
            a1.export_summary_equivalences(output_file_filtered_just_rbbhs_and_rbhs, filtered=True, unidirectional_blasts=False, coverage_threshold=coverage, identity_threshold=identity, verbose=False, quiet=True)

        else:
            a1.export_summary_equivalences(output_file, filtered=True, coverage_threshold=coverage, identity_threshold=identity, verbose=False, quiet=True)

    if not keep_intermediate:
        if os.path.exists(str(results_directory)):
            shutil.rmtree(str(results_directory))

if __name__ == "__main__":
    app()
