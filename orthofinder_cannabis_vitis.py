
import os

root = f"../../genomes_and_annotation/"

orthofinder = "../../apps/OrthoFinder_source/orthofinder.py"

root = f"{os.path.abspath(root)}/"

threads = 8

orthofinder = os.path.abspath(orthofinder)

protein_file_tag = "_proteins_g_id_main.fasta"

input_folder = f"{root}cannabis/aegis_output/features/"

input_tags = ["CBDRx", "Finola", "PurpleKush", "JamaicanLion"]

input_fastas = ["CBDRx_on_CBDRx", "Finola_on_Finola", "PurpleKush_on_PurpleKush", "JamaicanLion_on_JamaicanLion"]

target_folder = f"{root}grapevine/PN40024/aegis_output/features/"

target_tags = ["5.1"]

target_fastas = ["5.1_on_T2T_ref"]

output_folder = f"{root}aegis_output/annotation_transfers/"

os.system(f"mkdir -p {output_folder}")

for n, input_fasta in enumerate(input_fastas):
    for n2, target_fasta in enumerate(target_fastas):
        orthofinder_folder = f"{output_folder}orthofinder_{input_tags[n]}_{target_tags[n2]}/"

        os.system(f"mkdir -p {orthofinder_folder}")

        os.system(f"cp {input_folder}{input_fasta}{protein_file_tag} {orthofinder_folder}{input_fasta}{protein_file_tag}")
        os.system(f"cp {target_folder}{target_fasta}{protein_file_tag} {orthofinder_folder}{target_fasta}{protein_file_tag}")

        os.system(f"python3 {orthofinder} -t {threads} -a {threads} -f {orthofinder_folder}")

        os.system(f"rm {orthofinder_folder}{input_fasta}{protein_file_tag}")
        os.system(f"rm {orthofinder_folder}{target_fasta}{protein_file_tag}")
