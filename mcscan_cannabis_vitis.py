
import os

root = f"../../genomes_and_annotation/"

root = f"{os.path.abspath(root)}/"

threads = 8

input_folder = f"{root}cannabis/aegis_output/features/"

input_tags = ["CBDRx", "Finola", "PurpleKush", "JamaicanLion"]

input_fastas = ["CBDRx_on_CBDRx", "Finola_on_Finola", "PurpleKush_on_PurpleKush", "JamaicanLion_on_JamaicanLion"]

target_folder = f"{root}grapevine/PN40024/aegis_output/features/"

target_tags = ["5_1"]

target_tag_final_rename = ["5.1"]

target_fastas = ["5.1_on_T2T_ref"]

gff_input_folder = f"{root}cannabis/aegis_output/gffs/"

gff_target_folder = f"{root}grapevine/PN40024/aegis_output/gffs/"

input_gffs = ["CBDRx_on_CBDRx_aegis_fcounts.gff3", "Finola_on_Finola_aegis_fcounts.gff3", "PurpleKush_on_PurpleKush_aegis_fcounts.gff3", "JamaicanLion_on_JamaicanLion_aegis_fcounts.gff3"]

target_gffs = ["5.1_on_T2T_ref_aegis_fcounts_dapfit.gff3"]

output_folder = f"{root}aegis_output/annotation_transfers/"

current_path = os.getcwd()

os.system(f"mkdir -p {output_folder}")

for n, input_fasta in enumerate(input_fastas):
    for n2, target_fasta in enumerate(target_fastas):
        
        mcscan_folder = f"{output_folder}mcscan_{input_tags[n]}_{target_tag_final_rename[n2]}/"

        anchors = f"{mcscan_folder}{input_tags[n]}.{target_tags[n2]}.anchors"
        last_filtered = f"{mcscan_folder}{input_tags[n]}.{target_tags[n2]}.last.filtered"

        final_anchors = f"{mcscan_folder}{input_tags[n]}.{target_tag_final_rename[n2]}.anchors"
        final_last_filtered = f"{mcscan_folder}{input_tags[n]}.{target_tag_final_rename[n2]}.last.filtered"

        os.system(f"mkdir -p {mcscan_folder}")

        os.system(f"rm -f {final_last_filtered}")
        os.system(f"rm -f {final_anchors}")
        os.system(f"rm -f {last_filtered}")
        os.system(f"rm -f {anchors}")

        os.system(f"python3 -m jcvi.formats.fasta format {input_folder}{input_fasta}_CDSs_g_id_main.fasta {mcscan_folder}{input_tags[n]}.cds")

        os.system(f"python3 -m jcvi.formats.fasta format {target_folder}{target_fasta}_CDSs_g_id_main.fasta {mcscan_folder}{target_tags[n2]}.cds")
        
        os.system(f"python3 -m jcvi.formats.gff bed --type=mRNA --key=featurecounts_id --primary_only {gff_input_folder}{input_gffs[n]} -o {mcscan_folder}{input_tags[n]}.bed")

        os.system(f"python3 -m jcvi.formats.gff bed --type=mRNA --key=featurecounts_id --primary_only {gff_target_folder}{target_gffs[n2]} -o {mcscan_folder}{target_tags[n2]}.bed")

        os.chdir(mcscan_folder)

        os.system(f"python3 -m jcvi.compara.catalog ortholog {input_tags[n]} {target_tags[n2]} --no_strip_names")

        os.chdir(current_path)

        if target_tag_final_rename[n2] != target_tags[n2]:

            if not os.path.isfile(final_anchors):
                os.system(f"mv {anchors} {final_anchors}")
            if not os.path.isfile(final_last_filtered):
                os.system(f"mv {last_filtered} {final_last_filtered}")