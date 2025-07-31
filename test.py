from aegis.genome import Genome
from aegis.annotation import Annotation
import os
import re
import shutil
import subprocess
from pathlib import Path

def clean_fasta_headers(input_path, output_path):
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                line = re.sub(r"_t\d{3}_CDS\d*", "", line)
            outfile.write(line)

def run_docker_jcvi_format(pwd, input_file, output_file):
    subprocess.run([
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", pwd,
        "jcvi-image2", "python", "-m", "jcvi.formats.fasta", "format",
        input_file, output_file
    ], check=True)

def run_docker_jcvi_gff_to_bed(pwd, gff_file, output_bed):
    subprocess.run([
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", pwd,
        "jcvi-image2", "python", "-m", "jcvi.formats.gff", "bed",
        "--type=mRNA", "--key=Parent", "--primary_only",
        gff_file, "-o", output_bed
    ], check=True)

def run_ortholog_analysis(pwd, species1, species2):
    subprocess.run([
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", pwd,
        "jcvi-image2", "python", "-m", "jcvi.compara.catalog",
        "ortholog", species1, species2, "--no_strip_names"
    ], check=True)

def prepare_gff(pwd, original_gff_path, target_species_name):
    renamed_path = Path(pwd) / f"{target_species_name}.gff3"
    shutil.copyfile(original_gff_path, renamed_path)
    return str(renamed_path)


source_genome = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/T2T_ref_dapfit_organellefree_chr00.fasta'
target_genome = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_organellefree.fasta'

source_gff = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/5.1_on_T2T_ref_aegis_fcounts_dapfit_for_lifton.gff3'
target_gff = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_on_CBDRx_aegis_fcounts_for_lifton.gff3'

pwd = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/tmp'
os.makedirs(pwd, exist_ok=True)
threads = '5'


os.system(f'docker run --rm -v /media:/media -w {pwd} quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0 liftoff {source_genome} {target_genome} -g {target_gff} -o forward_liftoff.gff')
os.system(f'docker run --rm -v /media:/media -w {pwd} quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0 liftoff {target_genome} {source_genome} -g {source_gff} -o reverse_liftoff.gff')

os.system(f'docker run --rm -v /media:/media -w {pwd} antsanpaj/lifton lifton -g {target_gff} -o forward_lifton.gff3 -copies {source_genome} {target_genome}')
os.system(f'docker run --rm -v /media:/media -w {pwd} antsanpaj/lifton lifton -g {source_gff} -o forward_lifton.gff3 -copies {target_genome} {source_genome}')

genome = Genome('source_genome', source_genome)
annot = Annotation('source_annot', source_gff, genome)
annot.generate_sequences(genome)
annot.export_proteins(only_main=True, custom_path=pwd, unique_proteins_per_gene=True, used_id="gene", verbose=False)
annot.export_CDSs(only_main=True, custom_path=pwd, used_id="gene", verbose=False)

genome = Genome('target_genome', target_genome)
annot = Annotation('target_annot', target_gff, genome)
annot.generate_sequences(genome)
annot.export_proteins(only_main=True, custom_path=pwd, unique_proteins_per_gene=True, used_id="gene", verbose=False)
annot.export_CDSs(only_main=True, custom_path=pwd, used_id="gene", verbose=False)

os.system(f'diamond makedb -p {threads} --in {pwd}/features/source_annot_proteins_g_id_main.fasta -d {pwd}/source_diamond_db')
os.system(f'diamond blastp --threads {threads} --db {pwd}/source_diamond_db --ultra-sensitive --out {pwd}/source.diamond --outfmt 6 qseqid sseqid pident qcovhsp qlen slen length bitscore evalue --query {pwd}/features/target_annot_proteins_g_id_main.fasta --max-target-seqs 1 --evalue 0.00001 --max-hsps 1')

os.system(f'diamond makedb -p {threads} --in {pwd}/features/target_annot_proteins_g_id_main.fasta -d {pwd}/target_diamond_db')
os.system(f'diamond blastp --threads {threads} --db {pwd}/target_diamond_db --ultra-sensitive --out {pwd}/target.diamond --outfmt 6 qseqid sseqid pident qcovhsp qlen slen length bitscore evalue --query {pwd}/features/source_annot_proteins_g_id_main.fasta --max-target-seqs 1 --evalue 0.00001 --max-hsps 1')

fasta1 = f'{pwd}/features/source_annot_CDSs_g_id_main.fasta'
fasta2 = f'{pwd}/features/target_annot_CDSs_g_id_main.fasta'

species1 = Path(fasta1).stem
species2 = Path(fasta2).stem

raw1 = Path(pwd) / f"{species1}_raw.cds"
raw2 = Path(pwd) / f"{species2}_raw.cds"
final1 = Path(pwd) / f"{species1}.cds"
final2 = Path(pwd) / f"{species2}.cds"
bed1 = Path(pwd) / f"{species1}.bed"
bed2 = Path(pwd) / f"{species2}.bed"

clean_fasta_headers(fasta1, raw1)
clean_fasta_headers(fasta2, raw2)

run_docker_jcvi_format(pwd, str(raw1), str(final1))
run_docker_jcvi_format(pwd, str(raw2), str(final2))

raw1.unlink()
raw2.unlink()

gff1 = prepare_gff(pwd, source_gff, species1)
gff2 = prepare_gff(pwd, target_gff, species2)
run_docker_jcvi_gff_to_bed(pwd, gff1, str(bed1))
run_docker_jcvi_gff_to_bed(pwd, gff2, str(bed2))

run_ortholog_analysis(pwd, species1, species2)
run_ortholog_analysis(pwd, species2, species1)