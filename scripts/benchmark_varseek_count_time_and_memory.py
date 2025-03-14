import os
import random
import gget
import subprocess
import time
import logging

import numpy as np
import pandas as pd
import pysam

import varseek as vk
from varseek.utils import (
    is_program_installed,
    report_time_and_memory_of_script,
    run_command_with_error_logging,
    convert_mutation_cds_locations_to_cdna
)

logger = logging.getLogger(__name__)
logger.setLevel("INFO")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S")
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
output_dir = os.path.join(data_dir, "benchmarking_out_dir")
vk_sim_out_dir = os.path.join(output_dir, "vk_sim_out")
reference_out_dir = os.path.join(data_dir, "reference")

all_supported_tools_to_benchmark = {"varseek", "gatk_haplotypecaller", "gatk_mutect2", "strelka2", "varscan", "deepvariant"}
tools_that_require_star_alignment = {"gatk_haplotypecaller", "gatk_mutect2", "strelka2", "varscan", "deepvariant"}

### ARGUMENTS ###
number_of_reads_list = [0.001, 0.002]  # [1, 4, 16, 64, 256, 1024]  # number of reads, in millions
tools_to_benchmark = ["varseek", "gatk_haplotypecaller", "gatk_mutect2", "strelka2", "varscan", "deepvariant"]
dry_run = False  # only applies to the variant calling steps, not to the preparation (ie STAR, data downloads, etc)

read_length = 150
k = 51
w = 47
strand = None  # None for strand-agnostic (randomly-selected), "f" for forward, "r" for reverse, "both" for both - make sure this matches the reference genome (vk build command) - strand = True -> "f" or "r" here; strand = False -> None or "both" here - note that the strand is randomly selected per *transcript*, such that all drawn reads will come from the same strand no matter what
add_noise_sequencing_error = True
add_noise_base_quality = False
error_rate = 0.0001  # only if add_noise_sequencing_error=True
error_distribution = (0.85, 0.1, 0.05)  # sub, del, ins  # only if add_noise_sequencing_error=True
max_errors = float("inf")  # only if add_noise_sequencing_error=True
seq_id_column = "seq_ID"
var_column = "mutation_cdna"
threads = 16
random_seed = 42
qc_against_gene_matrix = False

# varseek ref parameters
vk_ref_out = os.path.join(data_dir, "vk_ref_out")
vk_ref_index_path = os.path.join(vk_ref_out, "vcrs_index.idx")  # for vk count
vk_ref_t2g_path = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")  # for vk count
dlist_reference_source = "t2t"

# normal reference genome
reference_genome_index_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "index.idx")  # can either already exist or will be created; only used if qc_against_gene_matrix=True
reference_genome_t2g_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "t2g.txt")  # can either already exist or will be created; only used if qc_against_gene_matrix=True

# varseek ref small k
w_small = 27
k_small = 31
vk_ref_small_k_out = os.path.join(data_dir, "vk_ref_small_k_out")
vk_ref_small_k_index_path = os.path.join(vk_ref_out, "vcrs_index.idx")  # for vk count
vk_ref_small_k_t2g_path = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")  # for vk count


cosmic_mutations_path = os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")  # for vk sim
reference_cdna_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.cdna.all.fa")  # for vk sim
reference_cds_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.cds.all.fa")  # for vk sim
reference_genome_fasta = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")  # can either already exist or will be downloaded
reference_genome_gtf = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.87.gtf")  # can either already exist or will be downloaded
genomes1000_vcf = os.path.join(reference_out_dir, "ensembl_grch37_release93", "1000GENOMES-phase_3.vcf")
star_genome_dir = os.path.join(reference_out_dir, "ensembl_grch37_release93", "star_reference")

seqtk = "seqtk"
STAR = "STAR"  # "/home/jmrich/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR"
java = "/home/jmrich/opt/jdk-17.0.12+7/bin/java"
picard_jar = "/home/jmrich/opt/picard.jar"
gatk = "/home/jmrich/opt/gatk-4.6.0.0/gatk"
STRELKA_INSTALL_PATH = "/home/jmrich/opt/strelka-2.9.10.centos6_x86_64"
python2_env = "python2_env"
VARSCAN_INSTALL_PATH = "/home/jmrich/opt/VarScan.v2.3.9.jar"

tmp_dir = "tmp"

deepvariant_model = os.path.join(tmp_dir, "deepvariant_model")
model_checkpoint_data_path = os.path.join(deepvariant_model, "model.ckpt.data-00000-of-00001")
model_checkpoint_example_path = os.path.join(deepvariant_model, "model.ckpt.example_info.json")
model_checkpoint_index_path = os.path.join(deepvariant_model, "model.ckpt.index")
model_checkpoint_meta_path = os.path.join(deepvariant_model, "model.ckpt.meta")
### ARGUMENTS ###

# # set random seeds
# random.seed(random_seed)
# np.random.seed(random_seed)

os.makedirs(tmp_dir)  # purposely not using exist_ok=True to ensure that the directory is non-existent  #* comment out for debugging to keep tmp_dir between runs
os.makedirs(deepvariant_model, exist_ok=True)

if not star_genome_dir:
    star_genome_dir = os.path.join(tmp_dir, "star_reference")

if strand is None or strand == "both":
    kb_count_strand = "unstranded"
elif strand == "f":
    kb_count_strand = "forward"
elif strand == "r":
    kb_count_strand = "reverse"

vk_count_script_path = os.path.join(script_dir, "run_varseek_count_for_benchmarking.py")
gatk_haplotypecaller_script_path = os.path.join(script_dir, "run_gatk_haplotypecaller_for_benchmarking.py")
gatk_mutect2_script_path = os.path.join(script_dir, "run_gatk_mutect2_for_benchmarking.py")
strelka_script_path = os.path.join(script_dir, "run_strelka_for_benchmarking.py")
varscan_script_path = os.path.join(script_dir, "run_varscan_for_benchmarking.py")
deepvariant_script_path = os.path.join(script_dir, "run_deepvariant_for_benchmarking.py")

star_output_file = os.path.join(output_dir, "STAR_time.txt")

# create synthetic reads
if k and w:
    if k <= w:
        raise ValueError("k must be greater than w")
    read_w = read_length - (k - w)  # note that this does not affect read length, just read *parent* length
else:
    read_w = read_length - 1


for tool in tools_to_benchmark:
    if tool not in all_supported_tools_to_benchmark:
        raise ValueError(f"Tool {tool} is not supported. Supported tools are: {all_supported_tools_to_benchmark}")

#* download COSMIC and sequences for vk sim if not already downloaded
# download cosmic and cdna
if not os.path.exists(reference_cdna_path):
    logger.info("Downloading cDNA")
    reference_cdna_dir = os.path.dirname(reference_cdna_path) if os.path.dirname(reference_cdna_path) else "."
    gget_ref_command = ["gget", "ref", "-w", "cdna", "-r", "93", "--out_dir", reference_cdna_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    subprocess.run(["gunzip", f"{reference_cdna_path}.gz"], check=True)
if not os.path.exists(reference_cds_path):
    logger.info("Downloading CDS")
    reference_cds_dir = os.path.dirname(reference_cds_path) if os.path.dirname(reference_cds_path) else "."
    gget_ref_command = ["gget", "ref", "-w", "cds", "-r", "93", "--out_dir", reference_cds_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    subprocess.run(["gunzip", f"{reference_cds_path}.gz"], check=True)

if not os.path.exists(cosmic_mutations_path):
    logger.info("Downloading COSMIC")
    gget.cosmic(
        None,
        grch_version=37,
        cosmic_version=101,
        out=reference_out_dir,
        mutation_class="cancer",
        download_cosmic=True,
    )

cosmic_mutations = pd.read_csv(cosmic_mutations_path)

if "mutation_cdna" not in cosmic_mutations.columns:
    logger.info("Converting CDS to cDNA in COSMIC")
    cosmic_mutations, _ = convert_mutation_cds_locations_to_cdna(input_csv_path=cosmic_mutations, output_csv_path=cosmic_mutations_path, cds_fasta_path=reference_cds_path, cdna_fasta_path=reference_cdna_path, logger=logger, verbose=True)

if "header" not in cosmic_mutations.columns:
    logger.info("Making header column in COSMIC")
    cosmic_mutations["header"] = cosmic_mutations["seq_ID"] + ":" + cosmic_mutations["mutation_cdna"]
    cosmic_mutations.to_csv(cosmic_mutations_path, index=False)

#* Make synthetic reads corresponding to the largest value in number_of_reads_list - if desired, I can replace this with real data
number_of_reads_max = int(max(number_of_reads_list) * 10**6)  # convert to millions

number_of_reads_per_variant_alt=100
number_of_reads_per_variant_ref=150
number_of_reads_per_variant_total = number_of_reads_per_variant_alt + number_of_reads_per_variant_ref

if number_of_reads_max > (len(cosmic_mutations) * number_of_reads_per_variant_total):
    raise ValueError("Max reads is too large. Either increase number_of_reads_per_variant_alt and/or number_of_reads_per_variant_ref, or choose a larger variant database.")

number_of_variants_to_sample = number_of_reads_max // number_of_reads_per_variant_total
fastq_output_path_max_reads = os.path.join(tmp_dir, f"reads_{number_of_reads_max}_fastq.fastq")

if not os.path.isfile(fastq_output_path_max_reads):
    logger.info(f"Building synthetic reads for {number_of_reads_max} reads")
    _ = vk.sim(
        variants = cosmic_mutations,
        reads_fastq_out = fastq_output_path_max_reads,
        number_of_variants_to_sample=number_of_variants_to_sample,
        strand=strand,
        number_of_reads_per_variant_alt=number_of_reads_per_variant_alt,
        number_of_reads_per_variant_ref=number_of_reads_per_variant_ref,
        read_length=read_length,
        seed=random_seed,
        add_noise_sequencing_error=add_noise_sequencing_error,
        add_noise_base_quality=add_noise_base_quality,
        error_rate=error_rate,
        error_distribution=error_distribution,
        max_errors=max_errors,
        with_replacement=False,
        gzip_reads_fastq_out=False,
        sequences=reference_cdna_path,
        seq_id_column=seq_id_column,
        var_column=var_column,
        header_column="header",
        variant_type_column=None,
        reference_out_dir=reference_out_dir,
        out=vk_sim_out_dir,
        k=k,
        w=w
    )

#* Download varseek index
if not os.path.exists(vk_ref_index_path) or not os.path.exists(vk_ref_t2g_path):
    vk.ref(variants="cosmic_cmc", sequences="cdna", w=w, k=k, out=vk_ref_out, dlist_reference_source=dlist_reference_source, download=True, index_out=vk_ref_index_path, t2g_out=vk_ref_t2g_path)
    # alternatively, to build from scratch: subprocess.run([os.path.join(script_dir, "run_vk_ref.py")], check=True)

if not os.path.exists(vk_ref_small_k_index_path) or not os.path.exists(vk_ref_small_k_t2g_path):
    try:
        vk.ref(variants="cosmic_cmc", sequences="cdna", w=w_small, k=k_small, out=vk_ref_out, dlist_reference_source=dlist_reference_source, download=True, index_out=vk_ref_small_k_index_path, t2g_out=vk_ref_small_k_t2g_path)
    except ValueError:
        logger.info(f"Cannot download vk ref index/t2g with w={w_small} and k={k_small}. Will skip this condition")

#* install seqtk if not installed
if not is_program_installed(seqtk):
    raise ValueError("seqtk is required to run this script. Please install seqtk and ensure that it is in your PATH.")
    # subprocess.run("git clone https://github.com/lh3/seqtk.git", shell=True, check=True)
    # subprocess.run("cd seqtk && make", shell=True, check=True)
    # seqtk = os.path.join(script_dir, "seqtk/seqtk")

#* Build normal genome reference (for vk clean in vk count) when qc_against_gene_matrix=True
if qc_against_gene_matrix and (not os.path.exists(reference_genome_index_path) or not os.path.exists(reference_genome_t2g_path)):  # download reference if does not exist
    if not os.path.exists(reference_genome_fasta) or not os.path.exists(reference_genome_gtf):
        reference_genome_out_dir = os.path.dirname(reference_genome_fasta) if os.path.dirname(reference_genome_fasta) else "."
        subprocess.run(["gget", "ref", "-w", "dna,gtf", "-r", "93", "--out_dir", reference_genome_out_dir, "-d", "human_grch37"], check=True)  # using grch37, ensembl 93 to agree with COSMIC
        subprocess.run(["gunzip", f"{reference_genome_fasta}.gz"], check=True)
        subprocess.run(["gunzip", f"{reference_genome_gtf}.gz"], check=True)
    reference_genome_f1 = os.path.join(reference_out_dir, "ensembl_grch37_release93", "f1.fasta")
    subprocess.run(["kb", "ref", "-t", str(threads), "-i", reference_genome_index_path, "-g", reference_genome_t2g_path, "-f1", reference_genome_f1, reference_genome_fasta, reference_genome_gtf], check=True)

#* Download/build necessary files for alternative variant calling tools
if any(tool in tools_that_require_star_alignment for tool in tools_to_benchmark):  # check if any tool in tools_to_benchmark requires STAR alignment
    #* Download reference genome information
    # reference_genome_fasta_url = "https://ftp.ensembl.org/pub/grch37/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
    # reference_genome_gtf_url = "https://ftp.ensembl.org/pub/grch37/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
    genomes1000_vcf_url = "https://ftp.ensembl.org/pub/grch37/release-93/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"

    reference_genome_out_dir = os.path.dirname(reference_genome_fasta) if os.path.dirname(reference_genome_fasta) else "."
    os.makedirs(reference_genome_out_dir, exist_ok=True)
    download_reference_genome_fasta_command = ["gget", "ref", "-w", "dna", "-r", "93", "--out_dir", reference_genome_out_dir, "-d", "human_grch37"]
    unzip_reference_genome_fasta_command = ["gunzip", f"{reference_genome_fasta}.gz"]

    reference_genome_out_dir = os.path.dirname(reference_genome_gtf) if os.path.dirname(reference_genome_gtf) else "."
    os.makedirs(reference_genome_out_dir, exist_ok=True)
    download_reference_genome_gtf_command = ["gget", "ref", "-w", "gtf", "-r", "93", "--out_dir", reference_genome_gtf, "-d", "human_grch37"]
    unzip_reference_genome_gtf_command = ["gunzip", f"{reference_genome_gtf}.gz"]

    if os.path.dirname(genomes1000_vcf):
        os.makedirs(os.path.dirname(genomes1000_vcf), exist_ok=True)
    download_1000_genomes_command = ["wget", "-O", f"{genomes1000_vcf}.gz", genomes1000_vcf_url]
    unzip_1000_genomes_command = ["gunzip", f"{genomes1000_vcf}.gz"]

    if not os.path.exists(reference_genome_fasta):
        run_command_with_error_logging(download_reference_genome_fasta_command)
        run_command_with_error_logging(unzip_reference_genome_fasta_command)

    if not os.path.exists(reference_genome_gtf):
        run_command_with_error_logging(download_reference_genome_gtf_command)
        run_command_with_error_logging(unzip_reference_genome_gtf_command)

    if not os.path.exists(genomes1000_vcf):
        run_command_with_error_logging(download_1000_genomes_command)
        run_command_with_error_logging(unzip_1000_genomes_command)

    read_length_minus_one = read_length - 1
    os.makedirs(star_genome_dir, exist_ok=True)

    if STAR != "STAR" and not os.path.exists(STAR):
        raise ValueError("STAR is required to run STAR. Please install STAR and, if installed from source, ensure that it is in your PATH.")
        # star_tarball = os.path.join(opt_dir, "2.7.11b.tar.gz")
        # subprocess.run(["wget", "-O", star_tarball, "https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz"], check=True)
        # subprocess.run(["tar", "-xzf", star_tarball, "-C", opt_dir], check=True)

    #* Build STAR index
    star_build_command = [
        STAR,
        "--runThreadN", str(threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", star_genome_dir,
        "--genomeFastaFiles", reference_genome_fasta,
        "--sjdbGTFfile", reference_genome_gtf,
        "--sjdbOverhang", str(read_length_minus_one),
    ]

    if not os.listdir(star_genome_dir):
        elapsed_time = run_command_with_error_logging(star_build_command, track_time=True)
        with open(star_output_file, "w", encoding="utf-8") as f:
            f.write(f"STAR build runtime: {elapsed_time[0]} minutes, {elapsed_time[1]} seconds\n")

    #* Index reference genome
    if not os.path.exists(f"{reference_genome_fasta}.fai"):
        start_time = time.time()
        _ = pysam.faidx(reference_genome_fasta)
        minutes, seconds = divmod(time.time() - start_time, 60)
        with open(star_output_file, "a", encoding="utf-8") as f:
            f.write(f"Genome indexing runtime: {minutes} minutes, {seconds} seconds\n")

    #* Index 1000 genomes standard variants
    index_feature_file_command = [
        gatk, "IndexFeatureFile",
        "-I", genomes1000_vcf
    ]

    if not os.path.exists(f"{genomes1000_vcf}.idx"):
        run_command_with_error_logging(index_feature_file_command)

if "deepvariant" in tools_to_benchmark:
    if not os.path.exists(model_checkpoint_data_path):
        subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.data-00000-of-00001 > {model_checkpoint_data_path}", check=True, shell=True)
    
    if not os.path.exists(model_checkpoint_example_path):
        subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.example_info.json > {model_checkpoint_example_path}", check=True, shell=True)

    if not os.path.exists(model_checkpoint_index_path):
        subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.index > {model_checkpoint_index_path}", check=True, shell=True)

    if not os.path.exists(model_checkpoint_meta_path):
        subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.meta > {model_checkpoint_meta_path}", check=True, shell=True)

if ("gatk_haplotypecaller" in tools_to_benchmark or "gatk_mutect2" in tools_to_benchmark):
    if not is_program_installed(java):
        raise ValueError("Java is required to run GATK. Please install Java and ensure that it is in your PATH.")
        # below are commands that will install it on x64_linux
        # wget https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.12%2B7/OpenJDK17U-jdk_x64_linux_hotspot_17.0.12_7.tar.gz
        # tar -xvf OpenJDK17U-jdk_x64_linux_hotspot_17.0.12_7.tar.gz
    if not os.path.exists(picard_jar):
        raise ValueError("Picard is required to run GATK. Please install Picard and ensure that it is in your PATH.")
        # subprocess.run(["wget", "https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar", "-O", picard_jar], check=True)
    if not os.path.exists(gatk):
        raise ValueError("GATK is required to run GATK. Please install GATK and ensure that it is in your PATH.")
        # os.chdir(opt_dir)
        # gatk_dir_name = os.path.join(opt_dir, "gatk-4.6.0.0")
        # subprocess.run(["wget", "https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip", "-O", "gatk-4.6.0.0.zip"], check=True)
        # subprocess.run(["unzip", "gatk-4.6.0.0.zip"], check=True)
        # os.environ['PATH'] = f"{gatk_dir_name}:{os.environ['PATH']}"


if "strelka2" in tools_to_benchmark:
    if not os.path.exists(STRELKA_INSTALL_PATH):
        raise ValueError("Strelka2 is required to run Strelka2. Please install Strelka2 and ensure that it is in your PATH.")
        # strelka_tarball = f"{STRELKA_INSTALL_PATH}.tar.bz2"
        # subprocess.run(["wget", "-O", strelka_tarball, "https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2"], check=True)
        # subprocess.run(["tar", "-xvjf", strelka_tarball, "-C", opt_dir], check=True)

    import shutil
    if not shutil.which("python2"):
        check_conda_environments_result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
        if python2_env not in check_conda_environments_result.stdout:
            raise FileNotFoundError(f"System python2 and conda environment {python2_env} not found. Please install system python2 or create a conda environment with python2 with `conda create -n {python2_env} python=2.7`, and supply it as an argument --python2_env {python2_env}.") 

if "varscan" in tools_to_benchmark and not os.path.exists(VARSCAN_INSTALL_PATH):
    raise ValueError("VarScan is required to run VarScan. Please install VarScan and ensure that it is in your PATH.")
    # subprocess.run(["wget", "-O", VARSCAN_INSTALL_PATH, "https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download"], check=True)

output_file = os.path.join(output_dir, "time_and_memory_benchmarking_report.txt")

#* Run variant calling tools
for number_of_reads in number_of_reads_list:
    number_of_reads = int(number_of_reads * 10**6)  # convert to millions
    fastq_output_path = os.path.join(tmp_dir, f"reads_{number_of_reads}_fastq.fastq")
    if number_of_reads != number_of_reads_max:
        seqtk_sample_command = f"{seqtk} sample -s {random_seed} {fastq_output_path_max_reads} {number_of_reads} > {fastq_output_path}"
        subprocess.run(seqtk_sample_command, shell=True, check=True)

    kb_count_reference_genome_out_dir = os.path.join(tmp_dir, f"kb_count_reference_genome_out_dir_{number_of_reads}")
    if not os.path.exists(kb_count_reference_genome_out_dir):
        # kb count, reference genome
        kb_count_standard_index_command = [
            "kb",
            "count",
            "-t",
            str(threads),
            "-i",
            reference_genome_index_path,
            "-g",
            reference_genome_t2g_path,
            "-x",
            "bulk",
            "--h5ad",
            "--parity",
            "single",
            "--strand",
            kb_count_strand,
            "-o",
            kb_count_reference_genome_out_dir,
            fastq_output_path
        ]

        logger.info(f"kb count, reference genome, {number_of_reads} reads")
        subprocess.run(kb_count_standard_index_command, check=True)
            
    #* Variant calling: varseek
    if "varseek" in tools_to_benchmark:
        logger.info(f"varseek, {number_of_reads} reads")
        script_title = f"varseek {number_of_reads} reads {threads} threads"
        vk_count_out_tmp = os.path.join(tmp_dir, f"vk_count_threads_{number_of_reads}_reads_out")
        argparse_flags = f"--index {vk_ref_index_path} --t2g {vk_ref_t2g_path} --technology bulk --threads {threads} -k {k} --out {vk_count_out_tmp} --kb_count_reference_genome_out_dir {kb_count_reference_genome_out_dir} --disable-clean --disable_summarize {fastq_output_path}"
        if dry_run:
            print(f"python3 {vk_count_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(vk_count_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    #* Variant calling: varseek with smaller k
    if "varseek" in tools_to_benchmark and os.path.isfile(vk_ref_small_k_index_path) and os.path.isfile(vk_ref_small_k_t2g_path):
        logger.info(f"varseek k={k_small}, {number_of_reads} reads")
        script_title = f"varseek_k={k_small} {number_of_reads} reads {threads} threads"
        vk_count_out_tmp = os.path.join(tmp_dir, f"vk_count_k{k_small}_reads_{number_of_reads}_out")
        argparse_flags = f"--index {vk_ref_small_k_index_path} --t2g {vk_ref_small_k_t2g_path} --technology bulk --threads {threads} -k {k_small} --out {vk_count_out_tmp} -k {k_small} --kb_count_reference_genome_out_dir {kb_count_reference_genome_out_dir} --disable-clean --disable_summarize {fastq_output_path}"
        if dry_run:
            print(f"python3 {vk_count_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(vk_count_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    if any(tool in tools_that_require_star_alignment for tool in tools_to_benchmark):
        #* STAR alignment
        star_alignment_dir = os.path.join(tmp_dir, f"star_alignment_dir_{number_of_reads}")
        out_file_name_prefix = f"{star_alignment_dir}/sample_"
        aligned_and_unmapped_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
        star_align_command = [
            STAR,
            "--runThreadN", str(threads),
            "--genomeDir", star_genome_dir,
            "--readFilesIn", fastq_output_path,
            "--sjdbOverhang", str(read_length_minus_one),
            "--outFileNamePrefix", out_file_name_prefix,
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMunmapped", "Within",
            "--outSAMmapqUnique", "60",
            "--twopassMode", "Basic"
        ]

        logger.info(f"STAR alignment, {number_of_reads} reads")
        if not os.path.isfile(aligned_and_unmapped_bam):
            elapsed_time = run_command_with_error_logging(star_align_command, track_time=True)
            with open(star_output_file, "a", encoding="utf-8") as f:
                f.write(f"STAR alignment runtime for {number_of_reads} reads: {elapsed_time[0]} minutes, {elapsed_time[1]} seconds\n")

        #* Index BAM file
        bam_index_file = f"{aligned_and_unmapped_bam}.bai"
        if not os.path.isfile(bam_index_file):
            start_time = time.time()
            _ = pysam.index(aligned_and_unmapped_bam)
            minutes, seconds = divmod(time.time() - start_time, 60)
            with open(star_output_file, "a", encoding="utf-8") as f:
                f.write(f"BAM indexing runtime for {number_of_reads} reads: {minutes} minutes, {seconds} seconds\n")
    
    if "gatk_haplotypecaller" in tools_to_benchmark:
        #* Variant calling: GATK HaplotypeCaller
        logger.info(f"GATK HaplotypeCaller, {number_of_reads} reads")
        script_title = f"gatk_haplotypecaller {number_of_reads} reads {threads} threads"
        gatk_parent_haplotypecaller = os.path.join(tmp_dir, f"gatk_haplotypecaller_{number_of_reads}_reads")
        argparse_flags = f"--synthetic_read_fastq {fastq_output_path} --reference_genome_fasta {reference_genome_fasta} --reference_genome_gtf {reference_genome_gtf} --genomes1000_vcf {genomes1000_vcf} --star_genome_dir {star_genome_dir} --aligned_and_unmapped_bam {aligned_and_unmapped_bam} --out {gatk_parent_haplotypecaller} --threads {threads} --read_length {read_length} --STAR {STAR} --java {java} --picard_jar {picard_jar} --gatk {gatk} --skip_accuracy_analysis"
        if dry_run:
            print(f"python3 {gatk_haplotypecaller_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(gatk_haplotypecaller_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    if "gatk_mutect2" in tools_to_benchmark:
        #* Variant calling: GATK Mutect2
        logger.info(f"GATK Mutect2, {number_of_reads} reads")
        script_title = f"gatk_mutect2 {number_of_reads} reads {threads} threads"
        gatk_parent_mutect2 = os.path.join(tmp_dir, f"gatk_mutect2_{number_of_reads}_reads")
        argparse_flags = f"--synthetic_read_fastq {fastq_output_path} --reference_genome_fasta {reference_genome_fasta} --reference_genome_gtf {reference_genome_gtf} --genomes1000_vcf {genomes1000_vcf} --star_genome_dir {star_genome_dir} --aligned_and_unmapped_bam {aligned_and_unmapped_bam} --out {gatk_parent_mutect2} --threads {threads} --read_length {read_length} --STAR {STAR} --java {java} --picard_jar {picard_jar} --gatk {gatk} --skip_accuracy_analysis"
        if dry_run:
            print(f"python3 {gatk_mutect2_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(gatk_mutect2_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    if "strelka2" in tools_to_benchmark:
        #* Variant calling: Strelka2
        logger.info(f"Strelka2, {number_of_reads} reads")
        script_title = f"strelka2 {number_of_reads} reads {threads} threads"
        strelka2_output_dir = os.path.join(tmp_dir, f"strelka2_simulated_data_dir_{number_of_reads}_reads")
        argparse_flags = f"--synthetic_read_fastq {fastq_output_path} --reference_genome_fasta {reference_genome_fasta} --reference_genome_gtf {reference_genome_gtf} --star_genome_dir {star_genome_dir} --aligned_and_unmapped_bam {aligned_and_unmapped_bam} --out {strelka2_output_dir} --threads {threads} --read_length {read_length} --STRELKA_INSTALL_PATH {STRELKA_INSTALL_PATH} --python2_env {python2_env} --skip_accuracy_analysis"
        if dry_run:
            print(f"python3 {strelka_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(strelka_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    if "varscan" in tools_to_benchmark:
        #* Variant calling: VarScan
        logger.info(f"VarScan, {number_of_reads} reads")
        script_title = f"varscan {number_of_reads} reads {threads} threads"
        varscan_output_dir = os.path.join(tmp_dir, f"varscan_simulated_data_dir_{number_of_reads}_reads")
        argparse_flags = f"--synthetic_read_fastq {fastq_output_path} --reference_genome_fasta {reference_genome_fasta} --reference_genome_gtf {reference_genome_gtf} --star_genome_dir {star_genome_dir} --aligned_and_unmapped_bam {aligned_and_unmapped_bam} --out {varscan_output_dir} -threads {threads} --read_length {read_length} --VARSCAN_INSTALL_PATH {VARSCAN_INSTALL_PATH} --skip_accuracy_analysis"
        if dry_run:
            print(f"python3 {varscan_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(varscan_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

    if "deepvariant" in tools_to_benchmark:
        #* Variant calling: Deepvariant
        logger.info(f"Deepvariant, {number_of_reads} reads")
        script_title = f"deepvariant {number_of_reads} reads {threads} threads"
        deepvariant_output_dir = os.path.join(tmp_dir, f"deepvariant_simulated_data_dir_{number_of_reads}_reads")
        argparse_flags = f"--synthetic_read_fastq {fastq_output_path} --reference_genome_fasta {reference_genome_fasta} --reference_genome_gtf {reference_genome_gtf} --star_genome_dir {star_genome_dir} --threads {threads} --read_length {read_length} --out {deepvariant_output_dir} --aligned_and_unmapped_bam {aligned_and_unmapped_bam} --model_dir {deepvariant_model} --threads {threads} --read_length {read_length} --skip_accuracy_analysis"
        if dry_run:
            print(f"python3 {deepvariant_script_path} {argparse_flags}")
        else:
            _ = report_time_and_memory_of_script(deepvariant_script_path, output_file = output_file, argparse_flags = argparse_flags, script_title = script_title)

# delete tmp directory
# os.system(f"rm -rf {tmp_dir}")  #!!! uncomment later to delete tmp directory
