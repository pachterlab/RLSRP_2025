import argparse
import os
import subprocess
import sys

import shutil

RLSRWP_2025_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #!!! erase
sys.path.append(RLSRWP_2025_dir)  #!!! erase
from RLSRWP_2025.seq_utils import perform_analysis

from varseek.utils import (
    run_command_with_error_logging,
    add_vcf_info_to_cosmic_tsv,
)

parser = argparse.ArgumentParser(description="Run Strelka2 on a set of reads and report the time and memory usage")

# Paths
parser.add_argument("--synthetic_read_fastq", help="Path to synthetic read FASTQ")
parser.add_argument("--reference_genome_fasta", help="Path to reference genome fasta")
parser.add_argument("--reference_genome_gtf", help="Path to reference genome GTF")
parser.add_argument("--star_genome_dir", default="", help="Path to star_genome_dir")
parser.add_argument("--aligned_bam", default="", help="Path to aligned_bam. If not provided, will be created")
parser.add_argument("--out", default="out", help="Path to out folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")
parser.add_argument("--STRELKA_INSTALL_PATH", default="STRELKA_INSTALL_PATH", help="Path to STRELKA_INSTALL_PATH parent")
parser.add_argument("--python2_env", default="python2_env", help="Conda environment name with python2")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--cosmic_cdna_info_csv", help="Path to COSMIC csv with cdna info from gget cosmic")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")


args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.out, "star_genome")
strelka2_output_dir = args.out
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = int(args.read_length) - 1
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_bam = args.aligned_bam

STAR = args.STAR
STRELKA_INSTALL_PATH = args.STRELKA_INSTALL_PATH
python2_env = args.python2_env

cosmic_tsv = args.cosmic_tsv
cosmic_cdna_info_csv = args.cosmic_cdna_info_csv
unique_mcrs_df_path = args.unique_mcrs_df_path

for name, path in {"STAR": STAR, "STRELKA_INSTALL_PATH": STRELKA_INSTALL_PATH}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

if not shutil.which("python2"):
    check_conda_environments_result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if python2_env not in check_conda_environments_result.stdout:
        raise FileNotFoundError(f"System python2 and conda environment {python2_env} not found. Please install system python2 or create a conda environment with python2 with `conda create -n {python2_env} python=2.7`, and supply it as an argument --python2_env {python2_env}.") 
    system_python2_installed = False
else:
    system_python2_installed = True

os.makedirs(strelka2_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)

alignment_folder = f"{strelka2_output_dir}/alignment"
os.makedirs(alignment_folder, exist_ok=True)
out_file_name_prefix = f"{alignment_folder}/sample_"
vcf_file = os.path.join(strelka2_output_dir, "results", "variants", "variants.vcf.gz")


# commented out, as these should already be done prior to running this script
#* STAR Build
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
    run_command_with_error_logging(star_build_command)

#* Reference genome index file
if not os.path.exists(f"{reference_genome_fasta}.fai"):
    import pysam
    _ = pysam.faidx(reference_genome_fasta)
    # commented out, as these should already be done prior to running this script

#* STAR Alignment
star_align_command = [
    STAR,
    "--runThreadN", str(threads),
    "--genomeDir", star_genome_dir,
    "--readFilesIn", synthetic_read_fastq,
    "--sjdbOverhang", str(read_length_minus_one),
    "--outFileNamePrefix", out_file_name_prefix,
    "--outSAMtype", "BAM", "SortedByCoordinate",
    "--outSAMmapqUnique", "60",
    "--twopassMode", "Basic"
]
if not os.path.exists(aligned_bam):
    aligned_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
    run_command_with_error_logging(star_align_command)

#* BAM index file creation
bam_index_file = f"{aligned_bam}.bai"
if not os.path.exists(bam_index_file):
    import pysam
    _ = pysam.index(aligned_bam)

#* Strelka2 variant calling configuration
strelka2_configure_command = [
    f"{STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py",
    "--bam", aligned_bam,
    "--referenceFasta", reference_genome_fasta,
    "--rna",
    "--runDir", strelka2_output_dir
]
if not system_python2_installed:
    strelka2_configure_command = ["conda", "run", "-n", python2_env] + strelka2_configure_command
if not os.path.exists(os.path.join(strelka2_output_dir, "runWorkflow.py")):
    run_command_with_error_logging(strelka2_configure_command)

#* Strelka2 variant calling configuration
strelka2_run_command = [
    f"{strelka2_output_dir}/runWorkflow.py",
    "-m", "local",
    "-j", str(threads)
]
if not system_python2_installed:
    strelka2_run_command = ["conda", "run", "-n", python2_env] + strelka2_run_command
run_command_with_error_logging(strelka2_run_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out = None, cosmic_cdna_info_csv = cosmic_cdna_info_csv, mutation_source = "cdna")
perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=strelka2_output_dir, package_name="strelka2")