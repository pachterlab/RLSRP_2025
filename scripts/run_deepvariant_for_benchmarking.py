import argparse
import os
import sys
import subprocess
import shutil

import pysam

RLSRWP_2025_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #!!! erase
sys.path.append(RLSRWP_2025_dir)  #!!! erase
from RLSRWP_2025.seq_utils import perform_analysis

from varseek.utils import (
    run_command_with_error_logging,
    add_vcf_info_to_cosmic_tsv,
)

parser = argparse.ArgumentParser(description="Run Deepvariant on a set of reads and report the time and memory usage")

# Paths
parser.add_argument("--synthetic_read_fastq", help="Path to synthetic read FASTQ")
parser.add_argument("--reference_genome_fasta", help="Path to reference genome fasta")
parser.add_argument("--reference_genome_gtf", help="Path to reference genome GTF")
parser.add_argument("--star_genome_dir", default="", help="Path to star_genome_dir")
parser.add_argument("--aligned_and_unmapped_bam", default="", help="Path to aligned_and_unmapped_bam. If not provided, will be created")
parser.add_argument("--model_dir", default="model", help="Path to model files. If not provided, will be created")
parser.add_argument("--out", default="out", help="Path to out folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--cosmic_cdna_info_csv", help="Path to COSMIC csv with cdna info from gget cosmic")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")

args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.out, "star_genome")
deepvariant_output_dir = args.out
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = int(args.read_length) - 1
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_and_unmapped_bam = args.aligned_and_unmapped_bam
deepvariant_model = args.model_dir

STAR = args.STAR

cosmic_tsv = args.cosmic_tsv
cosmic_cdna_info_csv = args.cosmic_cdna_info_csv
unique_mcrs_df_path = args.unique_mcrs_df_path

os.makedirs(deepvariant_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)

BIN_VERSION="1.4.0"
alignment_folder = f"{deepvariant_output_dir}/alignment"
out_file_name_prefix = f"{alignment_folder}/sample_"

intermediate_data_dir = os.path.join(deepvariant_output_dir, "data_intermediate")
os.makedirs(intermediate_data_dir, exist_ok=True)
x3_coverage_bed_file = os.path.join(intermediate_data_dir, "x3_coverage.bed")  # used in run_bedtools.sh

output_dir = os.path.join(deepvariant_output_dir, "output")
os.makedirs(output_dir, exist_ok=True)

deepvariant_vcf = os.path.join(output_dir, "out.vcf.gz")

intermediate_results = os.path.join(deepvariant_output_dir, "intermediate_results_dir")
os.makedirs(intermediate_results, exist_ok=True)

if reference_genome_fasta[0] == "/":
    reference_genome_fasta_directory = os.path.dirname(reference_genome_fasta)
    reference_genome_fasta_directory_for_docker = "/reference_dir"
    reference_genome_fasta_for_docker = f"{reference_genome_fasta_directory_for_docker}/{reference_genome_fasta}"
else:
    reference_genome_fasta_for_docker = reference_genome_fasta

for name, path in {"STAR": STAR}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

#* File download
model_checkpoint_data_path = os.path.join(deepvariant_model, "model.ckpt.data-00000-of-00001")
model_checkpoint_example_path = os.path.join(deepvariant_model, "model.ckpt.example_info.json")
model_checkpoint_index_path = os.path.join(deepvariant_model, "model.ckpt.index")
model_checkpoint_meta_path = os.path.join(deepvariant_model, "model.ckpt.meta")
if not os.path.exists(model_checkpoint_data_path):
    subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.data-00000-of-00001 > {model_checkpoint_data_path}", check=True, shell=True)
if not os.path.exists(model_checkpoint_example_path):
    subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.example_info.json > {model_checkpoint_example_path}", check=True, shell=True)
if not os.path.exists(model_checkpoint_index_path):
    subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.index > {model_checkpoint_index_path}", check=True, shell=True)
if not os.path.exists(model_checkpoint_meta_path):
    subprocess.run(f"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.meta > {model_checkpoint_meta_path}", check=True, shell=True)

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
    "--outSAMunmapped", "Within",
    "--outSAMmapqUnique", "60",
    "--twopassMode", "Basic"
]
if not os.path.exists(aligned_and_unmapped_bam):
    aligned_and_unmapped_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
    os.makedirs(alignment_folder, exist_ok=True)
    run_command_with_error_logging(star_align_command)

#* BAM index file creation
bam_index_file = f"{aligned_and_unmapped_bam}.bai"
if not os.path.exists(bam_index_file):
    _ = pysam.index(aligned_and_unmapped_bam)

#* Filtering by 3x coverage regions of BAM file
#!!! these both use sudo too (here and in the run_bedtools.sh file)
generate_3x_coverage_files = [
    "docker", "run",
    "-v", f"{os.getcwd()}:{os.getcwd()}",
    "-w", os.getcwd(),
    "-it", "quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1",
    "mosdepth",
    "--threads", threads,
    f"{intermediate_data_dir}/data_coverage",
    aligned_and_unmapped_bam
]
# run_command_with_error_logging(generate_3x_coverage_files)   #!!! uncomment and debug
# run_command_with_error_logging(["bash", "run_bedtools.sh", intermediate_data_dir])

#* DeepVariant variant calling
deepvariant_command = [
    "docker", "run",
    "-v", f"{os.getcwd()}:{os.getcwd()}",
    "-w", os.getcwd(),
    f"google/deepvariant:{BIN_VERSION}",
    "run_deepvariant",
    "--model_type=WES",
    f"--customized_model={deepvariant_model}/model.ckpt",
    f"--ref={reference_genome_fasta_for_docker}",
    f"--reads={aligned_and_unmapped_bam}",
    f"--output_vcf={deepvariant_vcf}",
    f"--num_shards={os.cpu_count()}",
    # f"--regions={x3_coverage_bed_file}",   #!!! uncomment when getting the portion above working
    "--make_examples_extra_args=split_skip_reads=true,channels=''",
    "--intermediate_results_dir", intermediate_results
]
if reference_genome_fasta_for_docker != reference_genome_fasta:
    deepvariant_command[4:4] = ["-v", f"{reference_genome_fasta_directory}:{reference_genome_fasta_directory_for_docker}"]
run_command_with_error_logging(deepvariant_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

vcf_file = deepvariant_vcf
cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out = None, cosmic_cdna_info_csv = cosmic_cdna_info_csv, mutation_source = "cdna")
perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=deepvariant_output_dir, package_name="deepvariant")