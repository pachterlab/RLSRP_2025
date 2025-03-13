import argparse
import os
import sys
import shutil

import pysam

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
deepvariant_output_dir = os.path.join(args.out, "deepvariant_simulated_data_dir")
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = args.read_length - 1
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_and_unmapped_bam = args.aligned_and_unmapped_bam

STAR = args.STAR

cosmic_tsv = args.cosmic_tsv
cosmic_cdna_info_csv = args.cosmic_cdna_info_csv
unique_mcrs_df_path = args.unique_mcrs_df_path

os.makedirs(deepvariant_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)

alignment_folder = f"{deepvariant_output_dir}/alignment"
out_file_name_prefix = f"{alignment_folder}/sample_"

deepvariant_vcf = os.path.join(deepvariant_output_dir, "results/variants/genome.vcf.gz")

intermediate_results = os.path.join(deepvariant_output_dir, "intermediate_results_dir")
os.makedirs(intermediate_results, exist_ok=True)

for name, path in {"STAR": STAR}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

#* Genome alignment with STAR
star_build_command = [
    STAR,
    "--runThreadN", str(threads),
    "--runMode", "genomeGenerate",
    "--genomeDir", star_genome_dir,
    "--genomeFastaFiles", reference_genome_fasta,
    "--sjdbGTFfile", reference_genome_gtf,
    "--sjdbOverhang", str(read_length_minus_one),
]

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

#* deepvariant variant calling
BIN_VERSION="1.4.0"

deepvariant_command = [
    "sudo", "docker", "run",
    "-v", f"{os.getcwd()}:{os.getcwd()}",
    "-w", os.getcwd(),
    f"google/deepvariant:{BIN_VERSION}",
    "run_deepvariant",
    "--model_type=WES",
    "--customized_model=model/model.ckpt",
    f"--ref={reference_genome_fasta}",
    f"--reads={aligned_and_unmapped_bam}",
    f"--output_vcf={deepvariant_vcf}",
    f"--num_shards={os.cpu_count()}",
    "--make_examples_extra_args=split_skip_reads=true,channel_list='BASE_CHANNELS'",
    "--intermediate_results_dir", intermediate_results
]

# commented out, as these should already be done prior to running this script
if not os.listdir(star_genome_dir):
    run_command_with_error_logging(star_build_command)

if not os.path.exists(f"{reference_genome_fasta}.fai"):
    _ = pysam.faidx(reference_genome_fasta)
# commented out, as these should already be done prior to running this script

if not os.path.exists(aligned_and_unmapped_bam):
    aligned_and_unmapped_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
    os.makedirs(alignment_folder, exist_ok=True)
    run_command_with_error_logging(star_align_command)

bam_index_file = f"{aligned_and_unmapped_bam}.bai"
if not os.path.exists(bam_index_file):
    _ = pysam.index(aligned_and_unmapped_bam)

run_command_with_error_logging(deepvariant_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

vcf_file = deepvariant_vcf
cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out = None, cosmic_cdna_info_csv = cosmic_cdna_info_csv, mutation_source = "cdna")
perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=deepvariant_output_dir, package_name="deepvariant")