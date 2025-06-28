import argparse
import os
import pandas as pd
import sys
import subprocess
import shutil
from pathlib import Path

import pysam

RLSRWP_2025_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #!!! erase
sys.path.append(RLSRWP_2025_dir)  #!!! erase
from RLSRWP_2025.seq_utils import perform_analysis, compare_two_vcfs_with_hap_py

scripts_dir = os.path.dirname(os.path.abspath(__file__))

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
parser.add_argument("--aligned_bam", default="", help="Path to aligned_bam. If not provided, will be created")
parser.add_argument("--model_dir", default="model", help="Path to model files. If not provided, will be created")
parser.add_argument("--out", default="out", help="Path to out folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--min_coverage", default=3, help="Min coverage. Set to 1 to forgo this filtering")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--cosmic_vcf", help="Path to COSMIC vcf")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")
parser.add_argument("--cosmic_version", default=101, help="COSMIC version. Default: 101")

args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.out, "star_genome")
deepvariant_output_dir = args.out
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = int(args.read_length) - 1
min_coverage = args.min_coverage
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_bam = args.aligned_bam
deepvariant_model = args.model_dir

STAR = args.STAR

cosmic_tsv = args.cosmic_tsv
unique_mcrs_df_path = args.unique_mcrs_df_path
cosmic_version = args.cosmic_version

os.makedirs(deepvariant_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)

BIN_VERSION="1.4.0"
alignment_folder = f"{deepvariant_output_dir}/alignment"
out_file_name_prefix = f"{alignment_folder}/sample_"
aligned_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam" if not aligned_bam else aligned_bam
aligned_bam = os.path.realpath(aligned_bam)  # ensure absolute path

output_dir = os.path.join(deepvariant_output_dir, "vcf_output")
os.makedirs(output_dir, exist_ok=True)

deepvariant_vcf = os.path.join(output_dir, "deepvariant_variants.vcf.gz")  # out.vcf.gz
deepvariant_vcf = os.path.realpath(deepvariant_vcf)  # ensure absolute path

plot_output_folder = f"{deepvariant_output_dir}/plots"
os.makedirs(plot_output_folder, exist_ok=True)

intermediate_results = os.path.join(deepvariant_output_dir, "intermediate_results_dir")
intermediate_results = os.path.realpath(intermediate_results)  # ensure absolute path
os.makedirs(intermediate_results, exist_ok=True)
above_min_coverage_bed_file = os.path.join(intermediate_results, "above_min_coverage.bed")  # used in run_bedtools.sh

reference_genome_fasta = os.path.realpath(reference_genome_fasta)  # ensure absolute path
deepvariant_model = os.path.realpath(deepvariant_model)  # ensure absolute path

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
if len(os.listdir(star_genome_dir)) == 0:
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
    "--outSAMmapqUnique", "60",
    "--twopassMode", "Basic"
]
if not os.path.exists(aligned_bam):
    os.makedirs(alignment_folder, exist_ok=True)
    run_command_with_error_logging(star_align_command)

#* BAM index file creation
bam_index_file = f"{aligned_bam}.bai"
if not os.path.exists(bam_index_file):
    _ = pysam.index(aligned_bam)

#* Filtering by 3x coverage regions of BAM file
generate_3x_coverage_files = [
    "docker", "run", "--rm",
    "-v", f"{os.getcwd()}:{os.getcwd()}",
    "-v", f"{intermediate_results}:{intermediate_results}",
    "-v", f"{os.path.dirname(aligned_bam)}:{os.path.dirname(aligned_bam)}",
    "-w", os.getcwd(),
    "-it", "quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1",
    "mosdepth",
    "--threads", threads,
    f"{intermediate_results}/data_coverage",
    aligned_bam
]
if min_coverage > 1 and not os.path.isfile(above_min_coverage_bed_file):
    if not os.path.exists(f"{intermediate_results}/data_coverage.per-base.bed.gz"):
        run_command_with_error_logging(generate_3x_coverage_files)
    run_command_with_error_logging(["bash", f"{scripts_dir}/run_bedtools.sh", intermediate_results, str(min_coverage)])


def find_independent_parents(paths):
    """
    Given a list of directory paths, returns only the top-level (parent) directories,
    removing any directories that are inside another.
    """
    cwd = os.getcwd()
    if cwd not in paths:
        paths.append(cwd)
    resolved_paths = {Path(p).resolve() for p in paths}  # Convert to absolute paths

    # Filter out child directories
    independent_parents = set()
    for path in resolved_paths:
        if not any(path.is_relative_to(other) for other in resolved_paths if path != other):
            independent_parents.add(str(path))

    return sorted(independent_parents)  # Sorted for consistency


#* DeepVariant variant calling
deepvariant_command = [
    "docker", "run", "--rm",
    # "-v", f"{os.getcwd()}:{os.getcwd()}",  # Added later
    "-w", os.getcwd(),
    f"google/deepvariant:{BIN_VERSION}",
    "run_deepvariant",
    "--model_type=WES",
    f"--customized_model={deepvariant_model}/model.ckpt",
    f"--ref={reference_genome_fasta}",
    f"--reads={aligned_bam}",
    f"--output_vcf={deepvariant_vcf}",
    f"--num_shards={os.cpu_count()}",
    "--make_examples_extra_args=split_skip_reads=true,channels=''",
    "--intermediate_results_dir", intermediate_results
]

directories_to_mount = find_independent_parents([os.path.dirname(reference_genome_fasta), os.path.dirname(deepvariant_model), os.path.dirname(aligned_bam), os.path.dirname(deepvariant_vcf), intermediate_results])
for directory in directories_to_mount:
    deepvariant_command[3:3] = ["-v", f"{directory}:{directory}"]

if min_coverage > 1 and os.path.isfile(above_min_coverage_bed_file):
    deepvariant_command.insert(-4, f"--regions={above_min_coverage_bed_file}")

if not os.path.isfile(deepvariant_vcf):
    run_command_with_error_logging(deepvariant_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

# cosmic_df_out = cosmic_tsv.replace(".tsv", "_vcf_info_for_fig2.csv")
# if not os.path.exists(cosmic_df_out):
#     cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out=cosmic_df_out, sequences="cdna", cosmic_version=cosmic_version)
# else:
#     cosmic_df = pd.read_csv(cosmic_df_out)

vcf_file = deepvariant_vcf
# perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=plot_output_folder, package_name="deepvariant", dp_column="default_DP")

package_name = "deepvariant"
cosmic_vcf = args.cosmic_vcf
happy_out = os.path.join(args.out, "hap_py_out", package_name)

cosmic_vcf = os.path.realpath(cosmic_vcf)
vcf_file = os.path.realpath(vcf_file)
reference_genome_fasta = os.path.realpath(reference_genome_fasta)
happy_out = os.path.realpath(happy_out)
compare_two_vcfs_with_hap_py(ground_truth_vcf=cosmic_vcf, test_vcf=vcf_file, reference_fasta=reference_genome_fasta, output_dir = happy_out, unique_mcrs_df = unique_mcrs_df_path, unique_mcrs_df_out = None, package_name = package_name, dry_run = False)