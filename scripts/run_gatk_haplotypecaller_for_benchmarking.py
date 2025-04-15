import argparse
import os
import sys
import shutil

import pysam
import pandas as pd

from varseek.utils import (
    run_command_with_error_logging,
    add_vcf_info_to_cosmic_tsv,
)

RLSRWP_2025_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #!!! erase
sys.path.append(RLSRWP_2025_dir)  #!!! erase
from RLSRWP_2025.seq_utils import perform_analysis

parser = argparse.ArgumentParser(description="Run GATK Haplotypecaller on a set of reads and report the time and memory usage")

# Paths
parser.add_argument("--synthetic_read_fastq", help="Path to synthetic read FASTQ")
parser.add_argument("--reference_genome_fasta", help="Path to reference genome fasta")
parser.add_argument("--reference_genome_gtf", help="Path to reference genome GTF")
parser.add_argument("--genomes1000_vcf", default="1000GENOMES-phase_3.vcf", help="Path to 1000 genomes vcf file")
parser.add_argument("--star_genome_dir", default="", help="Path to star_genome_dir")
parser.add_argument("--aligned_and_unmapped_bam", default="", help="Path to aligned_and_unmapped_bam. If not provided, will be created")
parser.add_argument("--out", default="out", help="Path to out folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--apply_mutation_filters", action="store_true", help="Use filtered vcf for accuracy analysis (otherwise use unfiltered)")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")
parser.add_argument("--java", default="java", help="Path to java executable")
parser.add_argument("--picard_jar", default="picard.jar", help="Path to picard.jar executable")
parser.add_argument("--gatk", default="gatk", help="Path to gatk executable")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")
parser.add_argument("--cosmic_version", default=101, help="COSMIC version. Default: 101")


args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.out, "star_genome")
gatk_parent = args.out
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
genomes1000_vcf = args.genomes1000_vcf
threads = args.threads
read_length_minus_one = int(args.read_length) - 1
apply_mutation_filters = args.apply_mutation_filters
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_and_unmapped_bam = args.aligned_and_unmapped_bam

STAR = args.STAR
java = args.java
picard_jar = args.picard_jar
gatk = args.gatk

cosmic_tsv = args.cosmic_tsv
unique_mcrs_df_path = args.unique_mcrs_df_path
cosmic_version = args.cosmic_version

for name, path in {"STAR": STAR, "java": java, "picard_jar": picard_jar, "gatk": gatk}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

if not os.environ.get('JAVA_HOME'):
    os.environ['JAVA_HOME'] = os.path.dirname(os.path.dirname(java))

os.environ['PATH'] = f"{os.environ['JAVA_HOME']}/bin:" + os.environ['PATH']

os.makedirs(star_genome_dir, exist_ok=True)

alignment_folder = f"{gatk_parent}/alignment"
os.makedirs(alignment_folder, exist_ok=True)

gatk_supporting_files = f"{gatk_parent}/supporting_files"
os.makedirs(gatk_supporting_files, exist_ok=True)

plot_output_folder = f"{gatk_parent}/plots"
os.makedirs(plot_output_folder, exist_ok=True)

out_file_name_prefix = f"{alignment_folder}/sample_"

vcf_folder = f"{gatk_parent}/vcfs"
haplotypecaller_folder = f"{vcf_folder}/haplotypecaller"
haplotypecaller_folder = f"{vcf_folder}/haplotypecaller"

os.makedirs(vcf_folder, exist_ok=True)
os.makedirs(haplotypecaller_folder, exist_ok=True)
os.makedirs(haplotypecaller_folder, exist_ok=True)

aligned_only_bam = f"{alignment_folder}/aligned_only.bam"
unmapped_bam = f"{alignment_folder}/unmapped.bam"
merged_bam = f"{alignment_folder}/merged.bam"

marked_duplicates_bam = f"{alignment_folder}/marked_duplicates.bam"
marked_dup_metrics_txt = f"{alignment_folder}/marked_dup_metrics.txt"

split_n_cigar_reads_bam = f"{alignment_folder}/split_n_cigar_reads.bam"
recal_data_table = f"{alignment_folder}/recal_data.table"
recalibrated_bam = f"{alignment_folder}/recalibrated.bam"
covariates_plot = f"{alignment_folder}/AnalyzeCovariates.pdf"
haplotypecaller_unfiltered_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_unfiltered.g.vcf.gz"

haplotypecaller_filtered_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_filtered.vcf.gz"
haplotypecaller_filtered_applied_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_filtered_applied.vcf.gz"

panel_of_normals_vcf = f"{gatk_supporting_files}/1000g_pon.hg38.vcf.gz"
panel_of_normals_vcf_filtered = f"{gatk_supporting_files}/1000g_pon.hg38_filtered.vcf.gz"
haplotypecaller_unfiltered_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_unfiltered.g.vcf.gz"
haplotypecaller_filtered_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_filtered.vcf.gz"
haplotypecaller_filtered_applied_vcf = f"{haplotypecaller_folder}/haplotypecaller_output_filtered_applied.vcf.gz"

reference_genome_dict = reference_genome_fasta.replace(".fa", ".dict")

# commented out, as these should already be done prior to running this script
reference_genome_fasta_url = "https://ftp.ensembl.org/pub/grch37/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
reference_genome_gtf_url = "https://ftp.ensembl.org/pub/grch37/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
genomes1000_vcf_url = "https://ftp.ensembl.org/pub/grch37/release-93/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"

download_reference_genome_fasta_command = ["wget", "-O", f"{reference_genome_fasta}.gz", reference_genome_fasta_url]
unzip_reference_genome_fasta_command = ["gunzip", f"{reference_genome_fasta}.gz"]

download_reference_genome_gtf_command = ["wget", "-O", f"{reference_genome_gtf}.gz", reference_genome_gtf_url]
unzip_reference_genome_gtf_command = ["gunzip", f"{reference_genome_gtf}.gz"]

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
    "--outSAMunmapped", "Within",
    "--outSAMmapqUnique", "60",
    "--twopassMode", "Basic"
]
if not aligned_and_unmapped_bam:
    aligned_and_unmapped_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
if not os.path.exists(aligned_and_unmapped_bam):
    run_command_with_error_logging(star_align_command)

#* FASTQ to SAM
fastq_to_sam_command = [
    java, "-jar", picard_jar, "FastqToSam",
    "-FASTQ", synthetic_read_fastq,
    "-OUTPUT", unmapped_bam,
    "-READ_GROUP_NAME", "rg1",
    "-SAMPLE_NAME", "sample1",
    "-LIBRARY_NAME", "lib1",
    "-PLATFORM_UNIT", "unit1",
    "-PLATFORM", "ILLUMINA",
    "-SEQUENCING_CENTER", "center1"
]
if not os.path.exists(unmapped_bam):
    run_command_with_error_logging(fastq_to_sam_command)

#* CreateSequenceDictionary
create_sequence_dict_command = [
    java, "-jar", picard_jar, "CreateSequenceDictionary",
    "-R", reference_genome_fasta,
    "-O", reference_genome_dict
]
if not os.path.exists(reference_genome_dict):
    run_command_with_error_logging(create_sequence_dict_command)

#* MergeBamAlignment
merge_bam_alignment_command = [
    java, "-jar", picard_jar, "MergeBamAlignment",
    "--ALIGNED_BAM", aligned_and_unmapped_bam,
    "--UNMAPPED_BAM", unmapped_bam,
    "--OUTPUT", merged_bam,
    "--REFERENCE_SEQUENCE", reference_genome_fasta,
    "--SORT_ORDER", "coordinate",
    "--INCLUDE_SECONDARY_ALIGNMENTS", "false",
    "--VALIDATION_STRINGENCY", "SILENT"
]
if not os.path.exists(merged_bam):
    run_command_with_error_logging(merge_bam_alignment_command)

#* MarkDuplicates
mark_duplicates_command = [
    java, "-jar", picard_jar, "MarkDuplicates",
    "--INPUT", merged_bam,
    "--OUTPUT", marked_duplicates_bam,
    "--METRICS_FILE", marked_dup_metrics_txt,
    "--CREATE_INDEX", "true",
    "--VALIDATION_STRINGENCY", "SILENT"
]
if not os.path.exists(marked_duplicates_bam):
    run_command_with_error_logging(mark_duplicates_command)

#* SplitNCigarReads
split_n_cigar_reads_command = [
    gatk, "SplitNCigarReads",
    "-R", reference_genome_fasta,
    "-I", marked_duplicates_bam,
    "-O", split_n_cigar_reads_bam
]
if not os.path.exists(split_n_cigar_reads_bam):
    run_command_with_error_logging(split_n_cigar_reads_command)

#* IndexFeatureFile
index_feature_file_command = [
    gatk, "IndexFeatureFile",
    "-I", genomes1000_vcf
]
if not os.path.exists(f"{genomes1000_vcf}.idx"):
    run_command_with_error_logging(index_feature_file_command)

#* BaseRecalibrator
base_recalibrator_command = [
    gatk, "BaseRecalibrator",
    "-I", split_n_cigar_reads_bam,
    "-R", reference_genome_fasta,
    "--use-original-qualities",
    "--known-sites", genomes1000_vcf,
    "-O", recal_data_table
]
if not os.path.exists(recal_data_table):
    run_command_with_error_logging(base_recalibrator_command)

#* ApplyBQSR
apply_bqsr_command = [
    gatk, "ApplyBQSR",
    "--add-output-sam-program-record",
    "-R", reference_genome_fasta,
    "-I", split_n_cigar_reads_bam,
    "--use-original-qualities",
    "--bqsr-recal-file", recal_data_table,
    "-O", recalibrated_bam
]
if not os.path.exists(recalibrated_bam):
    run_command_with_error_logging(apply_bqsr_command)

#* AnalyzeCovariates
analyze_covariates_command = [
    gatk, "AnalyzeCovariates",
    "-bqsr", recal_data_table,
    "-plots", covariates_plot
]
if not os.path.exists(covariates_plot):
    run_command_with_error_logging(analyze_covariates_command)

#* HaplotypeCaller
haplotypecaller_command = [
    gatk,
    "--java-options", "-Xmx4g",
    "HaplotypeCaller",
    "-R", reference_genome_fasta,
    "-I", recalibrated_bam,
    "-O", haplotypecaller_unfiltered_vcf,
    "--dont-use-soft-clipped-bases",
    "--disable-tool-default-read-filters",
    "--standard-min-confidence-threshold-for-calling", "10"
]
if not os.path.exists(haplotypecaller_unfiltered_vcf):
    run_command_with_error_logging(haplotypecaller_command)

#* VariantFiltration
variantfiltration_command = [
    gatk,
    "VariantFiltration",
    "-R", reference_genome_fasta,
    "-V", haplotypecaller_unfiltered_vcf,
    "-O", haplotypecaller_filtered_vcf,
    "--window", "35",
    "--cluster", "3",
    "--filter-name", "FS",
    "--filter", "FS > 30.0",
    "--filter-name", "QD",
    "--filter", "QD < 2.0"
]
if not os.path.exists(haplotypecaller_filtered_vcf):
    run_command_with_error_logging(variantfiltration_command)

#* SelectVariants
selectvariants_command = [
    gatk,
    "SelectVariants",
    "-V", haplotypecaller_filtered_vcf,
    "--exclude-filtered", "true",
    "-O", haplotypecaller_filtered_applied_vcf
]
if not os.path.exists(haplotypecaller_filtered_applied_vcf):
    run_command_with_error_logging(selectvariants_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

cosmic_df_out = cosmic_tsv.replace(".tsv", "_vcf_info_for_fig2.csv")
if not os.path.exists(cosmic_df_out):
    cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out=cosmic_df_out, sequences="cdna", cosmic_version=cosmic_version)
else:
    cosmic_df = pd.read_csv(cosmic_df_out)

vcf_file = haplotypecaller_filtered_applied_vcf if apply_mutation_filters else haplotypecaller_unfiltered_vcf
perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=plot_output_folder, package_name="gatk_haplotypecaller")