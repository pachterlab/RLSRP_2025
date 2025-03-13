import argparse
import os
import sys

import pysam
import pandas as pd
import shutil

from varseek.utils import (
    run_command_with_error_logging,
    add_vcf_info_to_cosmic_tsv,
    calculate_metrics,
    calculate_sensitivity_specificity,
    create_venn_diagram,
    draw_confusion_matrix,
    merge_gatk_and_cosmic,
    safe_literal_eval,
    vcf_to_dataframe,
)

parser = argparse.ArgumentParser(description="Run Strelka2 on a set of reads and report the time and memory usage")

# Paths
parser.add_argument("--synthetic_read_fastq", help="Path to synthetic read FASTQ")
parser.add_argument("--reference_genome_fasta", help="Path to reference genome fasta")
parser.add_argument("--reference_genome_gtf", help="Path to reference genome GTF")
parser.add_argument("--star_genome_dir", default="", help="Path to star_genome_dir")
parser.add_argument("--aligned_and_unmapped_bam", default="", help="Path to aligned_and_unmapped_bam. If not provided, will be created")
parser.add_argument("--tmp", default="tmp", help="Path to temp folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")
parser.add_argument("--STRELKA_INSTALL_PATH", default="STRELKA_INSTALL_PATH", help="Path to STRELKA_INSTALL_PATH parent")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--cosmic_cdna_info_csv", help="Path to COSMIC csv with cdna info from gget cosmic")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")


args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.tmp, "star_genome")
strelka2_output_dir = os.path.join(args.tmp, "strelka2_simulated_data_dir")
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = args.read_length - 1
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_and_unmapped_bam = args.aligned_and_unmapped_bam

STAR = args.STAR
STRELKA_INSTALL_PATH = args.STRELKA_INSTALL_PATH

cosmic_tsv = args.cosmic_tsv
cosmic_cdna_info_csv = args.cosmic_cdna_info_csv
unique_mcrs_df_path = args.unique_mcrs_df_path

for name, path in {"STAR": STAR, "STRELKA_INSTALL_PATH": STRELKA_INSTALL_PATH}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

os.makedirs(strelka2_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)

alignment_folder = f"{strelka2_output_dir}/alignment"
out_file_name_prefix = f"{alignment_folder}/sample_"

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

#* Strelka2 variant calling
strelka2_configure_command = [
    f"{STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py",
    "--bam", aligned_and_unmapped_bam,
    "--referenceFasta", reference_genome_fasta,
    "--rna",
    "--runDir", strelka2_output_dir
]

strelka2_run_command = [
    f"{strelka2_output_dir}/runWorkflow.py",
    "-m", "local",
    "-j", str(threads)
]


# commented out, as these should already be done prior to running this script
if not os.listdir(star_genome_dir):
    run_command_with_error_logging(star_build_command)

if not os.path.exists(f"{reference_genome_fasta}.fai"):
    _ = pysam.faidx(reference_genome_fasta)

if not os.path.exists(aligned_and_unmapped_bam):
    aligned_and_unmapped_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam"
    os.makedirs(alignment_folder, exist_ok=True)
    run_command_with_error_logging(star_align_command)

bam_index_file = f"{aligned_and_unmapped_bam}.bai"
if not os.path.exists(bam_index_file):
    _ = pysam.index(aligned_and_unmapped_bam)

run_command_with_error_logging(strelka2_configure_command)

run_command_with_error_logging(strelka2_run_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

#* Merging into COSMIC
# Convert VCF to DataFrame
strelka2_vcf = "output.vcf"  #!!!! change
df_strelka2 = vcf_to_dataframe(strelka2_vcf, additional_columns = True)
df_strelka2 = df_strelka2[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO_DP']].rename(columns={'INFO_DP': 'DP_strelka2'})

cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out = None, cosmic_cdna_info_csv = cosmic_cdna_info_csv, mutation_source = "cdna")

# load in unique_mcrs_df
unique_mcrs_df = pd.read_csv(unique_mcrs_df_path)
unique_mcrs_df["header_list"] = unique_mcrs_df["header_list"].apply(safe_literal_eval)

strelka2_cosmic_merged_df = merge_gatk_and_cosmic(df_strelka2, cosmic_df, exact_position=False)  # change exact_position to True to merge based on exact position as before
id_set_strelka2 = set(strelka2_cosmic_merged_df['ID'])

# Merge DP values into unique_mcrs_df
# Step 1: Remove rows with NaN values in 'ID' column
strelka2_cosmic_merged_df_for_merging = strelka2_cosmic_merged_df[['ID', 'DP_strelka2']].dropna(subset=['ID']).rename(columns={'ID': 'vcrs_header'})

# Step 2: Drop duplicates from 'ID' column
strelka2_cosmic_merged_df_for_merging = strelka2_cosmic_merged_df_for_merging.drop_duplicates(subset=['vcrs_header'])

# Step 3: Left merge with unique_mcrs_df
unique_mcrs_df = pd.merge(
    unique_mcrs_df,               # Left DataFrame
    strelka2_cosmic_merged_df_for_merging,         # Right DataFrame
    on='vcrs_header',
    how='left'
)

number_of_mutations_strelka2 = len(df_strelka2.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']))
number_of_cosmic_mutations_strelka2 = len(strelka2_cosmic_merged_df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']))

# unique_mcrs_df['header_list'] each contains a list of strings. I would like to make a new column unique_mcrs_df['mutation_detected_strelka2'] where each row is True if any value from the list unique_mcrs_df['vcrs_header'] is in the set id_set_mut  # keep in mind that my IDs are the mutation headers (ENST...), NOT mcrs headers or mcrs ids
unique_mcrs_df['mutation_detected_strelka2'] = unique_mcrs_df['header_list'].apply(
    lambda header_list: any(header in id_set_strelka2 for header in header_list)
)

# calculate expression error
unique_mcrs_df['mutation_expression_prediction_error'] = unique_mcrs_df['DP_strelka2'] - unique_mcrs_df['number_of_reads_mutant']  # positive means overpredicted, negative means underpredicted

unique_mcrs_df['TP'] = (unique_mcrs_df['included_in_synthetic_reads_mutant'] & unique_mcrs_df['mutation_detected_strelka2'])
unique_mcrs_df['FP'] = (~unique_mcrs_df['included_in_synthetic_reads_mutant'] & unique_mcrs_df['mutation_detected_strelka2'])
unique_mcrs_df['FN'] = (unique_mcrs_df['included_in_synthetic_reads_mutant'] & ~unique_mcrs_df['mutation_detected_strelka2'])
unique_mcrs_df['TN'] = (~unique_mcrs_df['included_in_synthetic_reads_mutant'] & ~unique_mcrs_df['mutation_detected_strelka2'])

strelka2_stat_path = f"{strelka2_output_dir}/reference_metrics_strelka2.txt"
metric_dictionary_reference = calculate_metrics(unique_mcrs_df, header_name = "vcrs_header", check_assertions = False, out = strelka2_stat_path)
draw_confusion_matrix(metric_dictionary_reference)

true_set = set(unique_mcrs_df.loc[unique_mcrs_df['included_in_synthetic_reads_mutant'], 'vcrs_header'])
positive_set = set(unique_mcrs_df.loc[unique_mcrs_df['mutation_detected_strelka2'], 'vcrs_header'])
create_venn_diagram(true_set, positive_set, TN = metric_dictionary_reference['TN'], out_path = f"{strelka2_output_dir}/venn_diagram_reference_cosmic_only_strelka2.png")

noncosmic_mutation_id_set = {f'strelka2_fp_{i}' for i in range(1, number_of_mutations_strelka2 - number_of_cosmic_mutations_strelka2 + 1)}

positive_set_including_noncosmic_mutations = positive_set.union(noncosmic_mutation_id_set)
false_positive_set = set(unique_mcrs_df.loc[unique_mcrs_df['FP'], 'vcrs_header'])
false_positive_set_including_noncosmic_mutations = false_positive_set.union(noncosmic_mutation_id_set)

FP_including_noncosmic = len(false_positive_set_including_noncosmic_mutations)
accuracy, sensitivity, specificity = calculate_sensitivity_specificity(metric_dictionary_reference['TP'], metric_dictionary_reference['TN'], FP_including_noncosmic, metric_dictionary_reference['FN'])

with open(strelka2_stat_path, "a", encoding="utf-8") as file:
    file.write(f"FP including non-cosmic: {FP_including_noncosmic}\n")
    file.write(f"accuracy including non-cosmic: {accuracy}\n")
    file.write(f"specificity including non-cosmic: {specificity}\n")



create_venn_diagram(true_set, positive_set_including_noncosmic_mutations, TN = metric_dictionary_reference['TN'], out_path = f"{strelka2_output_dir}/venn_diagram_reference_including_noncosmics_strelka2.png")

unique_mcrs_df.rename(columns={'TP': 'TP_strelka2', 'FP': 'FP_strelka2', 'TN': 'TN_strelka2', 'FN': 'FN_strelka2', 'mutation_expression_prediction_error': 'mutation_expression_prediction_error_strelka2'}, inplace=True)

unique_mcrs_df_out = unique_mcrs_df_path  # unique_mcrs_df_path.replace(".csv", "_with_strelka2.csv")
unique_mcrs_df.to_csv(unique_mcrs_df_out, index=False)