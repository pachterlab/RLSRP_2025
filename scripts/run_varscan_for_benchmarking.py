import argparse
import os
import subprocess
import sys

import pysam
import pandas as pd
import shutil

from varseek.utils import (
    run_command_with_error_logging,
    add_vcf_info_to_cosmic_tsv,
)

RLSRP_2025_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #!!! erase
sys.path.append(RLSRP_2025_dir)  #!!! erase
from RLSRP_2025.seq_utils import perform_analysis, compare_two_vcfs_with_hap_py

parser = argparse.ArgumentParser(description="Run VarScan on a set of reads and report the time and memory usage")

# Paths
parser.add_argument("--synthetic_read_fastq", help="Path to synthetic read FASTQ")
parser.add_argument("--reference_genome_fasta", help="Path to reference genome fasta")
parser.add_argument("--reference_genome_gtf", help="Path to reference genome GTF")
parser.add_argument("--star_genome_dir", default="", help="Path to star_genome_dir")
parser.add_argument("--aligned_bam", default="", help="Path to aligned_bam. If not provided, will be created")
parser.add_argument("--reference_genome_gtf_cleaned", default="", help="Path to reference_genome_gtf_cleaned. If not provided, it will not be used. If provided but does not exist, will be created")
parser.add_argument("--exons_bed", default="", help="Path to exons_bed. If not provided, it will not be used. If provided but does not exist, will be created")
parser.add_argument("--out", default="out", help="Path to out folder")

# Parameters
parser.add_argument("--threads", default=2, help="Number of threads")
parser.add_argument("--read_length", default=150, help="Read length")
parser.add_argument("--skip_accuracy_analysis", action="store_true", help="Skip accuracy analysis (beyond simple time and memory benchmarking)")

# Executables
parser.add_argument("--STAR", default="STAR", help="Path to STAR executable")
parser.add_argument("--VARSCAN_INSTALL_PATH", default="VARSCAN_INSTALL_PATH", help="Path to VARSCAN_INSTALL_PATH parent")

# Just for accuracy analysis
parser.add_argument("--cosmic_tsv", help="Path to COSMIC tsv")
parser.add_argument("--cosmic_vcf", help="Path to COSMIC vcf")
parser.add_argument("--unique_mcrs_df_path", help="Path to unique_mcrs_df_path from notebook 2")
parser.add_argument("--cosmic_version", default=101, help="COSMIC version. Default: 101")

args = parser.parse_args()

star_genome_dir = args.star_genome_dir if args.star_genome_dir else os.path.join(args.out, "star_genome")
varscan_output_dir = args.out
reference_genome_fasta = args.reference_genome_fasta
reference_genome_gtf = args.reference_genome_gtf
threads = args.threads
read_length_minus_one = int(args.read_length) - 1
skip_accuracy_analysis = args.skip_accuracy_analysis
synthetic_read_fastq = args.synthetic_read_fastq
aligned_bam = args.aligned_bam
reference_genome_gtf_cleaned = args.reference_genome_gtf_cleaned
exons_bed = args.exons_bed

STAR = args.STAR
VARSCAN_INSTALL_PATH = args.VARSCAN_INSTALL_PATH

cosmic_tsv = args.cosmic_tsv
unique_mcrs_df_path = args.unique_mcrs_df_path
cosmic_version = args.cosmic_version

for name, path in {"STAR": STAR, "VARSCAN_INSTALL_PATH": VARSCAN_INSTALL_PATH}.items():
    if not os.path.exists(path) and not shutil.which(path):
        raise FileNotFoundError(f"{name} not found or installed properly.")

os.makedirs(varscan_output_dir, exist_ok=True)
os.makedirs(star_genome_dir, exist_ok=True)


alignment_folder = f"{varscan_output_dir}/alignment"
out_file_name_prefix = f"{alignment_folder}/sample_"
aligned_bam = f"{out_file_name_prefix}Aligned.sortedByCoord.out.bam" if not aligned_bam else aligned_bam
data_pileup_file = f"{varscan_output_dir}/simulated_data.pileup"
vcf_file = f"{varscan_output_dir}/varscan_variants.vcf"  # variants.vcf
plot_output_folder = f"{varscan_output_dir}/plots"
os.makedirs(plot_output_folder, exist_ok=True)


# commented out, as these should already be done prior to running this script
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
if len(os.listdir(star_genome_dir)) == 0:
    run_command_with_error_logging(star_build_command)

#* Reference genome index file
if not os.path.exists(f"{reference_genome_fasta}.fai"):
    _ = pysam.faidx(reference_genome_fasta)
# commented out, as these should already be done prior to running this script

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

#* Exons BED file creation (iff exons_bed provided and does not exist)
if exons_bed and not os.path.isfile(exons_bed):
    if not reference_genome_gtf_cleaned:
        reference_genome_gtf_cleaned = "gtf_cleaned.gtf"
    if not os.path.isfile(reference_genome_gtf_cleaned):
        gffread_command = ["gffread", "-E", reference_genome_gtf, "-T", "-o", reference_genome_gtf_cleaned]  # gffread -E data/reference/ensembl_grch37_release93/Homo_sapiens.GRCh37.87.gtf -T -o gtf_cleaned.gtf
        try:
            run_command_with_error_logging(gffread_command)
        except Exception as e:
            reference_genome_gtf_cleaned = reference_genome_gtf

    exon_bed_command = f"awk '$3 == \"exon\" {{print $1, $4-1, $5}}' OFS='\t' {reference_genome_gtf_cleaned} > {exons_bed}"  # awk '$3 == "exon" {print $1, $4-1, $5}' OFS='\t' gtf_cleaned.gtf > exons.bed
    run_command_with_error_logging(exon_bed_command)

#* Samtools mpileup
samtools_mpileup_command = f"samtools mpileup -B -f {reference_genome_fasta} {aligned_bam} > {data_pileup_file}"
if exons_bed:
    samtools_mpileup_command = samtools_mpileup_command.replace(aligned_bam, f"-l {exons_bed} {aligned_bam}")
if not os.path.isfile(data_pileup_file):
    print("Running samtools mpileup")
    run_command_with_error_logging(samtools_mpileup_command)

#* Varscan variant calling
varscan_command = f"java -jar {VARSCAN_INSTALL_PATH} mpileup2cns {data_pileup_file} --output-vcf 1 --variants 1 --min-coverage 1 --min-reads2 1 --min-var-freq 0.01 --strand-filter 0 > {vcf_file}"
if not os.path.isfile(vcf_file):
    print("Running varscan")
    run_command_with_error_logging(varscan_command)

if skip_accuracy_analysis:
    print("Skipping accuracy analysis")
    sys.exit()

# cosmic_df_out = cosmic_tsv.replace(".tsv", "_vcf_info_for_fig2.csv")
# if not os.path.exists(cosmic_df_out):
#     cosmic_df = add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=reference_genome_fasta, cosmic_df_out=cosmic_df_out, sequences="cdna", cosmic_version=cosmic_version)
# else:
#     cosmic_df = pd.read_csv(cosmic_df_out)
# perform_analysis(vcf_file=vcf_file, unique_mcrs_df_path=unique_mcrs_df_path, cosmic_df=cosmic_df, plot_output_folder=plot_output_folder, package_name = "varscan", dp_column="INFO_ADP")


contigs_txt = os.path.join(varscan_output_dir, "contigs.txt")
header_txt = os.path.join(varscan_output_dir, "header.txt")
new_header_txt = os.path.join(varscan_output_dir, "new_header.txt")
output_vcf = os.path.join(varscan_output_dir, "varscan_variants_with_contigs.vcf")

if not os.path.exists(contigs_txt):
    with open(contigs_txt, "w") as f:
        subprocess.run(
            ["awk", "{print \"##contig=<ID=\"$1\",length=\"$2\">\"}", f"{reference_genome_fasta}.fai"],
            stdout=f,
            check=True
        )

if not os.path.exists(header_txt):
    with open(header_txt, "w") as f:
        subprocess.run(["bcftools", "view", "-h", vcf_file], stdout=f, check=True)

if not os.path.exists(new_header_txt):
    with open(header_txt) as f:
        lines = f.readlines()
    pre_chrom = [line for line in lines if not line.startswith("#CHROM")]
    chrom_line = [line for line in lines if line.startswith("#CHROM")]

    # Write new header with contigs before #CHROM
    with open(new_header_txt, "w") as out_f, open(contigs_txt) as f2:
        out_f.writelines(pre_chrom)
        out_f.writelines(f2.readlines())
        out_f.writelines(chrom_line)

if not os.path.exists(output_vcf):
    # Reheader the VCF file with the new contigs
    subprocess.run(["bcftools", "reheader", "-h", new_header_txt, "-o", output_vcf, vcf_file], check=True)

vcf_file = output_vcf  # Update vcf_file to the new VCF with contigs

package_name = "varscan"
cosmic_vcf = args.cosmic_vcf
happy_out = os.path.join(args.out, "hap_py_out", package_name)

cosmic_vcf = os.path.realpath(cosmic_vcf)
vcf_file = os.path.realpath(vcf_file)
reference_genome_fasta = os.path.realpath(reference_genome_fasta)
happy_out = os.path.realpath(happy_out)
compare_two_vcfs_with_hap_py(ground_truth_vcf=cosmic_vcf, test_vcf=vcf_file, reference_fasta=reference_genome_fasta, output_dir = happy_out, unique_mcrs_df = unique_mcrs_df_path, unique_mcrs_df_out = None, package_name = package_name, dry_run = False)