import os
import subprocess
import sys
import varseek as vk

# more on the dataset: 
#   PRJNA330719 (289.54GB for .sra files, 2.3T for .fastq files): https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_67eecb767ec86104c28549e7&o=acc_s%3Aa  # find the GEO here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465
#   PRJNA603103, PRJNA603104: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA603101&f=librarysource_s%3An%3Atranscriptomic%3Borganism_s%3An%3Ahomo%2520sapiens%3Bsource_name_sam_ss%3An%3Askin%3Ac%3Binstrument_s%3An&o=instrument_s%3Ad  # find the GEO here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144240
project_name_to_sra_link_dict = {
    "PRJNA330719": "https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_67eecb767ec86104c28549e7&o=acc_s%3Aa",
    "PRJNA603103": "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA603101&f=librarysource_s%3An%3Atranscriptomic%3Borganism_s%3An%3Ahomo%2520sapiens%3Bsource_name_sam_ss%3An%3Askin%3Bbioproject_s%3An%3Aprjna603103%3Ac&o=instrument_s%3Ad%3Bacc_s%3Aa",
    "PRJNA603104": "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA603101&f=librarysource_s%3An%3Atranscriptomic%3Borganism_s%3An%3Ahomo%2520sapiens%3Bsource_name_sam_ss%3An%3Askin%3Bbioproject_s%3An%3Aprjna603104%3Ac&o=instrument_s%3Ad%3Bacc_s%3Aa"
}

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
vk_count_out_dir = os.path.join(data_dir, "vk_count_out_fig1")
kb_count_reference_genome_dir = os.path.join(vk_count_out_dir, "kb_count_out_reference_genome")
reference_out_dir = os.path.join(data_dir, "reference")

project_name = "PRJNA330719"  # "PRJNA330719", "PRJNA603103", "PRJNA603104"
SRR_Acc_List_path = os.path.join(data_dir, "SRA", project_name, "SRR_Acc_List.txt")  # find the raw data here: 
download_only = False
overwrite_vk_count = False

# reference parameters
vk_ref_out = os.path.join(data_dir, "vk_ref_out")
vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
w=47
k=51
dlist_reference_source = "t2t"

# general vk count parameters
fastqs_dir = os.path.join(data_dir, "glioblastoma_smartseq_fastq_data")  # for PRJNA330719; change for each run
technology = "BULK"  # bulk is used for non-multiplexed SMARTSEQ2
parity = "paired"
threads = 16

# fastqpp
quality_control_fastqs = True
cut_front = True
cut_tail = True

# kb count, reference genome
reference_genome_index = os.path.join(reference_out_dir, "ensembl_grch37_release93", "index.idx")  # can either already exist or will be created; only used if qc_against_gene_matrix=True
reference_genome_t2g = os.path.join(reference_out_dir, "ensembl_grch37_release93", "t2g.txt")  # can either already exist or will be created; only used if qc_against_gene_matrix=True
reference_genome_fasta = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")  # can either already exist or will be downloaded; only used if qc_against_gene_matrix=True
reference_genome_gtf = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.87.gtf")  # can either already exist or will be downloaded; only used if qc_against_gene_matrix=True

# clean
qc_against_gene_matrix = False
save_vcf = True
vcf_out = os.path.join(vk_count_out_dir, "variants.vcf")
min_counts = 2
filter_cells_by_min_counts = False  # can be done later
filter_cells_by_min_genes = False  # can be done later
cpm_normalization = False  # False because I just want counts directly

# for qc_against_gene_matrix - same as from vk ref/build (not essential but speeds things up)
variants = None if not qc_against_gene_matrix else os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")
seq_id_column = "seq_ID"
var_column = "mutation_cdna"
gene_id_column = "gene_name"
variants_usecols = None  # all columns

# for making VCF
vcf_data_csv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_vcf_data.csv")
cosmic_tsv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37.tsv")
cosmic_reference_genome_fasta=os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
sequences="cdna"

# summarize

# maf
maf_out = os.path.join(vk_count_out_dir, "variants.maf")
vcf2maf_pl = "vcf2maf.pl"
genome_reference_fasta = os.path.join(data_dir, "reference", "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
vep_data = os.path.expanduser("~/.vep")
conda_path = os.path.expanduser("~/miniconda3")
conda_environment = "vep93"

# check for fastqs
if not os.path.isdir(fastqs_dir) or len(os.listdir(fastqs_dir)) == 0:
    if not os.path.isfile(SRR_Acc_List_path):
        raise ValueError(f"Must have {SRR_Acc_List_path} downloaded to run this script. Download the accession numbers list file from this link: {project_name_to_sra_link_dict[project_name]}")
    with open(SRR_Acc_List_path) as f:
        srr_list = f.read().split()
    try:
        print("Downloading fastq files (289.54GB for .sra files, 2.3T for .fastq files)")
        print("Running prefetch")
        prefetch_command = ["prefetch", "--progress", "--verbose"] + srr_list
        subprocess.run(prefetch_command, check=True)
        print("Downloading files")
        data_download_command = ["fasterq-dump", "--outdir", fastqs_dir, "--threads", str(threads), "--progress", "--verbose", "--split-files"] + srr_list
        subprocess.run(data_download_command, check=True)
        # print(" ".join(data_download_command))
        print(f"Fastq data downloaded successfully to {fastqs_dir}")
    except Exception as e:
        print(f"Error running fasterq-dump: {e}")
        raise  # re-raises the original exception
else:
    print(f"Fastq data already exists in {fastqs_dir}. If you want to re-download, please delete the directory first.")

if download_only:
    print("download_only=True. Exiting script.")
    sys.exit()

# check for VCRS reference files
if not os.path.isdir(vk_ref_out) or len(os.listdir(vk_ref_out)) == 0:
    vk.ref(variants="cosmic_cmc", sequences="cdna", w=w, k=k, out=vk_ref_out, dlist_reference_source=dlist_reference_source, download=True, index_out=vcrs_index, t2g_out=vcrs_t2g)
    # alternatively, to build from scratch: subprocess.run([os.path.join(script_dir, "run_vk_ref.py")], check=True)

# check for kb count reference genome files when needed in vk count (i.e., when qc_against_gene_matrix=True)
if qc_against_gene_matrix and (not os.path.exists(reference_genome_index) or not os.path.exists(reference_genome_t2g)):
    if not os.path.exists(reference_genome_fasta) or not os.path.exists(reference_genome_gtf):
        reference_genome_out_dir = os.path.dirname(reference_genome_fasta)
        subprocess.run(["gget", "ref", "-w", "dna,gtf", "-r", "93", "--out_dir", reference_genome_out_dir, "-d", "human_grch37"], check=True)  # using grch37, ensembl 93 to agree with COSMIC
    reference_genome_f1 = os.path.join(reference_out_dir, "ensembl_grch37_release93", "f1.fasta")
    subprocess.run(["kb", "ref", "-t", str(threads), "-i", reference_genome_index, "-g", reference_genome_t2g, "-f1", reference_genome_f1, reference_genome_fasta, reference_genome_gtf], check=True)

if save_vcf and not os.path.exists(vcf_data_csv):  # alternatively, I can do this in vk clean by passing in vcf_data_csv=vcf_data_csv, cosmic_tsv=cosmic_tsv, cosmic_reference_genome_fasta=cosmic_reference_genome_fasta, variants="cosmic_cmc", sequences="cdna", cosmic_version=101
    vk.utils.add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=cosmic_reference_genome_fasta, cosmic_df_out=vcf_data_csv, sequences=sequences, cosmic_version=101)

adata_cleaned_out = os.path.join(vk_count_out_dir, "adata_cleaned.h5ad")
if not os.path.exists(adata_cleaned_out) or overwrite_vk_count:
    vk_count_output_dict = vk.count(
        fastqs_dir,
        index=vcrs_index,
        t2g=vcrs_t2g,
        technology=technology,
        k=k,
        quality_control_fastqs=quality_control_fastqs,
        cut_front=cut_front,
        cut_tail=cut_tail,
        reference_genome_index=reference_genome_index,
        reference_genome_t2g=reference_genome_t2g,
        kb_count_reference_genome_out_dir=kb_count_reference_genome_dir,
        qc_against_gene_matrix=qc_against_gene_matrix,
        parity=parity,
        out=vk_count_out_dir,
        threads=threads,
        save_vcf=save_vcf,
        vcf_out=vcf_out,
        vcf_data_csv=vcf_data_csv,
        variants=variants,
        seq_id_column=seq_id_column,
        var_column=var_column,
        gene_id_column=gene_id_column,
        variants_usecols=variants_usecols,
        min_counts=min_counts,
        filter_cells_by_min_counts=filter_cells_by_min_counts,
        filter_cells_by_min_genes=filter_cells_by_min_genes,
        cpm_normalization=cpm_normalization,
        disable_summarize=True
    )

if save_vcf and maf_out is not None:
    if not os.path.exists(vcf2maf_pl) or not os.path.exists(f"{conda_path}/envs/{conda_environment}"):
        print(f"Please set up the vcf2maf command with scripts/vcf2maf_setup.sh")
        sys.exit()
    vcf2maf_command = ["conda", "run", "-n", conda_environment, "perl", vcf2maf_pl, "--input-vcf", vcf_out, "--output-maf", maf_out, "--ref-fasta", genome_reference_fasta, "--vep-path", f"{conda_path}/envs/{conda_environment}/bin", "--vep-data", vep_data, "--species", "homo_sapiens", "--ncbi-build", "GRCh37", "--cache-version", "93", "--retain-info", "AO,NS"]
    subprocess.run(vcf2maf_command, check=True)