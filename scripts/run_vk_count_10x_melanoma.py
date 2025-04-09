import os
import subprocess
import varseek as vk

# more on the dataset: https://www.10xgenomics.com/datasets/melanoma-tumor-derived-cells-v-2-2-standard-4-0-0

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
vk_count_out_dir = os.path.join(data_dir, "vk_count_out_fig1")
kb_count_reference_genome_dir = os.path.join(data_dir, "vk_count_out_reference_genome_fig1")
reference_out_dir = os.path.join(data_dir, "reference")

# reference parameters
vk_ref_out = os.path.join(data_dir, "vk_ref_out")
vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
w=47
k=51
dlist_reference_source = "t2t"

# general vk count parameters
fastqs_dir = os.path.join(data_dir, "sc5p_v2_hs_melanoma_10k_fastqs")
technology = "10xv2"
threads = 2

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
qc_against_gene_matrix = True
account_for_strand_bias = True
strand_bias_end = "5p"
save_vcf = True

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

# check for fastqs
if not os.path.isdir(fastqs_dir) or len(os.listdir(fastqs_dir)) == 0:
    subprocess.run(["curl", "-P", data_dir, "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/sc5p_v2_hs_melanoma_10k/sc5p_v2_hs_melanoma_10k_fastqs.tar"], check=True)
    subprocess.run(["tar", "-xvf", f"{fastqs_dir}.tar", "-C", fastqs_dir], check=True)

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
    account_for_strand_bias=account_for_strand_bias,
    strand_bias_end=strand_bias_end,
    out=vk_count_out_dir,
    threads=threads,
    save_vcf=save_vcf,
    vcf_data_csv=vcf_data_csv,
    variants=variants,
    seq_id_column=seq_id_column,
    var_column=var_column,
    gene_id_column=gene_id_column,
    variants_usecols=variants_usecols
)