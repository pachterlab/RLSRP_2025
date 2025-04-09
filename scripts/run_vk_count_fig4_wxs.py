import os
import subprocess
import json
import concurrent.futures
import pandas as pd
import varseek as vk

# more on the datasets:
    # CCLE:
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJNA523380
        # paper: https://www.nature.com/articles/s41586-019-1186-3
    # Geuvadis: 
        # RNA-seq sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
        # WXS sequencing data (IGSR): https://www.internationalgenome.org/data-portal/sample
        # paper: https://www.nature.com/articles/nature12531

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")

# Sequencing database parameters
sequencing_data_base = "geuvadis"  # "ccle" (Fig3) or "geuvadis" (Fig4) - ccle is ~11TB, geuvadis is ~2.3TB
number_of_threads_total = 32  # if too high (e.g., 64), then will not be able to download successfully (server error) - 8 seems like the sweet spot
number_of_threads_per_varseek_count_task = 32
max_retries = 5
download_only = True
delete_fastq_files = False
overwrite_vk_count = False
sequencing_data_out_base = os.path.join(data_dir, f"{sequencing_data_base}_data_base")
geuvadis_rna_json_file = os.path.join(sequencing_data_out_base, f"{sequencing_data_base}_metadata.json")
metadata_tsv_file = os.path.join(sequencing_data_out_base, "igsr_GRCh38.tsv")  # can be found at https://www.internationalgenome.org/data-portal/data-collection/grch38

# reference parameters
vk_ref_out = os.path.join(data_dir, "vk_ref_out")
vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
w=47
k=51
dlist_reference_source = "t2t"

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
save_vcf = False

# for qc_against_gene_matrix - same as from vk ref/build (not essential but speeds things up)
variants = None if not qc_against_gene_matrix else os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")
seq_id_column = "seq_ID"
var_column = "mutation_cdna"
gene_id_column = "gene_name"
variants_usecols = [seq_id_column, var_column, gene_id_column]
add_hgvs_breakdown_to_adata_var=False

# for making VCF
vcf_data_csv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_vcf_data.csv")
cosmic_tsv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37.tsv")
cosmic_reference_genome_fasta=os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
sequences="cdna"

# summarize

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
    subprocess.run(["kb", "ref", "-t", str(number_of_threads_total), "-i", reference_genome_index, "-g", reference_genome_t2g, "-f1", reference_genome_f1, reference_genome_fasta, reference_genome_gtf], check=True)

if not os.path.exists(metadata_tsv_file):
    raise FileNotFoundError(f"Metadata TSV file not found: {metadata_tsv_file}. Please download it from https://www.internationalgenome.org/data-portal/data-collection/grch38 --> 'Available Data' --> 'Download the list' and place it in this path.")

if not os.path.exists(geuvadis_rna_json_file):
    json_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"
    sequencing_metadata_download_command = ["wget", "-q", "-O", geuvadis_rna_json_file, json_url]
    subprocess.run(sequencing_metadata_download_command, check=True)

# Loop through json file and download fastqs
with open(geuvadis_rna_json_file, 'r', encoding="utf-8") as file:
    data = json.load(file)

rna_df = pd.DataFrame(data)
sample_titles_set = set(rna_df['sample_title'].tolist())

wxs_df = pd.read_csv(metadata_tsv_file, sep="\t", low_memory=False, usecols=["Sample", "Data type", "Analysis group", "url"])
wxs_df = wxs_df[
    (wxs_df["Sample"].isin(sample_titles_set)) &
    (wxs_df["Data type"] == "sequence") &
    (wxs_df["Analysis group"] == "Exome")
]
wxs_df = wxs_df.drop_duplicates(subset=["Sample"])
number_of_items = len(wxs_df)

if download_only:
    number_of_threads_per_varseek_count_task = 1

if number_of_threads_total > 10:
    print("WARNING: diminishing returns after 10 threads for downloading")

if save_vcf and not os.path.exists(vcf_data_csv):  # alternatively, I can do this in vk clean by passing in vcf_data_csv=vcf_data_csv, cosmic_tsv=cosmic_tsv, cosmic_reference_genome_fasta=cosmic_reference_genome_fasta, variants="cosmic_cmc", sequences="cdna", cosmic_version=101
    vk.utils.add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=cosmic_reference_genome_fasta, cosmic_df_out=vcf_data_csv, sequences=sequences, cosmic_version=101)

def download_sequencing_total(
    row,
    vcrs_index,
    vcrs_t2g,
    sequencing_data_out_base = ".",
    max_retries = 5,
    k=59,
    quality_control_fastqs=False,
    cut_front=False,
    cut_tail=False,
    reference_genome_index=False,
    reference_genome_t2g=False,
    qc_against_gene_matrix=False,
    save_vcf=False,
    vcf_data_csv=None,
    number_of_threads_per_varseek_count_task=2,
):
    sample = row['Sample']
    link = row['url']

    sample_out_folder = os.path.join(sequencing_data_out_base, sample)
    os.makedirs(sample_out_folder, exist_ok=True)

    failed_downloads = list()
    fastq_files = list()

    all_files_downloaded = True
    wxs_fastq_file = os.path.join(sample_out_folder, os.path.basename(link))
    tmp_fastq_file = os.path.join(sample_out_folder, "tmp_" + os.path.basename(link))  # Temporary file, just so I know files that have not finished downloading
    fastq_files.extend([wxs_fastq_file, tmp_fastq_file])

    if not os.path.exists(wxs_fastq_file):  # just here while I keep files permanently for debugging
        all_files_downloaded = False
        
        if os.path.exists(tmp_fastq_file):  # If an old tmp file exists, delete it to ensure a clean restart
            print(f"Removing incomplete file: {tmp_fastq_file}")
            os.remove(tmp_fastq_file)
        
        print(f"Downloading {link} to {wxs_fastq_file}")

        if not link.startswith(('ftp://', 'http://')):
            link = 'ftp://' + link

        # download_command = f"curl --connect-timeout 60 --speed-time 30 --speed-limit 10000 -o {wxs_fastq_file} {link}"
        download_command = f"wget -c --tries={max_retries} --retry-connrefused -O {tmp_fastq_file} {link}"  # --limit-rate=1m  (or some other rate)
        # download_command = f"aria2c -x 16 -d {sample_out_folder} -o {os.path.basename(link)} -c {link}"

        try:
            result = subprocess.run(download_command, shell=True, check=True)
            os.rename(tmp_fastq_file, wxs_fastq_file)  # If successful, rename to final filename
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {link} to {wxs_fastq_file}")
            print(e)
            failed_downloads.append(wxs_fastq_file)
    
    # count number of dirs in sequencing_data_out_base
    if not all_files_downloaded:
        if failed_downloads:
            print(f"Failed downloads in {sample}: {failed_downloads}")
        else:
            number_of_dirs = len(os.listdir(sequencing_data_out_base))
            print(f"Downloaded all files in {sample}, {number_of_dirs}/{number_of_items}")

    if download_only:
        return
    
    vcrs_metadata_df = os.path.join(sample_out_folder, "vcrs_metadata_df.csv")
    if os.path.isfile(vcrs_metadata_df):
        variants = None

    vk_count_out_dir = os.path.join(sample_out_folder, "vk_count_out")
    adata_cleaned_out = os.path.join(vk_count_out_dir, "adata_cleaned.h5ad")
    if not os.path.exists(adata_cleaned_out) or overwrite_vk_count:
        print(f"Running vk.count on {sample}")
        vk_count_out_dir = os.path.join(sample_out_folder, "vk_count_out")
        vk_count_output_dict = vk.count(
            sample_out_folder,
            index=vcrs_index,
            t2g=vcrs_t2g,
            technology="bulk",
            k=k,
            quality_control_fastqs=quality_control_fastqs,
            cut_front=cut_front,
            cut_tail=cut_tail,
            reference_genome_index=reference_genome_index,
            reference_genome_t2g=reference_genome_t2g,
            qc_against_gene_matrix=qc_against_gene_matrix,
            out=vk_count_out_dir,
            threads=number_of_threads_per_varseek_count_task,
            save_vcf=save_vcf,
            vcf_data_csv=vcf_data_csv,
            variants=variants,
            seq_id_column=seq_id_column,
            var_column=var_column,
            gene_id_column=gene_id_column,
            variants_usecols=variants_usecols,
            add_hgvs_breakdown_to_adata_var=add_hgvs_breakdown_to_adata_var,
            vcrs_metadata_df=vcrs_metadata_df
        )

        print(f"Finished vk.count on {sample}")

    if delete_fastq_files:
        for fastq_file in fastq_files:
            if os.path.isfile(fastq_file):
                os.remove(fastq_file)


number_of_tasks = number_of_threads_total / number_of_threads_per_varseek_count_task
with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks) as executor:
    futures = [
        executor.submit(
            download_sequencing_total,
            row=row,
            vcrs_index=vcrs_index,
            vcrs_t2g=vcrs_t2g,
            sequencing_data_out_base=sequencing_data_out_base,
            max_retries=max_retries,
            k=k,
            quality_control_fastqs=quality_control_fastqs,
            cut_front=cut_front,
            cut_tail=cut_tail,
            reference_genome_index=reference_genome_index,
            reference_genome_t2g=reference_genome_t2g,
            qc_against_gene_matrix=qc_against_gene_matrix,
            save_vcf=save_vcf,
            vcf_data_csv=vcf_data_csv,
            number_of_threads_per_varseek_count_task=number_of_threads_per_varseek_count_task,
        )
        for _, row in wxs_df.iterrows()
    ]

    concurrent.futures.wait(futures)