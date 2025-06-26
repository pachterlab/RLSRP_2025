import os
import subprocess
import json
import pandas as pd
import anndata as ad
import concurrent.futures
import varseek as vk
from RLSRWP_2025.constants import box_links_dict, ccle_glioblastoma_rnaseq_experiment_aliases, ccle_glioblastoma_wxs_experiment_aliases, ccle_glioblastoma_wgs_experiment_aliases

# more on the datasets:
    # CCLE:
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJNA523380
        # paper: https://www.nature.com/articles/s41586-019-1186-3

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")

# Sequencing database parameters
sequencing_data_base = "ccle"  # ~11TB
data_to_use = "rnaseq"  # only used when sequencing_data_base == "ccle"  # options: "rnaseq", "wxs", "wxs_with_corresponding_rnaseq_sample", "rnaseq_and_wxs", "wxs_and_rnaseq"
technology = "bulk"
parity = "paired"
number_of_threads_total = 16  # if too high (e.g., 64), then will not be able to download successfully (server error) - 8 seems like the sweet spot
number_of_threads_per_varseek_count_task = 16
max_retries = 5
download_only = False
delete_fastq_files = False
overwrite_vk_count = False
sequencing_data_out_base = os.path.join(data_dir, f"{sequencing_data_base}_data_base")
experiment_aliases_to_keep = ccle_glioblastoma_rnaseq_experiment_aliases | ccle_glioblastoma_wxs_experiment_aliases | ccle_glioblastoma_wgs_experiment_aliases  # os.path.join(data_dir, "ccle_data_base", f"{sequencing_data_base}_metadata_filtered_experiment_aliases.txt")  # None to use all

# output
adata_combined_path = os.path.join(sequencing_data_out_base, "adata_vcrs_glioblastoma.h5ad")

# reference parameters
vk_ref_out = os.path.join(data_dir, "vk_ref_out")
vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
w=47
k=51

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
save_vcf = False

# for qc_against_gene_matrix - same as from vk ref/build (not essential but speeds things up)
variants = None if not qc_against_gene_matrix else os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")
seq_id_column = "seq_ID"
var_column = "mutation_cdna"
gene_id_column = "gene_name"
variants_usecols = [seq_id_column, var_column, gene_id_column]
add_hgvs_breakdown_to_adata_var = False
vcrs_metadata_df = os.path.join(reference_out_dir, "cosmic", "vcrs_metadata_df.csv")  # if it doesn't exist then it's fine

# for making VCF
vcf_data_csv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_vcf_data.csv")
cosmic_tsv=os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37.tsv")
cosmic_reference_genome_fasta=os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
sequences="cdna"

# summarize

# check for VCRS reference files
if not os.path.isdir(vk_ref_out) or len(os.listdir(vk_ref_out)) == 0:
    index_link, t2g_link = box_links_dict["cosmic_index"], box_links_dict["cosmic_t2g"]
    vk.utils.download_box_url(index_link, output_file_name=vcrs_index)
    vk.utils.download_box_url(t2g_link, output_file_name=vcrs_t2g)
    # vk.ref(variants="cosmic_cmc", sequences="cdna", w=w, k=k, out=vk_ref_out, download=True, index_out=vcrs_index, t2g_out=vcrs_t2g)
    # alternatively, to build from scratch: subprocess.run([os.path.join(script_dir, "run_vk_ref_cosmic.py")], check=True)

# check for kb count reference genome files when needed in vk count (i.e., when qc_against_gene_matrix=True)
if qc_against_gene_matrix and (not os.path.exists(reference_genome_index) or not os.path.exists(reference_genome_t2g)):
    if not os.path.exists(reference_genome_fasta) or not os.path.exists(reference_genome_gtf):
        reference_genome_out_dir = os.path.dirname(reference_genome_fasta)
        subprocess.run(["gget", "ref", "-w", "dna,gtf", "-r", "93", "--out_dir", reference_genome_out_dir, "-d", "human_grch37"], check=True)  # using grch37, ensembl 93 to agree with COSMIC
    reference_genome_f1 = os.path.join(reference_out_dir, "ensembl_grch37_release93", "f1.fasta")
    subprocess.run(["kb", "ref", "-t", str(number_of_threads_total), "-i", reference_genome_index, "-g", reference_genome_t2g, "-f1", reference_genome_f1, reference_genome_fasta, reference_genome_gtf], check=True)

#* CCLE/Geuvadis
ccle_ena_project = "PRJNA523380"
geuvadis_ena_project = "PRJEB3366"
if sequencing_data_base == "ccle":
    ena_project = ccle_ena_project
else:
    raise ValueError("ccle")

json_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ena_project}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"

if experiment_aliases_to_keep is not None:
    if isinstance(experiment_aliases_to_keep, str) and os.path.isfile(experiment_aliases_to_keep) and experiment_aliases_to_keep.endswith(".txt"):
        with open(experiment_aliases_to_keep, 'r', encoding="utf-8") as file:
            experiment_aliases_to_keep = {line.strip() for line in file.readlines()}
    elif isinstance(experiment_aliases_to_keep, (set, list)):
        experiment_aliases_to_keep = set(experiment_aliases_to_keep)
    else:
        raise ValueError("experiment_aliases_to_keep must be a file path or a set/list of experiment aliases")

    # replace "-" with "_"
    experiment_aliases_to_keep = {alias.replace("-", "_") for alias in experiment_aliases_to_keep}

# metadata json
os.makedirs(sequencing_data_out_base, exist_ok=True)
json_path = os.path.join(sequencing_data_out_base, f"{sequencing_data_base}_metadata.json")
if not os.path.exists(json_path):
    sequencing_metadata_download_command = ["wget", "-q", "-O", json_path, json_url]
    subprocess.run(sequencing_metadata_download_command, check=True)

    if sequencing_data_base == "ccle":
        try:
            json_updated_path_tmp = os.path.join(sequencing_data_out_base, "ccle_metadata_updated_tmp.json")
            update_ccle_metadata_script = os.path.join(script_dir, "update_ccle_metadata.R")
            subprocess.run(["Rscript", update_ccle_metadata_script, json_path, json_updated_path_tmp], check=True)  # TODO: will try to install dependencies in conda environment - unsure if this works
            os.rename(json_updated_path_tmp, json_path)
        except Exception as e:
            print("Error in updating CCLE metadata:", e)

# Loop through json file and download fastqs
with open(json_path, 'r', encoding="utf-8") as file:
    data = json.load(file)
json_df = pd.DataFrame(data)
json_df['experiment_alias_underscores_only'] = json_df['experiment_alias'].str.replace("-", "_")

if sequencing_data_base == "ccle":
    rnaseq_data = [study for study in data if study['library_strategy'] == 'RNA-Seq']
    wxs_data = [study for study in data if study['library_strategy'] == 'WXS']
    rnaseq_study_accessions = {study['study_accession'] for study in rnaseq_data}
    wxs_data_with_corresponding_rnaseq_sample = [study for study in wxs_data if study['sample_accession'] in rnaseq_study_accessions]

    if data_to_use.lower() == "rnaseq":
        data_list_to_run = rnaseq_data
    elif data_to_use.lower() == "wxs":
        data_list_to_run = wxs_data
    elif data_to_use.lower() == "wxs_with_corresponding_rnaseq_sample":
        data_list_to_run = wxs_data_with_corresponding_rnaseq_sample
    elif data_to_use.lower() == "rnaseq_and_wxs" or data_to_use.lower() == "wxs_and_rnaseq":
        data_list_to_run = rnaseq_data + wxs_data
    else:
        raise ValueError("data_to_use must be one of 'rnaseq', 'wxs', 'wxs_with_corresponding_rnaseq_sample', 'rnaseq_and_wxs', or 'wxs_and_rnaseq'")

number_of_items = len(data_list_to_run)

if download_only:
    number_of_threads_per_varseek_count_task = 1

if number_of_threads_total > 10:
    print("WARNING: diminishing returns after 10 threads for downloading")

if save_vcf and not os.path.exists(vcf_data_csv):  # alternatively, I can do this in vk clean by passing in vcf_data_csv=vcf_data_csv, cosmic_tsv=cosmic_tsv, cosmic_reference_genome_fasta=cosmic_reference_genome_fasta, variants="cosmic_cmc", sequences="cdna", cosmic_version=101
    vk.utils.add_vcf_info_to_cosmic_tsv(cosmic_tsv=cosmic_tsv, reference_genome_fasta=cosmic_reference_genome_fasta, cosmic_df_out=vcf_data_csv, sequences=sequences, cosmic_version=101)

def download_sequencing_total(
    record,
    vcrs_index,
    vcrs_t2g,
    technology="bulk",
    parity="paired",
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
    variants=None,
    seq_id_column=None,
    var_column=None,
    gene_id_column=None,
    variants_usecols=None,
    add_hgvs_breakdown_to_adata_var=None,
    vcrs_metadata_df=None,
    number_of_threads_per_varseek_count_task=2,
):
    experiment_alias = record.get('experiment_alias')
    experiment_alias_underscores_only = experiment_alias.replace("-", "_")
    sample = experiment_alias_underscores_only  # f"{experiment_alias_underscores_only}___{sample_accession}___{experiment_accession}___{run_accession}"

    if experiment_aliases_to_keep is not None and sample not in experiment_aliases_to_keep:
        # print(f"Skipping {sample} as it is not in the list of experiment aliases to keep.")
        return

    fastq_ftp = record.get('fastq_ftp')
    fastq_links = fastq_ftp.split(';')

    sample_out_folder = os.path.join(sequencing_data_out_base, sample)
    os.makedirs(sample_out_folder, exist_ok=True)

    failed_downloads = list()
    fastq_files = list()

    all_files_downloaded = True
    for link in fastq_links:
        rnaseq_fastq_file = os.path.join(sample_out_folder, os.path.basename(link))
        tmp_fastq_file = os.path.join(sample_out_folder, "tmp_" + os.path.basename(link))  # Temporary file, just so I know files that have not finished downloading
        fastq_files.extend([rnaseq_fastq_file, tmp_fastq_file])

        if not os.path.exists(rnaseq_fastq_file):  # just here while I keep files permanently for debugging
            all_files_downloaded = False
            
            if os.path.exists(tmp_fastq_file):  # If an old tmp file exists, delete it to ensure a clean restart
                print(f"Removing incomplete file: {tmp_fastq_file}")
                os.remove(tmp_fastq_file)
            
            print(f"Downloading {link} to {rnaseq_fastq_file}")

            if not link.startswith(('ftp://', 'http://')):
                link = 'ftp://' + link

            # download_command = f"curl --connect-timeout 60 --speed-time 30 --speed-limit 10000 -o {rnaseq_fastq_file} {link}"
            download_command = f"wget -c --tries={max_retries} --retry-connrefused -O {tmp_fastq_file} {link}"  # --limit-rate=1m  (or some other rate)
            # download_command = f"aria2c -x 16 -d {sample_out_folder} -o {os.path.basename(link)} -c {link}"

            try:
                result = subprocess.run(download_command, shell=True, check=True)
                os.rename(tmp_fastq_file, rnaseq_fastq_file)  # If successful, rename to final filename
            except subprocess.CalledProcessError as e:
                print(f"Error downloading {link} to {rnaseq_fastq_file}")
                print(e)
                failed_downloads.append(rnaseq_fastq_file)
                continue
    
    # count number of dirs in sequencing_data_out_base
    if not all_files_downloaded:
        if failed_downloads:
            print(f"Failed downloads in {sample}: {failed_downloads}")
        else:
            number_of_dirs = len(os.listdir(sequencing_data_out_base))
            print(f"Downloaded all files in {sample}, {number_of_dirs}/{number_of_items}")

    if download_only:
        return
    
    if os.path.isfile(vcrs_metadata_df):
        variants = None
    
    vk_count_out_dir = os.path.join(sample_out_folder, "vk_count_out")
    adata_cleaned_out = os.path.join(vk_count_out_dir, "adata_cleaned.h5ad")
    if not os.path.exists(adata_cleaned_out) or overwrite_vk_count:
        print(f"Running vk count on {sample}")
        vk_count_output_dict = vk.count(
            sample_out_folder,
            index=vcrs_index,
            t2g=vcrs_t2g,
            technology=technology,
            parity=parity,
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
            vcrs_metadata_df=vcrs_metadata_df,
            disable_clean=True,
        )

        print(f"Finished vk.count on {sample}")
    
    if quality_control_fastqs:
        fastqs_quality_controlled_dir = os.path.join(vk_count_out_dir, "fastqs_quality_controlled")
        if os.path.exists(fastqs_quality_controlled_dir):
            for file in os.listdir(fastqs_quality_controlled_dir):
                if file.endswith(".fastq.gz") or file.endswith(".fastq") or file.endswith(".fq.gz") or file.endswith(".fq"):
                    fastq_files.append(os.path.join(fastqs_quality_controlled_dir, file))

    if delete_fastq_files:
        for fastq_file in fastq_files:
            os.remove(fastq_file)


number_of_tasks = number_of_threads_total / number_of_threads_per_varseek_count_task
# with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks) as executor:
#     futures = [
#         executor.submit(
#             download_sequencing_total,
#             record=record,
#             vcrs_index=vcrs_index,
#             vcrs_t2g=vcrs_t2g,
#             technology=technology,
#             parity=parity,
#             sequencing_data_out_base=sequencing_data_out_base,
#             max_retries=max_retries,
#             k=k,
#             quality_control_fastqs=quality_control_fastqs,
#             cut_front=cut_front,
#             cut_tail=cut_tail,
#             reference_genome_index=reference_genome_index,
#             reference_genome_t2g=reference_genome_t2g,
#             qc_against_gene_matrix=qc_against_gene_matrix,
#             save_vcf=save_vcf,
#             vcf_data_csv=vcf_data_csv,
#             variants=variants,
#             seq_id_column=seq_id_column,
#             var_column=var_column,
#             gene_id_column=gene_id_column,
#             variants_usecols=variants_usecols,
#             add_hgvs_breakdown_to_adata_var=add_hgvs_breakdown_to_adata_var,
#             vcrs_metadata_df=vcrs_metadata_df,
#             number_of_threads_per_varseek_count_task=number_of_threads_per_varseek_count_task,
#         )
#         for record in data_list_to_run
#     ]

#     concurrent.futures.wait(futures)


for record in data_list_to_run:
    download_sequencing_total(
        record=record,
        vcrs_index=vcrs_index,
        vcrs_t2g=vcrs_t2g,
        technology=technology,
        parity=parity,
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
        variants=variants,
        seq_id_column=seq_id_column,
        var_column=var_column,
        gene_id_column=gene_id_column,
        variants_usecols=variants_usecols,
        add_hgvs_breakdown_to_adata_var=add_hgvs_breakdown_to_adata_var,
        vcrs_metadata_df=vcrs_metadata_df,
        number_of_threads_per_varseek_count_task=number_of_threads_per_varseek_count_task,
    )
        

if not os.path.exists(adata_combined_path):
    adata_list = []
    for sample in os.listdir(sequencing_data_out_base):
        if experiment_aliases_to_keep is not None and sample not in experiment_aliases_to_keep:
            continue
        if os.path.exists(os.path.join(sequencing_data_out_base, sample, "vk_count_out")):
            adata_single_path = os.path.join(sequencing_data_out_base, sample, "vk_count_out", "kb_count_out_vcrs", "counts_unfiltered", "adata.h5ad")
            adata_single = ad.read_h5ad(adata_single_path)
            adata_single.obs["experiment_alias_underscores_only"] = sample
            adata_list.append(adata_single)

    adata = ad.concat(adata_list, join='outer')

    adata.obs = adata.obs.merge(
        json_df[['experiment_accession', 'library_strategy', 'primary_disease', 
                'subtype_disease', 'sex', 'age', 'lineage_subtype', 
                'Cellosaurus_NCIt_disease', 'lineage', 'experiment_alias_underscores_only']],
        on='experiment_alias_underscores_only', 
        how='left'
    )
    adata.obs.index = adata.obs.index.astype(str)

    adata.write_h5ad(adata_combined_path)
print(f"Combined adata saved to {adata_combined_path}")