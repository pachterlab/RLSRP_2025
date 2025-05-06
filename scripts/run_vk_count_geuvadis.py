import os
import subprocess
import json
import pandas as pd
import re
import anndata as ad
import concurrent.futures
import varseek as vk
from RLSRWP_2025.constants import box_links_dict

# more on the datasets:
    # Geuvadis: 
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
        # paper: https://www.nature.com/articles/nature12531

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")

# Sequencing database parameters
sequencing_data_base = "geuvadis"  # ~2.3TB
technology = "bulk"
parity = "paired"
number_of_threads_total = 32  # if too high (e.g., 64), then will not be able to download successfully (server error) - 8 seems like the sweet spot
number_of_threads_per_varseek_count_task = 32
max_retries = 5
download_only = False
delete_fastq_files = False
overwrite_vk_count = False
sequencing_data_out_base = os.path.join(data_dir, f"{sequencing_data_base}_data_base")
experiment_aliases_to_keep = {"E_GEUV_1:HG00377.1.M_120209_6"}  # {"E_GEUV_1:HG00377.1.M_120209_6"}  # None to use all

# reference parameters
vk_ref_out_parent = os.path.join(data_dir, "vk_ref_out_geuvadis")

sequences = os.path.join(reference_out_dir, "ensembl_grch37_release113", "Homo_sapiens.GRCh37.cdna.all.fa.gz")
w_and_k_list_of_dicts = [
    {"w": 27, "k": 31},
    {"w": 37, "k": 41},
    {"w": 47, "k": 51},
]

geuvadis_reference_files_dir = os.path.join(reference_out_dir, "geuvadis")
variants = os.path.join(geuvadis_reference_files_dir, "variants_transcriptome.parquet")
geuvadis_reference_variants_prefix = os.path.join(geuvadis_reference_files_dir, "1kg_phase1_all")
geuvadis_true_vcf = os.path.join(geuvadis_reference_files_dir, f"{geuvadis_reference_variants_prefix}.vcf.gz")
sample_metadata_tsv_file = os.path.join(sequencing_data_out_base, "sample_metadata.tsv")

# fastqpp
quality_control_fastqs = True
cut_front = True
cut_tail = True

# kb count, reference genome
reference_genome_index = os.path.join(reference_out_dir, "ensembl_grch37_release113", "index.idx")  # can either already exist or will be created; only used if qc_against_gene_matrix=True
reference_genome_t2g = os.path.join(reference_out_dir, "ensembl_grch37_release113", "t2g.txt")  # can either already exist or will be created; only used if qc_against_gene_matrix=True
reference_genome_fasta = os.path.join(reference_out_dir, "ensembl_grch37_release113", "Homo_sapiens.GRCh37.dna.primary_assembly.fa")  # can either already exist or will be downloaded; only used if qc_against_gene_matrix=True
reference_genome_gtf = os.path.join(reference_out_dir, "ensembl_grch37_release113", "Homo_sapiens.GRCh37.87.gtf")  # can either already exist or will be downloaded; only used if qc_against_gene_matrix=True

# clean
qc_against_gene_matrix = True
qc_against_gene_matrix_mistake_ratio = 0.5  # None for none
save_vcf = False

# for qc_against_gene_matrix - same as from vk ref/build (not essential but speeds things up)
seq_id_column = "transcript_ID"
var_column = "variant_cdna"
gene_id_column = "gene_name"
variants_usecols = [seq_id_column, var_column, gene_id_column]
add_hgvs_breakdown_to_adata_var = False

# summarize

if not os.path.exists(sample_metadata_tsv_file):
    raise FileNotFoundError(f"Metadata TSV file not found: {sample_metadata_tsv_file}. Please download it from https://www.internationalgenome.org/data-portal/sample --> 'Filter by data collection' --> 'Geuvadis' --> 'Download the list' and place it in this path.")

# check for VCRS reference files
# alternatively, to build from scratch: subprocess.run([os.path.join(script_dir, "run_vk_ref_geuvadis.py")], check=True)    
for w_and_k_dict in w_and_k_list_of_dicts:
    w, k = w_and_k_dict["w"], w_and_k_dict["k"]
    vk_ref_out = os.path.join(vk_ref_out_parent, f"w{w}_k{k}")
    vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
    vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")

    if not os.path.isdir(vk_ref_out) or len(os.listdir(vk_ref_out)) == 0:
        index_link, t2g_link = box_links_dict[f"geuvadis_index_w{w}_k{k}"], box_links_dict[f"geuvadis_t2g_w{w}_k{k}"]
        vk.utils.download_box_url(index_link, output_file_name=vcrs_index)
        vk.utils.download_box_url(t2g_link, output_file_name=vcrs_t2g)
        # vk.ref(variants="geuvadis", sequences="cdna", w=w, k=k, out=vk_ref_out, index_out=vcrs_index, t2g_out=vcrs_t2g, download=True)

# check for kb count reference genome files
if not os.path.exists(reference_genome_index) or not os.path.exists(reference_genome_t2g):
    if not os.path.exists(reference_genome_fasta) or not os.path.exists(reference_genome_gtf):
        reference_genome_out_dir = os.path.dirname(reference_genome_fasta)
        subprocess.run(["gget", "ref", "-w", "dna,gtf", "-r", "113", "--out_dir", reference_genome_out_dir, "-d", "human_grch37"], check=True)
    reference_genome_f1 = os.path.join(reference_out_dir, "ensembl_grch37_release93", "f1.fasta")
    subprocess.run(["kb", "ref", "-t", str(number_of_threads_total), "-i", reference_genome_index, "-g", reference_genome_t2g, "-f1", reference_genome_f1, reference_genome_fasta, reference_genome_gtf], check=True)

#* Geuvadis
geuvadis_ena_project = "PRJEB3366"
if sequencing_data_base == "geuvadis":
    ena_project = geuvadis_ena_project
else:
    raise ValueError("geuvadis")

json_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ena_project}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"

if experiment_aliases_to_keep is not None:
    if isinstance(experiment_aliases_to_keep, str) and os.path.isfile(experiment_aliases_to_keep):
        with open(experiment_aliases_to_keep, 'r', encoding="utf-8") as file:
            experiment_aliases_to_keep = {line.strip() for line in file.readlines()}

    experiment_aliases_to_keep = {re.sub(r"[-.:]", "_", experiment_alias) for experiment_alias in experiment_aliases_to_keep}  # replace dash, period, and colon with underscore

# metadata json
os.makedirs(sequencing_data_out_base, exist_ok=True)
json_path = os.path.join(sequencing_data_out_base, f"{sequencing_data_base}_metadata.json")
if not os.path.exists(json_path):
    sequencing_metadata_download_command = ["wget", "-q", "-O", json_path, json_url]
    subprocess.run(sequencing_metadata_download_command, check=True)

# Loop through json file and download fastqs
with open(json_path, 'r', encoding="utf-8") as file:
    data = json.load(file)
json_df = pd.DataFrame(data)
json_df['experiment_alias_underscores_only'] = json_df['experiment_alias'].str.replace(r"[-.:]", "_", regex=True)  # replace dash, period, and colon with underscore

sample_metadata_df = pd.read_csv(sample_metadata_tsv_file, sep="\t")
json_df = json_df.merge(sample_metadata_df, left_on='sample_title', right_on='Sample name', how='left')

data_list_to_run = data

number_of_items = len(data_list_to_run)

if download_only:
    number_of_threads_per_varseek_count_task = 1

if number_of_threads_total > 10:
    print("WARNING: diminishing returns after 10 threads for downloading")

def download_sequencing_total(
    record,
    vk_ref_out_parent,
    technology="bulk",
    parity="paired",
    sequencing_data_out_base = ".",
    max_retries = 5,
    w_and_k_dict=None,
    quality_control_fastqs=False,
    cut_front=False,
    cut_tail=False,
    reference_genome_index=False,
    reference_genome_t2g=False,
    qc_against_gene_matrix=False,
    qc_against_gene_matrix_mistake_ratio=None,
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
    experiment_alias_underscores_only = re.sub(r"[-.:]", "_", experiment_alias)
    sample = experiment_alias_underscores_only  # f"{experiment_alias_underscores_only}___{sample_accession}___{experiment_accession}___{run_accession}"

    if sample in {"E_GEUV_1_NA20527_1_M_111124_6"}:
        return  # bad sample - has 3 FASTQs

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
    
    if isinstance(vcrs_metadata_df, str) and os.path.isfile(vcrs_metadata_df):
        variants = None
    
    if quality_control_fastqs:
        quality_control_fastqs_out_dir = os.path.join(sample_out_folder, "fastqs_quality_controlled")
        vk.fastqpp(
            sample_out_folder,
            technology=technology,
            parity=parity,
            quality_control_fastqs=quality_control_fastqs,
            cut_front=cut_front,
            cut_tail=cut_tail,
            length=31,
            quality_control_fastqs_out_dir=quality_control_fastqs_out_dir,
            threads=number_of_threads_per_varseek_count_task,
        )
        fastq_dir = quality_control_fastqs_out_dir
        for file in os.listdir(fastq_dir):
            if file.endswith(".fastq.gz") or file.endswith(".fastq") or file.endswith(".fq.gz") or file.endswith(".fq"):
                fastq_files.append(os.path.join(fastq_dir, file))
    else:
        fastq_dir = sample_out_folder
    
    for w_and_k_dict in w_and_k_list_of_dicts:
        w, k = w_and_k_dict["w"], w_and_k_dict["k"]
        vk_ref_out = os.path.join(vk_ref_out_parent, f"w{w}_k{k}")
        vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
        vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
        
        vk_count_out_dir = os.path.join(sample_out_folder, f"vk_count_out_w{w}_k{k}")
        adata_cleaned_out = os.path.join(vk_count_out_dir, "adata_cleaned.h5ad")
        if not os.path.exists(adata_cleaned_out) or overwrite_vk_count:
            print(f"Running vk.count on {sample}")
            vk_count_output_dict = vk.count(
                fastq_dir,
                index=vcrs_index,
                t2g=vcrs_t2g,
                technology=technology,
                parity=parity,
                k=k,
                reference_genome_index=reference_genome_index,
                reference_genome_t2g=reference_genome_t2g,
                qc_against_gene_matrix=qc_against_gene_matrix,
                mistake_ratio=qc_against_gene_matrix_mistake_ratio,
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
                disable_fastqpp=True,
                # disable_clean=True,
                disable_summarize=True,
            )

            print(f"Finished vk.count on {sample}")
        
    kb_count_out_reference_genome = os.path.join(sample_out_folder, "kb_count_out_reference_genome")
    if not os.path.exists(kb_count_out_reference_genome) or len(os.listdir(kb_count_out_reference_genome)) == 0 or overwrite_vk_count:
        os.makedirs(kb_count_out_reference_genome, exist_ok=True)
        fastqs_final = [os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir) if (f.endswith(".fastq") or f.endswith(".fastq.gz") or f.endswith(".fq") or f.endswith(".fq.gz"))]
        kb_count_standard_index_command = [
            "kb",
            "count",
            "-t",
            str(number_of_threads_per_varseek_count_task),
            "-i",
            reference_genome_index,
            "-g",
            reference_genome_t2g,
            "-x",
            technology,
            "--h5ad",
            "--parity",
            parity,
            "-o",
            kb_count_out_reference_genome,
            "--overwrite",
        ] + fastqs_final
        print(f"Running kb count on {sample} against reference genome")
        subprocess.run(kb_count_standard_index_command, check=True)

    if delete_fastq_files:
        for fastq_file in fastq_files:
            os.remove(fastq_file)


number_of_tasks = number_of_threads_total / number_of_threads_per_varseek_count_task
with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks) as executor:
    futures = [
        executor.submit(
            download_sequencing_total,
            record=record,
            vk_ref_out_parent=vk_ref_out_parent,
            technology=technology,
            parity=parity,
            sequencing_data_out_base=sequencing_data_out_base,
            max_retries=max_retries,
            w_and_k_dict=w_and_k_dict,
            quality_control_fastqs=quality_control_fastqs,
            cut_front=cut_front,
            cut_tail=cut_tail,
            reference_genome_index=reference_genome_index,
            reference_genome_t2g=reference_genome_t2g,
            qc_against_gene_matrix=qc_against_gene_matrix,
            qc_against_gene_matrix_mistake_ratio=qc_against_gene_matrix_mistake_ratio,
            save_vcf=save_vcf,
            variants=variants,
            seq_id_column=seq_id_column,
            var_column=var_column,
            gene_id_column=gene_id_column,
            variants_usecols=variants_usecols,
            add_hgvs_breakdown_to_adata_var=add_hgvs_breakdown_to_adata_var,
            number_of_threads_per_varseek_count_task=number_of_threads_per_varseek_count_task,
        )
        for record in data_list_to_run
    ]

    concurrent.futures.wait(futures)

for w_and_k_dict in w_and_k_list_of_dicts:
    w, k = w_and_k_dict["w"], w_and_k_dict["k"]
    adata_combined_path_vcrs = os.path.join(sequencing_data_out_base, f"adata_vcrs_combined_w{w}_k{k}.h5ad")
    if not os.path.exists(adata_combined_path_vcrs):
        adata_vcrs_list = []
        for sample in os.listdir(sequencing_data_out_base):
            adata_vcrs_single_path = os.path.join(sequencing_data_out_base, sample, f"vk_count_out_w{w}_k{k}", "kb_count_out_vcrs", "counts_unfiltered", "adata.h5ad")
            if os.path.exists(adata_vcrs_single_path):
                adata_vcrs_single = ad.read_h5ad(adata_vcrs_single_path)
                adata_vcrs_single.obs["experiment_alias_underscores_only"] = sample
                adata_vcrs_list.append(adata_vcrs_single)

        adata_vcrs = ad.concat(adata_vcrs_list, join='outer')

        adata_vcrs.obs = adata_vcrs.obs.merge(
            json_df[['experiment_accession', 'library_strategy', 'sample_title', 
                    'Sex', 'Biosample ID', 'Population name', 'Superpopulation name', 
                    'experiment_alias_underscores_only']],
            on='experiment_alias_underscores_only', 
            how='left'
        )
        adata_vcrs.obs.index = adata_vcrs.obs.index.astype(str)

        adata_vcrs.write_h5ad(adata_combined_path_vcrs)
    print(f"Combined adata saved to {adata_combined_path_vcrs}")



adata_combined_path_reference_genome = os.path.join(sequencing_data_out_base, "adata_reference_genome_combined.h5ad")
if not os.path.exists(adata_combined_path_reference_genome):
    adata_reference_genome_list = []
    for sample in os.listdir(sequencing_data_out_base):
        adata_reference_genome_single_path = os.path.join(sequencing_data_out_base, sample, "kb_count_out_reference_genome", "counts_unfiltered", "adata.h5ad")
        if os.path.exists(adata_reference_genome_single_path):
            adata_reference_genome_single = ad.read_h5ad(adata_reference_genome_single_path)
            adata_reference_genome_single.obs["experiment_alias_underscores_only"] = sample
            adata_reference_genome_list.append(adata_reference_genome_single)

    adata_reference_genome = ad.concat(adata_reference_genome_list, join='outer')

    adata_reference_genome.obs = adata_reference_genome.obs.merge(
        json_df[['experiment_accession', 'library_strategy', 'sample_title', 
                'Sex', 'Biosample ID', 'Population name', 'Superpopulation name', 
                'experiment_alias_underscores_only']],
        on='experiment_alias_underscores_only', 
        how='left'
    )
    adata_reference_genome.obs.index = adata_reference_genome.obs.index.astype(str)

    adata_reference_genome.write_h5ad(adata_combined_path_reference_genome)
print(f"Combined adata saved to {adata_combined_path_reference_genome}")