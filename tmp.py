import os
import subprocess
import json
import concurrent.futures
import varseek as vk

# more on the datasets:
    # CCLE:
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJNA523380
        # paper: https://www.nature.com/articles/s41586-019-1186-3
    # Geuvadis: 
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
        # paper: https://www.nature.com/articles/nature12531

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")

# Sequencing database parameters
sequencing_data_base = "geuvadis"  # "ccle" (Fig3) or "geuvadis" (Fig4) - ccle is ~11TB, geuvadis is ~2.3TB
data_to_use = "rnaseq"  # only used when sequencing_data_base == "ccle"  # options: "rnaseq", "wxs", "wxs_with_corresponding_rnaseq_sample", "rnaseq_and_wxs", "wxs_and_rnaseq"
number_of_threads_total = 32  # if too high (e.g., 64), then will not be able to download successfully (server error) - 8 seems like the sweet spot
number_of_threads_per_varseek_count_task = 32
max_retries = 5
download_only = False
delete_fastq_files = False
sequencing_data_out_base = os.path.join(data_dir, f"{sequencing_data_base}_data_base")

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

#* CCLE/Geuvadis
ccle_ena_project = "PRJNA523380"
geuvadis_ena_project = "PRJEB3366"
if sequencing_data_base == "ccle":
    ena_project = ccle_ena_project
elif sequencing_data_base == "geuvadis":
    ena_project = geuvadis_ena_project
else:
    raise ValueError("ccle or geuvadis")

json_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ena_project}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"

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
elif sequencing_data_base == "geuvadis":
    data_list_to_run = data

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
    experiment_alias = record.get('experiment_alias')
    experiment_alias_underscores_only = experiment_alias.replace("-", "_")
    sample = experiment_alias_underscores_only  # f"{experiment_alias_underscores_only}___{sample_accession}___{experiment_accession}___{run_accession}"

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
        
        fastq_files.append(rnaseq_fastq_file)
    
    # count number of dirs in sequencing_data_out_base
    if not all_files_downloaded:
        if failed_downloads:
            print(f"Failed downloads in {sample}: {failed_downloads}")
        else:
            number_of_dirs = len(os.listdir(sequencing_data_out_base))
            print(f"Downloaded all files in {sample}, {number_of_dirs}/{number_of_items}")

    if download_only:
        return
    
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
    )

    print(f"Finished vk.count on {sample}")

    if delete_fastq_files:
        for fastq_file in fastq_files:
            os.remove(fastq_file)

    import sys  #!!! erase
    sys.exit()  #!!! erase


number_of_tasks = number_of_threads_total / number_of_threads_per_varseek_count_task
with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks) as executor:
    futures = [
        executor.submit(
            download_sequencing_total,
            record=record,
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
        for record in data_list_to_run
    ]

    concurrent.futures.wait(futures)