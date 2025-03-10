import argparse
import concurrent.futures
import json
import os
import subprocess


def check_for_successful_downloads(ccle_data_out_base, save_fastq_files = False):
    bad_samples = []

    for subfolder in os.listdir(ccle_data_out_base):
        subfolder_path = os.path.join(ccle_data_out_base, subfolder)
        
        if os.path.isdir(subfolder_path) and subfolder != "multiqc_total_data":
            print(f"Checking {subfolder}")
            sample_name = subfolder.split("___")[-1]

            fastqs_to_check = [f"{sample_name}_1", f"{sample_name}_2"]

            for fastq in fastqs_to_check:
                gzip_check_command = f"gzip -t {subfolder_path}/{fastq}.fastq.gz"
                try:
                    subprocess.run(gzip_check_command, shell=True, check=True)
                except subprocess.CalledProcessError:
                    bad_samples.append(subfolder)
                    break

    print(f"Samples with failed downloads: {bad_samples}")

    return bad_samples



def download_ccle_total(
    record,
    ccle_data_out_base = ".",
    max_retries = 5
):
    sample_accession = record.get('sample_accession')
    experiment_accession = record.get('experiment_accession')
    run_accession = record.get('run_accession')
    fastq_ftp = record.get('fastq_ftp')
    experiment_alias = record.get('experiment_alias')

    fastq_links = fastq_ftp.split(';')

    experiment_alias_underscores_only = experiment_alias.replace("-", "__")

    sample = f"{experiment_alias_underscores_only}___{sample_accession}___{experiment_accession}___{run_accession}"

    sample_out_folder = os.path.join(ccle_data_out_base, sample)

    failed_downloads = list()
    fastq_files = list()

    for link in fastq_links:
        rnaseq_fastq_file = os.path.join(sample_out_folder, os.path.basename(link))

        if not os.path.exists(rnaseq_fastq_file):  # just here while I keep files permanently for debugging
            print(f"Downloading {link} to {rnaseq_fastq_file}")

            if not link.startswith(('ftp://', 'http://')):
                link = 'ftp://' + link

            # download_command = f"curl --connect-timeout 60 --speed-time 30 --speed-limit 10000 -o {rnaseq_fastq_file} {link}"
            download_command = f"wget -c --tries={max_retries} --retry-connrefused -O {rnaseq_fastq_file} {link}"  # --limit-rate=1m  (or some other rate)
            # download_command = f"aria2c -x 16 -d {sample_out_folder} -o {os.path.basename(link)} -c {link}"

            try:
                result = subprocess.run(download_command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error downloading {link} to {rnaseq_fastq_file}")
                print(e)
                failed_downloads.append(rnaseq_fastq_file)
                continue
        
        fastq_files.append(rnaseq_fastq_file)
    
    print(f"Downloaded all files for run accession {run_accession}")
    print(f"Failed downloads: {failed_downloads}")
        

def main(args):
    # Paths and settings from arguments
    data_to_use = args.data_to_use
    json_path = args.json_path
    ccle_data_out_base = args.ccle_data_out_base
    number_of_tasks_total = args.number_of_tasks_total
    max_retries = args.max_retries

    if not os.path.exists(ccle_data_out_base):
        os.makedirs(ccle_data_out_base, exist_ok=True)

    if not json_path:
        json_path = "./ccle_metadata.json"

    if not os.path.exists(json_path):
        os.makedirs(os.path.dirname(json_path), exist_ok=True)
        ccle_metadata_download_command = f"wget -q -O {json_path} 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA523380&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0'"
        subprocess.run(ccle_metadata_download_command, shell=True, check=True)

        #!!! ad improvemenet R script

    # Loop through json file and download fastqs
    with open(json_path, 'r', encoding="utf-8") as file:
        data = json.load(file)

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

    with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks_total) as executor:
        futures = [
            executor.submit(
                download_ccle_total,
                record=record,
                ccle_data_out_base=ccle_data_out_base,
                max_retries=max_retries
            )
            for record in data_list_to_run
        ]

        concurrent.futures.wait(futures)

    # bad_samples = check_for_successful_downloads(ccle_data_out_base, save_fastq_files = save_fastq_files)

    with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_tasks_total) as executor:
        futures = [
            executor.submit(
                download_ccle_total,
                record=record,
                ccle_data_out_base=ccle_data_out_base,
                max_retries=max_retries
            )
            for record in data_list_to_run
        ]

        concurrent.futures.wait(futures)
    
    print("Program complete")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA-seq processing script with configurable parameters.")
    
    # Define arguments with argparse
    parser.add_argument("--data_to_use", type=str, default="rnaseq", help="Data to use (default: rnaseq) Options: 'rnaseq', 'wxs', 'wxs_with_corresponding_rnaseq_sample', 'rnaseq_and_wxs'.")
    parser.add_argument("--json_path", type=str, default="", help="Path to JSON file containing metadata.")
    parser.add_argument("--ccle_data_out_base", type=str, default=".", help="Base directory for kb output.")
    parser.add_argument("--number_of_tasks_total", type=int, default=4, help="Total number of concurrent tasks (default: 4).")
    parser.add_argument("--max_retries", type=int, default=20, help="Maximum number of retries for downloading files (default: 5).")

    # Parse arguments
    args = parser.parse_args()

    # Run main processing with parsed arguments
    main(args)


# swap --data_to_use rnaseq with wxs or wxs_with_corresponding_rnaseq_sample
#* see older version (before March 9 2025) for more args
# python3 download_and_align_ccle.py --mutation_index XXX --mutation_t2g XXX --standard_index XXX --standard_t2g XXX --json_path XXX --ccle_data_out_base XXX -k 59 --threads_per_task 1 --number_of_tasks_total 20 --data_to_use rnaseq