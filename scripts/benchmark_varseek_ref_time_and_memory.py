# import varseek as vk
import os
import logging

import pandas as pd

from varseek.utils import report_time_and_memory_of_script

logger = logging.getLogger(__name__)
logger.setLevel("INFO")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S")
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")

### ARGUMENTS ###
number_of_variants_list = [1, 4]  # [1, 4, 16, 64, 256, 1024, 4096]  # number of variants, in thousands
cosmic_mutations_path = os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")  # for vk sim
sequences_fasta_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.cds.all.fa")
gtf_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.87.gtf")
threads = 16

# only if cosmic_mutations_path does not exist
cosmic_version = "101"
grch = "37"
os.environ['COSMIC_EMAIL'] = 'your_email'  # replace with your email, or store in an environment variable COSMIC_EMAIL; only needed if cosmic_mutations_path does not exist
os.environ['COSMIC_PASSWORD'] = 'your_password'  # replace with your password, or store in an environment variable COSMIC_PASSWORD; only needed if cosmic_mutations_path does not exist

store_reference_files_in_tmp = True  # if True, will store reference files in tmp directory; if False, will store reference files in data directory
### ARGUMENTS ###

varseek_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(data_dir, "time_and_memory_benchmarking_out_dir_vk_ref")  #* change for each run
tmp_dir = "/data/benchmarking_vk_ref_tmp"  #!! replace with "tmp"

random_seed = 42
vk_ref_script_path = os.path.join(script_dir, "run_varseek_ref_for_benchmarking.py")

# download reference genome
if not os.path.isfile(sequences_fasta_path):
    import subprocess
    sequences_fasta_path = os.path.join(reference_out_dir, "Homo_sapiens.GRCh37.cds.all.fa")
    os.makedirs(reference_out_dir, exist_ok=True)
    
    if grch == "37":
        gget_ref_species = "human_grch37"
    else:
        gget_ref_species = "human"

    ref_sequence_download_command = f"gget ref -w cds -r 93 --out_dir {reference_out_dir} -d {gget_ref_species}"
    sequences_download_command_list = ref_sequence_download_command.split(" ")

    logger.info(f"Downloading reference sequences with {' '.join(sequences_download_command_list)}. Note that this requires curl >=7.73.0")
    subprocess.run(sequences_download_command_list, check=True)
    subprocess.run(["gunzip", f"{sequences_fasta_path}.gz"], check=True)

# download COSMIC
if not os.path.isfile(cosmic_mutations_path):
    import gget
    reference_out_dir_cosmic = os.path.dirname(os.path.dirname(cosmic_mutations_path))
    os.makedirs(os.path.dirname(cosmic_mutations_path), exist_ok=True)
    
    logger.info(f"Downloading COSMIC mutations with gget")
    gget.cosmic(
        None,
        grch_version=grch,
        cosmic_version=cosmic_version,
        out=reference_out_dir,
        cosmic_project="cancer",
        download_cosmic=True,
        gget_mutate=True,
        keep_genome_info=True,
        remove_duplicates=True,
        email=os.environ.get('COSMIC_EMAIL'),
        password=os.environ.get('COSMIC_PASSWORD'),
    )

logger.info("Loading in COSMIC mutations")
cosmic_df = pd.read_csv(cosmic_mutations_path, usecols=["seq_ID", "mutation"])  #!!! uncomment

os.makedirs(tmp_dir, exist_ok=True)

for number_of_variants in number_of_variants_list:
    number_of_variants *= 1000
    logger.info(f"Subsampling {number_of_variants} variants")
    tmp_dir_specific_run = os.path.join(tmp_dir, f"vk_ref_{number_of_variants}_variants")
    cosmic_df_subsampled_path = os.path.join(tmp_dir_specific_run, f"cosmic_subsampled_{number_of_variants}_variants.csv")
    if not os.path.exists(cosmic_df_subsampled_path):
        cosmic_df_subsampled = cosmic_df.sample(n=number_of_variants, random_state=random_seed)
        cosmic_df_subsampled.to_csv(cosmic_df_subsampled_path, index=False)
    
    output_file = os.path.join(output_dir, f"vk_ref_threads_{threads}_variants_{number_of_variants}_time_and_memory.txt")
    argparse_flags = f"--variants {cosmic_df_subsampled_path} --sequences {sequences_fasta_path} --out {tmp_dir_specific_run} --seq_id_column seq_ID --var_column mutation --reference_out_dir {reference_out_dir} --dlist_reference_source t2t --threads {threads}"  # make sure to provide all keyword args with two-dashes (ie full argument name)
    report_time_and_memory_of_script(vk_ref_script_path, output_file = output_file, argparse_flags = argparse_flags)

