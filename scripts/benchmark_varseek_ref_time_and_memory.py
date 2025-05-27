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
number_of_variants_list = [1, 4, 16, 64, 256, 1024, 4096]  # number of variants, in thousands
cosmic_mutations_path = os.path.join(reference_out_dir, "cosmic", "CancerMutationCensus_AllData_Tsv_v101_GRCh37", "CancerMutationCensus_AllData_v101_GRCh37_mutation_workflow.csv")  # for vk sim
sequences_fasta_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.cds.all.fa")
gtf_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.87.gtf")
threads = 16
w = 37
k = 41

# only if cosmic_mutations_path does not exist
cosmic_version = "101"
grch = "37"
os.environ['COSMIC_EMAIL'] = 'your_email'  # replace with your email, or store in an environment variable COSMIC_EMAIL; only needed if cosmic_mutations_path does not exist
os.environ['COSMIC_PASSWORD'] = 'your_password'  # replace with your password, or store in an environment variable COSMIC_PASSWORD; only needed if cosmic_mutations_path does not exist

store_reference_files_in_tmp = True  # if True, will store reference files in tmp directory; if False, will store reference files in data directory
### ARGUMENTS ###

varseek_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(data_dir, "Fig3")  #* change for each run
tmp_dir = "/data/benchmarking_vk_ref_tmp"  #!! replace with "tmp"

random_seed = 42
vk_ref_script_path = os.path.join(script_dir, "run_varseek_ref_for_benchmarking.py")
vk_count_script_path = os.path.join(script_dir, "run_varseek_count_for_benchmarking.py")

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

number_of_reads =  256_000_000
fastq_output_path = os.path.join(tmp_dir, f"reads_{number_of_reads}_fastq.fastq")
if not os.path.exists(fastq_output_path):   
    import varseek as vk
    reference_cdna_path = os.path.join(reference_out_dir, "ensembl_grch37_release93", "Homo_sapiens.GRCh37.cdna.all.fa")  # for vk sim
    number_of_reads_per_variant_alt=100
    number_of_reads_per_variant_ref=150
    number_of_reads_per_variant_total = number_of_reads_per_variant_alt + number_of_reads_per_variant_ref
    number_of_variants_to_sample = number_of_reads // number_of_reads_per_variant_total

    logger.info(f"Building synthetic reads for {number_of_reads} reads")
    _ = vk.sim(
        variants=cosmic_mutations_path,
        reads_fastq_out=fastq_output_path,
        number_of_variants_to_sample=number_of_variants_to_sample,
        strand=None,
        number_of_reads_per_variant_alt=number_of_reads_per_variant_alt,
        number_of_reads_per_variant_ref=number_of_reads_per_variant_ref,
        read_length=150,
        seed=random_seed,
        add_noise_sequencing_error=True,
        add_noise_base_quality=False,
        error_rate=0.0001,
        error_distribution=(0.85, 0.1, 0.05),
        max_errors=float("inf"),
        with_replacement=True,
        gzip_reads_fastq_out=False,
        sequences=reference_cdna_path,
        seq_id_column="seq_ID",
        var_column="mutation_cdna",
        variant_type_column=None,
        reference_out_dir=reference_out_dir,
        out=os.path.join(tmp_dir, "vk_sim_out"),
        k=k,
        w=w,
        make_dataframes=False
    )

logger.info("Loading in COSMIC mutations")
cosmic_df = pd.read_csv(cosmic_mutations_path, usecols=["seq_ID", "mutation"])

os.makedirs(tmp_dir, exist_ok=True)
output_file_vk_ref = os.path.join(output_dir, "time_and_memory_benchmarking_report_vk_ref.txt")
output_file_vk_count = os.path.join(output_dir, "time_and_memory_benchmarking_report_vk_count.txt")

for number_of_variants in number_of_variants_list:
    number_of_variants *= 1000
    logger.info(f"Subsampling {number_of_variants} variants")
    cosmic_df_subsampled_path = os.path.join(tmp_dir, f"cosmic_subsampled_{number_of_variants}_variants.csv")
    if not os.path.exists(cosmic_df_subsampled_path):
        cosmic_df_subsampled = cosmic_df.sample(n=number_of_variants, random_state=random_seed)
        cosmic_df_subsampled.to_csv(cosmic_df_subsampled_path, index=False)
    
    logger.info("varseek ref")
    tmp_dir_specific_run = os.path.join(tmp_dir, f"vk_ref_{number_of_variants}_variants")
    os.makedirs(tmp_dir_specific_run, exist_ok=True)
    script_title = f"vk ref {number_of_variants} variants {threads} threads"
    argparse_flags = f"--variants {cosmic_df_subsampled_path} --sequences {sequences_fasta_path} --out {tmp_dir_specific_run} --seq_id_column seq_ID --var_column mutation --reference_out_dir {reference_out_dir} --w {w} --k {k} --dlist_reference_source t2t --threads {threads}"  # make sure to provide all keyword args with two-dashes (ie full argument name)
    print(f"python3 {vk_ref_script_path} {argparse_flags}")
    _ = report_time_and_memory_of_script(vk_ref_script_path, output_file = output_file_vk_ref, argparse_flags = argparse_flags, script_title = script_title)

    vk_ref_index_path = os.path.join(tmp_dir_specific_run, "vcrs_index.idx")
    vk_ref_t2g_path = os.path.join(tmp_dir_specific_run, "vcrs_t2g_filtered.txt")

    logger.info("varseek count")
    tmp_dir_specific_run = os.path.join(tmp_dir, f"vk_count_{number_of_variants}_variants")
    os.makedirs(tmp_dir_specific_run, exist_ok=True)
    script_title = f"vk count {number_of_variants} variants {number_of_reads} reads {threads} threads"
    argparse_flags = f"--index {vk_ref_index_path} --t2g {vk_ref_t2g_path} --technology bulk --threads {threads} --k {k} --out {tmp_dir_specific_run} --disable_summarize --fastqs {fastq_output_path}"  # --disable_fastqpp
    print(f"python3 {vk_count_script_path} {argparse_flags}")
    _ = report_time_and_memory_of_script(vk_count_script_path, output_file = output_file_vk_count, argparse_flags = argparse_flags, script_title = script_title)

