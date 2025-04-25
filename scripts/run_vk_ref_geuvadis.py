import os
import subprocess
import numpy as np
import requests
import anndata as ad
import pickle
from scipy.sparse import lil_matrix, csr_matrix
from tqdm import tqdm
import pandas as pd
tqdm.pandas()
import json
import varseek as vk
from varseek.utils import is_program_installed, vcf_to_dataframe, assign_transcript_and_cds, add_variant_type, reverse_complement, write_to_vcf, add_variant_type_column_to_vcf_derived_df, add_variant_column_to_vcf_derived_df, convert_mutation_cds_locations_to_cdna
from varseek.constants import mutation_pattern
from RLSRWP_2025.seq_utils import make_transcript_df_from_gtf, convert_vcf_samples_to_anndata, keep_only_exons_in_vcf, chunks, add_to_hgvsc_dict, make_hgvsc_to_dbsnp_dict_from_vep_vcf

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


from pdb import set_trace as st

# more on the datasets:
    # Geuvadis: 
        # sequencing data (ENA): https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
        # paper: https://www.nature.com/articles/nature12531

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
reference_out_dir = os.path.join(data_dir, "reference")
geuvadis_reference_files_dir = os.path.join(reference_out_dir, "geuvadis")
ensembl_reference_files_dir = os.path.join(reference_out_dir, "ensembl_grch37_release113")

#* reference parameters
# # inputs
variants = os.path.join(geuvadis_reference_files_dir, "variants_transcriptome.parquet")
sequences = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.cdna.all.fa")
reference_genome_gtf = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.87.gtf")

# # parameters
w_and_k_list_of_dicts = [
    {"w": 27, "k": 31},
    {"w": 37, "k": 41},
    {"w": 47, "k": 51},
]
seq_id_column = "transcript_ID"
var_column = "variant_cdna"
dlist_reference_source = "t2t"
threads = 16
chunksize = None

# # output paths
vk_ref_out_parent = os.path.join(data_dir, "vk_ref_out_geuvadis")

# # misc
geuvadis_reference_variants_prefix = os.path.join(geuvadis_reference_files_dir, "1kg_phase1_all")
geuvadis_preliminary_vcf = f"{geuvadis_reference_variants_prefix}_preliminary.vcf.gz"
geuvadis_preliminary_vcf_exons_only = f"{geuvadis_reference_variants_prefix}_preliminary_exons.vcf.gz"
geuvadis_preliminary_vcf_exons_only_df_path = f"{geuvadis_reference_variants_prefix}_preliminary_vcf_exons.parquet"
geuvadis_true_vcf = f"{geuvadis_reference_variants_prefix}_true.vcf.gz"
geuvadis_genotype_preliminary_adata = os.path.join(geuvadis_reference_files_dir, "genotypes_adata_preliminary.h5ad")
geuvadis_genotype_true_adata = os.path.join(geuvadis_reference_files_dir, "genotypes_adata_true.h5ad")
gtf_df_path = os.path.join(geuvadis_reference_files_dir, "gtf_df.parquet")
transcript_df_path = os.path.join(geuvadis_reference_files_dir, "transcript_df.parquet")
save_vcf_samples = False
plink = "plink"
conda_path = os.path.expanduser("~/miniconda3")
conda_environment = "vep113"
vep_cache_dir = os.path.expanduser("~/.vep")
geuvadis_preliminary_vep_vcf = f"{geuvadis_reference_variants_prefix}_preliminary.vep.vcf.gz"
genome_fasta = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
cds_path = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.cds.all.fa")
hgvsc_to_dbsnp_pkl_path = os.path.join(geuvadis_reference_files_dir, "hgvsc_to_dbsnp_dict_tmp.pkl")

geuvadis_rna_json_file = os.path.join(geuvadis_reference_files_dir, "geuvadis_metadata.json")
total_entries = 517114  # zcat /home/jrich/Desktop/RLSRWP_2025/data/reference/geuvadis/1kg_phase1_all_preliminary_exons.vcf.gz | wc -l

#* Download geuvadis RNA metadata (for extracting sample names)
if not os.path.exists(geuvadis_rna_json_file):
    logger.info("Downloading Geuvadis RNA metadata")
    json_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"
    sequencing_metadata_download_command = ["wget", "-q", "-O", geuvadis_rna_json_file, json_url]
    if os.path.dirname(geuvadis_rna_json_file) and not os.path.exists(os.path.dirname(geuvadis_rna_json_file)):
        os.makedirs(os.path.dirname(geuvadis_rna_json_file), exist_ok=True)
    subprocess.run(sequencing_metadata_download_command, check=True)

#* download cDNA
os.makedirs(ensembl_reference_files_dir, exist_ok=True)
if not os.path.exists(sequences):
    logger.info("Downloading cDNA sequences")
    gget_ref_command = ["gget", "ref", "-w", "cdna", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    unzip_command = ["gunzip", f"{sequences}.gz"]
    subprocess.run(unzip_command, check=True)

#* download CDS
if not os.path.exists(cds_path):
    logger.info("Downloading CDS sequences")
    gget_ref_command = ["gget", "ref", "-w", "cds", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    unzip_command = ["gunzip", f"{sequences}.gz"]
    subprocess.run(unzip_command, check=True)

#* download gtf
if not os.path.exists(reference_genome_gtf):
    logger.info("Downloading GTF file")
    gget_ref_command = ["gget", "ref", "-w", "gtf", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    unzip_command = ["gunzip", f"{reference_genome_gtf}.gz"]
    subprocess.run(unzip_command, check=True)

#* download Geuvadis reference variants
if not os.path.exists(f"{geuvadis_reference_variants_prefix}.bim") or not os.path.exists(f"{geuvadis_reference_variants_prefix}.bed") or not os.path.exists(f"{geuvadis_reference_variants_prefix}.fam"):
    logger.info("Downloading Geuvadis reference variants")
    plink_file_download_command = ["wget", "-O", f"{geuvadis_reference_variants_prefix}.tar.gz", "https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz"]
    subprocess.run(plink_file_download_command, check=True)
    untar_command = ["tar", "-xzvf", f"{geuvadis_reference_variants_prefix}.tar.gz", "-C", geuvadis_reference_files_dir]
    subprocess.run(untar_command, check=True)

#* Convert plink variant files to preliminary VCF
if not os.path.exists(geuvadis_preliminary_vcf):  # plink pseudo-VCF doesn't exist
    logger.info("Converting plink files to preliminary VCF")
    if not is_program_installed(plink):
        raise ValueError(f"plink is not installed. Please install plink to run this script.")
    plink_to_preliminary_vcf_command = [plink, "--bfile", "/data/geuvadis_data_base/genome_variants_ground_truth/1kg_phase1_all", "--threads", str(threads), "--recode", "vcf", "bgz", "--out", geuvadis_reference_variants_prefix]
    subprocess.run(plink_to_preliminary_vcf_command, check=True)

if not os.path.exists(geuvadis_preliminary_vcf_exons_only):
    keep_only_exons_in_vcf(vcf=geuvadis_preliminary_vcf, gtf=reference_genome_gtf, vcf_out=geuvadis_preliminary_vcf_exons_only, total=total_entries)

#* Make VCF index file
if not os.path.exists(f"{geuvadis_preliminary_vcf_exons_only}.tbi"):
    logger.info("Creating VCF index")
    vcf_index_creation_command = ["tabix", "-p", "vcf", geuvadis_preliminary_vcf_exons_only]
    subprocess.run(vcf_index_creation_command, check=True)

if total_entries is None:
    logger.info("Counting total entries in VCF")
    total_entries = subprocess.run(f'zcat {geuvadis_preliminary_vcf_exons_only} | grep -v "^#" | wc -l', capture_output=True, text=True, shell=True).stdout.strip()
    logger.info(f"Total entries in VCF: {total_entries}")

#* make adata
if not os.path.exists(geuvadis_genotype_preliminary_adata):
    logger.info("Converting preliminary VCF to adata")
    with open(geuvadis_rna_json_file, 'r', encoding="utf-8") as file:
        data = json.load(file)

    rna_df = pd.DataFrame(data)
    sample_titles_set = set(rna_df['sample_title'].tolist())

    adata = convert_vcf_samples_to_anndata(geuvadis_preliminary_vcf_exons_only, sample_titles_set=sample_titles_set, adata_out=geuvadis_genotype_preliminary_adata, total=total_entries)
else:
    logger.info("Loading preliminary adata")
    adata = ad.read_h5ad(geuvadis_genotype_preliminary_adata)

#* Extract transcriptome variants using dbSNP IDs (dbSNP --> HGVSC) and GTF (HGVSC CDS --> cDNA)
if not os.path.exists(variants) or not os.path.exists(geuvadis_genotype_true_adata) or not os.path.exists(geuvadis_true_vcf):  # transcriptome variants don't exist
    # convert preliminary VCF to true VCF
    if not os.path.exists(geuvadis_preliminary_vcf_exons_only_df_path):
        logger.info("Making vcf dataframe")
        variants_plink_df = vcf_to_dataframe(geuvadis_preliminary_vcf_exons_only, additional_columns=False, explode_alt=False, filter_empty_alt=True, verbose=True, total=total_entries)
        variants_plink_df.to_parquet(geuvadis_preliminary_vcf_exons_only_df_path, index=False)
    else:
        logger.info("Loading vcf dataframe")
        variants_plink_df = pd.read_parquet(geuvadis_preliminary_vcf_exons_only_df_path)
    # variants_plink_df.rename(columns={"REF": "MINOR", "ALT": "MAJOR"}, inplace=True)
    variants_plink_df["must_swap_id"] = variants_plink_df["ID"].isna()
    variants_plink_df["ID"] = variants_plink_df["ID"].fillna(variants_plink_df.index.to_series().apply(lambda i: f"temp_{i}"))  #$ matches convert_vcf_samples_to_anndata

    if not os.path.exists(gtf_df_path) or not os.path.exists(transcript_df_path):
        logger.info("Making gtf and transcript dfs")
        gtf_df, transcript_df = make_transcript_df_from_gtf(reference_genome_gtf)
        gtf_df.to_parquet(gtf_df_path, index=False)
        transcript_df.to_parquet(transcript_df_path, index=False)
    else:
        logger.info("Loading gtf and transcript dfs")
        gtf_df = pd.read_parquet(gtf_df_path)
        transcript_df = pd.read_parquet(transcript_df_path)

    #* Make hgvsc_to_dbsnp_dict mapping
    # # TODO: account for swapping REF and ALT as needed (which I do below, but I would just need to move up)
    # # With local VEP
    # if not os.path.exists(genome_fasta):
    #     logger.info("Downloading DNA sequences")
    #     gget_ref_command = ["gget", "ref", "-w", "dna", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    #     subprocess.run(gget_ref_command, check=True)
    #     unzip_command = ["gunzip", f"{genome_fasta}.gz"]
    #     subprocess.run(unzip_command, check=True)
    # if not os.path.exists(f"{conda_path}/envs/{conda_environment}"):
    #     logger.info("Creating VEP conda environment")
    #     subprocess.run(f"conda create -y -n {conda_environment} && conda activate {conda_environment} && conda install -y -c conda-forge -c bioconda -c defaults ensembl-vep=113 htslib bcftools samtools ucsc-liftover perl-list-moreutils", check=True, shell=True)
    # vep_113_dir = os.path.join(vep_cache_dir, "homo_sapiens", "113_GRCh37")
    # if not os.path.exists(vep_113_dir):
    #     logger.info("Downloading VEP cache")
    #     subprocess.run(f"wget ftp://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz && mkdir -p {vep_cache_dir} && tar -xzf homo_sapiens_vep_113_GRCh37.tar.gz -C {vep_cache_dir}", check=True, shell=True)
    # if not os.path.exists(geuvadis_preliminary_vep_vcf):
    #     logger.info("Running VEP")
    #     subprocess.run(f"conda run -n {conda_environment} vep -i {geuvadis_preliminary_vcf_exons_only} --cache --offline --hgvs --vcf -o {geuvadis_preliminary_vep_vcf} --fasta {genome_fasta} --assembly GRCh37 --cache_version 113", check=True, shell=True)
    # # parse through geuvadis_preliminary_vep_vcf to get HGVSC to dbSNP mapping
    # logger.info("Making HGVSc:dbSNP mapping")
    # hgvsc_to_dbsnp_dict, ids_not_correctly_recorded = make_hgvsc_to_dbsnp_dict_from_vep_vcf(geuvadis_preliminary_vep_vcf, total_entries=total_entries)
    
    # With Ensembl REST API (Ensembl release 113 as of the date of running 4/8/2025)
    url = "https://grch37.rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json"}
    hgvsc_to_dbsnp_dict = {}
    ids_not_correctly_recorded = set()
    chunk_size = 200

    saved_chunk_number = 0
    if os.path.isfile(hgvsc_to_dbsnp_pkl_path):
        logger.info("Loading saved HGVSC to dbSNP mapping from pickle file")
        with open(hgvsc_to_dbsnp_pkl_path, "rb") as f:
            hgvsc_to_dbsnp_dict = pickle.load(f)
        saved_chunk_number = hgvsc_to_dbsnp_dict["chunk"]

    ids_list = variants_plink_df['ID'].tolist()
    num_chunks = (len(ids_list) + chunk_size - 1) // chunk_size
    hgvsc_to_dbsnp_dict_alread_completed = hgvsc_to_dbsnp_dict.get("chunk", 0) == num_chunks-1
    if not hgvsc_to_dbsnp_dict_alread_completed:
        logger.info("Getting HGVSC from dbSNP IDs. Warning: This may take a while.")
        for chunk_number, ids_chunk in enumerate(tqdm(chunks(ids_list, chunk_size), total=num_chunks, desc="Processing batches")):  #* go through 200 at a time
            if chunk_number <= saved_chunk_number:  #* skip already processed chunks
                continue
            
            # Create the payload with the list of IDs and additional parameters (e.g., hgvs)
            payload = {"ids": ids_chunk, "hgvs": 1}
            
            try:
                response = requests.post(url, headers=headers, data=json.dumps(payload), timeout=90)
                response.raise_for_status()  # Raise an exception for HTTP errors

                hgvsc_to_dbsnp_dict, ids_not_correctly_recorded = add_to_hgvsc_dict(response, hgvsc_to_dbsnp_dict, ids_not_correctly_recorded)
            except Exception as e:  #* go through 20 at a time
                # print("Outer exception")
                print(f"Entering subchunk for chunk {chunk_number}")
                # If an error occurs, log all ids in the current chunk as not recorded
                for ids_subchunk in list(chunks(ids_chunk, (chunk_size//10)+1)):
                    payload = {"ids": ids_subchunk, "hgvs": 1}
                    try:
                        response = requests.post(url, headers=headers, data=json.dumps(payload), timeout=60)
                        response.raise_for_status()  # Raise an exception for HTTP errors

                        hgvsc_to_dbsnp_dict, ids_not_correctly_recorded = add_to_hgvsc_dict(response, hgvsc_to_dbsnp_dict, ids_not_correctly_recorded)
                    except Exception as e:  #* go through 1 at a time
                        # print("Inner exception")
                        if (chunk_size//10)+1 == 1:
                            ids_not_correctly_recorded.add(ids_subchunk[0])  # both cases where there was no dbsnp ID at all, as well as any other causes of failure
                        else:
                            # If an error occurs, log all ids in the current chunk as not recorded
                            for dbsnp_id in ids_chunk:
                                # ids_not_correctly_recorded.add(dbsnp_id)
                                try:
                                    url_single = f"{url}/{dbsnp_id}?hgvs=1"
                                    response = requests.get(url_single, headers=headers, timeout=10)
                                    response.raise_for_status()  # Raise an exception for HTTP errors

                                    hgvsc_to_dbsnp_dict, ids_not_correctly_recorded = add_to_hgvsc_dict(response, hgvsc_to_dbsnp_dict, ids_not_correctly_recorded)
                                except Exception as e:  #* give up
                                    ids_not_correctly_recorded.add(dbsnp_id)  # both cases where there was no dbsnp ID at all, as well as any other causes of failure
            
            if chunk_number % 100 == 0:
                hgvsc_to_dbsnp_dict["chunk"] = chunk_number
                with open(hgvsc_to_dbsnp_pkl_path, "wb") as f:
                    pickle.dump(hgvsc_to_dbsnp_dict, f)

        hgvsc_to_dbsnp_dict["chunk"] = num_chunks-1
        with open(hgvsc_to_dbsnp_pkl_path, "wb") as f:
            pickle.dump(hgvsc_to_dbsnp_dict, f)
    hgvsc_to_dbsnp_dict.pop("chunk", None)

    variants_plink_df_subset = None
    if ids_not_correctly_recorded:
        pass  # can uncomment later once I address the TODO
    #     logger.info("Some dbSNP IDs were not correctly recorded, fetching HGVSC from GTF")
    #     variants_plink_df_subset = variants_plink_df.loc[variants_plink_df['ID'].isin(ids_not_correctly_recorded), ["CHROM", "POS", "ID", "REF", "ALT"]].copy()
    #     variants_plink_df_subset.rename(columns={"CHROM": "chromosome", "POS": "start_variant_position_genome"}, inplace=True)
    #     variants_plink_df_subset["variant_source"] = "genome"
    #     variants_plink_df_subset['end_variant_position_genome'] = variants_plink_df_subset['start_variant_position_genome'] + variants_plink_df_subset['ALT'].str.len() - 1
    #     new_cols = variants_plink_df_subset.progress_apply(lambda row: assign_transcript_and_cds(row, transcript_df, gtf_df), axis=1)  # adds columns ["transcript_ID", "start_variant_position", "end_variant_position"]
    #     variants_plink_df_subset = variants_plink_df_subset.merge(new_cols[["transcript_ID", "start_variant_position", "end_variant_position"]], left_index=True, right_index=True)
    #     variants_plink_df_subset.rename(columns={"chromosome": "CHROM", "start_variant_position_genome": "POS"}, inplace=True)
    #     del new_cols
        
    #     # drop rows that didn't map to a gene
    #     variants_plink_df_subset = variants_plink_df_subset[
    #         (variants_plink_df_subset["transcript_ID"] != "unknown") &
    #         (variants_plink_df_subset["start_variant_position"].notna()) &
    #         (variants_plink_df_subset["end_variant_position"].notna())
    #     ]
        
    #     # make HGVSC from this
    #     if len(variants_plink_df_subset) > 0:
    #         variants_plink_df_subset = variants_plink_df_subset.merge(
    #             transcript_df[["transcript_ID", "strand"]],
    #             on="transcript_ID",
    #             how="left"
    #         ).reset_index(drop=True)

    #         # TODO: account for swapping REF and ALT as needed (which I do below, but I would just need to move up) - look at CHROM and POS, do a lookup of this base; if base == REF[0], then keep; elif base == ALT[0], then swap; else give error (and likely exit whole function because it indicates systemic error)
                #* note that the above method only works for substitution/delins - cannot work for deletion or insertion because the REF[0] base will be correct regardless - I can look to see if the entirety of REF matches the reference genome for deletion, and I can look to see if the entirety of ALT matches the reference genome for insertion; if not a perfect match, then I know insertion; but if it is a perfect match, then I don't know whether deletion or duplication, so I'm still stumped
            
    #         # ensure REF and ALT reflect what is going on on the cDNA
    #         variants_plink_df_subset["pos_original"] = variants_plink_df_subset["POS"]  # just as a copy, in case I want it later
    #         variants_plink_df_subset["ref_original"] = variants_plink_df_subset["REF"]  # just as a copy, in case I want it later
    #         variants_plink_df_subset["alt_original"] = variants_plink_df_subset["ALT"]  # just as a copy, in case I want it later
    #         variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "REF"] = variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "REF"].apply(reverse_complement)
    #         variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "ALT"] = variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "ALT"].apply(reverse_complement)
    #         variants_plink_df_subset["POS"] = variants_plink_df_subset["start_variant_position"].astype(int)

            
    #         # create HGVSC info
    #         add_variant_type_column_to_vcf_derived_df(variants_plink_df_subset)
    #         add_variant_column_to_vcf_derived_df(variants_plink_df_subset, var_column="variant")
    #         variants_plink_df_subset['variant'] = variants_plink_df_subset['variant'].str.replace(r'^g\.', 'c.', regex=True)
    #         variants_plink_df_subset["variant_header"] = variants_plink_df_subset["transcript_ID"] + ":" + variants_plink_df_subset["variant"]
    #         variants_plink_df_subset = variants_plink_df_subset.dropna(subset=["variant_header"])
    #         variants_plink_df_subset = variants_plink_df_subset.rename(columns={"ID": "dbsnp_id"})
    # del variants_plink_df

    if not os.path.exists(variants):
        logger.info("Creating hgvs_df")
        hgvs_df = pd.DataFrame(
            [(v, k) for k, v in hgvsc_to_dbsnp_dict.items()],
            columns=['dbsnp_id', 'variant_header']
        )
        
        if variants_plink_df_subset is not None and len(variants_plink_df_subset) > 0:
            hgvs_df = pd.concat([hgvs_df, variants_plink_df_subset], ignore_index=True)  # add the ones not correctly recorded
            del variants_plink_df_subset

        hgvs_df[["transcript_ID", "variant"]] = hgvs_df["variant_header"].str.split(":", expand=True)
        hgvs_df["transcript_ID"] = hgvs_df["transcript_ID"].str.replace(r"\.\d+$", "", regex=True)
        hgvs_df[["nucleotide_positions", "actual_variant"]] = hgvs_df["variant"].str.extract(mutation_pattern)
        hgvs_df["after_ins_or_gt"] = hgvs_df["actual_variant"].str.extract(r"(?:ins|>)(.+)")
        hgvs_df['ALT_CDNA_TRUE'] = hgvs_df['actual_variant'].str.extract(r'(?:>|ins)([A-Za-z]+)')
        # split_positions = hgvs_df["nucleotide_positions"].str.split("_", expand=True)
        # hgvs_df["start_variant_position"] = split_positions[0]
        # if split_positions.shape[1] > 1:
        #     hgvs_df["end_variant_position"] = split_positions[1].fillna(split_positions[0])
        # else:
        #     hgvs_df["end_variant_position"] = hgvs_df["start_variant_position"]
        # hgvs_df.loc[hgvs_df["end_variant_position"].isna(), "end_variant_position"] = hgvs_df["start_variant_position"]
        # hgvs_df[["start_variant_position", "end_variant_position"]] = hgvs_df[["start_variant_position", "end_variant_position"]].astype(int)

        logger.info("Assigning transcript and CDS positions")
        hgvs_df = hgvs_df.merge(
            transcript_df[["transcript_ID", "utr_length_preceding_transcript", "strand"]],
            on="transcript_ID",
            how="left"
        ).reset_index(drop=True)
        del transcript_df

        # variant_cdna_NOT_nan_mask = hgvs_df[["start_variant_position", "end_variant_position", "utr_length_preceding_transcript"]].notna().all(axis=1)
        # hgvs_df.loc[variant_cdna_NOT_nan_mask, "start_variant_position_cdna"] = (hgvs_df.loc[variant_cdna_NOT_nan_mask, "start_variant_position"] + hgvs_df.loc[variant_cdna_NOT_nan_mask, "utr_length_preceding_transcript"])
        # hgvs_df.loc[variant_cdna_NOT_nan_mask, "end_variant_position_cdna"] = (hgvs_df.loc[variant_cdna_NOT_nan_mask, "end_variant_position"] + hgvs_df.loc[variant_cdna_NOT_nan_mask, "utr_length_preceding_transcript"])

        # hgvs_df["start_variant_position_cdna"] = hgvs_df["start_variant_position_cdna"].astype("Int64")
        # hgvs_df["end_variant_position_cdna"] = hgvs_df["end_variant_position_cdna"].astype("Int64")

        # same_pos = variant_cdna_NOT_nan_mask & (hgvs_df["start_variant_position_cdna"] == hgvs_df["end_variant_position_cdna"])
        # diff_pos = variant_cdna_NOT_nan_mask & (hgvs_df["start_variant_position_cdna"] != hgvs_df["end_variant_position_cdna"])
        # hgvs_df.loc[same_pos, "variant_cdna"] = ("c." + hgvs_df.loc[same_pos, "start_variant_position_cdna"].astype(str) + hgvs_df.loc[same_pos, "actual_variant"])
        # hgvs_df.loc[diff_pos, "variant_cdna"] = ("c." + hgvs_df.loc[diff_pos, "start_variant_position_cdna"].astype(str) + "_" + hgvs_df.loc[diff_pos, "end_variant_position_cdna"].astype(str) + hgvs_df.loc[diff_pos, "actual_variant"])

        # # add in the missing cDNA positions
        # hgvs_df_missing_cdna = hgvs_df[hgvs_df[["utr_length_preceding_transcript"]].isna().any(axis=1)].copy()
        hgvs_df_missing_cdna = hgvs_df.copy()
        hgvs_df_missing_cdna = hgvs_df_missing_cdna[["variant_header", "transcript_ID", "variant"]].rename(columns={"transcript_ID": "seq_ID", "variant": "mutation"})
        hgvs_df_missing_cdna, _ = convert_mutation_cds_locations_to_cdna(hgvs_df_missing_cdna, cdna_fasta_path=sequences, cds_fasta_path=cds_path, output_csv_path=None, verbose=True, strip_leading_Ns_cds=False)
        hgvs_df = hgvs_df.merge(hgvs_df_missing_cdna, on="variant_header", how="left", suffixes=("", "_missing_cdna"))
        # hgvs_df["variant_cdna"] = hgvs_df["mutation_cdna"].fillna(hgvs_df["mutation_cdna"])
        hgvs_df["variant_cdna"] = hgvs_df["mutation_cdna"]
        hgvs_df.drop(columns=["seq_ID", "mutation", "mutation_cdna"], inplace=True)
        del hgvs_df_missing_cdna

        hgvs_df = hgvs_df[~hgvs_df["variant_cdna"].str.startswith("c.nan", na=False)]

        # make variant_type based actual_variant
        add_variant_type(hgvs_df, var_column="variant")

        hgvs_df.to_parquet(variants, index=False)
    else:
        logger.info("Loading transcriptome variants into hgvs_df")
        hgvs_df = pd.read_parquet(variants)

    if not os.path.exists(geuvadis_genotype_true_adata) or not os.path.exists(geuvadis_true_vcf):
        # update adata and make true VCF
        genotype_swap_dict = {0: 2, 1: 1, 2: 0}

        # make variant_type_plink column that is substitution if REF and ALT 1, deletion if REF > 1 and ALT == 1, insertion if REF == 1 and ALT > 1, and delins if REF > 1 and ALT > 1
        conditions = [
            (variants_plink_df['REF'].str.len() == 1) & (variants_plink_df['ALT'].str.len() == 1),
            (variants_plink_df['REF'].str.len() > 1) & (variants_plink_df['ALT'].str.len() == 1),
            (variants_plink_df['REF'].str.len() == 1) & (variants_plink_df['ALT'].str.len() > 1),
            (variants_plink_df['REF'].str.len() > 1) & (variants_plink_df['ALT'].str.len() > 1),
        ]

        choices = ['substitution', 'deletion', 'insertion', 'delins']
        variants_plink_df['variant_type_plink'] = np.select(conditions, choices, default='unknown')

        # merge in strand from hgvs_df
        logger.info("Merging hgvs_df into variants_plink_df")
        number_of_variants_in_df_originally = len(variants_plink_df)

        variants_plink_df = variants_plink_df.merge(
            hgvs_df[["dbsnp_id", "strand"]].drop_duplicates(subset=["dbsnp_id"], keep="first").reset_index(drop=True),
            left_on="ID",
            right_on="dbsnp_id",
            how="left"
        ).reset_index(drop=True)
        variants_plink_df = variants_plink_df.dropna(subset=["dbsnp_id"])
        variants_plink_df = variants_plink_df.drop(columns="dbsnp_id")

        # make ALT_CDNA_PLINK which is reverse-complement of ALT for - strand and equal to ALT for + strand
        variants_plink_df["ALT_CDNA_PLINK"] = variants_plink_df["ALT"]
        variants_plink_df.loc[variants_plink_df["strand"] == "-", "ALT_CDNA_PLINK"] = variants_plink_df.loc[variants_plink_df["strand"] == "-", "ALT_CDNA_PLINK"].apply(reverse_complement)
        variants_plink_df["REF_CDNA_PLINK"] = variants_plink_df["REF"]
        variants_plink_df.loc[variants_plink_df["strand"] == "-", "REF_CDNA_PLINK"] = variants_plink_df.loc[variants_plink_df["strand"] == "-", "REF_CDNA_PLINK"].apply(reverse_complement)

        # optional - get rid of alternative splice transcript isoforms (and get rid of deletions that span different bases but have the same dbsnp_id, which should be zero but is still a nice sanity check to get for free)
        hgvs_df = hgvs_df.drop_duplicates(subset=["dbsnp_id", "ALT_CDNA_TRUE"], keep="first")
        
        variants_plink_df["ALT_CDNA_PLINK_no_first_base"] = variants_plink_df["ALT_CDNA_PLINK"].str[1:]
        variants_plink_df["REF_CDNA_PLINK_no_first_base"] = variants_plink_df["REF_CDNA_PLINK"].str[1:]

        # do my merges
        variants_plink_df["must_swap_genotype"] = False
        #* Case 1: Insertion, plink was correct
        logger.info("Case 1: Insertion, plink was correct")
        insertion_merged_correct = variants_plink_df[variants_plink_df["variant_type_plink"] == "insertion"].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on=["ID", "ALT_CDNA_PLINK_no_first_base"],
            right_on=["dbsnp_id", "ALT_CDNA_TRUE"],
            how="left"
        ).reset_index(drop=True)
        insertion_merged_correct = insertion_merged_correct.dropna(subset=["variant_header"])
        insertion_merged_correct = insertion_merged_correct.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        #* Case 2: Insertion, plink was incorrect
        logger.info("Case 2: Insertion, plink was incorrect")
        insertion_merged_wrong = variants_plink_df[variants_plink_df["variant_type_plink"] == "deletion"].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on=["ID", "REF_CDNA_PLINK_no_first_base"],
            right_on=["dbsnp_id", "ALT_CDNA_TRUE"],
            how="left"
        ).reset_index(drop=True)
        insertion_merged_wrong = insertion_merged_wrong.dropna(subset=["variant_header"])
        insertion_merged_wrong = insertion_merged_wrong.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        insertion_merged_wrong["must_swap_genotype"] = True
        #* Case 3: Deletion, plink was correct
        logger.info("Case 3: Deletion, plink was correct")
        deletion_merged_correct = variants_plink_df[variants_plink_df["variant_type_plink"] == "deletion"].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on="ID",
            right_on="dbsnp_id",
            how="left"
        ).reset_index(drop=True)
        deletion_merged_correct = deletion_merged_correct.dropna(subset=["variant_header"])
        deletion_merged_correct = deletion_merged_correct.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        #* Case 4: Deletion, plink was incorrect
        logger.info("Case 4: Deletion, plink was incorrect")
        deletion_merged_wrong = variants_plink_df[variants_plink_df["variant_type_plink"] == "insertion"].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on="ID",
            right_on="dbsnp_id",
            how="left"
        ).reset_index(drop=True)
        deletion_merged_wrong = deletion_merged_wrong.dropna(subset=["variant_header"])
        deletion_merged_wrong = deletion_merged_wrong.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        deletion_merged_wrong["must_swap_genotype"] = True
        #* Case 5: Substitution/delins, plink was correct
        logger.info("Case 5: Substitution/delins, plink was correct")
        substitution_merged_correct = variants_plink_df[(variants_plink_df["variant_type_plink"] == "substitution") | (variants_plink_df["variant_type_plink"] == "delins")].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on=["ID", "ALT_CDNA_PLINK"],
            right_on=["dbsnp_id", "ALT_CDNA_TRUE"],
            how="left"
        ).reset_index(drop=True)
        substitution_merged_correct = substitution_merged_correct.dropna(subset=["variant_header"])
        substitution_merged_correct = substitution_merged_correct.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        #* Case 6: Substitution/delins, plink was incorrect
        logger.info("Case 6: Substitution/delins, plink was incorrect")
        substitution_merged_wrong = variants_plink_df[(variants_plink_df["variant_type_plink"] == "substitution") | (variants_plink_df["variant_type_plink"] == "delins")].merge(
            hgvs_df[["variant_header", "ALT_CDNA_TRUE", "variant_type", "dbsnp_id"]],
            left_on=["ID", "REF_CDNA_PLINK"],
            right_on=["dbsnp_id", "ALT_CDNA_TRUE"],
            how="left"
        ).reset_index(drop=True)
        substitution_merged_wrong = substitution_merged_wrong.dropna(subset=["variant_header"])
        substitution_merged_wrong = substitution_merged_wrong.drop_duplicates(subset=["ID", "POS", "REF", "ALT"])
        substitution_merged_wrong["must_swap_genotype"] = True

        # variants_plink_df_copy = variants_plink_df.copy()
        variants_plink_df = pd.concat([insertion_merged_correct, insertion_merged_wrong, deletion_merged_correct, deletion_merged_wrong, substitution_merged_correct, substitution_merged_wrong], ignore_index=True)
        variants_plink_df = variants_plink_df.dropna(subset=["variant_header"])  # should be 0, but still a sanity check

        # optional - drop rows where dbsnp ID is duplicated, keeping only those where must_swap_genotype is False (because these represent the vast majority of cases, I trust that this is what Geuvadis meant)
        variants_plink_df = variants_plink_df.sort_values(
            by=["must_swap_genotype"], ascending=True
        )
        variants_plink_df = variants_plink_df.drop_duplicates(subset=["ID"], keep="first").reset_index(drop=True)
        assert len(variants_plink_df[variants_plink_df.duplicated(subset=["ID"], keep=False)]) == 0, "Duplicate IDs found in variants_plink_df after merging with hgvs_df and allegedly controlling for all duplicates"

        if len(variants_plink_df) < number_of_variants_in_df_originally:
            valid_ids = set(variants_plink_df["ID"])
            adata = adata[:, adata.var.index.isin(valid_ids)].copy()

        final_headers = set(variants_plink_df["variant_header"].tolist())
        hgvs_df = hgvs_df[hgvs_df["variant_header"].isin(final_headers)]
        hgvs_df.to_parquet(variants, index=False)
        
        # # old, bad code
        # # make variant_type_plink based on preliminary vcf
        # variants_plink_df["variant_type_plink"] = variants_plink_df.apply(
        #     lambda row: "sub_or_delins"
        #     if len(str(row["ALT"])) == len(str(row["REF"]))
        #     else ("insertion" if len(str(row["ALT"])) > len(str(row["REF"])) else "deletion"),
        #     axis=1
        # )

        # # convert duplication → insertion, and inversion → delins
        # variants_plink_df.loc[variants_plink_df["variant_type"] == "duplication", "variant_type"] = "insertion"
        # variants_plink_df.loc[variants_plink_df["variant_type"] == "inversion", "variant_type"] = "delins"
        
        # # make must_swap_genotype column as needed, i.e., when ALT_CDNA_PLINK != ALT_CDNA_TRUE (subs and delins), or where deletion and insertion are swapped
        # sub_or_delins_mask = variants_plink_df["variant_type_plink"] == "sub_or_delins"
        # variants_plink_df.loc[sub_or_delins_mask, "must_swap_genotype"] = variants_plink_df.loc[sub_or_delins_mask, "ALT_CDNA_PLINK"] != variants_plink_df.loc[sub_or_delins_mask, "ALT_CDNA_TRUE"]
        # variants_plink_df.loc[~sub_or_delins_mask, "must_swap_genotype"] = variants_plink_df.loc[~sub_or_delins_mask, "variant_type"] != variants_plink_df.loc[~sub_or_delins_mask, "variant_type_plink"]
        
        if not os.path.exists(geuvadis_genotype_true_adata):
            logger.info("Updating adata with true genotypes")
            # remove variants/cols from adata that aren't in variants_plink_df
            variant_ids_to_keep = set(variants_plink_df["ID"].unique())
            keep_mask = adata.var.index.isin(variant_ids_to_keep)
            adata = adata[:, keep_mask].copy()
            
            # Ensure adata.var.index is aligned with 'ID' from variants_plink_df
            swap_ids = variants_plink_df.loc[variants_plink_df["must_swap_genotype"], "ID"].values

            # Get the column indices in adata corresponding to those IDs
            swap_cols = adata.var.index.get_indexer(swap_ids)

            # Create a vectorized function for mapping
            swap_func = np.vectorize(genotype_swap_dict.get)

            # Apply the mapping column by column
            adata.X = adata.X.tolil()
            for col in tqdm(swap_cols, total=len(swap_cols), desc="Swapping genotypes"):
                if col != -1:  # skip if ID not found
                    adata.X[:, col] = swap_func(adata.X[:, col].toarray())
            adata.X = adata.X.tocsr()

            # add the HGVSC header to adata.var, and use HGVSC in place of dbsnp when dbsnp not available
            adata.var = adata.var.merge(variants_plink_df[["ID", "variant_header"]], how="left", left_index=True, right_on="ID")  #$ if I ever implement the variants_plink_df_subset uncommenting, then don't just merge by ID - instead, ensure that adata and variants_plink_df are aligned and then simply set    
            adata.write_h5ad(geuvadis_genotype_true_adata)
        elif os.path.exists(geuvadis_genotype_true_adata) and save_vcf_samples:
            adata = ad.read_h5ad(geuvadis_genotype_true_adata)

        if not os.path.exists(geuvadis_true_vcf):
            logger.info("Writing true VCF")
            variants_plink_df[["REF_original", "ALT_original"]] = variants_plink_df[["REF", "ALT"]]

            # swap REF and ALT as needed
            swap_mask = variants_plink_df["must_swap_genotype"]
            temp = variants_plink_df.loc[swap_mask, "REF"]
            variants_plink_df.loc[swap_mask, "REF"] = variants_plink_df.loc[swap_mask, "ALT"]
            variants_plink_df.loc[swap_mask, "ALT"] = temp
            del temp

            # write to a true VCF now
            adata_for_vcf = adata if save_vcf_samples else None
            write_to_vcf(variants_plink_df, output_file=geuvadis_true_vcf, save_vcf_samples=save_vcf_samples, adata=adata_for_vcf)

#* Run vk ref
for w_and_k_dict in w_and_k_list_of_dicts:
    w, k = w_and_k_dict["w"], w_and_k_dict["k"]
    logger.info(f"Running vk.ref with w={w} and k={k}")
    vk_ref_out = os.path.join(vk_ref_out_parent, f"w{w}_k{k}")
    vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
    vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")
    vk.ref(
        variants=variants,
        sequences=sequences,
        w=w,
        k=k,
        seq_id_column=seq_id_column,
        var_column=var_column,
        out=vk_ref_out,
        reference_out_dir=reference_out_dir,
        dlist_reference_source=dlist_reference_source,
        index_out=vcrs_index,
        t2g_out=vcrs_t2g,
        threads=threads,
        chunksize=chunksize,
        save_logs=False,
        verbose=True,
        merge_identical=False
    )