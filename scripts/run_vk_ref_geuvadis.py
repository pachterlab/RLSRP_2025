import os
import subprocess
import numpy as np
import requests
import anndata as ad
import pandas as pd
import json
import varseek as vk
from varseek.utils import is_program_installed, vcf_to_dataframe, assign_transcript_and_cds, add_variant_type, reverse_complement, write_to_vcf, add_variant_type_column_to_vcf_derived_df, add_variant_column_to_vcf_derived_df
from varseek.constants import mutation_pattern
from RLSRWP_2025.seq_utils import make_transcript_df_from_gtf, convert_vcf_samples_to_anndata

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
variants = os.path.join(geuvadis_reference_files_dir, "variants_transcriptome.csv")
sequences = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.cdna.all.fa")
reference_genome_gtf = os.path.join(ensembl_reference_files_dir, "Homo_sapiens.GRCh37.87.gtf")

# # parameters
w=47
k=51
seq_id_column = "transcript_ID"
var_column = "variant"
dlist_reference_source = "t2t"
threads = 4
chunksize = None

# # output paths
vk_ref_out = os.path.join(data_dir, "vk_ref_out_geuvadis")
vcrs_index = os.path.join(vk_ref_out, "vcrs_index.idx")
vcrs_t2g = os.path.join(vk_ref_out, "vcrs_t2g_filtered.txt")

# # misc
geuvadis_reference_variants_prefix = os.path.join(geuvadis_reference_files_dir, "1kg_phase1_all")
geuvadis_preliminary_vcf = f"{geuvadis_reference_variants_prefix}_preliminary.vcf.gz"
geuvadis_true_vcf = f"{geuvadis_reference_variants_prefix}_true.vcf.gz"
geuvadis_genotype_preliminary_adata = os.path.join(geuvadis_reference_files_dir, "genotypes_adata_preliminary.h5ad")
geuvadis_genotype_true_adata = os.path.join(geuvadis_reference_files_dir, "genotypes_adata_true.h5ad")
save_vcf_samples = False
plink = "plink"


geuvadis_rna_json_file = os.path.join(geuvadis_reference_files_dir, "geuvadis_metadata.json")

#* Download geuvadis RNA metadata (for extracting sample names)
if not os.path.exists(geuvadis_rna_json_file):
    json_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,scientific_name,library_strategy,experiment_title,experiment_alias,fastq_bytes,fastq_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"
    sequencing_metadata_download_command = ["wget", "-q", "-O", geuvadis_rna_json_file, json_url]
    if os.path.dirname(geuvadis_rna_json_file) and not os.path.exists(os.path.dirname(geuvadis_rna_json_file)):
        os.makedirs(os.path.dirname(geuvadis_rna_json_file), exist_ok=True)
    subprocess.run(sequencing_metadata_download_command, check=True)

#* download cDNA
os.makedirs(ensembl_reference_files_dir, exist_ok=True)
if not os.path.exists(sequences):
    gget_ref_command = ["gget", "ref", "-w", "cdna", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    unzip_command = ["gunzip", f"{sequences}.gz"]
    subprocess.run(unzip_command, check=True)

#* download gtf
if not os.path.exists(reference_genome_gtf):
    gget_ref_command = ["gget", "ref", "-w", "gtf", "-r", "113", "--out_dir", ensembl_reference_files_dir, "-d", "human_grch37"]
    subprocess.run(gget_ref_command, check=True)
    unzip_command = ["gunzip", f"{reference_genome_gtf}.gz"]
    subprocess.run(unzip_command, check=True)

#* download Geuvadis reference variants
if not os.path.exists(f"{geuvadis_reference_variants_prefix}.bim") or not os.path.exists(f"{geuvadis_reference_variants_prefix}.bed") or not os.path.exists(f"{geuvadis_reference_variants_prefix}.fam"):
    plink_file_download_command = ["wget", "-O", f"{geuvadis_reference_variants_prefix}.tar.gz", "https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz"]
    subprocess.run(plink_file_download_command, check=True)
    untar_command = ["tar", "-xzvf", f"{geuvadis_reference_variants_prefix}.tar.gz", "-C", geuvadis_reference_files_dir]
    subprocess.run(untar_command, check=True)

#* Convert plink variant files to preliminary VCF
if not os.path.exists(geuvadis_preliminary_vcf):  # plink pseudo-VCF doesn't exist
    if not is_program_installed(plink):
        raise ValueError(f"plink is not installed. Please install plink to run this script.")
    plink_to_preliminary_vcf_command = [plink, "--bfile", "/data/geuvadis_data_base/genome_variants_ground_truth/1kg_phase1_all", "--threads", str(threads), "--recode", "vcf", "bgz", "--out", geuvadis_reference_variants_prefix]
    subprocess.run(plink_to_preliminary_vcf_command, check=True)

#* make adata
if not os.path.exists(geuvadis_genotype_preliminary_adata):
    with open(geuvadis_rna_json_file, 'r', encoding="utf-8") as file:
        data = json.load(file)

    rna_df = pd.DataFrame(data)
    sample_titles_set = set(rna_df['sample_title'].tolist())

    adata = convert_vcf_samples_to_anndata(geuvadis_preliminary_vcf, sample_titles_set=sample_titles_set, adata_out=geuvadis_genotype_preliminary_adata)
else:
    adata = ad.read_h5ad(geuvadis_genotype_preliminary_adata)

#* Extract transcriptome variants using dbSNP IDs (dbSNP --> HGVSC) and GTF (HGVSC CDS --> cDNA)
if not os.path.exists(variants):  # transcriptome variants don't exist
    # convert preliminary VCF to true VCF
    variants_plink_df = vcf_to_dataframe(geuvadis_preliminary_vcf, additional_columns=False, explode_alt=False, filter_empty_alt=True, verbose=True)
    # variants_plink_df.rename(columns={"REF": "MINOR", "ALT": "MAJOR"}, inplace=True)
    variants_plink_df["must_swap_id"] = variants_plink_df["ID"].isna()
    variants_plink_df["ID"] = variants_plink_df["ID"].fillna(variants_plink_df.index.to_series().apply(lambda i: f"temp_{i}"))  #$ matches convert_vcf_samples_to_anndata

    gtf_df, transcript_df = make_transcript_df_from_gtf(reference_genome_gtf)
    
    # convert dbSNP IDs --> HGVSC with Ensembl REST API (Ensembl release 113 as of the date of running 4/8/2025)
    headers = {"Content-Type": "application/json"}
    hgvsc_to_dbsnp_dict = {}
    ids_not_correctly_recorded = set()
    for dbsnp_id in variants_plink_df['ID']:
        try:
            url = f"https://grch37.rest.ensembl.org/vep/human/id/{dbsnp_id}?hgvs=1"
            response = requests.get(url, headers=headers)
            hgvsc = None
            for entry in response.json():
                for tx in entry.get("transcript_consequences", []):
                    if "hgvsc" in tx:
                        hgvsc = tx["hgvsc"]
                        if ":n." in hgvsc or "-" in hgvsc or "+" in hgvsc:
                            continue
                        hgvsc_to_dbsnp_dict[hgvsc] = dbsnp_id
        except Exception as e:
            # print(f"Error fetching data for dbSNP ID {dbsnp_id}: {e}")
            ids_not_correctly_recorded.add(dbsnp_id)  # both cases where there was no dbsnp ID at all, as well as any other causes of failure
            continue
    if ids_not_correctly_recorded:
        variants_plink_df_subset = variants_plink_df.loc[variants_plink_df['ID'].isin(ids_not_correctly_recorded), ["CHROM", "POS", "ID", "REF", "ALT"]].copy()
        variants_plink_df_subset.rename(columns={"CHROM": "chromosome", "POS": "start_variant_position_genome"}, inplace=True)
        variants_plink_df_subset['end_variant_position_genome'] = variants_plink_df_subset['start_variant_position_genome'] + variants_plink_df_subset['ALT'].str.len() - 1
        variants_plink_df_subset = variants_plink_df_subset.apply(lambda row: assign_transcript_and_cds(row, transcript_df, gtf_df), axis=1)  # adds columns ["transcript_ID", "start_variant_position", "end_variant_position"]
        variants_plink_df_subset.rename(columns={"chromosome": "CHROM", "start_variant_position_genome": "POS"}, inplace=True)
        
        # drop rows that didn't map to a gene
        variants_plink_df_subset = variants_plink_df_subset[variants_plink_df_subset["transcript_ID"] != "unknown"]
        
        # make HGVSC from this
        if len(variants_plink_df_subset) > 0:
            variants_plink_df_subset = variants_plink_df_subset.merge(
                transcript_df[["transcript_ID", "strand"]],
                on="transcript_ID",
                how="left"
            ).reset_index(drop=True)
            
            # ensure REF and ALT reflect what is going on on the cDNA
            variants_plink_df_subset["ref_original"] = variants_plink_df_subset["REF"]  # just as a copy, in case I want it later
            variants_plink_df_subset["alt_original"] = variants_plink_df_subset["ALT"]  # just as a copy, in case I want it later
            variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "REF"] = variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "REF"].apply(reverse_complement)
            variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "ALT"] = variants_plink_df_subset.loc[variants_plink_df_subset["strand"] == "-", "ALT"].apply(reverse_complement)
            
            # create HGVSC info
            add_variant_type_column_to_vcf_derived_df(variants_plink_df_subset)
            add_variant_column_to_vcf_derived_df(variants_plink_df_subset, var_column="variant")
            variants_plink_df_subset['variant'] = variants_plink_df_subset['variant'].str.replace(r'^g\.', 'c.', regex=True)
            variants_plink_df_subset["variant_header"] = variants_plink_df_subset["transcript_ID"] + ":" + variants_plink_df_subset["variant"]
            variants_plink_df_subset = variants_plink_df_subset.dropna(subset=["variant_header"])
            variants_plink_df_subset = variants_plink_df_subset.rename(columns={"ID": "dbsnp_id"})
    # del variants_plink_df

    hgvs_df = pd.DataFrame(
        [(v, k) for k, v in hgvsc_to_dbsnp_dict.items()],
        columns=['dbsnp_id', 'variant_header']
    )
    if len(variants_plink_df_subset) > 0:
        hgvs_df = pd.concat([hgvs_df, variants_plink_df_subset], ignore_index=True)  # add the ones not correctly recorded
    del variants_plink_df_subset

    hgvs_df = hgvs_df.drop_duplicates(subset=['dbsnp_id'])  # optional - just if I don't want multiple transcripts per genomic mutation
    hgvs_df[["transcript_ID", "variant"]] = hgvs_df["variant_header"].str.split(":", expand=True)
    hgvs_df[["nucleotide_positions", "actual_variant"]] = hgvs_df[var_column].str.extract(mutation_pattern)
    split_positions = hgvs_df["nucleotide_positions"].str.split("_", expand=True)
    hgvs_df["start_variant_position"] = split_positions[0]
    if split_positions.shape[1] > 1:
        hgvs_df["end_variant_position"] = split_positions[1].fillna(split_positions[0])
    else:
        hgvs_df["end_variant_position"] = hgvs_df["start_variant_position"]
    hgvs_df.loc[hgvs_df["end_variant_position"].isna(), "end_variant_position"] = hgvs_df["start_variant_position"]
    hgvs_df[["start_variant_position", "end_variant_position"]] = hgvs_df[["start_variant_position", "end_variant_position"]].astype(int)

    hgvs_df = hgvs_df.merge(
        transcript_df[["transcript_ID", "utr_length_preceding_transcript", "strand"]],
        on="transcript_ID",
        how="left"
    ).reset_index(drop=True)
    del transcript_df

    hgvs_df["start_variant_position_cdna"] += hgvs_df["utr_length_preceding_transcript"]
    hgvs_df["end_variant_position_cdna"] += hgvs_df["utr_length_preceding_transcript"]

    hgvs_df.loc[hgvs_df["start_variant_position_cdna"] == hgvs_df["end_variant_position_cdna"], "variant_cdna"] = "c." + hgvs_df["start_variant_position_cdna"].astype(str) + hgvs_df["actual_variant"]
    hgvs_df.loc[hgvs_df["start_variant_position_cdna"] != hgvs_df["end_variant_position_cdna"], "variant_cdna"] = "c." + hgvs_df["start_variant_position_cdna"].astype(str) + "_" + hgvs_df["end_variant_position_cdna"].astype(str) + hgvs_df["actual_variant"]

    hgvs_df.to_parquet(variants, index=False)

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

    # merge in strand and actual_variant from hgvs_df
    variants_plink_df = variants_plink_df.merge(
        hgvs_df[["variant_header", "transcript_ID", "utr_length_preceding_transcript", "strand"]],
        on="transcript_ID",
        how="left"
    ).reset_index(drop=True)

    # make ALT_CDNA_PLINK which is reverse-complement of ALT for - strand and equal to ALT for + strand
    variants_plink_df["ALT_CDNA_PLINK"] = variants_plink_df["ALT"]
    variants_plink_df.loc[variants_plink_df["strand"] == "-", "ALT_CDNA_PLINK"] = variants_plink_df.loc[variants_plink_df["strand"] == "-", "ALT"].apply(reverse_complement)

    # make variant_type_true based actual_variant
    add_variant_type(variants_plink_df, var_column="actual_variant")

    # convert duplication → insertion, and inversion → delins
    variants_plink_df.loc[variants_plink_df["variant_type_true"] == "duplication", "variant_type_true"] = "insertion"
    variants_plink_df.loc[variants_plink_df["variant_type_true"] == "inversion", "variant_type_true"] = "delins"

    # make ALT_CDNA_TRUE column that extracts letters after >/ins
    variants_plink_df['ALT_CDNA_TRUE'] = variants_plink_df['actual_variant'].str.extract(r'(?:>|ins)([A-Za-z]+)')
    
    # make must_swap_genotype column as needed, i.e., when ALT_CDNA_PLINK != ALT_CDNA_TRUE:
    variants_plink_df["must_swap_genotype"] = variants_plink_df["ALT_CDNA_PLINK"] != variants_plink_df["ALT_CDNA_TRUE"]
    
    if not os.path.exists(geuvadis_genotype_true_adata):
        # remove variants/cols from adata that aren't in variants_plink_df
        variant_ids_to_keep = set(variants_plink_df["ID"].unique())
        keep_mask = adata.var.index.isin(variant_ids_to_keep)
        adata = adata[:, keep_mask].copy()
        
        # Ensure adata.var.index is aligned with 'ID' from variants_plink_df
        swap_ids = variants_plink_df.loc[variants_plink_df["must_swap_genotype"], "ID"]

        # Get the column indices in adata corresponding to those IDs
        swap_cols = adata.var.index.get_indexer(swap_ids)

        # Create a vectorized function for mapping
        swap_func = np.vectorize(genotype_swap_dict.get)

        # Apply the mapping column by column
        for col in swap_cols:
            if col != -1:  # skip if ID not found
                adata.X[:, col] = swap_func(adata.X[:, col])

        # add the HGVSC header to adata.var, and use HGVSC in place of dbsnp when dbsnp not available
        adata.var[["vcrs_header", "must_swap_id"]] = variants_plink_df[["vcrs_header", "must_swap_id"]]
        new_index = adata.var.index.to_series()
        mask = adata.var["must_swap_id"] == True
        new_index[mask] = adata.var.loc[mask, "vcrs_header"]
        adata.var.index = new_index
        adata.var.drop(columns=["must_swap_id"], inplace=True)
        
        adata.write_h5ad(geuvadis_genotype_true_adata)
    elif os.path.exists(geuvadis_genotype_true_adata) and save_vcf_samples:
        adata = ad.read_h5ad(geuvadis_genotype_true_adata)

    if not os.path.exists(geuvadis_true_vcf):
        variants_plink_df[["REF_original", "ALT_original"]] = variants_plink_df[["REF", "ALT"]]

        # swap REF and ALT as needed
        swap_mask = variants_plink_df["must_swap_genotype"]
        temp = variants_plink_df.loc[swap_mask, "REF"]
        variants_plink_df.loc[swap_mask, "REF"] = variants_plink_df.loc[swap_mask, "ALT"]
        variants_plink_df.loc[swap_mask, "ALT"] = temp
        del temp

        # write to a true VCF now
        adata_for_vcf = adata if save_vcf_samples else None
        write_to_vcf(variants_plink_df, output_vcf=geuvadis_true_vcf, save_vcf_samples=save_vcf_samples, adata=adata_for_vcf)

#* Run vk ref
vk.ref(
    variants=variants,
    sequences=sequences,
    w=w,
    k=k,
    seq_id_column=seq_id_column,
    var_column=var_column,
    out=vk_ref_out,
    dlist_reference_source=dlist_reference_source,
    index_out=vcrs_index,
    t2g_out=vcrs_t2g,
    threads=threads,
    chunksize=chunksize,
    save_logs=True,
    verbose=True,
)