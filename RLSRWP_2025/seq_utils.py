"""varseek sequencing utilities."""
import os
import scanpy as sc
import scipy.sparse as sp
import networkx as nx

import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm
import re
import scanpy as sc
from sklearn.metrics import euclidean_distances
from sklearn.neighbors import NearestNeighbors

from RLSRWP_2025.visualization_utils import (
    plot_ascending_bar_plot_of_cluster_distances,
    plot_jaccard_bar_plot,
    plot_knn_tissue_frequencies,
)

from varseek.utils import (
    calculate_metrics,
    calculate_sensitivity_specificity,
    create_venn_diagram,
    draw_confusion_matrix,
    safe_literal_eval,
    vcf_to_dataframe,
)

# tqdm.pandas()

def compute_intersection(set1, set2):
    return len(set1.intersection(set2))


def compute_jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union


# * TODO: this considers only the largest number of nearest neighbors, but not the ordinality, which is especially important because if I look for 15 nearest neighbors but a cancer type of interest only has 3 samples in my dataset then it will have a hard time winning out
def compute_knn(adata_unknown, adata_combined_ccle_rnaseq, out_dir=".", k=10):
    # Combine PCA embeddings from both datasets
    combined_pca_embeddings = np.vstack([adata_combined_ccle_rnaseq.obsm["X_pca"], adata_unknown.obsm["X_pca"]])

    # Fit the k-NN model on the reference dataset
    reference_embeddings = combined_pca_embeddings.obsm["X_pca"]
    knn_model = NearestNeighbors(n_neighbors=k, metric="euclidean")
    knn_model.fit(reference_embeddings)

    # Query the nearest neighbors for each point in `adata_unknown`
    unknown_embeddings = adata_unknown.obsm["X_pca"]
    distances, indices = knn_model.kneighbors(unknown_embeddings)

    tissue_counts = plot_knn_tissue_frequencies(
        indices,
        adata_combined_ccle_rnaseq,
        output_plot_file=f"{out_dir}/knn_tissue_frequencies.png",
    )

    plurality_tissue_counts = np.argmax(tissue_counts, axis=1)
    print("True cancer type: ", adata_unknown.uns["sample_name"])
    print("Predicted cancer type: ", plurality_tissue_counts)
    print(
        "Predicted cancer type and true cancer type match: ",
        adata_unknown.uns["sample_name"] == plurality_tissue_counts,
    )


def compute_cluster_centroid_distances(adata_unknown, adata_combined_ccle_rnaseq, out_dir="."):
    # Ensure PCA has been performed on `adata_combined_ccle_rnaseq` and `adata_unknown`
    # Get the PCA coordinates for both datasets
    reference_pca = adata_combined_ccle_rnaseq.obsm["X_pca"]
    unknown_pca = adata_unknown.obsm["X_pca"]

    # Calculate the centroids of each Leiden cluster
    cluster_centroids = {}
    for cluster in adata_combined_ccle_rnaseq.obs["tissue"].unique():
        # Select cells in the current cluster and compute the centroid
        cluster_cells = reference_pca[adata_combined_ccle_rnaseq.obs["tissue"] == cluster]
        centroid = cluster_cells.mean(axis=0)
        cluster_centroids[cluster] = centroid

    # Convert centroids to an array for easy distance calculation
    centroid_matrix = np.array(list(cluster_centroids.values()))

    # Calculate distances from the unknown point to each centroid
    distances = euclidean_distances(unknown_pca, centroid_matrix)

    sorted_distances = sorted(zip(cluster_centroids.keys(), distances[0]), key=lambda x: x[1])

    cluster_centroid_distance_file = f"{out_dir}/cluster_centroid_distances.txt"
    # Write the sorted distances to a file
    with open(cluster_centroid_distance_file, "w", encoding="utf-8") as f:
        for cluster, dist in sorted_distances:
            f.write(f"Distance from unknown sample to centroid of cluster {cluster}: {dist}\n")

    cluster_centroid_distance_plot_file = f"{out_dir}/cluster_centroid_distances.png"
    plot_ascending_bar_plot_of_cluster_distances(sorted_distances, output_plot_file=cluster_centroid_distance_plot_file)

    # Output the closest cluster
    closest_cluster, closest_distance = sorted_distances[0]
    print("True cancer type: ", adata_unknown.uns["sample_name"])
    print("Predicted cancer type and distance: ", closest_cluster, closest_distance)
    print(
        "Predicted cancer type and true cancer type match: ",
        adata_unknown.uns["sample_name"] == closest_cluster,
    )


def compute_jaccard_indices(
    adata_unknown,
    grouped_tissue_dir,
    jaccard_type="mutation",
    out_dir=".",
    jaccard_or_intersection="jaccard",
):
    if jaccard_type == "mutation":
        column_name = "header_with_gene_name"
        text_file = "sorted_mutations.txt"
    elif jaccard_type == "mutated_genes":
        column_name = "gene_name"
        text_file = "sorted_mutated_genes.txt"
    else:
        raise ValueError("Invalid jaccard_type. Must be 'mutation' or 'mutated_genes'")
    # Initialize a dictionary to store mutations for each tissue
    tissue_mutations = {}

    adata_unknown_unique_mutations = set(adata_unknown.var[column_name])

    # Loop through each subfolder in `grouped_tissue_dir`
    for tissue in os.listdir(grouped_tissue_dir):
        specific_tissue_dir = os.path.join(grouped_tissue_dir, tissue)
        sorted_mutations_file = os.path.join(specific_tissue_dir, text_file)

        # Initialize an empty set for mutations in this tissue
        mutations = set()

        # Check if the file exists in the current tissue folder
        if os.path.isfile(sorted_mutations_file):
            with open(sorted_mutations_file, "r", encoding="utf-8") as f:
                for line in f:
                    mutation_name = line.split()[0]  # Grab the first column (mutation name)
                    mutations.add(mutation_name)

        if jaccard_or_intersection == "jaccard":
            jaccard_index_mutations = compute_jaccard_index(mutations, adata_unknown_unique_mutations)
        elif jaccard_or_intersection == "intersection":
            jaccard_index_mutations = compute_intersection(mutations, adata_unknown_unique_mutations)

        tissue_mutations[tissue] = jaccard_index_mutations

    tissues = list(adata_unknown_unique_mutations.keys())
    jaccard_values = list(adata_unknown_unique_mutations.values())

    sorted_data = sorted(zip(jaccard_values, tissues), reverse=True)
    tissues, jaccard_values = zip(*sorted_data)

    plot_jaccard_bar_plot(
        tissues,
        jaccard_values,
        output_plot_file=f"{out_dir}/jaccard_bar_plot_{jaccard_type}.png",
    )

    # Find the tissue with the maximum Jaccard index
    max_tissue = max(adata_unknown_unique_mutations, key=adata_unknown_unique_mutations.get)
    max_jaccard = adata_unknown_unique_mutations[max_tissue]

    print("True cancer type: ", adata_unknown.uns["sample_name"])
    print("Predicted cancer type and jaccard for mutations: ", max_tissue, max_jaccard)
    print(
        "Predicted cancer type and true cancer type match: ",
        adata_unknown.uns["sample_name"] == max_tissue,
    )


# TODO: write this
def compute_mx_metric():
    pass
    # return a list of marker mutations for each cancer type, and provide these as input along with the unknown sample's mutation matrix to *mx assign* (might only work for single-cell)
    # - identify cancer-specific marker mutations (or genes, and just apply to all of it's mutations in my database) from the literature ("old")
    # - return a list of marker mutations for each cancer type from my cell lines ("new")
    # - merge these with *ec merge*
    # - run *ec index* on the merged marker list
    # - run *mx extract* on the output of *ec index*
    # - run *mx clean* on the output of *mx extract*
    # - run *mx assign* on the output of *mx clean*


def predict_cancer_type(
    adata_path,
    adata_combined_ccle_rnaseq,
    pca_components,
    grouped_tissue_dir,
    sample_name,
    sample_to_cancer_type,
    metric="knn",
    k=10,
    use_binary_matrix=False,
):
    sample_dir = os.path.dirname(os.path.dirname(adata_path))
    adata_unknown = sc.read_h5ad(adata_path)
    adata_unknown.uns["sample_name"] = sample_name
    adata_unknown.uns["cancer_type"] = sample_to_cancer_type[sample_name]

    if use_binary_matrix:
        adata_unknown.X = (adata_unknown.X > 0).astype(int)

    # Identify common genes between the two datasets
    common_genes = adata_unknown.var_names.intersection(adata_combined_ccle_rnaseq.var_names)

    # Subset adata_unknown to include only these common genes
    adata_unknown = adata_unknown[:, common_genes].copy()

    # Step 2: Center the new data using the mean of the original data
    # Note: adata_original.X.mean(axis=0) assumes both datasets have the same set of genes/variables
    mean_center = np.mean(adata_combined_ccle_rnaseq.X, axis=0)
    adata_unknown_centered = adata_unknown.X - mean_center

    # Step 3: Project the new data into the PCA space of the original data
    # This will give you the PCA coordinates for adata_unknown in the same space as adata_original
    adata_unknown.obsm["X_pca"] = adata_unknown_centered.dot(pca_components)

    # Optional: Store the explained variance for reference if needed
    adata_unknown.uns["pca"] = adata_combined_ccle_rnaseq.uns["pca"]  # Copy explained variance, etc.

    if metric == "knn":
        compute_knn(
            adata_unknown=adata_unknown,
            adata_combined_ccle_rnaseq=adata_combined_ccle_rnaseq,
            out_dir=sample_dir,
            k=k,
        )
    elif metric == "euclidean":
        compute_cluster_centroid_distances(
            adata_unknown=adata_unknown,
            adata_combined_ccle_rnaseq=adata_combined_ccle_rnaseq,
            out_dir=sample_dir,
        )
    elif metric == "jaccard":
        compute_jaccard_indices(
            adata_unknown=adata_unknown,
            grouped_tissue_dir=grouped_tissue_dir,
            out_dir=sample_dir,
        )
    elif metric == "intersection":
        compute_jaccard_indices(
            adata_unknown=adata_unknown,
            grouped_tissue_dir=grouped_tissue_dir,
            out_dir=sample_dir,
        )
    elif metric == "mx":
        compute_mx_metric()




def add_variant_type_gatk(df):
    # Define the conditions
    conditions = [(df["REF"].str.len() > 1) & (df["ALT"].str.len() == 1), (df["REF"].str.len() == 1) & (df["ALT"].str.len() > 1), (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1), (df["REF"].str.len() > 1) & (df["ALT"].str.len() > 1)]  # Deletion  # Insertion  # Substitution  # Delins

    # Define the corresponding mutation types
    variant_types = ["deletion", "insertion", "substitution", "delins"]

    # Apply the conditions and assign the values to the new column
    df["variant_type"] = np.select(conditions, variant_types, default="unknown")

    # For 'deletion', add 'deleted_bases' column with REF[1:]
    df.loc[df["variant_type"] == "deletion", "deleted_bases"] = df["REF"].str[1:]

    # For 'insertion', add 'inserted_bases' column with ALT[1:]
    df.loc[df["variant_type"] == "insertion", "inserted_bases"] = df["ALT"].str[1:]

    return df

def merge_gatk_and_cosmic(df_mut, cosmic_df, exact_position=False):
    df_mut = df_mut.copy()
    cosmic_df = cosmic_df.copy()

    if exact_position:
        # take the intersection of COSMIC and STAR dfs based on CHROM, POS, REF, ALT - but keep the ID from the COSMIC vcf
        # note that because there may be 2+ rows in COSMIC that share the same chrom, pos, ref, and alt (ie from alternatively spliced mutations), there may be some rows in mut_cosmic_merged_df that are duplicates for all but the ID column - but I don't want to drop these because I want to make sure I consider both headers in my set
        merged_df = pd.merge(df_mut, cosmic_df, on=["CHROM", "POS", "REF", "ALT"], how="inner", suffixes=("_df1", "_df2"))

        merged_df = merged_df.drop(columns=["ID_df1", "POS_df1"]).rename(columns={"ID_df2": "ID", "POS_df2": "POS"})
    else:
        if "variant_type" not in df_mut.columns:
            df_mut = add_variant_type_gatk(df_mut)
        if "variant_type" not in cosmic_df.columns:
            cosmic_df = add_variant_type_gatk(cosmic_df)

        # Split `df_mut` and `cosmic_df` by mutation type
        sub_delins_mut = df_mut[df_mut["variant_type"].isin(["substitution", "delins"])]
        del_mut = df_mut[df_mut["variant_type"] == "deletion"]
        ins_mut = df_mut[df_mut["variant_type"] == "insertion"]

        sub_delins_cosmic = cosmic_df[cosmic_df["variant_type"].isin(["substitution", "delins"])]
        del_cosmic = cosmic_df[cosmic_df["variant_type"] == "deletion"]
        ins_cosmic = cosmic_df[cosmic_df["variant_type"] == "insertion"]

        # 1. Merge substitution and delins
        sub_delins_merged = pd.merge(sub_delins_mut, sub_delins_cosmic, on=["CHROM", "POS", "REF", "ALT"], how="left", suffixes=("_df1", "_df2"))

        sub_delins_merged = sub_delins_merged.drop(columns=["ID_df1", "variant_type_df1", "variant_type_df2", "deleted_bases_df1", "deleted_bases_df2", "inserted_bases_df1", "inserted_bases_df2"]).rename(columns={"ID_df2": "ID"})

        # 2. Merge deletion
        del_merged = pd.merge(del_mut, del_cosmic, on=["CHROM", "deleted_bases"], how="left", suffixes=("_df1", "_df2"))

        # Filter rows where POS is within a certain range, setting the threshold dynamically based on 'sub' column - subs must be perfect, indels can be within 5
        del_merged = del_merged[abs(del_merged["POS_df1"] - del_merged["POS_df2"]) <= 5]

        del_merged = del_merged.drop(columns=["ID_df1", "POS_df1", "REF_df1", "ALT_df1", "variant_type_df1", "variant_type_df2", "inserted_bases_df1", "inserted_bases_df2", "deleted_bases"]).rename(columns={"ID_df2": "ID", "POS_df2": "POS", "REF_df2": "REF", "ALT_df2": "ALT"})

        # 3. Merge insertion
        ins_merged = pd.merge(ins_mut, ins_cosmic, on=["CHROM", "inserted_bases"], how="left", suffixes=("_df1", "_df2"))

        # Filter rows where POS is within a certain range, setting the threshold dynamically based on 'sub' column - subs must be perfect, indels can be within 5
        ins_merged = ins_merged[abs(ins_merged["POS_df1"] - ins_merged["POS_df2"]) <= 5]

        ins_merged = ins_merged.drop(columns=["ID_df1", "POS_df1", "REF_df1", "ALT_df1", "variant_type_df1", "variant_type_df2", "deleted_bases_df1", "deleted_bases_df2", "inserted_bases"]).rename(columns={"ID_df2": "ID", "POS_df2": "POS", "REF_df2": "REF", "ALT_df2": "ALT"})

        # Combine all results
        merged_df = pd.concat([sub_delins_merged, del_merged, ins_merged], ignore_index=True)

        return merged_df


def create_mutated_gene_count_matrix_from_mutation_count_matrix(adata, sum_strategy="total_reads", merge_strategy="all", use_binary_matrix=False):
    """
    This function takes a mutation count matrix and aggregates the counts for mutations belonging to the same gene. The function assumes that the AnnData object has the following columns in adata.var:
    - gene_name_set_string: a string containing a semi-colon separated list of gene names for each mutation
    - vcrs_id: a unique identifier for each mutation

    Parameters
    ----------
    adata : AnnData

    merge_strategy : str
        The strategy to use when merging mutations. The following options are available:
        - 'all': merge based on all genes matching (i.e., gene_name_set_string)
        - 'any': merge based on any genes mapping (i.e., any match in gene_name_set)
    sum_strategy: str
        The strategy for summing VCRSs - options:
        - 'total_reads': sum the total reads for each VCRS
        - 'unique_variants': sum the number of unique variants detected for a gene
    """

    if sum_strategy == "unique_variants":
        adata.X = (adata.X > 0).astype(int)  # convert to binary matrix
        count_column = "variant_count"
    else:
        count_column = "vcrs_count"

    if merge_strategy == "all":
        gene_column = "gene_name_set_string"
    elif merge_strategy == "any":  # TODO: untested for merge_strategy == "any"
        gene_column = "gene_name_set"
        gene_names = adata.var[gene_column]
        vcrs_ids = adata.var_names
        # Create a graph where each node is an vcrs_id
        graph = nx.Graph()
        for i, genes in enumerate(gene_names):
            for j in range(i + 1, len(gene_names)):
                if set(genes).intersection(gene_names[j]):
                    graph.add_edge(vcrs_ids[i], vcrs_ids[j])

        # Find connected components (each component is a group of columns to merge)
        components = list(nx.connected_components(graph))

        # Step 2: Create a mapping for new groups
        new_var = []
        group_mapping = {}
        for group_id, component in enumerate(components):
            # Combine gene names and vcrs_ids for the group
            group_genes = sorted(set.union(*(set(gene_names[vcrs_ids.tolist().index(mcrs)]) for mcrs in component)))
            group_vcrs_ids = sorted(component)

            # Use a representative name for the group
            group_name = ";".join(group_genes)
            for mcrs in component:
                group_mapping[mcrs] = group_name

            # Store new metadata
            new_var.append({"gene_name_set_string": group_name, "vcrs_id_list": group_vcrs_ids})

    # Step 1: Extract mutation-gene mappings
    gene_mapping = adata.var[gene_column]  # because I am using gene_name_set_string, this means that any merged mcrs's with different gene names will not be included in merging/summing
    vcrs_id_mapping = adata.var["vcrs_id"]

    # Step 2: Convert your data to a DataFrame for easier manipulation
    if sp.issparse(adata.X):
        data_df = pd.DataFrame.sparse.from_spmatrix(adata.X, index=adata.obs_names, columns=adata.var_names)
    else:
        data_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    # Step 3: Add gene mapping to the DataFrame for aggregation
    if merge_strategy == "all":
        data_df.columns = gene_mapping.values
    elif merge_strategy == "any":
        data_df.columns = [group_mapping[col] for col in data_df.columns]

    vcrs_id_df = pd.Series(vcrs_id_mapping.values, index=adata.var_names).groupby(gene_mapping).agg(list)

    # Step 4: Group by gene and sum across mutations belonging to the same gene
    data_gene_df = data_df.groupby(axis=1, level=0).sum()

    # Step 5: Convert the result back into an AnnData object
    adata_gene = sc.AnnData(data_gene_df, obs=adata.obs.copy())
    adata_gene.var_names = data_gene_df.columns  # Gene names

    adata_gene.var[gene_column] = adata_gene.var_names  # make this a column
    adata_gene.var["vcrs_id_list"] = vcrs_id_df.loc[data_gene_df.columns].values

    if use_binary_matrix:
        adata_gene.X = (adata_gene.X > 0).astype(int)

    adata_gene.var[count_column] = adata_gene.X.sum(axis=0).A1 if hasattr(adata_gene.X, "A1") else np.asarray(adata_gene.X.sum(axis=0)).flatten()

    return adata_gene


def perform_analysis(vcf_file, unique_mcrs_df_path, cosmic_df, plot_output_folder = "plots", unique_mcrs_df_path_out = None, package_name = "tool"):
    if not unique_mcrs_df_path_out:
        unique_mcrs_df_path_out = unique_mcrs_df_path

    #* Merging into COSMIC
    # Convert VCF to DataFrame
    df_tool = vcf_to_dataframe(vcf_file, additional_columns = True)
    df_tool = df_tool[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO_DP']].rename(columns={'INFO_DP': f'DP_{package_name}'})

    # load in unique_mcrs_df
    unique_mcrs_df = pd.read_csv(unique_mcrs_df_path)
    if "header_list" in unique_mcrs_df.columns:
        unique_mcrs_df["header_list"] = unique_mcrs_df["header_list"].apply(safe_literal_eval)
    else:
        unique_mcrs_df["header_list"] = unique_mcrs_df["vcrs_header"].str.split(";")

    tool_cosmic_merged_df = merge_gatk_and_cosmic(df_tool, cosmic_df, exact_position=False)  # change exact_position to True to merge based on exact position as before
    id_set_tool = set(tool_cosmic_merged_df['ID'])

    # Merge DP values into unique_mcrs_df
    # Step 1: Remove rows with NaN values in 'ID' column
    tool_cosmic_merged_df_for_merging = tool_cosmic_merged_df[['ID', f'DP_{package_name}']].dropna(subset=['ID']).rename(columns={'ID': 'vcrs_header'})

    # Step 2: Drop duplicates from 'ID' column
    tool_cosmic_merged_df_for_merging = tool_cosmic_merged_df_for_merging.drop_duplicates(subset=['vcrs_header'])

    # Step 3: Left merge with unique_mcrs_df
    unique_mcrs_df = pd.merge(
        unique_mcrs_df,               # Left DataFrame
        tool_cosmic_merged_df_for_merging,         # Right DataFrame
        on='vcrs_header',
        how='left'
    )

    number_of_mutations_tool = len(df_tool.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']))
    number_of_cosmic_mutations_tool = len(tool_cosmic_merged_df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']))

    # unique_mcrs_df['header_list'] each contains a list of strings. I would like to make a new column unique_mcrs_df[f'mutation_detected_{package_name}'] where each row is True if any value from the list unique_mcrs_df['vcrs_header'] is in the set id_set_mut  # keep in mind that my IDs are the mutation headers (ENST...), NOT mcrs headers or mcrs ids
    unique_mcrs_df[f'mutation_detected_{package_name}'] = unique_mcrs_df['header_list'].apply(
        lambda header_list: any(header in id_set_tool for header in header_list)
    )

    # calculate expression error
    unique_mcrs_df['mutation_expression_prediction_error'] = unique_mcrs_df[f'DP_{package_name}'] - unique_mcrs_df['number_of_reads_mutant']  # positive means overpredicted, negative means underpredicted

    unique_mcrs_df['TP'] = (unique_mcrs_df['included_in_synthetic_reads_mutant'] & unique_mcrs_df[f'mutation_detected_{package_name}'])
    unique_mcrs_df['FP'] = (~unique_mcrs_df['included_in_synthetic_reads_mutant'] & unique_mcrs_df[f'mutation_detected_{package_name}'])
    unique_mcrs_df['FN'] = (unique_mcrs_df['included_in_synthetic_reads_mutant'] & ~unique_mcrs_df[f'mutation_detected_{package_name}'])
    unique_mcrs_df['TN'] = (~unique_mcrs_df['included_in_synthetic_reads_mutant'] & ~unique_mcrs_df[f'mutation_detected_{package_name}'])

    tool_stat_path = f"{plot_output_folder}/reference_metrics_{package_name}.txt"
    metric_dictionary_reference = calculate_metrics(unique_mcrs_df, header_name = "vcrs_header", check_assertions = False, out = tool_stat_path)
    draw_confusion_matrix(metric_dictionary_reference)

    true_set = set(unique_mcrs_df.loc[unique_mcrs_df['included_in_synthetic_reads_mutant'], 'vcrs_header'])
    positive_set = set(unique_mcrs_df.loc[unique_mcrs_df[f'mutation_detected_{package_name}'], 'vcrs_header'])
    create_venn_diagram(true_set, positive_set, TN = metric_dictionary_reference['TN'], out_path = f"{plot_output_folder}/venn_diagram_reference_cosmic_only_{package_name}.png")

    noncosmic_mutation_id_set = {f'{package_name}_fp_{i}' for i in range(1, number_of_mutations_tool - number_of_cosmic_mutations_tool + 1)}

    positive_set_including_noncosmic_mutations = positive_set.union(noncosmic_mutation_id_set)
    false_positive_set = set(unique_mcrs_df.loc[unique_mcrs_df['FP'], 'vcrs_header'])
    false_positive_set_including_noncosmic_mutations = false_positive_set.union(noncosmic_mutation_id_set)

    FP_including_noncosmic = len(false_positive_set_including_noncosmic_mutations)
    accuracy, sensitivity, specificity = calculate_sensitivity_specificity(metric_dictionary_reference['TP'], metric_dictionary_reference['TN'], FP_including_noncosmic, metric_dictionary_reference['FN'])

    with open(tool_stat_path, "a", encoding="utf-8") as file:
        file.write(f"FP including non-cosmic: {FP_including_noncosmic}\n")
        file.write(f"accuracy including non-cosmic: {accuracy}\n")
        file.write(f"specificity including non-cosmic: {specificity}\n")

    create_venn_diagram(true_set, positive_set_including_noncosmic_mutations, TN = metric_dictionary_reference['TN'], out_path = f"{plot_output_folder}/venn_diagram_reference_including_noncosmics_{package_name}.png")

    unique_mcrs_df.rename(columns={'TP': f'TP_{package_name}', 'FP': f'FP_{package_name}', 'TN': f'TN_{package_name}', 'FN': f'FN_{package_name}', 'mutation_expression_prediction_error': f'mutation_expression_prediction_error_{package_name}'}, inplace=True)
    unique_mcrs_df.to_csv(unique_mcrs_df_path_out, index=False)


def make_transcript_df_from_gtf(gtf):
    if gtf is None:
        raise ValueError("gtf must be provided if variant_position_annotations is not 'cdna' or strand_bias_end is '3p'")
    if isinstance(gtf, str):
        gtf_cols = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
        gtf_df = pd.read_csv(gtf, sep="\t", comment="#", names=gtf_cols)
    elif isinstance(gtf, pd.DataFrame):
        gtf_df = gtf.copy()
    else:
        raise ValueError("gtf must be a path to a GTF file or a pandas DataFrame")
    gtf_df["region_length"] = gtf_df["end"] - gtf_df["start"] + 1  # this corresponds to the unspliced transcript length, NOT the spliced transcript/CDS length (which is what I care about)
    gtf_df["transcript_ID"] = gtf_df["attributes"].str.extract(r'transcript_id "([^"]+)"')
    transcript_df = gtf_df[gtf_df["feature"] == "transcript"].copy()
    transcript_df = transcript_df.drop_duplicates(subset="transcript_ID", keep="first")

    five_prime_lengths = []
    three_prime_lengths = []
    transcript_lengths = []
    for transcript_id in transcript_df["transcript_ID"]:
        specific_transcript_df = gtf_df.loc[gtf_df["transcript_ID"] == transcript_id].copy()
        five_sum = specific_transcript_df[(specific_transcript_df["feature"] == "five_prime_utr")]["region_length"].sum()
        three_sum = specific_transcript_df[(specific_transcript_df["feature"] == "three_prime_utr")]["region_length"].sum()
        exon_sum = specific_transcript_df[(specific_transcript_df["feature"] == "exon")]["region_length"].sum()
        transcript_length = five_sum + three_sum + exon_sum

        five_prime_lengths.append(five_sum)
        three_prime_lengths.append(three_sum)
        transcript_lengths.append(transcript_length)
    transcript_df["five_prime_length"] = five_prime_lengths
    transcript_df["three_prime_length"] = three_prime_lengths
    transcript_df["transcript_length"] = transcript_lengths

    transcript_df["utr_length_preceding_transcript"] = transcript_df.apply(
        lambda row: row["three_prime_utr_length"] if row["strand"] == "-" else row["five_prime_utr_length"],
        axis=1
    )

    return gtf_df, transcript_df


def convert_vcf_samples_to_anndata(vcf_path, adata_out=None, sample_titles_set=None, total=None):
    from cyvcf2 import VCF

    vcf = VCF(vcf_path)

    sample_names = vcf.samples
    n_samples = len(sample_names)

    if sample_titles_set:
        sample_indices_to_keep_set = {i for i, name in enumerate(sample_names) if name in sample_titles_set}
    else:
        sample_indices_to_keep_set = set(range(n_samples))

    # Store sparse matrix data
    data = []
    rows = []
    cols = []
    variant_ids = []
    has_id = []

    for var_idx, variant in tqdm(enumerate(vcf), desc="Processing VCF", unit="variants", total=total):
        genotypes = variant.genotypes  # List of [GT1, GT2, phased, ...]
        if not genotypes or len(genotypes) != n_samples:
            continue  # Skip malformed rows

        sample_idx_for_row = 0
        for sample_idx, gt in enumerate(genotypes):
            if sample_idx not in sample_indices_to_keep_set:
                continue
            sample_idx_for_row += 1
            
            if gt[0] == -1 or gt[1] == -1:
                continue  # missing
            
            gt_sum = gt[0] + gt[1]
            rows.append(sample_idx)
            cols.append(var_idx)
            data.append(gt_sum)

        if variant.ID:
            variant_ids.append(variant.ID)
            has_id.append(True)
        else:
            variant_ids.append(f"temp_{var_idx}")
            has_id.append(False)

    # Create sparse matrix (samples x variants)
    X = sp.csr_matrix((data, (rows, cols)), shape=(n_samples, len(variant_ids)))

    # Build obs and var
    obs = pd.DataFrame(index=sample_names)
    var = pd.DataFrame(index=variant_ids)
    var["has_id"] = has_id

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var)
    if adata_out is None:
        adata_out = re.sub(r'\.vcf(?:\.gz)?$', '.h5ad', vcf_path)
    adata.write_h5ad(adata_out)
    return adata
