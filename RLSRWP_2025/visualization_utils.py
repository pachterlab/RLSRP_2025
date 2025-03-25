import json
import os
from collections import Counter, defaultdict

import scanpy as sc

import re
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import seaborn as sns
from rich.console import Console
from matplotlib.ticker import FuncFormatter, LogLocator
from scipy import stats
from scipy.stats import t, ttest_rel
from statsmodels.stats.contingency_tables import mcnemar

console = Console()

# Set global settings
plt.rcParams.update(
    {
        "savefig.dpi": 450,  # Set resolution to 450 dpi
        "font.family": "DejaVu Sans",  # Set font to Arial  # TODO: replace with Arial for Nature
        "pdf.fonttype": 42,  # Embed fonts as TrueType (keeps text editable)
        "ps.fonttype": 42,  # Same for PostScript files
        "savefig.format": "pdf",  # Default save format as PNG
        "savefig.bbox": "tight",  # Adjust bounding box to fit tightly
        "figure.facecolor": "white",  # Set figure background to white (common for RGB)
        "savefig.transparent": False,  # Disable transparency
    }
)

color_map_10 = plt.get_cmap("tab10").colors  # Default color map with 10 colors

color_map_20_original = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"]  # plotly category 20

color_map_20 = ["#f08925", "#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"]  # modified to swap 1 and 2 (orange first), and replaced the orange with varseek orange

SAVE_PDF_GLOBAL = (os.getenv("VARSEEK_SAVE_PDF") == "TRUE")
DPI = 450

def plot_scree(adata, output_plot_file=None):
    variance_explained = adata.uns["pca"]["variance_ratio"]
    num_components = len(variance_explained)

    # Plot the scree plot
    plt.figure(figsize=(10, 5))
    plt.plot(
        np.arange(1, len(variance_explained) + 1),
        variance_explained,
        marker="o",
        linestyle="-",
    )
    plt.xticks(ticks=np.arange(1, num_components + 1))
    plt.xlabel("Principal Component")
    plt.ylabel("Variance Explained")
    plt.title("Scree Plot")
    if output_plot_file:
        os.makedirs(os.path.dirname(output_plot_file), exist_ok=True)
        plt.savefig(output_plot_file, format="png", dpi=DPI)
        if SAVE_PDF_GLOBAL:
            plt.savefig(output_plot_file.replace(".png", ".pdf"), format="pdf", dpi=DPI)

    plt.show()
    plt.close()


def plot_loading_contributions(adata, PC_index=0, top_genes_stats=100, top_genes_plot=10, output_stats_file=None, output_plot_file=None, show=False):
    # Get PCA loadings for the selected component
    loadings = adata.varm["PCs"][:, PC_index]

    # Find indices of top genes by absolute loading values
    top_gene_indices_stats = np.argsort(np.abs(loadings))[::-1][:top_genes_stats]
    top_gene_names_stats = adata.var_names[top_gene_indices_stats]
    top_gene_loadings_stats = loadings[top_gene_indices_stats]

    if output_stats_file:
        os.makedirs(os.path.dirname(output_stats_file), exist_ok=True)
        with open(output_stats_file, "w", encoding="utf-8") as f:
            for gene, loading in zip(top_gene_names_stats, top_gene_loadings_stats):
                f.write(f"{gene} {loading}\n")

    top_gene_indices_plot = np.argsort(np.abs(loadings))[::-1][:top_genes_plot]
    top_gene_names_plot = adata.var_names[top_gene_indices_plot]
    top_gene_loadings_plot = loadings[top_gene_indices_plot]

    # Plot as a horizontal bar chart
    plt.figure(figsize=(8, 6))
    plt.barh(top_gene_names_plot, top_gene_loadings_plot, color="skyblue")
    plt.xlabel("Contribution to PC1")
    plt.ylabel("Gene")
    plt.title("Top Gene Contributions to PC1")
    plt.gca().invert_yaxis()  # Invert Y-axis for descending order
    if output_plot_file:
        os.makedirs(os.path.dirname(output_plot_file), exist_ok=True)
        plt.savefig(output_plot_file, format="png", dpi=DPI)
        if SAVE_PDF_GLOBAL:
            plt.savefig(output_plot_file.replace(".png", ".pdf"), format="pdf", dpi=DPI)
    if show:
        plt.show()
    plt.close()


def find_resolution_for_target_clusters(adata, target_clusters, tolerance=3, max_iters=10):
    # Initial bounds for resolution
    lower, upper = 0.1, 10
    if max_iters <= 0:
        raise ValueError("max_iters must be a positive integer")
    for i in range(max_iters):
        # Take the midpoint as the next test resolution
        adata_copy = adata.copy()
        resolution = (lower + upper) / 2
        sc.tl.leiden(adata_copy, resolution=resolution)
        num_clusters = adata_copy.obs["leiden"].nunique()

        print(f"Iteration {i + 1}: Resolution = {resolution}, Clusters = {num_clusters}")

        # Check if the number of clusters is within the tolerance of the target
        if abs(num_clusters - target_clusters) <= tolerance:
            return adata_copy, resolution, num_clusters

        # Update bounds based on whether we have too many or too few clusters
        if num_clusters < target_clusters:
            lower = resolution
        else:
            upper = resolution

    return (
        adata_copy,
        resolution,
        num_clusters,
    )  # Return last tested resolution if exact match not found


def plot_contingency_table(adata, column1="tissue", column2="leiden", output_plot_file=None):
    # Create a contingency table (counts of cells in each combination of tissue and leiden cluster)
    contingency_table = pd.crosstab(adata.obs[column1], adata.obs[column2])

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(contingency_table, annot=True, cmap="Blues", fmt="d")
    plt.xlabel("Leiden Cluster")
    plt.ylabel("Tissue Type")
    plt.title(f"Heatmap of Agreement Between {column1} and {column2}")
    if output_plot_file:
        plt.savefig(output_plot_file)
    plt.show()
    plt.close()


def plot_knn_tissue_frequencies(indices, adata_combined_ccle_rnaseq, output_plot_file=None):
    # Split the 'obs_names' to extract the tissue component
    neighbor_tissues = [adata_combined_ccle_rnaseq.obs_names[idx].split("_")[-1] for neighbors in indices for idx in neighbors]  # Extract tissue from "experiment_tissue"

    # Count occurrences of each tissue in the nearest neighbors
    tissue_counts = Counter(neighbor_tissues)

    # Plot the tissue frequencies
    plt.figure(figsize=(10, 6))
    plt.bar(tissue_counts.keys(), tissue_counts.values(), color="skyblue")
    plt.xlabel("Tissue Type")
    plt.ylabel("Frequency")
    plt.title("Frequency of Each Tissue in Nearest Neighbors")
    plt.xticks(rotation=45)
    if output_plot_file:
        plt.savefig(output_plot_file, format="png", dpi=DPI)
        if SAVE_PDF_GLOBAL:
            plt.savefig(output_plot_file.replace(".png", ".pdf"), format="pdf", dpi=DPI)
    plt.show()
    plt.close()

    return tissue_counts


def plot_ascending_bar_plot_of_cluster_distances(sorted_distances, output_plot_file=None):
    # Separate clusters and distances for plotting
    clusters_sorted, distances_sorted = zip(*sorted_distances)

    plt.figure(figsize=(10, 6))
    plt.bar(clusters_sorted, distances_sorted, color="skyblue")
    plt.xlabel("Cluster")
    plt.ylabel("Distance to Unknown Sample")
    plt.title("Distance from Unknown Sample to Each Cluster Centroid (Ascending Order)")
    plt.xticks(rotation=45)
    if output_plot_file:
        plt.savefig(output_plot_file, format="png", dpi=DPI)
        if SAVE_PDF_GLOBAL:
            plt.savefig(output_plot_file.replace(".png", ".pdf"), format="pdf", dpi=DPI)

    plt.show()
    plt.close()


def plot_jaccard_bar_plot(tissues, jaccard_values, output_plot_file=None):
    sorted_data = sorted(zip(jaccard_values, tissues), reverse=True)
    sorted_jaccard_values, sorted_tissues = zip(*sorted_data)

    plt.figure(figsize=(10, 6))
    plt.bar(sorted_tissues, sorted_jaccard_values, color="skyblue")
    plt.xlabel("Tissue")
    plt.ylabel("Jaccard Index")
    plt.title("Jaccard Index for Each Tissue")
    plt.xticks(rotation=45)
    if output_plot_file:
        plt.savefig(output_plot_file, format="png", dpi=DPI)
        if SAVE_PDF_GLOBAL:
            plt.savefig(output_plot_file.replace(".png", ".pdf"), format="pdf", dpi=DPI)

    plt.show()
    plt.close()




def plot_overall_metrics(metric_dict_collection, primary_metrics=("accuracy", "sensitivity", "specificity"), display_numbers=False, unique_mcrs_df=None, show_p_values=False, bonferroni=True, output_file=None, show=True, output_file_p_values=None, filter_real_negatives=False):
    if not isinstance(primary_metrics, (str, list, tuple)):
        raise ValueError("Primary metrics must be a string, list, or tuple.")

    if unique_mcrs_df is not None:
        unique_mcrs_df = unique_mcrs_df.copy()

    if not filter_real_negatives and "expression_error" in primary_metrics:
        print("Warning: filtering real negatives is recommended when using expression error as a primary metric, but this setting is not currently enabled. Recommended: filter_real_negatives = True.")

    if filter_real_negatives:
        unique_mcrs_df = unique_mcrs_df[unique_mcrs_df["included_in_synthetic_reads_mutant"] == True]

    if isinstance(primary_metrics, str):
        primary_metrics = [primary_metrics]
    elif isinstance(primary_metrics, tuple):
        primary_metrics = list(primary_metrics)

    # Extract keys and values
    groups = list(metric_dict_collection.keys())  # Outer keys: varseek, mutect2, haplotypecaller
    colors = color_map_20[: len(groups)]

    # Prepare data
    x_primary = np.arange(len(primary_metrics))  # Positions for the metrics on the x-axis
    bar_width = 0.25  # Width of each primary_metrics
    offsets = np.arange(len(groups)) * bar_width - (len(groups) - 1) * bar_width / 2  # Centered offsets

    # Create the plot
    fig, ax1 = plt.subplots(figsize=(8, 6))
    y_values_primary_total = []  # y_values_primary_total holds all metrics across all tools - I could make this a nested dict if desired
    for i, group in enumerate(groups):
        y_values_primary = [metric_dict_collection[group][metric] for metric in primary_metrics]  # y_values_primary just holds all metrics for a specific tool
        bars_primary = ax1.bar(x_primary + offsets[i], y_values_primary, bar_width, label=group, color=colors[i])

        # Add value annotations for the primary metrics
        if display_numbers:
            for specific_bar, value in zip(bars_primary, y_values_primary):
                ax1.text(specific_bar.get_x() + specific_bar.get_width() / 2, specific_bar.get_height(), f"{value:.3f}", ha="center", va="bottom", fontsize=10)

        y_values_primary_total.extend(y_values_primary)

    # Customize the plot
    # ax1.set_xlabel("Metrics", fontsize=12)
    # ax1.set_title("Comparison of Metrics Across Tools", fontsize=14)
    if "accuracy" in primary_metrics or "sensitivity" in primary_metrics or "specificity" in primary_metrics:
        ax1.set_ylim(0, 1.05)
    ax1.set_xticks(x_primary)
    ax1.set_xticklabels(primary_metrics, fontsize=12)

    metric_to_tool_to_p_value_dict_of_dicts = {}
    margins_of_error = {}
    if show_p_values or output_file_p_values:
        for metric in primary_metrics:
            margins_of_error[metric] = {}
            if metric in {"accuracy", "sensitivity", "specificity"}:
                tool_to_p_value_dict_aggregate = calculate_mcnemar(unique_mcrs_df, tools=groups, metric=metric)  # don't pass output file here because I want output file to include all p-values; and don't do bonferroni here for a similar reason
            elif metric in {"mean_magnitude_expression_error"}:
                tool_to_p_value_dict_aggregate = calculate_paired_t_test(unique_mcrs_df, column_root="mutation_expression_prediction_error", tools=groups, take_absolute_value=True)
                for group in groups:  # calculate 95% confidence intervals
                    mean_value, margin_of_error = compute_95_confidence_interval_margin_of_error(unique_mcrs_df[f"mutation_expression_prediction_error_{group}"], take_absolute_value=True)
                    margins_of_error[metric][group] = (mean_value, margin_of_error)
            elif metric in {"mean_expression_error"}:
                tool_to_p_value_dict_aggregate = calculate_paired_t_test(unique_mcrs_df, column_root="mutation_expression_prediction_error", tools=groups, take_absolute_value=False)
                for group in groups:  # calculate 95% confidence intervals
                    mean_value, margin_of_error = compute_95_confidence_interval_margin_of_error(unique_mcrs_df[f"mutation_expression_prediction_error_{group}"], take_absolute_value=False)
                    margins_of_error[metric][group] = (mean_value, margin_of_error)
            else:
                raise ValueError(f"Invalid metric for p-value calculation: {metric}. Valid options are 'accuracy', 'sensitivity', 'specificity', 'mean_magnitude_expression_error', 'mean_expression_error'")

            metric_to_tool_to_p_value_dict_of_dicts[metric] = tool_to_p_value_dict_aggregate

        if bonferroni:
            n_tests = count_leaves(metric_to_tool_to_p_value_dict_of_dicts)  # counts the number of leaves in the nested dictionary - makes it generalizable  # (len(groups)-1) * len(primary_metrics)
            for metric in primary_metrics:
                for group in groups:
                    if metric in metric_to_tool_to_p_value_dict_of_dicts and group in metric_to_tool_to_p_value_dict_of_dicts[metric]:
                        metric_to_tool_to_p_value_dict_of_dicts[metric][group] = min(metric_to_tool_to_p_value_dict_of_dicts[metric][group] * n_tests, 1.0)

        # Save to a file
        if output_file_p_values:
            with open(output_file_p_values, "w", encoding="utf-8") as f:
                json.dump(metric_to_tool_to_p_value_dict_of_dicts, f, indent=4)

        # # toy p-values for accuracy, sensitivity, and specificity
        # metric_to_tool_to_p_value_dict_of_dicts = {"accuracy": {
        #     "gatk_mutect2": 0.0001,
        #     "gatk_haplotypecaller": 0.04
        # },
        # "sensitivity": {
        #     "gatk_mutect2": 0.007,
        #     "gatk_haplotypecaller": 0.99
        # },
        # "specificity": {
        #     "gatk_mutect2": 0.02,
        #     "gatk_haplotypecaller": 0.99
        # }}

        # # toy p-values for mean_magnitude_expression_error
        # metric_to_tool_to_p_value_dict_of_dicts = {"mean_magnitude_expression_error": {
        #     "gatk_mutect2": 0.0001,
        #     "gatk_haplotypecaller": 0.04
        # }}

        if show_p_values:
            for i, metric in enumerate(primary_metrics):
                if metric in metric_to_tool_to_p_value_dict_of_dicts:
                    number_of_p_values_in_this_cluster = 0
                    for j, group in enumerate(groups):
                        # 95% confidence intervals
                        if metric in margins_of_error and group in margins_of_error[metric]:
                            # Calculate error values
                            mean_value, margin_of_error = margins_of_error[metric][group]
                            if margin_of_error != 0:
                                yerr = [margin_of_error, margin_of_error]
                                x_value = x_primary[i] + offsets[j]

                                # Plot the point with the confidence interval
                                ax1.errorbar(x_value, mean_value, yerr=np.array([yerr]).T, fmt="", capsize=5, label="Mean with 95% CI", color="black")  # Transpose to match dimensions  # Marker for the point  # Adds caps to the error bars

                        # p-values
                        if group in metric_to_tool_to_p_value_dict_of_dicts[metric]:
                            p_value = metric_to_tool_to_p_value_dict_of_dicts[metric][group]
                            if p_value >= 0.05:  # * increase these values to show more p-values for debugging
                                continue
                            elif p_value < 0.05 and p_value >= 0.01:
                                symbol = "*"
                            elif p_value < 0.01 and p_value >= 0.001:
                                symbol = "**"
                            else:
                                symbol = "***"

                            start_x = x_primary[i] + offsets[0]  # assuming varseek is first element
                            end_x = x_primary[i] + offsets[j]

                            if metric in {"accuracy", "sensitivity", "specificity"}:
                                y_start = max(y_values_primary_total) + (number_of_p_values_in_this_cluster * 0.05) + 0.1
                            else:
                                y_start = (max(y_values_primary_total) + (number_of_p_values_in_this_cluster * 1.7)) * 1.08  # 1.7 (left constant) adjusts based on other bars; 1.08 (right constant) adjusts to make sure it doesn't hit the top bar

                            y_end = y_start

                            ax1.plot([start_x, start_x, end_x, end_x], [y_start, y_end, y_end, y_start], lw=1.5, c="k")  # plot the bar
                            ax1.text((start_x + end_x) * 0.5, y_end, symbol, ha="center", va="bottom", color="k")  # plot the asterisk(s)

                            number_of_p_values_in_this_cluster += 1

    # ax1.legend(title="Tools", loc="upper left", bbox_to_anchor=(1.05, 1))
    ax1.grid(axis="y", linestyle="--", alpha=0.7)

    # Show the plot
    plt.tight_layout()
    if output_file:
        if os.path.dirname(output_file):
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()


def calculate_grouped_metric(grouped_df, y_metric, tool):
    TP_column = f"TP_{tool}"
    FP_column = f"FP_{tool}"
    FN_column = f"FN_{tool}"
    TN_column = f"TN_{tool}"
    mutation_expression_prediction_error_column = f"mutation_expression_prediction_error_{tool}"
    y_metric_output_column = f"{y_metric}_{tool}"

    if y_metric == "accuracy":
        grouped_df[y_metric_output_column] = (grouped_df[TP_column] + grouped_df[TN_column]) / (grouped_df[TP_column] + grouped_df[TN_column] + grouped_df[FP_column] + grouped_df[FN_column])
    elif y_metric == "sensitivity":
        grouped_df[y_metric_output_column] = grouped_df[TP_column] / (grouped_df[TP_column] + grouped_df[FN_column])
        grouped_df.loc[(grouped_df[TP_column] + grouped_df[FN_column]) == 0, y_metric] = 1.0
    elif y_metric == "specificity":
        grouped_df[y_metric_output_column] = grouped_df[TN_column] / (grouped_df[TN_column] + grouped_df[FP_column])
        grouped_df.loc[(grouped_df[TN_column] + grouped_df[FP_column]) == 0, y_metric] = 1.0
    elif y_metric in {"mean_magnitude_expression_error", "mean_expression_error"}:
        grouped_df[y_metric_output_column] = grouped_df[mutation_expression_prediction_error_column] / grouped_df["number_of_elements_in_the_group"]
    else:
        raise ValueError(f"Invalid y_metric: {y_metric}. Valid options are 'accuracy', 'sensitivity', 'specificity', and 'mutation_expression_prediction_error'")

    return grouped_df


def create_stratified_metric_line_plot(unique_mcrs_df, x_stratification, y_metric, tools, bins=None, keep_strict_bins=False, show_p_values=False, show_confidence_intervals=False, bonferroni=True, output_file=None, show=True, output_file_p_values=None, filter_real_negatives=False):
    if x_stratification not in unique_mcrs_df.columns:
        raise ValueError(f"Invalid x_stratification: {x_stratification}. Valid options are {unique_mcrs_df.columns.tolist()}")

    # removes unnecessary columns for the function
    columns_to_keep_for_function = list(set([x_stratification, "included_in_synthetic_reads_mutant", "number_of_reads_mutant"]))
    for tool in tools:
        columns_to_keep_for_function.extend([f"TP_{tool}", f"TN_{tool}", f"FP_{tool}", f"FN_{tool}", f"mutation_expression_prediction_error_{tool}", f"mutation_detected_{tool}", f"DP_{tool}"])
    for column in columns_to_keep_for_function:
        if column not in unique_mcrs_df.columns:
            columns_to_keep_for_function.remove(column)

    unique_mcrs_df = unique_mcrs_df.loc[:, columns_to_keep_for_function].copy()  # make a copy to avoid modifying the original DataFrame
    # unique_mcrs_df = unique_mcrs_df.copy()  # make a copy to avoid modifying the original DataFrame

    if "expression_error" in y_metric and not filter_real_negatives and x_stratification not in {"number_of_reads_mutant"}:
        print("Warning: filtering real negatives is recommended when using expression error as a primary metric and stratifying by something other than number_of_reads_mutant, but this setting is not currently enabled. Recommended: filter_real_negatives = True.")

    if y_metric == "sensitivity" or filter_real_negatives:
        unique_mcrs_df = unique_mcrs_df[(unique_mcrs_df["included_in_synthetic_reads_mutant"] == True) & (unique_mcrs_df["number_of_reads_mutant"] > 0)]
    elif y_metric == "specificity":
        unique_mcrs_df = unique_mcrs_df[(unique_mcrs_df["included_in_synthetic_reads_mutant"] == False) & (unique_mcrs_df["number_of_reads_mutant"] == 0)]

    if keep_strict_bins:
        unique_mcrs_df = unique_mcrs_df[unique_mcrs_df[x_stratification].astype(int).isin(bins)]

    # Prepare for plotting
    plt.figure(figsize=(10, 6))

    x_values_raw = sorted(list(unique_mcrs_df[x_stratification].unique()))

    x_stratification_original = x_stratification  # to label x-axis, as x_stratification may be changed to "bin"

    if bins and not keep_strict_bins:  # remember bins are left-inclusive and right-exclusive
        if bins[-2] > x_values_raw[-1]:
            raise ValueError(f"Invalid bins: {bins}. The 2nd to last bin value {bins[-2]} is greater than the maximum value in the data {x_values_raw[-1]}")
        # list comprehension to assign labels list
        labels = [f"({bins[i]}, {bins[i+1]}]" for i in range(len(bins) - 1)]  # eg bins [0, 0.25, 0.5, 0.75, 1] --> labels ["[0, 0.25)", "[0.25, 0.5)", "[0.5, 0.75)", "[0.75, 1)"]

        # replace "inf" with true start and end values
        if "-inf" in labels[0]:
            labels[0] = labels[0].replace("-inf", str(x_values_raw[0]))
        if "inf" in labels[-1]:
            labels[-1] = labels[-1].replace("inf", str(x_values_raw[-1]))

        number_of_rows_before_filtering = len(unique_mcrs_df)
        # remove rows lower than lower bound or higher than upper bound
        unique_mcrs_df = unique_mcrs_df[(unique_mcrs_df[x_stratification] >= bins[0]) & (unique_mcrs_df[x_stratification] < bins[-1])]
        number_of_rows_after_filtering = len(unique_mcrs_df)

        if number_of_rows_before_filtering != number_of_rows_after_filtering:
            print(f"Removed {number_of_rows_before_filtering - number_of_rows_after_filtering} rows due to binning.")

        # Assign bins to a new column
        unique_mcrs_df["bin"] = pd.cut(unique_mcrs_df[x_stratification], bins=bins, labels=labels, right=True, include_lowest=False)

        x_values = labels
        x_stratification = "bin"
    else:
        x_values = x_values_raw

    if not keep_strict_bins:
        x_indices = range(len(x_values))
    else:
        x_indices = x_values

    if y_metric == "mutation_expression_prediction_error":
        for tool in tools:
            if f"mutation_expression_prediction_error_{tool}" not in unique_mcrs_df.columns:
                raise ValueError(f"mutation_expression_prediction_error_{tool} not in unique_mcrs_df.columns")

    # created grouped_df
    if y_metric == "mean_magnitude_expression_error":  # calculate sum of magnitudes for this one column
        aggregation_functions = {}
        for tool in tools:
            # Group by tumor purity and calculate sensitivity
            aggregation_functions[f"mutation_expression_prediction_error_{tool}"] = lambda x: x.abs().sum()  # Sum of absolute values

        # Use the default sum for all other columns
        grouped_df = unique_mcrs_df.groupby(x_stratification).agg({col: aggregation_functions.get(col, "sum") for col in unique_mcrs_df.columns if col != x_stratification})  # sum is the default aggregation function
    else:  # including if y_metric == mean_expression_error:  # calculate sum for all columns
        grouped_df = unique_mcrs_df.groupby(x_stratification).sum(numeric_only=True)

    grouped_df["number_of_elements_in_the_group"] = unique_mcrs_df.groupby(x_stratification).size()

    # redundant code for calculating y-max (because I need this for setting p-value asterisk height)
    if y_metric in {"accuracy", "sensitivity", "specificity"}:
        custom_y_limit = 1.05
    else:
        custom_y_limit = 0
        for i, tool in enumerate(tools):
            grouped_df = calculate_grouped_metric(grouped_df, y_metric, tool)
            y_metric_tool_specific = f"{y_metric}_{tool}"
            custom_y_limit = max(custom_y_limit, grouped_df[y_metric_tool_specific].max())

    number_of_valid_p_values = 0
    nested_dict = lambda: defaultdict(nested_dict)
    metric_to_tool_to_p_value_dict_of_dicts = nested_dict()

    for i, tool in enumerate(tools):
        grouped_df = calculate_grouped_metric(grouped_df, y_metric, tool)
        y_metric_tool_specific = f"{y_metric}_{tool}"  # matches column created by calculate_grouped_metric - try not changing this name if possible

        # Plot sensitivity as a function of tumor purity
        plt.plot(x_indices, grouped_df[y_metric_tool_specific], label=tool, marker="o", color=color_map_20[i])  # use grouped_df.index to get the x-axis values and plot numerically (vs converting to categorical)

        if (show_p_values or output_file_p_values) and tool != "varseek":  # because varseek is the reference tool
            p_value_list = []
            confidence_intervals_list = []
            for x_value in x_values:
                filtered_unique_mcrs_df_for_p_value = unique_mcrs_df.loc[unique_mcrs_df[x_stratification] == x_value]
                if len(filtered_unique_mcrs_df_for_p_value) > 1:
                    number_of_valid_p_values += 1
                    if y_metric in {"accuracy", "sensitivity", "specificity"}:  # Mcnemar
                        p_value = calculate_individual_mcnemar(filtered_unique_mcrs_df_for_p_value, "mutation_detected_varseek", f"mutation_detected_{tool}")
                        margin_of_error = 0

                    elif y_metric in {"mean_magnitude_expression_error"}:  # paired t-test
                        p_value = calculate_individual_paired_t_test(filtered_unique_mcrs_df_for_p_value, column1="mutation_expression_prediction_error_varseek", column2=f"mutation_expression_prediction_error_{tool}", take_absolute_value=True)
                        _, margin_of_error = compute_95_confidence_interval_margin_of_error(filtered_unique_mcrs_df_for_p_value[f"mutation_expression_prediction_error_{tool}"], take_absolute_value=True)

                    elif y_metric in {"mean_expression_error"}:
                        p_value = calculate_individual_paired_t_test(filtered_unique_mcrs_df_for_p_value, column1="mutation_expression_prediction_error_varseek", column2=f"mutation_expression_prediction_error_{tool}", take_absolute_value=False)
                        _, margin_of_error = compute_95_confidence_interval_margin_of_error(filtered_unique_mcrs_df_for_p_value[f"mutation_expression_prediction_error_{tool}"], take_absolute_value=False)

                    else:
                        raise ValueError(f"Invalid metric for p-value calculation: {y_metric}. Valid options are 'accuracy', 'sensitivity', 'specificity', 'mean_magnitude_expression_error', 'mean_expression_error'")

                    p_value_list.append(p_value)
                    confidence_intervals_list.append(margin_of_error)
                else:
                    p_value_list.append(1.0)

            # bonferroni
            if bonferroni:
                n_tests = number_of_valid_p_values * (len(tools) - 1)  # because varseek is the reference tool
                p_value_list = [min((p_value * n_tests), 1.0) for p_value in p_value_list]

            # Plot '*' above points where p-value < 0.05
            if show_p_values:
                for x, y, p_value, margin_of_error in zip(x_indices, list(grouped_df[y_metric_tool_specific]), p_value_list, confidence_intervals_list):
                    if show_confidence_intervals and margin_of_error != 0:
                        # confidence interval errors
                        yerr = [margin_of_error, margin_of_error]

                        # Plot the point with the confidence interval
                        plt.errorbar(x, y, yerr=np.array([yerr]).T, fmt="", capsize=5, label="Mean with 95% CI", color=color_map_20[i])  # Transpose to match dimensions  # Marker for the point  # Adds caps to the error bars

                    # p-values
                    metric_to_tool_to_p_value_dict_of_dicts[str(x)][tool] = p_value
                    if p_value >= 0.05:  # * increase these values to show more p-values for debugging
                        continue
                    elif p_value < 0.05 and p_value >= 0.01:
                        symbol = "*"
                    elif p_value < 0.01 and p_value >= 0.001:
                        symbol = "**"
                    else:
                        symbol = "***"
                    plt.text(x, y + (custom_y_limit * 0.01), symbol, color=color_map_20[i], fontsize=12, ha="center")  # Slightly above the point

    # # Set x-axis to log2 scale
    # if log:  # log can be False (default, not log) or True (defaults to 2) or int (log base)
    #     if log is True:
    #         log = 2
    #     plt.xscale("log", base=log)

    if output_file_p_values:
        with open(output_file_p_values, "w", encoding="utf-8") as f:
            json.dump(metric_to_tool_to_p_value_dict_of_dicts, f, indent=4)

    # Customize plot
    if y_metric in {"accuracy", "sensitivity", "specificity"}:  # "accuracy" in primary_metrics or "sensitivity" in primary_metrics or "specificity" in primary_metrics:
        plt.ylim(0, custom_y_limit)

    x_values_new = []
    for x_value in x_values:  # add the number of elements in each stratification below the x-axis label
        number_of_elements = grouped_df.loc[grouped_df.index == x_value]["number_of_elements_in_the_group"].iloc[0]
        x_values_new.append(f"{x_value}\n(n={number_of_elements})")
    x_values = x_values_new

    plt.xticks(ticks=x_indices, labels=x_values)
    plt.xlabel(x_stratification_original, fontsize=12)
    plt.ylabel(y_metric, fontsize=12)
    # plt.legend(title="Tools")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Show the plot
    plt.tight_layout()
    if output_file:
        if os.path.dirname(output_file):
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()


def create_benchmarking_legend(tools, outfile=None, show=True):
    # Define colors
    colors = color_map_20[: len(tools)]

    # Create a figure for the legend
    fig, ax = plt.subplots(figsize=(0.1, 0.1))  # Adjust size as needed to change whitespace margins

    # Create proxy artists (invisible items for the legend)
    proxies = [plt.Line2D([0], [0], color=color, lw=4) for color in colors]

    # Add the legend to the figure
    ax.legend(proxies, tools, title="Legend", loc="center")
    ax.axis("off")  # Turn off the axes

    # Show the legend-only figure
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()


def write_p_values_to_file(tool_to_p_value_dict, out_file):
    if out_file.endswith(".txt"):
        with open(out_file, "w", encoding="utf-8") as f:
            for tool, p_value in tool_to_p_value_dict.items():
                f.write(f"{tool}: {p_value}\n")
    elif out_file.endswith(".json"):
        # Save to a file
        with open(out_file, "w", encoding="utf-8") as f:
            json.dump(tool_to_p_value_dict, f, indent=4)
    else:
        raise ValueError(f"Invalid file extension: {out_file}. Accepted extensions are '.txt' and '.json'.")


def calculate_individual_mcnemar(df, column1, column2):
    # Counting the true/false values
    contingency_table = pd.crosstab(df[column1], df[column2])

    # McNemar's test
    exact = False if (contingency_table.values >= 25).all() else True
    result = mcnemar(contingency_table, exact=exact)
    return result.pvalue


# see the function above for simple mncnemar, as the following function was written specifically for figure 2c
def calculate_mcnemar(unique_mcrs_df, tools, metric="accuracy", out_file=None, bonferroni=False, do_sensitivity_specificity_filtering=True):
    # unique_mcrs_df = unique_mcrs_df.copy()  # make a copy to avoid modifying the original DataFrame

    # Filtering the DataFrame
    if do_sensitivity_specificity_filtering:
        if metric == "sensitivity":
            filtered_df = unique_mcrs_df[unique_mcrs_df["included_in_synthetic_reads_mutant"] == True]
        elif metric == "specificity":
            filtered_df = unique_mcrs_df[unique_mcrs_df["included_in_synthetic_reads_mutant"] == False]
        elif metric == "accuracy":
            filtered_df = unique_mcrs_df
        else:
            raise ValueError(f"Invalid metric: {metric}. Accepted values are 'sensitivity', 'specificity', and 'accuracy'.")
    else:
        filtered_df = unique_mcrs_df

    tool_to_p_value_dict = {}
    for tool in tools:
        if tool == "varseek":
            continue  # Skip varseek as it is the reference tool

        tool_to_p_value_dict[tool] = calculate_individual_mcnemar(filtered_df, "mutation_detected_varseek", f"mutation_detected_{tool}")

    # Bonferroni correction
    if bonferroni:
        n_tests = count_leaves(tool_to_p_value_dict)
        for tool in tools:
            tool_to_p_value_dict[tool] = min(tool_to_p_value_dict[tool] * n_tests, 1.0)

    if out_file:
        write_p_values_to_file(tool_to_p_value_dict, out_file)

    return tool_to_p_value_dict


def calculate_individual_paired_t_test(df, column1, column2, tails=2, larger_column_expected=None, take_absolute_value=False):
    df = df.copy()  # make a copy to avoid modifying the original DataFrame

    if take_absolute_value:
        df[column1] = df[column1].abs()
        df[column2] = df[column2].abs()

    t_stat, p_value = ttest_rel(df[column1], df[column2])  # if one-tailed, then the sign of t_stat indicates if 1st arg > 2nd arg

    if tails == 1:
        if not larger_column_expected:
            raise ValueError("larger_column_expected must be provided when tails == 1")
        if larger_column_expected != column1 and larger_column_expected != column2:
            raise ValueError("larger_column_expected must be one of the two columns passed in")

        if t_stat < 0 and larger_column_expected == column1:  # corresponds to varseek being passed in first above
            p_value /= 2
        elif t_stat > 0 and larger_column_expected != column2:
            p_value /= 2
        else:
            p_value = 1.0

    return p_value


# make sure I've already computed the difference between predicted expression and true expression for each tool (i.e., DP_TOOL - number_of_reads_mutant)
# two-tailed - to make one-tailed, divide by 2 and make sure that
def calculate_paired_t_test(unique_mcrs_df, column_root, tools, take_absolute_value=False, out_file=None, bonferroni=False, tails=2, larger_column_expected="varseek"):
    unique_mcrs_df = unique_mcrs_df.copy()  # make a copy to avoid modifying the original DataFrame

    tool_to_p_value_dict = {}
    for tool in tools:
        if tool == "varseek":
            continue  # Skip varseek as it is the reference tool

        if larger_column_expected == "varseek":
            larger_column_expected = f"{column_root}_varseek"
        elif larger_column_expected in tools:
            larger_column_expected = f"{column_root}_{tool}"

        p_value = calculate_individual_paired_t_test(unique_mcrs_df, column1=f"{column_root}_varseek", column2=f"{column_root}_{tool}", take_absolute_value=take_absolute_value, tails=tails, larger_column_expected=larger_column_expected)

        tool_to_p_value_dict[tool] = p_value

    # Bonferroni correction
    if bonferroni:
        n_tests = count_leaves(tool_to_p_value_dict)
        for tool in tools:
            tool_to_p_value_dict[tool] = min(tool_to_p_value_dict[tool] * n_tests, 1.0)

    if out_file:
        write_p_values_to_file(tool_to_p_value_dict, out_file)

    return tool_to_p_value_dict


def print_json(json_file):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    print(json.dumps(data, indent=4))


def count_leaves(d):
    """
    Recursively counts the number of leaves in a nested dictionary.

    Args:
        d (dict): The dictionary to traverse.

    Returns:
        int: The count of leaves.
    """
    if not isinstance(d, dict):  # Base case: not a dictionary
        return 1
    return sum(count_leaves(value) for value in d.values())


def compute_95_confidence_interval_margin_of_error(values, take_absolute_value=False):
    # values = unique_mcrs_df[column]
    if take_absolute_value:
        values = values.abs()

    # Step 1: Compute mean
    mean = np.mean(values)

    # Step 2: Compute Standard Error of the Mean (SEM)
    sem = np.std(values, ddof=1) / np.sqrt(len(values))

    # Step 3: Degrees of freedom
    df = len(values) - 1

    # Step 4: Critical t-value for 95% confidence
    t_critical = t.ppf(0.975, df)  # Two-tailed 95% confidence

    # Step 5: Margin of error
    margin_of_error = t_critical * sem

    # # Step 6: Compute confidence interval
    # ci_lower = mean - margin_of_error
    # ci_upper = mean + margin_of_error

    return mean, margin_of_error  # , ci_lower, ci_upper


def create_stratified_metric_bar_plot_updated(unique_mcrs_df, x_stratification, y_metric, tools, display_numbers=False, show_p_values=False, show_confidence_intervals=True, bonferroni=True, output_file=None, show=True, output_file_p_values=None, filter_real_negatives=False):
    if x_stratification not in unique_mcrs_df.columns:
        raise ValueError(f"Invalid x_stratification: {x_stratification}. Valid options are {unique_mcrs_df.columns.tolist()}")

    # removes unnecessary columns for the function
    columns_to_keep_for_function = list(set([x_stratification, "included_in_synthetic_reads_mutant", "number_of_reads_mutant"]))
    for tool in tools:
        columns_to_keep_for_function.extend([f"TP_{tool}", f"TN_{tool}", f"FP_{tool}", f"FN_{tool}", f"mutation_expression_prediction_error_{tool}", f"mutation_detected_{tool}", f"DP_{tool}"])

    columns_to_keep_for_function = [column for column in columns_to_keep_for_function if column in unique_mcrs_df.columns]

    unique_mcrs_df = unique_mcrs_df.loc[:, columns_to_keep_for_function].copy()  # make a copy to avoid modifying the original DataFrame
    # unique_mcrs_df = unique_mcrs_df.copy()  # make a copy to avoid modifying the original DataFrame

    if x_stratification == "vcrs_variant_type":
        # remove any values in "vcrs_variant_type" equal to "mixed"
        unique_mcrs_df = unique_mcrs_df[unique_mcrs_df["vcrs_variant_type"] != "mixed"]

    if "expression_error" in y_metric and not filter_real_negatives and x_stratification not in {"number_of_reads_mutant"}:
        print("Warning: filtering real negatives is recommended when using expression error as a primary metric and stratifying by something other than number_of_reads_mutant, but this setting is not currently enabled. Recommended: filter_real_negatives = True.")

    if y_metric == "sensitivity" or filter_real_negatives:
        unique_mcrs_df = unique_mcrs_df[unique_mcrs_df["included_in_synthetic_reads_mutant"] == True]
    elif y_metric == "specificity":
        unique_mcrs_df = unique_mcrs_df[unique_mcrs_df["included_in_synthetic_reads_mutant"] == False]

    # Prepare for plotting
    plt.figure(figsize=(10, 6))

    if y_metric == "mutation_expression_prediction_error":
        for tool in tools:
            if f"mutation_expression_prediction_error_{tool}" not in unique_mcrs_df.columns:
                raise ValueError(f"mutation_expression_prediction_error_{tool} not in unique_mcrs_df.columns")

    # created grouped_df
    if y_metric == "mean_magnitude_expression_error":  # calculate sum of magnitudes for this one column
        aggregation_functions = {}
        for tool in tools:
            # Group by tumor purity and calculate sensitivity
            aggregation_functions[f"mutation_expression_prediction_error_{tool}"] = lambda x: x.abs().sum()  # Sum of absolute values

        # Use the default sum for all other columns
        grouped_df = unique_mcrs_df.groupby(x_stratification).agg({col: aggregation_functions.get(col, "sum") for col in unique_mcrs_df.columns if col != x_stratification})  # sum is the default aggregation function
    else:  # including if y_metric == mean_expression_error:  # calculate sum for all columns
        grouped_df = unique_mcrs_df.groupby(x_stratification).sum(numeric_only=True)

    grouped_df["number_of_elements_in_the_group"] = unique_mcrs_df.groupby(x_stratification).size()

    # redundant code for calculating y-max (because I need this for setting p-value asterisk height)
    if y_metric in {"accuracy", "sensitivity", "specificity"}:
        custom_y_limit = 1.05
    else:
        custom_y_limit = 0
        for i, tool in enumerate(tools):
            grouped_df = calculate_grouped_metric(grouped_df, y_metric, tool)
            y_metric_tool_specific = f"{y_metric}_{tool}"
            custom_y_limit = max(custom_y_limit, grouped_df[y_metric_tool_specific].max())

    number_of_valid_p_values = 0
    nested_dict = lambda: defaultdict(nested_dict)
    stratification_to_tool_to_p_value_dict_of_dicts = nested_dict()
    stratification_to_tool_to_error_dict_of_dicts = nested_dict()

    # Prepare data
    bar_names = unique_mcrs_df[x_stratification].unique()
    x_primary = np.arange(len(bar_names))  # Positions for the metrics on the x-axis
    bar_width = 0.25  # Width of each primary_metrics
    offsets = np.arange(len(tools)) * bar_width - (len(tools) - 1) * bar_width / 2  # Centered offsets

    # Create the plot
    y_values_primary_total = []  # y_values_primary_total holds all metrics across all tools - I could make this a nested dict if desired
    for i, tool in enumerate(tools):
        grouped_df = calculate_grouped_metric(grouped_df, y_metric, tool)
        y_metric_tool_specific = f"{y_metric}_{tool}"  # matches column created by calculate_grouped_metric - try not changing this name if possible
        y_values_primary = [grouped_df.loc[grouped_df.index == bar_name, y_metric_tool_specific][0] for bar_name in bar_names]  # y_values_primary just holds all metrics for a specific tool
        bars_primary = plt.bar(x_primary + offsets[i], y_values_primary, bar_width, label=tool, color=color_map_20[i])

        # Add value annotations for the primary metrics
        if display_numbers:
            for specific_bar, value in zip(bars_primary, y_values_primary):
                plt.text(specific_bar.get_x() + specific_bar.get_width() / 2, specific_bar.get_height(), f"{value:.3f}", ha="center", va="bottom", fontsize=10)

        y_values_primary_total.extend(y_values_primary)

        if (show_p_values or output_file_p_values) and tool != "varseek":  # because varseek is the reference tool
            p_value_list = []
            confidence_intervals_list = []
            for x_value in bar_names:
                filtered_unique_mcrs_df_for_p_value = unique_mcrs_df.loc[unique_mcrs_df[x_stratification] == x_value]
                if len(filtered_unique_mcrs_df_for_p_value) > 1:
                    number_of_valid_p_values += 1
                    if y_metric in {"accuracy", "sensitivity", "specificity"}:  # Mcnemar
                        p_value = calculate_individual_mcnemar(filtered_unique_mcrs_df_for_p_value, "mutation_detected_varseek", f"mutation_detected_{tool}")
                        mean_value, margin_of_error = 0, 0

                    elif y_metric in {"mean_magnitude_expression_error"}:  # paired t-test
                        p_value = calculate_individual_paired_t_test(filtered_unique_mcrs_df_for_p_value, column1="mutation_expression_prediction_error_varseek", column2=f"mutation_expression_prediction_error_{tool}", take_absolute_value=True)
                        mean_value, margin_of_error = compute_95_confidence_interval_margin_of_error(filtered_unique_mcrs_df_for_p_value[f"mutation_expression_prediction_error_{tool}"], take_absolute_value=True)

                    elif y_metric in {"mean_expression_error"}:
                        p_value = calculate_individual_paired_t_test(filtered_unique_mcrs_df_for_p_value, column1="mutation_expression_prediction_error_varseek", column2=f"mutation_expression_prediction_error_{tool}", take_absolute_value=False)
                        mean_value, margin_of_error = compute_95_confidence_interval_margin_of_error(filtered_unique_mcrs_df_for_p_value[f"mutation_expression_prediction_error_{tool}"], take_absolute_value=False)

                    else:
                        raise ValueError(f"Invalid metric for p-value calculation: {y_metric}. Valid options are 'accuracy', 'sensitivity', 'specificity', 'mean_magnitude_expression_error', 'mean_expression_error'")
                else:
                    p_value = 1.0
                    mean_value, margin_of_error = 0, 0

                p_value_list.append(p_value)
                confidence_intervals_list.append(margin_of_error)
                stratification_to_tool_to_p_value_dict_of_dicts[x_value][tool] = p_value
                stratification_to_tool_to_error_dict_of_dicts[x_value][tool] = (mean_value, margin_of_error)

            # bonferroni
            if bonferroni:
                n_tests = number_of_valid_p_values * (len(tools) - 1)  # because varseek is the reference tool
                p_value_list = [min((p_value * n_tests), 1.0) for p_value in p_value_list]
                for x_value in bar_names:
                    stratification_to_tool_to_p_value_dict_of_dicts[x_value][tool] = min((stratification_to_tool_to_p_value_dict_of_dicts[x_value][tool] * n_tests), 1.0)

    # Plot '*' above points where p-value < 0.05
    if show_p_values:
        for i, bar_name in enumerate(bar_names):
            if bar_name in stratification_to_tool_to_p_value_dict_of_dicts:
                number_of_p_values_in_this_cluster = 0
                for j, tool in enumerate(tools):
                    # 95% confidence intervals
                    if bar_name in stratification_to_tool_to_error_dict_of_dicts and tool in stratification_to_tool_to_error_dict_of_dicts[bar_name]:
                        # Calculate error values
                        mean_value, margin_of_error = stratification_to_tool_to_error_dict_of_dicts[bar_name][tool]
                        if margin_of_error != 0:
                            yerr = [margin_of_error, margin_of_error]
                            x_value = x_primary[i] + offsets[j]

                            # Plot the point with the confidence interval
                            plt.errorbar(x_value, mean_value, yerr=np.array([yerr]).T, fmt="", capsize=5, label="Mean with 95% CI", color="black")  # Transpose to match dimensions  # Marker for the point  # Adds caps to the error bars

                    # p-values
                    if tool in stratification_to_tool_to_p_value_dict_of_dicts[bar_name]:
                        p_value = stratification_to_tool_to_p_value_dict_of_dicts[bar_name][tool]

                        if p_value >= 0.05:  # * increase these values to show more p-values for debugging
                            continue
                        elif 0.01 <= p_value < 0.05:
                            symbol = "*"
                        elif 0.001 <= p_value < 0.01:
                            symbol = "**"
                        else:  # p_value < 0.001:
                            symbol = "***"

                        start_x = x_primary[i] + offsets[0]  # assuming varseek is first element
                        end_x = x_primary[i] + offsets[j]

                        if y_metric in {"accuracy", "sensitivity", "specificity"}:
                            y_start = max(y_values_primary_total) + (number_of_p_values_in_this_cluster * 0.05) + 0.05
                        else:
                            y_start = (max(y_values_primary_total) + (number_of_p_values_in_this_cluster * 1.7)) * 1.08  # 1.7 (left constant) adjusts based on other bars; 1.08 (right constant) adjusts to make sure it doesn't hit the top bar

                        y_end = y_start

                        plt.plot([start_x, start_x, end_x, end_x], [y_start, y_end, y_end, y_start], lw=1.5, c="k")  # plot the bar
                        plt.text((start_x + end_x) * 0.5, y_end, symbol, ha="center", va="bottom", color="k")  # plot the asterisk(s)

                        number_of_p_values_in_this_cluster += 1

    if output_file_p_values:
        with open(output_file_p_values, "w", encoding="utf-8") as f:
            json.dump(stratification_to_tool_to_p_value_dict_of_dicts, f, indent=4)

    if y_metric in {"accuracy", "sensitivity", "specificity"}:  # "accuracy" in primary_metrics or "sensitivity" in primary_metrics or "specificity" in primary_metrics:
        plt.ylim(0, custom_y_limit)

    x_values_new = []
    for x_value in bar_names:  # add the number of elements in each stratification below the x-axis label
        number_of_elements = grouped_df.loc[grouped_df.index == x_value]["number_of_elements_in_the_group"].iloc[0]
        x_values_new.append(f"{x_value}\n(n={number_of_elements})")
    bar_names = x_values_new

    plt.xticks(ticks=x_primary, labels=bar_names, fontsize=12)
    plt.ylabel(y_metric, fontsize=12)

    # Show the plot
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    if output_file:
        if os.path.dirname(output_file):
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()


def plot_frequency_histogram(unique_mcrs_df, column_base, tools, fraction=False, output_file=None, show=True):
    """
    Plots a histogram of the mutation expression prediction errors for each tool.
    """
    errors_dict = {}

    plt.figure(figsize=(10, 6))
    for index, tool in enumerate(tools):
        errors_dict[tool] = unique_mcrs_df.loc[unique_mcrs_df[f"FP_{tool}"], f"{column_base}_{tool}"]
        if fraction:
            total_count = len(errors_dict[tool])  # Total number of errors for this tool
            weights = [1 / total_count] * total_count  # Fractional weights for each error
            y_axis_label = "Fraction of FPs"
        else:
            weights = None  # No weights for absolute counts
            y_axis_label = "Number of FPs"
        plt.hist(errors_dict[tool], bins=30, alpha=0.6, label=tool, color=color_map_20[index], weights=weights)

    # Add labels, legend, and title
    plt.xscale("log", base=2)

    # Customize ticks to show all powers of 2
    log_locator = LogLocator(base=2.0, subs=[], numticks=30)  # `subs=[]` means only major ticks are shown
    log_formatter = FuncFormatter(lambda x, _: f"{int(x)}" if x >= 1 else "")

    ax = plt.gca()
    ax.xaxis.set_major_locator(log_locator)
    ax.xaxis.set_major_formatter(log_formatter)

    plt.xlabel("Counts Detected")
    plt.ylabel(y_axis_label)
    plt.title("Histogram of Counts Detected for FPs")
    plt.legend(loc="upper right")

    # Show the plot
    plt.tight_layout()

    if output_file:
        if os.path.dirname(output_file):
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()


def plot_time_and_memory_benchmarking(df, metric_name, units=None, log_y=False, x_col="Reads", y_col="Value_total", y_thresholds=None, output_file=None, show=True):
    df = df.copy()  # make a copy to avoid modifying the original DataFrame
    df = df.loc[df["Metric"] == metric_name]
       
    tools = df["Tool"].unique()
    tool_colors = {tool: color_map_20[i % len(color_map_20)] for i, tool in enumerate(tools)}

    df[x_col] /= 1_000_000  # Convert to millions
    x_ticks = sorted(df[x_col].unique())
    
    plt.figure(figsize=(6, 4))  # Adjust figure size
    for tool in tools:
        subset = df.loc[df["Tool"] == tool]
        plt.plot(subset[x_col], subset[y_col], marker="o", linestyle="-", color=tool_colors[tool], label=tool)

    ylabel = metric_name
    if units:
        ylabel += f" ({units})"
    
    plt.xlabel("Number of Reads (millions)")
    plt.xticks(x_ticks)  # Set only the unique "Reads" values as x-ticks  (eg the only ticks marked are 1, 4, 16, 64, etc)
    plt.ylabel(ylabel)
    if log_y:
        plt.yscale("log")
        min_value = 10 ** math.floor(math.log10(min(df[y_col])))
        max_value = 10 ** math.ceil(math.log10(max(df[y_col])))
        y_ticks = np.logspace(math.log10(min_value), math.log10(max_value), num=int(math.log10(max_value) - math.log10(min_value)) + 1, base=10)
        plt.yticks(y_ticks)
        plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))  # integers instead of scientific notation
    if y_thresholds is not None:
        for key, value in y_thresholds.items():
            plt.axhline(y=value, color='gray', linestyle="--", label=key)
            plt.text(
                x=plt.xlim()[1],  # Right side of the plot
                y=value, 
                s=key, 
                color='gray', 
                ha='right', 
                va='bottom'  # Position text above the line
            )
    # plt.legend(loc="upper left", fontsize=8)
    # plt.title(f"{metric_name} vs Number of Reads")
    # plt.ticklabel_format(style="plain", axis="x")  # remove the 1e7 in the bottom right
    if output_file:
        if os.path.dirname(output_file):
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()
