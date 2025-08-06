import os
import subprocess
import requests
import pandas as pd
from bs4 import BeautifulSoup

def find_key_with_trace(data, target_key, path=""):
    """
    Recursively search for a key in a nested JSON object (dict or list) and return its full path and values.

    Parameters:
        data (dict or list): The JSON object to search.
        target_key (str): The key to find.
        path (str): Tracks the current path in the JSON object.

    Returns:
        list: A list of tuples containing the path to the key and the associated value.
    """
    found = []

    if isinstance(data, dict):
        for key, value in data.items():
            current_path = f"{path}.{key}" if path else key
            if key == target_key:
                found.append((current_path, value))
            # Recursively search if the value is a nested structure
            if isinstance(value, (dict, list)):
                found.extend(find_key_with_trace(value, target_key, current_path))
    elif isinstance(data, list):
        for index, item in enumerate(data):
            current_path = f"{path}[{index}]"
            # Recursively search each item
            found.extend(find_key_with_trace(item, target_key, current_path))

    return found


def get_experiment_links(search_url, base_url):
    response = requests.get(search_url, timeout=(10, 90))
    response.raise_for_status()  # Ensure the request was successful
    soup = BeautifulSoup(response.text, "html.parser")
    # Find all experiment links
    links = [base_url + link["href"] for link in soup.find_all("a", href=True) if "/experiments/" in link["href"]]
    return list(set(links))  # Remove duplicates if any


def make_entex_df():
    print("making entex df")
    base_url = "https://www.encodeproject.org"
    search_url = "https://www.encodeproject.org/search/?type=Experiment&status=released&internal_tags=ENTEx&files.file_type=fastq&control_type!=*&limit=200&assay_title=total+RNA-seq"
    experiment_links = get_experiment_links(search_url=search_url, base_url=base_url)

    entex_list_of_dicts = []

    for experiment_link in experiment_links:
        entex_experiment_dict = {}
        experiment_id = experiment_link.split("/")[-2]
        json_data_url = f"{experiment_link}/?format=json"

        response = requests.get(json_data_url, timeout=(10, 90))
        response.raise_for_status()
        experiment_data = response.json()

        # donor_values = find_key_with_trace(experiment_data, 'donor')  # to find the trace in the json - returns as list of key_path, value pairs

        description = experiment_data["description"]
        tissue = experiment_data["biosample_ontology"]["term_name"]
        organ_slims = experiment_data["biosample_ontology"]["organ_slims"]
        age = experiment_data["replicates"][0]["library"]["biosample"]["donor"]["age"]  # same as donor_values[0][1]['age']
        sex = experiment_data["replicates"][0]["library"]["biosample"]["donor"]["sex"]  # same as donor_values[0][1]['sex']

        # Extract FASTQ file links
        fastq_files_metadata = [file for file in experiment_data.get("files", []) if file["file_type"] == "fastq" and file["replicate"]["technical_replicate_number"] == 1]

        if len(fastq_files_metadata) > 2:
            date_created_to_match = fastq_files_metadata[0]["date_created"].split("T")[0]
            fastq_files_metadata = [metadata for metadata in fastq_files_metadata if metadata["date_created"].split("T")[0] == date_created_to_match]

        fastq_files = [file["href"] for file in fastq_files_metadata]
        read_lengths = [file["read_length"] for file in fastq_files_metadata]
        paired_ends = [file["paired_end"] for file in fastq_files_metadata]
        date_created = [file["date_created"] for file in fastq_files_metadata]

        # reorder fastq_files according to paired_ends
        fastq_files = [fastq_files[int(i) - 1] for i in paired_ends]

        if len(set(read_lengths)) != 1:
            print(f"skipping experiment {experiment_id} due to having {len(set(read_lengths))} read lengths: {set(read_lengths)}")

        read_length = read_lengths[0]

        fastq_links = [base_url + link for link in fastq_files]

        if len(fastq_links) != 2:
            print(f"skipping experiment {experiment_id} due to having {len(fastq_links)} fastq links: {fastq_links}")

        entex_experiment_dict["experiment_id"] = experiment_id
        entex_experiment_dict["description"] = description
        entex_experiment_dict["tissue"] = tissue
        entex_experiment_dict["organ_slims"] = organ_slims
        entex_experiment_dict["age"] = age
        entex_experiment_dict["sex"] = sex
        entex_experiment_dict["read_length"] = read_length
        entex_experiment_dict["fastq_link_pair_1"] = fastq_links[0]
        entex_experiment_dict["fastq_link_pair_2"] = fastq_links[1]

        # experiment_data["biosample_ontology"]["term_name"]

        entex_list_of_dicts.append(entex_experiment_dict)

    entex_df = pd.DataFrame(entex_list_of_dicts)

    return entex_df


def download_entex_fastq_links(entex_df, tissue=None, data_download_base="."):
    print("downloading fastq files")
    if tissue:
        entex_df_tissue_selection = entex_df.loc[entex_df["tissue"] == tissue].reset_index(drop=True)
    else:
        entex_df_tissue_selection = entex_df

    # iterate through rows of entex_df_tissue_selection
    number_of_samples = len(entex_df_tissue_selection)
    for index, row in entex_df_tissue_selection.iterrows():
        tissue_underscore_separated = row["tissue"].replace(" ", "_")
        sample_base_dir = os.path.join(data_download_base, tissue_underscore_separated, row["experiment_id"])
        for i in [1, 2]:
            pair_dir = f"{sample_base_dir}/pair{i}"
            os.makedirs(pair_dir, exist_ok=True)

            link = row[f"fastq_link_pair_{i}"]

            download_command = f"wget -c --tries=20 --retry-connrefused -P {pair_dir} {link}"

            if not os.path.exists(f"{pair_dir}/{link.split('/')[-1]}"):
                try:
                    print(f"Downloading sample {index + 1}/{number_of_samples}, pair {i}/2 to {pair_dir}")
                    result = subprocess.run(download_command, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error downloading {link} to {pair_dir}")
                    print(e)
            else:
                print(f"File {pair_dir}/{link.split('/')[-1]} already exists, skipping download")