import os
import json

cancer_types_to_keep = ["NSCLC", "melanoma", "colorectal_adenocarcinoma", "ovary_adenocarcinoma", "Ductal Adenocarcinoma", "exocrine", "gastric_adenocarcinoma", "upper_aerodigestive_squamous", "hepatocellular_carcinoma", "bladder_carcinoma", "renal_cell_carcinoma"]  # tissues to keep  # # cancer_types_to_keep = ["LUNG", "SKIN", "LARGE_INTESTINE", "OVARY", "PANCREAS", "STOMACH", "UPPER_AERODIGESTIVE_TRACT", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LIVER", "KIDNEY"]
number_to_keep = 5  # Maximum number of records to keep per strategy

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
sequencing_data_out_base = os.path.join(data_dir, "ccle_data_base")
json_path = os.path.join(sequencing_data_out_base, f"ccle_metadata.json")
output_txt_path = json_path.replace(".json", "_filtered_experiment_aliases.txt")

with open(json_path, 'r', encoding="utf-8") as file:
    data = json.load(file)

filtered_data = []

# Filter records by library strategies
for cancer_type in cancer_types_to_keep:
    # Get records matching the current strategy
    strategy_records = [
        study for study in data 
        if study['library_strategy'] == 'RNA-Seq' and 
        'subtype_disease' in study and 
        'lineage_subtype' in study and 
        (study['subtype_disease'] == cancer_type or study['lineage_subtype'] == cancer_type)
    ]
    
    # Limit the number of records to `number_to_keep` or however many are available
    filtered_data.extend(strategy_records[:number_to_keep])

experimant_aliases = []
for study in filtered_data:
    experiment_alias = study.get('experiment_alias')
    experiment_alias_underscores_only = experiment_alias.replace("-", "_")
    experimant_aliases.append(experiment_alias_underscores_only)

# write to text file
with open(output_txt_path, 'w', encoding="utf-8") as file:
    for alias in experimant_aliases:
        file.write(f"{alias}\n")

print(f"Filtered experiment aliases to {output_txt_path}")