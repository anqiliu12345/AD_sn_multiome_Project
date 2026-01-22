### Merge cisTopic objects

import pandas as pd
import pickle
import os
import argparse

# Get cell type from command-line arguments
parser = argparse.ArgumentParser(description="Process cell type.")
parser.add_argument('cell_type', type=str, help="The cell type to process.")
args = parser.parse_args()
cell_type = args.cell_type

# Load Donor_IDs
file = "/work/aliu10/AD_sc_Project/sea_AD_data/MTG/multiome/sample_donor_ID_matching.csv"
donor_data = pd.read_csv(file, header=None, sep="\t", names=["Fragment_File", "Donor_ID"]).iloc[:, 1]

# Load cisTopic objects
cistopic_obj_list = []
for sample_id in donor_data.values.tolist():
    pickle_path = os.path.join("/work/aliu10/AD_sc_Project/results", sample_id, f"cistopic_obj_{cell_type}.pkl")
    
    if not os.path.isfile(pickle_path):
        print(f"File not found for {sample_id}, skipping: {pickle_path}")
        continue
        
    with open(pickle_path, "rb") as file:
        cistopic_obj = pickle.load(file)
        cistopic_obj_list.append(cistopic_obj)

# Merge cisTopic objects
cistopic_obj = cistopic_obj_list[0]
cistopic_obj.merge(cistopic_obj_list[1:])

# Save object
pickle.dump(
    cistopic_obj,
    open(os.path.join("/work/aliu10/AD_sc_Project/results/all", f"cistopic_obj_{cell_type}.pkl"), "wb")
)

