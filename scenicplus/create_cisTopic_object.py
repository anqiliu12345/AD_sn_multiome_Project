### Creating a cisTopic object for each sample using consensus peaks

import os
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
import pandas as pd
import scanpy as sc
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
import pickle
import scrublet as scr
import argparse

# Define the path to fragments and data directories
data_path = "/work/aliu10/AD_sc_Project/sea_AD_data/MTG/multiome/processed_data"

# Get donor_id from command-line arguments
parser = argparse.ArgumentParser(description="Process donor ID.")
parser.add_argument('donor_id', type=str, help="The donor ID to process.")
args = parser.parse_args()
sample_id = args.donor_id

# Construct file paths
fragment_file_path = os.path.join(data_path, f"{sample_id}_filtered_fragments.tsv.gz")

# Load the h5ad file and append it to the list
print(f"Load donor {sample_id}...")
fragments_dict = {sample_id: fragment_file_path}

metadata_file = os.path.join(data_path, "metadata", f"metadata_{sample_id}.csv")
cell_data = pd.read_csv(metadata_file, index_col=0) 

def define_cell_type(row):
    if row['Subclass'] in ['Oligodendrocyte', 'Astrocyte', 'Microglia-PVM', 'OPC', 'Endothelial', 'VLMC']:
        return row['Subclass']
    elif row['Class'] == 'Neuronal: Glutamatergic':
        return 'Excitatory'
    elif row['Class'] == 'Neuronal: GABAergic':
        return 'Inhibitory'

cell_data['CellType'] = cell_data.apply(define_cell_type, axis=1)

####################################
out_dir = f"/work/aliu10/AD_sc_Project/results/{sample_id}"

# Creating a cisTopic object for each sample
path_to_blacklist = "/work/aliu10/AD_sc_Project/scenicplus/pycisTopic/blacklist/hg38-blacklist.v2.bed"

cell_types = cell_data["CellType"].unique()

for cell_type in cell_types:
    print(f"Processing cell type: {cell_type}")
    
    # Subset metadata
    cell_data_subset = cell_data[cell_data["CellType"] == cell_type]
    
    # Construct path to celltype-specific consensus regions
    path_to_regions = os.path.join("/work/aliu10/AD_sc_Project/results/all/", f"consensus_regions_{cell_type}.bed")
    
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,  
        path_to_blacklist = path_to_blacklist,
        metrics = None,
        valid_bc = list(cell_data_subset.index),
        n_cpu = 64,
        project = sample_id
    )

    # Add metadata
    cistopic_obj.add_cell_data(cell_data_subset, split_pattern='___')

    # Save object
    output_path = os.path.join("/work/aliu10/AD_sc_Project/results", sample_id, f"cistopic_obj_{cell_type}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(cistopic_obj, f)
