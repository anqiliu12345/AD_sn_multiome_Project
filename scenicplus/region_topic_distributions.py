import pickle
import os
from pycisTopic.clust_vis import (find_clusters, run_umap, run_tsne)
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import topic_annotation
from pycisTopic.diff_features import (impute_accessibility, normalize_scores, find_highly_variable_features, find_diff_features)
import numpy as np
from pycisTopic.utils import region_names_to_coordinates
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity
import pandas as pd
import argparse

# Get cell type from command-line arguments
parser = argparse.ArgumentParser(description="Process cell type.")
parser.add_argument('cell_type', type=str, help="The cell type to process.")
args = parser.parse_args()
cell_type = args.cell_type
##############################
# Load object
pickle_path = os.path.join("/work/aliu10/AD_sc_Project/results/all", f"cistopic_obj_model_{cell_type}.pkl")
with open(pickle_path, "rb") as file:
    cistopic_obj = pickle.load(file)
    
cistopic_obj.projections = {"cell": {}, "region": {}}

##############################
# Topic binarization & QC
## Binarize region-topic distributions
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=False, num_columns=5
)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=False, num_columns=5
)

## Binarize cell-topic distributions (for subsequent usages)
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=False,
    num_columns=5, nbins=100)

os.makedirs(os.path.join("/work/aliu10/AD_sc_Project/results", "region_sets", cell_type), exist_ok = True)
os.makedirs(os.path.join("/work/aliu10/AD_sc_Project/results", "region_sets", cell_type, "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join("/work/aliu10/AD_sc_Project/results", "region_sets", cell_type, "Topics_otsu"), exist_ok = True)

## Save region-topic distributions table
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("/work/aliu10/AD_sc_Project/results", "region_sets", cell_type, "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("/work/aliu10/AD_sc_Project/results", "region_sets", cell_type, "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

## Save object
pickle.dump(
    cistopic_obj,
    open(os.path.join("/work/aliu10/AD_sc_Project/results/all", f"cistopic_obj_regions_{cell_type}.pkl"), "wb")
    )

