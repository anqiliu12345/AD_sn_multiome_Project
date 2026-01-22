### Peak calling (and quality control)

import os
import pandas as pd
import scanpy as sc
import pycisTopic
import argparse

##############################
# Prepare public data
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes = chromsizes[~chromsizes["Chromosome"].str.contains("chrUn|chrX|chrY|chrM")] # Remove chr X, Y, M and unknown

######################################################################
# Define the path to fragments and data directories
data_path = "/work/aliu10/AD_sc_Project/sea_AD_data/MTG/multiome/processed_data"

# Get donor_id from command-line arguments
parser = argparse.ArgumentParser(description="Process donor ID.")
parser.add_argument('donor_id', type=str, help="The donor ID to process.")
args = parser.parse_args()
donor_id = args.donor_id

# Construct file paths
fragment_file_path = os.path.join(data_path, f"{donor_id}_filtered_fragments.tsv.gz")

print(f"Load donor {donor_id}...")
fragments_dict = {donor_id: fragment_file_path}

# Prepare output directory for each donor
out_dir = f"/work/aliu10/AD_sc_Project/results/{donor_id}"
os.makedirs(out_dir, exist_ok=True)

metadata_file = os.path.join(data_path, "metadata", f"metadata_{donor_id}.csv")
cell_data = pd.read_csv(metadata_file, index_col=0) 

def define_cell_type(row):
    if row['Subclass'] in ['Oligodendrocyte', 'Astrocyte', 'Microglia-PVM', 'OPC', 'Endothelial', 'VLMC']:
        return row['Subclass']
    elif row['Class'] == 'Neuronal: Glutamatergic':
        return 'Excitatory'
    elif row['Class'] == 'Neuronal: GABAergic':
        return 'Inhibitory'

cell_data['CellType'] = cell_data.apply(define_cell_type, axis=1)

######################################################################
# Pseudobulk files
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

os.makedirs("tmp", exist_ok = True)

bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "CellType",
    sample_id_col = "donor_id",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 64, 
    normalize_bigwig = True,
    temp_dir = "/tmp/"
)

bed_paths = {key.replace('/', '_').replace(' ', '_'): value for key, value in bed_paths.items()}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")


######################################################################
# Peak calling
bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

# Call peaks for each pseudobulk fragments.tsv.gz file.
from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)
os.makedirs("/tmp/", exist_ok=True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    skip_empty_peaks=True,  # Allow samples with no peaks to be skipped
    genome_size = 'hs',
    n_cpu = 64, 
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = "/tmp/"
)
