import os
import scanpy as sc
import pandas as pd
import anndata as ad
import flowsig as fs
import numpy as np

case_group = "AD_female"
control_group = "CT_female"
compare_groups = "female_AD_vs_CT"

# Load the two AnnData files
adata1 = sc.read_h5ad("/work/aliu10/AD_sc_Project/results/" + case_group + "/adata_" + case_group + ".h5ad")
adata2 = sc.read_h5ad("/work/aliu10/AD_sc_Project/results/" + control_group + "/adata_" + control_group + ".h5ad")
# Concatenate along observations (cells)
adata_merged = adata1.concatenate(adata2, batch_key="Status", batch_categories=[case_group, control_group])

condition_key = 'Status'

adata = adata_merged
# Both the index name and column names in adata.var cannot have the same names

data_dir = '/work/aliu10/AD_sc_Project/results/CCI/'
cellchat_IFNg = pd.read_csv(data_dir + 'seaAD_communications_case_' + compare_groups + '.csv')
cellchat_Ctrl = pd.read_csv(data_dir + 'seaAD_communications_control_' + compare_groups + '.csv')

cellchat_output_key = 'cellchat_output'
# Make sure your keys for these align with their condition labels
adata.uns[cellchat_output_key] = {'Ctrl': cellchat_Ctrl,
                                  'IFNg': cellchat_IFNg}

adata.obs_names = adata.obs_names.str.replace(r'-(AD|CT)_(male|female)$', '', regex=True)

###################################################################################################################
cell_types = ['Astrocyte', 'Microglia-PVM', 'Excitatory', 'Inhibitory', 'Oligodendrocyte', 'OPC']

# -------- Loop over each cell type -------- #
for cell_type in cell_types:
    print(f"Running FlowSig for cell type: {cell_type}")

    # Step 1: Load Gene × GEM (kME) matrices
    gene_gem_1 = pd.read_csv("/work/aliu10/AD_sc_Project/results/hdWGCNA/GEM/GeneByGEM_" + cell_type + "_" + compare_groups + "_case.csv", index_col=0)
    gene_gem_2 = pd.read_csv("/work/aliu10/AD_sc_Project/results/hdWGCNA/GEM/GeneByGEM_" + cell_type + "_" + compare_groups + "_control.csv", index_col=0)
    # Take absolute value to reflect "importance" regardless of direction
    gene_gem_1 = gene_gem_1.abs()
    gene_gem_2 = gene_gem_2.abs()
    union_gems = gene_gem_1.columns.union(gene_gem_2.columns)

    W_1 = gene_gem_1
    W_2 = gene_gem_2

    # Step 2: Normalize W to V (L2 normalization by column)
    from sklearn.preprocessing import normalize

    V_1 = normalize(W_1.values, axis=0)  # shape: genes × gems
    V_2 = normalize(W_2.values, axis=0)

    # Step 3: Load Cell x GEM (hME) for both conditions
    cell_gem_1 = pd.read_csv("/work/aliu10/AD_sc_Project/results/hdWGCNA/GEM/CellByGEM_" + cell_type + "_" + compare_groups + "_case.csv", index_col=0)
    cell_gem_2 = pd.read_csv("/work/aliu10/AD_sc_Project/results/hdWGCNA/GEM/CellByGEM_" + cell_type + "_" + compare_groups + "_control.csv", index_col=0)
    # Take absolute value to reflect "importance" regardless of direction
    cell_gem_1 = cell_gem_1.abs()
    cell_gem_2 = cell_gem_2.abs()

    # Reindex both matrices to align with full set
    H_1 = cell_gem_1.reindex(columns = gene_gem_1.columns)
    H_2 = cell_gem_2.reindex(columns = gene_gem_2.columns)

    # Step 4: Prepare X_gem matrix (cells × union_gems)
    H_combined = pd.concat([H_1, H_2], axis=1)
    idx = adata.obs_names.intersection(H_combined.index)
    adata_updated = adata[adata.obs_names.isin(idx)].copy()
    H_combined = H_combined.loc[idx]

    # Scale so that the GEM memberships sum to 1 per cell (flowsig idea)
    gem_sum = H_combined.sum(axis=0)
    H_combined = H_combined / gem_sum

    # Save to AnnData object
    adata_updated.obsm["X_gem"] = H_combined
    adata_updated.uns["pyliger_info"] = {
        case_group: {
            "H": H_1,
            "W": W_1.values,
            "V": V_1
        },
        control_group: {
            "H": H_2,
            "W": W_2.values,
            "V": V_2
        },
        "vars": W_1.index.to_numpy(),
        "n_gems": H_combined.shape[1],
        "gems": H_combined.columns.to_numpy()
    }
    ###################################################################################################################
    # We first construct the potential cellular flows from the cellchat output, i.e., separate the inflows from the outflows.
    config = fs.pp.FlowSigConfig(
        gem_expr_key='X_gem',       
        scale_gem_expr=True,
        flowsig_network_key = 'flowsig_network',
        flowsig_expr_key = 'X_flow'
    )

    fs.pp.construct_flows_from_cellchat(
        adata_updated,
        cellchat_output_key,
        model_organism = 'human',
        tfs_to_use = None,
        config = config,
        construction = 'v2'
    )

    # Then we subset for "differentially flowing" variables, using a Mann-Whitney U test on the inflow and outflow expressions, separately.
    fs.pp.determine_informative_variables(adata_updated,  
                                        config=config,
                                        spatial = False,
                                        condition_key = condition_key,
                                        control=control_group,
                                        qval_threshold = 0.05,
                                        logfc_threshold = 0.1, # 0.5
                                        construction = 'v2')

    # Now we are ready to learn the network
    fs.tl.learn_intercellular_flows(adata_updated,
                            config=config,
                            condition_key = condition_key,
                            control = control_group, 
                            alpha_ci = 1e-3,
                            alpha_inv = 1e-3,
                            use_spatial = False,
                            n_jobs = 32,
                            n_bootstraps = 100)


    # Now we do post-learning validation to reorient the network and remove low-quality edges.
    # This part is key for reducing false positives 
    fs.tl.apply_biological_flow(adata_updated,
                            flowsig_network_key = 'flowsig_network',
                            adjacency_key = 'adjacency',
                            validated_key = 'validated')

    edge_threshold = 0.5

    fs.tl.filter_low_confidence_edges(adata_updated,
                                    edge_threshold = edge_threshold,
                                    flowsig_network_key = 'flowsig_network',
                                    adjacency_key = 'adjacency',
                                    filtered_key = 'filtered')


    output_dir = '/work/aliu10/AD_sc_Project/results/FlowSig/' + compare_groups + '_GEM/' + cell_type
    os.makedirs(output_dir, exist_ok=True)  
    adata_updated.write(output_dir + '/seaAD.h5ad', compression='gzip')
