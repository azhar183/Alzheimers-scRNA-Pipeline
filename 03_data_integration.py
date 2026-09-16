import scanpy as sc
import pandas as pd

# 1. Setup
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sc.settings.figdir = './figures_integrated/' # Naya folder merged graphs ke liye

print("1. Loading processed datasets...")
# Load APOE3 and APOE4 processed objects
adata_apoe3 = sc.read_h5ad('sample_APOE3_processed.h5ad')
adata_apoe4 = sc.read_h5ad('sample_APOE4_processed.h5ad')

# Identify which cells belong to which patient
adata_apoe3.obs['Condition'] = 'APOE3_Normal'
adata_apoe4.obs['Condition'] = 'APOE4_HighRisk'

print("2. Merging Datasets...")
# Combine them into one big dataset
adata_merged = sc.concat([adata_apoe3, adata_apoe4], label='dataset')

print(f"Merge Complete! Total Cells: {adata_merged.n_obs}, Total Genes: {adata_merged.n_vars}")

# 3. Save the merged dataset
adata_merged.write('integrated_APOE3_APOE4.h5ad')
print("✅ Merged data saved as 'integrated_APOE3_APOE4.h5ad'")