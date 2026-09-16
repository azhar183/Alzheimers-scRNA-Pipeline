import scanpy as sc
import pandas as pd

# 1. Setup
sc.settings.verbosity = 3
sc.set_figure_params(dpi=100)
sc.settings.figdir = './figures_comparison/' 

print("1. Loading Merged Data...")
adata = sc.read_h5ad('integrated_APOE3_APOE4.h5ad')

print("2. Re-computing Joint UMAP...")
# Dono patients ke cells ko ek sath mila kar naya UMAP ban banana
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)

print("3. Saving Comparison Graphs...")
# Graph 1: Condition (APOE3 vs APOE4)
# Graph 2: Leiden (Cell Types)
# wspace=0.5 dono graphs ke beech extra space add kar dega
sc.pl.umap(adata, color=['Condition', 'leiden'], wspace=0.5, save='_patients_vs_clusters.png', show=False)

print("4. Finding Alzheimer's Disease Genes (DGE)...")
# Computer ko bolna ki dono conditions ko compare kare
sc.tl.rank_genes_groups(adata, groupby='Condition', method='wilcoxon')

print("5. Saving Genes to Excel/CSV...")
# Results ko ek proper table mein convert karke CSV mein save karna
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

df = pd.DataFrame(
    {group + '_' + key: result[key][group]
    for group in groups for key in ['names', 'pvals', 'logfoldchanges']}
)
df.to_csv('Alzheimers_Disease_Genes.csv', index=False)

print("✅ Analysis Complete! You are now a Bioinformatician!")