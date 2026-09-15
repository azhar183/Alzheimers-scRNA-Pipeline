import scanpy as sc

# 1. Setup
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sc.settings.figdir = './figures_apoe4/' 

print("1. Loading APOE4 Data...")
adata = sc.read_10x_mtx('sample_APOE4/', var_names='gene_symbols', cache=True)

# 2. QC & Filtering
print("2. Filtering low-quality cells...")
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :] 

# 3. Normalization & Log Transform
print("3. Normalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# 4. Math: PCA, Neighbors, and UMAP
print("4. Running PCA & UMAP (Please wait 10-20 seconds)...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)

# 5. Save the Output
print("5. Saving UMAP plot and Data...")
sc.pl.umap(adata, color=['leiden'], save='_apoe4_clusters.png', show=False)
adata.write('sample_APOE4_processed.h5ad')

print("Pipeline 100% Complete! Check the 'figures_apoe4' folder.")