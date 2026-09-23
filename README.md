🧬 Alzheimer's Disease scRNA-seq Pipeline (APOE3 vs APOE4)
An end-to-end automated Single-Cell RNA Sequencing (scRNA-seq) bioinformatics pipeline built in Python. This project investigates the cellular and molecular differences between a normal brain (APOE3) and a high-risk Alzheimer's brain (APOE4) using computational biology and unsupervised machine learning.

🎯 Project Overview
The goal of this project is to process raw 10x Genomics scRNA-seq data and build a modular, reproducible pipeline to identify disease-specific biomarkers.

Key Biological Findings:

Distinct Cell States: Joint UMAP visualization reveals that APOE4 (Alzheimer's) cells form entirely distinct clusters compared to APOE3 (Normal) cells, indicating massive transcriptomic shifts.

Biomarker Discovery: Successfully performed Differential Gene Expression (DGE) using the Wilcoxon rank-sum test to extract the Top 10 disease-driving genes.

🛠️ Technologies & Tools Used
Language: Python 3

Bioinformatics Core: Scanpy, AnnData

Data Science: Pandas, NumPy

Machine Learning: PCA, UMAP, Leiden Clustering

Environment: VS Code, Jupyter Notebooks, Conda, Git

📂 Pipeline Architecture
This project transitions from exploratory Jupyter notebooks to professional, automated Python scripts for scalable computation:

01 & 02 - Data QC & Preprocessing:
Loads raw matrices, filters dead cells/mitochondrial genes, normalizes total counts, and extracts highly variable genes for both APOE3 and APOE4 samples individually.

03_data_integration.py:
Combines the processed patients' datasets into a single highly-dimensional matrix (17,148 cells and 2,590 genes). Handles unique barcode generation.

04_comparative_analysis.py:
Runs PCA, computes the neighborhood graph, and performs Leiden clustering. Executes DGE (Wilcoxon test) to mathematically rank genes comparing APOE3 vs. APOE4, exporting results to CSV.

05_final_analysis.ipynb:
Interactive notebook for high-quality data visualization (Dot Plots, Joint UMAPs) ready for research publications.

🚀 How to Run This Project
1. Setup the Environment:
Ensure you have Conda installed, then create and activate the environment:
conda create -n scanpy_env python=3.11
conda activate scanpy_env

2. Install Dependencies:
pip install scanpy pandas numpy jupyter igraph leidenalg

3. Execute the Pipeline:
Run the scripts in numerical order:
python 03_data_integration.py
python 04_comparative_analysis.py

(Note: Raw .h5ad data files are not uploaded due to GitHub size limits. You will need the raw GSE datasets to run the preprocessing scripts).

📊 Results Visualization
All generated graphs (UMAPs, Violin Plots for QC, and DGE Dot Plots) are automatically saved in the figures_comparison/ directory. The Differential Gene Expression results are exported as Alzheimers_Disease_Genes.csv.

Author: [Apna Naam Yahan Likhein]
Status: Completed
