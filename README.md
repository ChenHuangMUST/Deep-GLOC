# Deep-GLOC
Knowledge-guided graph framework for disentangling cholesterol metabolic processes and deriving the free cholesterol loading score (FCLS) from bulk and single-cell transcriptomic data.

**Deep-GLOC** is a knowledge-guided graph framework for disentangling cholesterol metabolic processes from transcriptomic data and deriving the **free cholesterol loading score (FCLS)** in bulk and single-cell settings.

Cholesterol homeostasis is jointly regulated by multiple tightly connected processes, including **biosynthesis**, **uptake**, **esterification**, and **excretion**. Conventional transcriptomic analyses often collapse these related processes into a single signature or score, which limits biological interpretation.
Deep-GLOC was developed to address this problem. Starting from curated process-specific seed genes, the framework constructs a liver-relevant knowledge-augmented graph by integrating transcriptomic similarity and protein–protein interaction information, and then expands process-specific signatures through graph-based learning. Based on the resulting signatures, Deep-GLOC further quantifies cholesterol loading states using the **free cholesterol loading score (FCLS)**.


**Contents**

 - runICA/: This folder contains the scripts for running ICA.
 - getimmIC.R: Filtering for independent components related to immunology.
 - getiDICmatrix.R: Calculation of TIME-driver independent components profile by integrating signature components and somatic mutation data.
 - getiDICss.R: Development of iDIC-based risk score system.
 - getMLcombine.R：A function that develops robust predictive models with excellent performance, using patient prognosis as the outcome.
 - plot.heatmap.R：Plot the heatmap of the C-index of all combined prediction models.
 - compare.performance.R：Calculating the performance measurements of all combined prediction models.- Knowledge-guided graph modeling of cholesterol metabolism
- Process disentanglement for four core cholesterol-related programs:
  - Biosynthesis
  - Uptake
  - Esterification
  - Excretion
- Expansion of process-specific gene signatures from curated seed genes
- Calculation of composite indices including:
  - FCLS
- Support for both bulk transcriptomic and single-cell RNA-seq data
- Benchmarking against alternative graph-based methods such as DeepWalk and random walk with restart (RWR)

The overall Deep-GLOC workflow consists of the following steps:

1. Curate seed genes for each cholesterol metabolic process  
2. Define a liver-expressed gene universe  
3. Construct a liver knowledge-augmented network by integrating:
   - transcriptomic similarity
   - protein–protein interaction information  
4. Train the Deep-GLOC model to learn process-aware graph representations  
5. Expand process-specific gene signatures  
6. Compute process scores and derive FCLS in bulk and single-cell datasets  
7. Perform downstream validation, benchmarking, and biological interpretation  

**Installation**

#R
install.packages(c(
  "Seurat",
  "GSVA",
  "ggplot2",
  "dplyr",
  "data.table",
  "Matrix",
  "pheatmap",
  "survival",
  "survminer"
))
#Python
conda create -n deepgloc python=3.10
conda activate deepgloc
pip install numpy pandas scipy scikit-learn networkx torch

**Contact**

For any questions or feedback, please reach out to Chen Huang at chuang@must.edu.mo.


