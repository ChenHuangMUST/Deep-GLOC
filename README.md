# Deep-GLOC

**Deep-GLOC** is a knowledge-guided graph framework for disentangling
cholesterol metabolic processes from transcriptomic data and deriving the
**free cholesterol loading score (FCLS)** in bulk and single-cell settings.

Cholesterol homeostasis is jointly regulated by multiple connected processes,
including **biosynthesis**, **uptake**, **esterification**, and **excretion**.
Conventional transcriptomic analyses often collapse these related processes
into a single signature or score, which limits biological interpretation.

Deep-GLOC was developed to address this problem. Starting from curated
process-specific seed genes, the framework constructs a liver-relevant
knowledge-augmented graph by integrating transcriptomic similarity and
protein-protein interaction information, expands process-specific signatures
through graph-based learning, and quantifies cholesterol loading states using
FCLS.

## Repository Structure

```text
Deep-GLOC/
├─ get_coexpression_network/
│  ├─ 01_get_G_liver.R
│  └─ 02_get_meta_graph.R
├─ run_GAT/
│  ├─ 01_prepare_augmented_graph_for_gat.ipynb
│  ├─ 02_train_final_model.ipynb
│  ├─ 03_generate_masked_seed_splits.ipynb
│  └─ 04_run_masked_seed_benchmark.ipynb
├─ demo/
│  ├─ run_demo.py
│  ├─ config_demo.json
│  └─ test_data/
├─ requirements.txt
├─ environment.yml
└─ README.md
```

## Contents

- Knowledge-guided graph modeling of cholesterol metabolism.
- Process disentanglement for four core cholesterol-related programs:
  Biosynthesis, Uptake, Esterification, and Excretion.
- Expansion of process-specific gene signatures from curated seed genes.
- Calculation of composite indices including Influx, Removal, and FCLS.
- Support for both bulk transcriptomic and single-cell RNA-seq analyses.
- Benchmarking against alternative graph-based methods such as DeepWalk and
  random walk with restart (RWR).

## Workflow

The overall Deep-GLOC workflow consists of the following steps:

1. Curate seed genes for each cholesterol metabolic process.
2. Define a liver-expressed gene universe.
3. Construct a liver knowledge-augmented network by integrating transcriptomic
   similarity and protein-protein interaction information.
4. Train the Deep-GLOC model to learn process-aware graph representations.
5. Expand process-specific gene signatures.
6. Compute process scores and derive FCLS in bulk and single-cell datasets.
7. Perform downstream validation, benchmarking, and biological interpretation.

## Installation

Create the Python environment:

```bash
conda env create -f environment.yml
conda activate deepgloc
```

Alternatively:

```bash
conda create -n deepgloc python=3.10
conda activate deepgloc
pip install -r requirements.txt
```

Install the main R packages used by the co-expression and downstream analysis
scripts:

```r
install.packages(c(
  "Seurat",
  "GSVA",
  "ggplot2",
  "dplyr",
  "data.table",
  "Matrix",
  "pheatmap",
  "survival",
  "survminer",
  "BiocManager"
))

BiocManager::install(c(
  "edgeR",
  "EnsDb.Hsapiens.v86",
  "AnnotationDbi",
  "biomaRt"
))
```

## Quick Demo

A small self-contained demo is provided in `demo/`. It uses toy data to verify
the expected input format and graph-based process scoring workflow.

Run from the repository root:

```bash
python demo/run_demo.py
```

The demo writes:

- `demo/results/demo_gene_process_scores.tsv`
- `demo/results/demo_sample_fcls.tsv`
- `demo/results/demo_expansion_candidates.tsv`

See `demo/README.md` for details.

## Main Analysis

The full study workflow is organized in two parts:

1. `get_coexpression_network/`: define the liver gene universe and construct
   the meta co-expression network.
2. `run_GAT/`: prepare the augmented graph, train the graph attention model,
   generate masked-seed benchmark splits, and run Deep-GLOC benchmarking.

The full analysis requires large public transcriptomic and interaction datasets
that are not included in this repository.

Example for defining the liver-expressed gene universe:

```bash
Rscript get_coexpression_network/01_get_G_liver.R \
  --input=/path/to/gene_reads_v11_liver.gct \
  --outdir=/path/to/output/G_liver
```

Example for building the meta co-expression graph:

```bash
Rscript get_coexpression_network/02_get_meta_graph.R \
  --g-liver-file=/path/to/G_liver_protein_coding_with_symbol.txt \
  --expr-dir=/path/to/harmonized_logCPM_files \
  --core-gene-file=/path/to/core_gene_all.csv \
  --outdir=/path/to/output/meta_graph
```

The full GAT notebooks in `run_GAT/` are provided as the study workflow used
for model preparation, training, and benchmarking. Before running them, update
the input and output paths to match your local project structure.

## Data Note

The files in `demo/test_data/` are synthetic toy data for software testing and
format demonstration only. They do not contain individual-level human data and
should not be used for biological interpretation.

## Full Reproducibility Data

The full manuscript-level reproducibility package contains large processed
inputs and outputs that are not suitable for direct GitHub upload. See
`reproducibility/` for the recommended data release layout, file manifest, data
availability template, and upload checklist.

After depositing the full package in Zenodo, Figshare, OSF, or an institutional
repository, add the DOI or stable URL here:

```text
Full reproducibility package: TODO
```

## Contact

For questions or feedback, please contact Chen Huang at chuang@must.edu.mo.
