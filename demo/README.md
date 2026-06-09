# Deep-GLOC Demo

This folder contains a minimal, self-contained demo for checking the Deep-GLOC
input format and graph-based cholesterol process scoring workflow.

The demo uses a small toy graph and synthetic expression matrix. It is designed
to run quickly on a normal laptop and does not reproduce the full trained GAT
model or the full study results.

## Files

- `run_demo.py`: one-command demo runner.
- `config_demo.json`: paths and demo parameters.
- `test_data/node_table_demo.tsv`: gene node table.
- `test_data/edge_list_demo.tsv`: weighted gene-gene edge list.
- `test_data/core_gene_demo.csv`: process-specific seed genes.
- `test_data/expression_demo.tsv`: toy bulk expression matrix.
- `expected_output/`: expected result files generated from the demo inputs.
- `results/`: generated outputs when running the demo.

## Install

From the repository root:

```bash
conda create -n deepgloc python=3.10
conda activate deepgloc
pip install -r requirements.txt
```

The demo itself only requires `numpy`, `pandas`, and `networkx`.

## Run

From the repository root:

```bash
python demo/run_demo.py
```

Or from inside this folder:

```bash
python run_demo.py
```

## Outputs

Running the demo writes the following files to `demo/results/`:

- `demo_gene_process_scores.tsv`: propagated gene scores for the four
  cholesterol-related processes.
- `demo_sample_fcls.tsv`: sample-level Biosynthesis, Uptake, Esterification,
  Excretion, Influx, Removal, and FCLS scores.
- `demo_expansion_candidates.tsv`: top non-seed candidate genes for each
  process in the toy graph.

The expected outputs in `expected_output/` can be used to confirm that the demo
ran successfully.

## Input Format

`node_table_demo.tsv`

```text
gene_id    gene_symbol
```

`edge_list_demo.tsv`

```text
source    target    weight    source_type
```

`core_gene_demo.csv`

```text
gene_id,gene_symbol,process
```

`expression_demo.tsv`

```text
gene_id    gene_symbol    Sample_1    Sample_2    ...
```

Supported process labels are `Biosynthesis`, `Uptake`, `Esterification`, and
`Excretion`.
