#!/usr/bin/env python
"""Run a small Deep-GLOC demo on toy cholesterol-metabolism data.

The demo checks the expected tabular inputs, performs a lightweight
graph-based propagation from process-specific seed genes, and calculates
sample-level influx, removal, and FCLS scores. It is intended for format and
workflow validation, not for reproducing the full trained GAT model.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


PROCESSES = ["Biosynthesis", "Uptake", "Esterification", "Excretion"]


def resolve_path(base_dir: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base_dir / path


def load_config(config_file: Path) -> dict:
    with config_file.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config["_base_dir"] = config_file.parent
    return config


def read_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base_dir = config["_base_dir"]
    node_df = pd.read_csv(resolve_path(base_dir, config["node_file"]), sep="\t")
    edge_df = pd.read_csv(resolve_path(base_dir, config["edge_file"]), sep="\t")
    core_df = pd.read_csv(resolve_path(base_dir, config["core_gene_file"]))
    expr_df = pd.read_csv(resolve_path(base_dir, config["expression_file"]), sep="\t")

    required = {
        "node_table": (node_df, {"gene_id", "gene_symbol"}),
        "edge_list": (edge_df, {"source", "target", "weight"}),
        "core_gene": (core_df, {"gene_id", "gene_symbol", "process"}),
        "expression": (expr_df, {"gene_id", "gene_symbol"}),
    }
    for name, (df, columns) in required.items():
        missing = columns.difference(df.columns)
        if missing:
            raise ValueError(f"{name} is missing required columns: {sorted(missing)}")

    unknown_process = sorted(set(core_df["process"]) - set(PROCESSES))
    if unknown_process:
        raise ValueError(f"Unknown process labels in core gene file: {unknown_process}")

    return node_df, edge_df, core_df, expr_df


def build_graph(node_df: pd.DataFrame, edge_df: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    for row in node_df.itertuples(index=False):
        graph.add_node(row.gene_id, gene_symbol=row.gene_symbol)

    for row in edge_df.itertuples(index=False):
        if row.source in graph and row.target in graph:
            graph.add_edge(row.source, row.target, weight=float(row.weight))

    if graph.number_of_edges() == 0:
        raise ValueError("The demo graph has no valid edges after matching node IDs.")
    return graph


def propagate_process_scores(
    graph: nx.Graph,
    core_df: pd.DataFrame,
    alpha: float,
    iterations: int,
) -> pd.DataFrame:
    genes = list(graph.nodes)
    index = {gene: i for i, gene in enumerate(genes)}
    scores = pd.DataFrame(0.0, index=genes, columns=PROCESSES)

    adjacency = nx.to_numpy_array(graph, nodelist=genes, weight="weight", dtype=float)
    degree = adjacency.sum(axis=1)
    transition = np.divide(
        adjacency,
        degree[:, None],
        out=np.zeros_like(adjacency),
        where=degree[:, None] > 0,
    )

    for process in PROCESSES:
        seed_genes = core_df.loc[core_df["process"] == process, "gene_id"]
        seed_vector = np.zeros(len(genes), dtype=float)
        for gene in seed_genes:
            if gene in index:
                seed_vector[index[gene]] = 1.0

        if seed_vector.sum() == 0:
            continue

        seed_vector = seed_vector / seed_vector.sum()
        propagated = seed_vector.copy()
        for _ in range(iterations):
            propagated = alpha * transition.T.dot(propagated) + (1 - alpha) * seed_vector

        if propagated.max() > 0:
            propagated = propagated / propagated.max()
        scores[process] = propagated

    scores.index.name = "gene_id"
    return scores.reset_index()


def score_samples(
    expr_df: pd.DataFrame,
    gene_scores: pd.DataFrame,
    node_df: pd.DataFrame,
) -> pd.DataFrame:
    expr = expr_df.set_index("gene_id").drop(columns=["gene_symbol"])
    expr = expr.apply(pd.to_numeric)
    z_expr = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0)
    z_expr = z_expr.fillna(0.0)

    weights = gene_scores.set_index("gene_id")[PROCESSES]
    common = z_expr.index.intersection(weights.index)
    process_scores = {}
    for process in PROCESSES:
        w = weights.loc[common, process]
        if w.sum() == 0:
            process_scores[process] = pd.Series(0.0, index=z_expr.columns)
        else:
            process_scores[process] = z_expr.loc[common].mul(w, axis=0).sum(axis=0) / w.sum()

    sample_df = pd.DataFrame(process_scores)
    sample_df["Influx"] = sample_df["Biosynthesis"] + sample_df["Uptake"]
    sample_df["Removal"] = sample_df["Esterification"] + sample_df["Excretion"]
    sample_df["FCLS"] = sample_df["Influx"] - sample_df["Removal"]
    sample_df.index.name = "sample_id"
    return sample_df.reset_index()


def select_candidates(
    gene_scores: pd.DataFrame,
    core_df: pd.DataFrame,
    node_df: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:
    symbols = node_df.set_index("gene_id")["gene_symbol"]
    seed_genes = set(core_df["gene_id"])
    records = []

    for process in PROCESSES:
        ranked = gene_scores.loc[~gene_scores["gene_id"].isin(seed_genes), ["gene_id", process]]
        ranked = ranked.sort_values(process, ascending=False).head(top_n)
        for rank, row in enumerate(ranked.itertuples(index=False), start=1):
            records.append(
                {
                    "process": process,
                    "rank": rank,
                    "gene_id": row.gene_id,
                    "gene_symbol": symbols.get(row.gene_id, row.gene_id),
                    "demo_score": getattr(row, process),
                }
            )

    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Deep-GLOC toy demo.")
    parser.add_argument(
        "--config",
        default=Path(__file__).with_name("config_demo.json"),
        type=Path,
        help="Path to the demo JSON config file.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    node_df, edge_df, core_df, expr_df = read_inputs(config)
    graph = build_graph(node_df, edge_df)

    gene_scores = propagate_process_scores(
        graph,
        core_df,
        alpha=float(config.get("alpha", 0.85)),
        iterations=int(config.get("iterations", 30)),
    )
    gene_scores = gene_scores.merge(node_df, on="gene_id", how="left")
    gene_scores = gene_scores[["gene_id", "gene_symbol"] + PROCESSES]

    sample_scores = score_samples(expr_df, gene_scores, node_df)
    candidates = select_candidates(
        gene_scores,
        core_df,
        node_df,
        top_n=int(config.get("top_candidates_per_process", 3)),
    )

    output_dir = resolve_path(config["_base_dir"], config.get("output_dir", "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_scores.to_csv(output_dir / "demo_gene_process_scores.tsv", sep="\t", index=False)
    sample_scores.to_csv(output_dir / "demo_sample_fcls.tsv", sep="\t", index=False)
    candidates.to_csv(output_dir / "demo_expansion_candidates.tsv", sep="\t", index=False)

    print("Deep-GLOC demo completed.")
    print(f"Nodes: {graph.number_of_nodes()}  Edges: {graph.number_of_edges()}")
    print(f"Results written to: {output_dir}")
    print(sample_scores[["sample_id", "Influx", "Removal", "FCLS"]].to_string(index=False))


if __name__ == "__main__":
    main()
