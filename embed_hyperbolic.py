#!/usr/bin/env python3

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import networkx as nx


# ── Argument parsing ──────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Hyperbolic PSO embedding of STRING LCC.")
    p.add_argument("--string-dir", default=os.path.join("data", "STRING"),
                   help="Directory containing the decompressed STRING links file.")
    p.add_argument("--score-min", type=int, default=900,
                   help="Minimum combined_score threshold (default: 900).")
    p.add_argument("--out-csv", default=None,
                   help="Output CSV path. Defaults to <string-dir>/embedding.csv.")
    p.add_argument("--epochs", type=int, default=20,
                   help="Number of gradient-descent epochs (default: 20).")
    p.add_argument("--seed", type=int, default=1,
                   help="Random seed for reproducibility (default: 1).")
    p.add_argument("--temperature", type=float, default=0.6,
                   help="Temperature parameter T for the PSO model (default: 0.6).")
    p.add_argument("--lr", type=float, default=0.01,
                   help="Angular learning rate (default: 0.01).")
    p.add_argument("--neg-samples", type=int, default=5,
                   help="Negative samples per positive edge per epoch (default: 5).")
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────
def find_string_links_file(dir_path):
    pattern = os.path.join(dir_path, "*protein.links*.txt*")
    candidates = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
    candidates = [c for c in candidates if not c.lower().endswith(".gz")]
    if not candidates:
        raise FileNotFoundError(
            f"No decompressed STRING links file found in: {dir_path}"
        )
    return os.path.normpath(candidates[0])


def extract_string_metadata(file_path):
    fname = os.path.basename(file_path)
    species_m = re.match(r"^(\d+)\.", fname)
    if not species_m:
        raise ValueError(f"Cannot extract species ID from filename: {fname}")
    version_m = re.search(r"v([0-9.]+)", fname)
    return {
        "species": species_m.group(1),
        "version": version_m.group(1) if version_m else None,
    }


def hyperbolic_distance(r, theta, i, j):
    ri, rj = r[i], r[j]
    dt = np.pi - abs(np.pi - abs(theta[i] - theta[j]))
    arg = np.cosh(ri) * np.cosh(rj) - np.sinh(ri) * np.sinh(rj) * np.cos(dt)
    return np.arccosh(max(arg, 1.0))


# ── Main ──────────────────────────────────────────────────────
def main():
    args = parse_args()
    np.random.seed(args.seed)

    string_file = find_string_links_file(args.string_dir)
    meta = extract_string_metadata(string_file)
    out_csv = args.out_csv or os.path.join(args.string_dir, "embedding.csv")

    df = pd.read_csv(string_file, sep=r"\s+")
    df = df[df["combined_score"] >= args.score_min].copy()

    tax_pat = rf"^{re.escape(meta['species'])}\."
    df["protein1"] = df["protein1"].astype(str).str.replace(tax_pat, "", regex=True)
    df["protein2"] = df["protein2"].astype(str).str.replace(tax_pat, "", regex=True)

    edges = df[["protein1", "protein2"]].drop_duplicates()
    edges = edges[edges["protein1"] != edges["protein2"]]

    G = nx.from_pandas_edgelist(edges, "protein1", "protein2")
    lcc = max(nx.connected_components(G), key=len)
    G = G.subgraph(lcc).copy()

    # ── Radial coordinates (PSO model) ────────────────────────
    N = G.number_of_nodes()
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}

    deg = np.array([G.degree(n) for n in nodes], dtype=float)
    order = np.argsort(-deg, kind="mergesort")
    rank = np.empty_like(order)
    rank[order] = np.arange(N) + 1  # 1-based rank

    r = 2.0 * np.log(rank)
    r = (r - r.min()) / np.percentile(r, 95) * 8.0

    # ── Angular coordinate optimisation ───────────────────────
    theta = np.random.uniform(0.0, 2.0 * np.pi, size=N)
    adj = {node_idx[n]: [node_idx[m] for m in G.neighbors(n)] for n in nodes}
    all_idx = np.arange(N)

    for epoch in range(args.epochs):
        total_ll = 0.0
        for i in range(N):
            for j in adj[i]:
                d = hyperbolic_distance(r, theta, i, j)
                p = 1.0 / (1.0 + np.exp((d - 2.0 * np.log(N)) / (2.0 * args.temperature)))
                theta[i] -= args.lr * (1.0 - p) * np.sign(np.sin(theta[i] - theta[j]))
                total_ll += np.log(p + 1e-12)

            negs = np.random.choice(all_idx, args.neg_samples, replace=False)
            for j in negs:
                if j in adj[i] or j == i:
                    continue
                d = hyperbolic_distance(r, theta, i, j)
                p = 1.0 / (1.0 + np.exp((d - 2.0 * np.log(N)) / (2.0 * args.temperature)))
                theta[i] -= args.lr * (-p) * np.sign(np.sin(theta[i] - theta[j]))
                total_ll += np.log(1.0 - p + 1e-12)

    theta = np.mod(theta, 2.0 * np.pi)

    # ── Convert to Poincaré disk ──────────────────────────────
    rho = np.tanh(r / 2.0)
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)

    # ── Save ──────────────────────────────────────────────────
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    coords = pd.DataFrame({
        "node":  nodes,
        "r_hyp": r,
        "theta": theta,
        "rho":   rho,
        "x":     x,
        "y":     y,
        "degree": deg.astype(int),
    })
    coords.to_csv(out_csv, index=False)


if __name__ == "__main__":
    main()
