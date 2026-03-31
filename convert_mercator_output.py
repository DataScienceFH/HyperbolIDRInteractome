#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Convert Mercator .inf_coord to CSV.")
    p.add_argument(
        "--inf-coord",
        default=os.path.join("data", "STRING", "STRING_9606_H2.inf_coord"),
        help="Path to the Mercator .inf_coord output file.",
    )
    p.add_argument(
        "--out-csv",
        default=os.path.join("data", "STRING", "embedding.csv"),
        help="Path for the output CSV.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.inf_coord):
        raise FileNotFoundError(f"Mercator output not found: {args.inf_coord}")

    df = pd.read_csv(
        args.inf_coord,
        sep=r"\s+",
        comment="#",
        header=None,
        engine="python",
    )

    if df.shape[1] < 4:
        raise ValueError(
            f"Expected at least 4 columns in .inf_coord, found {df.shape[1]}"
        )

    # Mercator column order: node, kappa, theta, r_hyp, [...]
    col_names = ["node", "kappa", "theta", "r_hyp"] + [
        f"extra_{i}" for i in range(4, df.shape[1])
    ]
    df.columns = col_names

    out = pd.DataFrame({
        "node":  df["node"].astype(str),
        "r_hyp": pd.to_numeric(df["r_hyp"], errors="coerce"),
        "theta": pd.to_numeric(df["theta"], errors="coerce"),
    })

    out["rho"] = np.tanh(out["r_hyp"] / 2.0)
    out["x"]   = out["rho"] * np.cos(out["theta"])
    out["y"]   = out["rho"] * np.sin(out["theta"])

    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    out.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
