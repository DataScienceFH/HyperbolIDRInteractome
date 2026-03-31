#!/usr/bin/env python3

import argparse
import csv
import os


def parse_args():
    p = argparse.ArgumentParser(description="Extract FuzDrop LLPS scores.")
    p.add_argument("--data-dir", required=True,
                   help="Root directory containing per-protein FuzDrop result subdirectories.")
    p.add_argument("--out-csv", default="llps_scores.csv",
                   help="Output CSV path (semicolon-delimited, default: llps_scores.csv).")
    return p.parse_args()


def extract_llps_score(path):
    """
    Return the LLPS propensity score from a protein_res.txt file, or None if not found.
    """
    target = "Liquid-liquid phase separation propensity p(LLPS)"
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if target in line:
                try:
                    value_str = line.split("=", 1)[1].strip().split()[0]
                    return float(value_str)
                except (IndexError, ValueError):
                    return None
    return None


def collect_scores(data_dir, out_csv):
    rows = []
    n_found = n_missing = 0

    for root, _dirs, files in os.walk(data_dir):
        if "protein_res.txt" not in files:
            continue
        txt_path = os.path.join(root, "protein_res.txt")
        string_id = os.path.basename(root)
        score = extract_llps_score(txt_path)
        if score is not None:
            rows.append((string_id, score))
            n_found += 1
        else:
            n_missing += 1

    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["string_id", "llps_score"])
        writer.writerows(rows)



def main():
    args = parse_args()
    if not os.path.isdir(args.data_dir):
        raise NotADirectoryError(f"FuzDrop data directory not found: {args.data_dir}")
    collect_scores(args.data_dir, args.out_csv)


if __name__ == "__main__":
    main()
