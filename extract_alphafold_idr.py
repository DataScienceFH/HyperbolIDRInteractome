#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Extract IDR percentages from AlphaFold CIF files.")
    p.add_argument("--cif-dir", required=True,
                   help="Directory containing AlphaFold .cif files.")
    p.add_argument("--map-csv", required=True,
                   help="Semicolon-delimited STRING_id;UniProt_id mapping CSV.")
    p.add_argument("--out-csv", required=True,
                   help="Output CSV path (semicolon-delimited).")
    p.add_argument("--cutoff", type=float, default=50.0,
                   help="pLDDT threshold below which a residue is considered disordered (default: 50.0).")
    p.add_argument("--cif-version", type=int, default=6,
                   help="AlphaFold model version number in the CIF filename (default: 6).")
    return p.parse_args()


def load_mapping(path: Path):
    """Load STRING → UniProt mapping from a semicolon-delimited CSV."""
    mapping = []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=";")
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            string_id, uniprot_id = row[0].strip(), row[1].strip()
            if string_id and uniprot_id:
                mapping.append((string_id, uniprot_id))
    return mapping


def disorder_pct_from_cif(path: Path, cutoff: float):
    """
    Parse an AlphaFold mmCIF file and return the percentage of residues
    with pLDDT < cutoff, or None if the relevant loop is absent.

    The function looks for the _ma_qa_metric_local loop (7 columns) and
    extracts column index 4 (metric_value = per-residue pLDDT).
    """
    text = path.read_text(encoding="utf-8")
    target_cols = [
        "_ma_qa_metric_local.label_asym_id",
        "_ma_qa_metric_local.label_comp_id",
        "_ma_qa_metric_local.label_seq_id",
        "_ma_qa_metric_local.metric_id",
        "_ma_qa_metric_local.metric_value",
        "_ma_qa_metric_local.model_id",
        "_ma_qa_metric_local.ordinal_id",
    ]

    vals = []
    in_loop = False
    cols = []

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")

        if line.startswith("loop_"):
            in_loop = True
            cols = []
            continue

        if in_loop and line.startswith("_"):
            cols.append(line.strip())
            continue

        if in_loop and line and not line.startswith("_"):
            parts = line.split()
            if len(cols) >= 7 and cols[:7] == target_cols and len(parts) >= 5:
                try:
                    vals.append(float(parts[4]))
                except ValueError:
                    pass
            continue

        if in_loop and line.startswith("data_"):
            in_loop = False

    if not vals:
        return None

    disordered = sum(1 for v in vals if v < cutoff)
    return 100.0 * disordered / len(vals)


def main():
    args = parse_args()
    cif_dir = Path(args.cif_dir)
    map_csv = Path(args.map_csv)
    out_csv = Path(args.out_csv)

    if not cif_dir.is_dir():
        raise NotADirectoryError(f"CIF directory not found: {cif_dir}")
    if not map_csv.is_file():
        raise FileNotFoundError(f"Mapping CSV not found: {map_csv}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    mapping = load_mapping(map_csv)

    found = missing = 0
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["STRING", "UniProt", "disorder_pct"])

        for string_id, uniprot_id in mapping:
            cif_name = f"AF-{uniprot_id}-F1-model_v{args.cif_version}.cif"
            cif_path = cif_dir / cif_name

            if cif_path.exists():
                found += 1
                pct = disorder_pct_from_cif(cif_path, args.cutoff)
                writer.writerow([
                    string_id,
                    uniprot_id,
                    f"{pct:.6f}" if pct is not None else "",
                ])
            else:
                missing += 1
                writer.writerow([string_id, uniprot_id, ""])



if __name__ == "__main__":
    main()
