# HyperbolIDRInteractome

A modular analysis pipeline for studying **intrinsically disordered regions (IDRs)** in the human proteome through the lens of **hyperbolic network geometry**. The pipeline embeds the STRING human protein–protein interaction (PPI) network in two-dimensional hyperbolic space (H²) and integrates IDR content, liquid–liquid phase separation (LLPS) propensity, network topology metrics, and community structure to test whether the position of a protein in hyperbolic space is a predictor of its disorder properties and interaction behaviour.

---

## Scientific Background

Scale-free networks like PPI networks have a latent hyperbolic geometry: highly connected hub proteins occupy the centre of the Poincaré disk, while peripheral proteins reside at the boundary. Hyperbolic distance between two proteins therefore summarises both their structural similarity and their tendency to interact. Intrinsically disordered regions – sequence stretches lacking a fixed 3D fold – are enriched in hub proteins and drive LLPS-mediated condensate formation. This pipeline asks whether H² radial position, angular proximity, and community membership can predict IDR content, LLPS propensity, and functional annotation.

---

## Repository Structure

```
.
├── download_string.R            # 01 – Download STRING v12.0 human PPI data
├── build_string_network.R       # 05 – Filter, build igraph LCC, save GraphML
├── prepare_mercator_input.R     # 06 – Export edge list for Mercator
├── embed_hyperbolic.py          # 07 – Self-contained PSO hyperbolic embedding (H²)
├── convert_mercator_output.py   # 08 – Convert Mercator .inf_coord → CSV + Poincaré coords
├── extract_alphafold_idr.py     # 09 – Extract IDR % from AlphaFold mmCIF files
├── extract_fuzdrop_llps.py      # 10 – Extract LLPS scores from FuzDrop result directories
├── ensp_hgnc_mapping.R          # 11 – Map ENSP → HGNC / ENSG / Entrez via BioMart
├── hyperbolic_distances.R       # 12 – All-pairs H² distances and angular differences
├── topology_metrics.R           # 13 – Degree, betweenness, closeness, clustering, PageRank
├── louvain_communities.R        # 15 – Louvain community detection
├── walktrap_communities.R       # 16 – Walktrap community detection (step sweep)
├── ora_communities.R            # 17 – GO + Reactome ORA per community
├── kmer_enrichment.R            # 18 – k-mer extraction (k = 3–12) with PWMs
└── kmer_motif_enrichment.R      # 19 – Fisher ORA of top k-mers against ELM + PROSITE
```

Script numbers reflect the original pipeline ordering; gaps indicate steps that are upstream data-preparation tasks handled outside this repository (e.g. downloading AlphaFold CIF archives and running FuzDrop).

---

## Data Sources

| Dataset | Version | Source |
|---|---|---|
| STRING human PPI | v12.0 | https://stringdb-downloads.org |
| AlphaFold human proteome | v4 / v6 CIF | https://ftp.ebi.ac.uk/pub/databases/alphafold |
| FuzDrop LLPS scores | web tool | https://fuzdrop.chem.elte.hu |
| Mercator hyperbolic embedder | latest | https://mercator.unizar.es |
| Ensembl BioMart | current | https://www.ensembl.org |
| ELM linear motif database | — | http://elm.eu.org |
| PROSITE motif database | — | https://prosite.expasy.org |

---

## Pipeline Overview

### Data acquisition and network construction

| Script | Input | Output |
|---|---|---|
| `download_string.R` | — | `data/STRING/9606.protein.links.full.v12.0.txt` |
| `build_string_network.R` | STRING links file | `data/STRING/*.graphml` (LCC, score ≥ 900) |
| `prepare_mercator_input.R` | GraphML | `data/STRING/STRING_9606_LCC_for_mercator.txt` |

### Hyperbolic embedding

Two alternative routes are provided:

**Route A – Mercator** (recommended, higher quality)
```bash
python -c "import mercator; mercator.embed(
    'data/STRING/STRING_9606_LCC_for_mercator.txt',
    output_name='data/STRING/STRING_9606_H2',
    seed=1, validation_mode=True, screen_mode=True)"
python convert_mercator_output.py
```

**Route B – self-contained PSO** (no external dependency)
```bash
python embed_hyperbolic.py
```

Both routes produce `data/STRING/embedding.csv` with columns `node, r_hyp, theta, rho, x, y`.

### IDR and LLPS data extraction

```bash
python extract_alphafold_idr.py \
    --cif-dir  data/AlphaFold/cif \
    --map-csv  data/AlphaFold/string_to_uniprot.csv \
    --out-csv  data/AlphaFold/alphafold_idr.csv

python extract_fuzdrop_llps.py \
    --data-dir data/FuzDrop \
    --out-csv  data/FuzDrop/llps_scores.csv
```

AlphaFold IDR is defined as the fraction of residues with pLDDT < 50.

### Identifier mapping and graph metrics

```r
Rscript ensp_hgnc_mapping.R     # ENSP → HGNC / Entrez (BioMart)
Rscript hyperbolic_distances.R  # all-pairs H² distance and angular difference matrices
Rscript topology_metrics.R      # node-level metrics added to GraphML
```

### Community detection and functional enrichment

```r
Rscript louvain_communities.R   # Louvain; output: data/STRING/Louvain_communities.csv
Rscript walktrap_communities.R  # Walktrap step sweep (2–50); exports best + user steps
Rscript ora_communities.R       # GO BP/MF/CC + Reactome ORA per Walktrap community
```

### k-mer analysis

Operates on per-group residue-state sequence files produced by upstream analyses (AIUPred / AlphaFold × DN / DM / ON / OM state classes).

```r
Rscript kmer_enrichment.R        # top-100 k-mers (k = 3–12) + PWMs per group
Rscript kmer_motif_enrichment.R  # Fisher ORA of top k-mers vs. ELM + PROSITE motifs
```

---

## Requirements

### R (≥ 4.3)

All packages are auto-installed on first run if absent.

| Package | Source |
|---|---|
| `data.table` | CRAN |
| `igraph` | CRAN |
| `R.utils` | CRAN |
| `openxlsx` | CRAN |
| `ggplot2` | CRAN |
| `biomaRt` | Bioconductor |
| `clusterProfiler` | Bioconductor |
| `org.Hs.eg.db` | Bioconductor |
| `ReactomePA` | Bioconductor |

### Python (≥ 3.9)

```
pip install numpy pandas networkx
```

The AlphaFold and FuzDrop extraction scripts use the standard library only.

### External tools (optional)

- **Mercator** – for Route A embedding: `pip install mercator-embedding`

---

## Suggested Execution Order

```
download_string.R
build_string_network.R
prepare_mercator_input.R
→ embed_hyperbolic.py  OR  (mercator + convert_mercator_output.py)
extract_alphafold_idr.py
extract_fuzdrop_llps.py
ensp_hgnc_mapping.R
hyperbolic_distances.R
topology_metrics.R
louvain_communities.R
walktrap_communities.R
ora_communities.R
kmer_enrichment.R
kmer_motif_enrichment.R
```

All R scripts are run from the repository root with `Rscript <script>.R`. All Python scripts accept `--help` for CLI options.

---

## Output Data Layout

```
data/
└── STRING/
    ├── 9606.protein.links.full.v12.0.txt
    ├── *.graphml                          (LCC with topology attributes)
    ├── STRING_9606_LCC_for_mercator.txt
    ├── embedding.csv
    ├── distance_matrix.rds
    ├── angular_difference_matrix.rds
    ├── STRING_ENSP_to_HGNC.csv
    ├── STRING_9606_TopologyMetrics.{csv,xlsx}
    ├── Louvain_communities.csv
    ├── steps_summary.csv
    └── Walktrap_steps<N>_Communities.csv
data/
└── AlphaFold/
    └── alphafold_idr.csv
data/
└── FuzDrop/
    └── llps_scores.csv
ORA_CommunityWise/
    ├── community_gene_mapping.csv
    ├── ORA_community_*.csv
    └── ORA_community_*.xlsx
Results_*/
    └── .../kmer_analysis/
        ├── top100_kmers_*.csv
        ├── top100_kmers_*_pwms_*.csv
        └── motif_enrichment/
            ├── motif_enrichment_*.csv
            └── motif_enrichment_all_groups.xlsx
```

---

## Citation

If you use this pipeline, please cite the associated manuscript (in preparation) and the underlying tools:

- **STRING** – Szklarczyk et al. (2023) *Nucleic Acids Research*
- **Mercator** – Papadopoulos et al. (2023) *Nature Communications*
- **AlphaFold** – Jumper et al. (2021) *Nature*
- **FuzDrop** – Hardenberg et al. (2022) *PNAS*
- **clusterProfiler** – Wu et al. (2021) *The Innovation*
- **ReactomePA** – Yu & He (2016) *Molecular BioSystems*

---

## Author

Frank Hause · `frank.hause@uvic.cat`
