# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table",  quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",      quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt",     quietly = TRUE)) BiocManager::install("biomaRt", ask = FALSE, update = FALSE)

# ── Settings ──────────────────────────────────────────────────
string_dir <- file.path("data", "STRING")
out_path   <- file.path(string_dir, "STRING_ENSP_to_HGNC.csv")

# ── Helpers ───────────────────────────────────────────────────
find_graphml <- function(dir_path, prefer = "ENSP_unweighted_LCC") {
  candidates <- list.files(dir_path, pattern = "\\.graphml$", full.names = TRUE)
  if (length(candidates) == 0L) stop("No GraphML files found in: ", dir_path)
  if (!is.null(prefer)) {
    pref <- candidates[grepl(prefer, basename(candidates), fixed = TRUE)]
    if (length(pref) > 0L) candidates <- pref
  }
  if (length(candidates) > 1L)
    candidates <- candidates[which.max(file.info(candidates)$mtime)]
  normalizePath(candidates, winslash = "/", mustWork = TRUE)
}

use_ensembl_safe <- function(dataset = "hsapiens_gene_ensembl") {
  hosts <- c(
    "https://www.ensembl.org",
    "https://useast.ensembl.org",
    "https://asia.ensembl.org"
  )
  for (h in hosts) {
    mart <- try(biomaRt::useMart("ensembl", dataset = dataset, host = h), silent = TRUE)
    if (!inherits(mart, "try-error")) {
      return(mart)
    }
  }
  stop("Could not connect to any Ensembl host via biomaRt.")
}

# ── Load graph + extract ENSP IDs ─────────────────────────────
graph_path <- find_graphml(string_dir)
g    <- igraph::read_graph(graph_path, format = "graphml")
ensp <- unique(igraph::V(g)$name)
ensp <- ensp[!is.na(ensp) & ensp != ""]

# ── BioMart query ─────────────────────────────────────────────
ensembl <- use_ensembl_safe()

raw <- biomaRt::getBM(
  attributes = c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters    = "ensembl_peptide_id",
  values     = ensp,
  mart       = ensembl
)

if (nrow(raw) == 0L) stop("biomaRt returned zero rows.")

dt <- data.table::as.data.table(raw)
data.table::setnames(
  dt,
  old = c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  new = c("ENSP", "ENSG", "HGNC", "EntrezID")
)

dt[, HGNC     := ifelse(HGNC == "",     NA_character_, HGNC)]
dt[, EntrezID := as.character(EntrezID)]
dt[EntrezID %in% c("", "NA"), EntrezID := NA_character_]
dt <- dt[!is.na(ENSP) & ENSP != ""]
dt <- unique(dt, by = c("ENSP", "HGNC", "ENSG", "EntrezID"))

# ── Save ───────────────────────────────────────────────────────
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(dt, out_path)
