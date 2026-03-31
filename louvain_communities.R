# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",     quietly = TRUE)) install.packages("igraph")

# ── Settings ──────────────────────────────────────────────────
string_dir <- file.path("data", "STRING")
out_file   <- file.path(string_dir, "Louvain_communities.csv")

# ── Helpers ───────────────────────────────────────────────────
find_graphml <- function(dir_path, prefer = "ENSP_unweighted_LCC") {
  candidates <- list.files(dir_path, pattern = "\\.graphml$", full.names = TRUE)
  if (!is.null(prefer)) {
    pref <- candidates[grepl(prefer, basename(candidates), fixed = TRUE)]
    if (length(pref) > 0L) candidates <- pref
  }
  if (length(candidates) > 1L)
    candidates <- candidates[which.max(file.info(candidates)$mtime)]
  normalizePath(candidates, winslash = "/", mustWork = TRUE)
}

# ── Load graph ─────────────────────────────────────────────────
graph_path <- find_graphml(string_dir)
g <- igraph::read_graph(graph_path, format = "graphml")

# ── Louvain community detection ────────────────────────────────
comm       <- igraph::cluster_louvain(g)
membership <- igraph::membership(comm)
modularity <- igraph::modularity(comm)

# ── Save ───────────────────────────────────────────────────────
out_dt <- data.table::data.table(
  node      = names(membership),
  community = as.integer(membership)
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out_dt, out_file)
