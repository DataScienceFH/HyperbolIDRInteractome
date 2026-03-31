# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",     quietly = TRUE)) install.packages("igraph")

# ── Settings ──────────────────────────────────────────────────
string_dir <- file.path("data", "STRING")
out_file   <- file.path(string_dir, "STRING_9606_LCC_for_mercator.txt")

# ── Locate GraphML ─────────────────────────────────────────────
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

graph_path <- find_graphml(string_dir)
g <- igraph::read_graph(graph_path, format = "graphml")

# ── Export edge list ───────────────────────────────────────────
dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)

el <- igraph::as_data_frame(g, what = "edges")
data.table::fwrite(
  el[, c("from", "to")],
  file      = out_file,
  sep       = " ",
  col.names = FALSE
)


