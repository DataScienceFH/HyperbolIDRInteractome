# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",     quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("openxlsx",   quietly = TRUE)) install.packages("openxlsx")

# ── Settings ──────────────────────────────────────────────────
string_dir <- file.path("data", "STRING")

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

derive_prefix <- function(graph_path, suffix) {
  base <- tools::file_path_sans_ext(basename(graph_path))
  m    <- regexpr("tax([0-9]+)", base, perl = TRUE)
  if (m[1] > 0L) {
    species <- substr(base, m[1] + 3L, m[1] + attr(m, "match.length") - 1L)
    return(paste0("STRING_", species, "_", suffix))
  }
  paste0(base, "_", suffix)
}

# ── Load graph ─────────────────────────────────────────────────
graph_path <- find_graphml(string_dir)
g          <- igraph::read_graph(graph_path, format = "graphml")
prefix     <- derive_prefix(graph_path, "TopologyMetrics")
label      <- tools::file_path_sans_ext(basename(graph_path))

message("Graph: ", igraph::vcount(g), " nodes | ", igraph::ecount(g), " edges")

# ── Compute metrics on LCC ────────────────────────────────────
comp  <- igraph::components(g)
giant <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))

deg <- igraph::degree(giant, mode = "all", loops = FALSE, normalized = FALSE)

bet <- igraph::betweenness(giant, directed = FALSE, normalized = TRUE)

clo <- igraph::closeness(giant, mode = "all", normalized = TRUE)

clu <- igraph::transitivity(giant, type = "localundirected", isolates = "zero")

pgr <- igraph::page_rank(giant, directed = FALSE)$vector

metrics_dt <- data.table::data.table(
  node                   = igraph::V(giant)$name,
  degree                 = as.numeric(deg),
  betweenness            = as.numeric(bet),
  closeness              = as.numeric(clo),
  clustering_coefficient = as.numeric(clu),
  pagerank               = as.numeric(pgr)
)

# ── Attach to graph + overwrite GraphML ───────────────────────
vnames <- igraph::V(g)$name
idx    <- match(vnames, metrics_dt$node)

igraph::V(g)$Degree                 <- ifelse(is.na(idx), NA_real_, metrics_dt$degree[idx])
igraph::V(g)$Betweenness            <- ifelse(is.na(idx), NA_real_, metrics_dt$betweenness[idx])
igraph::V(g)$Closeness              <- ifelse(is.na(idx), NA_real_, metrics_dt$closeness[idx])
igraph::V(g)$ClusteringCoefficient  <- ifelse(is.na(idx), NA_real_, metrics_dt$clustering_coefficient[idx])
igraph::V(g)$PageRank               <- ifelse(is.na(idx), NA_real_, metrics_dt$pagerank[idx])

igraph::write_graph(g, graph_path, format = "graphml")

# ── Save table + workbook ─────────────────────────────────────
csv_path  <- file.path(string_dir, paste0(prefix, ".csv"))
xlsx_path <- file.path(string_dir, paste0(prefix, ".xlsx"))

data.table::fwrite(metrics_dt, csv_path)
wb    <- openxlsx::createWorkbook()
sheet <- substr(prefix, 1L, 31L)
openxlsx::addWorksheet(wb, sheet)
openxlsx::writeData(wb, sheet, metrics_dt)
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

