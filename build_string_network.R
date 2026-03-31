# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",     quietly = TRUE)) install.packages("igraph")

# ── Settings ──────────────────────────────────────────────────
string_dir <- file.path("data", "STRING")
score_min  <- 900L

# ── Helpers ───────────────────────────────────────────────────
find_string_links_file <- function(dir_path) {
  candidates <- list.files(dir_path, pattern = "protein\\.links", full.names = TRUE)
  candidates <- candidates[!grepl("\\.gz$", candidates, ignore.case = TRUE) &
                             grepl("\\.txt$", candidates, ignore.case = TRUE)]
  if (length(candidates) == 0L)
    stop("No decompressed STRING links file found in: ", dir_path)
  if (length(candidates) > 1L)
    candidates <- candidates[which.max(file.info(candidates)$mtime)]
  normalizePath(candidates, winslash = "/", mustWork = TRUE)
}

extract_string_metadata <- function(file_path) {
  fname   <- basename(file_path)
  species <- sub("\\..*$", "", fname)
  if (!grepl("^[0-9]+$", species))
    stop("Cannot extract numeric species ID from filename: ", fname)
  vm <- regexpr("v[0-9.]+", fname)
  version <- if (vm[1] > 0L)
    substr(fname, vm[1] + 1L, vm[1] + attr(vm, "match.length") - 1L)
  else NA_character_
  list(species = species, version = version)
}

# ── Locate input ───────────────────────────────────────────────
string_file <- find_string_links_file(string_dir)
meta        <- extract_string_metadata(string_file)

out_prefix <- paste(
  "STRING",
  paste0("tax", meta$species),
  if (!is.na(meta$version)) paste0("v", meta$version) else NULL,
  paste0("score", score_min),
  "ENSP_unweighted_LCC",
  sep = "_"
)
out_graphml <- file.path(string_dir, paste0(out_prefix, ".graphml"))

# ── Load + filter STRING ───────────────────────────────────────
dt <- data.table::fread(string_file)

required <- c("protein1", "protein2", "combined_score")
missing  <- setdiff(required, colnames(dt))
if (length(missing) > 0L) stop("Missing columns: ", paste(missing, collapse = ", "))

dt <- dt[combined_score >= score_min]

tax_prefix <- paste0("^", meta$species, "\\.")
dt[, protein1 := sub(tax_prefix, "", protein1)]
dt[, protein2 := sub(tax_prefix, "", protein2)]

edges <- unique(dt[, .(protein1, protein2)])
edges <- edges[protein1 != protein2]

# ── Build igraph ───────────────────────────────────────────────
g <- igraph::graph_from_data_frame(
  d        = data.frame(from = edges$protein1, to = edges$protein2),
  directed = FALSE
)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# ── Largest connected component ────────────────────────────────
comp   <- igraph::components(g)
lcc_id <- which.max(comp$csize)
g_lcc  <- igraph::induced_subgraph(g, vids = which(comp$membership == lcc_id))

# ── Save ───────────────────────────────────────────────────────
igraph::write_graph(g_lcc, file = out_graphml, format = "graphml")
