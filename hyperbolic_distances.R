# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

# ── Settings ──────────────────────────────────────────────────
embedding_path   <- file.path("data", "STRING", "embedding.csv")
out_dir          <- file.path("data", "STRING")

write_dense_csv  <- FALSE  # set TRUE to export CSV.GZ (very large)
write_angular_csv <- FALSE

# ── Load embedding ─────────────────────────────────────────────
if (!file.exists(embedding_path)) stop("Embedding file not found: ", embedding_path)

emb <- data.table::fread(embedding_path)
required <- c("node", "r_hyp", "theta")
missing  <- setdiff(required, colnames(emb))
if (length(missing) > 0L) stop("Missing columns: ", paste(missing, collapse = ", "))

emb <- emb[!is.na(node) & !is.na(r_hyp) & !is.na(theta)]
if (anyDuplicated(emb$node)) {
  emb <- emb[!duplicated(node)]
}

nodes <- emb$node
r_hyp <- as.numeric(emb$r_hyp)
theta <- as.numeric(emb$theta) %% (2 * pi)
n     <- length(nodes)

# ── Angular difference matrix ──────────────────────────────────
theta_diff   <- abs(outer(theta, theta, FUN = "-"))
angular_diff <- pi - abs(pi - theta_diff)
diag(angular_diff) <- 0
rownames(angular_diff) <- colnames(angular_diff) <- nodes

# ── Hyperbolic distance matrix ─────────────────────────────────
h2_dist      <- acosh(arg)
diag(h2_dist) <- 0
rownames(h2_dist) <- colnames(h2_dist) <- nodes

# ── Save RDS (primary format) ──────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_dist <- file.path(out_dir, "distance_matrix.rds")
rds_ang  <- file.path(out_dir, "angular_difference_matrix.rds")
saveRDS(h2_dist,      file = rds_dist)
saveRDS(angular_diff, file = rds_ang)

# ── Optional CSV.GZ export ────────────────────────────────────
if (isTRUE(write_dense_csv)) {
  csv_dist <- file.path(out_dir, "STRING_9606_H2_distance_matrix.csv.gz")
  data.table::fwrite(data.table::as.data.table(h2_dist, keep.rownames = "node"), file = csv_dist)
}

if (isTRUE(write_angular_csv)) {
  csv_ang <- file.path(out_dir, "STRING_9606_H2_angular_difference_matrix.csv.gz")
  data.table::fwrite(data.table::as.data.table(angular_diff, keep.rownames = "node"), file = csv_ang)
}

