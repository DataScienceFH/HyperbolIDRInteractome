# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")

# ── Settings ──────────────────────────────────────────────────
species        <- 9606     # NCBI taxonomy ID (human)
string_version <- "12.0"
out_dir        <- file.path("data", "STRING")
overwrite      <- FALSE    # set TRUE to force re-download

# ── Download + decompress ─────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gz_name  <- paste0(species, ".protein.links.full.v", string_version, ".txt.gz")
txt_name <- sub("\\.gz$", "", gz_name)
gz_path  <- file.path(out_dir, gz_name)
txt_path <- file.path(out_dir, txt_name)

url <- paste0(
  "https://stringdb-downloads.org/download/protein.links.full.v",
  string_version, "/",
  species, ".protein.links.full.v", string_version, ".txt.gz"
)

if (!file.exists(gz_path) || isTRUE(overwrite)) {
  utils::download.file(url = url, destfile = gz_path, mode = "wb", quiet = FALSE)
}

if (!file.exists(txt_path) || isTRUE(overwrite)) {
  R.utils::gunzip(filename = gz_path, destname = txt_path, overwrite = isTRUE(overwrite))
}

