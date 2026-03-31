# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

# ── Settings ──────────────────────────────────────────────────
seq_dir <- file.path(
  "Results_ABYSS_H4", "Community_IDR_MBM_Correlation",
  "Sector_Residue_State_Strings", "group_sequence_textfiles"
)
out_dir  <- file.path(seq_dir, "kmer_analysis")
k_values <- 3L:12L
top_n    <- 100L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seq_files <- list.files(
  seq_dir,
  pattern   = "^(AIUPred|AlphaFold)_(DN|DM|ON|OM)_sequences\\.txt$",
  full.names = TRUE
)
if (length(seq_files) == 0L)
  stop("No sequence text files found in: ", seq_dir)

# ── Helpers ───────────────────────────────────────────────────
parse_file_meta <- function(fp) {
  bn  <- basename(fp)
  m   <- regexec("^(AIUPred|AlphaFold)_(DN|DM|ON|OM)_sequences\\.txt$", bn)
  hit <- regmatches(bn, m)[[1]]
  if (length(hit) != 3L)
    return(data.table::data.table(source = NA_character_, state = NA_character_))
  data.table::data.table(source = hit[2], state = hit[3])
}

read_seq_lines <- function(fp) {
  x <- readLines(fp, warn = FALSE, encoding = "UTF-8")
  x <- toupper(trimws(x))
  x[!is.na(x) & x != ""]
}

extract_kmers <- function(seqs, k) {
  if (length(seqs) == 0L) return(character(0))
  out <- vector("list", length(seqs))
  for (i in seq_along(seqs)) {
    s <- seqs[i]
    n <- nchar(s)
    if (is.na(n) || n < k) {
      out[[i]] <- character(0)
    } else {
      starts   <- seq_len(n - k + 1L)
      out[[i]] <- substring(s, starts, starts + k - 1L)
    }
  }
  unlist(out, use.names = FALSE)
}

build_pwm_for_kmer <- function(kmer, alphabet) {
  aa <- strsplit(kmer, "", fixed = TRUE)[[1]]
  k  <- length(aa)
  if (k == 0L) return(data.table::data.table())
  data.table::rbindlist(lapply(seq_len(k), function(pos) {
    data.table::data.table(
      position = pos,
      residue  = alphabet,
      prob     = as.numeric(alphabet == aa[pos])
    )
  }), use.names = TRUE)
}

# ── Main loop ─────────────────────────────────────────────────
all_top <- list()
all_pwm <- list()

for (fp in seq_files) {
  meta <- parse_file_meta(fp)
  seqs <- read_seq_lines(fp)

  aa_alphabet <- sort(unique(unlist(strsplit(paste(seqs, collapse = ""), "", fixed = TRUE),
                                   use.names = FALSE)))
  if (length(aa_alphabet) == 0L)
    aa_alphabet <- c("A","C","D","E","F","G","H","I","K","L",
                     "M","N","P","Q","R","S","T","V","W","Y","X")

  for (k in k_values) {
    km <- extract_kmers(seqs, k)
    if (length(km) == 0L) next

    km_tab <- data.table::as.data.table(table(km))
    data.table::setnames(km_tab, c("kmer", "count"))
    km_tab[, count := as.integer(count)]
    data.table::setorder(km_tab, -count, kmer)
    km_top <- km_tab[seq_len(min(top_n, .N))]
    km_top[, `:=`(
      source     = meta$source[1],
      state      = meta$state[1],
      input_file = basename(fp),
      k          = as.integer(k),
      rank       = seq_len(.N)
    )]
    data.table::setcolorder(
      km_top,
      c("source", "state", "input_file", "k", "rank", "kmer", "count")
    )
    all_top[[length(all_top) + 1L]] <- km_top

    pwm_top <- data.table::rbindlist(lapply(seq_len(nrow(km_top)), function(i) {
      pw <- build_pwm_for_kmer(km_top$kmer[i], alphabet = aa_alphabet)
      if (nrow(pw) == 0L) return(NULL)
      pw[, `:=`(
        source     = km_top$source[i],
        state      = km_top$state[i],
        input_file = km_top$input_file[i],
        k          = km_top$k[i],
        rank       = km_top$rank[i],
        kmer       = km_top$kmer[i],
        count      = km_top$count[i]
      )]
      data.table::setcolorder(
        pw,
        c("source", "state", "input_file", "k", "rank", "kmer", "count",
          "position", "residue", "prob")
      )
      pw
    }), use.names = TRUE, fill = TRUE)

    if (nrow(pwm_top) > 0L)
      all_pwm[[length(all_pwm) + 1L]] <- pwm_top
  }
}

# ── Combine ───────────────────────────────────────────────────
top_dt <- if (length(all_top) > 0L)
  data.table::rbindlist(all_top, use.names = TRUE, fill = TRUE) else
  data.table::data.table()
pwm_dt <- if (length(all_pwm) > 0L)
  data.table::rbindlist(all_pwm, use.names = TRUE, fill = TRUE) else
  data.table::data.table()

if (nrow(top_dt) == 0L)
  stop("No k-mers were generated. Check sequence files in: ", seq_dir)

# ── Save combined ─────────────────────────────────────────────
data.table::fwrite(top_dt, file.path(out_dir, "top100_kmers_k3_to_k12_all_files.csv"))
if (nrow(pwm_dt) > 0L)
  data.table::fwrite(pwm_dt, file.path(out_dir, "top100_kmers_k3_to_k12_pwms_long.csv"))

# ── Save per source/state ─────────────────────────────────────
for (src in unique(top_dt$source)) {
  for (st in unique(top_dt$state[top_dt$source == src])) {
    slug    <- paste0(src, "_", st)
    sub_top <- top_dt[source == src & state == st]
    if (nrow(sub_top) == 0L) next
    data.table::fwrite(
      sub_top,
      file.path(out_dir, paste0("top100_kmers_k3_to_k12_", slug, ".csv"))
    )
    if (nrow(pwm_dt) > 0L) {
      sub_pwm <- pwm_dt[source == src & state == st]
      if (nrow(sub_pwm) > 0L)
        data.table::fwrite(
          sub_pwm,
          file.path(out_dir, paste0("top100_kmers_k3_to_k12_pwms_", slug, ".csv"))
        )
    }
  }
}
