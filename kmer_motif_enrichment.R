# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("openxlsx",   quietly = TRUE)) install.packages("openxlsx")

# ── Settings ──────────────────────────────────────────────────
kmer_dir      <- file.path(
  "Results_ABYSS_H4", "Community_IDR_MBM_Correlation",
  "Sector_Residue_State_Strings", "group_sequence_textfiles", "kmer_analysis"
)
seq_dir       <- file.path(
  "Results_ABYSS_H4", "Community_IDR_MBM_Correlation",
  "Sector_Residue_State_Strings", "group_sequence_textfiles"
)
elm_path      <- "ELM_Regex.tsv"
prosite_path  <- "PROSITE_Regex.dat"
k_values      <- 3L:12L
out_dir       <- file.path(kmer_dir, "motif_enrichment")

for (fp in c(kmer_dir, seq_dir, elm_path, prosite_path)) {
  if (!file.exists(fp) && !dir.exists(fp))
    stop("Required input not found: ", fp)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ───────────────────────────────────────────────────
read_seq_file <- function(fp) {
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

prosite_token_to_regex <- function(tok) {
  if (tok == "") return("")
  lead_anchor <- ""
  tail_anchor <- ""
  if (startsWith(tok, "^")) {
    lead_anchor <- "^"; tok <- substring(tok, 2L)
  }
  if (endsWith(tok, "$")) {
    tail_anchor <- "$"; tok <- substring(tok, 1L, nchar(tok) - 1L)
  }
  rep_m   <- regexec("^(.+)\\((\\d+)(,(\\d+))?\\)$", tok)
  rep_hit <- regmatches(tok, rep_m)[[1]]
  if (length(rep_hit) > 0L) {
    core     <- rep_hit[2]
    n1       <- as.integer(rep_hit[3])
    n2       <- if (length(rep_hit) >= 5L && !is.na(rep_hit[5]) && rep_hit[5] != "")
      as.integer(rep_hit[5]) else NA_integer_
    core_rgx <- prosite_token_to_regex(core)
    core_rgx <- sub("^\\^", "", core_rgx)
    core_rgx <- sub("\\$$", "", core_rgx)
    quant    <- if (is.na(n2)) paste0("{", n1, "}") else paste0("{", n1, ",", n2, "}")
    return(paste0(lead_anchor, "(?:", core_rgx, ")", quant, tail_anchor))
  }
  core <- if (tok == "x") {
    "."
  } else if (grepl("^\\[[A-Z]+\\]$", tok)) {
    tok
  } else if (grepl("^\\{[A-Z]+\\}$", tok)) {
    paste0("[^", substring(tok, 2L, nchar(tok) - 1L), "]")
  } else if (grepl("^[A-Z]$", tok)) {
    tok
  } else {
    tok
  }
  paste0(lead_anchor, core, tail_anchor)
}

prosite_to_regex <- function(pa) {
  p    <- gsub("\\s+", "", pa)
  p    <- gsub("\\.$", "", p)
  p    <- gsub("<", "^", p, fixed = TRUE)
  p    <- gsub(">", "$", p, fixed = TRUE)
  toks <- strsplit(p, "-", fixed = TRUE)[[1]]
  toks <- toks[toks != ""]
  if (length(toks) == 0L) return(NA_character_)
  rgx  <- vapply(toks, prosite_token_to_regex, character(1L))
  paste0(rgx, collapse = "")
}

read_prosite_regex <- function(path) {
  lines   <- readLines(path, warn = FALSE, encoding = "UTF-8")
  entries <- list()
  cur_id  <- NA_character_
  cur_ac  <- NA_character_
  cur_de  <- NA_character_
  cur_pa  <- character(0)

  flush_entry <- function() {
    if (!is.na(cur_id) && length(cur_pa) > 0L) {
      pa_join <- gsub("\\s+", "", paste(cur_pa, collapse = ""))
      entries[[length(entries) + 1L]] <<- data.table::data.table(
        motif_db    = "PROSITE",
        motif_id    = cur_id,
        accession   = cur_ac,
        motif_name  = cur_de,
        raw_pattern = pa_join,
        regex       = prosite_to_regex(pa_join)
      )
    }
  }

  for (ln in lines) {
    if (startsWith(ln, "//")) {
      flush_entry()
      cur_id <- NA_character_; cur_ac <- NA_character_
      cur_de <- NA_character_; cur_pa <- character(0)
    } else if (startsWith(ln, "ID")) {
      cur_id <- sub(";.*$", "", trimws(sub("^ID\\s+", "", ln)))
    } else if (startsWith(ln, "AC")) {
      cur_ac <- sub(";.*$", "", trimws(sub("^AC\\s+", "", ln)))
    } else if (startsWith(ln, "DE")) {
      cur_de <- trimws(sub("^DE\\s+", "", ln))
    } else if (startsWith(ln, "PA")) {
      cur_pa <- c(cur_pa, gsub(";", "", trimws(sub("^PA\\s+", "", ln))))
    }
  }

  if (length(entries) == 0L) return(data.table::data.table())
  out <- data.table::rbindlist(entries, use.names = TRUE, fill = TRUE)
  out[, regex_ok := vapply(regex, function(r) {
    ok <- TRUE
    tryCatch(grepl(r, "AAAA", perl = TRUE), error = function(e) { ok <<- FALSE })
    ok
  }, logical(1L))]
  out[regex_ok == TRUE][, regex_ok := NULL]
}

read_elm_regex <- function(path) {
  dt <- data.table::fread(path, sep = "\t", quote = "\"", comment.char = "#")
  data.table::setnames(dt, gsub('"', "", names(dt), fixed = TRUE))
  out <- dt[, .(
    motif_db    = "ELM",
    motif_id    = as.character(ELMIdentifier),
    accession   = as.character(Accession),
    motif_name  = as.character(FunctionalSiteName),
    raw_pattern = as.character(Regex),
    regex       = as.character(Regex)
  )]
  out[, regex_ok := vapply(regex, function(r) {
    ok <- TRUE
    tryCatch(grepl(r, "AAAA", perl = TRUE), error = function(e) { ok <<- FALSE })
    ok
  }, logical(1L))]
  out[regex_ok == TRUE][, regex_ok := NULL]
}

# ── Load motif databases ──────────────────────────────────────
elm_dt     <- read_elm_regex(elm_path)
prosite_dt <- read_prosite_regex(prosite_path)
motif_dt   <- data.table::rbindlist(list(elm_dt, prosite_dt), use.names = TRUE, fill = TRUE)
motif_dt   <- motif_dt[!is.na(regex) & regex != ""]
motif_dt   <- unique(motif_dt[, .(motif_db, motif_id, accession, motif_name, raw_pattern, regex)])

if (nrow(motif_dt) == 0L)
  stop("No valid motifs loaded from ELM / PROSITE.")

# ── Top k-mer files (produced by 18_kmer_enrichment.R) ────────
top_files <- list.files(
  kmer_dir,
  pattern    = "^top100_kmers_k3_to_k12_(AIUPred|AlphaFold)_(DN|DM|ON|OM)\\.csv$",
  full.names = TRUE
)
if (length(top_files) == 0L)
  stop("No per-group top-kmer files found in: ", kmer_dir,
       "\nRun 18_kmer_enrichment.R first.")

parse_top_meta <- function(fp) {
  bn  <- basename(fp)
  m   <- regexec("^top100_kmers_k3_to_k12_(AIUPred|AlphaFold)_(DN|DM|ON|OM)\\.csv$", bn)
  hit <- regmatches(bn, m)[[1]]
  if (length(hit) != 3L)
    return(data.table::data.table(source = NA_character_, state = NA_character_))
  data.table::data.table(source = hit[2], state = hit[3])
}

# ── ORA loop ──────────────────────────────────────────────────
all_res_by_group <- list()

for (tf in top_files) {
  meta <- parse_top_meta(tf)
  src  <- meta$source[1]
  st   <- meta$state[1]

  top_dt <- data.table::fread(tf)
  if (!all(c("kmer", "k", "count") %in% names(top_dt))) {
    warning("Skipping file with unexpected schema: ", tf)
    next
  }

  fg_kmers <- unique(toupper(as.character(top_dt$kmer)))

  seq_fp <- file.path(seq_dir, paste0(src, "_", st, "_sequences.txt"))
  if (!file.exists(seq_fp)) {
    warning("Sequence file missing for ", src, "_", st, ": ", seq_fp)
    next
  }

  seqs     <- read_seq_file(seq_fp)
  all_km   <- unlist(lapply(k_values, function(k) extract_kmers(seqs, k)), use.names = FALSE)
  universe <- unique(toupper(all_km))
  bg_kmers <- setdiff(universe, fg_kmers)

  if (length(fg_kmers) == 0L || length(universe) == 0L) {
    warning("Skipping file with empty foreground/universe: ", tf)
    next
  }

  n_motif  <- nrow(motif_dt)
  res_list <- vector("list", n_motif)

  for (i in seq_len(n_motif)) {
    r      <- motif_dt$regex[i]
    fg_hit <- sum(grepl(r, fg_kmers, perl = TRUE))
    bg_hit <- if (length(bg_kmers) > 0L) sum(grepl(r, bg_kmers, perl = TRUE)) else 0L

    a <- as.integer(fg_hit)
    b <- as.integer(length(fg_kmers) - fg_hit)
    c <- as.integer(bg_hit)
    d <- as.integer(length(bg_kmers) - bg_hit)

    p_val  <- 1
    or_val <- NA_real_
    if ((a + c) > 0L) {
      ft <- tryCatch(
        stats::fisher.test(
          matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE),
          alternative = "greater"
        ),
        error = function(e) NULL
      )
      if (!is.null(ft)) {
        p_val  <- ft$p.value
        or_val <- unname(ft$estimate)
      }
    }

    res_list[[i]] <- data.table::data.table(
      source       = src,
      state        = st,
      top_file     = basename(tf),
      seq_file     = basename(seq_fp),
      motif_db     = motif_dt$motif_db[i],
      motif_id     = motif_dt$motif_id[i],
      accession    = motif_dt$accession[i],
      motif_name   = motif_dt$motif_name[i],
      regex        = motif_dt$regex[i],
      n_fg_unique  = length(fg_kmers),
      n_bg_unique  = length(bg_kmers),
      n_fg_match   = a,
      n_bg_match   = c,
      fg_prop      = a / length(fg_kmers),
      bg_prop      = if (length(bg_kmers) > 0L) c / length(bg_kmers) else NA_real_,
      odds_ratio   = or_val,
      p_value      = p_val
    )
  }

  res_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  res_dt[, p_adj_bh         := stats::p.adjust(p_value, method = "BH")]
  res_dt[, significant_bh_0_05 := p_adj_bh < 0.05]
  data.table::setorder(res_dt, p_adj_bh, p_value, -fg_prop, motif_db, motif_id)

  grp_key <- paste0(src, "_", st)
  all_res_by_group[[grp_key]] <- res_dt

  data.table::fwrite(
    res_dt,
    file.path(out_dir, paste0("motif_enrichment_", src, "_", st, ".csv"))
  )
}

# ── Save combined XLSX ────────────────────────────────────────
if (length(all_res_by_group) == 0L)
  stop("No enrichment results generated. Check inputs.")

all_res_dt <- data.table::rbindlist(
  all_res_by_group,
  use.names = TRUE, fill = TRUE, idcol = "group_key"
)
data.table::setorder(all_res_dt, source, state, p_adj_bh, p_value, -fg_prop, motif_db, motif_id)

xlsx_out <- file.path(out_dir, "motif_enrichment_all_groups.xlsx")
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "all_results")
openxlsx::writeData(wb, sheet = "all_results", x = all_res_dt)

for (grp in names(all_res_by_group)) {
  sheet_nm <- substr(paste0("group_", grp), 1L, 31L)
  openxlsx::addWorksheet(wb, sheet_nm)
  openxlsx::writeData(wb, sheet = sheet_nm, x = all_res_by_group[[grp]])
}

openxlsx::saveWorkbook(wb, file = xlsx_out, overwrite = TRUE)
