# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table",  quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("openxlsx",    quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("clusterProfiler", "org.Hs.eg.db", "biomaRt", "ReactomePA")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

# ── Settings ──────────────────────────────────────────────────
input_file      <- file.path("data", "STRING",
                              "STRING_9606_Walktrap_steps18_Communities.csv")
out_dir         <- "ORA_CommunityWise"
go_ontologies   <- c("BP", "MF", "CC")
p_cutoff        <- 0.05
q_cutoff        <- 0.20
min_gs_size     <- 10
max_gs_size     <- 500

# ── Helpers ───────────────────────────────────────────────────
use_ensembl_safe <- function(dataset = "hsapiens_gene_ensembl") {
  hosts <- c("www.ensembl.org", "useast.ensembl.org", "asia.ensembl.org")
  for (host in hosts) {
    mart <- tryCatch(
      biomaRt::useEnsembl(biomart = "genes", dataset = dataset, host = host),
      error = function(e) NULL
    )
    if (!is.null(mart)) {
      return(mart)
    }
  }
  stop("Could not connect to any Ensembl host.")
}

# ── Load community file ────────────────────────────────────────
if (!file.exists(input_file)) {
  stop("Community file not found: ", input_file,
       "\nRun 16_walktrap_communities.R first.")
}

comm_dt <- data.table::fread(input_file)
comm_dt[, node      := as.character(node)]
comm_dt[, community := as.character(community)]

# ── ENSP → Entrez mapping via BioMart ─────────────────────────
mart <- use_ensembl_safe()

raw_map <- biomaRt::getBM(
  attributes = c("ensembl_peptide_id", "ensembl_gene_id",
                 "external_gene_name", "entrezgene_id"),
  filters    = "ensembl_peptide_id",
  values     = unique(comm_dt$node),
  mart       = mart
)

map_dt <- data.table::as.data.table(raw_map)
data.table::setnames(
  map_dt,
  old = c("ensembl_peptide_id", "ensembl_gene_id", "external_gene_name", "entrezgene_id"),
  new = c("node", "ensembl_gene_id", "gene_symbol", "entrez_id")
)
map_dt <- map_dt[!is.na(entrez_id) & entrez_id != ""]
map_dt[, entrez_id := as.character(entrez_id)]
map_dt <- unique(map_dt, by = c("node", "entrez_id"))

comm_map_dt <- merge(comm_dt, map_dt, by = "node", all.x = TRUE, all.y = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(comm_map_dt, file.path(out_dir, "community_gene_mapping.csv"))

universe_entrez <- unique(comm_map_dt[!is.na(entrez_id), entrez_id])
if (length(universe_entrez) == 0L) stop("No Entrez IDs could be mapped.")

# ── ORA per community ─────────────────────────────────────────
community_ids <- sort(unique(comm_map_dt$community))
all_results   <- list()
summary_list  <- list()

for (cid in community_ids) {
  sub       <- comm_map_dt[community == cid]
  entrez    <- unique(sub[!is.na(entrez_id), entrez_id])
  symbols   <- unique(sub[!is.na(gene_symbol) & gene_symbol != "", gene_symbol])
  n_nodes   <- nrow(sub)
  n_entrez  <- length(entrez)
  n_symbols <- length(symbols)

  if (n_nodes < min_gs_size || n_nodes > max_gs_size || n_entrez < min_gs_size) {
    status <- if (n_nodes < min_gs_size || n_nodes > max_gs_size)
      "skipped_size_filter" else "skipped_too_few_mapped_genes"
    summary_list[[length(summary_list) + 1L]] <- data.table::data.table(
      community = cid, n_nodes = n_nodes, n_gene_symbols = n_symbols,
      n_entrez = n_entrez, n_sig_GO_BP = 0L, n_sig_GO_MF = 0L,
      n_sig_GO_CC = 0L, n_sig_REAC = 0L, n_sig_total = 0L, status = status
    )
    next
  }

  term_counts <- list(GO_BP = 0L, GO_MF = 0L, GO_CC = 0L, REAC = 0L)

  # GO ORA
  for (ont in go_ontologies) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene          = entrez,
        universe      = universe_entrez,
        OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cutoff,
        qvalueCutoff  = q_cutoff,
        minGSSize     = min_gs_size,
        maxGSSize     = max_gs_size,
        readable      = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0L) {
      go_dt <- data.table::as.data.table(as.data.frame(ego))
      go_dt[, `:=`(database = paste0("GO_", ont), community = cid,
                   n_nodes = n_nodes, n_gene_symbols = n_symbols, n_entrez = n_entrez)]
      all_results[[length(all_results) + 1L]] <- go_dt
      term_counts[[paste0("GO_", ont)]] <- nrow(go_dt)
    }
  }

  # Reactome ORA
  ereac <- tryCatch(
    ReactomePA::enrichPathway(
      gene          = entrez,
      universe      = universe_entrez,
      organism      = "human",
      pvalueCutoff  = p_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff  = q_cutoff,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      readable      = TRUE
    ),
    error = function(e) NULL
  )
  if (!is.null(ereac) && nrow(as.data.frame(ereac)) > 0L) {
    reac_dt <- data.table::as.data.table(as.data.frame(ereac))
    reac_dt[, `:=`(database = "REAC", community = cid,
                   n_nodes = n_nodes, n_gene_symbols = n_symbols, n_entrez = n_entrez)]
    all_results[[length(all_results) + 1L]] <- reac_dt
    term_counts[["REAC"]] <- nrow(reac_dt)
  }

  total_terms <- Reduce(`+`, term_counts)
  status      <- if (total_terms > 0L) "annotated" else "no_significant_terms"
  summary_list[[length(summary_list) + 1L]] <- data.table::data.table(
    community      = cid,
    n_nodes        = n_nodes,
    n_gene_symbols = n_symbols,
    n_entrez       = n_entrez,
    n_sig_GO_BP    = term_counts$GO_BP,
    n_sig_GO_MF    = term_counts$GO_MF,
    n_sig_GO_CC    = term_counts$GO_CC,
    n_sig_REAC     = term_counts$REAC,
    n_sig_total    = total_terms,
    status         = status
  )
}

# ── Save results ───────────────────────────────────────────────
summary_dt <- data.table::rbindlist(summary_list, fill = TRUE)
data.table::fwrite(
  summary_dt,
  file.path(out_dir, "ORA_community_summary_GO_BP_GO_MF_GO_CC_REAC.csv")
)

combined_dt <- data.table::data.table()
if (length(all_results) > 0L) {
  combined_dt  <- data.table::rbindlist(all_results, fill = TRUE)
  front_cols   <- c("database", "community", "n_nodes", "n_gene_symbols", "n_entrez")
  data.table::setcolorder(combined_dt, c(front_cols, setdiff(names(combined_dt), front_cols)))
  data.table::fwrite(combined_dt,
    file.path(out_dir, "ORA_community_all_combined_GO_BP_GO_MF_GO_CC_REAC.csv"))
}

# XLSX workbook
xlsx_out <- file.path(out_dir, "ORA_community_enrichment_GO_BP_GO_MF_GO_CC_REAC.xlsx")
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "community_gene_mapping"); openxlsx::writeData(wb, "community_gene_mapping", comm_map_dt)
openxlsx::addWorksheet(wb, "summary");               openxlsx::writeData(wb, "summary", summary_dt)
openxlsx::addWorksheet(wb, "all_combined");          openxlsx::writeData(wb, "all_combined", combined_dt)
openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
