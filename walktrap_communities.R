# ── Packages ──────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("igraph",     quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggplot2",    quietly = TRUE)) install.packages("ggplot2")

# ── Settings ──────────────────────────────────────────────────
string_dir           <- file.path("data", "STRING")
figures_dir          <- "figures"
walktrap_steps_range <- seq(2L, 50L, by = 1L)  # walk-length sweep
walktrap_user_steps  <- 18L                     # step count exported as community file

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

# ── Output paths ──────────────────────────────────────────────
dir.create(string_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv  <- file.path(string_dir, "steps_summary.csv")
user_csv     <- file.path(string_dir,
                          sprintf("Walktrap_steps%d_Communities.csv",
                                  walktrap_user_steps))
pdf_path     <- file.path(figures_dir, "Walktrap_modularity_vs_steps.pdf")

# ── Load graph ─────────────────────────────────────────────────
graph_path <- find_graphml(string_dir)
g <- igraph::read_graph(graph_path, format = "graphml")

# ── Walk-length sweep ─────────────────────────────────────────
steps_range <- sort(unique(as.integer(walktrap_steps_range)))

stats_df     <- data.frame(steps = steps_range, n_communities = NA_integer_, modularity = NA_real_)
memberships  <- vector("list", length(steps_range))

for (ii in seq_along(steps_range)) {
  s   <- steps_range[ii]
  wt  <- igraph::cluster_walktrap(g, steps = s)
  mem <- igraph::membership(wt)
  memberships[[ii]]           <- mem
  stats_df$n_communities[ii]  <- length(unique(mem))
  stats_df$modularity[ii]     <- igraph::modularity(g, mem)
}

utils::write.csv(stats_df, summary_csv, row.names = FALSE)

# ── Best steps by modularity ──────────────────────────────────
best_idx   <- which.max(stats_df$modularity)
best_steps <- stats_df$steps[best_idx]

# ── User-defined step export ──────────────────────────────────
wt_user       <- igraph::cluster_walktrap(g, steps = walktrap_user_steps)
mem_user      <- igraph::membership(wt_user)
mod_user      <- igraph::modularity(g, mem_user)

utils::write.csv(
  data.frame(node = names(mem_user), community = as.integer(mem_user)),
  user_csv,
  row.names = FALSE
)

# ── Modularity vs. steps plot ─────────────────────────────────
mod_range  <- range(stats_df$modularity,    na.rm = TRUE)
comm_range <- range(stats_df$n_communities, na.rm = TRUE)

scale_factor <- if (diff(comm_range) > 0)
  diff(mod_range) / diff(comm_range)
else 1

stats_df$n_communities_scaled <-
  mod_range[1] + (stats_df$n_communities - comm_range[1]) * scale_factor

grDevices::pdf(pdf_path, width = 7.5, height = 4.8)

p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = steps)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mod_range[1], ymax = modularity),
    fill = "#751D6C", alpha = 0.45
  ) +
  ggplot2::geom_line(ggplot2::aes(y = modularity),
                     color = "#751D6C", linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(y = modularity),
                      color = "black", size = 2.6) +
  ggplot2::geom_line(ggplot2::aes(y = n_communities_scaled),
                     color = "#163646", linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(y = n_communities_scaled),
                      color = "#163646", size = 2.4) +
  ggplot2::geom_vline(xintercept = best_steps, linetype = 2,
                      linewidth = 1.25, color = "#DAA520") +
  ggplot2::annotate("text", x = best_steps,
                    y = max(stats_df$modularity, na.rm = TRUE),
                    label = paste0("best steps = ", best_steps),
                    vjust = -0.6, hjust = 0.5, size = 3.5) +
  ggplot2::scale_y_continuous(
    name = "Modularity",
    sec.axis = ggplot2::sec_axis(
      trans = ~ (. - mod_range[1]) / scale_factor + comm_range[1],
      name  = "Number of communities"
    )
  ) +
  ggplot2::scale_x_continuous(breaks = steps_range) +
  ggplot2::labs(title = "Modularity vs random-walk steps",
                x     = "Walk length (steps)") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

print(p)
grDevices::dev.off()
