# Directorio de trabajo dinámico (establece el WD en la ubicación del script)
if (!interactive()) {
  args <- commandArgs()
  script_path <- sub("--file=", "", grep("--file=", args, value = TRUE))
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(script_path)))
  }
} else if (requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
}

# ==============================================================================
# SCRIPT 11 - VALIDACIÓN MOLECULAR DE CONSENSO: TRIPLE CORE
#             (SUBNÚCLEO COMPLEMENTARIO)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(stringr)
})

# --- CONFIGURACIÓN DE FIRMAS DE ML ---

# 1. TRIPLE CORE PRINCIPAL (46 genes de consenso ML comunes de los 3 targets)
triple_core <- character(0)
tc_path <- "biological_analysis/gene_lists/01_genes_core_masters.txt"
if (file.exists(tc_path)) {
  triple_core <- readLines(tc_path) |> str_trim()
  triple_core <- triple_core[triple_core != ""]
}

# 2. SUBNÚCLEO COMPLEMENTARIO (16 genes del meta-análisis extra 3^2)
#    Adicional o complementario
subnucleo <- character(0)
sub_path <- "biological_analysis/gene_lists/08_genes_subnucleus.txt"
if (file.exists(sub_path)) {
  subnucleo <- readLines(sub_path) |> str_trim()
  subnucleo <- subnucleo[subnucleo != ""]
}

datasets <- c("braaksc", "ceradsc", "cogdx")

comparisons_list <- list(
  braaksc = list(
    control = c(0, 1, 2),
    case = c(5, 6),
    name = "Braak"
  ),
  ceradsc = list(
    control = c(4),
    case = c(1),
    name = "CERAD"
  ),
  cogdx = list(
    control = c(1),
    case = c(4, 5),
    name = "COGDX"
  )
)

# Definir directorios de salida
out_dir <- "boxplot_triple_core_validation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sub_out_dir <- file.path(out_dir, "Complementary_Subnucleus")
dir.create(sub_out_dir, showWarnings = FALSE, recursive = TRUE)

# --- FUNCIÓN AUXILIAR DE PLOTTING PREMIUM ---
generate_boxplots_grid <- function(
  df_plot, gene_list, title_text, save_path, columns_count = 4
) {
  # Filtrar solo genes que existen en las columnas del dataframe
  valid_genes <- intersect(gene_list, colnames(df_plot))

  if (length(valid_genes) == 0) {
    cat(
      "    [ADVERTENCIA] Ninguno de los genes especificados",
      "existe en las columnas del dataset.\n"
    )
    return(NULL)
  }

  # Seleccionar solo las columnas necesarias sin variables globales
  df_sub <- df_plot |>
    dplyr::select("Group", tidyselect::all_of(valid_genes))

  # Pivotar a largo
  df_long <- df_sub |>
    tidyr::pivot_longer(
      cols = -"Group",
      names_to = "Gene",
      values_to = "Expression"
    )

  # Dimensionar dinámicamente según número de genes
  n_genes <- length(valid_genes)
  n_rows <- ceiling(n_genes / columns_count)
  h <- 2.2 * n_rows + 1.6
  w <- 2.8 * columns_count

  # Crear ggplot con estética premium de publicación
  p <- ggplot(
    df_long,
    aes(x = .data$Group, y = .data$Expression, fill = .data$Group)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.75,
      color = "#2C3E50",
      width = 0.6,
      size = 0.7
    ) +
    geom_jitter(
      width = 0.2,
      alpha = 0.25,
      size = 0.6,
      aes(color = .data$Group)
    ) +
    facet_wrap(~ Gene, scales = "free_y", ncol = columns_count) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      label.x = 1.4
    ) +
    scale_fill_manual(
      values = c("Control" = "#00A896", "Case" = "#E63946")
    ) +
    scale_color_manual(
      values = c("Control" = "#028090", "Case" = "#D90429")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold", size = 15, color = "#2C3E50", hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10, color = "#7F8C8D", hjust = 0.5, margin = margin(b = 12)
      ),
      strip.background = element_rect(fill = "#EDF2F4", color = NA),
      strip.text = element_text(
        face = "bold", size = 10, color = "#2C3E50"
      ),
      panel.spacing = unit(1.0, "lines"),
      panel.border = element_rect(
        color = "#E5E7E9", fill = NA, size = 0.8
      ),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold", size = 10, color = "#2C3E50"
      ),
      legend.text = element_text(size = 9, color = "#2C3E50"),
      axis.title.y = element_text(
        face = "bold", size = 11, margin = margin(r = 8), color = "#2C3E50"
      ),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 9, color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#F0F3F4")
    ) +
    labs(
      title = title_text,
      subtitle = paste(
        "Comparison between Alzheimer's extreme phenotypes.",
        "Wilcoxon Rank Sum Test significance."
      ),
      y = "Normalized Expression Level",
      fill = "Clinical Group",
      color = "Clinical Group"
    )

  ggsave(save_path, plot = p, width = w, height = h, dpi = 300)
}

# --- EJECUCIÓN DEL PIPELINE ---
for (ds in datasets) {
  df_path <- file.path("data", paste0("genes-", ds, ".csv"))
  if (!file.exists(df_path)) {
    cat(sprintf("  [ADVERTENCIA] No existe el archivo: %s\n", df_path))
    next
  }

  cat(sprintf(
    "\n>>> VALIDANDO MOLECULARMENTE EL DATASET: %s\n", toupper(ds)
  ))
  df <- read_csv(df_path, show_col_types = FALSE)
  cfg <- comparisons_list[[ds]]

  # Filtrar grupos de contraste y etiquetar
  df_plot <- df |>
    filter(.data[[ds]] %in% c(cfg$control, cfg$case)) |>
    mutate(
      Group = factor(
        ifelse(.data[[ds]] %in% cfg$control, "Control", "Case"),
        levels = c("Control", "Case")
      )
    )

  # 1. Generar Figuras Principales: Triple Core (46 genes)
  #    - Segmentado por partes en la raíz de salida
  if (length(triple_core) > 0) {
    cat("    - Graficando Triple Core Principal (Consenso ML)... ")
    n_total <- length(triple_core)
    chunk_size <- 16
    n_chunks <- ceiling(n_total / chunk_size)

    for (i in seq_len(n_chunks)) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_total)
      chunk_genes <- triple_core[start_idx:end_idx]

      save_path_tc <- file.path(
        out_dir,
        paste0("boxplots_triple_core_part", i, "_", ds, ".png")
      )

      title_tc <- sprintf(
        "Expression Profiles: Triple Core Consenso (%s - Part %d of %d)",
        cfg$name, i, n_chunks
      )

      generate_boxplots_grid(
        df_plot, chunk_genes, title_tc, save_path_tc, columns_count = 4
      )
    }
    cat(sprintf("[ÉXITO] (Generados %d bloques principales)\n", n_chunks))
  }

  # 2. Generar Figura Extra/Complementaria: Subnúcleo Irreducible (16 genes)
  #    en subcarpeta dedicada
  cat("    - Graficando Subnúcleo Irreducible Complementario (Extra)...")
  save_path_sub <- file.path(
    sub_out_dir, paste0("boxplots_subnucleus_complementary_", ds, ".png")
  )
  title_sub <- paste(
    "Complementary Profiles: Subnúcleo Irreducible (Extra -",
    cfg$name, ")"
  )
  generate_boxplots_grid(
    df_plot, subnucleo, title_sub, save_path_sub, columns_count = 4
  )
  cat(" [ÉXITO]\n")
}

cat("\nProceso de validación de consenso molecular completado con éxito.\n")
