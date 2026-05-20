# Directorio de trabajo dinámico (establece el WD en la ubicación del script)
if (!interactive()) {
  script_path <- sub(
    "--file=", "", grep("--file=", commandArgs(), value = TRUE)
  )
  if (length(script_path) > 0) setwd(dirname(normalizePath(script_path)))
} else if (requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
}
# ==============================================================================
# SCRIPT 8 - EXPRESIÓN DIFERENCIAL Y VALIDACIÓN DE ML
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggrepel)
  library(ggplot2)
  library(tidyverse)
  library(rlang)
})

# --- CONFIGURACIÓN DE RUTAS ---
path_genes <- list(
  braaksc = "data/genes-braaksc.csv",
  ceradsc = "data/genes-ceradsc.csv",
  cogdx   = "data/genes-cogdx.csv"
)

path_core_genes <- "biological_analysis/gene_lists/01_genes_core_masters.txt"

# Definición de contrastes biológicos
comparisons <- list(
  braaksc = list(
    control = c(0, 1, 2), case = c(5, 6), name = "Braak Low vs High"
  ),
  ceradsc = list(
    control = c(4), case = c(1), name = "No Plaques vs Frequent"
  ),
  cogdx = list(
    control = c(1), case = c(4, 5), name = "Healthy vs Alzheimer"
  )
)

cols_meta <- c(
  "projid", "individualID", "SampleID", "braaksc", "ceradsc", "cogdx"
)

# Asegurar estructura de carpetas
dir.create(
  "differential_expression/ml_validation",
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  "differential_expression/full_range",
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  "differential_expression/significant_ranked_genes",
  showWarnings = FALSE, recursive = TRUE
)

# ==============================================================================
# FUNCIÓN PARA CARGAR GENES SELECCIONADOS POR ML
# ==============================================================================
get_ml_genes <- function(ds_key) {
  if (!dir.exists("best_k_selection")) {
    return(character(0))
  }
  matches <- list.dirs("best_k_selection", recursive = FALSE, full.names = TRUE)
  winner_dir <- matches[grepl(ds_key, matches) & grepl("winner$", matches)][1]

  if (is.na(winner_dir) || length(winner_dir) == 0) {
    return(character(0))
  }

  path_genes <- file.path(winner_dir, "selected_genes.txt")
  if (!file.exists(path_genes)) {
    return(character(0))
  }

  lines <- readLines(path_genes) |>
    stringr::str_trim()
  lines[lines != ""]
}

# ==============================================================================
# FUNCIÓN PRINCIPAL DE VALIDACIÓN
# ==============================================================================
validate_ml_genes <- function(ds_key) {
  cat(paste("\n  Procesando validación para:", toupper(ds_key), "\n"))
  cfg <- comparisons[[ds_key]]
  ml_genes <- get_ml_genes(ds_key)

  if (length(ml_genes) == 0) {
    cat(
      "    Advertencia: No se hallaron genes de ML",
      "para este conjunto de datos.\n"
    )
    return(NULL)
  }

  # 1. Cargar conjunto de datos (dataset)
  df <- readr::read_csv(path_genes[[ds_key]], show_col_types = FALSE)

  # 2. Filtrar muestras
  df_filt <- df |>
    dplyr::filter(!!rlang::sym(ds_key) %in% c(cfg$control, cfg$case)) |>
    dplyr::mutate(
      Group = factor(
        ifelse(!!rlang::sym(ds_key) %in% cfg$control, "Control", "Case"),
        levels = c("Control", "Case")
      )
    )

  # 3. Preparar matriz
  gene_cols <- setdiff(colnames(df_filt), c(cols_meta, "Group"))
  expr_mat <- t(as.matrix(df_filt[, gene_cols]))

  # 4. Limma
  design <- model.matrix(~Group, data = df_filt)
  colnames(design) <- c("Intercept", "Case_vs_Control")

  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  # 5. Resultados
  res <- topTable(fit, coef = "Case_vs_Control", number = Inf) |>
    tibble::rownames_to_column("Symbol") |>
    dplyr::mutate(
      ML_Selected   = .data$Symbol %in% ml_genes,
      Significant   = .data$adj.P.Val < 0.05 & abs(.data$logFC) > 0.5,
      Validated     = .data$ML_Selected & .data$Significant
    )

  # --- Métricas ---
  n_ml <- length(ml_genes)
  n_validated <- sum(res$Validated)
  n_degs_total <- sum(res$Significant)
  sensitivity <- round((n_validated / max(n_degs_total, 1)) * 100, 1)

  cat(sprintf(
    "    Validados (ML + DEA): %d/%d (%.1f%%)\n",
    n_validated, n_ml, (n_validated / n_ml) * 100
  ))
  cat(sprintf(
    "    Cobertura de DEGs: %.1f%% del total de DEGs capturados.\n", sensitivity
  ))

  # --- GUARDAR ---
  path_complete <- file.path(
    "differential_expression", paste0("DEA_", ds_key, "_complete.csv")
  )
  path_validated <- file.path(
    "differential_expression", "ml_validation",
    paste0("validated_", ds_key, ".csv")
  )

  readr::write_csv(res, path_complete)
  res |>
    dplyr::filter(.data$Validated == TRUE) |>
    readr::write_csv(path_validated)

  # --- GENES SIGNIFICATIVOS ORDENADOS (CSV y TXT) ---
  res_sig_ranked <- res |>
    dplyr::filter(.data$Significant == TRUE) |>
    dplyr::arrange(dplyr::desc(abs(.data$t)))

  path_sig_csv <- file.path(
    "differential_expression", "significant_ranked_genes",
    paste0("significant_ranked_", ds_key, ".csv")
  )
  path_sig_txt <- file.path(
    "differential_expression", "significant_ranked_genes",
    paste0("symbols_significant_", ds_key, ".txt")
  )

  readr::write_csv(res_sig_ranked, path_sig_csv)
  writeLines(res_sig_ranked$Symbol, path_sig_txt)

  # --------------------------------------------------------

  # --- GRÁFICOS DE VOLCÁN (VOLCANO PLOTS) ---
  core_genes <- if (file.exists(path_core_genes)) {
    tmp_genes <- readLines(path_core_genes) |> stringr::str_trim()
    tmp_genes[tmp_genes != ""]
  } else {
    character(0)
  }

  plot_data <- res |>
    dplyr::mutate(
      Category = factor(dplyr::case_when(
        .data$Validated ~ "Validated (ML + DEA)",
        .data$ML_Selected & !.data$Significant ~ "ML Only",
        .data$Significant & !.data$ML_Selected ~ "DEA Only",
        TRUE ~ "Not Significant"
      ), levels = c(
        "Not Significant", "DEA Only", "ML Only", "Validated (ML + DEA)"
      ))
    ) |>
    dplyr::arrange(.data$Category) |>
    dplyr::filter(!is.na(.data$adj.P.Val), .data$adj.P.Val > 0)

  # Colores y tamaños mejorados
  colors <- c(
    "Validated (ML + DEA)" = "#E64A19",
    "ML Only" = "#FFB300",
    "DEA Only" = "#1E88E5",
    "Not Significant" = "#CFD8DC"
  )
  sizes <- c(
    "Validated (ML + DEA)" = 2.8, "ML Only" = 1.4,
    "DEA Only" = 1, "Not Significant" = 0.6
  )
  alphas <- c(
    "Validated (ML + DEA)" = 0.9, "ML Only" = 0.6,
    "DEA Only" = 0.4, "Not Significant" = 0.15
  )

  # Base del gráfico común
  base_p <- ggplot(
    plot_data,
    aes(
      x = .data$logFC, y = -log10(.data$adj.P.Val),
      color = .data$Category, size = .data$Category, alpha = .data$Category
    )
  ) +
    geom_point() +
    geom_hline(
      yintercept = -log10(0.05), linetype = "dashed",
      color = "gray40", alpha = 0.6
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5), linetype = "dashed",
      color = "gray40", alpha = 0.6
    ) +
    scale_color_manual(values = colors) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) +
    theme_minimal() +
    theme(
      text = element_text(color = "#263238"),
      plot.title = element_text(
        face = "bold", size = 16, margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = 11, color = "gray30", margin = margin(b = 20)
      ),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#F0F0F0")
    ) +
    labs(
      title = paste("Biological Validation:", cfg$name),
      subtitle = sprintf(
        "Cobertura de DEGs: %.1f%% | Validated Genes: %d",
        sensitivity, n_validated
      ),
      x = "Log2 Fold Change", y = "-Log10 FDR"
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 4)),
      size = "none",
      alpha = "none"
    )

  # 1. Gráfico completo sin zoom
  top_15_validated <- head(
    dplyr::filter(plot_data, .data$Validated == TRUE), 15
  )

  p_full <- base_p +
    geom_text_repel(
      data = top_15_validated,
      aes(label = .data$Symbol), size = 3, fontface = "bold"
    )

  path_full <- file.path(
    "differential_expression", "full_range",
    paste0("volcano_", ds_key, ".png")
  )
  ggsave(path_full, p_full, width = 9, height = 7, dpi = 300)

  # 2. Gráfico mejorado con zoom y etiquetas de genes centrales (Core Genes)
  limit_x <- max(
    abs(quantile(plot_data$logFC, c(0.01, 0.99), na.rm = TRUE))
  )
  limit_x <- min(max(limit_x, 1.5), 5)

  # Unimos top 15 y core_genes validados en un único set sin duplicados
  zoom_labeled_genes <- plot_data |>
    dplyr::filter(
      .data$Validated == TRUE,
      (.data$Symbol %in% top_15_validated$Symbol) |
        (.data$Symbol %in% core_genes)
    ) |>
    dplyr::distinct(.data$Symbol, .keep_all = TRUE)

  p_zoom <- base_p +
    coord_cartesian(xlim = c(-limit_x, limit_x)) +
    geom_text_repel(
      data = zoom_labeled_genes,
      aes(label = .data$Symbol),
      size = 3.5, fontface = "bold",
      box.padding = 0.6, point.padding = 0.4,
      force = 15, segment.color = "gray30", segment.alpha = 0.5,
      max.overlaps = 50, min.segment.length = 0
    )

  path_zoom <- file.path(
    "differential_expression",
    paste0("volcano_", ds_key, ".png")
  )
  ggsave(
    path_zoom,
    p_zoom,
    width = 9, height = 7, dpi = 300
  )

  cat(sprintf("    Archivos guardados en 'differential_expression/'\n"))
}

# --- EJECUCIÓN ---
for (ds in names(comparisons)) {
  validate_ml_genes(ds)
}

cat("\n  Proceso completado. Resultados en 'differential_expression/'.\n")
