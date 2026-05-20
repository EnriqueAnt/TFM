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
# SCRIPT 12 - ANÁLISIS DE CORRELACIÓN INTER-ÓMICA MULTIÓMICA (GEN-METABOLITO)
# ==============================================================================
# Este script realiza la correlación estadística entre las firmas génicas
# seleccionadas por machine learning (consenso de Triple Core y complemento
# de Subnucleus) y los perfiles metabolómicos de ROSMAP. Identifica
# interacciones multiómicas significativas, calcula el rho de Spearman,
# p-valores nominales y FDR de Benjamini-Hochberg, y genera visualizaciones
# de grado de publicación (heatmaps y dotplots).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gtable)
})

# --- CONFIGURACIÓN DE FIRMAS DE ML ---

# 1. TRIPLE CORE PRINCIPAL (Consenso principal de 46 genes ML comunes)
triple_core <- character(0)
tc_path <- "biological_analysis/gene_lists/01_genes_core_masters.txt"
if (file.exists(tc_path)) {
  triple_core <- readLines(tc_path) |> str_trim()
  triple_core <- triple_core[triple_core != ""]
}

# 2. SUBNÚCLEO COMPLEMENTARIO (16 genes del meta-análisis extra 3^2)
subnucleo <- character(0)
sub_path <- "biological_analysis/gene_lists/08_genes_subnucleus.txt"
if (file.exists(sub_path)) {
  subnucleo <- readLines(sub_path) |> str_trim()
  subnucleo <- subnucleo[subnucleo != ""]
}

# Resolución robusta de rutas para los archivos integrados de Arrangements
genes_file <- "../Arrangements/genes_clinical_integrated.csv"
metab_file <- "../Arrangements/metabolites_clinical_integrated.csv"

if (!file.exists(genes_file)) {
  genes_file <- "Arrangements/genes_clinical_integrated.csv"
  metab_file <- "Arrangements/metabolites_clinical_integrated.csv"
}

if (!file.exists(genes_file)) {
  genes_file <- "genes_clinical_integrated.csv"
  metab_file <- "metabolites_clinical_integrated.csv"
}

# La carpeta de salida es hermana de results_enrichment, no hija
out_dir <- "results_multiomics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- FUNCIÓN CENTRAL PARA CORRELACIÓN MULTIÓMICA ---
run_correlation_analysis <- function(genes_df,
                                     metab_df,
                                     gene_list,
                                     prefix_name,
                                     base_out_dir,
                                     is_complementary = FALSE) {

  # Filtrar la lista de genes para conservar solo los presentes en el dataset
  valid_genes <- intersect(gene_list, colnames(genes_df))

  if (length(valid_genes) == 0) {
    cat(sprintf(
      "  [ADVERTENCIA] No genes from '%s' exist in the dataset.\n",
      prefix_name
    ))
    return(NULL)
  }

  cat(sprintf(
    "\n>>> RUNNING MULTI-OMICS CORRELATION FOR: %s (%d genes)\n",
    prefix_name, length(valid_genes)
  ))
  cat(sprintf(
    ">>> EJECUTANDO CORRELACIÓN MULTIÓMICA PARA: %s (%d genes)\n",
    prefix_name, length(valid_genes)
  ))

  # Alinear muestras comunes por projid
  common_samples <- intersect(genes_df$projid, metab_df$projid)
  cat(sprintf(
    "    - Common samples aligned (projid): %d\n",
    length(common_samples)
  ))

  if (length(common_samples) == 0) {
    cat("    [ERROR] No common samples found between datasets. Skipping.\n")
    return(NULL)
  }

  # Subconjunto de datos de expresión génica
  genes_subset <- genes_df |>
    dplyr::filter(.data$projid %in% common_samples) |>
    dplyr::select(tidyselect::all_of(c("projid", valid_genes))) |>
    dplyr::arrange(.data$projid)

  # Control de calidad: Filtrar genes con varianza cero y al menos 10
  # muestras válidas
  genes_data <- genes_subset |> dplyr::select(-"projid")
  genes_keep <- colnames(genes_data)[sapply(
    genes_data,
    function(x) {
      v <- var(x, na.rm = TRUE)
      n_valid <- sum(!is.na(x))
      !is.na(v) && v > 0 && n_valid >= 10
    }
  )]

  if (length(genes_keep) == 0) {
    cat(paste(
      "    [ERROR] No valid genes with expression",
      "variance and >= 10 samples found. Skipping.\n"
    ))
    return(NULL)
  }
  genes_subset <- genes_subset |>
    dplyr::select(tidyselect::all_of(c("projid", genes_keep)))

  # Identificar columnas de metabolitos reales (excluir lotes y covariables)
  # Conservando ÚNICAMENTE las columnas que comienzan con 'X' y dígitos
  all_cols <- colnames(metab_df)
  metab_cols_raw <- all_cols[grepl("^X\\d+$", all_cols)]

  # Asegurar que sean numéricas en el dataframe
  metab_cols <- metab_cols_raw[sapply(metab_df[metab_cols_raw], is.numeric)]

  if (length(metab_cols) == 0) {
    cat("    [ERROR] No genuine metabolite columns found. Skipping.\n")
    return(NULL)
  }

  # Subconjunto de datos de metabolitos
  metab_subset <- metab_df |>
    dplyr::filter(.data$projid %in% common_samples) |>
    dplyr::select(tidyselect::all_of(c("projid", metab_cols))) |>
    dplyr::arrange(.data$projid)

  # Control de calidad: Filtrar metabolitos con varianza cero y al menos 10
  # muestras válidas
  metab_data <- metab_subset |> dplyr::select(-"projid")
  metab_keep <- colnames(metab_data)[sapply(
    metab_data,
    function(x) {
      v <- var(x, na.rm = TRUE)
      n_valid <- sum(!is.na(x))
      !is.na(v) && v > 0 && n_valid >= 10
    }
  )]

  if (length(metab_keep) == 0) {
    cat(paste(
      "    [ERROR] No valid metabolites with variance and >= 10",
      "samples found. Skipping.\n"
    ))
    return(NULL)
  }
  metab_subset <- metab_subset |>
    dplyr::select(tidyselect::all_of(c("projid", metab_keep)))

  # --- CALCULAR LA MATRIZ DE SPEARMAN Y SIGNIFICACIÓN ESTADÍSTICA ---
  cat("    - Calculating Spearman rho coefficients...\n")
  cor_mat <- cor(
    genes_subset |> dplyr::select(-"projid"),
    metab_subset |> dplyr::select(-"projid"),
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  # Matriz de indicadores de presencia (1 si no es NA, 0 si es NA) para
  # el cálculo exacto de N por par
  g_present <- !is.na(as.matrix(genes_subset |> dplyr::select(-"projid")))
  m_present <- !is.na(as.matrix(metab_subset |> dplyr::select(-"projid")))

  # Convertir a numérico (TRUE -> 1, FALSE -> 0)
  g_present <- g_present + 0
  m_present <- m_present + 0

  # Multiplicación matricial para obtener el N exacto de observaciones
  # completas por par
  n_mat_size <- t(g_present) %*% m_present

  cat(paste(
    "    - Computing statistical significance using Student's t",
    "approximation and direct survival function (lower.tail = FALSE)...\n"
  ))
  # Cálculo vectorial de p-valores usando la aproximación t de Student
  # con el N específico de cada par
  t_mat <- cor_mat * sqrt((n_mat_size - 2) / (1 - cor_mat^2 + 1e-16))
  # Usar lower.tail = FALSE en pt() para evitar cancelaciones
  # numéricas que producen p-valores de exactamente 0
  p_mat <- 2 * pt(abs(t_mat), df = n_mat_size - 2, lower.tail = FALSE)

  # Calcular la corrección FDR de Benjamini-Hochberg en toda la matriz
  p_adj_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    ncol = ncol(p_mat)
  )
  rownames(p_adj_mat) <- rownames(p_mat)
  colnames(p_adj_mat) <- colnames(p_mat)

  # Resolver ruta de salida
  out_path <- if (is_complementary) {
    file.path(base_out_dir, "Complementary_Subnucleus")
  } else {
    base_out_dir
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  # --- GUARDAR LA MATRIZ DE CORRELACIÓN Y LAS TABLAS EN FORMATO LARGO ---

  # 1. Guardar en formato de matriz ancha estándar
  cor_df_wide <- as.data.frame(cor_mat) |> tibble::rownames_to_column("Gene")
  csv_file_wide <- sprintf(
    "correlations_%s_metabolites_matrix.csv",
    tolower(str_replace_all(prefix_name, " ", "_"))
  )
  write_csv(cor_df_wide, file.path(out_path, csv_file_wide))

  # 2. Guardar en formato largo detallado (Gene, Metabolite, r, pvalue, FDR)
  cor_long <- as.data.frame.table(cor_mat, responseName = "r") |>
    dplyr::rename(Gene = "Var1", Metabolite = "Var2")
  p_long <- as.data.frame.table(p_mat, responseName = "pvalue") |>
    dplyr::rename(Gene = "Var1", Metabolite = "Var2")
  padj_long <- as.data.frame.table(p_adj_mat, responseName = "FDR") |>
    dplyr::rename(Gene = "Var1", Metabolite = "Var2")

  # Unir y ordenar por significación
  cor_results <- cor_long |>
    dplyr::left_join(p_long, by = c("Gene", "Metabolite")) |>
    dplyr::left_join(padj_long, by = c("Gene", "Metabolite")) |>
    dplyr::mutate(
      Gene = as.character(.data$Gene),
      Metabolite = as.character(.data$Metabolite)
    ) |>
    dplyr::arrange(.data$FDR, .data$pvalue)

  csv_file_full <- sprintf(
    "correlations_%s_metabolites_full.csv",
    tolower(str_replace_all(prefix_name, " ", "_"))
  )
  write_csv(cor_results, file.path(out_path, csv_file_full))
  cat(sprintf(
    "    [SUCCESS] Full correlation table saved in: %s\n",
    file.path(out_path, csv_file_full)
  ))

  # 3. Guardar interacciones significativas (|r| >= 0.20 y pvalue < 0.05)
  significant_corrs <- cor_results |>
    dplyr::filter(.data$pvalue < 0.05 & abs(.data$r) >= 0.20)

  csv_file_sig <- sprintf(
    "correlations_%s_metabolites_significant.csv",
    tolower(str_replace_all(prefix_name, " ", "_"))
  )
  write_csv(significant_corrs, file.path(out_path, csv_file_sig))
  cat(sprintf(
    "    [SUCCESS] Significant correlations table saved (N = %d): %s\n",
    nrow(significant_corrs), file.path(out_path, csv_file_sig)
  ))

  # --- SELECCIONAR METABOLITOS PARA EL GRÁFICO DE HEATMAP ---
  # Filtrar metabolitos con al menos una correlación moderada y significativa
  sig_mask <- apply(abs(cor_mat) >= 0.25 & p_mat < 0.05, 2, any)
  top_metabs <- names(sig_mask[sig_mask])

  # Alternativas (fallbacks) para la representación gráfica
  max_cor_by_metab <- apply(abs(cor_mat), 2, max, na.rm = TRUE)
  if (length(top_metabs) == 0) {
    cat("    - No metabolites with |r| > 0.25 & p < 0.05. Using top 40.\n")
    top_metabs <- names(sort(max_cor_by_metab, decreasing = TRUE)[
      seq_len(min(40, length(max_cor_by_metab)))
    ])
  } else if (length(top_metabs) > 40) {
    cat(sprintf(
      "    - Found %d significant metabolites. Selecting top 40 by p-value.\n",
      length(top_metabs)
    ))
    min_p_by_metab <- apply(
      p_mat[, top_metabs, drop = FALSE],
      2, min, na.rm = TRUE
    )
    top_metabs <- names(sort(min_p_by_metab)[
      seq_len(min(40, length(min_p_by_metab)))
    ])
  } else {
    cat(sprintf(
      "    - Selected %d significant metabolites for visual rendering.\n",
      length(top_metabs)
    ))
  }

  if (length(top_metabs) == 0) {
    cat("    [WARNING] No numeric metabolites available for visualization.\n")
    return(NULL)
  }

  # --- GENERAR MAPA DE CALOR (HEATMAP) PREMIUM ---
  mat_plot <- cor_mat[, top_metabs, drop = FALSE]
  colnames(mat_plot) <- str_trunc(colnames(mat_plot), 40)

  # Dimensionamiento dinámico basado en el tamaño de los datos
  w <- 5 + (ncol(mat_plot) * 0.22)
  w <- min(max(w, 8), 16)
  h <- 4 + (nrow(mat_plot) * 0.22)
  h <- min(max(h, 6), 14)

  heatmap_file_name <- sprintf(
    "heatmap_%s_metabolites.png",
    tolower(str_replace_all(prefix_name, " ", "_"))
  )
  heatmap_path <- file.path(out_path, heatmap_file_name)

  # Solo mostrar números si el gráfico es legible
  display_cell_nums <- (nrow(mat_plot) * ncol(mat_plot) <= 300)

  tryCatch({
    png(
      heatmap_path,
      width = w,
      height = h,
      units = "in",
      res = 300
    )
    p <- pheatmap(
      mat_plot,
      # Paleta RdYlBu para publicaciones
      color = colorRampPalette(
        c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
          "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
      )(100),
      breaks = seq(-0.6, 0.6, length.out = 101),
      cluster_rows = (nrow(mat_plot) >= 2),
      cluster_cols = (ncol(mat_plot) >= 2),
      display_numbers = display_cell_nums,
      number_format = "%.2f",
      fontsize_number = 6.5,
      number_color = "black",
      border_color = "white",
      lwd = 0.5,
      main = sprintf(
        "Spearman Correlation: %s vs Top Metabolites",
        prefix_name
      )
    )

    # Añadir el título de la leyenda dynamically encima de la barra de color
    gt <- p[["gtable"]]
    leg_idx <- which(gt$layout$name == "legend")
    if (length(leg_idx) > 0) {
      leg_row <- gt$layout$t[leg_idx]
      leg_col <- gt$layout$l[leg_idx]

      # Crear el grob del título, centrado y desplazado hacia arriba
      title_grob <- grid::textGrob(
        "Spearman r",
        y = unit(1, "npc") + unit(8, "bigpts"),
        just = c("center", "bottom"),
        gp = grid::gpar(fontsize = 8, fontface = "bold")
      )

      # Desactivar clipping en la celda de la leyenda
      gt$layout$clip[leg_idx] <- "off"

      # Añadir a la celda del legend
      gt <- gtable::gtable_add_grob(
        gt,
        title_grob,
        t = leg_row, l = leg_col
      )
    }

    grid.newpage()
    grid.draw(gt)
    dev.off()
    cat(sprintf(
      "    [SUCCESS] Correlation Heatmap saved in:\n    %s\n",
      heatmap_path
    ))
  }, error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    cat(sprintf(
      "    [ERROR] Heatmap generation failed: %s\n",
      e$message
    ))
  })

  # --- GENERAR GRÁFICO DE BURBUJAS (DOTPLOT) ELEGANTE ---
  # Muestra simultáneamente el tamaño del efecto y la significación
  top_dots <- significant_corrs |>
    dplyr::arrange(.data$pvalue) |>
    head(40)

  if (nrow(top_dots) > 0) {
    dotplot_file_name <- sprintf(
      "dotplot_%s_metabolites.png",
      tolower(str_replace_all(prefix_name, " ", "_"))
    )
    dotplot_path <- file.path(out_path, dotplot_file_name)

    top_dots <- top_dots |>
      dplyr::mutate(
        # Evitar Inf al aplicar -log10 sobre p-valores que son exactamente 0
        log10p = -log10(pmax(.data$pvalue, 1e-300)),
        Significance = ifelse(
          .data$FDR < 0.05,
          "FDR < 0.05",
          "Nominal p < 0.05"
        )
      )

    # Orden jerárquico determinista libre de promedios para factores
    ordered_metabolites <- top_dots |>
      dplyr::group_by(.data$Metabolite) |>
      dplyr::slice_min(order_by = .data$pvalue, n = 1, with_ties = FALSE) |>
      dplyr::arrange(.data$r) |>
      dplyr::pull(.data$Metabolite)

    ordered_genes <- top_dots |>
      dplyr::group_by(.data$Gene) |>
      dplyr::slice_min(order_by = .data$pvalue, n = 1, with_ties = FALSE) |>
      dplyr::arrange(.data$r) |>
      dplyr::pull(.data$Gene)

    top_dots <- top_dots |>
      dplyr::mutate(
        Metabolite = factor(.data$Metabolite, levels = ordered_metabolites),
        Gene = factor(.data$Gene, levels = ordered_genes)
      )

    # Gráfico de burbujas bellamente ordenado (Spearman r con límites c(-1, 1))
    p_dot <- ggplot(
      top_dots,
      aes(
        x = .data$Metabolite,
        y = .data$Gene,
        size = .data$log10p
      )
    ) +
      geom_point(
        alpha = 0.85, stroke = 0.5, shape = 21,
        aes(fill = .data$r), color = "#4B5563"
      ) +
      scale_fill_gradient2(
        low = "#313695", mid = "#FFFFFF", high = "#A50026",
        midpoint = 0, name = "Spearman r", limits = c(-1, 1)
      ) +
      scale_size_continuous(name = "-log10(p-value)", range = c(3, 8)) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(
          angle = 45, hjust = 1, vjust = 1,
          face = "bold", size = 8
        ),
        axis.text.y = element_text(face = "bold", size = 8),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "#4B5563"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        x = "Metabolite (Serum/Brain ID)",
        y = "Gene Symbol",
        title = sprintf("Top Spearman Correlations: %s", prefix_name),
        subtitle = "Top 40 significant multi-omic interactions (p < 0.05)"
      )

    # Guardar dotplot
    ggsave(dotplot_path, plot = p_dot, width = 9, height = 7, dpi = 300)
    cat(sprintf(
      "    [SUCCESS] Correlation Bubble Plot saved in: %s\n",
      dotplot_path
    ))
  }
}

# --- EJECUCIÓN DEL PIPELINE ---
if (!file.exists(genes_file) || !file.exists(metab_file)) {
  cat("  [ERROR] Integrated Arrangements files could not be found.\n")
  cat(sprintf("  Searched paths: %s and %s\n", genes_file, metab_file))
} else {
  cat("\nLoading integrated multi-omics datasets...\n")
  genes_df <- read_csv(genes_file, show_col_types = FALSE)
  metab_df <- read_csv(metab_file, show_col_types = FALSE)

  # 1. Análisis Principal: FIRMA DE CONSENSO TRIPLE CORE
  if (length(triple_core) > 0) {
    run_correlation_analysis(
      genes_df = genes_df,
      metab_df = metab_df,
      gene_list = triple_core,
      prefix_name = "Triple Core",
      base_out_dir = out_dir,
      is_complementary = FALSE
    )
  } else {
    cat("  [WARNING] Triple Core gene list is empty or file not found.\n")
  }

  # 2. Análisis Complementario: META-FIRMA IRREDUCIBLE DEL SUBNÚCLEO
  if (length(subnucleo) > 0) {
    run_correlation_analysis(
      genes_df = genes_df,
      metab_df = metab_df,
      gene_list = subnucleo,
      prefix_name = "Subnucleus",
      base_out_dir = out_dir,
      is_complementary = TRUE
    )
  }
}

cat("\n[COMPLETE] Multi-omics inter-omic correlation finished successfully.\n")
