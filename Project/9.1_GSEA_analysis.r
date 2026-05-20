# Directorio de trabajo dinámico (establece el WD en la ubicación del script)
if (!interactive()) {
  script_path <- sub(
    "--file=", "",
    grep("--file=", commandArgs(), value = TRUE)
  )
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(script_path)))
  }
} else if (requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
}

# =====================================================================
# SCRIPT: ANÁLISIS DE GSEA (ACADEMIC-GRADE Y ROBUSTO)
# =====================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyverse)
  library(fgsea)
  library(data.table)
  library(enrichplot)
  library(dplyr)
  library(patchwork) # Para combinar plots en 2x2
  library(aplot)
  library(stringr)
  library(readr)
})

# Establecer semilla para reproducibilidad (fgsea usa permutaciones aleatorias)
set.seed(42)

# Directorio base
base_dir_gsea <- "results_enrichment/GSEA"
dir.create(
  file.path(base_dir_gsea, "Plots"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_gsea, "Significative"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_gsea, "validated_GSEA_Global"),
  showWarnings = FALSE, recursive = TRUE
)

# Parámetros globales de corte para GSEA
nes_threshold <- 2.0
q_cutoff <- 0.10 # FDR estricto al 10%
# Umbral relajado de simulación para colapsado robusto
pvalue_simulation_cutoff <- 0.25

# --- FUNCIÓN PARA GENERAR EL COMPOSITE 2X2 ---
generate_composite_gsea <- function(obj, prefix, ds_key, plot_dir) {
  cat(
    paste(
      "    - Entrando en generate_composite_gsea para",
      prefix, "-", ds_key, "\n"
    )
  )
  df <- as.data.frame(obj)
  if (nrow(df) < 1) {
    cat("    - [ADVERTENCIA] No hay rutas significativas para graficar.\n")
    return(NULL)
  }

  # 1. Seleccionar IDs siguiendo las reglas del usuario
  # El objeto 'obj' ya viene filtrado por significancia en el script principal
  n_total <- nrow(df)
  top_n <- min(3, n_total)
  top_3_df <- df[1:top_n, ]

  ids_to_plot <- as.character(top_3_df$ID)
  ranks <- seq_len(top_n)

  if (n_total >= 3) {
    all_pos <- all(top_3_df$NES > 0)
    all_neg <- all(top_3_df$NES < 0)

    # Buscar la 4ª ruta (si existe)
    found_4th <- FALSE
    if (all_pos) {
      # Si los 3 son positivos, buscamos la mejor suprimida (NES < 0)
      # significativa
      opp_idx <- which(df$NES < 0)[1]
      if (!is.na(opp_idx)) {
        ids_to_plot <- c(ids_to_plot, as.character(df$ID[opp_idx]))
        ranks <- c(ranks, opp_idx)
        found_4th <- TRUE
        cat(
          sprintf(
            "    - Regla: Top 3 Pos. Neg. Sig. en rango: %d\n",
            opp_idx
          )
        )
      }
    } else if (all_neg) {
      # Si los 3 son negativos, buscamos la mejor enriquecida (NES > 0)
      # significativa
      opp_idx <- which(df$NES > 0)[1]
      if (!is.na(opp_idx)) {
        ids_to_plot <- c(ids_to_plot, as.character(df$ID[opp_idx]))
        ranks <- c(ranks, opp_idx)
        found_4th <- TRUE
        cat(
          sprintf(
            "    - Regla: Top 3 Neg. Pos. Sig. en rango: %d\n",
            opp_idx
          )
        )
      }
    }

    # Si no se encontró por la regla anterior o era variada,
    # tomamos la top 4 si existe
    if (!found_4th && n_total >= 4) {
      ids_to_plot <- c(ids_to_plot, as.character(df$ID[4]))
      ranks <- c(ranks, 4)
      cat("    - Regla: Top 4 estándar.\n")
    }
  }

  # 2. Crear los plots individuales
  plot_list <- list()
  for (i in seq_along(ids_to_plot)) {
    id <- ids_to_plot[i]
    rank <- ranks[i]

    # Verificar si el ID existe en el objeto
    if (!id %in% obj@result$ID) {
      cat(paste("    - [ADVERTENCIA] ID no está en el resultado:", id, "\n"))
      next
    }

    row <- df[df$ID == id, ]
    desc <- if (nrow(row) > 0 && !is.na(row$Description[1])) {
      row$Description[1]
    } else {
      id
    }

    # Leyenda con p.adj, FDR (qvalue) y NES
    stats_text <- paste0(
      "Rank: #", rank, " | p.adj: ",
      formatC(row$p.adjust[1], format = "e", digits = 2),
      " | FDR: ", formatC(row$qvalue[1], format = "e", digits = 2),
      "\nNES: ", round(row$NES[1], 2)
    )

    p_item <- tryCatch(
      {
        gseaplot2(
          obj,
          geneSetID = id,
          title = stringr::str_trunc(desc, 50),
          pvalue_table = FALSE,
          base_size = 9
        )
      },
      error = function(e) {
        cat(paste(
          "    - [ERROR] gseaplot2 falló para", id, ":", e$message, "\n"
        ))
        NULL
      }
    )

    if (!is.null(p_item)) {
      # Inyectar la leyenda y personalizar subplots internos del gglist
      tryCatch(
        {
          # El último subplot del panel es ideal para llevar el caption/leyenda
          last_idx <- length(p_item)
          p_item[[last_idx]] <- p_item[[last_idx]] +
            labs(caption = stats_text) +
            theme(
              plot.caption = element_text(
                hjust = 0.5, size = 8,
                face = "italic", color = "darkblue"
              )
            )

          # El primer subplot (arriba) lleva el título
          p_item[[1]] <- p_item[[1]] +
            theme(plot.title = element_text(size = 10, face = "bold"))

          plot_list[[length(plot_list) + 1]] <- p_item
          cat(paste("    - Gráfico añadido con éxito para el ID:", id, "\n"))
        },
        error = function(e) {
          cat(
            sprintf(
              "    - [ERROR] Personalización falló para %s: %s\n",
              id, e$message
            )
          )
        }
      )
    }
  }

  # 3. Combinar en 2x2
  cat(paste("    - Se encontraron", length(plot_list), "gráficos.\n"))
  if (length(plot_list) == 0) {
    return(NULL)
  }

  tryCatch(
    {
      cat("    - Generando gráfico compuesto...\n")
      # aplot::plot_list combina con éxito paneles complejos del tipo gglist
      composite <- aplot::plot_list(gglist = plot_list, ncol = 2)

      filename <- file.path(
        plot_dir,
        paste0("composite_GSEA_", prefix, "_", ds_key, ".png")
      )
      ggsave(filename, plot = composite, width = 14, height = 11, dpi = 300)
      cat(paste("    - [ÉXITO] Gráfico compuesto guardado:", filename, "\n"))
    },
    error = function(e) {
      cat(
        paste(
          "    - [ERROR] aplot::plot_list o ggsave falló:",
          e$message, "\n"
        )
      )
    }
  )
}

# --- FUNCIONES DE SOPORTE PARA GRÁFICOS ---
create_gsea_barplot <- function(obj, prefix, ds_key, base_dir) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0) {
    return()
  }

  # Tomamos top 20 por valor absoluto de NES
  df_plot <- df |>
    arrange(desc(abs(.data$NES))) |>
    head(20)

  p <- ggplot(
    df_plot,
    aes(
      x = reorder(.data$Description, .data$NES),
      y = .data$NES,
      fill = .data$NES > 0
    )
  ) +
    geom_col(width = 0.8, color = "black", linewidth = 0.25) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray40", linewidth = 0.5
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#d63031", "FALSE" = "#0984e3"),
      labels = c(
        "TRUE" = "Enriched (NES > 0)",
        "FALSE" = "Suppressed (NES < 0)"
      )
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste("Análisis de Enriquecimiento GSEA:", prefix, "-", ds_key),
      subtitle = "Top 20 rutas biológicas ordenadas por NES absoluto",
      x = "Ruta Biológica (Pathway)", y = "Normalized Enrichment Score (NES)",
      fill = "Dirección"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(
        face = "italic", size = 10,
        hjust = 0.5, color = "gray30"
      ),
      axis.text.y = element_text(size = 9, face = "bold", color = "black"),
      axis.text.x = element_text(size = 9, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  filename <- file.path(
    base_dir, "Plots",
    paste0("fold_enrichment_barplot_", prefix, "_", ds_key, ".png")
  )
  ggsave(filename, plot = p, width = 11, height = 9, dpi = 300)
}

# --- FUNCIONES DE SOPORTE ---
generate_plots_gsea <- function(obj, prefix, ds_key, base_dir,
                                fold_change_vector = NULL) {
  if (is.null(obj) || nrow(obj@result) == 0) {
    return()
  }
  plot_dir <- file.path(base_dir, "Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Dotplot
  try(
    {
      p1 <- dotplot(obj, showCategory = 15, split = ".sign") +
        facet_grid(. ~ .sign) +
        theme_bw()
      ggsave(
        file.path(plot_dir, paste0("dotplot_", prefix, "_", ds_key, ".png")),
        plot = p1, width = 10, height = 8
      )
    },
    silent = FALSE
  )

  # Emapplot (Corregido para usar similitud Jaccard "JC"
  # en todas las categorías)
  # Esto soluciona el fallo silencioso del método "Wang" que requería semData.
  try(
    {
      if (prefix != "GLOBAL" && nrow(obj@result) >= 2) {
        obj_sim <- pairwise_termsim(obj, method = "JC")
        p2 <- emapplot(obj_sim, showCategory = 20) +
          ggtitle(paste("Mapa de Enriquecimiento GSEA:", prefix, "-", ds_key))
        ggsave(
          file.path(plot_dir, paste0("emapplot_", prefix, "_", ds_key, ".png")),
          plot = p2, width = 12, height = 10
        )
      }
    },
    silent = FALSE
  )

  # Cnetplot (Coloreado premium de genes utilizando sus t-statistics)
  try(
    {
      p3 <- cnetplot(obj, showCategory = 5, foldChange = fold_change_vector) +
        ggtitle(paste("Red Gen-Ruta GSEA:", prefix, "-", ds_key))
      ggsave(
        file.path(plot_dir, paste0("cnetplot_", prefix, "_", ds_key, ".png")),
        plot = p3, width = 12, height = 10
      )
    },
    silent = FALSE
  )

  # GENERAR COMPOSITE 2X2
  try(
    {
      generate_composite_gsea(obj, prefix, ds_key, plot_dir)
    },
    silent = FALSE
  )

  # Barplot de NES Premium
  try(
    {
      create_gsea_barplot(obj, prefix, ds_key, base_dir)
    },
    silent = FALSE
  )
}

update_gsea_result <- function(obj, df) {
  df_df <- as.data.frame(df)
  if (nrow(df_df) > 0) {
    rownames(df_df) <- df_df$ID
  }
  obj@result <- df_df
  obj
}

save_and_plot_gsea <- function(obj, prefix, ds_key, val_genes,
                               fold_change_vector = NULL) {
  if (is.null(obj)) {
    return(NULL)
  }

  df <- as.data.frame(obj)
  if (nrow(df) == 0) {
    return(NULL)
  }

  cat(paste("    - GSEA", prefix, ":", nrow(df), "rutas encontradas.\n"))
  df$Category <- prefix

  # Identificar genes validados en el Core Enrichment
  df <- df |>
    mutate(
      gene_list = strsplit(as.character(.data$core_enrichment), "/")
    ) |>
    rowwise() |>
    mutate(
      Validated_Genes_Found = paste(
        intersect(.data$gene_list, val_genes),
        collapse = "/"
      ),
      n_Validated = length(intersect(.data$gene_list, val_genes))
    ) |>
    ungroup() |>
    dplyr::select(-.data$gene_list)

  # Actualizamos el objeto
  obj <- update_gsea_result(obj, df)

  # 1. Guardar en Significativas (Significative)
  write_csv(
    df,
    file.path(
      base_dir_gsea, "Significative",
      paste0("gsea_", prefix, "_", ds_key, ".csv")
    )
  )
  generate_plots_gsea(obj, prefix, ds_key, base_dir_gsea, fold_change_vector)

  # 2. Guardar en Validadas (validated_GSEA_Global)
  df_val <- df |> filter(.data$n_Validated > 0)
  if (nrow(df_val) > 0) {
    write_csv(
      df_val,
      file.path(
        base_dir_gsea, "validated_GSEA_Global",
        paste0("validated_gsea_", prefix, "_", ds_key, ".csv")
      )
    )
    obj_val <- obj
    obj_val <- update_gsea_result(obj_val, df_val)
    generate_plots_gsea(
      obj_val, prefix, ds_key,
      file.path(base_dir_gsea, "validated_GSEA_Global"),
      fold_change_vector
    )
  }

  df
}

# --- PROCESAMIENTO POR DATASET ---
for (ds in c("braaksc", "ceradsc", "cogdx")) {
  cat(paste("\n>>> PROCESANDO GSEA:", toupper(ds), "\n"))
  path_dea <- paste0("differential_expression/DEA_", ds, "_complete.csv")

  if (!file.exists(path_dea)) {
    cat(paste("Error: Archivo no encontrado:", path_dea, "\n"))
    next
  }

  dea_results <- read_csv(path_dea, show_col_types = FALSE)
  is_val <- dea_results$Validated == TRUE | dea_results$Validated == "TRUE"
  val_genes <- dea_results$Symbol[is_val & !is.na(dea_results$Symbol)]

  cat("    - Mapeando símbolos a Entrez...\n")
  entrez_map <- suppressMessages(
    bitr(
      dea_results$Symbol,
      fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  rank_df <- dea_results |>
    inner_join(entrez_map, by = c("Symbol" = "SYMBOL")) |>
    filter(!is.na(.data$t), !is.na(.data$ENTREZID)) |>
    dplyr::select(.data$ENTREZID, .data$t) |>
    group_by(.data$ENTREZID) |>
    summarize(t = mean(.data$t)) |>
    arrange(desc(.data$t))

  ranked_list <- rank_df$t
  names(ranked_list) <- rank_df$ENTREZID

  if (length(ranked_list) == 0) {
    next
  }

  # Crear vector de t-statistics con símbolos para colorear cnetplot
  symbols_t_df <- dea_results |>
    filter(!is.na(Symbol), !is.na(t)) |>
    group_by(Symbol) |>
    summarize(t = mean(t))
  fold_change_vector <- symbols_t_df$t
  names(fold_change_vector) <- symbols_t_df$Symbol

  res_bp <- NULL
  res_mf <- NULL
  res_kegg <- NULL
  res_bp_obj <- NULL
  res_mf_obj <- NULL
  res_kegg_obj <- NULL

  # GSEA BP
  cat("    - Ejecutando GSEA GO BP...\n")
  gsea_bp <- gseGO(
    ranked_list,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    # Simulación relajada para colapsado robusto
    pvalueCutoff  = pvalue_simulation_cutoff,
    pAdjustMethod = "BH",
    minGSSize     = 15, # Filtro biológico mínimo
    maxGSSize     = 500, # Filtro biológico máximo
    eps           = 0, # Precisión exacta de p-valores
    verbose       = FALSE
  )
  if (!is.null(gsea_bp) && nrow(as.data.frame(gsea_bp)) > 0) {
    gs <- geneInCategory(gsea_bp)

    # CORRECCIÓN DEFENSIVA: Asegurar que todos los genes de las rutas estén
    # contenidos estrictamente en ranked_list para evitar caídas
    # en collapsePathways
    gs <- lapply(gs, function(g) intersect(g, names(ranked_list)))

    fg <- as.data.table(as.data.frame(gsea_bp))
    setnames(
      fg,
      old = c("ID", "pvalue", "enrichmentScore"),
      new = c("pathway", "pval", "ES")
    )
    main <- collapsePathways(fg[order(pval)], gs, ranked_list)
    filtered_bp <- gsea_bp@result |>
      filter(
        .data$ID %in% main$mainPathways,
        abs(.data$NES) > nes_threshold,
        .data$qvalue < q_cutoff
      )
    gsea_bp <- update_gsea_result(gsea_bp, filtered_bp)

    if (nrow(filtered_bp) > 0) {
      res_bp_obj <- setReadable(gsea_bp, org.Hs.eg.db, "ENTREZID")
      res_bp <- save_and_plot_gsea(
        res_bp_obj, "GO_BP", ds, val_genes,
        fold_change_vector
      )
    }
  }

  # GSEA MF
  cat("    - Ejecutando GSEA GO MF...\n")
  gsea_mf <- gseGO(
    ranked_list,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pvalueCutoff  = pvalue_simulation_cutoff, # Simulación relajada
    pAdjustMethod = "BH",
    minGSSize     = 15,
    maxGSSize     = 500,
    eps           = 0,
    verbose       = FALSE
  )
  if (!is.null(gsea_mf) && nrow(as.data.frame(gsea_mf)) > 0) {
    gs_m <- geneInCategory(gsea_mf)

    # CORRECCIÓN DEFENSIVA: Asegurar que todos los genes de las rutas estén
    # contenidos estrictamente en ranked_list para evitar caídas
    # en collapsePathways
    gs_m <- lapply(gs_m, function(g) intersect(g, names(ranked_list)))

    fg_m <- as.data.table(as.data.frame(gsea_mf))
    setnames(
      fg_m,
      old = c("ID", "pvalue", "enrichmentScore"),
      new = c("pathway", "pval", "ES")
    )
    main_m <- collapsePathways(fg_m[order(pval)], gs_m, ranked_list)
    filtered_mf <- gsea_mf@result |>
      filter(
        .data$ID %in% main_m$mainPathways,
        abs(.data$NES) > nes_threshold,
        .data$qvalue < q_cutoff
      )
    gsea_mf <- update_gsea_result(gsea_mf, filtered_mf)

    if (nrow(filtered_mf) > 0) {
      res_mf_obj <- setReadable(gsea_mf, org.Hs.eg.db, "ENTREZID")
      res_mf <- save_and_plot_gsea(
        res_mf_obj, "GO_MF", ds, val_genes,
        fold_change_vector
      )
    }
  }

  # GSEA KEGG
  cat("    - Ejecutando GSEA KEGG...\n")
  gsea_kegg <- gseKEGG(
    ranked_list,
    organism      = "hsa",
    pvalueCutoff  = pvalue_simulation_cutoff, # Simulación relajada
    pAdjustMethod = "BH",
    minGSSize     = 15,
    maxGSSize     = 500,
    eps           = 0,
    verbose       = FALSE
  )
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    gs_k <- geneInCategory(gsea_kegg)

    # CORRECCIÓN DEFENSIVA: Asegurar que todos los genes de las rutas estén
    # contenidos estrictamente en ranked_list para evitar caídas
    # en collapsePathways
    gs_k <- lapply(gs_k, function(g) intersect(g, names(ranked_list)))

    fg_k <- as.data.table(as.data.frame(gsea_kegg))
    setnames(
      fg_k,
      old = c("ID", "pvalue", "enrichmentScore"),
      new = c("pathway", "pval", "ES")
    )
    main_k <- collapsePathways(fg_k[order(pval)], gs_k, ranked_list)
    filtered_kegg <- gsea_kegg@result |>
      filter(
        .data$ID %in% main_k$mainPathways,
        abs(.data$NES) > nes_threshold,
        .data$qvalue < q_cutoff
      )
    gsea_kegg <- update_gsea_result(gsea_kegg, filtered_kegg)

    if (nrow(filtered_kegg) > 0) {
      res_kegg_obj <- setReadable(gsea_kegg, org.Hs.eg.db, "ENTREZID")
      res_kegg <- save_and_plot_gsea(
        res_kegg_obj, "KEGG", ds, val_genes,
        fold_change_vector
      )
    }
  }

  # --- GENERACIÓN DE ARCHIVOS GLOBALES ---
  res_list <- list(res_bp, res_mf, res_kegg)
  valid_dfs <- res_list[
    sapply(res_list, function(x) !is.null(x) && nrow(x) > 0)
  ]

  objs_list <- list(res_bp_obj, res_mf_obj, res_kegg_obj)
  valid_objs <- objs_list[
    sapply(objs_list, function(x) !is.null(x) && nrow(as.data.frame(x)) > 0)
  ]

  if (length(valid_dfs) > 0) {
    global_df <- bind_rows(valid_dfs)
    write_csv(
      global_df,
      file.path(
        base_dir_gsea, "Significative",
        paste0("gsea_GLOBAL_", ds, ".csv")
      )
    )

    df_val_glob <- global_df |>
      filter(.data$n_Validated > 0)
    if (nrow(df_val_glob) > 0) {
      write_csv(
        df_val_glob,
        file.path(
          base_dir_gsea, "validated_GSEA_Global",
          paste0("validated_gsea_GLOBAL_", ds, ".csv")
        )
      )
    }
  }
}

cat("\nReorganización de GSEA y gráficos compuestos completados.\n")
