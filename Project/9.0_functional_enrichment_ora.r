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
# ANÁLISIS DE SOBRE-REPRESENTACIÓN (ORA) PRINCIPAL - functional_enrichment_ora.r
# =====================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyverse)
  library(enrichplot)
  library(dplyr)
  library(readr)
})

# Directorio base
base_dir_ora <- "results_enrichment/ORA"
dir.create(
  file.path(base_dir_ora, "Plots"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_ora, "Significative"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_ora, "validated_ORA_Global"),
  showWarnings = FALSE, recursive = TRUE
)

# Nuevo Directorio base para genes de ML Validados
base_dir_ml <- "results_enrichment/ORA/ORA_validated"
dir.create(
  file.path(base_dir_ml, "Plots"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_ml, "Significative"),
  showWarnings = FALSE, recursive = TRUE
)

# --- FUNCIONES DE SOPORTE ---
symbols_to_entrez <- function(symbols) {
  symbols <- as.character(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != "" & symbols != "NA"]
  if (length(symbols) == 0) {
    character(0)
  } else {
    mapping <- suppressMessages(
      bitr(
        symbols,
        fromType = "SYMBOL", toType = "ENTREZID",
        OrgDb = org.Hs.eg.db, drop = TRUE
      )
    )
    unique(mapping$ENTREZID)
  }
}

generate_plots_ora <- function(obj, prefix, ds_key, base_dir) {
  if (is.null(obj) || nrow(obj@result) == 0) {
    return()
  }
  plot_dir <- file.path(base_dir, "Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # 3. Suite de Gráficos (Usamos try() independientes para evitar caídas)

  # Dotplot
  try({
    p1 <- dotplot(obj, showCategory = 15) +
      ggtitle(paste(prefix, "-", ds_key))
    ggsave(
      file.path(plot_dir, paste0("dotplot_", prefix, "_", ds_key, ".png")),
      plot = p1, width = 10, height = 8
    )
  }, silent = TRUE)

  # Emapplot (requiere similitud y >= 2 categorías.
  # Evitamos el plot GLOBAL por solicitud del usuario;
  # si es KEGG, usamos similitud Jaccard JC)
  try({
    if (prefix != "GLOBAL" && nrow(obj@result) >= 2) {
      method_sim <- if (prefix == "KEGG") "JC" else "Wang"
      obj_sim <- pairwise_termsim(obj, method = method_sim)
      p2 <- emapplot(obj_sim, showCategory = 20) +
        ggtitle(paste("Enrichment Map", prefix, "-", ds_key))
      ggsave(
        file.path(plot_dir, paste0("emapplot_", prefix, "_", ds_key, ".png")),
        plot = p2, width = 12, height = 10
      )
    }
  }, silent = TRUE)

  # Cnetplot (Red de relaciones entre genes y rutas)
  try({
    p3 <- cnetplot(obj, showCategory = 5) +
      ggtitle(paste("Gene-Pathway Network", prefix, "-", ds_key))
    ggsave(
      file.path(plot_dir, paste0("cnetplot_", prefix, "_", ds_key, ".png")),
      plot = p3, width = 12, height = 10
    )
  }, silent = TRUE)

  # 4. Nuevos Gráficos: Barplot de Fold Enrichment
  try({
    create_ora_barplot(obj, prefix, ds_key, base_dir)
  }, silent = TRUE)
}

save_and_plot_ora <- function(
  obj, prefix, ds_key, val_genes, base_dir = base_dir_ora
) {
  if (is.null(obj)) {
    return(NULL)
  }

  df <- as.data.frame(obj)
  if (nrow(df) == 0) {
    return(NULL)
  }

  cat(paste("    -", prefix, ":", nrow(df), "rutas encontradas.\n"))
  df$Category <- prefix

  # Identificar genes validados
  df <- df |>
    mutate(gene_list = strsplit(as.character(.data$geneID), "/")) |>
    rowwise() |>
    mutate(
      Validated_Genes_Found = paste(
        intersect(.data$gene_list, val_genes),
        collapse = "/"
      ),
      n_Validated = length(intersect(.data$gene_list, val_genes))
    ) |>
    ungroup() |>
    dplyr::select(-"gene_list")

  # Actualizamos el objeto con el dataframe anotado
  obj@result <- df

  # 1. Guardar en Significativas (Significative)
  write_csv(
    df,
    file.path(
      base_dir, "Significative",
      paste0("ora_", prefix, "_", ds_key, ".csv")
    )
  )
  generate_plots_ora(obj, prefix, ds_key, base_dir)

  # 2. Guardar y graficar en Valificados (validated_ORA_Global)
  df_val <- df |> filter(.data$n_Validated > 0)
  if (nrow(df_val) > 0) {
    # Solo duplicamos en la carpeta validated_ORA_Global
    # para la ejecución del DEA general
    if (base_dir == base_dir_ora) {
      write_csv(
        df_val,
        file.path(
          base_dir, "validated_ORA_Global",
          paste0("validated_ora_", prefix, "_", ds_key, ".csv")
        )
      )
      obj_val <- obj
      obj_val@result <- df_val
      generate_plots_ora(
        obj_val, prefix, ds_key,
        file.path(base_dir, "validated_ORA_Global")
      )
    }
  }

  df
}

# --- FUNCIONES DE SOPORTE PARA GRÁFICOS NUEVOS ---
create_ora_barplot <- function(obj, prefix, ds_key, base_dir) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0) {
    return()
  }

  # Analizar proporciones para calcular el enriquecimiento (Fold Enrichment)
  parse_ratio <- function(x) {
    sapply(strsplit(as.character(x), "/"), function(y) {
      if (length(y) < 2) {
        NA
      } else {
        as.numeric(y[1]) / as.numeric(y[2])
      }
    })
  }

  df$FoldEnrichment <- parse_ratio(df$GeneRatio) / parse_ratio(df$BgRatio)

  # Tomamos top 20 por Fold Enrichment
  df_plot <- df |>
    arrange(desc(.data$FoldEnrichment)) |>
    head(20)

  p <- ggplot(
    df_plot,
    aes(
      x = reorder(.data$Description, .data$FoldEnrichment),
      y = .data$FoldEnrichment,
      fill = .data$FoldEnrichment
    )
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme_minimal() +
    labs(
      title = paste("Fold Enrichment:", prefix, "-", ds_key),
      subtitle = "Top 20 routes by Fold Enrichment",
      x = "Pathway", y = "Fold Enrichment"
    ) +
    theme(
      axis.text.y = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 14, face = "bold")
    )

  filename <- file.path(
    base_dir, "Plots",
    paste0("fold_enrichment_barplot_", prefix, "_", ds_key, ".png")
  )
  ggsave(filename, plot = p, width = 11, height = 9, dpi = 300)
}

# --- DEFINICIÓN EL UNIVERSO DE GENES DE FONDO ---
cat("\n>>> DEFINIENDO EL UNIVERSO DE GENES DE FONDO DESDE DATA...\n")
# Leemos únicamente las columnas de los tres archivos de entrada
# (sin cargar datos) para extraer la unión de genes iniciales.
raw_cols_braak <- colnames(
  read_csv("data/genes-braaksc.csv", n_max = 0, show_col_types = FALSE)
)
raw_cols_cerad <- colnames(
  read_csv("data/genes-ceradsc.csv", n_max = 0, show_col_types = FALSE)
)
raw_cols_cogdx <- colnames(
  read_csv("data/genes-cogdx.csv", n_max = 0, show_col_types = FALSE)
)

meta_cols <- c(
  "projid", "individualID", "SampleID",
  "braaksc", "ceradsc", "cogdx"
)

primary_symbols <- unique(c(
  setdiff(raw_cols_braak, meta_cols),
  setdiff(raw_cols_cerad, meta_cols),
  setdiff(raw_cols_cogdx, meta_cols)
))

entrez_universe <- symbols_to_entrez(primary_symbols)
cat(
  sprintf(
    paste0(
      "    - Universo de genes de partida definido (unión de 3 archivos): ",
      "%d símbolos a %d Entrez.\n"
    ),
    length(primary_symbols), length(entrez_universe)
  )
)

# --- PROCESAMIENTO POR DATASET ---
for (ds in c("braaksc", "ceradsc", "cogdx")) {
  cat(paste("\n>>> PROCESANDO ORA:", toupper(ds), "\n"))
  path_dea <- paste0("differential_expression/DEA_", ds, "_complete.csv")

  if (!file.exists(path_dea)) {
    cat(paste("Error: Archivo no encontrado:", path_dea, "\n"))
    next
  }

  dea_data <- read_csv(path_dea, show_col_types = FALSE)

  is_sig <- dea_data$Significant == TRUE | dea_data$Significant == "TRUE"
  is_val <- dea_data$Validated == TRUE | dea_data$Validated == "TRUE"

  val_genes <- dea_data$Symbol[is_val & !is.na(dea_data$Symbol)]
  degs_symbols <- dea_data$Symbol[is_sig & !is.na(dea_data$Symbol)]

  cat(paste("    - Total de DEGs encontrados:", length(degs_symbols), "\n"))
  cat(paste("    - Total de genes validados:", length(val_genes), "\n"))

  entrez_degs <- symbols_to_entrez(degs_symbols)
  cat(paste("    - DEGs mapeados (Entrez):", length(entrez_degs), "\n"))

  if (length(entrez_degs) == 0) {
    cat("    - Omitir: No hay DEGs mapeados a Entrez.\n")
    next
  }

  # NOTA SOBRE EL UNIVERSO DE GENES:
  # Restringimos el universo de ORA exactamente a los 12,556 genes iniciales
  # medidos en la carpeta data/ (transcriptoma de partida de ROSMAP).
  # Esto asegura que el análisis hipergeométrico se realice sobre el trasfondo
  # correcto de la muestra.

  # GO BP
  cat("    - Ejecutando GO BP...\n")
  res_bp_obj <- enrichGO(
    gene = entrez_degs,
    universe = entrez_universe,
    OrgDb = org.Hs.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.10,
    qvalueCutoff = 0.25, readable = TRUE
  )
  res_bp <- save_and_plot_ora(res_bp_obj, "GO_BP", ds, val_genes)

  # GO MF
  cat("    - Ejecutando GO MF...\n")
  res_mf_obj <- enrichGO(
    gene = entrez_degs,
    universe = entrez_universe,
    OrgDb = org.Hs.eg.db, ont = "MF",
    pAdjustMethod = "BH", pvalueCutoff = 0.10,
    qvalueCutoff = 0.25, readable = TRUE
  )
  res_mf <- save_and_plot_ora(res_mf_obj, "GO_MF", ds, val_genes)

  # KEGG
  cat("    - Ejecutando KEGG...\n")
  res_kegg_obj <- enrichKEGG(
    gene = entrez_degs,
    universe = entrez_universe,
    organism = "hsa", pAdjustMethod = "BH",
    pvalueCutoff = 0.10, qvalueCutoff = 0.25
  )
  if (!is.null(res_kegg_obj) && nrow(as.data.frame(res_kegg_obj)) > 0) {
    res_kegg_obj <- setReadable(
      res_kegg_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
    )
    res_kegg_df <- save_and_plot_ora(res_kegg_obj, "KEGG", ds, val_genes)
  } else {
    res_kegg_df <- NULL
  }

  # --- GENERACIÓN DE ARCHIVOS GLOBALES ---
  res_list <- list(res_bp, res_mf, res_kegg_df)
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
        base_dir_ora, "Significative",
        paste0("ora_GLOBAL_", ds, ".csv")
      )
    )

    df_val_glob <- global_df |>
      filter(.data$n_Validated > 0)

    if (nrow(df_val_glob) > 0) {
      write_csv(
        df_val_glob,
        file.path(
          base_dir_ora, "validated_ORA_Global",
          paste0("validated_ora_GLOBAL_", ds, ".csv")
        )
      )
    }
  } else {
    cat(paste("    - Sin resultados significativos de ORA para", ds, "\n"))
  }

  # =====================================================================
  # ORA INDEPENDIENTE SOBRE GENES VALIDADOS (ML + DEA)
  # =====================================================================
  cat("    - Iniciando ORA independiente sobre genes validados (ML + DEA)...\n")
  entrez_ml_val <- symbols_to_entrez(val_genes)

  if (length(entrez_ml_val) > 0) {
    # GO BP para validados
    cat("      * Ejecutando GO BP para genes validados...\n")
    res_bp_obj_ml <- enrichGO(
      gene = entrez_ml_val,
      universe = entrez_universe,
      OrgDb = org.Hs.eg.db, ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 0.10,
      qvalueCutoff = 0.25, readable = TRUE
    )
    res_bp_ml <- save_and_plot_ora(
      res_bp_obj_ml, "GO_BP", ds, val_genes,
      base_dir = base_dir_ml
    )

    # GO MF para validados
    cat("      * Ejecutando GO MF para genes validados...\n")
    res_mf_obj_ml <- enrichGO(
      gene = entrez_ml_val,
      universe = entrez_universe,
      OrgDb = org.Hs.eg.db, ont = "MF",
      pAdjustMethod = "BH", pvalueCutoff = 0.10,
      qvalueCutoff = 0.25, readable = TRUE
    )
    res_mf_ml <- save_and_plot_ora(
      res_mf_obj_ml, "GO_MF", ds, val_genes,
      base_dir = base_dir_ml
    )

    # KEGG para validados
    cat("      * Ejecutando KEGG para genes validados...\n")
    res_kegg_obj_ml <- enrichKEGG(
      gene = entrez_ml_val,
      universe = entrez_universe,
      organism = "hsa", pAdjustMethod = "BH",
      pvalueCutoff = 0.10, qvalueCutoff = 0.25
    )
    if (
      !is.null(res_kegg_obj_ml) &&
        nrow(as.data.frame(res_kegg_obj_ml)) > 0
    ) {
      res_kegg_obj_ml <- setReadable(
        res_kegg_obj_ml, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
      )
      res_kegg_df_ml <- save_and_plot_ora(
        res_kegg_obj_ml, "KEGG", ds, val_genes,
        base_dir = base_dir_ml
      )
    } else {
      res_kegg_df_ml <- NULL
    }

    # --- GENERACIÓN DE ARCHIVOS GLOBALES PARA ML ---
    res_list_ml <- list(res_bp_ml, res_mf_ml, res_kegg_df_ml)
    valid_dfs_ml <- res_list_ml[
      sapply(res_list_ml, function(x) !is.null(x) && nrow(x) > 0)
    ]

    objs_list_ml <- list(res_bp_obj_ml, res_mf_obj_ml, res_kegg_obj_ml)
    valid_objs_ml <- objs_list_ml[
      sapply(
        objs_list_ml,
        function(x) !is.null(x) && nrow(as.data.frame(x)) > 0
      )
    ]

    if (length(valid_dfs_ml) > 0) {
      global_df_ml <- bind_rows(valid_dfs_ml)
      write_csv(
        global_df_ml,
        file.path(
          base_dir_ml, "Significative",
          paste0("ora_GLOBAL_", ds, ".csv")
        )
      )

      cat(paste(
        "      * [ÉXITO] ORA independiente para genes validados guardado",
        "en 'results_enrichment/ORA/ORA_validated'.\n"
      ))
    } else {
      cat(paste(
        "      * [INFO] Sin rutas significativas en el ORA independiente",
        "de genes validados.\n"
      ))
    }
  } else {
    cat(paste(
      "      * [INFO] No hay genes validados mapeados a Entrez",
      "para ejecutar ORA independiente.\n"
    ))
  }
}

cat("\nEjecución del script maestro de ORA completada.\n")
