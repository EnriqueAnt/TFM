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

# ====================================================================
# SCRIPT 10 - ENTORNO DE CONJUNTOS DE GENES (GENE SET LANDSCAPE)
# ====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(stringr)
  library(clusterProfiler)
  library(enrichplot)
  library(readr)
})

base_dir_landscape <- "results_enrichment/Landscape"
dir.create(
  file.path(base_dir_landscape, "Data"),
  showWarnings = FALSE, recursive = TRUE
)
dir.create(
  file.path(base_dir_landscape, "Plots"),
  showWarnings = FALSE, recursive = TRUE
)

safe_read <- function(path) {
  if (file.exists(path)) {
    lines <- readLines(path) |> str_trim()
    lines[lines != ""]
  } else {
    character(0)
  }
}

# 1. Definir los genes de interés
subnucleo <- safe_read(
  "biological_analysis/gene_lists/08_genes_subnucleus.txt"
)

triple_core <- safe_read(
  "biological_analysis/gene_lists/01_genes_core_masters.txt"
)
braak_unique <- safe_read(
  "biological_analysis/gene_lists/02_genes_unique_braak.txt"
)
cerad_unique <- safe_read(
  "biological_analysis/gene_lists/03_genes_unique_cerad.txt"
)
cogdx_unique <- safe_read(
  "biological_analysis/gene_lists/04_genes_unique_cogdx.txt"
)

path_resil <- safe_read(
  "biological_analysis/gene_lists/05_genes_pathology_resilience.txt"
)
tau_cog <- safe_read(
  "biological_analysis/gene_lists/06_genes_tau_cognitive.txt"
)
amy_cog <- safe_read(
  "biological_analysis/gene_lists/07_genes_amyloid_cognitive.txt"
)

intersection_all <- unique(c(path_resil, tau_cog, amy_cog))

# Inicializar lista para acumular rutas de todos los datasets
# para el reporte global
global_pathways_list <- list()

# Procesamiento por conjunto de datos (dataset)
for (ds in c("braaksc", "ceradsc", "cogdx")) {
  cat(paste("\n>>> PROCESANDO ENTORNO (LANDSCAPE):", toupper(ds), "\n"))

  # Cargar genes significativos de DEA para este dataset
  dea_file <- paste0("differential_expression/DEA_", ds, "_complete.csv")
  sig_dea_genes <- character(0)
  if (file.exists(dea_file)) {
    dea_data <- read_csv(dea_file, show_col_types = FALSE)
    is_sig <- dea_data$Significant == TRUE | dea_data$Significant == "TRUE"
    sig_dea_genes <- dea_data$Symbol[is_sig]
    sig_dea_genes <- sig_dea_genes[!is.na(sig_dea_genes)]
  }

  if (ds == "braaksc") {
    unique_list <- braak_unique
  } else if (ds == "ceradsc") {
    unique_list <- cerad_unique
  } else {
    unique_list <- cogdx_unique
  }

  # Buscar archivos de ORA y GSEA
  ora_files <- list.files(
    "results_enrichment/ORA/Significative",
    pattern = paste0("ora_.*_", ds, "\\.csv"),
    full.names = TRUE
  )
  ora_files <- ora_files[!grepl("GLOBAL", ora_files)]

  gsea_files <- list.files(
    "results_enrichment/GSEA/Significative",
    pattern = paste0("gsea_.*_", ds, "\\.csv"),
    full.names = TRUE
  )
  gsea_files <- gsea_files[!grepl("GLOBAL", gsea_files)]

  if (length(ora_files) == 0 && length(gsea_files) == 0) {
    cat("    - No se encontraron archivos de enriquecimiento. Omitiendo.\n")
    next
  }

  ora_df <- data.frame()
  if (length(ora_files) > 0) {
    ora_df <- bind_rows(
      lapply(ora_files, read_csv, show_col_types = FALSE)
    )
    if (nrow(ora_df) > 0) {
      ora_df <- ora_df |>
        dplyr::select(
          any_of(c("ID", "Description", "Category", "geneID"))
        ) |>
        dplyr::rename(Genes = .data$geneID) |>
        mutate(Source = "ORA")
    }
  }

  gsea_df <- data.frame()
  if (length(gsea_files) > 0) {
    gsea_df <- bind_rows(
      lapply(gsea_files, read_csv, show_col_types = FALSE)
    )
    if (nrow(gsea_df) > 0) {
      gsea_df <- gsea_df |>
        dplyr::select(
          any_of(c("ID", "Description", "Category", "core_enrichment"))
        ) |>
        dplyr::rename(Genes = .data$core_enrichment) |>
        mutate(Source = "GSEA")
    }
  }

  all_df <- bind_rows(ora_df, gsea_df)
  if (nrow(all_df) == 0) {
    cat("    - Los archivos de enriquecimiento están vacíos. Omitiendo.\n")
    next
  }

  # Unificar por ID
  cat("    - Unificando rutas...\n")
  unified_df <- all_df |>
    group_by(.data$ID, .data$Description, .data$Category) |>
    summarize(
      ORA_Presence = ifelse(any(.data$Source == "ORA"), 1, 0),
      GSEA_Presence = ifelse(any(.data$Source == "GSEA"), 1, 0),
      All_Genes = paste(
        unique(unlist(strsplit(Genes[!is.na(Genes)], "/"))),
        collapse = "/"
      ),
      .groups = "drop"
    )

  # Calcular conteos e intersecciones específicos
  cat("    - Calculando intersecciones de genes específicas...\n")
  unified_df <- unified_df |>
    rowwise() |>
    mutate(
      gene_list = list(unlist(strsplit(.data$All_Genes, "/"))),
      Triple_Core = length(intersect(.data$gene_list, triple_core)),
      Interseccion_Braaksc_Ceradsc = length(
        intersect(.data$gene_list, path_resil)
      ),
      Interseccion_Braaksc_Cogdx = length(
        intersect(.data$gene_list, tau_cog)
      ),
      Interseccion_Ceradsc_Cogdx = length(
        intersect(.data$gene_list, amy_cog)
      ),
      Unique_Braaksc = length(intersect(.data$gene_list, braak_unique)),
      Unique_Ceradsc = length(intersect(.data$gene_list, cerad_unique)),
      Unique_Cogdx = length(intersect(.data$gene_list, cogdx_unique)),
      Subnucleus = length(intersect(.data$gene_list, subnucleo)),
      Total_ML_Genes = .data$Triple_Core +
        .data$Interseccion_Braaksc_Ceradsc +
        .data$Interseccion_Braaksc_Cogdx +
        .data$Interseccion_Ceradsc_Cogdx +
        .data$Unique_Braaksc +
        .data$Unique_Ceradsc +
        .data$Unique_Cogdx +
        .data$Subnucleus,
      Triple_Core_Genes = paste(
        intersect(.data$gene_list, triple_core),
        collapse = "/"
      ),
      Interseccion_Braaksc_Ceradsc_Genes = paste(
        intersect(.data$gene_list, path_resil),
        collapse = "/"
      ),
      Interseccion_Braaksc_Cogdx_Genes = paste(
        intersect(.data$gene_list, tau_cog),
        collapse = "/"
      ),
      Interseccion_Ceradsc_Cogdx_Genes = paste(
        intersect(.data$gene_list, amy_cog),
        collapse = "/"
      ),
      Unique_Braaksc_Genes = paste(
        intersect(.data$gene_list, braak_unique),
        collapse = "/"
      ),
      Unique_Ceradsc_Genes = paste(
        intersect(.data$gene_list, cerad_unique),
        collapse = "/"
      ),
      Unique_Cogdx_Genes = paste(
        intersect(.data$gene_list, cogdx_unique),
        collapse = "/"
      ),
      Subnucleus_Genes = paste(
        intersect(.data$gene_list, subnucleo),
        collapse = "/"
      )
    ) |>
    ungroup() |>
    dplyr::select(-.data$gene_list) |>
    arrange(desc(.data$Total_ML_Genes))

  # Crear una versión con cabeceras claras para Excel en Español
  unified_csv_df <- unified_df |>
    dplyr::rename(
      Descripcion = .data$Description,
      Categoria = .data$Category,
      ORA = .data$ORA_Presence,
      GSEA = .data$GSEA_Presence,
      Todos_los_genes = .data$All_Genes,
      `Triple core` = .data$Triple_Core,
      `Interseccion braaksc-ceradsc` = .data$Interseccion_Braaksc_Ceradsc,
      `Interseccion braaksc-cogdx` = .data$Interseccion_Braaksc_Cogdx,
      `Interseccion ceradsc-cogdx` = .data$Interseccion_Ceradsc_Cogdx,
      `Unique braaksc` = .data$Unique_Braaksc,
      `Unique ceradsc` = .data$Unique_Ceradsc,
      `Unique cogdx` = .data$Unique_Cogdx,
      `Subnucleus` = .data$Subnucleus,
      `Total ML Genes` = .data$Total_ML_Genes,
      `Triple core genes` = .data$Triple_Core_Genes,
      `Interseccion braaksc-ceradsc genes` =
        .data$Interseccion_Braaksc_Ceradsc_Genes,
      `Interseccion braaksc-cogdx genes` =
        .data$Interseccion_Braaksc_Cogdx_Genes,
      `Interseccion ceradsc-cogdx genes` =
        .data$Interseccion_Ceradsc_Cogdx_Genes,
      `Unique braaksc genes` = .data$Unique_Braaksc_Genes,
      `Unique ceradsc genes` = .data$Unique_Ceradsc_Genes,
      `Unique cogdx genes` = .data$Unique_Cogdx_Genes,
      `Subnucleus genes` = .data$Subnucleus_Genes
    )

  # Exportar CSV
  write_excel_csv(
    unified_csv_df,
    file.path(base_dir_landscape, "Data", paste0("landscape_", ds, ".csv"))
  )
  cat(
    paste(
      "    - Datos de entorno guardados en",
      file.path(base_dir_landscape, "Data", paste0("landscape_", ds, ".csv")),
      "\n"
    )
  )

  # Acumular rutas para el reporte global por fragmento
  global_pathways_list[[ds]] <- unified_df |> mutate(Dataset = ds)

  # Mapa de calor (Heatmap)
  cat("    - Generando mapas de calor fragmentados...\n")
  out_dir <- file.path(
    base_dir_landscape, "Plots", ds, "functional_heatmap_global"
  )
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Seleccionar solo rutas con al menos 1 gen validado
  plot_df <- unified_df |> filter(.data$Total_ML_Genes > 0)

  if (nrow(plot_df) == 0) {
    cat("    - No se encontraron rutas con genes de ML para graficar.\n")
  } else {
    n_total <- nrow(plot_df)
    chunk_size <- 20
    n_chunks <- ceiling(n_total / chunk_size)

    cat(
      paste(
        "    - Total de rutas significativas con genes de ML:",
        n_total, "\n"
      )
    )
    cat(paste("    - Generando", n_chunks, "bloques de heatmaps...\n"))

    mat_cols <- c(
      "Triple_Core", "Interseccion_Braaksc_Ceradsc",
      "Interseccion_Braaksc_Cogdx", "Interseccion_Ceradsc_Cogdx",
      "Unique_Braaksc", "Unique_Ceradsc", "Unique_Cogdx", "Subnucleus"
    )

    for (i in seq_len(n_chunks)) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx   <- min(i * chunk_size, n_total)

      chunk_df <- plot_df[start_idx:end_idx, ]

      mat_data <- chunk_df |>
        dplyr::select(all_of(mat_cols)) |>
        as.matrix()

      # Renombrar columnas para que salgan perfectas en el Heatmap
      colnames(mat_data) <- c(
        "Triple core", "Interseccion braaksc-ceradsc",
        "Interseccion braaksc-cogdx", "Interseccion ceradsc-cogdx",
        "Unique braaksc", "Unique ceradsc", "Unique cogdx", "Subnucleus"
      )

      # Rownames = Description truncado
      labels <- str_trunc(chunk_df$Description, 50)
      rownames(mat_data) <- make.unique(labels)

      # Anotaciones
      anno_row <- data.frame(
        Database = chunk_df$Category,
        ORA = factor(
          ifelse(chunk_df$ORA_Presence == 1, "Presente", "Ausente"),
          levels = c("Ausente", "Presente")
        ),
        GSEA = factor(
          ifelse(chunk_df$GSEA_Presence == 1, "Presente", "Ausente"),
          levels = c("Ausente", "Presente")
        )
      )
      rownames(anno_row) <- rownames(mat_data)

      # Paleta de colores Premium
      anno_colors <- list(
        Database = c(
          GO_BP = "#2ECC71", GO_MF = "#3498DB", KEGG = "#E74C3C"
        ),
        ORA = c(Ausente = "#F2F4F4", Presente = "#F39C12"),
        GSEA = c(Ausente = "#F2F4F4", Presente = "#8E44AD")
      )

      filename <- file.path(
        out_dir,
        paste0("heatmap_landscape_chunk", i, ".png")
      )

      # Ajuste de tamaño
      h <- 4 + (nrow(mat_data) * 0.3)

      tryCatch({
        pheatmap(
          mat_data,
          color = colorRampPalette(
            c("#FFFFFF", "#EBF5FB", "#AED6F1", "#5DADE2", "#2E86C1")
          )(50),
          cluster_rows = (nrow(mat_data) >= 2),
          cluster_cols = FALSE,
          annotation_row = anno_row,
          annotation_colors = anno_colors,
          display_numbers = TRUE,
          number_format = "%.0f",
          fontsize_number = 10,
          number_color = "black",
          main = paste(
            "Functional Landscape:", toupper(ds), "(Part", i, ")"
          ),
          filename = filename,
          width = 11, height = h
        )
      }, error = function(e) {
        cat(
          sprintf(
            "      [ERROR] Falló generación del mapa de calor: %s\n",
            e$message
          )
        )
      })
    }
  }

  # 4.5 Generar Reportes CSV Específicos por Dataset
  cat("    - Exportando archivos CSV específicos del entorno...\n")
  tryCatch({
    # A. CSV: Rutas de GSEA que NO tienen genes significativos en DEA
    gsea_no_dea_df <- unified_df |>
      filter(.data$GSEA_Presence == 1) |>
      rowwise() |>
      mutate(
        gene_list = list(unlist(strsplit(.data$All_Genes, "/"))),
        n_DEA_Sig_Genes = length(intersect(.data$gene_list, sig_dea_genes))
      ) |>
      ungroup() |>
      filter(.data$n_DEA_Sig_Genes == 0) |>
      dplyr::select(
        .data$ID, .data$Description, .data$Category,
        .data$ORA_Presence, .data$GSEA_Presence, .data$All_Genes
      )

    gsea_no_dea_filename <- file.path(
      base_dir_landscape, "Data",
      paste0("gsea_no_dea_sig_pathways_", ds, ".csv")
    )
    gsea_no_dea_csv_df <- gsea_no_dea_df |>
      dplyr::rename(
        Descripcion = .data$Description,
        Categoria = .data$Category,
        ORA = .data$ORA_Presence,
        GSEA = .data$GSEA_Presence,
        Todos_los_genes = .data$All_Genes
      )
    write_excel_csv(gsea_no_dea_csv_df, gsea_no_dea_filename)

    # B. CSV: Rutas con genes del Triple Core
    triple_core_paths_df <- unified_df |>
      filter(.data$Triple_Core > 0) |>
      dplyr::select(
        .data$ID, .data$Description, .data$Category,
        .data$ORA_Presence, .data$GSEA_Presence, .data$Triple_Core,
        .data$Triple_Core_Genes, .data$Total_ML_Genes
      )

    triple_core_paths_filename <- file.path(
      base_dir_landscape, "Data",
      paste0("pathways_with_triple_core_", ds, ".csv")
    )
    triple_core_paths_csv_df <- triple_core_paths_df |>
      dplyr::rename(
        Descripcion = .data$Description,
        Categoria = .data$Category,
        ORA = .data$ORA_Presence,
        GSEA = .data$GSEA_Presence,
        `Triple core` = .data$Triple_Core,
        `Triple core genes` = .data$Triple_Core_Genes,
        `Total ML Genes` = .data$Total_ML_Genes
      )
    write_excel_csv(triple_core_paths_csv_df, triple_core_paths_filename)

    # C. CSV: Rutas con genes del Subnúcleo Irreducible
    subnucleus_paths_df <- unified_df |>
      filter(.data$Subnucleus > 0) |>
      dplyr::select(
        .data$ID, .data$Description, .data$Category,
        .data$ORA_Presence, .data$GSEA_Presence, .data$Subnucleus,
        .data$Subnucleus_Genes, .data$Total_ML_Genes
      )

    subnucleus_paths_filename <- file.path(
      base_dir_landscape, "Data",
      paste0("pathways_with_subnucleus_", ds, ".csv")
    )
    subnucleus_paths_csv_df <- subnucleus_paths_df |>
      dplyr::rename(
        Descripcion = .data$Description,
        Categoria = .data$Category,
        ORA = .data$ORA_Presence,
        GSEA = .data$GSEA_Presence,
        `Subnucleus` = .data$Subnucleus,
        `Subnucleus genes` = .data$Subnucleus_Genes,
        `Total ML Genes` = .data$Total_ML_Genes
      )
    write_excel_csv(subnucleus_paths_csv_df, subnucleus_paths_filename)

    cat("      [ÉXITO] Se exportaron 3 archivos CSV específicos.\n")
  }, error = function(e) {
    cat(
      sprintf(
        "      [ERROR] Falló exportación de CSV específicos: %s\n",
        e$message
      )
    )
  })

  # 4.6 Generar Reporte de Trazabilidad TXT por Dataset
  cat("    - Generando el reporte de trazabilidad biológica (.txt)...\n")
  tryCatch({
    report_filename <- file.path(
      base_dir_landscape, "Data", paste0("landscape_report_", ds, ".txt")
    )
    report_conn <- file(report_filename, "w")
    writeLines(
      "====================================================================",
      report_conn
    )
    writeLines(
      paste(
        paste(
          "BIOLOGICAL LANDSCAPE REPORT (TRAZABILIDAD DE",
          "MACHINE LEARNING Y RUTAS):"
        ),
        toupper(ds)
      ),
      report_conn
    )
    writeLines(
      "====================================================================",
      report_conn
    )
    writeLines(
      paste("Total unified pathways:", nrow(unified_df)),
      report_conn
    )
    writeLines(
      paste(
        "Pathways detected by ORA only:",
        sum(unified_df$ORA_Presence == 1 & unified_df$GSEA_Presence == 0)
      ),
      report_conn
    )
    writeLines(
      paste(
        "Pathways detected by GSEA only:",
        sum(unified_df$ORA_Presence == 0 & unified_df$GSEA_Presence == 1)
      ),
      report_conn
    )
    writeLines(
      paste(
        "Pathways detected by Both (ORA & GSEA):",
        sum(unified_df$ORA_Presence == 1 & unified_df$GSEA_Presence == 1)
      ),
      report_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      report_conn
    )
    writeLines("ML Gene Set Overlaps in Pathways:", report_conn)
    writeLines(
      paste(
        "- Pathways with Triple Core genes:",
        sum(unified_df$Triple_Core > 0)
      ),
      report_conn
    )
    writeLines(
      paste(
        "- Pathways with Subnucleus genes:",
        sum(unified_df$Subnucleus > 0)
      ),
      report_conn
    )

    if (ds == "braaksc") {
      writeLines(
        paste(
          "- Pathways with unique genes for Braaksc (Solo Tau):",
          sum(unified_df$Unique_Braaksc > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
          sum(unified_df$Interseccion_Braaksc_Ceradsc > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
          sum(unified_df$Interseccion_Braaksc_Cogdx > 0)
        ),
        report_conn
      )
    } else if (ds == "ceradsc") {
      writeLines(
        paste(
          "- Pathways with unique genes for Ceradsc (Solo Amiloide):",
          sum(unified_df$Unique_Ceradsc > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
          sum(unified_df$Interseccion_Braaksc_Ceradsc > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
          sum(unified_df$Interseccion_Ceradsc_Cogdx > 0)
        ),
        report_conn
      )
    } else if (ds == "cogdx") {
      writeLines(
        paste(
          "- Pathways with unique genes for Cogdx (Clínica):",
          sum(unified_df$Unique_Cogdx > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
          sum(unified_df$Interseccion_Braaksc_Cogdx > 0)
        ),
        report_conn
      )
      writeLines(
        paste(
          "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
          sum(unified_df$Interseccion_Ceradsc_Cogdx > 0)
        ),
        report_conn
      )
    }

    # PÁRRAFO GSEA SIN GENES DEA SIGNIFICATIVOS
    gsea_no_dea_count <- if (exists("gsea_no_dea_df")) {
      nrow(gsea_no_dea_df)
    } else {
      0
    }
    writeLines(
      "--------------------------------------------------------------------",
      report_conn
    )
    writeLines(
      "ANÁLISIS DE RUTAS GSEA EXCLUSIVAS (SIN GENES SIGNIFICATIVOS EN DEA):",
      report_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      report_conn
    )
    writeLines(
      paste(
        "El análisis tradicional de sobre-representación (ORA) requiere",
        "que los genes pasen un umbral estricto de significación en DEA.",
        "Sin embargo, GSEA funciona ordenando todos los genes del genoma,",
        "lo que permite detectar cambios sutiles y coordinados incluso",
        "si ningún gen por sí solo es estadísticamente significativo en DEA."
      ),
      report_conn
    )
    writeLines(
      paste(
        "Para el dataset", toupper(ds), ", se han identificado un total",
        "de", gsea_no_dea_count, "rutas enriquecidas por GSEA que",
        "contienen exactamente CERO genes significativos en DEA dentro",
        "de su 'core enrichment'."
      ),
      report_conn
    )
    writeLines(
      paste(
        "Este hallazgo metodológico demuestra la tremenda capacidad de GSEA",
        "para revelar mecanismos biológicos coordinados (por ejemplo,",
        "cascadas metabólicas completas con pequeñas desregulaciones en",
        "cada paso) que habrían permanecido completamente invisibles si",
        "solo se hubieran estudiado los genes significativos del DEA."
      ),
      report_conn
    )
    writeLines(
      paste(
        "El desglose de estas rutas se ha guardado en el archivo CSV:",
        paste0("gsea_no_dea_sig_pathways_", ds, ".csv")
      ),
      report_conn
    )

    if (gsea_no_dea_count > 0) {
      writeLines(
        "\nTop 10 rutas de GSEA sin genes significativos en DEA:",
        report_conn
      )
      top_no_dea <- head(gsea_no_dea_df, 10)
      for (k in seq_len(nrow(top_no_dea))) {
        writeLines(
          paste0(
            "  - ", top_no_dea$Description[k], " (",
            top_no_dea$ID[k], ") [", top_no_dea$Category[k], "]"
          ),
          report_conn
        )
      }
    }

    # RESUMEN DETALLADO DATABASE POR DATABASE
    writeLines(
      "--------------------------------------------------------------------",
      report_conn
    )
    writeLines(
      paste(
        "RESUMEN DETALLADO RUTA POR RUTA Y DATABASE POR",
        "DATABASE (CON GENES DE ML):"
      ),
      report_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      report_conn
    )

    categories <- c("GO_BP", "GO_MF", "KEGG")
    for (cat_name in categories) {
      cat_df <- unified_df |> filter(.data$Category == cat_name)
      writeLines(
        paste(
          "\n=== DATABASE:", cat_name,
          paste0("(", nrow(cat_df), " rutas unificadas) ===")
        ),
        report_conn
      )
      if (nrow(cat_df) == 0) {
        writeLines(
          "  No se encontraron rutas significativas en esta base de datos.",
          report_conn
        )
        next
      }

      for (r in seq_len(nrow(cat_df))) {
        pathway_line <- paste0(
          "  * ", cat_df$Description[r], " (", cat_df$ID[r], ")"
        )
        if (cat_df$Total_ML_Genes[r] > 0) {
          genes_str <- c()
          if (cat_df$Triple_Core[r] > 0) {
            genes_str <- c(
              genes_str,
              paste0(
                "Triple Core [", cat_df$Triple_Core[r], "]: ",
                cat_df$Triple_Core_Genes[r]
              )
            )
          }
          if (cat_df$Subnucleus[r] > 0) {
            genes_str <- c(
              genes_str,
              paste0(
                "Subnúcleo [", cat_df$Subnucleus[r], "]: ",
                cat_df$Subnucleus_Genes[r]
              )
            )
          }
          if (cat_df$Unique_Braaksc[r] > 0) {
            genes_str <- c(
              genes_str,
              paste0(
                "Unique Braak [", cat_df$Unique_Braaksc[r], "]: ",
                cat_df$Unique_Braaksc_Genes[r]
              )
            )
          }
          if (cat_df$Unique_Ceradsc[r] > 0) {
            genes_str <- c(
              genes_str,
              paste0(
                "Unique Cerad [", cat_df$Unique_Ceradsc[r], "]: ",
                cat_df$Unique_Ceradsc_Genes[r]
              )
            )
          }
          if (cat_df$Unique_Cogdx[r] > 0) {
            genes_str <- c(
              genes_str,
              paste0(
                "Unique Cogdx [", cat_df$Unique_Cogdx[r], "]: ",
                cat_df$Unique_Cogdx_Genes[r]
              )
            )
          }

          writeLines(
            paste0(
              pathway_line, " | Genes ML: ",
              paste(genes_str, collapse = " ; ")
            ),
            report_conn
          )
        } else {
          writeLines(
            paste0(pathway_line, " | Sin genes del conjunto de ML"),
            report_conn
          )
        }
      }
    }
    close(report_conn)
    cat(
      paste(
        "      [ÉXITO] Reporte de trazabilidad guardado en:",
        report_filename, "\n"
      )
    )
  }, error = function(e) {
    cat(
      sprintf(
        "      [ERROR] Falló el reporte de trazabilidad biológica: %s\n",
        e$message
      )
    )
  })
}

# --------------------------------------------------------------------
# GENERACIÓN DE REPORTES DE RUTAS SIGNIFICATIVAS GLOBALES
# --------------------------------------------------------------------
cat("\n>>> GENERANDO REPORTES BIOLÓGICOS GLOBALES DE RUTAS...\n")
tryCatch({
  global_pathways_all <- bind_rows(global_pathways_list)

  if (nrow(global_pathways_all) > 0) {
    # Agrupar las rutas comunes que aparecen en múltiples datasets
    global_pathways_grouped <- global_pathways_all |>
      group_by(.data$ID, .data$Description, .data$Category) |>
      summarize(
        Datasets_Detected = paste(
          unique(toupper(.data$Dataset)),
          collapse = ", "
        ),
        Triple_Core_Count = max(.data$Triple_Core, na.rm = TRUE),
        Subnucleus_Count = max(.data$Subnucleus, na.rm = TRUE),
        Unique_Braaksc_Count = max(.data$Unique_Braaksc, na.rm = TRUE),
        Unique_Ceradsc_Count = max(.data$Unique_Ceradsc, na.rm = TRUE),
        Unique_Cogdx_Count = max(.data$Unique_Cogdx, na.rm = TRUE),
        Int_Braak_Cerad_Count = max(
          .data$Interseccion_Braaksc_Ceradsc,
          na.rm = TRUE
        ),
        Int_Braak_Cogdx_Count = max(
          .data$Interseccion_Braaksc_Cogdx,
          na.rm = TRUE
        ),
        Int_Cerad_Cogdx_Count = max(
          .data$Interseccion_Ceradsc_Cogdx,
          na.rm = TRUE
        ),
        Triple_Core_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Triple_Core_Genes[
                  !is.na(Triple_Core_Genes) & Triple_Core_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Subnucleus_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Subnucleus_Genes[
                  !is.na(Subnucleus_Genes) & Subnucleus_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Unique_Braak_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Unique_Braaksc_Genes[
                  !is.na(Unique_Braaksc_Genes) & Unique_Braaksc_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Unique_Cerad_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Unique_Ceradsc_Genes[
                  !is.na(Unique_Ceradsc_Genes) & Unique_Ceradsc_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Unique_Cogdx_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Unique_Cogdx_Genes[
                  !is.na(Unique_Cogdx_Genes) & Unique_Cogdx_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Int_Braak_Cerad_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Interseccion_Braaksc_Ceradsc_Genes[
                  !is.na(Interseccion_Braaksc_Ceradsc_Genes) &
                    Interseccion_Braaksc_Ceradsc_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Int_Braak_Cogdx_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Interseccion_Braaksc_Cogdx_Genes[
                  !is.na(Interseccion_Braaksc_Cogdx_Genes) &
                    Interseccion_Braaksc_Cogdx_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        Int_Cerad_Cogdx_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                Interseccion_Ceradsc_Cogdx_Genes[
                  !is.na(Interseccion_Ceradsc_Cogdx_Genes) &
                    Interseccion_Ceradsc_Cogdx_Genes != ""
                ],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        # Consolidador de TODOS los genes de ML en esta ruta
        All_ML_Genes_List = paste(
          unique(
            unlist(
              strsplit(
                All_Genes[!is.na(All_Genes) & All_Genes != ""],
                "/"
              )
            )
          ),
          collapse = "/"
        ),
        .groups = "drop"
      ) |>
      rowwise() |>
      mutate(
        Total_ML_Genes_Count = length(
          unique(unlist(strsplit(.data$All_ML_Genes_List, "/")))
        )
      ) |>
      ungroup()

    # Calcular estadísticas globales unificadas para las cabeceras
    total_global_pathways <- nrow(global_pathways_grouped)

    pathway_methods <- global_pathways_all |>
      group_by(.data$ID) |>
      summarize(
        ORA_Glob = any(.data$ORA_Presence == 1),
        GSEA_Glob = any(.data$GSEA_Presence == 1)
      )

    ora_only_global <- sum(
      pathway_methods$ORA_Glob & !pathway_methods$GSEA_Glob
    )
    gsea_only_global <- sum(
      pathway_methods$GSEA_Glob & !pathway_methods$ORA_Glob
    )
    both_global <- sum(
      pathway_methods$ORA_Glob & pathway_methods$GSEA_Glob
    )

    tc_global_count <- sum(global_pathways_grouped$Triple_Core_Count > 0)
    sub_global_count <- sum(global_pathways_grouped$Subnucleus_Count > 0)
    uniq_br_global_count <- sum(
      global_pathways_grouped$Unique_Braaksc_Count > 0
    )
    uniq_ce_global_count <- sum(
      global_pathways_grouped$Unique_Ceradsc_Count > 0
    )
    uniq_co_global_count <- sum(
      global_pathways_grouped$Unique_Cogdx_Count > 0
    )
    int_bc_global_count <- sum(
      global_pathways_grouped$Int_Braak_Cerad_Count > 0
    )
    int_bg_global_count <- sum(
      global_pathways_grouped$Int_Braak_Cogdx_Count > 0
    )
    int_cg_global_count <- sum(
      global_pathways_grouped$Int_Cerad_Cogdx_Count > 0
    )

    # ==========================================================================
    # ARCHIVO 1: REPORTE GLOBAL POR FRAGMENTOS (con cabecera unificada)
    # ==========================================================================
    global_pathways_filename <- file.path(
      base_dir_landscape, "Data",
      "landscape_global_pathways_by_biological_fragment.txt"
    )
    global_pathways_conn <- file(global_pathways_filename, "w")

    writeLines(
      "====================================================================",
      global_pathways_conn
    )
    writeLines(
      "REPORTE GLOBAL DE RUTAS SIGNIFICATIVAS POR FRAGMENTO BIOLÓGICO DE ML",
      global_pathways_conn
    )
    writeLines(
      "====================================================================",
      global_pathways_conn
    )
    writeLines(
      paste("Total unified pathways:", total_global_pathways),
      global_pathways_conn
    )
    writeLines(
      paste("Pathways detected by ORA only:", ora_only_global),
      global_pathways_conn
    )
    writeLines(
      paste("Pathways detected by GSEA only:", gsea_only_global),
      global_pathways_conn
    )
    writeLines(
      paste("Pathways detected by Both (ORA & GSEA):", both_global),
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_pathways_conn
    )
    writeLines(
      "ML Gene Set Overlaps in Pathways (Global):",
      global_pathways_conn
    )
    writeLines(
      paste("- Pathways with Triple Core genes:", tc_global_count),
      global_pathways_conn
    )
    writeLines(
      paste("- Pathways with Subnucleus genes:", sub_global_count),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Braaksc (Solo Tau):",
        uniq_br_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Ceradsc (Solo Amiloide):",
        uniq_ce_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Cogdx (Clínica):",
        uniq_co_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
        int_bc_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
        int_bg_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
        int_cg_global_count
      ),
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------\n",
      global_pathways_conn
    )

    # 1. TRIPLE CORE
    tc_paths <- global_pathways_grouped |>
      filter(.data$Triple_Core_Count > 0) |>
      arrange(desc(.data$Triple_Core_Count))
    writeLines(
      paste(
        "1. FRAGMENTO TRIPLE CORE (Rutas biológicas de herencia sistémica",
        "afectadas por componentes del Triple Core de ML)"
      ),
      global_pathways_conn
    )
    writeLines(
      paste(
        "   Total de rutas que contienen genes del Triple Core:",
        nrow(tc_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(tc_paths) > 0) {
      for (r in seq_len(nrow(tc_paths))) {
        writeLines(
          paste0(
            "  * ", tc_paths$Description[r], " (", tc_paths$ID[r], ") [",
            tc_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "    Datasets: ", tc_paths$Datasets_Detected[r],
            " | Genes [", tc_paths$Triple_Core_Count[r], "]: ",
            tc_paths$Triple_Core_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "  No se encontraron rutas que contengan genes del Triple Core.",
        global_pathways_conn
      )
    }
    writeLines(
      "\n--------------------------------------------------------------------",
      global_pathways_conn
    )

    # 2. SUBNUCLEUS
    sub_paths <- global_pathways_grouped |>
      filter(.data$Subnucleus_Count > 0) |>
      arrange(desc(.data$Subnucleus_Count))
    writeLines(
      "2. FRAGMENTO SUBNÚCLEO IRREDUCIBLE (16 genes biomarcadores clave)",
      global_pathways_conn
    )
    writeLines(
      paste(
        "   Total de rutas que contienen genes del Subnúcleo Irreducible:",
        nrow(sub_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(sub_paths) > 0) {
      for (r in seq_len(nrow(sub_paths))) {
        writeLines(
          paste0(
            "  * ", sub_paths$Description[r], " (", sub_paths$ID[r], ") [",
            sub_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "    Datasets: ", sub_paths$Datasets_Detected[r],
            " | Genes [", sub_paths$Subnucleus_Count[r], "]: ",
            sub_paths$Subnucleus_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "  No se encontraron rutas que contengan genes del Subnúcleo.",
        global_pathways_conn
      )
    }
    writeLines(
      "\n--------------------------------------------------------------------",
      global_pathways_conn
    )

    # 3. INTERSECCIONES
    writeLines(
      "3. FRAGMENTOS DE INTERSECCIÓN ENTRE DOS DATASETS (2 vs 1)",
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_pathways_conn
    )

    # A. Braaksc & Ceradsc (Resiliencia)
    res_paths <- global_pathways_grouped |>
      filter(.data$Int_Braak_Cerad_Count > 0) |>
      arrange(desc(.data$Int_Braak_Cerad_Count))
    writeLines(
      "   A. Intersección Braaksc y Ceradsc (Resiliencia Cognitiva)",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes de Resiliencia:",
        nrow(res_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(res_paths) > 0) {
      for (r in seq_len(nrow(res_paths))) {
        writeLines(
          paste0(
            "    * ", res_paths$Description[r], " (", res_paths$ID[r], ") [",
            res_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", res_paths$Datasets_Detected[r],
            " | Genes [", res_paths$Int_Braak_Cerad_Count[r], "]: ",
            res_paths$Int_Braak_Cerad_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes de esta intersección.",
        global_pathways_conn
      )
    }
    writeLines("", global_pathways_conn)

    # B. Braaksc & Cogdx (Tau-Cogdx)
    tau_paths <- global_pathways_grouped |>
      filter(.data$Int_Braak_Cogdx_Count > 0) |>
      arrange(desc(.data$Int_Braak_Cogdx_Count))
    writeLines(
      "   B. Intersección Braaksc y Cogdx (Disociación Tau-Cogdx)",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes de Tau-Cogdx:",
        nrow(tau_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(tau_paths) > 0) {
      for (r in seq_len(nrow(tau_paths))) {
        writeLines(
          paste0(
            "    * ", tau_paths$Description[r], " (", tau_paths$ID[r], ") [",
            tau_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", tau_paths$Datasets_Detected[r],
            " | Genes [", tau_paths$Int_Braak_Cogdx_Count[r], "]: ",
            tau_paths$Int_Braak_Cogdx_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes de esta intersección.",
        global_pathways_conn
      )
    }
    writeLines("", global_pathways_conn)

    # C. Ceradsc & Cogdx (Amiloide-Cogdx)
    amy_paths <- global_pathways_grouped |>
      filter(.data$Int_Cerad_Cogdx_Count > 0) |>
      arrange(desc(.data$Int_Cerad_Cogdx_Count))
    writeLines(
      "   C. Intersección Ceradsc y Cogdx (Disociación Amiloide-Cogdx)",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes de Amiloide-Cogdx:",
        nrow(amy_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(amy_paths) > 0) {
      for (r in seq_len(nrow(amy_paths))) {
        writeLines(
          paste0(
            "    * ", amy_paths$Description[r], " (", amy_paths$ID[r], ") [",
            amy_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", amy_paths$Datasets_Detected[r],
            " | Genes [", amy_paths$Int_Cerad_Cogdx_Count[r], "]: ",
            amy_paths$Int_Cerad_Cogdx_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes de esta intersección.",
        global_pathways_conn
      )
    }
    writeLines(
      "\n--------------------------------------------------------------------",
      global_pathways_conn
    )

    # 4. FRAGMENTOS EXCLUSIVOS
    writeLines(
      "4. FRAGMENTOS DE GENES EXCLUSIVOS (Unique por Dataset)",
      global_pathways_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_pathways_conn
    )

    # A. Braaksc Unique
    uniq_br_paths <- global_pathways_grouped |>
      filter(.data$Unique_Braaksc_Count > 0) |>
      arrange(desc(.data$Unique_Braaksc_Count))
    writeLines(
      "   A. Exclusivos de Braaksc (Solo Tau):",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes exclusivos de Braaksc:",
        nrow(uniq_br_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(uniq_br_paths) > 0) {
      for (r in seq_len(nrow(uniq_br_paths))) {
        writeLines(
          paste0(
            "    * ", uniq_br_paths$Description[r], " (",
            uniq_br_paths$ID[r], ") [", uniq_br_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", uniq_br_paths$Datasets_Detected[r],
            " | Genes [", uniq_br_paths$Unique_Braaksc_Count[r], "]: ",
            uniq_br_paths$Unique_Braak_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes exclusivos de Braaksc.",
        global_pathways_conn
      )
    }
    writeLines("", global_pathways_conn)

    # B. Ceradsc Unique
    uniq_ce_paths <- global_pathways_grouped |>
      filter(.data$Unique_Ceradsc_Count > 0) |>
      arrange(desc(.data$Unique_Ceradsc_Count))
    writeLines(
      "   B. Exclusivos de Ceradsc (Solo Amiloide):",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes exclusivos de Ceradsc:",
        nrow(uniq_ce_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(uniq_ce_paths) > 0) {
      for (r in seq_len(nrow(uniq_ce_paths))) {
        writeLines(
          paste0(
            "    * ", uniq_ce_paths$Description[r], " (",
            uniq_ce_paths$ID[r], ") [", uniq_ce_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", uniq_ce_paths$Datasets_Detected[r],
            " | Genes [", uniq_ce_paths$Unique_Ceradsc_Count[r], "]: ",
            uniq_ce_paths$Unique_Cerad_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes exclusivos de Ceradsc.",
        global_pathways_conn
      )
    }
    writeLines("", global_pathways_conn)

    # C. Cogdx Unique
    uniq_co_paths <- global_pathways_grouped |>
      filter(.data$Unique_Cogdx_Count > 0) |>
      arrange(desc(.data$Unique_Cogdx_Count))
    writeLines(
      "   C. Exclusivos de Cogdx (Clínica):",
      global_pathways_conn
    )
    writeLines(
      paste(
        "      Total de rutas que contienen genes exclusivos de Cogdx:",
        nrow(uniq_co_paths)
      ),
      global_pathways_conn
    )
    writeLines(
      "      ------------------------------------------------------------",
      global_pathways_conn
    )
    if (nrow(uniq_co_paths) > 0) {
      for (r in seq_len(nrow(uniq_co_paths))) {
        writeLines(
          paste0(
            "    * ", uniq_co_paths$Description[r], " (",
            uniq_co_paths$ID[r], ") [", uniq_co_paths$Category[r], "]"
          ),
          global_pathways_conn
        )
        writeLines(
          paste0(
            "      Datasets: ", uniq_co_paths$Datasets_Detected[r],
            " | Genes [", uniq_co_paths$Unique_Cogdx_Count[r], "]: ",
            uniq_co_paths$Unique_Cogdx_Genes_List[r]
          ),
          global_pathways_conn
        )
      }
    } else {
      writeLines(
        "      No se encontraron rutas con genes exclusivos de Cogdx.",
        global_pathways_conn
      )
    }

    close(global_pathways_conn)
    cat(
      paste(
        "      [ÉXITO] Reporte global de rutas guardado:",
        global_pathways_filename, "\n"
      )
    )

    # ==========================================================================
    # ARCHIVO 2: REPORTE GLOBAL INTEGRADO EN LISTA ÚNICA
    # ==========================================================================
    global_unique_filename <- file.path(
      base_dir_landscape, "Data",
      "landscape_global_pathways_unique_list.txt"
    )
    global_unique_conn <- file(global_unique_filename, "w")

    writeLines(
      "====================================================================",
      global_unique_conn
    )
    writeLines(
      "REPORTE GLOBAL DE RUTAS SIGNIFICATIVAS UNIFICADAS DE ML (LISTA ÚNICA)",
      global_unique_conn
    )
    writeLines(
      "====================================================================",
      global_unique_conn
    )
    writeLines(
      paste("Total unified pathways:", total_global_pathways),
      global_unique_conn
    )
    writeLines(
      paste("Pathways detected by ORA only:", ora_only_global),
      global_unique_conn
    )
    writeLines(
      paste("Pathways detected by GSEA only:", gsea_only_global),
      global_unique_conn
    )
    writeLines(
      paste("Pathways detected by Both (ORA & GSEA):", both_global),
      global_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_unique_conn
    )
    writeLines(
      "ML Gene Set Overlaps in Pathways (Global):",
      global_unique_conn
    )
    writeLines(
      paste("- Pathways with Triple Core genes:", tc_global_count),
      global_unique_conn
    )
    writeLines(
      paste("- Pathways with Subnucleus genes:", sub_global_count),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Braaksc (Solo Tau):",
        uniq_br_global_count
      ),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Ceradsc (Solo Amiloide):",
        uniq_ce_global_count
      ),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Cogdx (Clínica):",
        uniq_co_global_count
      ),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
        int_bc_global_count
      ),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
        int_bg_global_count
      ),
      global_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
        int_cg_global_count
      ),
      global_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      global_unique_conn
    )
    writeLines(
      paste(
        "Este reporte consolida de forma unificada todas las rutas",
        "biológicas significativas encontradas para Braaksc, Ceradsc y",
        "Cogdx. Cada ruta se muestra una sola vez con la lista total",
        "integrada de todos sus genes solapantes del conjunto de ML."
      ),
      global_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------\n",
      global_unique_conn
    )

    # Ordenar las rutas por número de genes de ML solapantes desc
    global_pathways_sorted <- global_pathways_grouped |>
      arrange(desc(.data$Total_ML_Genes_Count), .data$Description)

    for (r in seq_len(nrow(global_pathways_sorted))) {
      writeLines(
        paste0(
          "  * ", global_pathways_sorted$Description[r], " (",
          global_pathways_sorted$ID[r], ") [",
          global_pathways_sorted$Category[r], "]"
        ),
        global_unique_conn
      )
      writeLines(
        paste0(
          "    Datasets: ", global_pathways_sorted$Datasets_Detected[r],
          " | Total ML Genes [",
          global_pathways_sorted$Total_ML_Genes_Count[r], "]: ",
          global_pathways_sorted$All_ML_Genes_List[r]
        ),
        global_unique_conn
      )
    }

    close(global_unique_conn)
    cat(
      paste(
        "      [ÉXITO] Reporte global de rutas unificadas guardado:",
        global_unique_filename, "\n"
      )
    )

    # ==========================================================================
    # GENERACIÓN DE REPORTES INVERSOS: GENES -> RUTAS
    # ==========================================================================
    cat("\n>>> GENERANDO REPORTES BIOLÓGICOS GLOBALES DE GENES A RUTAS...\n")

    # Expandir rutas a nivel de gen
    genes_to_pathways <- global_pathways_all |>
      dplyr::select(
        .data$ID, .data$Description, .data$Category,
        .data$Dataset, .data$All_Genes
      ) |>
      separate_rows(.data$All_Genes, sep = "/") |>
      dplyr::rename(Gene = .data$All_Genes) |>
      filter(!is.na(.data$Gene) & .data$Gene != "")

    # Agrupar por Gen y Ruta
    gene_pathway_grouped <- genes_to_pathways |>
      group_by(.data$Gene, .data$ID, .data$Description, .data$Category) |>
      summarize(
        Datasets_Detected = paste(
          unique(toupper(.data$Dataset)),
          collapse = ", "
        ),
        .groups = "drop"
      )

    # Calcular el total de rutas por gen
    gene_pathway_counts <- gene_pathway_grouped |>
      group_by(.data$Gene) |>
      summarize(Total_Pathways = n(), .groups = "drop")

    # Unir conteos y clasificar biológicamente cada gen
    gene_pathway_merged <- gene_pathway_grouped |>
      left_join(gene_pathway_counts, by = "Gene") |>
      rowwise() |>
      mutate(
        Fragment = case_when(
          .data$Gene %in% subnucleo ~ "SUBNÚCLEO IRREDUCIBLE",
          .data$Gene %in% triple_core ~ "TRIPLE CORE",
          .data$Gene %in% path_resil ~
            "BRAAKSC-CERADSC INTERSECTION (RESILIENCIA COGNITIVA)",
          .data$Gene %in% tau_cog ~
            "BRAAKSC-COGDX INTERSECTION (TAU-COGNITIVO)",
          .data$Gene %in% amy_cog ~
            "CERADSC-COGDX INTERSECTION (AMILOIDE-COGNITIVO)",
          .data$Gene %in% braak_unique ~ "EXCLUSIVO BRAAKSC (SOLO TAU)",
          .data$Gene %in% cerad_unique ~
            "EXCLUSIVO CERADSC (SOLO AMILOIDE)",
          .data$Gene %in% cogdx_unique ~ "EXCLUSIVO COGDX (CLÍNICA)",
          TRUE ~ "OTRO GEN DE ML"
        )
      ) |>
      ungroup()

    # Obtener la lista única de genes con sus conteos y fragmento
    unique_genes_list <- gene_pathway_merged |>
      dplyr::select(.data$Gene, .data$Total_Pathways, .data$Fragment) |>
      distinct() |>
      arrange(desc(.data$Total_Pathways), .data$Gene)

    # ==========================================================================
    # ARCHIVO 3: GENES -> RUTAS POR FRAGMENTO BIOLÓGICO
    # ==========================================================================
    genes_fragment_filename <- file.path(
      base_dir_landscape, "Data",
      "landscape_global_genes_by_biological_fragment.txt"
    )
    genes_fragment_conn <- file(genes_fragment_filename, "w")

    writeLines(
      "====================================================================",
      genes_fragment_conn
    )
    writeLines(
      "REPORTE GLOBAL DE GENES DE ML Y SUS RUTAS ASOCIADAS POR FRAGMENTO",
      genes_fragment_conn
    )
    writeLines(
      "====================================================================",
      genes_fragment_conn
    )
    writeLines(
      paste("Total unified pathways:", total_global_pathways),
      genes_fragment_conn
    )
    writeLines(
      paste("Pathways detected by ORA only:", ora_only_global),
      genes_fragment_conn
    )
    writeLines(
      paste("Pathways detected by GSEA only:", gsea_only_global),
      genes_fragment_conn
    )
    writeLines(
      paste("Pathways detected by Both (ORA & GSEA):", both_global),
      genes_fragment_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_fragment_conn
    )
    writeLines(
      "ML Gene Set Overlaps in Pathways (Global):",
      genes_fragment_conn
    )
    writeLines(
      paste("- Pathways with Triple Core genes:", tc_global_count),
      genes_fragment_conn
    )
    writeLines(
      paste("- Pathways with Subnucleus genes:", sub_global_count),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Braaksc (Solo Tau):",
        uniq_br_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Ceradsc (Solo Amiloide):",
        uniq_ce_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Cogdx (Clínica):",
        uniq_co_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
        int_bc_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
        int_bg_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      paste(
        "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
        int_cg_global_count
      ),
      genes_fragment_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_fragment_conn
    )
    writeLines(
      paste(
        "Este reporte organiza todos los genes del conjunto de Machine",
        "Learning que tienen solapamiento con al menos una ruta biológica",
        "significativa, clasificados por su correspondiente",
        "fragmento biológico."
      ),
      genes_fragment_conn
    )
    writeLines(
      "--------------------------------------------------------------------\n",
      genes_fragment_conn
    )

    write_gene_section <- function(genes_sub, section_title, conn) {
      writeLines(section_title, conn)
      writeLines(
        paste(rep("-", nchar(section_title)), collapse = ""),
        conn
      )
      if (nrow(genes_sub) > 0) {
        for (g in seq_len(nrow(genes_sub))) {
          gene_name <- genes_sub$Gene[g]
          writeLines(paste0("  * GEN: ", gene_name), conn)
          writeLines(
            paste0(
              "    Total de rutas asociadas: ",
              genes_sub$Total_Pathways[g]
            ),
            conn
          )
          writeLines("    Rutas:", conn)

          # Filtrar rutas de este gen
          gene_paths <- gene_pathway_merged |>
            filter(.data$Gene == gene_name) |>
            arrange(.data$Description)

          for (p in seq_len(nrow(gene_paths))) {
            writeLines(
              paste0(
                "      - ", gene_paths$Description[p], " (",
                gene_paths$ID[p], ") [", gene_paths$Category[p],
                "] | Datasets: ", gene_paths$Datasets_Detected[p]
              ),
              conn
            )
          }
          writeLines("", conn)
        }
      } else {
        writeLines(
          "  No hay genes con rutas significativas en este fragmento.",
          conn
        )
        writeLines("", conn)
      }
      writeLines(
        "------------------------------------------------------------\n",
        conn
      )
    }

    # A. Subnúcleo
    sub_genes <- unique_genes_list |>
      filter(.data$Fragment == "SUBNÚCLEO IRREDUCIBLE")
    write_gene_section(
      sub_genes,
      "1. FRAGMENTO SUBNÚCLEO IRREDUCIBLE (Biomarcadores Clave)",
      genes_fragment_conn
    )

    # B. Triple Core
    tc_genes <- unique_genes_list |>
      filter(.data$Fragment == "TRIPLE CORE")
    write_gene_section(
      tc_genes,
      "2. FRAGMENTO TRIPLE CORE (Genes comunes en 3 Datasets)",
      genes_fragment_conn
    )

    # C. Intersecciones dobles
    writeLines(
      "3. FRAGMENTOS DE INTERSECCIÓN ENTRE DOS DATASETS (Disociaciones)",
      genes_fragment_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_fragment_conn
    )

    bc_genes <- unique_genes_list |>
      filter(
        .data$Fragment ==
          "BRAAKSC-CERADSC INTERSECTION (RESILIENCIA COGNITIVA)"
      )
    write_gene_section(
      bc_genes,
      "   A. Intersección Braaksc y Ceradsc (Resiliencia Cognitiva)",
      genes_fragment_conn
    )

    bg_genes <- unique_genes_list |>
      filter(
        .data$Fragment == "BRAAKSC-COGDX INTERSECTION (TAU-COGNITIVO)"
      )
    write_gene_section(
      bg_genes,
      "   B. Intersección Braaksc y Cogdx (Tau-Cogdx)",
      genes_fragment_conn
    )

    cg_genes <- unique_genes_list |>
      filter(
        .data$Fragment == "CERADSC-COGDX INTERSECTION (AMILOIDE-COGNITIVO)"
      )
    write_gene_section(
      cg_genes,
      "   C. Intersección Ceradsc y Cogdx (Amiloide-Cogdx)",
      genes_fragment_conn
    )

    # D. Exclusivos
    writeLines(
      "4. FRAGMENTOS DE GENES EXCLUSIVOS (Unique por Dataset)",
      genes_fragment_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_fragment_conn
    )

    uniq_br_g <- unique_genes_list |>
      filter(.data$Fragment == "EXCLUSIVO BRAAKSC (SOLO TAU)")
    write_gene_section(
      uniq_br_g,
      "   A. Exclusivos de Braaksc (Solo Tau)",
      genes_fragment_conn
    )

    uniq_ce_g <- unique_genes_list |>
      filter(.data$Fragment == "EXCLUSIVO CERADSC (SOLO AMILOIDE)")
    write_gene_section(
      uniq_ce_g,
      "   B. Exclusivos de Ceradsc (Solo Amiloide)",
      genes_fragment_conn
    )

    uniq_co_g <- unique_genes_list |>
      filter(.data$Fragment == "EXCLUSIVO COGDX (CLÍNICA)")
    write_gene_section(
      uniq_co_g,
      "   C. Exclusivos de Cogdx (Clínica)",
      genes_fragment_conn
    )

    close(genes_fragment_conn)
    cat(
      paste(
        "      [ÉXITO] Reporte de genes por fragmento guardado:",
        genes_fragment_filename, "\n"
      )
    )

    # ==========================================================================
    # ARCHIVO 4: GENES -> RUTAS EN LISTA ÚNICA INTEGRADA
    # ==========================================================================
    genes_unique_filename <- file.path(
      base_dir_landscape, "Data",
      "landscape_global_genes_unique_list.txt"
    )
    genes_unique_conn <- file(genes_unique_filename, "w")

    writeLines(
      "====================================================================",
      genes_unique_conn
    )
    writeLines(
      "REPORTE GLOBAL DE GENES DE ML Y SUS RUTAS ASOCIADAS (LISTA ÚNICA)",
      genes_unique_conn
    )
    writeLines(
      "====================================================================",
      genes_unique_conn
    )
    writeLines(
      paste("Total unified pathways:", total_global_pathways),
      genes_unique_conn
    )
    writeLines(
      paste("Pathways detected by ORA only:", ora_only_global),
      genes_unique_conn
    )
    writeLines(
      paste("Pathways detected by GSEA only:", gsea_only_global),
      genes_unique_conn
    )
    writeLines(
      paste("Pathways detected by Both (ORA & GSEA):", both_global),
      genes_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_unique_conn
    )
    writeLines(
      "ML Gene Set Overlaps in Pathways (Global):",
      genes_unique_conn
    )
    writeLines(
      paste("- Pathways with Triple Core genes:", tc_global_count),
      genes_unique_conn
    )
    writeLines(
      paste("- Pathways with Subnucleus genes:", sub_global_count),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Braaksc (Solo Tau):",
        uniq_br_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Ceradsc (Solo Amiloide):",
        uniq_ce_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with unique genes for Cogdx (Clínica):",
        uniq_co_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Ceradsc intersection (Resiliencia):",
        int_bc_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Braaksc-Cogdx intersection (Tau-Cogdx):",
        int_bg_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      paste(
        "- Pathways with Ceradsc-Cogdx intersection (Amiloide-Cogdx):",
        int_cg_global_count
      ),
      genes_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      genes_unique_conn
    )
    writeLines(
      paste(
        "Este reporte consolida en una lista integrada única todos los genes",
        "del conjunto de Machine Learning ordenados descendentemente por su",
        "nivel de conectividad funcional, mostrando para cada gen la totalidad",
        "de sus rutas biológicas asociadas."
      ),
      genes_unique_conn
    )
    writeLines(
      "--------------------------------------------------------------------\n",
      genes_unique_conn
    )

    for (g in seq_len(nrow(unique_genes_list))) {
      gene_name <- unique_genes_list$Gene[g]
      writeLines(paste0("  * GEN: ", gene_name), genes_unique_conn)
      writeLines(
        paste0("    Fragmento biológico: ", unique_genes_list$Fragment[g]),
        genes_unique_conn
      )
      writeLines(
        paste0(
          "    Total de rutas asociadas: ",
          unique_genes_list$Total_Pathways[g]
        ),
        genes_unique_conn
      )
      writeLines("    Rutas:", genes_unique_conn)

      # Filtrar rutas de este gen
      gene_paths <- gene_pathway_merged |>
        filter(.data$Gene == gene_name) |>
        arrange(.data$Description)

      for (p in seq_len(nrow(gene_paths))) {
        writeLines(
          paste0(
            "      - ", gene_paths$Description[p], " (",
            gene_paths$ID[p], ") [", gene_paths$Category[p],
            "] | Datasets: ", gene_paths$Datasets_Detected[p]
          ),
          genes_unique_conn
        )
      }
      writeLines("", genes_unique_conn)
    }

    close(genes_unique_conn)
    cat(
      paste(
        "      [ÉXITO] Reporte global de genes unificados guardado:",
        genes_unique_filename, "\n"
      )
    )
  }
}, error = function(e) {
  cat(
    sprintf(
      "      [ERROR] Falló generación reportes globales: %s\n",
      e$message
    )
  )
})

# --------------------------------------------------------------------
# COMPARACIÓN GLOBAL: NÚCLEOS ML VS DEA (FUERA DEL BUCLE DE DATASETS)
# --------------------------------------------------------------------
cat("\n>>> REALIZANDO COMPARACIÓN GLOBAL EXHAUSTIVA DE NÚCLEOS ML VS DEA...\n")
tryCatch({
  # Cargar significativos de DEA para los 3 datasets
  sig_dea_braak <- character(0)
  sig_dea_cerad <- character(0)
  sig_dea_cogdx <- character(0)

  if (file.exists("differential_expression/DEA_braaksc_complete.csv")) {
    df <- read_csv(
      "differential_expression/DEA_braaksc_complete.csv",
      show_col_types = FALSE
    )
    is_sig <- df$Significant == TRUE | df$Significant == "TRUE"
    sig_dea_braak <- df$Symbol[is_sig] |> na.omit() |> as.character()
  }
  if (file.exists("differential_expression/DEA_ceradsc_complete.csv")) {
    df <- read_csv(
      "differential_expression/DEA_ceradsc_complete.csv",
      show_col_types = FALSE
    )
    is_sig <- df$Significant == TRUE | df$Significant == "TRUE"
    sig_dea_cerad <- df$Symbol[is_sig] |> na.omit() |> as.character()
  }
  if (file.exists("differential_expression/DEA_cogdx_complete.csv")) {
    df <- read_csv(
      "differential_expression/DEA_cogdx_complete.csv",
      show_col_types = FALSE
    )
    is_sig <- df$Significant == TRUE | df$Significant == "TRUE"
    sig_dea_cogdx <- df$Symbol[is_sig] |> na.omit() |> as.character()
  }

  # Triple core de DEA: Genes significativos en los 3 datasets a la vez
  triple_core_dea <- intersect(
    intersect(sig_dea_braak, sig_dea_cerad),
    sig_dea_cogdx
  )

  # Comparaciones
  triple_core_ml <- triple_core
  subnucleus_ml  <- subnucleo

  # Genes que están en ML pero NO en DEA Core
  ml_not_in_dea_core <- setdiff(triple_core_ml, triple_core_dea)
  subnucleo_not_in_dea_core <- setdiff(subnucleus_ml, triple_core_dea)

  # 1. Crear Reporte Global (.txt)
  global_report_filename <- file.path(
    base_dir_landscape, "Data",
    "landscape_global_ml_vs_dea_comparison.txt"
  )
  global_conn <- file(global_report_filename, "w")
  writeLines(
    "====================================================================",
    global_conn
  )
  writeLines(
    "GLOBAL BIOLOGICAL DISCOVERY: MACHINE LEARNING (AI/ML) VS DEA CORES",
    global_conn
  )
  writeLines(
    "====================================================================",
    global_conn
  )
  writeLines(
    paste(
      "Este reporte analiza la diferencia entre los genes clave seleccionados",
      "por algoritmos de Inteligencia Artificial / Machine Learning (IA/ML)",
      "y los detectados por el análisis estadístico tradicional de",
      "Expresión Diferencial (DEA) en los 3 datasets."
    ),
    global_conn
  )
  writeLines(
    "--------------------------------------------------------------------",
    global_conn
  )
  writeLines(
    "1. TRIPLE CORE COMPARISON (Genes de ML vs Genes comunes de DEA):",
    global_conn
  )
  writeLines(
    paste(
      "   - Total Triple Core ML (Core Masters de ML):",
      length(triple_core_ml)
    ),
    global_conn
  )
  writeLines(
    paste(
      "   - Total Triple Core DEA (Comunes en los 3 DEAs):",
      length(triple_core_dea)
    ),
    global_conn
  )
  writeLines(
    paste(
      paste(
        "   - Genes del Triple Core ML que NO están en",
        "DEA Core (Descubrimientos):"
      ),
      length(ml_not_in_dea_core)
    ),
    global_conn
  )
  writeLines(
    paste("     ->", paste(ml_not_in_dea_core, collapse = ", ")),
    global_conn
  )
  writeLines(
    "--------------------------------------------------------------------",
    global_conn
  )
  writeLines(
    paste(
      "2. SUBNUCLEUS COMPARISON (Subnácleo Irreducible vs",
      "Genes comunes de DEA):"
    ),
    global_conn
  )
  writeLines(
    paste(
      "   - Total genes en el Subnúcleo Irreducible de ML:",
      length(subnucleus_ml)
    ),
    global_conn
  )
  writeLines(
    paste(
      "   - Genes del Subnúcleo de ML que NO están en el Triple Core DEA:",
      length(subnucleo_not_in_dea_core)
    ),
    global_conn
  )
  writeLines(
    paste("     ->", paste(subnucleo_not_in_dea_core, collapse = ", ")),
    global_conn
  )
  writeLines(
    "--------------------------------------------------------------------",
    global_conn
  )
  writeLines(
    "3. ANÁLISIS DE SIGNIFICACIÓN DE DEA PARA EL SUBNÚCLEO IRREDUCIBLE:",
    global_conn
  )
  for (g in subnucleus_ml) {
    in_braak <- g %in% sig_dea_braak
    in_cerad <- g %in% sig_dea_cerad
    in_cogdx <- g %in% sig_dea_cogdx
    writeLines(
      paste0(
        "  * ", g, " : ",
        "Significant in Braaksc (", ifelse(in_braak, "YES", "NO"), ") | ",
        "Ceradsc (", ifelse(in_cerad, "YES", "NO"), ") | ",
        "Cogdx (", ifelse(in_cogdx, "YES", "NO"), ")"
      ),
      global_conn
    )
  }
  close(global_conn)

  # 2. Crear CSV de Trazabilidad Global ML vs DEA
  global_comparison_df <- data.frame(
    Gene_Symbol = unique(c(triple_core_ml, subnucleus_ml))
  ) |>
    mutate(
      Is_In_Triple_Core_ML = .data$Gene_Symbol %in% triple_core_ml,
      Is_In_Subnucleus_ML  = .data$Gene_Symbol %in% subnucleus_ml,
      Sig_in_DEA_Braaksc   = .data$Gene_Symbol %in% sig_dea_braak,
      Sig_in_DEA_Ceradsc   = .data$Gene_Symbol %in% sig_dea_cerad,
      Sig_in_DEA_Cogdx     = .data$Gene_Symbol %in% sig_dea_cogdx,
      Is_In_Triple_Core_DEA = .data$Gene_Symbol %in% triple_core_dea,
      Missed_by_DEA_Core   = !.data$Is_In_Triple_Core_DEA
    )

  global_csv_filename <- file.path(
    base_dir_landscape, "Data", "landscape_global_ml_vs_dea_comparison.csv"
  )
  global_comparison_csv_df <- global_comparison_df |>
    dplyr::rename(
      Gene = .data$Gene_Symbol,
      `Is In Triple Core ML` = .data$Is_In_Triple_Core_ML,
      `Is In Subnucleus ML` = .data$Is_In_Subnucleus_ML,
      `Sig in DEA Braaksc` = .data$Sig_in_DEA_Braaksc,
      `Sig in DEA Ceradsc` = .data$Sig_in_DEA_Ceradsc,
      `Sig in DEA Cogdx` = .data$Sig_in_DEA_Cogdx,
      `Is In Triple Core DEA` = .data$Is_In_Triple_Core_DEA,
      `Missed by DEA Core` = .data$Missed_by_DEA_Core
    )
  write_excel_csv(global_comparison_csv_df, global_csv_filename)

  cat(
    paste(
      "      [ÉXITO] Reporte global de ML vs DEA guardado:",
      global_report_filename, "\n"
    )
  )
  cat(
    paste(
      "      [ÉXITO] CSV global de ML vs DEA guardado:",
      global_csv_filename, "\n"
    )
  )
}, error = function(e) {
  cat(
    sprintf(
      "      [ERROR] Falló la comparación global de ML vs DEA: %s\n",
      e$message
    )
  )
})

# ====================================================================
# REPORTE DE RUTAS GSEA SIGNIFICATIVAS CON COMPARATIVA CRUZADA DE DEGs Y RANGOS
# ====================================================================
cat(
  paste(
    "\n>>> GENERANDO REPORTE DETALLADO DE RUTAS GSEA",
    "SIGNIFICATIVAS SIN DEGs...\n"
  )
)
tryCatch({
  # Cargar archivos DEA completos y añadir rankings
  load_dea_r <- function(path) {
    if (!file.exists(path)) {
      NULL
    } else {
      df <- read_csv(path, show_col_types = FALSE)
      df$Rank_PVal <- seq_len(nrow(df))
      df
    }
  }

  dea_b <- load_dea_r("differential_expression/DEA_braaksc_complete.csv")
  dea_c <- load_dea_r("differential_expression/DEA_ceradsc_complete.csv")
  dea_g <- load_dea_r("differential_expression/DEA_cogdx_complete.csv")

  get_sig_metrics_r <- function(df) {
    if (is.null(df)) {
      list(tot = 0, sig = 0, max_r = 0)
    } else {
      tot <- nrow(df)
      sig_df <- df[df$Significant == TRUE | df$Significant == "TRUE", ]
      sig <- nrow(sig_df)
      max_r <- if (sig > 0) max(sig_df$Rank_PVal, na.rm = TRUE) else 0
      list(tot = tot, sig = sig, max_r = max_r)
    }
  }

  met_b <- get_sig_metrics_r(dea_b)
  met_c <- get_sig_metrics_r(dea_c)
  met_g <- get_sig_metrics_r(dea_g)

  report_filename <- file.path(
    base_dir_landscape, "Data", "landscape_gsea_no_dea_genes.txt"
  )
  conn <- file(report_filename, "w", encoding = "UTF-8")

  writeLines(
    "====================================================================",
    conn
  )
  writeLines(
    paste(
      "REPORTE DE RUTAS GSEA SIGNIFICATIVAS CON COMPARATIVA",
      "CRUZADA DE DEGs Y RANGOS"
    ),
    conn
  )
  writeLines(
    "====================================================================",
    conn
  )
  writeLines(
    paste(
      "Este reporte agrupa y describe las rutas identificadas como",
      "significativas en GSEA pero cuyos genes de ML solapantes NO",
      "alcanzan la significancia en DEA."
    ),
    conn
  )
  writeLines(
    paste(
      "Para cada gen solapante, se detalla su comportamiento completo",
      "(Expresión, FDR, Rango y Significancia) cruzado en los tres",
      "datasets (BRAAKSC, CERADSC, COGDX), así como la última posición",
      "(punto de corte) de significancia para cada uno."
    ),
    conn
  )
  writeLines(
    "====================================================================\n",
    conn
  )

  writeLines(
    "--------------------------------------------------------------------",
    conn
  )
  writeLines(
    "LÍMITES Y CRITERIOS DE SELECCIÓN PARA GENES SIGNIFICATIVOS (DEA):",
    conn
  )
  writeLines(
    "--------------------------------------------------------------------",
    conn
  )
  writeLines("- FDR (p-adj / adj.P.Val) límite: < 0.05", conn)
  writeLines(
    "- Log2 Fold Change (Log2FC) límite: valor absoluto > 0.5",
    conn
  )
  writeLines(
    paste(
      "- Límite en Posición de Ranking: No hay un límite",
      "de ranking para ser DEG;"
    ),
    conn
  )
  writeLines(
    paste(
      "  todos los genes que cumplen FDR < 0.05 y",
      "abs(logFC) > 0.5 son significativos."
    ),
    conn
  )
  writeLines(
    "--------------------------------------------------------------------",
    conn
  )
  writeLines(
    paste(
      "CONTEXTO HISTÓRICO Y PUNTOS DE CORTE",
      "(ÚLTIMO DEG SIGNIFICATIVO):"
    ),
    conn
  )
  writeLines(
    sprintf(
      paste(
        "- BRAAKSC: %d DEGs significativos de %d genes.",
        "Rango del útimo DEG: #%d"
      ),
      met_b$sig, met_b$tot, met_b$max_r
    ),
    conn
  )
  writeLines(
    sprintf(
      paste(
        "- CERADSC: %d DEGs significativos de %d genes.",
        "Rango del útimo DEG: #%d"
      ),
      met_c$sig, met_c$tot, met_c$max_r
    ),
    conn
  )
  writeLines(
    sprintf(
      paste(
        "- COGDX:   %d DEGs significativos de %d genes.",
        "Rango del útimo DEG: #%d"
      ),
      met_g$sig, met_g$tot, met_g$max_r
    ),
    conn
  )
  writeLines(
    paste(
      "*(Los genes están ordenados en el DEA por p-valor",
      "original de forma ascendente)"
    ),
    conn
  )
  writeLines(
    "--------------------------------------------------------------------\n",
    conn
  )

  datasets_info <- list(
    list(
      name = "BRAAKSC",
      csv = file.path(
        base_dir_landscape, "Data",
        "gsea_no_dea_sig_pathways_braaksc.csv"
      ),
      dea = dea_b
    ),
    list(
      name = "CERADSC",
      csv = file.path(
        base_dir_landscape, "Data",
        "gsea_no_dea_sig_pathways_ceradsc.csv"
      ),
      dea = dea_c
    ),
    list(
      name = "COGDX",
      csv = file.path(
        base_dir_landscape, "Data",
        "gsea_no_dea_sig_pathways_cogdx.csv"
      ),
      dea = dea_g
    )
  )

  format_gene_details_r <- function(sym, ds_name, df, max_rank, conn) {
    if (is.null(df)) {
      writeLines(
        sprintf("      * [ERROR] No se pudo cargar el dataset %s.", ds_name),
        conn
      )
    } else {
      g_row <- df[which(df$Symbol == sym), ]
      if (nrow(g_row) == 0) {
        writeLines(
          sprintf(
            "      * ¿Es DEG Significativo en %s? NO (Gen no detectado)",
            ds_name
          ),
          conn
        )
      } else {
        logfc <- g_row$logFC[1]
        padj <- g_row$adj.P.Val[1]
        pval <- g_row$P.Value[1]
        t_stat <- g_row$t[1]
        rank_pval <- g_row$Rank_PVal[1]
        sig_status <- g_row$Significant[1]

        sig_txt <- if (sig_status == TRUE || sig_status == "TRUE") {
          "SÍ"
        } else {
          "NO"
        }

        writeLines(
          sprintf("      * ¿Es DEG Significativo en %s? %s", ds_name, sig_txt),
          conn
        )
        writeLines(
          sprintf(
            "        - Rango por p-valor: #%d de %d (Último DEG: #%d)",
            rank_pval, nrow(df), max_rank
          ),
          conn
        )
        writeLines(sprintf("        - Log2 Fold Change: %.6f", logfc), conn)
        writeLines(sprintf("        - FDR (adj. p-valor): %.4e", padj), conn)
        writeLines(sprintf("        - p-valor original: %.4e", pval), conn)
        writeLines(
          sprintf("        - Estadístico t de Limma: %.4f", t_stat),
          conn
        )
      }
    }
  }

  for (ds_info in datasets_info) {
    writeLines(
      "====================================================================",
      conn
    )
    writeLines(sprintf("DATASET DE ORIGEN: %s", ds_info$name), conn)
    writeLines(
      "====================================================================",
      conn
    )

    if (!file.exists(ds_info$csv)) {
      writeLines("Archivo CSV de rutas no encontrado.\n", conn)
      next
    }

    df_pathways <- read_csv(ds_info$csv, show_col_types = FALSE)
    if (nrow(df_pathways) == 0) {
      writeLines(
        paste0(
          "No se detectaron rutas significativas en esta ",
          "condición GSEA sin DEG.\n"
        ),
        conn
      )
      next
    }

    writeLines(
      sprintf("Total de rutas en esta condición: %d", nrow(df_pathways)),
      conn
    )
    writeLines(
      "--------------------------------------------------------------------",
      conn
    )

    for (idx in seq_len(nrow(df_pathways))) {
      path_id <- df_pathways$ID[idx]
      desc <- df_pathways$Descripcion[idx]
      cat_db <- df_pathways$Categoria[idx]
      genes_str <- df_pathways$Todos_los_genes[idx]

      writeLines(
        sprintf("\n* RUTA: %s (%s) [%s]", desc, path_id, cat_db),
        conn
      )
      writeLines(
        "  Detalle de Expresión y Significancia Cruzada de los Genes:",
        conn
      )

      if (is.na(genes_str) || genes_str == "") {
        writeLines("  (Ningún gen listado)", conn)
        next
      }

      genes_list <- unlist(strsplit(genes_str, "/"))

      # Ordenar genes según su rango en el dataset de origen
      gene_ranks <- data.frame(
        Symbol = genes_list,
        Rank = 999999,
        stringsAsFactors = FALSE
      )
      if (!is.null(ds_info$dea)) {
        for (g_idx in seq_len(nrow(gene_ranks))) {
          g <- gene_ranks$Symbol[g_idx]
          g_row <- ds_info$dea[
            which(ds_info$dea$Symbol == g),
          ]
          if (nrow(g_row) > 0) {
            gene_ranks$Rank[g_idx] <- g_row$Rank_PVal[1]
          }
        }
      }
      gene_ranks <- gene_ranks[order(gene_ranks$Rank), ]
      sorted_genes <- gene_ranks$Symbol

      for (g in sorted_genes) {
        writeLines(sprintf("\n    - GEN: %s", g), conn)
        format_gene_details_r(g, "BRAAKSC", dea_b, met_b$max_r, conn)
        format_gene_details_r(g, "CERADSC", dea_c, met_c$max_r, conn)
        format_gene_details_r(g, "COGDX", dea_g, met_g$max_r, conn)
      }
    }
    writeLines("\n", conn)
  }

  close(conn)
  cat(
    paste(
      "      [ÉXITO] Reporte de rutas GSEA sin DEGs guardado:",
      report_filename, "\n"
    )
  )
}, error = function(e) {
  cat(
    sprintf(
      "      [ERROR] Falló generación reporte de rutas GSEA sin DEGs: %s\n",
      e$message
    )
  )
})

cat("\nEjecución del script de entorno de conjuntos de genes completada.\n")
