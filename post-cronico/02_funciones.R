# ============================================================
# 02_funciones.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# ============================================================

preparar_deseq <- function(dds_obj, grupo_ref, grupo_comp) {
  seleccion <- !is.na(colData(dds_obj)$grupo_biologico) &
    colData(dds_obj)$grupo_biologico %in% c(grupo_ref, grupo_comp)
  dds_sub <- dds_obj[, seleccion]

  if (ncol(dds_sub) == 0) {
    stop("No hay muestras para la comparación ", grupo_comp, " vs ", grupo_ref, ".")
  }

  dds_sub$grupo_biologico <- factor(colData(dds_sub)$grupo_biologico, levels = c(grupo_ref, grupo_comp))
  dds_sub$dataset <- droplevels(factor(colData(dds_sub)$dataset))

  n_por_grupo <- table(colData(dds_sub)$grupo_biologico)
  if (any(n_por_grupo == 0)) {
    stop(
      "No se puede calcular ", grupo_comp, " vs ", grupo_ref,
      ": uno de los dos grupos no tiene muestras."
    )
  }
  if (any(n_por_grupo < 2)) {
    warning("Hay menos de dos muestras en uno de los grupos; interpretar con mucha cautela.")
  }

  # Diseño: se prueban modelos de más a menos completo y se usa el primero
  # con rango completo. Cada covariable entra solo si tiene >1 nivel y sin NA.
  datos_modelo <- as.data.frame(colData(dds_sub))
  covariable_valida <- function(v) {
    v %in% names(datos_modelo) && !any(is.na(datos_modelo[[v]])) &&
      length(unique(datos_modelo[[v]])) > 1
  }
  candidatos <- list(c("dataset", "tejido"), c("tejido"), c("dataset"), character(0))
  if (!isTRUE(config$usar_tejido_en_diseno)) candidatos <- list(c("dataset"), character(0))

  formula_final <- ~ grupo_biologico
  for (covs in candidatos) {
    if (!all(vapply(covs, covariable_valida, logical(1)))) next
    for (v in covs) datos_modelo[[v]] <- droplevels(factor(datos_modelo[[v]]))
    f <- as.formula(paste("~", paste(c(covs, "grupo_biologico"), collapse = " + ")))
    mm <- model.matrix(f, data = datos_modelo)
    if (qr(mm)$rank == ncol(mm)) {
      for (v in covs) colData(dds_sub)[[v]] <- datos_modelo[[v]]
      formula_final <- f
      break
    }
  }
  design(dds_sub) <- formula_final
  cat("Diseño usado:", deparse(formula_final), "\n")
  if (!exists("disenos_usados", envir = .GlobalEnv)) assign("disenos_usados", list(), envir = .GlobalEnv)
  disenos_usados[[paste(grupo_comp, "vs", grupo_ref)]] <<- deparse(formula_final)

  return(dds_sub)
}

ejecutar_deseq <- function(dds_obj, grupo_ref, grupo_comp, nombre_comparacion, anotacion_genes) {
  cat("\nEjecutando DESeq2:", nombre_comparacion, "\n")
  dds_sub <- preparar_deseq(dds_obj, grupo_ref, grupo_comp)
  print(table(colData(dds_sub)$dataset, colData(dds_sub)$grupo_biologico))
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("grupo_biologico", grupo_comp, grupo_ref))
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  res_df <- res_df %>%
    mutate(ensembl_id_limpio = sub("\\..*", "", ensembl_id), comparacion = nombre_comparacion) %>%
    left_join(anotacion_genes, by = "ensembl_id_limpio") %>%
    mutate(
      simbolo = ifelse(!is.na(simbolo_mouse) & simbolo_mouse != "", simbolo_mouse, ensembl_id_limpio),
      inflamatorio = ifelse(is.na(inflamatorio), "No", inflamatorio),
      fold_change = 2^log2FoldChange,
      cambio = case_when(
        padj < config$padj_significativo & log2FoldChange > 1 ~ "Sube",
        padj < config$padj_significativo & log2FoldChange < -1 ~ "Baja",
        TRUE ~ "No significativo"
      )
    )
  return(res_df)
}

convertir_entrez <- function(res_df) {
  cat("Convirtiendo ENSEMBL a ENTREZ...\n")
  conversion <- bitr(
    unique(res_df$ensembl_id_limpio),
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Mm.eg.db
  )
  res_entrez <- res_df %>%
    left_join(conversion, by = c("ensembl_id_limpio" = "ENSEMBL")) %>%
    mutate(simbolo = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, simbolo))
  return(res_entrez)
}

crear_ranking <- function(res_entrez) {
  ranking <- res_entrez %>%
    filter(!is.na(stat), !is.na(ENTREZID)) %>%
    arrange(desc(abs(stat))) %>%
    group_by(ENTREZID) %>%
    slice_head(n = 1) %>%
    ungroup()
  gene_list <- ranking$stat
  names(gene_list) <- ranking$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

ejecutar_GO <- function(gene_list) {
  cat("Ejecutando GSEA GO BP...\n")
  set.seed(config$semilla)
  gseGO(
    geneList = gene_list,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = config$minGSSize,
    maxGSSize = config$maxGSSize,
    pvalueCutoff = config$gsea_pvalueCutoff,
    verbose = FALSE
  )
}

ejecutar_KEGG <- function(gene_list) {
  cat("Ejecutando GSEA KEGG...\n")
  set.seed(config$semilla)
  gseKEGG(
    geneList = gene_list,
    organism = "mmu",
    minGSSize = config$minGSSize,
    maxGSSize = config$maxGSSize,
    pvalueCutoff = config$gsea_pvalueCutoff,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

formatear_gsea <- function(gsea_obj, nombre_comparacion, base) {
  df <- as.data.frame(gsea_obj)
  if (nrow(df) == 0) return(data.frame())
  df %>%
    mutate(
      comparacion = nombre_comparacion,
      base = base,
      direccion = ifelse(NES > 0, "Enriquecido en grupo comparado", "Enriquecido en grupo referencia"),
      n_leading_edge = sapply(strsplit(core_enrichment, "/"), length)
    ) %>%
    dplyr::select(
      comparacion, base, ID, Description, NES, pvalue, p.adjust, qvalue,
      setSize, n_leading_edge, direccion, core_enrichment
    )
}

traducir_core_enrichment <- function(gsea_df, res_entrez) {
  if (nrow(gsea_df) == 0) return(gsea_df)
  mapa <- res_entrez %>%
    filter(!is.na(ENTREZID)) %>%
    dplyr::select(ENTREZID, simbolo) %>%
    distinct()
  traducir <- function(x) {
    ids <- unlist(strsplit(x, "/"))
    simbolos <- mapa$simbolo[match(ids, mapa$ENTREZID)]
    simbolos <- unique(na.omit(simbolos))
    paste(simbolos, collapse = ", ")
  }
  gsea_df %>%
    mutate(Genes = sapply(core_enrichment, traducir))
}

extraer_leading_edge <- function(gsea_df, res_entrez) {
  if (nrow(gsea_df) == 0) return(data.frame())
  res_simple <- res_entrez %>%
    dplyr::select(
      comparacion, ensembl_id_limpio, ENTREZID, simbolo, biotipo,
      inflamatorio, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    ) %>%
    distinct()
  gsea_df %>%
    dplyr::select(comparacion, base, ID, Description, NES, p.adjust, direccion, core_enrichment) %>%
    separate_rows(core_enrichment, sep = "/") %>%
    rename(ENTREZID = core_enrichment, Proceso = Description, padj_proceso = p.adjust) %>%
    left_join(res_simple, by = c("comparacion", "ENTREZID")) %>%
    transmute(
      comparacion, base, ID, Proceso, Gen = simbolo, ensembl_id_limpio, ENTREZID,
      log2FC = log2FoldChange, stat, pvalue, padj, NES, padj_proceso,
      direccion, inflamatorio, biotipo
    )
}

resumir_genes_en_procesos <- function(tabla_leading) {
  if (nrow(tabla_leading) == 0) return(data.frame())
  tabla_leading %>%
    group_by(comparacion, Gen, ensembl_id_limpio, ENTREZID) %>%
    summarise(
      log2FC = first(log2FC),
      stat = first(stat),
      pvalue = first(pvalue),
      padj = first(padj),
      inflamatorio = first(inflamatorio),
      GO = paste(unique(Proceso[base == "GO"]), collapse = " | "),
      KEGG = paste(unique(Proceso[base == "KEGG"]), collapse = " | "),
      n_GO = n_distinct(Proceso[base == "GO"]),
      n_KEGG = n_distinct(Proceso[base == "KEGG"]),
      n_total_procesos = n_distinct(paste(base, Proceso)),
      .groups = "drop"
    ) %>%
    arrange(desc(n_total_procesos), desc(abs(log2FC)))
}

preparar_deseq_exportar <- function(res_entrez) {
  res_entrez %>%
    dplyr::select(
      comparacion, simbolo, ensembl_id, ensembl_id_limpio, ENTREZID,
      biotipo, inflamatorio, simbolo_humano, baseMean,
      log2FoldChange, lfcSE, stat, pvalue, padj, fold_change, cambio
    ) %>%
    arrange(padj)
}

preparar_gsea_exportar <- function(gsea_df) {
  if (nrow(gsea_df) == 0) return(data.frame())
  gsea_df %>%
    dplyr::select(
      comparacion, base, ID, Description, NES, pvalue, p.adjust,
      qvalue, setSize, n_leading_edge, direccion, Genes
    ) %>%
    arrange(p.adjust)
}

