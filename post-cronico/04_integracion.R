# ============================================================
# 04_integracion.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# Integración de resultados para exportación final
# ============================================================

cat("\nIntegrando resultados...\n")

tabla_resumen <- bind_rows(lapply(resultados_analisis, function(x) x$resumen))
tabla_deseq <- bind_rows(lapply(resultados_analisis, function(x) x$deseq))
tabla_go <- bind_rows(lapply(resultados_analisis, function(x) x$go))
tabla_kegg <- bind_rows(lapply(resultados_analisis, function(x) x$kegg))
tabla_leading_edge <- bind_rows(lapply(resultados_analisis, function(x) x$leading_edge))
tabla_genes_en_procesos <- bind_rows(lapply(resultados_analisis, function(x) x$genes_en_procesos))

tabla_genes_inflamatorios <- if (nrow(tabla_genes_en_procesos) > 0) {
  tabla_genes_en_procesos %>%
    filter(inflamatorio == "Sí") %>%
    arrange(desc(n_total_procesos), desc(abs(log2FC)))
} else {
  data.frame()
}

resumir_procesos_principales <- function(tabla, nombre_columna) {
  if (nrow(tabla) == 0) {
    salida <- data.frame(comparacion = character(), stringsAsFactors = FALSE)
    salida[[nombre_columna]] <- character()
    return(salida)
  }

  tabla %>%
    group_by(comparacion) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    summarise(
      valor = paste(Description, collapse = " | "),
      .groups = "drop"
    ) %>%
    rename(!!nombre_columna := valor)
}

procesos_principales_go <- resumir_procesos_principales(tabla_go, "principales_GO")
procesos_principales_kegg <- resumir_procesos_principales(tabla_kegg, "principales_KEGG")

tabla_resumen <- tabla_resumen %>%
  left_join(procesos_principales_go, by = "comparacion") %>%
  left_join(procesos_principales_kegg, by = "comparacion")

tabla_frecuencia_genes <- if (nrow(tabla_leading_edge) > 0) {
  tabla_leading_edge %>%
    group_by(comparacion, Gen, ensembl_id_limpio, ENTREZID, inflamatorio) %>%
    summarise(
      n_procesos_total = n_distinct(paste(base, Proceso)),
      n_GO = n_distinct(Proceso[base == "GO"]),
      n_KEGG = n_distinct(Proceso[base == "KEGG"]),
      log2FC = first(log2FC),
      stat = first(stat),
      padj = first(padj),
      procesos = paste(unique(Proceso), collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(comparacion, desc(n_procesos_total), desc(abs(log2FC)))
} else {
  data.frame()
}

diccionario <- data.frame(
  columna = c(
    "comparacion", "baseMean", "log2FoldChange", "log2FC", "lfcSE",
    "stat", "pvalue", "padj", "NES", "p.adjust", "qvalue",
    "setSize", "n_leading_edge", "core_enrichment", "Genes",
    "GO", "KEGG", "n_GO", "n_KEGG", "n_total_procesos",
    "inflamatorio", "direccion"
  ),
  significado = c(
    "Nombre de la comparación biológica realizada.",
    "Promedio normalizado de expresión del gen en todas las muestras.",
    "Cambio de expresión en escala log2 estimado por DESeq2.",
    "Misma información que log2FoldChange, usada en tablas resumidas.",
    "Error estándar del log2FoldChange.",
    "Estadístico de Wald de DESeq2. Se usa como ranking para GSEA.",
    "Valor p sin ajustar.",
    "Valor p ajustado por múltiples comparaciones.",
    "Normalized Enrichment Score. Indica fuerza y dirección del enriquecimiento.",
    "Valor p ajustado del análisis GSEA.",
    "Valor q reportado por GSEA.",
    "Número de genes de la ruta o proceso.",
    "Número de genes del leading edge.",
    "Genes ENTREZ que explican el enriquecimiento.",
    "Genes del leading edge traducidos a símbolos oficiales.",
    "Procesos GO en los que participa el gen.",
    "Rutas KEGG en las que participa el gen.",
    "Número de procesos GO donde aparece el gen.",
    "Número de rutas KEGG donde aparece el gen.",
    "Número total de procesos/rutas donde aparece el gen.",
    "Indica si el gen pertenece a la lista inflamatoria de literatura.",
    "Indica si la ruta está enriquecida en el grupo comparado o en el grupo referencia."
  ),
  stringsAsFactors = FALSE
)

cat("\nTablas integradas correctamente:\n")
cat("- tabla_resumen\n")
cat("- tabla_deseq\n")
cat("- tabla_go\n")
cat("- tabla_kegg\n")
cat("- tabla_genes_en_procesos\n")
cat("- tabla_genes_inflamatorios\n")
cat("- tabla_frecuencia_genes\n")
cat("- diccionario\n")

