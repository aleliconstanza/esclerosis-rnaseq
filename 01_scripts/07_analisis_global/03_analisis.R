# ============================================================
# 03_analisis.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# Ejecución de DESeq2, ranking, GSEA GO y GSEA KEGG
# ============================================================

cat("\nIniciando análisis principal...\n")

resultados_analisis <- list()

for (nombre_comp in names(comparaciones)) {
  cat("\n============================================================\n")
  cat("Comparación:", nombre_comp, "\n")
  cat("============================================================\n")

  grupo_ref <- comparaciones[[nombre_comp]]$grupo_ref
  grupo_comp <- comparaciones[[nombre_comp]]$grupo_comp

  res_deseq <- ejecutar_deseq(
    dds_obj = dds,
    grupo_ref = grupo_ref,
    grupo_comp = grupo_comp,
    nombre_comparacion = nombre_comp,
    anotacion_genes = anotacion_genes
  )

  res_entrez <- convertir_entrez(res_deseq)
  gene_list <- crear_ranking(res_entrez)

  cat("Genes usados en ranking GSEA:", length(gene_list), "\n")

  if (length(gene_list) < config$minGSSize) {
    stop("El ranking de ", nombre_comp, " no contiene suficientes genes para GSEA.")
  }

  gsea_go_obj <- ejecutar_GO(gene_list)
  gsea_kegg_obj <- ejecutar_KEGG(gene_list)

  go_df <- formatear_gsea(
    gsea_obj = gsea_go_obj,
    nombre_comparacion = nombre_comp,
    base = "GO"
  )

  kegg_df <- formatear_gsea(
    gsea_obj = gsea_kegg_obj,
    nombre_comparacion = nombre_comp,
    base = "KEGG"
  )

  go_df <- traducir_core_enrichment(go_df, res_entrez)
  kegg_df <- traducir_core_enrichment(kegg_df, res_entrez)

  leading_go <- extraer_leading_edge(go_df, res_entrez)
  leading_kegg <- extraer_leading_edge(kegg_df, res_entrez)

  leading_total <- bind_rows(leading_go, leading_kegg)
  genes_procesos <- resumir_genes_en_procesos(leading_total)

  deseq_exportar <- preparar_deseq_exportar(res_entrez)
  go_exportar <- preparar_gsea_exportar(go_df)
  kegg_exportar <- preparar_gsea_exportar(kegg_df)

  resumen_comp <- data.frame(
    comparacion = nombre_comp,
    grupo_referencia = grupo_ref,
    grupo_comparado = grupo_comp,
    genes_analizados = nrow(deseq_exportar),
    genes_ranking_gsea = length(gene_list),
    genes_diferenciales = sum(deseq_exportar$cambio != "No significativo", na.rm = TRUE),
    genes_sube = sum(deseq_exportar$cambio == "Sube", na.rm = TRUE),
    genes_baja = sum(deseq_exportar$cambio == "Baja", na.rm = TRUE),
    procesos_GO = nrow(go_exportar),
    rutas_KEGG = nrow(kegg_exportar),
    genes_en_procesos = nrow(genes_procesos),
    genes_inflamatorios_en_procesos = sum(genes_procesos$inflamatorio == "Sí", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  resultados_analisis[[nombre_comp]] <- list(
    comparacion = nombre_comp,
    grupo_ref = grupo_ref,
    grupo_comp = grupo_comp,
    deseq = deseq_exportar,
    go = go_exportar,
    kegg = kegg_exportar,
    leading_edge = leading_total,
    genes_en_procesos = genes_procesos,
    resumen = resumen_comp,
    objetos_gsea = list(
      go = gsea_go_obj,
      kegg = gsea_kegg_obj
    ),
    res_entrez = res_entrez
  )

  cat("\nResumen de la comparación:\n")
  print(resumen_comp)
}

cat("\nAnálisis principal finalizado correctamente.\n")

