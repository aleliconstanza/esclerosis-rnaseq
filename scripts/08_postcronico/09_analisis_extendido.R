# ============================================================
# 09_analisis_extendido.R
# DESeq2 + GSEA GO/KEGG para las comparaciones extendidas
# (tejido y sexo). Reutiliza SIN MODIFICAR las funciones de
# 02_funciones.R: preparar_deseq, ejecutar_deseq, convertir_entrez,
# crear_ranking, ejecutar_GO, ejecutar_KEGG, formatear_gsea,
# traducir_core_enrichment, extraer_leading_edge,
# resumir_genes_en_procesos, preparar_deseq_exportar,
# preparar_gsea_exportar.
# Requiere en la sesion: dds (extendido, con grupo_ext/grupo_sexo),
# comparaciones_extendidas (08), anotacion_genes (01), y
# 02_funciones.R ya sourceado.
# ============================================================

cat("\nIniciando analisis extendido (tejido / sexo)...\n")

origen_grupo <- setNames(
  c(rep("grupo_ext", length(comparaciones_tejido)),
    rep("grupo_sexo", length(comparaciones_sexo))),
  names(comparaciones_extendidas)
)

resultados_analisis_ext <- list()

for (nombre_comp in names(comparaciones_extendidas)) {
  cat("\n============================================================\n")
  cat("Comparacion:", nombre_comp, "\n")
  cat("============================================================\n")

  grupo_ref <- comparaciones_extendidas[[nombre_comp]]$grupo_ref
  grupo_comp <- comparaciones_extendidas[[nombre_comp]]$grupo_comp

  # Se usa la columna correspondiente (tejido o sexo) como
  # grupo_biologico SOLO para esta iteracion, sobre una copia
  dds_iter <- dds
  colData(dds_iter)$grupo_biologico <- colData(dds_iter)[[origen_grupo[[nombre_comp]]]]

  res_deseq <- tryCatch(
    ejecutar_deseq(
      dds_obj = dds_iter,
      grupo_ref = grupo_ref,
      grupo_comp = grupo_comp,
      nombre_comparacion = nombre_comp,
      anotacion_genes = anotacion_genes
    ),
    error = function(e) {
      message("Se omitio ", nombre_comp, ": ", e$message)
      NULL
    }
  )
  if (is.null(res_deseq)) next

  res_entrez <- convertir_entrez(res_deseq)
  gene_list <- crear_ranking(res_entrez)

  cat("Genes usados en ranking GSEA:", length(gene_list), "\n")

  if (length(gene_list) < config$minGSSize) {
    message("Ranking insuficiente para GSEA en ", nombre_comp, ", se omite.")
    next
  }

  gsea_go_obj <- ejecutar_GO(gene_list)
  gsea_kegg_obj <- ejecutar_KEGG(gene_list)

  go_df <- formatear_gsea(gsea_go_obj, nombre_comp, "GO")
  kegg_df <- formatear_gsea(gsea_kegg_obj, nombre_comp, "KEGG")

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

  resultados_analisis_ext[[nombre_comp]] <- list(
    comparacion = nombre_comp,
    grupo_ref = grupo_ref,
    grupo_comp = grupo_comp,
    deseq = deseq_exportar,
    go = go_exportar,
    kegg = kegg_exportar,
    leading_edge = leading_total,
    genes_en_procesos = genes_procesos,
    resumen = resumen_comp,
    objetos_gsea = list(go = gsea_go_obj, kegg = gsea_kegg_obj),
    res_entrez = res_entrez
  )

  cat("\nResumen de la comparacion:\n")
  print(resumen_comp)
}

cat("\nAnalisis extendido finalizado. Comparaciones calculadas:",
    length(resultados_analisis_ext), "de", length(comparaciones_extendidas), "\n")
