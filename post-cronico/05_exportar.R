# ============================================================
# 05_exportar.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# Exportación de Excel, tablas y figuras
# ============================================================

cat("\nExportando resultados...\n")

ruta_excel <- file.path(config$carpeta_excel, config$archivo_excel)

wb <- createWorkbook()
agregar_hoja <- function(wb, nombre, tabla) {
  addWorksheet(wb, nombre)

  if (ncol(tabla) == 0) {
    writeData(wb, nombre, data.frame(Mensaje = "Sin resultados para esta tabla."))
    return(invisible(NULL))
  }

  writeData(wb, nombre, tabla)
  freezePane(wb, nombre, firstRow = TRUE)
  setColWidths(wb, nombre, cols = seq_len(ncol(tabla)), widths = "auto")
}

agregar_hoja(wb, "Resumen", tabla_resumen)
agregar_hoja(wb, "DESeq2", tabla_deseq)
agregar_hoja(wb, "GO", tabla_go)
agregar_hoja(wb, "KEGG", tabla_kegg)
agregar_hoja(wb, "Genes_en_procesos", tabla_genes_en_procesos)
agregar_hoja(wb, "Genes_inflamatorios", tabla_genes_inflamatorios)
agregar_hoja(wb, "Frecuencia_genes", tabla_frecuencia_genes)
agregar_hoja(wb, "Leading_edge", tabla_leading_edge)
agregar_hoja(wb, "Diccionario", diccionario)

saveWorkbook(wb, ruta_excel, overwrite = TRUE)

write.csv2(tabla_resumen, file.path(config$carpeta_tablas, "Resumen.csv"), row.names = FALSE)
write.csv2(tabla_deseq, file.path(config$carpeta_tablas, "DESeq2.csv"), row.names = FALSE)
write.csv2(tabla_go, file.path(config$carpeta_tablas, "GO.csv"), row.names = FALSE)
write.csv2(tabla_kegg, file.path(config$carpeta_tablas, "KEGG.csv"), row.names = FALSE)
write.csv2(tabla_genes_en_procesos, file.path(config$carpeta_tablas, "Genes_en_procesos.csv"), row.names = FALSE)
write.csv2(tabla_genes_inflamatorios, file.path(config$carpeta_tablas, "Genes_inflamatorios.csv"), row.names = FALSE)
write.csv2(tabla_frecuencia_genes, file.path(config$carpeta_tablas, "Frecuencia_genes.csv"), row.names = FALSE)

obtener_subtitulo <- function(nombre_comp, texto_normal) {
  if (nombre_comp == "Cronico_vs_Agudo") {
    return("Comparación entre fase aguda y fase crónica")
  }
  texto_normal
}

guardar_dotplot_manual <- function(tabla, nombre_comp, base) {
  if (nrow(tabla) == 0) return(NULL)
  df <- tabla %>%
    filter(comparacion == nombre_comp) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste(base, "-", nombre_comp),
      subtitle = obtener_subtitulo(nombre_comp, "GSEA usando todos los genes ordenados por stat"),
      x = "NES",
      y = "Proceso / ruta",
      size = "Tamaño",
      color = "p ajustado"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(config$carpeta_figuras, paste0("dotplot_", base, "_", nombre_comp, ".png")),
    plot = p, width = 11, height = 8, dpi = 300
  )
}

guardar_barplot_genes <- function(nombre_comp) {
  if (nrow(tabla_frecuencia_genes) == 0) return(NULL)
  df <- tabla_frecuencia_genes %>%
    filter(comparacion == nombre_comp) %>%
    arrange(desc(n_procesos_total)) %>%
    slice_head(n = 20)
  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = reorder(Gen, n_procesos_total), y = n_procesos_total)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Genes presentes en mayor número de rutas -", nombre_comp),
      subtitle = obtener_subtitulo(nombre_comp, NULL),
      x = "Gen",
      y = "Número de procesos/rutas"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(config$carpeta_figuras, paste0("barplot_genes_", nombre_comp, ".png")),
    plot = p, width = 10, height = 7, dpi = 300
  )
}

guardar_volcano <- function(nombre_comp) {
  if (nrow(tabla_deseq) == 0) return(NULL)
  df <- tabla_deseq %>%
    filter(comparacion == nombre_comp,
           !is.na(padj), !is.na(log2FoldChange))
  if (nrow(df) == 0) return(NULL)

  n_etiquetas <- 15

  # ---- Función interna que arma el gráfico dado un set de etiquetas ----
  construir_volcano <- function(etiquetas, sufijo_titulo, sufijo_archivo) {
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = cambio)) +
      geom_point(alpha = 0.6, size = 1.4) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      ggrepel::geom_text_repel(
        data = etiquetas, aes(label = simbolo),
        size = 3, max.overlaps = 30, show.legend = FALSE
      ) +
      scale_color_manual(
        values = c("Sube" = "#D6604D", "Baja" = "#4393C3", "No significativo" = "grey80"),
        name = "Cambio"
      ) +
      labs(
        title = paste0("Volcano plot - ", nombre_comp, sufijo_titulo),
        subtitle = obtener_subtitulo(
          nombre_comp,
          paste0("Sube: ", sum(df$cambio == "Sube"), " | Baja: ", sum(df$cambio == "Baja"))
        ),
        x = "log2 Fold Change", y = "-log10(padj)"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")

    ggsave(
      filename = file.path(config$carpeta_figuras,
                            paste0("volcano_", nombre_comp, sufijo_archivo, ".png")),
      plot = p, width = 9, height = 7, dpi = 300
    )
  }

  # ---- Versión 1: solo genes de literatura (la que ya tenías) ----
  etiquetas_literatura <- df %>%
    filter(cambio %in% c("Sube", "Baja"), inflamatorio == "Sí") %>%
    arrange(padj) %>%
    group_by(cambio) %>%
    slice_head(n = n_etiquetas) %>%
    ungroup()

  construir_volcano(
    etiquetas_literatura,
    sufijo_titulo  = " (genes de literatura)",
    sufijo_archivo = "_literatura"
  )

  # ---- Versión 2: los más significativos, sin filtrar por literatura ----
  etiquetas_significativos <- df %>%
    filter(cambio %in% c("Sube", "Baja")) %>%
    arrange(padj) %>%
    group_by(cambio) %>%
    slice_head(n = n_etiquetas) %>%
    ungroup()

  construir_volcano(
    etiquetas_significativos,
    sufijo_titulo  = " (genes más significativos)",
    sufijo_archivo = "_significativos"
  )
}

guardar_figuras_enrichplot <- function(nombre_comp) {
  go_obj <- resultados_analisis[[nombre_comp]]$objetos_gsea$go
  kegg_obj <- resultados_analisis[[nombre_comp]]$objetos_gsea$kegg

  tryCatch({
    if (nrow(as.data.frame(go_obj)) > 0) {
      png(file.path(config$carpeta_figuras, paste0("cnetplot_GO_", nombre_comp, ".png")), width = 3000, height = 2200, res = 300)
      p_cnet_go <- cnetplot(go_obj, showCategory = 8, foldChange = NULL, circular = FALSE) +
        labs(subtitle = obtener_subtitulo(nombre_comp, NULL))
      print(p_cnet_go)
      dev.off()

      png(file.path(config$carpeta_figuras, paste0("emapplot_GO_", nombre_comp, ".png")), width = 3000, height = 2200, res = 300)
      p_emap_go <- emapplot(pairwise_termsim(go_obj), showCategory = 15) +
        labs(subtitle = obtener_subtitulo(nombre_comp, NULL))
      print(p_emap_go)
      dev.off()
    }
  }, error = function(e) {
    message("No se pudo generar cnetplot/emapplot GO para ", nombre_comp, ": ", e$message)
    if (dev.cur() != 1) dev.off()
  })

  tryCatch({
    if (nrow(as.data.frame(kegg_obj)) > 0) {
      png(file.path(config$carpeta_figuras, paste0("cnetplot_KEGG_", nombre_comp, ".png")), width = 3000, height = 2200, res = 300)
      p_cnet_kegg <- cnetplot(kegg_obj, showCategory = 8, foldChange = NULL, circular = FALSE) +
        labs(subtitle = obtener_subtitulo(nombre_comp, NULL))
      print(p_cnet_kegg)
      dev.off()

      png(file.path(config$carpeta_figuras, paste0("emapplot_KEGG_", nombre_comp, ".png")), width = 3000, height = 2200, res = 300)
      p_emap_kegg <- emapplot(pairwise_termsim(kegg_obj), showCategory = 15) +
        labs(subtitle = obtener_subtitulo(nombre_comp, NULL))
      print(p_emap_kegg)
      dev.off()
    }
  }, error = function(e) {
    message("No se pudo generar cnetplot/emapplot KEGG para ", nombre_comp, ": ", e$message)
    if (dev.cur() != 1) dev.off()
  })
}

for (nombre_comp in names(resultados_analisis)) {
  guardar_dotplot_manual(tabla_go, nombre_comp, "GO")
  guardar_dotplot_manual(tabla_kegg, nombre_comp, "KEGG")
  guardar_barplot_genes(nombre_comp)
  guardar_volcano(nombre_comp)          # <- línea nueva
  guardar_figuras_enrichplot(nombre_comp)
}

cat("\nExportación finalizada correctamente.\n")
cat("Excel generado en:", ruta_excel, "\n")
cat("Tablas CSV guardadas en:", config$carpeta_tablas, "\n")
cat("Figuras guardadas en:", config$carpeta_figuras, "\n")


