# ============================================================
# 10_integracion_exportar_extendido.R
# Integra y exporta (Excel + CSV + figuras) los resultados de las
# comparaciones extendidas (09_analisis_extendido.R). Mismo formato
# que 04_integracion.R + 05_exportar.R, pero en archivos separados
# para no sobreescribir el analisis oficial (16w/20w/24w combinados).
# ============================================================

cat("\nIntegrando resultados extendidos...\n")

tabla_resumen <- bind_rows(lapply(resultados_analisis_ext, function(x) x$resumen))
tabla_deseq <- bind_rows(lapply(resultados_analisis_ext, function(x) x$deseq))
tabla_go <- bind_rows(lapply(resultados_analisis_ext, function(x) x$go))
tabla_kegg <- bind_rows(lapply(resultados_analisis_ext, function(x) x$kegg))
tabla_leading_edge <- bind_rows(lapply(resultados_analisis_ext, function(x) x$leading_edge))
tabla_genes_en_procesos <- bind_rows(lapply(resultados_analisis_ext, function(x) x$genes_en_procesos))

tabla_genes_inflamatorios <- if (nrow(tabla_genes_en_procesos) > 0) {
  tabla_genes_en_procesos %>%
    filter(inflamatorio == "Sí") %>%
    arrange(desc(n_total_procesos), desc(abs(log2FC)))
} else {
  data.frame()
}

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

# --- Excel + CSV (mismo patron que 05_exportar.R) ---
ruta_excel_ext <- file.path(config$carpeta_excel, "Analisis_GSEA_KEGG_extendido.xlsx")

wb <- createWorkbook()
agregar_hoja <- function(wb, nombre, tabla) {
  addWorksheet(wb, nombre)
  if (ncol(tabla) == 0 || nrow(tabla) == 0) {
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

saveWorkbook(wb, ruta_excel_ext, overwrite = TRUE)

write.csv2(tabla_resumen, file.path(config$carpeta_tablas, "Resumen_extendido.csv"), row.names = FALSE)
write.csv2(tabla_deseq, file.path(config$carpeta_tablas, "DESeq2_extendido.csv"), row.names = FALSE)
write.csv2(tabla_go, file.path(config$carpeta_tablas, "GO_extendido.csv"), row.names = FALSE)
write.csv2(tabla_kegg, file.path(config$carpeta_tablas, "KEGG_extendido.csv"), row.names = FALSE)
write.csv2(tabla_frecuencia_genes, file.path(config$carpeta_tablas, "Frecuencia_genes_extendido.csv"), row.names = FALSE)

cat("\nExcel extendido guardado en:", ruta_excel_ext, "\n")

# --- Figuras: volcano, dotplot, barplot (mismas funciones/estetica que 05_exportar.R) ---

guardar_dotplot_manual <- function(tabla, nombre_comp, base) {
  if (nrow(tabla) == 0) return(NULL)
  df <- tabla %>% filter(comparacion == nombre_comp) %>% arrange(p.adjust) %>% slice_head(n = 20)
  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste(base, "-", nombre_comp),
      subtitle = "GSEA usando todos los genes ordenados por stat",
      x = "NES", y = "Proceso / ruta", size = "Tamano", color = "p ajustado"
    ) +
    theme_bw()

  ggsave(file.path(config$carpeta_figuras, paste0("dotplot_", base, "_", nombre_comp, ".png")),
         p, width = 11, height = 8, dpi = 300)
}

guardar_barplot_genes <- function(nombre_comp) {
  if (nrow(tabla_frecuencia_genes) == 0) return(NULL)
  df <- tabla_frecuencia_genes %>% filter(comparacion == nombre_comp) %>%
    arrange(desc(n_procesos_total)) %>% slice_head(n = 20)
  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = reorder(Gen, n_procesos_total), y = n_procesos_total)) +
    geom_col() + coord_flip() +
    labs(
      title = paste("Genes presentes en mayor numero de rutas -", nombre_comp),
      x = "Gen", y = "Numero de procesos/rutas"
    ) +
    theme_bw()

  ggsave(file.path(config$carpeta_figuras, paste0("barplot_genes_", nombre_comp, ".png")),
         p, width = 10, height = 7, dpi = 300)
}

guardar_volcano <- function(nombre_comp) {
  if (nrow(tabla_deseq) == 0) return(NULL)
  df <- tabla_deseq %>% filter(comparacion == nombre_comp, !is.na(padj), !is.na(log2FoldChange))
  if (nrow(df) == 0) return(NULL)

  n_etiquetas <- 15

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
        subtitle = paste0("Sube: ", sum(df$cambio == "Sube"), " | Baja: ", sum(df$cambio == "Baja")),
        x = "log2 Fold Change", y = "-log10(padj)"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")

    ggsave(
      file.path(config$carpeta_figuras, paste0("volcano_", nombre_comp, sufijo_archivo, ".png")),
      p, width = 9, height = 7, dpi = 300
    )
  }

  etiquetas_literatura <- df %>%
    filter(cambio %in% c("Sube", "Baja"), inflamatorio == "Sí") %>%
    arrange(padj) %>% group_by(cambio) %>% slice_head(n = n_etiquetas) %>% ungroup()
  construir_volcano(etiquetas_literatura, " (genes de literatura)", "_literatura")

  etiquetas_significativos <- df %>%
    filter(cambio %in% c("Sube", "Baja")) %>%
    arrange(padj) %>% group_by(cambio) %>% slice_head(n = n_etiquetas) %>% ungroup()
  construir_volcano(etiquetas_significativos, " (genes mas significativos)", "_significativos")
}

for (nombre_comp in names(resultados_analisis_ext)) {
  guardar_dotplot_manual(tabla_go, nombre_comp, "GO")
  guardar_dotplot_manual(tabla_kegg, nombre_comp, "KEGG")
  guardar_barplot_genes(nombre_comp)
  guardar_volcano(nombre_comp)
}

cat("\nFiguras extendidas (volcano, dotplot, barplot) guardadas en:", config$carpeta_figuras, "\n")
cat("Compara visualmente Corteza/Hipocampo/CuerpoCalloso por ventana temporal,\n")
cat("y la version pareada por sexo (24w solo hembras) contra el resultado oficial.\n")
