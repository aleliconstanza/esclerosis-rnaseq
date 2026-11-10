# ============================================================
# figuras_poster.R
# Figuras del póster "Déficit mitocondrial persistente y debilitamiento
# del bloqueo diferenciativo de OPCs durante la remielinización
# post-crónica inducida por cuprizona".
#
# Uso (desde 01_scripts/09_figuras_poster/):
#   Rscript figuras_poster.R
#   Rscript figuras_poster/figuras_poster.R <tablas> <metadata_muestras> <metadata_sexo> <salida>
#
# Lee las tablas que exporta el pipeline (05 y 10) y guarda cada
# figura en PNG (300 dpi) y PDF, con el tamaño del espacio en el póster.
# Solo necesita ggplot2, dplyr, tidyr, stringr y scales.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
ruta_tablas   <- if (length(args) >= 1) args[1] else "../../03_resultados/08_postcronico/tablas"
ruta_metadata <- if (length(args) >= 2) args[2] else "../../02_intermedios/06_normalizacion/metadata_muestras.csv"
ruta_sexo     <- if (length(args) >= 3) args[3] else "../00_config/metadata/metadata_sexo.csv"
ruta_salida   <- if (length(args) >= 4) args[4] else "../../03_resultados/09_figuras_poster"
ruta_categorias <- "../00_config/genes_literatura/categorias_genes_curados.csv"
# Si aún no hay corrida definitiva del post-crónico, usa la corrida anterior (07_analisis_global)
if (length(args) == 0 && !file.exists(file.path(ruta_tablas, "KEGG.csv"))) {
  ruta_tablas <- "../../03_resultados/07_analisis_global/tablas"
  cat("AVISO: no hay tablas en 03_resultados/08_postcronico; uso", ruta_tablas, "\n")
}
dir.create(ruta_salida, showWarnings = FALSE, recursive = TRUE)

# ---------------- Estilo común ----------------
NARANJA <- "#EA580C"; OXIDO <- "#7C2D12"; AZUL <- "#1D4E89"; GRIS <- "#6B7280"
PADJ <- 0.05

tema <- theme_minimal(base_size = 24, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.text = element_text(color = "grey15"),
    axis.title = element_text(color = OXIDO, face = "bold"),
    strip.text = element_text(color = OXIDO, face = "bold", hjust = 0),
    legend.title = element_text(size = 18, color = "grey20"),
    legend.text = element_text(size = 16),
    plot.title = element_text(color = OXIDO, face = "bold", size = 26),
    plot.margin = margin(10, 15, 10, 10)
  )
escala_nes <- function(lim) scale_fill_gradient2(low = AZUL, mid = "white", high = NARANJA,
                                                 midpoint = 0, limits = c(-lim, lim), oob = squish, name = "NES")

guardar <- function(p, nombre, ancho, alto) {
  ggsave(file.path(ruta_salida, paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = 300, bg = "white")
  ggsave(file.path(ruta_salida, paste0(nombre, ".pdf")), p, width = ancho, height = alto, bg = "white",
         device = if (capabilities("cairo")) cairo_pdf else "pdf")
  cat("  ✔", nombre, "\n")
}

# Comparaciones y etiquetas
VENTANAS <- c(
  Cronico_vs_Control           = "12 sem\n(crónico)",
  Postcronico16w_vs_Control24w = "16 sem",
  Postcronico20w_vs_Control24w = "20 sem",
  Postcronico24w_vs_Control24w = "24 sem"
)
REMIE <- names(VENTANAS)[2:4]

leer <- function(archivo) {
  ruta <- file.path(ruta_tablas, archivo)
  if (!file.exists(ruta)) stop("No se encontró ", ruta)
  read.csv2(ruta, stringsAsFactors = FALSE, check.names = FALSE)
}

cat("\nLeyendo tablas desde", ruta_tablas, "...\n")
kegg <- leer("KEGG.csv"); go <- leer("GO.csv")
gsea <- bind_rows(kegg, go) %>% filter(comparacion %in% names(VENTANAS))
deseq_cols <- c("comparacion", "simbolo", "inflamatorio", "log2FoldChange", "lfcSE", "padj")
deseq <- leer("DESeq2.csv")[, deseq_cols] %>% filter(comparacion %in% names(VENTANAS))
resumen <- leer("Resumen.csv")

etiquetar <- function(df) df %>%
  mutate(ventana = factor(VENTANAS[comparacion], levels = VENTANAS))

dotplot_rutas <- function(df, terminos, titulo_y = NULL, lim = 3) {
  datos <- df %>%
    inner_join(terminos, by = c("base", "Description")) %>%
    etiquetar() %>%
    mutate(sig = p.adjust < PADJ,
           mlog = pmin(-log10(p.adjust), 10),
           etiqueta = factor(etiqueta, levels = rev(unique(terminos$etiqueta))))
  ggplot(datos, aes(ventana, etiqueta)) +
    geom_point(aes(size = mlog, fill = NES, shape = sig), color = "grey25", stroke = 0.8) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 4), labels = c(`TRUE` = "padj < 0,05", `FALSE` = "n.s."),
                       name = NULL, drop = FALSE) +
    scale_size_continuous(range = c(3, 15), limits = c(0, 10), name = expression(-log[10]~padj), breaks = c(2, 5, 10)) +
    escala_nes(lim) +
    labs(x = NULL, y = titulo_y) + tema +
    guides(shape = guide_legend(order = 3, override.aes = list(size = 7, fill = NARANJA)),
           size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey70")),
           fill = guide_colorbar(order = 1, barheight = unit(4, "cm")))
}

# ============================================================
# Fig. 1 — Diseño experimental
# ============================================================
cat("\nFig. 1: diseño experimental\n")
if (file.exists(ruta_metadata)) {
  meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
  sexo <- if (file.exists(ruta_sexo)) read.csv(ruta_sexo, stringsAsFactors = FALSE) else data.frame(muestra = character(), sexo = character())
  meta <- meta %>%
    filter(dataset != "db8") %>%
    mutate(g = str_remove(tolower(grupo), "(cx|cc|hp)$"),
           punto = case_when(g == "des12w" ~ "Crónico", g == "rem16w" ~ "16", g == "rem20w" ~ "20",
                             g == "rem24w" ~ "24", g == "control24w" ~ "Control 24", TRUE ~ NA_character_)) %>%
    filter(!is.na(punto)) %>% left_join(sexo, by = "muestra")
  n_tab <- meta %>% group_by(punto) %>%
    summarise(n = n(), H = sum(sexo == "F", na.rm = TRUE), M = sum(sexo == "M", na.rm = TRUE), .groups = "drop") %>%
    mutate(lab = ifelse(H + M > 0, sprintf("n = %d\n(%dH / %dM)", n, H, M), sprintf("n = %d", n)))
  lab_n <- setNames(n_tab$lab, n_tab$punto)

  puntos <- data.frame(
    x = c(12, 16, 20, 24), y = 1,
    nombre = c("Crónico", "16 sem", "20 sem", "24 sem"),
    lab = lab_n[c("Crónico", "16", "20", "24")]
  )
  p1 <- ggplot() +
    annotate("rect", xmin = 0, xmax = 12, ymin = 0.85, ymax = 1.15, fill = NARANJA, alpha = 0.9) +
    annotate("rect", xmin = 12, xmax = 24.6, ymin = 0.85, ymax = 1.15, fill = "#FED7AA") +
    annotate("text", x = 6, y = 1, label = "Cuprizona 0,2 %", color = "white", fontface = "bold", size = 9) +
    annotate("text", x = 12, y = 1.55, label = "Retiro de cuprizona", color = OXIDO, fontface = "bold", size = 7) +
    annotate("segment", x = 12, xend = 12, y = 1.45, yend = 1.18, color = OXIDO, linewidth = 1.2,
             arrow = arrow(length = unit(0.25, "cm"))) +
    geom_point(data = puntos, aes(x, y), shape = 21, size = 11, fill = "white", color = OXIDO, stroke = 2.5) +
    geom_text(data = puntos, aes(x, 0.55, label = nombre), fontface = "bold", color = OXIDO, size = 7.5) +
    geom_text(data = puntos, aes(x, 0.22, label = lab), color = "grey20", size = 6, lineheight = 0.9) +
    annotate("label", x = 24, y = 1.6, label = paste0("Control 24 sem\n", lab_n["Control 24"]),
             fill = "grey95", color = "grey20", size = 6, label.size = 0, lineheight = 0.9) +
    scale_x_continuous(breaks = seq(0, 24, 4), labels = function(b) paste(b, "sem"), limits = c(-0.5, 26)) +
    scale_y_continuous(limits = c(0, 1.85)) +
    labs(x = NULL, y = NULL) + tema +
    theme(axis.text.y = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(size = 18))
  guardar(p1, "Fig1_diseno_experimental", 16.8, 3.9)
} else {
  cat("  (omitida: no se encontró", ruta_metadata, ")\n")
}

# ============================================================
# Fig. 2 — Rutas mitocondriales a lo largo del tiempo (figura central)
# ============================================================
cat("Fig. 2: rutas mitocondriales en el tiempo\n")
terminos_mito <- tribble(
  ~base, ~Description, ~etiqueta, ~eje,
  "KEGG", "Oxidative phosphorylation", "Oxidative phosphorylation", "Bioenergética\nmitocondrial",
  "GO", "respiratory electron transport chain", "Respiratory electron transport chain", "Bioenergética\nmitocondrial",
  "GO", "mitochondrial ATP synthesis coupled electron transport", "Mitochondrial ATP synthesis", "Bioenergética\nmitocondrial",
  "KEGG", "Thermogenesis", "Thermogenesis", "Bioenergética\nmitocondrial",
  "KEGG", "Ribosome", "Ribosome", "Bioenergética\nmitocondrial",
  "KEGG", "Chemical carcinogenesis - reactive oxygen species", "Reactive oxygen species", "Bioenergética\nmitocondrial",
  "KEGG", "Pathways of neurodegeneration - multiple diseases", "Neurodegeneration (multiple)", "Neurodegene-\nración",
  "KEGG", "Parkinson disease", "Parkinson disease", "Neurodegene-\nración",
  "KEGG", "Huntington disease", "Huntington disease", "Neurodegene-\nración",
  "KEGG", "Alzheimer disease", "Alzheimer disease", "Neurodegene-\nración",
  "KEGG", "Amyotrophic lateral sclerosis", "Amyotrophic lateral sclerosis", "Neurodegene-\nración",
  "KEGG", "Prion disease", "Prion disease", "Neurodegene-\nración"
)
p2 <- dotplot_rutas(gsea, terminos_mito) +
  facet_grid(eje ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = OXIDO, linewidth = 1) +
  geom_text(data = data.frame(eje = terminos_mito$eje[1]), aes(x = 1.55, y = Inf), label = "  retiro →",
            inherit.aes = FALSE, hjust = 0, vjust = 1.4, color = OXIDO, size = 6.5, fontface = "italic") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 18))
guardar(p2, "Fig2_rutas_mitocondriales_tiempo", 21, 6.2)

# ============================================================
# Fig. 3 — DEGs por ventana
# ============================================================
cat("Fig. 3: DEGs por ventana\n")
degs <- resumen %>%
  filter(comparacion %in% REMIE) %>%
  transmute(comparacion, Suben = as.numeric(genes_sube), Bajan = -as.numeric(genes_baja)) %>%
  pivot_longer(c(Suben, Bajan), names_to = "dir", values_to = "n") %>%
  etiquetar() %>% mutate(ventana = droplevels(ventana))
lim3 <- max(abs(degs$n)) * 1.25
p3 <- ggplot(degs, aes(ventana, n, fill = dir)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_text(aes(label = abs(n), vjust = ifelse(n > 0, -0.4, 1.3)), size = 8, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = c(Suben = NARANJA, Bajan = AZUL), breaks = c("Suben", "Bajan"), name = NULL) +
  scale_y_continuous(labels = function(b) abs(b), limits = c(-lim3, lim3)) +
  labs(x = NULL, y = "N.º de DEGs") + tema +
  theme(legend.position = "top")
guardar(p3, "Fig3_DEGs_por_ventana", 12.6, 6.2)

# ============================================================
# Fig. 4 — Rutas positivas específicas por etapa
# ============================================================
cat("Fig. 4: rutas positivas por etapa\n")
terminos_pos <- tribble(
  ~base, ~Description, ~etiqueta,
  "KEGG", "Integrin signaling", "Integrin signaling",
  "GO", "actin cytoskeleton organization", "Actin cytoskeleton org.",
  "KEGG", "Autophagy - animal", "Autophagy (KEGG)",
  "GO", "autophagosome organization", "Autophagosome org.",
  "GO", "microtubule cytoskeleton organization", "Microtubule cytoskeleton org.",
  "GO", "extracellular matrix organization", "Extracellular matrix org."
)
p4 <- dotplot_rutas(gsea %>% filter(comparacion %in% REMIE), terminos_pos, lim = 2.5) +
  theme(legend.position = "bottom", legend.box = "horizontal", axis.text.y = element_text(size = 20),
        legend.title = element_text(size = 14), legend.text = element_text(size = 13)) +
  guides(fill = guide_colorbar(order = 1, barwidth = unit(4, "cm"), barheight = unit(0.5, "cm")),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey70")),
         shape = "none")
guardar(p4, "Fig4_rutas_positivas_por_etapa", 11, 4.8)

# ============================================================
# Fig. 5 — Genes curados (solo los significativos en ≥1 ventana)
# ============================================================
cat("Fig. 5: genes curados\n")
categorias <- read.csv(ruta_categorias, stringsAsFactors = FALSE)
curados <- deseq %>%
  filter(comparacion %in% REMIE, inflamatorio == "Sí") %>%
  left_join(categorias, by = "simbolo") %>%
  group_by(simbolo) %>% filter(any(padj < PADJ, na.rm = TRUE)) %>% ungroup() %>%
  etiquetar() %>% mutate(ventana = droplevels(ventana),
                         estrella = ifelse(!is.na(padj) & padj < PADJ, "*", ""),
                         categoria = recode(coalesce(categoria, "Otros"),
                                            "Matriz extracelular" = "MEC", "Astrogliosis" = "Astrogliosis",
                                            "Linaje oligodendroglial" = "OPC / OL",
                                            "Inflamación y microglía" = "Inflam.", "Neuronal / otros" = "Otros"),
                         categoria = factor(categoria, levels = c("MEC", "Astrogliosis", "OPC / OL", "Inflam.", "Otros")))
cat("  genes curados significativos:", n_distinct(curados$simbolo), "\n")
lim5 <- max(1, ceiling(max(abs(curados$log2FoldChange), na.rm = TRUE)))
curados <- curados %>% mutate(ventana_r = factor(ventana, levels = rev(levels(ventana))))
p5 <- ggplot(curados, aes(simbolo, ventana_r)) +
  geom_tile(aes(fill = log2FoldChange), color = "white", linewidth = 1) +
  geom_text(aes(label = estrella), size = 9, vjust = 0.75, color = "grey10") +
  facet_grid(. ~ categoria, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = AZUL, mid = "white", high = NARANJA, midpoint = 0,
                       limits = c(-lim5, lim5), oob = squish, name = "log2FC") +
  labs(x = NULL, y = NULL) + tema +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, face = "italic", size = 18),
        strip.text = element_text(size = 15, hjust = 0.5), panel.grid = element_blank(),
        legend.position = "right")
guardar(p5, "Fig5_genes_curados", 11, 4.8)

# ============================================================
# Fig. 6 — Hes5 frente a genes de mielina (recuperación disociada)
#   Mbp/Plp1/Mog caen en la fase crónica y vuelven al nivel control desde
#   16 sem; Hes5 solo baja a 24 sem. Punto relleno = padj < 0,05.
#   Si existe el análisis solo hembras, se agrega como rombo azul en 24 sem.
# ============================================================
cat("Fig. 6: Hes5 y genes de mielina\n")
genes6 <- c("Hes5", "Mbp", "Plp1", "Mog")
col6 <- c(Hes5 = NARANJA, Mbp = "#5B7FA6", Plp1 = "#2F4F7F", Mog = "#8FA9C8")
d6 <- deseq %>% filter(simbolo %in% genes6) %>% etiquetar() %>%
  mutate(x = as.numeric(ventana), sig = !is.na(padj) & padj < PADJ,
         simbolo = factor(simbolo, levels = genes6))

hembras <- NULL
ruta_ext <- file.path(ruta_tablas, "DESeq2_extendido.csv")
if (file.exists(ruta_ext)) {
  hembras <- read.csv2(ruta_ext, stringsAsFactors = FALSE)[, deseq_cols] %>%
    filter(simbolo == "Hes5", comparacion == "Postcronico24w_hembras_vs_Control24w_hembras") %>%
    mutate(x = 4.22, sig = !is.na(padj) & padj < PADJ)
  if (nrow(hembras) == 0) hembras <- NULL
}

# Etiquetas al final de cada línea, separadas para que no se encimen
etiq6 <- d6 %>% filter(x == max(x)) %>% arrange(log2FoldChange) %>%
  mutate(y_lab = log2FoldChange)
sep <- 0.17
for (i in seq_len(nrow(etiq6))[-1]) etiq6$y_lab[i] <- max(etiq6$y_lab[i], etiq6$y_lab[i - 1] + sep)
p6 <- ggplot(d6, aes(x, log2FoldChange, color = simbolo, group = simbolo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = GRIS) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = OXIDO) +
  annotate("text", x = 1.55, y = Inf, label = " retiro →", hjust = 0, vjust = 1.3,
           color = OXIDO, size = 5.5, fontface = "italic") +
  geom_line(data = filter(d6, x >= 2), aes(linewidth = simbolo == "Hes5")) +
  geom_errorbar(data = filter(d6, simbolo == "Hes5"),
                aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.08, linewidth = 1) +
  geom_point(aes(fill = ifelse(sig, as.character(simbolo), "ns")), shape = 21, size = 5.5, stroke = 1.6) +
  geom_text(data = etiq6, aes(x = x + 0.12, y = y_lab, label = simbolo), hjust = 0, size = 6.5,
            fontface = "bold.italic", show.legend = FALSE) +
  scale_linewidth_manual(values = c(`TRUE` = 2, `FALSE` = 1), guide = "none") +
  scale_color_manual(values = col6, guide = "none") +
  scale_fill_manual(values = c(col6, ns = "white"), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = levels(d6$ventana), expand = expansion(add = c(0.3, 0.7))) +
  labs(x = NULL, y = expression(bold(log[2]*FC~"vs control"))) + tema +
  theme(axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 18))
if (!is.null(hembras)) {
  p6 <- p6 +
    geom_errorbar(data = hembras, aes(x = x, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
                  inherit.aes = FALSE, width = 0.06, color = AZUL, linewidth = 1) +
    geom_point(data = hembras, aes(x = x, y = log2FoldChange), inherit.aes = FALSE,
               shape = 23, size = 5.5, stroke = 1.6, color = AZUL, fill = ifelse(hembras$sig, AZUL, "white"))
}
guardar(p6, "Fig6_Hes5_y_mielina", 11, 4.8)

# ============================================================
# Fig. S1 — Eje Tgfb1–Jag1–Notch1–Hes5 y marcadores de mielinización (apoyo)
# ============================================================
cat("Fig. S1: eje Notch y marcadores de mielinización\n")
eje <- c("Tgfb1", "Jag1", "Notch1", "Hes1", "Hes5", "Id2", "Id4", "Pdgfra", "Gpr17",
         "Olig2", "Sox10", "Myrf", "Mbp", "Plp1", "Mog")
s1 <- deseq %>% filter(simbolo %in% eje) %>% etiquetar() %>%
  mutate(simbolo = factor(simbolo, levels = rev(eje)),
         estrella = ifelse(!is.na(padj) & padj < PADJ, "*", ""),
         bloque = ifelse(simbolo %in% c("Tgfb1", "Jag1", "Notch1", "Hes1", "Hes5", "Id2", "Id4"),
                         "Señalización inhibitoria", "Diferenciación y mielina"),
         bloque = factor(bloque, levels = c("Señalización inhibitoria", "Diferenciación y mielina")))
pS1 <- ggplot(s1, aes(ventana, simbolo)) +
  geom_tile(aes(fill = log2FoldChange), color = "white", linewidth = 1) +
  geom_text(aes(label = estrella), size = 9, vjust = 0.75) +
  facet_grid(bloque ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = AZUL, mid = "white", high = NARANJA, midpoint = 0, limits = c(-2, 2),
                       oob = squish, name = "log2FC") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = OXIDO) +
  labs(x = NULL, y = NULL) + tema +
  theme(axis.text.y = element_text(face = "italic"), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, size = 16))
guardar(pS1, "FigS1_eje_Notch_mielina", 11, 8)

cat("\nFiguras guardadas en", ruta_salida, "\n")
