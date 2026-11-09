# esclerosis-rnaseq

Transcriptómica de la remielinización post-crónica en el modelo de cuprizona (esclerosis múltiple progresiva).
**Pregunta:** ¿el microambiente lesional crónico se mantiene estático o evoluciona tras retirar la cuprizona?
Foco en el déficit mitocondrial persistente y en *Hes5* como freno de la diferenciación de OPCs.

**Responsable:** Alelí Galaz Meléndez · **Institución:** USACH

## Flujo de ejecución

| Paso | Script(s) | Qué hace |
|---|---|---|
| 1. Descarga | `descargar_sra.sh` | Descarga las lecturas de SRA |
| 2. Control de calidad | `correr_fastqc.sh`, `fastp_db*.sh`, `correr_trimming.sh`, `organizar_trimmed.sh`, `MultiQC.sh` | QC y recorte con fastp; reportes MultiQC |
| 3. Alineamiento | `alinear_single_end.sh`, `alinear_paired_end_db9.sh` | Alineamiento con STAR |
| 4. Conteo | `contar_genes.sh` | Matrices de conteo (`Conteos/`) |
| 5. Preparación | `preparar_deseq2.R` | Une conteos + metadata, filtra, crea `dds_preparado.rds` |
| 6. Normalización | `Normalizacion/normalizar.R`, `Normalizacion/corregido_batch.R` | VST, PCA, heatmaps; guarda `dds_normalizado.rds` |
| 7. Análisis post-crónico | `post-cronico/main.R` (00–05) | DESeq2 + GSEA GO/KEGG, 6 comparaciones |
| 8. Extensión | `post-cronico/main_extendido.R` (06–10) | Comparaciones por tejido y 24 sem solo hembras |
| 9. Figuras del póster | `post-cronico/figuras_poster/figuras_poster.R` | Figuras 1–6 y S1 (PNG 300 dpi + PDF) |
| Reportes | `reporte_comparativo.sh`, `reporte_pesos.sh` | Resúmenes de tamaños y QC |

Pasos 7–9 en un solo comando (desde la raíz):

```bash
cp <ruta>/dds_normalizado.rds post-cronico/datos/
bash correr_postcronico.sh            # análisis + figuras
bash correr_postcronico.sh figuras    # solo figuras
```

`Analisis_KEGG/` es la versión anterior del mismo pipeline (00–05 idénticos); la versión oficial es `post-cronico/`.

## Decisiones de análisis

- DEG: padj < 0,05 y |log2FC| > 1. GSEA: ranking por estadístico de Wald, GO BP y KEGG, minGSSize 5, maxGSSize 1000, semilla 123.
- Diseño: se usa el modelo más completo con rango completo entre `~ dataset + tejido + grupo`, `~ tejido + grupo`, `~ dataset + grupo` y `~ grupo`. El diseño usado en cada comparación queda en `resultados/tablas/disenos_usados.csv`.
- Sexo: Postcrónico 24 sem tiene solo hembras (9H/0M) frente a Control 24 sem (10H/6M); por eso se corre un análisis de sensibilidad solo con hembras (`datos/metadata_sexo.csv`, tomado de `Metadata/metadata5.csv`).

## Datos

Los datos crudos, intermedios e índices no están en Git. Ver [DATOS.md](DATOS.md).
