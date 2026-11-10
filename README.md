# esclerosis-rnaseq

Transcriptómica de la remielinización post-crónica en el modelo de cuprizona (esclerosis múltiple progresiva).
Responsable: Alelí Galaz Meléndez · USACH

## Estructura

| Carpeta | Contenido |
|---|---|
| `scripts/00_config/` | Metadata, lista de genes curados y `decisiones.md` |
| `scripts/01_…09_/` | Código de cada etapa, numerado en orden de ejecución |
| `intermedios/01_…08_/` | Salidas intermedias: reportes fastp y MultiQC, logs de STAR, conteos, anotación, objetos de normalización |
| `resultados/06_…09_/` | Tablas y figuras finales |
| `disco/` | Acceso directo a `/data/Database` (FASTQ, BAM, índices, `.rds`; no está en Git) |

| Etapa | Scripts | Intermedios | Resultados |
|---|---|---|---|
| 01 Descarga | descarga SRA, pesos | reporte de pesos | |
| 02 Control de calidad | fastp, MultiQC | fastp/, multiqc/ | |
| 03 Alineamiento | STAR SE y PE | star_logs/ | |
| 04 Conteo | featureCounts | matrices y resúmenes | |
| 05 Anotación | anotación y genes curados | genes anotados | |
| 06 Normalización | preparar, normalizar, batch | metadata, size factors, `.rds` (enlaces) | PCA y heatmaps |
| 07 Análisis global | DESeq2 + GSEA, 6 comparaciones | | tablas y figuras |
| 08 Post-crónico | DESeq2 + GSEA + tejido y sexo | logs | tablas y figuras |
| 09 Figuras del póster | figuras_poster.R | | Fig. 1–6 y S1 |

## Cómo correr

```bash
bash scripts/correr_postcronico.sh            # análisis post-crónico + figuras
bash scripts/correr_postcronico.sh figuras    # solo figuras
```

## Pendiente
- Los scripts de las etapas 01–06 todavía usan rutas absolutas de `/data/Database`. Para reproducibilidad completa hay que pasarlos a rutas relativas (`disco/`, `intermedios/`).
