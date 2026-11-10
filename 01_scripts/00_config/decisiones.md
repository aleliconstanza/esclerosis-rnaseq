# Decisiones del análisis

Registro de qué se decidió y por qué. Actualizar cada vez que cambie algo.

## Datos
- Datasets: db5 (GSE212147, corteza), db5_cchp (GSE212147, cuerpo calloso + hipocampo), db8 (GSE124266), db9 (GSE136863, paired-end), db11 (GSE278137).
- db8 excluido del análisis conjunto: duplicación ≥ 80 %, baja asignación de lecturas, mezcla microglía / no-microglía y fracaso de la corrección de batch.
- db9: se descarta 5.5w. db8/db9/db11: solo WT.

## Muestras excluidas (control de calidad)
| Muestra | Motivo |
|---|---|
| SRR21224515 | %PF 83,1 % |
| SRR21224508 | %PF 84,8 % |
| SRR21224518 | GC 57,9 % y adaptador 5,6 % |
| SRR21224519, SRR21224520 | rem 6w quedaba con 2 réplicas |
| SRR10067436, SRR10067437 | tamaño de inserto anómalo |

## Parámetros
- fastp: Phred ≥ 20; longitud mínima según plataforma (20, 30 o 36 pb).
- STAR: GRCm39 Ensembl 111; índices sjdbOverhang 79 (< 80 pb) y 159 (80–160 pb).
- Filtro de genes: rowSums(counts) ≥ 10.
- DESeq2: DEG = padj < 0,05 y |log2FC| > 1. Diseño: el más completo con rango completo entre ~dataset + tejido + grupo, ~tejido + grupo, ~dataset + grupo, ~grupo.
- GSEA: ranking por estadístico de Wald; GO BP y KEGG; minGSSize 5, maxGSSize 1000; semilla 123.

## Post-crónico
- Control común: Control 24 sem (misma edad que 24 sem).
- 24 sem tiene 9 hembras y 0 machos → análisis de sensibilidad solo con hembras (metadata/metadata_sexo.csv).
