#!/bin/bash

# Ruta al archivo de anotación de ratón que usaste para el índice
GTF_FILE="/data/Database/referencia/Mus_musculus.GRCm39.111.gtf"

# Carpeta donde se guardarán las matrices de conteo finales
OUTPUT_DIR="/data/Database/Conteos"
mkdir -p "$OUTPUT_DIR"

echo "=========================================================="
echo " INICIANDO CUANTIFICACIÓN DE EXPRESIÓN GÉNICAS (featureCounts)"
echo "=========================================================="

# ---------------------------------------------------------
# PARTE A: Muestras Single-End (db5, db8, db5_cchp, db11)
# ---------------------------------------------------------
echo "--> Procesando carpetas SINGLE-END..."

# Buscamos todos los archivos .bam de las carpetas SE
BAM_SE=$(find /data/Database/Alineamiento/db5/ /data/Database/Alineamiento/db8/ /data/Database/Alineamiento/db5_cchp/ /data/Database/Alineamiento/db11/ -name "*_Aligned.sortedByCoord.out.bam" 2>/dev/null)

if [ -z "$BAM_SE" ]; then
    echo "No se encontraron archivos BAM para Single-End."
else
    # Ejecutar conteo SE usando 12 hilos (-T 12)
    featureCounts -T 12 \
                   -a "$GTF_FILE" \
                   -t exon \
                   -g gene_id \
                   -o "$OUTPUT_DIR/matriz_counts_single_end.txt" \
                   $BAM_SE
fi

# ---------------------------------------------------------
# PARTE B: Muestras Paired-End (db9)
# ---------------------------------------------------------
echo "--> Procesando carpeta PAIRED-END (db9)..."

BAM_PE=$(find /data/Database/Alineamiento/db9/ -name "*_Aligned.sortedByCoord.out.bam" 2>/dev/null)

if [ -z "$BAM_PE" ]; then
    echo "No se encontraron archivos BAM para Paired-End (db9)."
else
    # Ejecutar conteo PE (-p indica que son lecturas pareadas)
    featureCounts -p \
                   -T 12 \
                   -a "$GTF_FILE" \
                   -t exon \
                   -g gene_id \
                   -o "$OUTPUT_DIR/matriz_counts_paired_end_db9.txt" \
                   $BAM_PE
fi

echo "=========================================================="
echo " ¡PROCESO TERMINADO CON ÉXITO!"
echo " Revisar matrices de conteo en: $OUTPUT_DIR"
echo "=========================================================="
