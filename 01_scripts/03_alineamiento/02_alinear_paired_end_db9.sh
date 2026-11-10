#!/bin/bash

# Configuración para db9
INDEX_DIR="star_index_80"
DB="db9"
BASE_PATH="/data/Database/Trimmed/${DB}/*"

echo "=========================================================="
echo " Procesando PAIRED-END para: $DB (Índice: $INDEX_DIR)"
echo "=========================================================="

# Buscar archivos R1 que terminen en _1.fastq en cualquier timepoint
for R1 in $BASE_PATH/*_1.fastq; do
    
    [ -e "$R1" ] || continue
    
    # Identificar el archivo R2 correspondiente reemplazando _1.fastq por _2.fastq
    R2="${R1/_1.fastq/_2.fastq}"
    
    # Extraer nombre base de la muestra omitiendo el _1.fastq
    SAMPLE=$(basename "$R1" _1.fastq)
    
    # Extraer el nombre del timepoint
    TIMEPOINT=$(basename $(dirname "$R1"))
    
    # Carpeta de salida ordenada
    OUTPUT_DIR="/data/Database/Alineamiento/${DB}/${TIMEPOINT}"
    mkdir -p "$OUTPUT_DIR"
    
    echo "--> Alineando PE: Muestra $SAMPLE ($TIMEPOINT)"
    
    # Ejecutar star pasándole la pareja de archivos R1 y R2
    STAR --runThreadN 12 \
         --genomeDir "$INDEX_DIR" \
         --readFilesIn "$R1" "$R2" \
         --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate
         
done

echo "¡Alineamiento Paired-End de db9 completado con éxito!"
