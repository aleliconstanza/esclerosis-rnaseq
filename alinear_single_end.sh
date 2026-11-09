#!/bin/bash

# Lista de todas las carpetas que son Single-End
CARPETAS=("db5" "db8" "db5_cchp" "db11")

for DB in "${CARPETAS[@]}"; do
    
    # ASIGNACIÓN AUTOMÁTICA DEL ÍNDICE SEGÚN LA CARPETA
    if [ "$DB" == "db5" ] || [ "$DB" == "db8" ]; then
        INDEX_DIR="star_index_80"
    else
        INDEX_DIR="star_index_160"
    fi
    
    # Nueva ruta: Entra a la base de datos y a cualquier carpeta de timepoint (*)
    BASE_PATH="/data/Database/Trimmed/${DB}/*"
    
    echo "=========================================================="
    echo " Procesando SINGLE-END para: $DB (Índice: $INDEX_DIR)"
    echo "=========================================================="
    
    # Buscamos archivos .fastq descomprimidos directamente en el timepoint
    for R1 in $BASE_PATH/*.fastq; do
        
        [ -e "$R1" ] || continue
        
        # Extraer nombre de la muestra (quitando .fastq)
        SAMPLE=$(basename "$R1" .fastq)
        
        # Extraer el nombre del timepoint (es el nombre de la carpeta contenedora del archivo)
        TIMEPOINT=$(basename $(dirname "$R1"))
        
        # Carpeta de salida ordenada
        OUTPUT_DIR="/data/Database/Alineamiento/${DB}/${TIMEPOINT}"
        mkdir -p "$OUTPUT_DIR"
        
        echo "--> Alineando SE: Muestra $SAMPLE ($TIMEPOINT)"
        
        # Ejecutar star
        STAR --runThreadN 12 \
             --genomeDir "$INDEX_DIR" \
             --readFilesIn "$R1" \
             --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE}_" \
             --outSAMtype BAM SortedByCoordinate
             
    done
done

echo "¡Alineamiento Single-End completado con éxito!"
