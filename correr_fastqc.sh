#!/bin/bash

# Este script busca archivos .fastq (descomprimidos) ignorando la carpeta 'extra'
# y ejecuta FastQC guardando los resultados ordenadamente.

find . -type d -name "extra" -prune -o -name "*.fastq" -print | while read -r archivo; do
    # 1. Obtener la ruta de la carpeta donde está el archivo actual
    carpeta_origen=$(dirname "$archivo")
    
    # 2. Crear la carpeta 'fastqc_results' dentro de esa ubicación si no existe
    mkdir -p "$carpeta_origen/fastqc_results"
    
    echo "Procesando: $archivo"
    
    # 3. Llamado a tu ejecutable de FastQC
    /home/nvaldes/FastQC/fastqc "$archivo" -o "$carpeta_origen/fastqc_results" -t 4
done

