#!/bin/bash

# Lista de tus archivos de metadatos
ARCHIVOS_METADATA=("metadata5.csv" "metadata8.csv" "metadata9.csv" "metadata11.csv")

echo "=== Iniciando automatización de descargas SRA ==="

# Bucle para recorrer cada archivo metadata
for archivo in "${ARCHIVOS_METADATA[@]}"; do
    if [ -f "$archivo" ]; then
        echo "--------------------------------------------------"
        echo " Leyendo códigos del archivo: $archivo"
        echo "--------------------------------------------------"
        
        # Extrae la primera columna, filtra solo líneas que tengan SRR y recorre cada ID
        cut -d',' -f1 "$archivo" | grep -E '^SRR[0-9]+' | while read -r srr_id; do
            echo "📥 Descargando: $srr_id ..."
            
            # Ejecuta tu comando patrón
            fasterq-dump --progress "$srr_id"
            
            echo "✅ Finalizado: $srr_id"
            echo ""
        done
    else
        echo "⚠️ Advertencia: El archivo $archivo no se encuentra en esta carpeta."
    fi
done

echo "=== ¡Todas las descargas han finalizado con éxito! ==="
