#!/bin/bash
# fastp - Grupo 2: DB8 + DB5 Completo (Procesa todo lo que encuentre en ambas carpetas)

# Ruta absoluta exacta al binario de tu entorno bio
FASTP_EXEC="/home/nvaldes/miniconda3/envs/bio/bin/fastp"

DIR_CRUDOS="/data/Database/Datoscrudos"
DIR_RESULTADOS="/data/Database/Trimming"

echo "🚀 Iniciando fastp para la totalidad de DB8 y DB5 Regular (Single-End)..."

# 1. Buscar todas las subcarpetas dentro de db5 y db8 de forma directa
find "$DIR_CRUDOS/db5" "$DIR_CRUDOS/db8" -type d 2>/dev/null | while read -r subcarpeta; do
    
    # Verificar si la subcarpeta contiene archivos .fastq antes de procesar
    if ls "$subcarpeta"/*.fastq >/dev/null 2>&1; then
        
        # Detectar dinámicamente si la carpeta madre es db5 o db8
        db_nombre=$(echo "$subcarpeta" | sed "s|^$DIR_CRUDOS/||" | cut -d'/' -f1)
        
        # Extraer el nombre del timepoint (ej: control0w, 3w, 5w)
        timepoint=$(basename "$subcarpeta")
        
        # Estructurar la salida separada por base de datos y timepoint
        OUT_FINAL="$DIR_RESULTADOS/$db_nombre/$timepoint"
        
        mkdir -p "$OUT_FINAL/trimmed_reads" "$OUT_FINAL/fastqc_trimming"
        
        # 2. Procesar la totalidad de archivos .fastq en este timepoint
        for archivo in "$subcarpeta"/*.fastq; do
            [ -f "$archivo" ] || continue
            nombre_base=$(basename "$archivo" .fastq)
            
            echo "🧬 Procesando SE en $db_nombre/$timepoint: $nombre_base"
            
            $FASTP_EXEC \
              --in1 "$archivo" \
              --out1 "$OUT_FINAL/trimmed_reads/${nombre_base}_trimmed.fastq" \
              --qualified_quality_phred 20 \
              --length_required 20 \
              --trim_poly_g \
              --thread 8 \
              --json "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.json" \
              --html "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.html"
        done
    fi
done
echo "✅ Grupo 2 finalizado. Datos de db5 y db8 procesados por completo."
