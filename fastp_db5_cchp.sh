#!/bin/bash
# fastp - Grupo 3: DB5 Especial (Base de datos independiente db5_cchp)

# Ruta absoluta exacta al binario de tu entorno bio
FASTP_EXEC="/home/nvaldes/miniconda3/envs/bio/bin/fastp"

DIR_CRUDOS="/data/Database/Datoscrudos/db5_cchp"
DIR_RESULTADOS="/data/Database/Trimming/db5_cchp"

echo "🚀 Iniciando fastp para la Base de Datos Independiente: db5_cchp..."

# 1. Buscar todas las subcarpetas dentro de db5_cchp (ej: control24wcc)
find "$DIR_CRUDOS" -type d | while read -r subcarpeta; do
    
    # Verificar si la subcarpeta contiene archivos .fastq antes de procesar
    if ls "$subcarpeta"/*.fastq >/dev/null 2>&1; then
        
        # Extraer el nombre del timepoint (ej: control24wcc)
        timepoint=$(basename "$subcarpeta")
        OUT_FINAL="$DIR_RESULTADOS/$timepoint"
        
        # Crear la estructura limpia dentro de Trimming/db5_cchp/
        mkdir -p "$OUT_FINAL/trimmed_reads" "$OUT_FINAL/fastqc_trimming"
        
        # 2. Procesar la totalidad de los archivos .fastq dentro de esa carpeta
        for archivo in "$subcarpeta"/*.fastq; do
            [ -f "$archivo" ] || continue
            nombre_base=$(basename "$archivo" .fastq)
            
            echo "🧬 Procesando en db5_cchp/$timepoint: $nombre_base"
            
            $FASTP_EXEC \
              --in1 "$archivo" \
              --out1 "$OUT_FINAL/trimmed_reads/${nombre_base}_trimmed.fastq" \
              --qualified_quality_phred 20 \
              --length_required 30 \
              --trim_poly_g --poly_g_min_len 8 \
              --thread 8 \
              --json "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.json" \
              --html "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.html"
        done
    fi
done
echo "✅ Análisis de la base de datos db5_cchp finalizado por completo."
