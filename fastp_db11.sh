#!/bin/bash
# fastp - Grupo 4: DB11 Completo (La totalidad de archivos)

# Ruta absoluta exacta al binario de tu entorno bio
FASTP_EXEC="/home/nvaldes/miniconda3/envs/bio/bin/fastp"

DIR_CRUDOS="/data/Database/Datoscrudos/db11"
DIR_RESULTADOS="/data/Database/Trimming/db11"

echo "🚀 Iniciando fastp para la totalidad de Database 11 (Single-End 151pb)..."

find "$DIR_CRUDOS" -type d | while read -r subcarpeta; do
    if ls "$subcarpeta"/*.fastq >/dev/null 2>&1; then
        
        # Extraer DIRECTAMENTE el nombre del timepoint (ej: control0w, 3w, etc)
        timepoint=$(basename "$subcarpeta")
        OUT_FINAL="$DIR_RESULTADOS/$timepoint"
        
        mkdir -p "$OUT_FINAL/trimmed_reads" "$OUT_FINAL/fastqc_trimming"
        
        for archivo in "$subcarpeta"/*.fastq; do
            [ -f "$archivo" ] || continue
            nombre_base=$(basename "$archivo" .fastq)
            
            echo "🧬 Procesando DB11 en $timepoint: $nombre_base"
            
            # ¡CORREGIDO AQUÍ!: Ahora usamos la variable con el signo $ delante
            $FASTP_EXEC \
              --in1 "$archivo" \
              --out1 "$OUT_FINAL/trimmed_reads/${nombre_base}_trimmed.fastq" \
              --qualified_quality_phred 20 \
              --length_required 36 \
              --trim_poly_g --poly_g_min_len 8 \
              --thread 8 \
              --json "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.json" \
              --html "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.html"
        done
    fi
done
echo "✅ Grupo 4 finalizado por completo."
