#!/bin/bash
# fastp - Grupo 1: Database 9 (Paired-End)

# Ruta absoluta exacta al binario de tu entorno bio
FASTP_EXEC="/home/nvaldes/miniconda3/envs/bio/bin/fastp"

DIR_CRUDOS="/data/Database/Datoscrudos/db9"
DIR_RESULTADOS="/data/Database/Trimming/db9"

echo "🚀 Iniciando fastp para la totalidad de Database 9 (Paired-End)..."

find "$DIR_CRUDOS" -type d | while read -r subcarpeta; do
    for r1_file in "$subcarpeta"/*_1.fastq "$subcarpeta"/*_R1*.fastq; do
        [ -e "$r1_file" ] || continue
        
        if [[ "$r1_file" == *"_1.fastq" ]]; then
            r2_file="${r1_file/_1.fastq/_2.fastq}"
        else
            r2_file="${r1_file/_R1/_R2}"
        fi
        
        timepoint=$(basename "$subcarpeta")
        OUT_FINAL="$DIR_RESULTADOS/$timepoint"
        
        mkdir -p "$OUT_FINAL/trimmed_reads" "$OUT_FINAL/fastqc_trimming"
        nombre_base=$(basename "$r1_file" | sed -E 's/(_1|_R1).*//')
        
        echo "🧬 Procesando par en $timepoint: $nombre_base"
        
        $FASTP_EXEC \
          --in1 "$r1_file" \
          --in2 "$r2_file" \
          --out1 "$OUT_FINAL/trimmed_reads/${nombre_base}_trim_1.fastq" \
          --out2 "$OUT_FINAL/trimmed_reads/${nombre_base}_trim_2.fastq" \
          --qualified_quality_phred 20 \
          --length_required 20 \
          --trim_poly_g \
          --thread 8 \
          --json "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.json" \
          --html "$OUT_FINAL/fastqc_trimming/${nombre_base}_fastp.html"
    done
done
echo "✅ Grupo 1 finalizado."
