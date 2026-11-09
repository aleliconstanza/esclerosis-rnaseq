#!/bin/bash

# Agregamos la ruta de tu FastQC al PATH del sistema para que Trim Galore lo encuentre solo
export PATH="/home/nvaldes/FastQC:$PATH"

echo "=== Iniciando Flujo de Trimming y Control de Calidad Post-Trim (Corregido) ==="

# Buscar carpetas contenedoras que tengan archivos .fastq (excluyendo la carpeta 'extra')
find . -type d -name "extra" -prune -o -name "*.fastq" -print | while read -r archivo; do
    subcarpeta=$(dirname "$archivo")
    
    if [ -f "$subcarpeta/.procesado_trim" ]; then
        continue
    fi
    
    echo "------------------------------------------------------------"
    echo "🚀 PROCESANDO CARPETA: $subcarpeta"
    echo "------------------------------------------------------------"

    mkdir -p "$subcarpeta/trimmed_reads"
    mkdir -p "$subcarpeta/fastqc_trimming"
    mkdir -p "$subcarpeta/multiqc_after"

    echo "✂️  Iniciando Trim Galore..."
    
    # Se eliminó --path_to_fastqc ya que ahora se detecta automáticamente
    trim_galore --fastqc --quality 20 --length 20 --cores 4 \
                --output_dir "$subcarpeta/trimmed_reads" \
                "$subcarpeta"/*.fastq

    echo "📂 Organizando reportes FastQC post-trimming..."
    find "$subcarpeta/trimmed_reads" -name "*_fastqc.html" -exec mv {} "$subcarpeta/fastqc_trimming/" \;
    find "$subcarpeta/trimmed_reads" -name "*_fastqc.zip" -exec mv {} "$subcarpeta/fastqc_trimming/" \;

    echo "📉 Generando MultiQC final..."
    multiqc "$subcarpeta/fastqc_trimming" -o "$subcarpeta/multiqc_after"

    touch "$subcarpeta/.procesado_trim"
    echo "✅ Carpeta finalizada con éxito."
    echo ""

done

find . -name ".procesado_trim" -delete
echo "=== ¡Todo el dataset ha sido recortado y analizado con éxito! ==="

