#!/bin/bash
# MultiQC - Reportes individuales por Timepoint y Globales por Set

MULTIQC_EXEC="multiqc"

DIR_RESULTADOS="/data/Database/Trimming"
DIR_MULTIQC="/data/Database/MultiQC"
DIR_GLOBAL="$DIR_MULTIQC/Database"
DIR_GLOBAL_GLOBAL="$DIR_MULTIQC/Global"

echo "🚀 Iniciando generación de reportes MultiQC..."

# ==============================================================================
# PARTE 1: Reportes por Base de Datos y Timepoint (Local)
# ==============================================================================
echo "📂 Creando estructura de carpetas y reportes locales..."

# Busca directorios a nivel "timepoint" (Ej: Trimming/db5cchp/T1)
find "$DIR_RESULTADOS" -mindepth 2 -maxdepth 2 -type d | while read -r carpeta_trimming; do
    
    timepoint=$(basename "$carpeta_trimming")              
    db_nombre=$(basename "$(dirname "$carpeta_trimming")") 
    
    if ls "$carpeta_trimming/fastqc_trimming"/*_fastp.json >/dev/null 2>&1; then
        
        OUT_LOCAL="$DIR_MULTIQC/$db_nombre/$timepoint"
        mkdir -p "$OUT_LOCAL"
        
        echo "📊 Generando MultiQC local para: $db_nombre -> $timepoint"
        
        $MULTIQC_EXEC "$carpeta_trimming/fastqc_trimming" \
            -o "$OUT_LOCAL" \
            -n "MultiQC_${db_nombre}_${timepoint}" \
            -q
    fi
done

echo "✅ Reportes locales finalizados."
echo "---------------------------------------------------"

# ==============================================================================
# PARTE 2: Reportes Globales por Set
# ==============================================================================
echo "🌍 Generando reportes globales..."

# Definimos una lista con los nombres exactos de tus 5 carpetas base
BASES_DE_DATOS=("db5" "db5_cchp" "db8" "db9" "db11")

# Iteramos sobre cada base de datos para generar su reporte global
for db in "${BASES_DE_DATOS[@]}"; do
    if [ -d "$DIR_RESULTADOS/$db" ]; then
        echo "📈 Procesando Global: $db"
        
        # Creamos la carpeta global para esta base de datos específica
        mkdir -p "$DIR_GLOBAL/$db"
        
        # MultiQC escanea toda la carpeta base (incluyendo todos sus timepoints)
        $MULTIQC_EXEC "$DIR_RESULTADOS/$db" \
            -o "$DIR_GLOBAL/$db" \
            -n "MultiQC_Global_$db" \
            -q
    else
        echo "⚠️  Aviso: No se encontró el directorio $DIR_RESULTADOS/$db. Saltando..."
    fi
done

echo "✅ Reportes globales finalizados."
echo "🎉 Todo el flujo de MultiQC ha terminado exitosamente."

# PARTE 3: Reporte Máster (Comparativo de todo tu proyecto usando solo fastp)
# ==============================================================================
echo "🔄 Generando Reporte Máster: Consolidado de todo el proyecto..."

mkdir -p "$DIR_GLOBAL_GLOBAL"

echo "📊 Analizando todos los JSON de fastp en Trimming..."

# Al pasarle la raíz de Trimming, MultiQC juntará absolutamente todas las muestras 
# de todas tus bases de datos en un solo HTML, mostrando sus curvas antes/después automáticamente.
$MULTIQC_EXEC "$DIR_RESULTADOS" \
    -o "$DIR_GLOBAL_GLOBAL" \
    -n "MultiQC_Master_Proyecto_Completo" \
    -q

echo "✅ Reporte Máster finalizado."
echo "🎉 Todo el flujo de MultiQC ha terminado exitosamente!"
echo "📂 Encuentra tu reporte de todo el proyecto en: $DIR_GLOBAL_GLOBAL/MultiQC_Master_Proyecto_Completo.html"

