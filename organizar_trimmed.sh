#!/bin/bash
# Script para clonar Trimming a Trimmed simplificando rutas y aplicando Blacklist de SRRs

DIR_ORIGEN="/data/Database/Trimming"
DIR_DESTINO="/data/Database/Trimmed"

echo "🚀 Iniciando clonación limpia en: $DIR_DESTINO..."

# 1. Buscar todas las carpetas "trimmed_reads" dentro de la estructura actual
find "$DIR_ORIGEN/db8" -type d -name "trimmed_reads" | while read -r carpeta_trimmed; do
    
    # Extraer la ruta intermedia (ej: db5/control0w o db11/des5w)
    ruta_relativa=$(echo "$carpeta_trimmed" | sed "s|^$DIR_ORIGEN/||" | sed 's|/trimmed_reads||')
    NUEVO_DIR="$DIR_DESTINO/$ruta_relativa"
    
    if ls "$carpeta_trimmed"/*.fastq >/dev/null 2>&1; then
        
        # Crear la carpeta destino del timepoint solo si no existe
        mkdir -p "$NUEVO_DIR"
        
        # 2. Procesar los archivos .fastq dentro de la subcarpeta
        for archivo_fastq in "$carpeta_trimmed"/*.fastq; do
            [ -f "$archivo_fastq" ] || continue
            
            nombre_archivo=$(basename "$archivo_fastq")
            
            # Convertir el nombre a minúsculas para evitar problemas de "SRR" vs "srr"
            nombre_minusculas=$(echo "$nombre_archivo" | tr '[:upper:]' '[:lower:]')
            
            # 3. FILTRO / BLACKLIST: Verificar si el archivo contiene algún SRR descartado
            case "$nombre_minusculas" in
                # --- DATA 5 (Calidad y GC) ---
                *srr21224515*)
                    echo "❌ Excluido (db5 - Calidad): $nombre_archivo"
                    continue
                    ;;
                *srr21224508*)
                    echo "❌ Excluido (db5 - Calidad): $nombre_archivo"
                    continue
                    ;;
                *srr21224518*)
                    echo "❌ Excluido (db5 - Contenido GC): $nombre_archivo"
                    continue
                    ;;
                
                # --- DATA 9 (Insert Size) ---
                *srr10067436*)
                    echo "❌ Excluido (db9 - Tamaño Inserto): $nombre_archivo"
                    continue
                    ;;
                *srr10067437*)
                    echo "❌ Excluido (db9 - Tamaño Inserto): $nombre_archivo"
                    continue
                    ;;
                    
                # Si no está en la lista negra, se procede a copiar
                *)
                    echo "  ➡️ Copiando a $ruta_relativa: $nombre_archivo"
                    cp "$archivo_fastq" "$NUEVO_DIR/$nombre_archivo"
                    ;;
            esac
        done
    fi
done

# 4. Limpieza de directorios vacíos (en caso de que un timepoint se haya quedado sin archivos)
find "$DIR_DESTINO" -type d -empty -delete

echo "✅ ¡Filtro y reestructuración completados con éxito!"
echo "Los datos limpios para el alineamiento están en: $DIR_DESTINO"
