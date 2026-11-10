#!/bin/bash
# Uso: bash reporte_pesos.sh
# Genera: reporte_pesos.csv

DIR_CRUDOS="/data/Database/Datoscrudos"
REPORTE="/data/Database/Resultados/reporte_pesos.csv"

echo "Base_Datos,Set,Nombre_Archivo,Tamaño_GB,Longitud_bp" > "$REPORTE"

find "$DIR_CRUDOS" \
    \( -path "$DIR_CRUDOS/db5/*" \
       -o -path "$DIR_CRUDOS/db5_cchp/*" \
       -o -path "$DIR_CRUDOS/db8/*" \
       -o -path "$DIR_CRUDOS/db9/*" \
       -o -path "$DIR_CRUDOS/db11/*" \) \
    -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read -r archivo; do

    ruta_relativa=$(echo "$archivo" | sed "s|^$DIR_CRUDOS/||")
    db_nombre=$(echo "$ruta_relativa" | cut -d'/' -f1)
    set_nombre=$(echo "$ruta_relativa" | cut -d'/' -f2-)

    peso_gb=$(stat -c%s "$archivo" | awk '{printf "%.4f", $1/1073741824}')

    if [[ "$archivo" == *.gz ]]; then
        longitud=$(zcat "$archivo" | sed -n '2p' | tr -d '\n' | wc -c)
    else
        longitud=$(sed -n '2p' "$archivo" | tr -d '\n' | wc -c)
    fi

    echo "$db_nombre,$set_nombre,$(basename "$archivo"),$peso_gb,$longitud"

done >> "$REPORTE"

echo "Listo: $REPORTE"
