#!/usr/bin/env bash
# ============================================================
# correr_postcronico.sh
# Corre el análisis post-crónico completo y las figuras del póster.
#
# Uso (desde la raíz del repositorio):
#   bash correr_postcronico.sh            # todo
#   bash correr_postcronico.sh figuras    # solo figuras (usa resultados existentes)
#
# Requisito: post-cronico/datos/dds_normalizado.rds
#   (lo genera Normalizacion/normalizar.R; no está en Git por tamaño)
# ============================================================
set -euo pipefail

RAIZ="$(cd "$(dirname "$0")" && pwd)"
cd "$RAIZ/post-cronico"
mkdir -p logs
FECHA="$(date +%Y%m%d_%H%M)"
export LC_ALL="${LC_ALL:-C.UTF-8}"

if [[ "${1:-todo}" != "figuras" ]]; then
  if [[ ! -f datos/dds_normalizado.rds ]]; then
    echo "ERROR: falta post-cronico/datos/dds_normalizado.rds"
    echo "Cópialo con:  cp <ruta>/dds_normalizado.rds post-cronico/datos/"
    exit 1
  fi
  [[ -f datos/metadata_sexo.csv ]] || echo "AVISO: falta datos/metadata_sexo.csv (no correrá la comparación solo hembras)"

  echo "[1/3] Pipeline oficial (main.R)... log: post-cronico/logs/main_$FECHA.log"
  Rscript main.R 2>&1 | tee "logs/main_$FECHA.log"

  echo "[2/3] Pipeline extendido: tejido y solo hembras (main_extendido.R)..."
  Rscript main_extendido.R 2>&1 | tee "logs/extendido_$FECHA.log"
fi

echo "[3/3] Figuras del póster..."
Rscript figuras_poster/figuras_poster.R resultados/tablas ../Normalizacion/metadata_muestras.csv \
  datos/metadata_sexo.csv figuras_poster/salida 2>&1 | tee "logs/figuras_$FECHA.log"

echo
echo "Listo."
echo "  Tablas:   post-cronico/resultados/tablas/"
echo "  Diseños:  post-cronico/resultados/tablas/disenos_usados.csv"
echo "  Figuras:  post-cronico/figuras_poster/salida/"
