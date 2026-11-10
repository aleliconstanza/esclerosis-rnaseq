#!/usr/bin/env bash
# ============================================================
# correr_postcronico.sh
# Corre el análisis post-crónico (08) y las figuras del póster (09).
#
# Uso (desde cualquier lugar):
#   bash scripts/correr_postcronico.sh            # análisis + figuras
#   bash scripts/correr_postcronico.sh figuras    # solo figuras
# ============================================================
set -euo pipefail
RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
LOGS="$RAIZ/intermedios/08_postcronico/logs"; mkdir -p "$LOGS"
FECHA="$(date +%Y%m%d_%H%M)"
export LC_ALL="${LC_ALL:-C.UTF-8}"

if [[ "${1:-todo}" != "figuras" ]]; then
  DDS="$RAIZ/intermedios/06_normalizacion/dds_normalizado.rds"
  [[ -e "$DDS" ]] || { echo "ERROR: falta $DDS (enlace al disco)"; exit 1; }
  cd "$RAIZ/scripts/08_postcronico"
  echo "[1/3] Pipeline oficial (main.R)";        Rscript main.R           2>&1 | tee "$LOGS/main_$FECHA.log"
  echo "[2/3] Tejido y solo hembras";             Rscript main_extendido.R 2>&1 | tee "$LOGS/extendido_$FECHA.log"
fi
echo "[3/3] Figuras del póster"
cd "$RAIZ/scripts/09_figuras_poster"
Rscript figuras_poster.R 2>&1 | tee "$LOGS/figuras_$FECHA.log"
echo
echo "Tablas:  resultados/08_postcronico/tablas/"
echo "Figuras: resultados/09_figuras_poster/"
