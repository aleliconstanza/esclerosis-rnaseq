#!/usr/bin/env bash
# ============================================================
# organizar_repo.sh
# Copia lo liviano de /data/Database (scripts, reportes, tablas,
# metadata) al repositorio en tu cuenta, y deja un acceso directo
# al disco para lo pesado (FASTQ, BAM, índices STAR, .rds).
#
# Uso:
#   bash organizar_repo.sh            # SIMULA: muestra qué copiaría, no toca nada
#   bash organizar_repo.sh ejecutar   # copia de verdad
#
# Nunca modifica /data/Database: solo lee desde ahí.
# Nunca sobrescribe archivos que ya existen en el repo.
# ============================================================
set -uo pipefail

ORIGEN="${ORIGEN:-/data/Database}"
REPO="${REPO:-$HOME/esclerosis-rnaseq}"
MAX="50M"                       # tope por archivo (GitHub rechaza >100 MB)
MODO="${1:-simular}"

[[ -d "$ORIGEN" ]] || { echo "No existe $ORIGEN"; exit 1; }
[[ -d "$REPO/.git" ]] || { echo "No existe el repo en $REPO (clónalo primero)"; exit 1; }

if [[ "$MODO" == "ejecutar" ]]; then SIM=""; else SIM="--dry-run"; echo ">>> MODO SIMULACIÓN (nada se copia)"; fi

# Qué se copia: solo tipos livianos. Qué nunca se copia: carpetas y archivos pesados.
BASE=(
  --exclude='Datos_Crudos/' --exclude='Trimmed/' --exclude='referencia/'
  --exclude='star_index_80/' --exclude='star_index_160/' --exclude='.git/'
  --exclude='*.fastq' --exclude='*.fastq.gz' --exclude='*.fq' --exclude='*.fq.gz'
  --exclude='*.bam' --exclude='*.bai' --exclude='*.sam' --exclude='*.rds' --exclude='*.RData'
  --include='*/'
  --include='*.sh' --include='*.R' --include='*.py' --include='*.md' --include='*.txt'
  --include='*.csv' --include='*.tsv' --include='*.xlsx' --include='*.json'
  --include='*.html' --include='*.pdf' --include='*.png' --include='*.jpg' --include='*.svg'
  --include='*.out' --include='*.summary' --include='*.log'
  --exclude='*'
)
FILTROS=(--exclude='Normalización/' "${BASE[@]}")
OPC=(-a --prune-empty-dirs --ignore-existing --max-size="$MAX" --itemize-changes $SIM)

echo
echo "[1/4] Copiando archivos livianos de $ORIGEN → $REPO"
rsync "${OPC[@]}" "${FILTROS[@]}" "$ORIGEN/" "$REPO/" | grep '^>f' | sed 's/^>f[^ ]* /  + /'

echo
echo "[2/4] Normalización/ (con tilde en el disco) → Normalizacion/ (sin tilde en el repo)"
if [[ -d "$ORIGEN/Normalización" ]]; then
  rsync "${OPC[@]}" "${BASE[@]}" "$ORIGEN/Normalización/" "$REPO/Normalizacion/" \
    | grep '^>f' | sed 's/^>f[^ ]* /  + /'
fi

echo
echo "[3/4] Archivos que existen en ambos lados pero son DISTINTOS (no se tocaron; revísalos):"
rsync -rcn --existing --itemize-changes "${FILTROS[@]}" "$ORIGEN/" "$REPO/" | grep '^>f' | sed 's/^>f[^ ]* /  ≠ /'
[[ -d "$ORIGEN/Normalización" ]] && rsync -rcn --existing --itemize-changes "${BASE[@]}" "$ORIGEN/Normalización/" "$REPO/Normalizacion/" | grep '^>f' | sed 's/^>f[^ ]* /  ≠ Normalizacion\//'

echo
echo "[4/4] Accesos directos a lo pesado"
if [[ "$MODO" == "ejecutar" ]]; then
  ln -sfn "$ORIGEN" "$REPO/disco"
  echo "  disco → $ORIGEN   (ej.: ls disco/Datos_Crudos)"
  DDS="$(find "$ORIGEN" -name dds_normalizado.rds 2>/dev/null | head -1)"
  if [[ -n "$DDS" ]]; then
    mkdir -p "$REPO/post-cronico/datos"
    ln -sfn "$DDS" "$REPO/post-cronico/datos/dds_normalizado.rds"
    echo "  post-cronico/datos/dds_normalizado.rds → $DDS"
  else
    echo "  (no encontré dds_normalizado.rds en el disco)"
  fi
  for linea in "/disco" "*.rds" "post-cronico/datos/dds_normalizado.rds"; do
    grep -qxF "$linea" "$REPO/.gitignore" 2>/dev/null || echo "$linea" >> "$REPO/.gitignore"
  done
  echo "  .gitignore actualizado (el acceso directo y los .rds no se suben a GitHub)"
else
  echo "  (se crearán al ejecutar: disco → $ORIGEN y el enlace a dds_normalizado.rds)"
fi

echo
echo "Tamaño actual del repo (sin .git ni disco):"
du -sh --exclude=.git --exclude=disco "$REPO" 2>/dev/null
echo
[[ "$MODO" == "ejecutar" ]] && echo "Listo. Revisa con:  cd $REPO && git status" \
                            || echo "Si la lista se ve bien, corre:  bash organizar_repo.sh ejecutar"
