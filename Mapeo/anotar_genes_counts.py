#!/usr/bin/env python3
"""
Anota todos los genes presentes en la matriz de conteos usando el GTF
de referencia, e identifica cuáles coinciden con la lista de literatura
(genes_literatura_mapeado.csv generado por mapear_genes_literatura.R).

Uso:
    python3 anotar_genes_counts.py <gtf> <counts_csv> <literatura_csv> <salida_dir>

Ejemplo:
    python3 anotar_genes_counts.py \
        /data/Database/referencia/Mus_musculus.GRCm39.111.gtf \
        /data/Database/Resultados/matriz_conteos_filtrada.csv \
        genes_literatura_mapeado.csv \
        /data/Database/Resultados/anotacion
"""

import sys
import os
import re
import csv


def parsear_gtf(ruta_gtf):
    """Extrae gene_id -> gene_name desde las líneas 'gene' del GTF."""
    mapa_id_simbolo = {}
    mapa_id_biotipo = {}
    with open(ruta_gtf, "r") as f:
        for linea in f:
            if linea.startswith("#"):
                continue
            campos = linea.rstrip("\n").split("\t")
            if len(campos) < 9 or campos[2] != "gene":
                continue
            atributos = campos[8]

            m_id = re.search(r'gene_id "([^"]+)"', atributos)
            m_name = re.search(r'gene_name "([^"]+)"', atributos)
            m_biotipo = re.search(r'gene_biotype "([^"]+)"', atributos)

            if not m_id:
                continue
            gene_id = m_id.group(1)
            gene_name = m_name.group(1) if m_name else None
            biotipo = m_biotipo.group(1) if m_biotipo else None

            mapa_id_simbolo[gene_id] = gene_name
            mapa_id_biotipo[gene_id] = biotipo

    return mapa_id_simbolo, mapa_id_biotipo


def leer_ids_counts(ruta_counts):
    """Lee solo la primera columna (IDs de gen) de la matriz de conteos."""
    ids = []
    with open(ruta_counts, "r") as f:
        lector = csv.reader(f)
        encabezado = next(lector)
        for fila in lector:
            if not fila:
                continue
            gene_id = fila[0].strip('"')
            ids.append(gene_id)
    return ids


def leer_literatura(ruta_lit):
    """Lee genes_literatura_mapeado.csv -> dict ensembl_id -> simbolo_humano."""
    mapa = {}
    if not os.path.exists(ruta_lit):
        print(f"  AVISO: no se encontró {ruta_lit}, se omite el cruce con literatura.")
        return mapa
    with open(ruta_lit, "r") as f:
        lector = csv.DictReader(f)
        for fila in lector:
            ens = fila.get("ensembl_mouse", "").strip()
            if ens and ens.upper() != "NA":
                mapa[ens] = fila.get("simbolo_humano", "")
    return mapa


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)

    ruta_gtf, ruta_counts, ruta_lit, dir_salida = sys.argv[1:5]
    os.makedirs(dir_salida, exist_ok=True)

    print(f"[1/4] Parseando GTF: {ruta_gtf}")
    id_a_simbolo, id_a_biotipo = parsear_gtf(ruta_gtf)
    print(f"      Genes en GTF: {len(id_a_simbolo)}")

    print(f"[2/4] Leyendo IDs de la matriz de conteos: {ruta_counts}")
    ids_counts = leer_ids_counts(ruta_counts)
    print(f"      Genes en counts: {len(ids_counts)}")

    print(f"[3/4] Cruzando con lista de literatura: {ruta_lit}")
    id_a_lit = leer_literatura(ruta_lit)
    print(f"      Genes de literatura con ID válido: {len(id_a_lit)}")

    print("[4/4] Construyendo tabla de anotación...")
    salida_tabla = os.path.join(dir_salida, "genes_counts_anotados.csv")
    salida_sin_simbolo = os.path.join(dir_salida, "genes_sin_simbolo.txt")
    salida_literatura_presente = os.path.join(dir_salida, "genes_literatura_en_counts.csv")

    n_con_simbolo = 0
    n_sin_simbolo = 0
    n_no_en_gtf = 0
    n_literatura_presente = 0

    with open(salida_tabla, "w", newline="") as f_tab, \
         open(salida_sin_simbolo, "w") as f_sin, \
         open(salida_literatura_presente, "w", newline="") as f_lit:

        escritor = csv.writer(f_tab)
        escritor.writerow(["ensembl_id", "simbolo", "biotipo", "en_gtf", "en_literatura", "simbolo_humano_lit"])

        escritor_lit = csv.writer(f_lit)
        escritor_lit.writerow(["ensembl_id", "simbolo", "simbolo_humano_lit"])

        for gene_id in ids_counts:
            en_gtf = gene_id in id_a_simbolo
            simbolo = id_a_simbolo.get(gene_id) if en_gtf else None
            biotipo = id_a_biotipo.get(gene_id) if en_gtf else None
            en_lit = gene_id in id_a_lit
            simbolo_humano = id_a_lit.get(gene_id, "")

            if not en_gtf:
                n_no_en_gtf += 1
            elif simbolo:
                n_con_simbolo += 1
            else:
                n_sin_simbolo += 1
                f_sin.write(gene_id + "\n")

            if en_lit:
                n_literatura_presente += 1
                escritor_lit.writerow([gene_id, simbolo, simbolo_humano])

            escritor.writerow([gene_id, simbolo, biotipo, en_gtf, en_lit, simbolo_humano])

    print("\n============================================================")
    print(f"Total genes en counts:              {len(ids_counts)}")
    print(f"  Con símbolo identificado:          {n_con_simbolo}")
    print(f"  En GTF pero SIN símbolo (gene_name vacío): {n_sin_simbolo}")
    print(f"  NO encontrados en el GTF:           {n_no_en_gtf}")
    print(f"Genes de literatura presentes en counts: {n_literatura_presente} / {len(id_a_lit)}")
    print("============================================================")
    print(f"\nArchivos generados en {dir_salida}:")
    print(f"  - genes_counts_anotados.csv      (tabla completa)")
    print(f"  - genes_sin_simbolo.txt          (IDs en GTF sin gene_name, revisar manualmente)")
    print(f"  - genes_literatura_en_counts.csv (subset: genes de literatura presentes en tus datos)")


if __name__ == "__main__":
    main()
