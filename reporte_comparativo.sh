#!/bin/bash

#===============================================================================
# RECOPILACIÓN EXTRA Y REPORTE COMPARATIVO GLOBAL
# ==============================================================================

DIR_RESULTADOS="/data/Database/Resultados"
PYTHON_SCRIPT="/data/Database/comparador_qc.py"

echo "======================================================================"
echo "📊 GENERANDO REPORTE COMPARATIVO CSV (CRUDO VS TRIMMED)"
echo "======================================================================"

# El script escribe dinámicamente el archivo de python con las rutas corregidas a 'resultados'
cat << 'EOF' > "$PYTHON_SCRIPT"
import pandas as pd
import os

ruta_crudos = "/data/Database/Resultados/CENTRAL_MULTIQC/Crudos_Global/multiqc_data/multiqc_fastqc.txt"
ruta_trimmed = "/data/Database/Resultados/CENTRAL_MULTIQC/Trimmed_Global/multiqc_data/multiqc_fastqc.txt"
ruta_salida = "/data/Database/Resultados/reporte_evolucion_QC.csv"

if not os.path.exists(ruta_crudos) or not os.path.exists(ruta_trimmed):
    print("❌ Error: No se encontraron los archivos de datos globales de MultiQC.")
    exit()

# Las 3 métricas de interés seleccionadas para el análisis
metricas = ['per_base_sequence_quality', 'adapter_content', 'sequence_length_distribution']

# Cargar datos seleccionando las columnas clave
df_crudos = pd.read_csv(ruta_crudos, sep='\t')[['Sample'] + metricas]
df_trimmed = pd.read_csv(ruta_trimmed, sep='\t')[['Sample'] + metricas]

# Renombrar columnas para identificar origen y destino
df_crudos.columns = ['Sample', 'Calidad_Base_Orig', 'Adaptadores_Orig', 'Longitud_Dist_Orig']
df_trimmed.columns = ['Sample', 'Calidad_Base_Final', 'Adaptadores_Final', 'Longitud_Dist_Final']

# Limpieza de nombres de muestras generados por Trim Galore para poder emparejarlas
df_trimmed['Sample'] = df_trimmed['Sample'].str.replace(r'_val_[12]$', '', regex=True)
df_trimmed['Sample'] = df_trimmed['Sample'].str.replace(r'_trimmed$', '', regex=True)

# Unir crudos y procesados en una sola tabla matriz
df_comparativo = pd.merge(df_crudos, df_trimmed, on='Sample', how='inner')

# Diccionario de peso matemático para los semáforos de FastQC
valores_qc = {'pass': 2, 'warn': 1, 'fail': 0}

# ------------------------------------------------------------------------------
# FUNCIONES LÓGICAS DE DECISIÓN
# ------------------------------------------------------------------------------

# 1. Función para Calidad por Base y Adaptadores (Donde subir o mantenerse arriba es lo ideal)
def evaluar_mejora_estandar(orig_col, final_col):
    def logica(fila):
        val_orig = valores_qc.get(fila[orig_col], -1)
        val_final = valores_qc.get(fila[final_col], -1)
        if val_final > val_orig:
            return 'Mejora'
        elif val_final == val_orig:
            return 'Igual'
        else:
            return 'Empeora'
    return logica

# 2. Función para Distribución de Longitud (Donde pasar a 'warn' es lo normal tras recortar)
def evaluar_longitud(fila):
    orig = fila['Longitud_Dist_Orig']
    final = fila['Longitud_Dist_Final']
    
    if orig == 'pass' and final == 'warn':
        return 'Recortado Exitoso'
    elif orig == 'pass' and final == 'pass':
        return 'No requirió recortes'
    elif orig == 'warn' and final == 'warn':
        return 'No requirió recortes'
    elif final == 'fail':
        return 'Pérdida Crítica'
    else:
        return f'Cambio de {orig} a {final}'

# Aplicar las funciones lógicas a las columnas correspondientes
df_comparativo['Estado_Calidad_Base'] = df_comparativo.apply(evaluar_mejora_estandar('Calidad_Base_Orig', 'Calidad_Base_Final'), axis=1)
df_comparativo['Estado_Adaptadores'] = df_comparativo.apply(evaluar_mejora_estandar('Adaptadores_Orig', 'Adaptadores_Final'), axis=1)
df_comparativo['Estado_Distribucion_Longitud'] = df_comparativo.apply(evaluar_longitud, axis=1)

# Ordenar las columnas para que el archivo final sea cómodo de leer
columnas_ordenadas = [
    'Sample', 
    'Calidad_Base_Orig', 'Calidad_Base_Final', 'Estado_Calidad_Base',
    'Adaptadores_Orig', 'Adaptadores_Final', 'Estado_Adaptadores',
    'Longitud_Dist_Orig', 'Longitud_Dist_Final', 'Estado_Distribucion_Longitud'
]
df_comparativo = df_comparativo[columnas_ordenadas]

# Guardar a archivo final CSV
df_comparativo.to_csv(ruta_salida, index=False, encoding='utf-8')
print(f"✅ Reporte comparativo multi-métrica generado con éxito en: {ruta_salida}")
EOF

# Ejecutar el script generado pasándole el entorno de python3
python3 "$PYTHON_SCRIPT"

# Limpieza del script temporal de Python para dejar el directorio impecable
rm -f "$PYTHON_SCRIPT"

echo "¡ARCHIVO COMPARATIVO GENERADO EXITOSAMENTE!"
