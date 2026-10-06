# -*- coding: utf-8 -*-
"""Configuración centralizada de los scripts de Python del TFG.

Reúne en un único lugar las rutas, la paleta de colores, los parámetros del
ECG y los parámetros estadísticos que antes estaban repartidos (y a veces
duplicados) entre ``extraer_ecg.py`` y ``eda_analysis.py``.

IMPORTANTE
----------
Importar desde aquí **no cambia ningún resultado**: los valores son
exactamente los mismos que ya tenían los scripts. Este módulo solo centraliza
constantes; no altera la lógica, ni las figuras, ni las tablas.

La parte en R (``eda_funcional_superv.R``) mantiene sus propias constantes
equivalentes (``PINK``, ``BLUE``, factor ``1.5`` del fbplot, ``NHARM = 8``).
R no puede importar este módulo de Python; si cambias aquí un color o un
parámetro que también use el R, refléjalo a mano allí.
"""
import os

import numpy as np

# ---------------------------------------------------------------------------
# Rutas
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMAGES_DIR = os.path.join(ROOT, "IMAGES")
# datos/ agrupa los CSV originales del CHUAC.
DATOS_DIR = os.path.join(ROOT, "datos")
RESULTADOS_ECG_DIR = os.path.join(ROOT, "RESULTADOS_ECG")
# salidas/ es la carpeta que se sube a Overleaf (imaxes/eda,
# imaxes/curva y contido/tablas_eda).
SALIDAS_DIR = os.path.join(ROOT, "salidas")
FIG_DIR = os.path.join(SALIDAS_DIR, "figuras")
TAB_DIR = os.path.join(SALIDAS_DIR, "tablas")

# ---------------------------------------------------------------------------
# CSV de entrada
# ---------------------------------------------------------------------------
CSV_CLINICAL = os.path.join(DATOS_DIR, "BC_cardiotox_clinical_variables.csv")
CSV_FUNCTIONAL = os.path.join(DATOS_DIR, "BC_cardiotox_functional_variable.csv")
CSV_SEP = ";"
CSV_DECIMAL = ","
N_PATIENTS_EXPECTED = 531   # nº de pacientes del conjunto clínico (validación)

# ---------------------------------------------------------------------------
# Paleta de colores (plantilla UDC): 0 = no CTRCD, 1 = CTRCD
# ---------------------------------------------------------------------------
PINK = "#C5006E"
BLUE = "#4C9BD4"
GRAY = "#7a7a7a"
PAL = {0: BLUE, 1: PINK}

# ---------------------------------------------------------------------------
# Parámetros del ECG (extraer_ecg.py)
# ---------------------------------------------------------------------------
NPTS = 1001                            # nº de puntos de la curva discretizada
VERDE_BAJO = np.array([38, 100, 70])   # umbral HSV inferior del verde de la traza
VERDE_ALTO = np.array([89, 255, 255])  # umbral HSV superior del verde de la traza
PROM_QRS = 0.10      # prominencia mínima de los picos, fracción del rango
DIST_QRS = 0.75      # separación mínima entre picos, fracción de la anchura
DEFAULT_CROP = (0.42, 0.50, 0.24, 0.50)  # recorte (p_sup, p_inf, p_izq, p_der) por defecto

# ---------------------------------------------------------------------------
# Parámetros estadísticos y de representación (eda_analysis.py)
# ---------------------------------------------------------------------------
ALPHA = 0.05          # nivel de significación para marcar diferencias
IQR_WHISKER = 1.5     # factor de la regla del rango intercuartílico (atípicos)
HIST_BINS_CONT = 25   # nº de barras en los histogramas de variables continuas
HIST_BINS_TIME = 30   # nº de barras en el histograma del tiempo de seguimiento
PCA_COMPONENTS = 2    # componentes retenidas en el ACP
FIG_DPI = 120         # resolución de las figuras
