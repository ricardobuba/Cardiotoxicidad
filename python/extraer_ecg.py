# -*- coding: utf-8 -*-
"""
Extracción de la curva del ECG a partir de las imágenes de ecocardiografía.

Pipeline (Capítulo de extracción de la curva):
  1. Recorte de la región de interés, ajustado a mano imagen a imagen
     (config_recortes, en config_recortes.py).
  2. Segmentación del verde en HSV -> máscara binaria.
  3. Mediana por columna (huecos rellenados por interpolación lineal entre
     columnas vecinas) -> señal unidimensional; inversión del eje de píxel.
  4. Unificación de la polaridad: en parte de las imágenes el complejo QRS es
     una deflexión negativa. Se voltea la señal QUE SE GUARDA (no solo la que
     se usa para detectar), de modo que el pico R sea siempre positivo.
  5. Detección de los dos complejos QRS: prominencia mínima del 10 % del rango
     y separación mínima del 75 % de la anchura del recorte.
  6. El tramo entre ambos es un ciclo cardíaco -> spline penalizado, con el
     parámetro de suavizado λ elegido por validación cruzada generalizada ->
     1001 puntos equiespaciados en [0, 1].
  7. Centrado de cada curva. La AMPLITUD NO SE NORMALIZA: las diferencias de
     amplitud entre pacientes son potencialmente informativas.

Salidas (en RESULTADOS_ECG/):
  senales_ecg_full.csv    -> curvas válidas   (filename, num, CTRCD, time, p0..p1000)
  senales_malas_full.csv  -> señales descartadas

Uso:
  python python/extraer_ecg.py              # extracción
  python python/extraer_ecg.py --figuras    # extracción + figuras del capítulo
"""
import os, sys
import numpy as np
import pandas as pd
import cv2
from scipy.signal import find_peaks
from scipy.interpolate import make_smoothing_spline

# Los recortes manuales viven ahora en su propio módulo (config_recortes.py),
# copia verbatim de la celda de procesado2.ipynb. Antes se leían haciendo
# exec() de una celda del notebook, una dependencia frágil ya eliminada.
from config_recortes import config_recortes

# Rutas y parámetros del ECG centralizados en config.py (mismos valores).
from config import (
    ROOT, IMAGES_DIR as INPUT, RESULTADOS_ECG_DIR as OUT, FIG_DIR as FIGS,
    NPTS, VERDE_BAJO, VERDE_ALTO, PROM_QRS, DIST_QRS, DEFAULT_CROP,
    CSV_CLINICAL, CSV_SEP, CSV_DECIMAL,
)
os.makedirs(OUT, exist_ok=True)


# --- 1. Recortes manuales (definidos en config_recortes.py) ----------------
print("config_recortes cargado:", len(config_recortes), "entradas")

clin = pd.read_csv(CSV_CLINICAL, sep=CSV_SEP, decimal=CSV_DECIMAL)


# --- 2. Pipeline ----------------------------------------------------------
def mascara(img_bgr: np.ndarray) -> np.ndarray:
    """Segmenta el verde de la traza del ECG en una imagen BGR.

    Convierte a HSV, suaviza con un filtro gaussiano 5x5, umbraliza en el
    rango de verdes y limpia con un filtro de mediana.

    Args:
        img_bgr: Imagen (o recorte) en formato BGR de OpenCV.

    Returns:
        Máscara binaria (uint8) con los píxeles de la traza a 255.
    """
    hsv = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2HSV)
    blur = cv2.GaussianBlur(hsv, (5, 5), 0)
    return cv2.medianBlur(cv2.inRange(blur, VERDE_BAJO, VERDE_ALTO), 3)


def traza(mask: np.ndarray) -> np.ndarray:
    """Mediana de los píxeles activos en cada columna; los huecos (columnas
    sin píxeles activos) se rellenan por interpolación lineal entre las
    columnas vecinas con dato. Después se invierte el eje de píxel."""
    signal = np.full(mask.shape[1], np.nan)
    for col in range(mask.shape[1]):
        idx = np.where(mask[:, col] > 0)[0]
        if len(idx) > 0:
            signal[col] = np.median(idx)
    xs = np.arange(mask.shape[1])
    valido = ~np.isnan(signal)
    if valido.any():
        # Interpolación lineal: cada hueco toma el valor interpolado entre las
        # columnas vecinas con dato (para un hueco de una sola columna, la media
        # de la anterior y la siguiente); los extremos sin dato se completan con
        # el valor válido más cercano.
        signal = np.interp(xs, xs[valido], signal[valido])
    else:
        signal = np.zeros(mask.shape[1], dtype=float)
    return np.max(signal) - signal


def procesar(num: int):
    """Extrae la curva del ECG de la imagen ``num`` siguiendo el pipeline.

    Args:
        num: Número de imagen/paciente (1..270).

    Returns:
        Un ``dict`` con la curva discretizada y los productos intermedios
        (señal cruda, picos, recorte, etc.), o ``None`` si la imagen no
        existe.
    """
    img = cv2.imread(os.path.join(INPUT, f"image ({num}).png"))
    if img is None:
        return None
    p_sup, p_inf, p_izq, p_der = config_recortes.get(num, DEFAULT_CROP)
    alto, ancho = img.shape[:2]
    crop = img[int(alto * p_sup):int(alto * p_inf),
               int(ancho * p_izq):int(ancho * p_der)]

    signal = traza(mascara(crop))

    # Unificación de la polaridad
    sc = signal - np.mean(signal)
    invertida = abs(np.max(sc)) <= abs(np.min(sc))
    if invertida:
        signal = -signal
        sc = -sc

    largo = len(signal)
    rango = np.max(sc) - np.min(sc)
    picos, _ = find_peaks(sc, prominence=rango * PROM_QRS, distance=largo * DIST_QRS)

    valida = len(picos) >= 2
    rec = signal[picos[0]:picos[-1]] if valida else signal
    x_old = np.linspace(0, 1, len(rec))
    x_new = np.linspace(0, 1, NPTS)
    # Spline penalizado con lambda elegido por validacion cruzada generalizada
    # (Seccion de metodologia). NO usar UnivariateSpline con su s por defecto:
    # scipy toma s = len(x), que permite 1 px de error cuadratico medio y, sobre
    # una traza de ~12 px de amplitud, aplana el QRS (perdida mediana del 17 %
    # de la amplitud, y hasta el 87 % en las trazas mas pequenas).
    curva = make_smoothing_spline(x_old, rec)(x_new)
    curva = curva - np.mean(curva)      # centrado; la amplitud NO se toca

    return dict(num=num, curva=curva, valida=valida, picos=picos,
                signal=signal, invertida=invertida,
                recorte=(p_sup, p_inf, p_izq, p_der), crop=crop)


# --- 3. Extracción --------------------------------------------------------
filenames = sorted([f for f in os.listdir(INPUT) if f.lower().endswith(".png")],
                   key=lambda x: int(x.split("(")[1].split(")")[0]))
nums = [int(f.split("(")[1].split(")")[0]) for f in filenames]

buenas, malas, n_invertidas = [], [], 0
for num in nums:
    r = procesar(num)
    if r is None:
        continue
    n_invertidas += int(r["invertida"])
    (buenas if r["valida"] else malas).append(r)

print(f"Señales buenas: {len(buenas)} | malas: {len(malas)}")
print(f"Curvas con QRS negativo (polaridad corregida): {n_invertidas} de {len(nums)}")


# --- 4. Guardado ----------------------------------------------------------
def to_df(lst):
    rows = []
    for r in lst:
        row = {"filename": f"image ({r['num']}).png", "num": r["num"]}
        c = clin.iloc[r["num"] - 1]
        row["CTRCD"] = int(c["CTRCD"])
        row["time"] = int(c["time"])
        for j, v in enumerate(r["curva"]):
            row[f"p{j}"] = v
        rows.append(row)
    return pd.DataFrame(rows)


# Copia de seguridad de la versión previa a la corrección de polaridad.
for f in ("senales_ecg_full.csv", "senales_malas_full.csv"):
    src = os.path.join(OUT, f)
    bak = os.path.join(OUT, f.replace(".csv", "_SIN_CORREGIR.csv"))
    if os.path.exists(src) and not os.path.exists(bak):
        os.replace(src, bak)
        print(f"[backup] {f} -> {os.path.basename(bak)}")

df_b, df_m = to_df(buenas), to_df(malas)
df_b.to_csv(os.path.join(OUT, "senales_ecg_full.csv"), index=False)
df_m.to_csv(os.path.join(OUT, "senales_malas_full.csv"), index=False)
print("Guardado senales_ecg_full.csv", df_b.shape)
print("CTRCD en buenas:", df_b["CTRCD"].value_counts().to_dict())


# --- 5. Figuras del capítulo (opcional) -----------------------------------
if "--figuras" in sys.argv:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIGS, exist_ok=True)
    GREEN, PINK, BLUE = "#1a7f37", "#c0397a", "#1f4e79"
    xg = np.linspace(0, 1, NPTS)

    # (1) Etapas del pipeline sobre una imagen de ejemplo
    ej = buenas[0]
    img = cv2.imread(os.path.join(INPUT, f"image ({ej['num']}).png"))
    ps, pi, pz, pd_ = ej["recorte"]
    h, w = img.shape[:2]
    fig = plt.figure(figsize=(11, 7.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1, 1], hspace=.5, wspace=.22)
    a = fig.add_subplot(gs[0, :]); a.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
    a.add_patch(plt.Rectangle((w * pz, h * ps), w * (pd_ - pz), h * (pi - ps),
                              edgecolor="red", facecolor="none", lw=1.8))
    a.set_title("(a) Imagen original con la región de interés delimitada", fontsize=10); a.axis("off")
    a = fig.add_subplot(gs[1, 0]); a.imshow(cv2.cvtColor(ej["crop"], cv2.COLOR_BGR2RGB), aspect="auto")
    a.set_title("(b) Recorte de la región de interés", fontsize=10); a.axis("off")
    a = fig.add_subplot(gs[1, 1]); a.imshow(mascara(ej["crop"]), cmap="gray", aspect="auto")
    a.set_title("(c) Máscara binaria tras la segmentación en HSV", fontsize=10); a.axis("off")
    a = fig.add_subplot(gs[2, 0]); a.plot(ej["signal"], color=GREEN, lw=1)
    a.plot(ej["picos"], ej["signal"][ej["picos"]], "v", color=PINK, ms=8)
    a.set_title("(d) Traza extraída y complejos QRS detectados", fontsize=10)
    a.set_xlabel("columna de píxel", fontsize=8); a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[2, 1]); a.plot(xg, ej["curva"], color=BLUE, lw=1.2)
    a.set_title(f"(e) Ciclo normalizado y discretizado en {NPTS} puntos", fontsize=10)
    a.set_xlabel("tiempo normalizado", fontsize=8); a.tick_params(labelsize=7)
    fig.savefig(os.path.join(FIGS, "curva_pipeline.pdf"), bbox_inches="tight"); plt.close(fig)

    # (2) Necesidad del recorte manual
    full = mascara(img[int(h * ps):int(h * pi), :])
    ys, xs = np.where(full > 0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].imshow(mascara(img), cmap="gray")
    marc = (ys < h * ps) | (ys > h * pi)
    yy, xx = np.where(mascara(img) > 0)
    fuera = (yy < h * ps) | (yy > h * pi)
    ax[0].scatter(xx[fuera], yy[fuera], s=14, facecolors="none", edgecolors=PINK, lw=.9,
                  label=f"marcadores del ecógrafo ({int(fuera.sum())} px)")
    ax[0].axhspan(h * ps, h * pi, color="#e07b00", alpha=.15)
    ax[0].add_patch(plt.Rectangle((w * pz, h * ps), w * (pd_ - pz), h * (pi - ps),
                                  edgecolor="red", facecolor="none", lw=2))
    ax[0].legend(fontsize=8, loc="lower left")
    ax[0].set_title("(a) Segmentación del verde sobre la imagen completa.\n"
                    "La banda naranja es la traza; los círculos, marcadores del\n"
                    "equipo en el mismo verde. En rojo, el recorte manual.", fontsize=9)
    ax[0].axis("off")
    ax[1].imshow(mascara(ej["crop"]), cmap="gray", aspect="auto")
    ax[1].set_title("(b) Máscara dentro del recorte: sin marcadores y con\n"
                    "exactamente dos complejos QRS consecutivos.", fontsize=9)
    ax[1].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "curva_roi.pdf"), bbox_inches="tight"); plt.close(fig)

    # (3) Muestra de trazas con sus complejos QRS
    rs = np.random.RandomState(7)
    sel = rs.choice(len(buenas), 9, replace=False)
    fig, axes = plt.subplots(3, 3, figsize=(12, 6.5))
    for a, i in zip(axes.ravel(), sel):
        r = buenas[i]
        a.plot(r["signal"], color=GREEN, lw=1)
        a.plot(r["picos"], r["signal"][r["picos"]], "v", color=PINK, ms=8)
        a.set_title(f"Paciente {r['num']}", fontsize=9); a.tick_params(labelsize=7)
        a.set_xlabel("columna de píxel", fontsize=7)
    fig.suptitle("Trazas del ECG extraídas de nueve pacientes, con los dos complejos QRS detectados",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(os.path.join(FIGS, "curva_muestra.pdf"), bbox_inches="tight"); plt.close(fig)

    # (4) Discretización: el spline suaviza el escalonado
    rec = ej["signal"][ej["picos"][0]:ej["picos"][-1]]
    rec_c = rec - rec.mean()
    x_old = np.linspace(0, 1, len(rec_c))
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 2, hspace=.45, wspace=.25)
    a = fig.add_subplot(gs[0, :])
    a.plot(ej["signal"], color=GREEN, lw=1.2, label="traza cruda (un valor por columna)")
    a.plot(ej["picos"], ej["signal"][ej["picos"]], "v", color=PINK, ms=9, label="complejos QRS")
    a.axvspan(ej["picos"][0], ej["picos"][-1], color=BLUE, alpha=.10, label="tramo conservado")
    a.set_title(f"(a) Traza extraída ({len(ej['signal'])} columnas); se conserva el tramo "
                f"entre los dos QRS ({len(rec)} columnas)", fontsize=10)
    a.legend(fontsize=8); a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[1, 0])
    a.plot(x_old, rec_c, "o-", color=GREEN, ms=2.5, lw=.8, label=f"puntos crudos ({len(rec_c)})")
    a.plot(xg, ej["curva"], color=PINK, lw=1.6, label=f"spline en {NPTS} puntos")
    a.set_title("(b) Reparametrización a [0,1] y suavizado", fontsize=10)
    a.legend(fontsize=8); a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[1, 1])
    m1 = (xg >= .35) & (xg <= .55); m2 = (x_old >= .35) & (x_old <= .55)
    a.plot(x_old[m2], rec_c[m2], "o-", color=GREEN, ms=5, lw=.9, label="puntos crudos")
    a.plot(xg[m1], ej["curva"][m1], color=PINK, lw=1.8, label="spline")
    a.set_title("(c) Detalle: el spline suaviza el escalonado que\nintroduce la mediana por columna",
                fontsize=10)
    a.legend(fontsize=8); a.tick_params(labelsize=7)
    fig.savefig(os.path.join(FIGS, "curva_interpolacion.pdf"), bbox_inches="tight"); plt.close(fig)

    # (5) Efecto de la unificación de la polaridad y el centrado
    E = np.array([r["curva"] for r in buenas])
    E0 = np.array([r["curva"] * (-1 if r["invertida"] else 1) for r in buenas])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for a, X, t in [(ax[0], E0, "(a) Sin unificar la polaridad: se mezclan\nlas dos orientaciones del complejo QRS"),
                    (ax[1], E,  "(b) Con la polaridad unificada y centrado\n(la amplitud no se modifica)")]:
        for i in range(len(X)):
            a.plot(xg, X[i], color="gray", lw=.4, alpha=.28)
        a.plot(xg, X.mean(0), color=PINK, lw=2.2, label="Curva media")
        a.axhline(0, color="k", lw=.6, ls=":")
        a.set_title(t, fontsize=10); a.set_xlabel("Fase del ciclo cardíaco", fontsize=9)
        a.set_ylabel("Amplitud (px)", fontsize=9); a.legend(fontsize=8); a.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "curva_polaridad.pdf"), bbox_inches="tight"); plt.close(fig)

    # (6) Ejemplos de señales descartadas
    fig, axes = plt.subplots(2, 3, figsize=(11, 5.2))
    for a, r in zip(axes.ravel(), malas[:6]):
        a.plot(xg, r["curva"], color="firebrick", lw=1)
        a.set_title(f"Paciente {r['num']}", fontsize=9); a.tick_params(labelsize=6)
    fig.suptitle("Ejemplos de trazas descartadas: no se detectan dos complejos QRS "
                 "que delimiten un ciclo completo", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(os.path.join(FIGS, "curva_descartes.pdf"), bbox_inches="tight"); plt.close(fig)

    # (7) Un descarte explicado: la traza es legible, pero el recorte no encuadra
    #     los dos QRS lo bastante separados como para superar el criterio.
    cand = []
    for r in malas:
        s = r["signal"]
        sc = s - s.mean()
        rg = sc.max() - sc.min()
        p, pr = find_peaks(sc, prominence=rg * 0.5, distance=25)   # solo QRS nítidos
        # Se exige además que el descarte sea claro (no un empate al filo del
        # criterio), para que el ejemplo ilustre y no confunda.
        if len(p) == 2 and (p[1] - p[0]) < 0.90 * DIST_QRS * len(s):
            cand.append((r, p, len(s), pr["prominences"].min() / rg))
    if cand:
        r, p, L, _ = max(cand, key=lambda t: t[3])   # el de QRS más nítidos
        s, sep, req = r["signal"], int(p[1] - p[0]), DIST_QRS * L
        fig, a = plt.subplots(figsize=(11, 4.6))
        a.plot(s, color=GREEN, lw=1.3)
        a.plot(p, s[p], "o", ms=11, mfc="none", mec=BLUE, mew=2,
               label="complejos QRS nítidos (2)")
        dom = p[int(np.argmax(s[p]))]
        a.plot(dom, s[dom], "*", ms=20, color=PINK,
               label="pico dominante: el algoritmo parte de aquí")
        a.axvspan(max(0, dom - req), min(L, dom + req), color=PINK, alpha=.13,
                  label=f"zona de exclusión (±{req:.0f} px = 75 % del recorte)")
        a.annotate("", xy=(p[0], s[p].min() * .3), xytext=(p[1], s[p].min() * .3),
                   arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.4))
        a.text((p[0] + p[1]) / 2, s[p].min() * .38,
               f"separación real: {sep} px  (hacen falta {req:.0f} px)",
               ha="center", fontsize=9, color=BLUE)
        a.set_xlim(0, L)
        a.set_title(f"Paciente {r['num']} (descartada). Los dos complejos QRS son inequívocos, "
                    f"pero están a {sep} px\ny el criterio exige {req:.0f} px: el segundo cae "
                    f"dentro de la zona de exclusión del primero.", fontsize=10)
        a.set_xlabel("columna de píxel"); a.set_ylabel("amplitud (px)")
        a.legend(fontsize=8, loc="lower right", framealpha=.95); a.tick_params(labelsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, "curva_descarte_ejemplo.pdf"), bbox_inches="tight")
        plt.close(fig)
        print(f"[figura] descarte explicado: paciente {r['num']} "
              f"(separación {sep} px, exigida {req:.0f} px)")

    print(f"[figuras] generadas en {FIGS}")
