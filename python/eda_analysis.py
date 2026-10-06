# -*- coding: utf-8 -*-
"""Análisis exploratorio y descriptivo del conjunto de datos *BC_cardiotox*.

Genera todas las figuras (``imaxes/eda``) y las tablas LaTeX
(``contido/tablas_eda``) que se utilizan en el capítulo de Análisis
descriptivo de la memoria.

La parte de Kaplan--Meier y el análisis funcional (curvas ECG + FPCA) se
generan aparte en R (``eda_funcional_superv.R``, con *survival*/*survminer*
y ``fda::pca.fd``). Aquí solo se producen las figuras y tablas del EDA
"clásico".

Uso
---
    python python/eda_analysis.py

Salidas
-------
``salidas/figuras/*.pdf`` y ``salidas/tablas/*.tex``.

Nota
----
Este refactor conserva **exactamente** las mismas salidas numéricas y
gráficas que la versión previa: solo cambian la organización del código, la
documentación, la validación de entrada y el registro de mensajes.
"""
from __future__ import annotations

import logging
import os
import re
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

log = logging.getLogger("eda_analysis")


# ----------------------------------------------------------------------
# 0. Configuración
# ----------------------------------------------------------------------
# Rutas, paleta y parámetros centralizados en config.py (mismos valores).
from config import (  # noqa: E402
    ROOT,  # noqa: F401
    FIG_DIR as FIG, TAB_DIR as TAB,
    CSV_CLINICAL, CSV_FUNCTIONAL, CSV_SEP, CSV_DECIMAL,
    N_PATIENTS_EXPECTED,
    PINK, BLUE, GRAY, PAL,  # noqa: F401
    ALPHA, IQR_WHISKER, HIST_BINS_CONT, HIST_BINS_TIME,
    PCA_COMPONENTS, FIG_DPI,
)

sns.set_theme(style="whitegrid", context="paper")
plt.rcParams.update({"figure.dpi": FIG_DPI, "savefig.bbox": "tight",
                     "axes.titlesize": 11, "font.size": 10})

# Grupos de variables
CONT = ["age", "weight", "height", "heart_rate", "LVEF", "PWT",
        "LAd", "LVDd", "LVSd"]
BIN = ["heart_rhythm", "AC", "antiHER2", "HTA", "DL", "DM", "smoker",
       "exsmoker", "ACprev", "antiHER2prev", "RTprev", "CIprev",
       "ICMprev", "ARRprev", "VALVprev", "cxvalv"]
TARGET, TIME = "CTRCD", "time"

# Etiquetas legibles (es)
LAB = {
    "age": "Edad (años)", "weight": "Peso (kg)", "height": "Estatura (cm)",
    # OJO: heart_rate es la FRECUENCIA cardíaca; el ritmo cardíaco es
    # heart_rhythm, otra variable del mismo fichero.
    "heart_rate": "Frecuencia cardíaca (lpm)", "LVEF": "FEVI (%)",
    "PWT": "PWT (cm)", "LAd": "LAd (cm)", "LVDd": "LVDd (cm)",
    "LVSd": "LVSd (cm)", "time": "Tiempo de seguimiento (días)",
    "heart_rhythm": "Fibrilación auricular", "AC": "Antraciclinas",
    "antiHER2": "Terapia anti-HER2", "HTA": "Hipertensión",
    "DL": "Dislipidemia", "DM": "Diabetes mellitus", "smoker": "Fumadora",
    "exsmoker": "Exfumadora", "ACprev": "Antraciclinas previas",
    "antiHER2prev": "Anti-HER2 previa", "RTprev": "Radioterapia previa",
    "CIprev": "Insuf. cardíaca previa", "ICMprev": "Miocardiopatía isq. previa",
    "ARRprev": "Arritmia previa", "VALVprev": "Valvulopatía previa",
    "cxvalv": "Cirugía valvular",
}


# ----------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------
_DEC = re.compile(r"(?<=[0-9])\.(?=[0-9])")


def coma_decimal(s: str) -> str:
    """Sustituye el punto decimal por ``{,}`` en un texto destinado a LaTeX.

    Se emplea ``{,}`` y no una coma suelta porque en modo matemático esta
    última se compone como signo de puntuación. El patrón solo actúa entre
    dígitos, de modo que no toca ningún comando de LaTeX.

    Args:
        s: Texto de una celda, cabecera, pie o *caption*.

    Returns:
        El texto con los puntos decimales sustituidos.
    """
    return _DEC.sub("{,}", s)


def texsafe(s: object) -> str:
    """Escapa los caracteres especiales de LaTeX de una cadena.

    Args:
        s: Valor cualquiera; se convierte a ``str`` antes de escapar.

    Returns:
        La cadena con ``%``, ``&`` y ``_`` escapados para LaTeX.
    """
    return (str(s).replace("%", r"\%").replace("&", r"\&")
            .replace("_", r"\_"))


def savefig(fig: plt.Figure, name: str) -> None:
    """Guarda una figura en ``FIG`` como PDF y la cierra.

    Args:
        fig: Figura de Matplotlib a guardar.
        name: Nombre del fichero de salida (con extensión).
    """
    path = os.path.join(FIG, name)
    fig.savefig(path)
    plt.close(fig)
    log.info("  [fig] %s", name)


def write_table(name: str, header: Sequence[str], rows: Sequence[Sequence[object]],
                caption: str, label: str, colspec: str | None = None,
                note: str | None = None) -> None:
    """Escribe una tabla LaTeX con el estilo de la plantilla UDC.

    Args:
        name: Nombre del fichero ``.tex`` de salida.
        header: Encabezados de columna.
        rows: Filas de la tabla (cada fila, una secuencia de celdas).
        caption: Texto del ``\\caption``.
        label: Etiqueta del ``\\label``.
        colspec: Especificación de columnas de ``tabular``. Si es ``None``,
            se usa ``l`` para la primera columna y ``c`` para el resto.
        note: Nota opcional en letra pequeña bajo la tabla.
    """
    ncol = len(header)
    if colspec is None:
        colspec = "l" + "c" * (ncol - 1)
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \rowcolors{2}{white}{udcgray!25}")
    # Escala la tabla al ancho de texto solo si la excede (las estrechas
    # conservan su tamaño natural y quedan centradas).
    lines.append(r"  \resizebox{\ifdim\width>\linewidth\linewidth\else"
                 r"\width\fi}{!}{%")
    lines.append(r"  \begin{tabular}{%s}" % colspec)
    lines.append(r"    \rowcolor{udcpink!25}")
    lines.append("    " + " & ".join(r"\textbf{%s}" % coma_decimal(texsafe(h))
                                     for h in header) + r" \\ \hline")
    for r in rows:
        lines.append("    " + " & ".join(coma_decimal(str(c)) for c in r)
                     + r" \\")
    lines.append(r"  \end{tabular}}")
    if note:
        lines.append(r"  \\[2pt] {\footnotesize %s}" % coma_decimal(note))
    lines.append(r"  \caption{%s}" % coma_decimal(caption))
    lines.append(r"  \label{%s}" % label)
    lines.append(r"\end{table}")
    with open(os.path.join(TAB, name), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log.info("  [tab] %s", name)


def load_clinical() -> pd.DataFrame:
    """Carga el CSV clínico y comprueba su integridad básica.

    Returns:
        El ``DataFrame`` con las variables clínicas.

    Raises:
        FileNotFoundError: Si no se encuentra el CSV clínico.
        ValueError: Si faltan columnas esperadas en el CSV.
    """
    if not os.path.isfile(CSV_CLINICAL):
        raise FileNotFoundError(
            f"No se encuentra el CSV clínico: {CSV_CLINICAL}")
    df = pd.read_csv(CSV_CLINICAL, sep=CSV_SEP, decimal=CSV_DECIMAL)
    faltan = [c for c in (CONT + BIN + [TARGET, TIME]) if c not in df.columns]
    if faltan:
        raise ValueError(
            "Faltan columnas esperadas en el CSV clínico: " + ", ".join(faltan))
    if len(df) != N_PATIENTS_EXPECTED:
        log.warning("El CSV clínico tiene %d filas (se esperaban %d).",
                    len(df), N_PATIENTS_EXPECTED)
    return df


def main() -> None:
    """Ejecuta el análisis exploratorio completo y genera figuras y tablas."""
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Carga y limpieza
    # ------------------------------------------------------------------
    log.info("== Carga y limpieza ==")
    df = load_clinical()
    if not os.path.isfile(CSV_FUNCTIONAL):
        raise FileNotFoundError(
            f"No se encuentra el CSV funcional: {CSV_FUNCTIONAL}")
    fn = pd.read_csv(CSV_FUNCTIONAL, sep=CSV_SEP, decimal=CSV_DECIMAL)
    log.info("clinical: %s | functional: %s", df.shape, fn.shape)

    N = len(df)
    log.info("N = %d | eventos CTRCD = %d (%.1f%%)",
             N, int(df[TARGET].sum()), 100 * df[TARGET].mean())

    # ------------------------------------------------------------------
    # 2. Valores faltantes
    # ------------------------------------------------------------------
    log.info("== Valores faltantes ==")
    miss = df.isna().sum()
    miss = miss[miss > 0].sort_values(ascending=False)
    miss_pct = (miss / N * 100).round(1)
    log.info("\n%s", pd.concat([miss, miss_pct], axis=1, keys=["n", "%"]))

    fig, ax = plt.subplots(figsize=(7, 5))
    order = miss_pct.index
    ax.barh([LAB.get(v, v) for v in order], miss_pct.values, color=PINK)
    ax.invert_yaxis()
    ax.set_xlabel("Porcentaje de valores faltantes (%)")
    ax.set_title("Valores faltantes por variable")
    for i, v in enumerate(miss_pct.values):
        ax.text(v + 0.1, i, f"{v}%", va="center", fontsize=8)
    savefig(fig, "eda_missing.pdf")

    # Faltantes desglosados por clase de cardiotoxicidad (sugerencia de la
    # tutora): permite ver si los faltantes se concentran en un grupo, p. ej.
    # si los faltantes de "exfumadora" corresponden solo a pacientes sin CTRCD.
    miss0 = df.loc[df[TARGET] == 0].isna().sum()
    miss1 = df.loc[df[TARGET] == 1].isna().sum()
    N0 = int((df[TARGET] == 0).sum())
    N1 = int((df[TARGET] == 1).sum())
    rows = [(texsafe(LAB.get(v, v)),
             int(miss0[v]), int(miss1[v]),
             int(miss[v]), f"{miss_pct[v]:.1f}")
            for v in order]
    write_table("tab_missing.tex",
                ["Variable", "Sin CTRCD", "Con CTRCD", "Total", "% total"],
                rows,
                "Variables con valores faltantes en el conjunto clínico, "
                "desglosadas según el estado de cardiotoxicidad. El porcentaje "
                "se calcula sobre el total de pacientes.",
                "tab:missing",
                note=(f"Tamaños de grupo: sin CTRCD $n={N0}$, "
                      f"con CTRCD $n={N1}$."))

    # ------------------------------------------------------------------
    # 3. Variable objetivo y tiempo de seguimiento
    # ------------------------------------------------------------------
    log.info("== Variable objetivo / tiempo ==")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    vc = df[TARGET].value_counts().sort_index()
    bars = axes[0].bar(["No CTRCD\n(censurado)", "CTRCD\n(evento)"], vc.values,
                       color=[PAL[0], PAL[1]])
    axes[0].set_ylabel("Número de pacientes")
    axes[0].set_title("Distribución de la cardiotoxicidad (CTRCD)")
    for b, v in zip(bars, vc.values):
        axes[0].text(b.get_x() + b.get_width()/2, v + 5,
                     f"{v}\n({100*v/N:.1f}%)", ha="center", fontsize=9)

    for k in [0, 1]:
        axes[1].hist(df.loc[df[TARGET] == k, TIME], bins=HIST_BINS_TIME, alpha=0.65,
                     color=PAL[k], label=("Evento" if k else "Censurado"))
    axes[1].set_xlabel("Tiempo de seguimiento (días)")
    axes[1].set_ylabel("Frecuencia")
    axes[1].set_title("Tiempo hasta evento / censura")
    axes[1].legend()
    savefig(fig, "eda_target_time.pdf")

    # ------------------------------------------------------------------
    # 4. Univariante: variables continuas
    # ------------------------------------------------------------------
    log.info("== Univariante continuas ==")
    desc_rows = []
    for v in CONT + [TIME]:
        s = df[v].dropna()
        desc_rows.append((
            texsafe(LAB[v]), len(s), f"{s.mean():.2f}", f"{s.std():.2f}",
            f"{s.min():.1f}", f"{s.quantile(.25):.1f}", f"{s.median():.1f}",
            f"{s.quantile(.75):.1f}", f"{s.max():.1f}",
            f"{stats.skew(s):.2f}"))
    write_table("tab_desc_cont.tex",
                ["Variable", "n", "Media", "DE", "Mín", "Q1", "Mediana",
                 "Q3", "Máx", "Asim."],
                desc_rows,
                "Estadísticos descriptivos de las variables continuas.",
                "tab:desc_cont")

    # histogramas
    fig, axes = plt.subplots(3, 3, figsize=(12, 9))
    for ax, v in zip(axes.ravel(), CONT):
        s = df[v].dropna()
        ax.hist(s, bins=HIST_BINS_CONT, color=PAL[0], edgecolor="white", alpha=0.85)
        ax.axvline(s.mean(), color=PINK, ls="--", lw=1.3, label="Media")
        ax.axvline(s.median(), color="black", ls=":", lw=1.3, label="Mediana")
        ax.set_title(LAB[v])
        ax.set_ylabel("Frec.")
    axes.ravel()[0].legend(fontsize=7)
    fig.suptitle("Distribución de las variables continuas", y=1.01,
                 fontsize=13)
    fig.tight_layout()
    savefig(fig, "eda_hist_cont.pdf")

    # boxplots
    fig, axes = plt.subplots(3, 3, figsize=(12, 9))
    for ax, v in zip(axes.ravel(), CONT):
        ax.boxplot(df[v].dropna(), vert=True, patch_artist=True,
                   boxprops=dict(facecolor=PAL[0], alpha=0.7),
                   medianprops=dict(color=PINK, lw=2))
        ax.set_title(LAB[v])
        ax.set_xticks([])
    fig.suptitle("Diagramas de caja de las variables continuas (detección de "
                 "valores atípicos)", y=1.01, fontsize=13)
    fig.tight_layout()
    savefig(fig, "eda_box_cont.pdf")

    # conteo de outliers por regla IQR
    log.info("Outliers (regla %.1f*IQR):", IQR_WHISKER)
    for v in CONT:
        s = df[v].dropna()
        q1, q3 = s.quantile(.25), s.quantile(.75)
        iqr = q3 - q1
        lo, hi = q1 - IQR_WHISKER*iqr, q3 + IQR_WHISKER*iqr
        nout = ((s < lo) | (s > hi)).sum()
        log.info("  %-12s: %d outliers", v, nout)

    # ------------------------------------------------------------------
    # 5. Univariante: variables binarias
    # ------------------------------------------------------------------
    log.info("== Univariante binarias ==")
    bin_rows = []
    for v in BIN:
        s = df[v].dropna()
        n1 = int(s.sum())
        bin_rows.append((texsafe(LAB[v]), len(s), n1, f"{100*n1/len(s):.1f}"))
    write_table("tab_freq_bin.tex",
                ["Variable", "n", "n (factor presente)", "Prevalencia %"],
                bin_rows,
                "Frecuencia de las variables binarias (factores de riesgo y "
                "tratamientos).",
                "tab:freq_bin")

    fig, ax = plt.subplots(figsize=(8, 6))
    prev = pd.Series({v: 100*df[v].dropna().mean() for v in BIN}
                     ).sort_values()
    ax.barh([LAB[v] for v in prev.index], prev.values, color=PINK)
    ax.set_xlabel("Prevalencia (%)")
    ax.set_title("Prevalencia de las variables binarias")
    for i, v in enumerate(prev.values):
        ax.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=8)
    savefig(fig, "eda_prev_bin.pdf")

    # ------------------------------------------------------------------
    # 6. Bivariante: continuas vs CTRCD
    # ------------------------------------------------------------------
    log.info("== Bivariante continuas vs CTRCD ==")
    biv_rows = []
    for v in CONT:
        g0 = df.loc[df[TARGET] == 0, v].dropna()
        g1 = df.loc[df[TARGET] == 1, v].dropna()
        u, p = stats.mannwhitneyu(g0, g1, alternative="two-sided")
        star = "*" if p < ALPHA else ""

        # Mediana e IQR además de media y DT: el contraste es de rangos, así
        # que el resumen que le corresponde es el robusto. La media se conserva
        # porque es la que cita la prosa del capítulo.
        def _med(g):
            return (f"{g.median():.2f} "
                    f"[{g.quantile(0.25):.2f}--{g.quantile(0.75):.2f}]")

        biv_rows.append((
            texsafe(LAB[v]),
            _med(g0),
            f"{g0.mean():.2f} $\\pm$ {g0.std():.2f}",
            _med(g1),
            f"{g1.mean():.2f} $\\pm$ {g1.std():.2f}",
            (f"{p:.3f}" if p >= 0.001 else "$<$0.001") + star))
        log.info("  %-12s p=%.4f", v, p)
    write_table("tab_biv_cont.tex",
                ["Variable",
                 "No CTRCD: mediana [Q1--Q3]", "No CTRCD: media$\\pm$DT",
                 "CTRCD: mediana [Q1--Q3]", "CTRCD: media$\\pm$DT",
                 "p-valor"],
                biv_rows,
                "Comparación de las variables continuas según el estado de "
                "cardiotoxicidad (prueba U de Mann--Whitney). El contraste "
                "compara distribuciones por rangos, por lo que el resumen que "
                "le corresponde es la mediana con su rango intercuartílico; "
                "la media y la desviación típica se incluyen para facilitar "
                "la comparación con la literatura clínica.",
                "tab:biv_cont",
                colspec="lccccc",
                note="* indica diferencia significativa al nivel $\\alpha=0.05$.")

    # boxplots por grupo
    fig, axes = plt.subplots(3, 3, figsize=(12, 9))
    for ax, v in zip(axes.ravel(), CONT):
        data = [df.loc[df[TARGET] == 0, v].dropna(),
                df.loc[df[TARGET] == 1, v].dropna()]
        bp = ax.boxplot(data, patch_artist=True, labels=["No", "Sí"])
        for patch, k in zip(bp["boxes"], [0, 1]):
            patch.set_facecolor(PAL[k]); patch.set_alpha(0.7)
        for med in bp["medians"]:
            med.set_color("black")
        ax.set_title(LAB[v])
        ax.set_xlabel("CTRCD")
    fig.suptitle("Variables continuas según el estado de cardiotoxicidad",
                 y=1.01, fontsize=13)
    fig.tight_layout()
    savefig(fig, "eda_box_by_target.pdf")

    # ------------------------------------------------------------------
    # 7. Bivariante: binarias vs CTRCD (tasa de evento, OR, test)
    # ------------------------------------------------------------------
    log.info("== Bivariante binarias vs CTRCD ==")
    binbiv_rows = []
    or_data = []
    for v in BIN:
        sub = df[[v, TARGET]].dropna()
        a = ((sub[v] == 1) & (sub[TARGET] == 1)).sum()
        b = ((sub[v] == 1) & (sub[TARGET] == 0)).sum()
        c = ((sub[v] == 0) & (sub[TARGET] == 1)).sum()
        d = ((sub[v] == 0) & (sub[TARGET] == 0)).sum()
        table = np.array([[a, b], [c, d]])
        # tasa de evento en cada grupo
        r1 = 100*a/(a+b) if (a+b) else np.nan
        r0 = 100*c/(c+d) if (c+d) else np.nan
        # OR con corrección de Haldane--Anscombe si hay ceros. La fila se
        # marca: con una casilla nula el OR crudo es 0 o infinito y lo que se
        # tabula es una cota estabilizada, no una estimación del efecto.
        corregido = bool((table == 0).any())
        aa, bb, cc, dd = (table + 0.5).ravel() if corregido else table.ravel()
        OR = (aa*dd)/(bb*cc)
        # test: Fisher si alguna esperada < 5
        try:
            _, p_chi, _, expected = stats.chi2_contingency(table)
            use_fisher = (expected < 5).any()
        except Exception:
            use_fisher = True
        if use_fisher:
            _, p = stats.fisher_exact(table)
            test = "F"
        else:
            p = p_chi
            test = "$\\chi^2$"
        star = "*" if p < ALPHA else ""
        binbiv_rows.append((
            texsafe(LAB[v]),
            f"{r0:.1f}", f"{r1:.1f}",
            f"{OR:.2f}" + ("$^{\\dagger}$" if corregido else ""),
            (f"{p:.3f}" if p >= 0.001 else "$<$0.001") + star))
        or_data.append((LAB[v], OR, p))
        log.info("  %-12s tasa0=%.1f%% tasa1=%.1f%% OR=%.2f p=%.4f (%s)",
                 v, r0, r1, OR, p, test)
    write_table("tab_biv_bin.tex",
                ["Variable", "Tasa sin factor (%)", "Tasa con factor (%)", "OR",
                 "p-valor"],
                binbiv_rows,
                "Asociación de las variables binarias con la cardiotoxicidad: "
                "tasa de evento por grupo, razón de momios (OR) y contraste de "
                "independencia. Las columnas indican la tasa de cardiotoxicidad "
                "observada en dos subgrupos distintos, pacientes sin el "
                "antecedente y pacientes con él, por lo que no tienen por qué "
                "sumar 100.",
                "tab:biv_bin",
                colspec="lcccc",
                note="OR: razón de momios (odds ratio). Contraste de Fisher o "
                     "$\\chi^2$ según las frecuencias esperadas. "
                     "* significativo a $\\alpha=0.05$. "
                     "$^{\\dagger}$ alguna casilla de la tabla es nula: el OR "
                     "se ha calculado con la corrección de continuidad de "
                     "Haldane--Anscombe (se suma 0.5 a las cuatro casillas), "
                     "de modo que debe leerse como una cota estabilizada y no "
                     "como una estimación del efecto.")

    # ------------------------------------------------------------------
    # 7b. Bivariante: continuas frente al tiempo hasta el evento
    # ------------------------------------------------------------------
    # Las variables binarias frente al tiempo hasta el evento se analizan en
    # R (eda_funcional_superv.R) mediante curvas de Kaplan-Meier por subgrupo
    # y el contraste log-rank. Aqui solo se trata la parte continua: el
    # tiempo registrado solo es el tiempo "hasta el evento" real para quienes
    # lo sufrieron (CTRCD=1); en el resto es tiempo de censura, por lo que la
    # correlacion se restringe a los pacientes con evento.
    log.info("== Continuas vs tiempo hasta el evento (solo CTRCD=1) ==")
    sub_evt = df[df[TARGET] == 1]
    n_evt = len(sub_evt)
    time_rows = []
    for v in CONT:
        x = sub_evt[v]
        y = sub_evt[TIME]
        mask = x.notna() & y.notna()
        rho, p = stats.spearmanr(x[mask], y[mask])
        time_rows.append((texsafe(LAB[v]), f"{rho:.3f}", f"{p:.3f}"))
        log.info("  %-12s rho=%.3f p=%.4f (n=%d)", v, rho, p, mask.sum())
    write_table("tab_biv_time.tex",
                ["Variable", "Correlación de Spearman", "p-valor"],
                time_rows,
                "Correlación de Spearman entre las variables continuas y el "
                "tiempo hasta el evento, calculada entre los pacientes que "
                "desarrollaron cardiotoxicidad ($n={}$).".format(n_evt),
                "tab:biv_time")

    # ------------------------------------------------------------------
    # 8. Multivariante: matriz de correlación
    # ------------------------------------------------------------------
    log.info("== Correlación ==")
    # Matriz combinada: Spearman bajo la diagonal y Pearson por encima, de
    # modo que se muestran ambos coeficientes en una sola figura.
    corr_s = df[CONT].corr(method="spearman")
    corr_p = df[CONT].corr(method="pearson")
    # Se construye sobre un array NumPy escribible (compatible con pandas 2.x
    # y 3.x; en pandas 3 el ``.values`` de un DataFrame es de solo lectura).
    combined = corr_s.to_numpy(copy=True)
    iu = np.triu_indices_from(combined, k=1)   # triángulo superior (Pearson)
    combined[iu] = corr_p.to_numpy()[iu]
    combined = pd.DataFrame(combined, index=corr_s.index,
                            columns=corr_s.columns)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sns.heatmap(combined, annot=True, fmt=".2f", cmap="RdBu_r",
                center=0, vmin=-1, vmax=1, square=True,
                cbar_kws={"label": "Coeficiente de correlación"},
                xticklabels=[LAB[v] for v in CONT],
                yticklabels=[LAB[v] for v in CONT], ax=ax,
                annot_kws={"size": 8})
    ax.set_title("Matriz de correlación de las variables continuas\n"
                 "(Pearson sobre la diagonal, Spearman bajo la diagonal)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    savefig(fig, "eda_corr.pdf")

    # ------------------------------------------------------------------
    # 9. Multivariante: PCA
    # ------------------------------------------------------------------
    log.info("== PCA ==")
    Xc = df[CONT].dropna()
    yc = df.loc[Xc.index, TARGET]
    Z = StandardScaler().fit_transform(Xc)
    pca = PCA(n_components=PCA_COMPONENTS).fit(Z)
    P = pca.transform(Z)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for k in [0, 1]:
        m = (yc == k).values
        ax.scatter(P[m, 0], P[m, 1], s=22, alpha=0.6, color=PAL[k],
                   label=("CTRCD" if k else "No CTRCD"))
    ax.set_xlabel(f"CP1 ({100*pca.explained_variance_ratio_[0]:.1f}% var.)")
    ax.set_ylabel(f"CP2 ({100*pca.explained_variance_ratio_[1]:.1f}% var.)")
    ax.set_title("Análisis de componentes principales (variables continuas)")
    ax.legend()
    savefig(fig, "eda_pca.pdf")
    log.info("  var. explicada CP1+CP2 = %.1f%%",
             100*pca.explained_variance_ratio_[:2].sum())

    # Coeficientes (loadings) de cada variable en las dos primeras componentes,
    # para poder interpretar que combinacion de variables representa cada CP
    # (en vez de solo proyectar a las pacientes en el plano CP1-CP2).
    loadings = pca.components_.T.copy()  # (n_vars, n_components)
    for k in range(loadings.shape[1]):
        j = int(np.argmax(np.abs(loadings[:, k])))
        if loadings[j, k] < 0:
            loadings[:, k] *= -1
    pca_rows = [(texsafe(LAB[v]), f"{loadings[i, 0]:.3f}",
                 f"{loadings[i, 1]:.3f}") for i, v in enumerate(CONT)]
    write_table(
        "tab_pca_loadings.tex",
        ["Variable",
         f"CP1 ({100*pca.explained_variance_ratio_[0]:.1f}%)",
         f"CP2 ({100*pca.explained_variance_ratio_[1]:.1f}%)"],
        pca_rows,
        "Coeficientes (\\textit{loadings}) de las variables continuas en las "
        "dos primeras componentes principales del ACP.",
        "tab:pca_loadings")
    for i, v in enumerate(CONT):
        log.info("  %-25s %7.3f %7.3f", LAB[v], loadings[i, 0], loadings[i, 1])

    log.info("Figuras EDA clasico OK en %s", FIG)
    log.info("Recuerda generar la parte funcional/supervivencia con "
             "eda_funcional_superv.R")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()
