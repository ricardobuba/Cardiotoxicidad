# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Parte EN R del analisis exploratorio de la memoria (enfoque hibrido):
#   * Kaplan-Meier (global y por subgrupos)  -> paquetes survival / survminer
#   * Datos funcionales del ECG: curvas y FPCA -> paquete fda (pca.fd).
#
# El FPCA reutiliza el enfoque del script proporcionado por el Departamento de
# Estadistica de la UDC ("Codigo_aplicacion_smcure_sicure.R"): se construye el
# objeto funcional con fd() y se aplica pca.fd(..., nharm = 8). Ese script
# ORIGINAL no se modifica; aqui se adapta su metodologia a la senal ECG
# extraida de las imagenes (en lugar del ciclo TDI del dataset original).
# El resto de figuras y tablas del EDA las genera eda_analysis.py.
#
# La carga de datos, la deteccion de atipicos funcionales y la FPCA NO se
# programan aqui: viven en tfg_comun.R, que comparten este script y
# modelizacion_curacion.R. Asi el capitulo de analisis descriptivo y el de
# modelizacion parten exactamente del mismo subconjunto de curvas y de los
# mismos scores, y no pueden desincronizarse.
#
# Genera (en salidas/figuras/, que luego se sube a imaxes/eda):
#   eda_km_global.pdf, eda_km_groups.pdf,
#   eda_ecg_sample.pdf, eda_ecg_mean.pdf,
#   eda_fpca_var.pdf, eda_fpca_components.pdf, eda_fpca_scores.pdf
# y la tabla LaTeX salidas/tablas/tab_fpca.tex (-> contido/tablas_eda)
# ---------------------------------------------------------------------------

# Se ejecuta desde la carpeta del proyecto (Rscript R/eda_funcional_superv.R).
# Sin rutas absolutas: el guion tiene que funcionar en cualquier maquina.
if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(survminer))

# Rutas (FIG, TAB), colores (PINK, BLUE, PAL), NHARM, FACTOR_FB y los helpers
# openpdf()/done()/escribir_tabla_tex() los define tfg_comun.R.

cat("== Preparacion de datos (tfg_comun.R) ==\n")
D <- preparar_datos()

# ----------------------------------------------------------------------
# 1. Kaplan-Meier (datos clinicos completos, 531 pacientes)
# ----------------------------------------------------------------------
cat("== Kaplan-Meier ==\n")
clin <- D$clin
cat("  clinicas:", nrow(clin), "pacientes |", sum(clin$CTRCD), "eventos\n")

# (a) Curva global
fit_g <- survfit(Surv(time, CTRCD) ~ 1, data = clin)
g_glob <- ggsurvplot(fit_g, data = clin, conf.int = TRUE,
                     palette = PINK, ylim = c(0, 1),
                     xlab = "Tiempo (dias)",
                     ylab = "Prob. de no sufrir cardiotoxicidad",
                     title = "Curva de Kaplan-Meier (cohorte global)",
                     legend = "none", ggtheme = theme_bw())
openpdf("eda_km_global.pdf", 7.5, 5)
print(g_glob$plot)
done("eda_km_global.pdf")

# (a bis) Test de Maller-Zhou de follow-up suficiente. El plateau de la curva
# KM global sugiere una fraccion de pacientes que no desarrollara el evento
# (curadas); este contraste lo confirma formalmente y justifica los modelos de
# curacion. Se usa npcure::testmz (Maller y Zhou, 1992), del que la directora
# es coautora. maller_zhou() de tfg_comun.R reproduce la misma formula y sirve
# de comprobacion cruzada.
mz_f <- maller_zhou(clin$time, clin$CTRCD)
cat(sprintf("  Maller-Zhou (formula): N = %d | p = %.4g | intervalo = (%.0f, %.0f]\n",
            mz_f$statistic, mz_f$pvalue, mz_f$interval[1], mz_f$interval[2]))
if (requireNamespace("npcure", quietly = TRUE)) {
  mz <- npcure::testmz(time, CTRCD, clin)
  cat(sprintf("  Maller-Zhou (npcure): estadistico = %s | p-valor = %.4g\n",
              format(mz$aux$statistic), mz$pvalue))
  if (abs(mz$pvalue - mz_f$pvalue) > 1e-8)
    cat("  [AVISO] npcure y la formula discrepan: revisar maller_zhou().\n")
} else {
  cat("  [aviso] Instala 'npcure' para ejecutar el test de Maller-Zhou:\n",
      "          install.packages(\"npcure\")\n")
}

# (b) Curvas por subgrupos seleccionados, con p-valor log-rank
LAB <- c(AC = "Antraciclinas", antiHER2 = "Terapia anti-HER2",
         HTA = "Hipertension", heart_rhythm = "Fibrilacion auricular")
sel <- names(LAB)
plots <- list()
for (v in sel) {
  sub <- clin[!is.na(clin[[v]]), c("time", "CTRCD", v)]
  sub$grp <- factor(sub[[v]])
  fit <- survfit(Surv(time, CTRCD) ~ grp, data = sub)
  lr  <- survdiff(Surv(time, CTRCD) ~ grp, data = sub)
  pval <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  cat(sprintf("  log-rank %-12s p=%.4f\n", v, pval))
  gg <- ggsurvplot(fit, data = sub, conf.int = FALSE,
                   palette = c(BLUE, PINK), ylim = c(0.55, 1.0),
                   legend.labs = c(paste0(LAB[v], " = 0"),
                                   paste0(LAB[v], " = 1")),
                   legend.title = "", xlab = "Tiempo (dias)",
                   ylab = "Prob. supervivencia",
                   title = sprintf("%s  (log-rank p=%.3f)", LAB[v], pval),
                   ggtheme = theme_bw())
  plots[[v]] <- gg
}
g_groups <- arrange_ggsurvplots(plots, print = FALSE, ncol = 2, nrow = 2)
# OJO: print(g_groups) sobre un device pdf() deja una primera pagina EN BLANCO
# (y \includegraphics mostraria esa pagina vacia). ggsave escribe una sola pagina.
ggplot2::ggsave(file.path(FIG, "eda_km_groups.pdf"), g_groups, width = 11, height = 8)
cat("  [fig] eda_km_groups.pdf\n")

# ----------------------------------------------------------------------
# 2. Datos funcionales del ECG: curvas crudas y atipicos
# ----------------------------------------------------------------------
# Las curvas, la matriz E_all y la deteccion de atipicos (functional boxplot
# sobre la Modified Band Depth, Sun & Genton 2011, con el mismo factor 1.5 que
# la regla 1.5*IQR de la Seccion 4.2) las calcula preparar_datos().
cat("== Datos funcionales (ECG extraido) ==\n")
E_all  <- D$E_all
yE_all <- D$ecg$CTRCD
xg <- seq(0, 1, length.out = ncol(E_all))  # fase del ciclo cardiaco (0-1)
cat("  curvas ECG (total):", nrow(E_all), " (No CTRCD =", sum(yE_all == 0),
    ", CTRCD =", sum(yE_all == 1), ")\n")

is_outlier <- D$atipicos
cat("  atipicos funcionales detectados (fbplot, factor", FACTOR_FB, "):",
    sum(is_outlier), "de", nrow(E_all), "\n")
cat("  ficheros:", paste(D$ecg$filename[is_outlier], collapse = ", "), "\n")

E  <- D$E
yE <- yE_all[!is_outlier]
idx0 <- which(yE == 0); idx1 <- which(yE == 1)
cat("  curvas ECG (tras excluir atipicos):", nrow(E), " (No CTRCD =",
    length(idx0), ", CTRCD =", length(idx1), ")\n")

# (a) Todas las curvas por grupo + media del grupo (correccion de la tutora:
# se representan TODAS las curvas de cada grupo, no una muestra aleatoria).
openpdf("eda_ecg_sample.pdf", 12, 4.5)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
# Eje Y mas estrecho: percentiles 1-99 (tras excluir atipicos), en lugar
# del rango completo, para apreciar mejor la variacion de cada curva.
ylim <- quantile(E, probs = c(0.01, 0.99))
for (k in 0:1) {
  idx <- if (k == 0) idx0 else idx1
  col <- PAL[as.character(k)]
  plot(NA, xlim = c(0, 1), ylim = ylim, xlab = "Fase del ciclo cardiaco",
       ylab = if (k == 0) "Amplitud de la senal ECG" else "",
       main = sprintf("%s  (n=%d)", if (k) "CTRCD" else "No CTRCD",
                      length(idx)))
  # transparencia y grosor menores porque ahora se dibujan todas las curvas
  for (i in idx) lines(xg, E[i, ], col = adjustcolor(col, 0.12), lwd = 0.5)
  lines(xg, colMeans(E[idx, , drop = FALSE]), col = "black", lwd = 2)
  legend("topright", "Media del grupo", lwd = 2, col = "black", bty = "n")
}
par(op); done("eda_ecg_sample.pdf")

# (b) Curva media +/- 1 DE por grupo
openpdf("eda_ecg_mean.pdf", 9, 5)
op <- par(mar = c(4, 4, 3, 1))
mu0 <- colMeans(E[idx0, ]); sd0 <- apply(E[idx0, ], 2, sd)
mu1 <- colMeans(E[idx1, ]); sd1 <- apply(E[idx1, ], 2, sd)
plot(NA, xlim = c(0, 1), ylim = range(c(mu0 - sd0, mu0 + sd0, mu1 - sd1, mu1 + sd1)),
     xlab = "Fase del ciclo cardiaco", ylab = "Amplitud de la senal ECG",
     main = "Curva ECG media por grupo (+/- 1 desviacion tipica)")
polygon(c(xg, rev(xg)), c(mu0 - sd0, rev(mu0 + sd0)),
        col = adjustcolor(BLUE, 0.18), border = NA)
polygon(c(xg, rev(xg)), c(mu1 - sd1, rev(mu1 + sd1)),
        col = adjustcolor(PINK, 0.18), border = NA)
lines(xg, mu0, col = BLUE, lwd = 2); lines(xg, mu1, col = PINK, lwd = 2)
legend("topright", c("No CTRCD", "CTRCD"), lwd = 2, col = c(BLUE, PINK),
       bty = "n")
par(op); done("eda_ecg_mean.pdf")

# ----------------------------------------------------------------------
# 3. FPCA de la senal ECG (calculada en tfg_comun.R con fda::pca.fd)
# ----------------------------------------------------------------------
cat("== FPCA del ECG (fda::pca.fd) ==\n")
fun_pca <- D$pca
vr <- fun_pca$varprop                    # proporcion de varianza por componente
scoresF <- fun_pca$scores
cat("  var. explicada por componente:", round(100 * vr, 1), "\n")
cat("  var. acumulada:", round(100 * cumsum(vr), 1), "\n")

rng <- fun_pca$harmonics$basis$rangeval
gg  <- seq(rng[1], rng[2], length.out = ncol(E))
xgp <- seq(0, 1, length.out = ncol(E))
meanv <- as.vector(eval.fd(gg, fun_pca$meanfd))
harmv <- eval.fd(gg, fun_pca$harmonics)  # grid x NHARM

# (a) Modos de variacion de las 3 primeras componentes
openpdf("eda_fpca_components.pdf", 13, 4)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (j in 1:3) {
  delta <- 2 * sqrt(fun_pca$values[j])
  yl <- range(c(meanv, meanv + delta * harmv[, j], meanv - delta * harmv[, j]))
  plot(xgp, meanv, type = "l", col = "black", lwd = 2, ylim = yl,
       xlab = "Fase del ciclo cardiaco",
       ylab = if (j == 1) "Amplitud ECG" else "",
       main = sprintf("Componente %d  (%.1f%% var.)", j, 100 * vr[j]))
  lines(xgp, meanv + delta * harmv[, j], col = PINK, lwd = 1.4, lty = 2)
  lines(xgp, meanv - delta * harmv[, j], col = BLUE, lwd = 1.4, lty = 3)
  if (j == 1)
    legend("topright", c("Media", "Media + CP", "Media - CP"),
           col = c("black", PINK, BLUE), lty = c(1, 2, 3),
           lwd = c(2, 1.4, 1.4), bty = "n", cex = 0.9)
}
par(op); done("eda_fpca_components.pdf")

# (b) Varianza explicada (scree)
openpdf("eda_fpca_var.pdf", 7, 4.5)
op <- par(mar = c(4, 4, 3, 1))
comps <- 1:NHARM
bp <- barplot(100 * vr, names.arg = comps, col = adjustcolor(BLUE, 0.8),
              border = NA, ylim = c(0, 100),
              xlab = "Componente principal funcional",
              ylab = "Varianza explicada (%)",
              main = "Varianza explicada por la FPCA del ECG")
lines(bp, 100 * cumsum(vr), type = "o", col = PINK, pch = 19, lwd = 2)
abline(h = 90, col = "gray", lty = 2)
legend("right", c("Individual", "Acumulada"),
       col = c(BLUE, PINK), pch = c(15, 19), bty = "n")
par(op); done("eda_fpca_var.pdf")

# (c) Scores FPC1 vs FPC2 por grupo
openpdf("eda_fpca_scores.pdf", 7, 5.5)
op <- par(mar = c(4, 4, 3, 1))
plot(scoresF[, 1], scoresF[, 2], type = "n",
     xlab = sprintf("Score FPC1 (%.1f%% var.)", 100 * vr[1]),
     ylab = sprintf("Score FPC2 (%.1f%% var.)", 100 * vr[2]),
     main = "Scores de las dos primeras componentes funcionales del ECG")
points(scoresF[idx0, 1], scoresF[idx0, 2], col = adjustcolor(BLUE, 0.6),
       pch = 19)
points(scoresF[idx1, 1], scoresF[idx1, 2], col = adjustcolor(PINK, 0.6),
       pch = 19)
legend("topright", c("No CTRCD", "CTRCD"), col = c(BLUE, PINK), pch = 19,
       bty = "n")
par(op); done("eda_fpca_scores.pdf")

# ----------------------------------------------------------------------
# 4. Tabla LaTeX de varianza explicada
# ----------------------------------------------------------------------
cumv <- cumsum(vr)
escribir_tabla_tex("tab_fpca.tex", "lcc",
  "\\textbf{Componente} & \\textbf{Var. explicada (\\%)} & \\textbf{Var. acumulada (\\%)}",
  sprintf("FPC%d & %.1f & %.1f \\\\", 1:NHARM, 100 * vr, 100 * cumv),
  paste0("Varianza explicada por las ocho primeras componentes ",
         "principales funcionales (FPCA) de la se\u00f1al ECG."),
  "tab:fpca")

cat("\nTODO OK. Figuras en", FIG, "\n")
