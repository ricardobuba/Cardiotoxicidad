# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# modelizacion_curacion.R - Capitulo de modelizacion del TFG.
#
# Sustituye a aplicacion_modelos_curacion.R, ampliandolo a todos los modelos
# que se comparan en el Capitulo "Analisis de supervivencia y modelizacion":
#
#   (0) Contraste de Maller-Zhou de seguimiento suficiente.
#   (1) Cox clasico (referencia): con las 9 variables clinicas y con la unica
#       que resulta relevante, mas la comprobacion de riesgos proporcionales.
#   (2) Cribado univariante con modelos de curacion de las 9 clinicas.
#   (3) Modelos de mixtura de curacion (smcure, Sy y Taylor 2000):
#         M1 clinicas       : heart_rate + weight
#         M2 funcional      : FPC1 + FPC2   (aportacion propia: la senal ECG)
#         M3 mixto          : heart_rate + FPC1 + FPC2
#   (4) [el modelo single-index vive en modelizacion_sicure.R, que lo cachea]
#   (5) Evaluacion: C-index validado por validacion cruzada repetida.
#
# Toda la preparacion de datos (atipicos funcionales + FPCA + emparejamiento
# con las clinicas) se reutiliza de tfg_comun.R, que comparte con
# eda_funcional_superv.R, de modo que se parte exactamente del mismo
# subconjunto de curvas y de los mismos scores que en el capitulo de EDA.
#
# Salidas: salidas/tablas/tab_mod_{cox,screening,curacion,cindex}.tex
#          salidas/figuras/mod_{prob_curacion,superv_predicha,cindex}.png
#          RESULTADOS_ECG/modelos_curacion.rds
#
# Requiere: install.packages(c("smcure", "npcure"))
#           (ademas de readr, fda, survival, ya usados en el resto del TFG)
# ---------------------------------------------------------------------------

# Se ejecuta desde la carpeta del proyecto (Rscript R/modelizacion_curacion.R).
# Sin rutas absolutas: el guion tiene que funcionar en cualquier maquina.
if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

# ---------------------------------------------------------------------------
# 0. Datos y contraste de seguimiento suficiente
# ---------------------------------------------------------------------------
cat("== Preparacion de datos ==\n")
D  <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta

cat("\n== Contraste de Maller-Zhou (seguimiento suficiente) ==\n")
mz <- maller_zhou(D$clin$time, D$clin$CTRCD)
cat(sprintf("  cohorte completa (n=%d): N=%d, delta=%.0f dias, intervalo=(%.0f, %.0f], p=%.4g\n",
            nrow(D$clin), mz$statistic, mz$delta, mz$interval[1], mz$interval[2], mz$pvalue))
if (requireNamespace("npcure", quietly = TRUE)) {
  mzn <- npcure::testmz(time, CTRCD, D$clin)
  cat(sprintf("  npcure::testmz -> estadistico=%s, p=%.4g\n",
              format(mzn$aux$statistic), mzn$pvalue))
  if (abs(mzn$pvalue - mz$pvalue) > 1e-8)
    cat("  [AVISO] npcure y la formula de Maller-Zhou (1992) discrepan: revisar.\n")
} else cat("  [aviso] instala 'npcure' para contrastar con npcure::testmz\n")

# ---------------------------------------------------------------------------
# 1. Referencia: modelo de Cox clasico
# ---------------------------------------------------------------------------
cat("\n== (1) Cox clasico ==\n")
cox9  <- coxph(as.formula(paste("Surv(tt, dd) ~", paste(VARS_CLIN, collapse = " + "))),
               data = dat)
coxhr <- coxph(Surv(tt, dd) ~ heart_rate, data = dat)
print(summary(cox9)$coefficients)
zph <- cox.zph(cox9)
cat("\n  cox.zph GLOBAL: chisq =", round(zph$table["GLOBAL", "chisq"], 2),
    " gl =", zph$table["GLOBAL", "df"], " p =", round(zph$table["GLOBAL", "p"], 4), "\n")

ic9 <- summary(cox9)$conf.int      # columnas: exp(coef), exp(-coef), lower .95, upper .95
filas <- sprintf("%s & %s & %s & [%s; %s] & %s \\\\",
                 ETIQ[VARS_CLIN],
                 fmt_n(coef(cox9), 3), fmt_n(exp(coef(cox9)), 3),
                 fmt_n(ic9[, "lower .95"], 3), fmt_n(ic9[, "upper .95"], 3),
                 fmt_p(summary(cox9)$coefficients[, "Pr(>|z|)"]))
escribir_tabla_tex("tab_mod_cox.tex", "lcccc",
  "\\textbf{Variable} & \\textbf{$\\hat\\beta$} & \\textbf{HR} & \\textbf{IC 95\\,\\%} & \\textbf{$p$}",
  filas,
  sprintf("Modelo de Cox con las nueve variables clínicas sobre el subconjunto con imagen ($n=%d$, %d eventos). Variables estandarizadas, de modo que el HR corresponde a un incremento de una desviación típica. Contraste global de riesgos proporcionales: $p=%.3f$.",
          length(tt), sum(dd), zph$table["GLOBAL", "p"]),
  "tab:mod-cox")

# ---------------------------------------------------------------------------
# 2. Cribado univariante con modelos de curacion
# ---------------------------------------------------------------------------
cat("\n== (2) Cribado univariante (modelo de curacion, una covariable) ==\n")
scr <- do.call(rbind, lapply(VARS_CLIN, function(v) {
  f <- ajustar_curacion(v, dat, tt, dd)
  data.frame(var = v, b = f$b[2], b_p = f$b_pvalue[2],
             beta = f$beta[1], beta_p = f$beta_pvalue[1], row.names = NULL)
}))
scr <- scr[order(scr$b_p), ]
print(scr, digits = 3)

filas <- sprintf("%s & %s & %s & %s & %s \\\\", ETIQ[scr$var],
                 fmt_n(scr$b, 3), fmt_p(scr$b_p),
                 fmt_n(scr$beta, 3), fmt_p(scr$beta_p))
escribir_tabla_tex("tab_mod_screening.tex", "lcccc",
  paste("\\textbf{Variable} & \\textbf{$\\hat b$ (incid.)} & \\textbf{$p$}",
        "& \\textbf{$\\hat\\beta$ (laten.)} & \\textbf{$p$}"),
  filas,
  sprintf("Cribado univariante de las variables clínicas mediante modelos de mixtura de curación ajustados de una en una ($n=%d$, %d eventos), ordenadas por el $p$-valor de la componente de incidencia.",
          length(tt), sum(dd)),
  "tab:mod-screening")

# ---------------------------------------------------------------------------
# 3. Modelos de mixtura de curacion
# ---------------------------------------------------------------------------
COVS <- list(M1 = c("heart_rate", "weight"),
             M2 = c("FPC1", "FPC2"),
             M3 = c("heart_rate", "FPC1", "FPC2"))
fits <- list()
for (m in names(COVS)) {
  cat("\n== (3)", m, ":", paste(COVS[[m]], collapse = " + "), "==\n")
  fits[[m]] <- ajustar_curacion(COVS[[m]], dat, tt, dd)
  printsmcure(fits[[m]], Var = TRUE)
  cat(sprintf("  fraccion de curacion media estimada: %.3f\n",
              fraccion_curacion(fits[[m]], COVS[[m]], dat)))
}

# Tabla conjunta de los tres modelos
filas <- c()
NOM <- c(M1 = "M1: clínicas", M2 = "M2: funcional (ECG)", M3 = "M3: mixto")
for (m in names(COVS)) {
  f <- fits[[m]]
  filas <- c(filas, sprintf("\\multicolumn{6}{l}{\\textbf{%s}} \\\\", NOM[m]))
  for (i in seq_along(f$b))
    filas <- c(filas, sprintf("\\quad %s & Incidencia & %s & %s & [%s; %s] & %s \\\\",
      ifelse(is.na(ETIQ[f$bnm[i]]), f$bnm[i], ETIQ[f$bnm[i]]),
      fmt_n(f$b[i], 3), fmt_n(exp(f$b[i]), 3),
      fmt_n(exp(f$b[i] - 1.96 * f$b_sd[i]), 3),
      fmt_n(exp(f$b[i] + 1.96 * f$b_sd[i]), 3),
      fmt_p(f$b_pvalue[i])))
  for (i in seq_along(f$beta))
    filas <- c(filas, sprintf("\\quad %s & Latencia & %s & %s & [%s; %s] & %s \\\\",
      ETIQ[f$betanm[i]], fmt_n(f$beta[i], 3), fmt_n(exp(f$beta[i]), 3),
      fmt_n(exp(f$beta[i] - 1.96 * f$beta_sd[i]), 3),
      fmt_n(exp(f$beta[i] + 1.96 * f$beta_sd[i]), 3),
      fmt_p(f$beta_pvalue[i])))
}
escribir_tabla_tex("tab_mod_curacion.tex", "llcccc",
  "\\textbf{Covariable} & \\textbf{Componente} & \\textbf{Coef.} & \\textbf{OR / HR} & \\textbf{IC 95\\,\\%} & \\textbf{$p$}",
  filas,
  sprintf("Modelos de mixtura de curación ajustados con \\textit{smcure} ($n=%d$, %d eventos). En la incidencia el coeficiente es el de una regresión logística sobre la probabilidad de ser susceptible (se muestra el \\textit{odds ratio}); en la latencia, el de un modelo de Cox (se muestra el \\textit{hazard ratio}). Errores estándar por bootstrap ($B=%d$); el intervalo de confianza al 95\\,\\%% es $\\exp(\\hat\\theta \\pm 1{,}96\\,\\mathrm{EE})$, calculado sobre la escala del coeficiente. Covariables estandarizadas.",
          length(tt), sum(dd), NBOOT),
  "tab:mod-curacion")

# ---------------------------------------------------------------------------
# 4. Modelo single-index con la senal funcional completa
# ---------------------------------------------------------------------------
# El ajuste de sicure.vf NO se hace aqui: vive en modelizacion_sicure.R, que es
# el guion que produce todas las cifras y la figura del apartado 5.4.2 y que
# ademas cachea el resultado (tres ajustes de ~9 min cada uno). Antes se
# ajustaba tambien en este guion, con otra semilla, y el resultado no lo usaba
# ninguna tabla ni figura: eran ~9 minutos por nada y dos fuentes posibles para
# el mismo numero.
#
#   Rscript R/modelizacion_sicure.R

# ---------------------------------------------------------------------------
# 5. Evaluacion: C-index por validacion cruzada repetida
# ---------------------------------------------------------------------------
cat("\n== (5) C-index validado (validacion cruzada repetida) ==\n")
K <- 5; R <- 20
t0 <- median(tt[dd == 1])   # horizonte de prediccion: mediana del tiempo al evento
res  <- validacion_cruzada(dat, tt, dd, COVS, t0, K = K, R = R)
cind <- tabla_cindex(res, length(tt), sum(dd), t0, K = K, R = R)
print(cind, digits = 3)

# ---------------------------------------------------------------------------
# 6. Figuras
# ---------------------------------------------------------------------------
cat("\n== (6) Figuras ==\n")
figuras_modelo_m3(fits$M3, COVS$M3, dat, tt, dd)
figura_cindex(res)

# ---------------------------------------------------------------------------
# 7. Guardado
# ---------------------------------------------------------------------------
saveRDS(list(cox9 = cox9, coxhr = coxhr, zph = zph, screening = scr,
             curacion = fits, cindex = res, mz = mz,
             covs = COVS, t0 = t0),
        file.path(ROOT, "RESULTADOS_ECG", "modelos_curacion.rds"))
cat("\nTODO OK. Resultados en RESULTADOS_ECG/modelos_curacion.rds\n")
