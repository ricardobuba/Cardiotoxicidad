# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# cindex_v2.R - Recalculo del indice de concordancia para los modelos
# seleccionados por eliminacion hacia atras (correcciones A1 y A2).
#
# La tabla tab_mod_cindex.tex corresponde a los modelos ANTERIORES
# (M1 = heart_rate + weight, M2 = FPC1 + FPC2, M3 = heart_rate + FPC1 + FPC2)
# y por tanto ya no describe nada. Este guion la rehace con:
#
#   M1: incidencia heart_rate      | latencia LVEF
#   M2: incidencia FPC4            | latencia FPC5
#   M3: incidencia heart_rate      | latencia FPC6
#
# (los conjuntos se leen de los .rds, no se escriben a mano)
#
# La validacion cruzada de tfg_comun.R no sirve tal cual, porque asume que
# incidencia y latencia llevan las MISMAS covariables: su predecir_curacion()
# lo comprueba con stopifnot(identical(fit$betanm, fit$bnm[-1])). Aqui se
# generaliza a conjuntos distintos:
#
#   S(t | z) = (1 - pi(z_inc)) + pi(z_inc) * S0(t)^exp(beta' x_lat)
#
# Esta seccion es la que la tutora llama "5.6 Comparacion de resultados".
#
# DEPENDE de cache/backward_clinicas.rds y cache/modelos_v2.rds
#
# Uso:  Rscript R/cindex_v2.R
# Salidas: salidas/tablas/tab_mod_cindex_v2.tex
#          salidas/figuras/mod_cindex_v2.png
#          cache/cindex_v2.rds
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))
source(file.path(ROOT, "R", "curacion_lib.R"))
suppressPackageStartupMessages(library(smcure))

# tfg_comun.R solo etiqueta FPC1..FPC3; esto tolera FPC4..FPC8 sin tocarlo.
ETIQ_O <- function(v) ifelse(is.na(ETIQ[v]), v, ETIQ[v])

K <- 5; R <- 20; SEMILLA <- 20260810

f1 <- file.path(ROOT, "cache", "backward_clinicas.rds")
f2 <- file.path(ROOT, "cache", "modelos_v2.rds")
if (!file.exists(f1) || !file.exists(f2))
  stop("Faltan los .rds. Ejecuta antes seleccion_backward.R y modelizacion_funcional.R")
A1 <- readRDS(f1); V2 <- readRDS(f2)

MOD <- list(
  M1 = list(inc = A1$covs_inc,    lat = A1$covs_lat),
  M2 = list(inc = V2$M2$inc,      lat = V2$M2$lat),
  M3 = list(inc = V2$M3$inc,      lat = V2$M3$lat))

cat("== Modelos a comparar ==\n")
for (m in names(MOD))
  cat(sprintf("  %s: incidencia = %-28s latencia = %s\n", m,
              paste(MOD[[m]]$inc, collapse = " + "),
              paste(MOD[[m]]$lat, collapse = " + ")))

cat("\n== Preparacion de datos ==\n")
D <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta
t0 <- median(tt[dd == 1])
cat(sprintf("  horizonte t0 = %.0f dias\n", t0))

# --- Supervivencia poblacional predicha con incidencia y latencia distintas
predecir_curacion2 <- function(fit, newdata, t0, time_train, s_train) {
  s0 <- s0_en(time_train, s_train, t0)
  Zi <- cbind(1, as.matrix(newdata[, fit$covs_inc, drop = FALSE]))
  Xl <- as.matrix(newdata[, fit$covs_lat, drop = FALSE])
  stopifnot(ncol(Zi) == length(fit$b), ncol(Xl) == length(fit$beta))
  pi <- 1 / (1 + exp(-drop(Zi %*% fit$b)))
  (1 - pi) + pi * s0^exp(drop(Xl %*% fit$beta))
}

cindex <- function(score, t, d)
  as.numeric(concordance(Surv(t, d) ~ score)$concordance)

nombres <- c("Cox9", "CoxHR", names(MOD))
res <- setNames(vector("list", length(nombres)), nombres)
fallos <- 0L; descartadas <- 0L
t_ini <- Sys.time()
set.seed(SEMILLA)
for (r in seq_len(R)) {
  fold <- sample(rep(1:K, length.out = length(tt)))
  sc <- setNames(lapply(nombres, function(z) numeric(length(tt))), nombres)
  ok_rep <- TRUE
  for (k in 1:K) {
    tr <- fold != k; te <- !tr
    c9 <- coxph(as.formula(paste("Surv(tt[tr], dd[tr]) ~",
                                 paste(VARS_CLIN, collapse = " + "))),
                data = dat[tr, ])
    sc$Cox9[te]  <- -predict(c9, newdata = dat[te, ], type = "lp")
    c1 <- coxph(Surv(tt[tr], dd[tr]) ~ heart_rate, data = dat[tr, ])
    sc$CoxHR[te] <- -predict(c1, newdata = dat[te, ], type = "lp")
    for (m in names(MOD)) {
      f <- tryCatch(suppressMessages(suppressWarnings(
             ajustar_curacion2(MOD[[m]]$inc, MOD[[m]]$lat, dat[tr, ],
                               tt[tr], dd[tr], Var = FALSE))),
           error = function(e) NULL)
      if (is.null(f) || fallo(f)) { fallos <- fallos + 1L; ok_rep <- FALSE; next }
      sc[[m]][te] <- predecir_curacion2(f, dat[te, , drop = FALSE], t0, tt[tr], f$s)
    }
  }
  if (!ok_rep) { descartadas <- descartadas + 1L; next }
  for (nm in nombres) res[[nm]] <- c(res[[nm]], cindex(sc[[nm]], tt, dd))
  cat(sprintf("  repeticion %2d de %d  (%.1f min acumulados)\n", r, R,
              as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))
}
if (fallos > 0)
  cat("\n  [AVISO]", fallos, "ajustes no convergieron;", descartadas,
      "repeticiones descartadas.\n")
val <- sapply(res, length)
if (any(val < 3)) stop("demasiadas repeticiones descartadas: revisar")

cat("\n== Resultados ==\n")
cind <- data.frame(modelo = names(res), C = sapply(res, mean),
                   sd = sapply(res, sd), n = val, row.names = NULL)
print(cind, digits = 3)

NOM <- c(Cox9 = "Cox (9 clínicas)", CoxHR = "Cox (frecuencia cardíaca)",
         M1 = "Curación M1 (clínicas)", M2 = "Curación M2 (ECG)",
         M3 = "Curación M3 (mixto)")
NOMA <- c(Cox9 = "Cox (9 clinicas)", CoxHR = "Cox (frec. cardiaca)",
          M1 = "Curacion M1 (clinicas)", M2 = "Curacion M2 (ECG)",
          M3 = "Curacion M3 (mixto)")
det <- sapply(names(MOD), function(m)
  sprintf("%s: incidencia %s, latencia %s", m,
          paste(ETIQ_O(MOD[[m]]$inc), collapse = " + "),
          paste(ETIQ_O(MOD[[m]]$lat), collapse = " + ")))

escribir_tabla_tex("tab_mod_cindex_v2.tex", "lcc",
  "\\textbf{Modelo} & \\textbf{C-index} & \\textbf{DT}",
  sprintf("%s & %s & %s \\\\", NOM[cind$modelo], fmt_n(cind$C, 3), fmt_n(cind$sd, 3)),
  sprintf(paste("Índice de concordancia estimado por validación cruzada de %d particiones repetida %d veces",
                "($n=%d$, %d eventos) sobre los modelos seleccionados por eliminación hacia atrás.",
                "Para los modelos de curación la puntuación de riesgo es la supervivencia poblacional",
                "predicha en $t_0=%.0f$ días, la mediana del tiempo hasta el evento. La \\gls{fpca} y el",
                "horizonte $t_0$ se calculan una única vez sobre toda la muestra, por ser pasos no",
                "supervisados. La desviación típica recoge la variabilidad entre repeticiones, no la",
                "incertidumbre de la estimación."),
          K, R, length(tt), sum(dd), t0),
  "tab:mod-cindex",
  nota = paste(det, collapse = ". "))

figura("mod_cindex_v2.png", 9, 5, {
  op <- par(mar = c(4.5, 11, 3, 1.5), mgp = c(2.6, 0.8, 0))
  boxplot(rev(res), names = rev(NOMA[names(res)]), horizontal = TRUE, las = 1,
          col = adjustcolor(BLUE, 0.6), border = GRAY,
          xlab = "C-index (validacion cruzada)",
          main = "Discriminacion de los modelos")
  abline(v = 0.5, col = PINK, lty = 2, lwd = 2)
  par(op)
})

saveRDS(list(res = res, cind = cind, MOD = MOD, t0 = t0, K = K, R = R),
        file.path(ROOT, "cache", "cindex_v2.rds"))
cat(sprintf("\nTiempo total: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))
cat("TODO OK -> cache/cindex_v2.rds\n")
