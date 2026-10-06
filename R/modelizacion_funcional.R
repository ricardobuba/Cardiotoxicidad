# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# modelizacion_funcional.R - Correccion A2 de la tutora (10/09/2026):
#
#   "Esto no es cierto, se puede trabajar con mas covariables. Coge, por
#    ejemplo, 8, para llegar a un 90% y ajusta el modelo de la misma forma
#    que en la seccion anterior"
#
# Rehace M2 y M3 con las OCHO primeras componentes principales funcionales
# (el mismo 90 % de varianza que ya retiene el modelo single-index de la
# Seccion 5.4.2, con lo que desaparece la incoherencia de usar 2 componentes
# en un apartado y 8 en el siguiente).
#
# "De la misma forma que en la seccion anterior" se interpreta como el
# PROCEDIMIENTO COMPLETO, backward incluido:
#
#   M2 : arranca con FPC1..FPC8 en incidencia y en latencia, y se aplica la
#        misma eliminacion hacia atras desacoplada que en A1.
#   M3 : arranca con las clinicas que sobrevivieron en A1 MAS FPC1..FPC8, en
#        las dos componentes, y se repite la eliminacion.
#
# Aplicar en 5.4 un criterio de seleccion distinto del de 5.3 volveria a abrir
# justo el frente que la correccion A1 cierra.
#
# DEPENDE de cache/backward_clinicas.rds -> ejecuta antes seleccion_backward.R
#
# NO modifica tfg_comun.R, modelizacion_curacion.R ni seleccion_backward.R.
#
# Uso:   Rscript R/modelizacion_funcional.R
# Salidas:
#   salidas/tablas/tab_mod_backward_m2.tex
#   salidas/tablas/tab_mod_backward_m3.tex
#   salidas/tablas/tab_mod_curacion_v2.tex   (M1, M2 y M3 juntos)
#   cache/modelos_v2.rds
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

B_BUSQUEDA    <- 100
B_FINAL       <- NBOOT
ALPHA_SALIDA  <- 0.05
SEED_BUSQUEDA <- 20260910
source(file.path(ROOT, "R", "curacion_lib.R"))

K_FUN     <- 8
VARS_FPC  <- paste0("FPC", seq_len(K_FUN))

# tfg_comun.R solo etiqueta FPC1..FPC3; se completa aqui sin tocarlo.
ETIQ[VARS_FPC] <- VARS_FPC

BW <- file.path(ROOT, "cache", "backward_clinicas.rds")
if (!file.exists(BW))
  stop("Falta ", BW, ".\n  Ejecuta antes:  Rscript R/seleccion_backward.R")
A1 <- readRDS(BW)

cat("== Preparacion de datos ==\n")
D <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta
stopifnot(all(VARS_FPC %in% names(dat)))

cat("\n== Varianza funcional retenida ==\n")
vp <- D$pca$varprop
cat(sprintf("  componentes 1..%d: %s\n", K_FUN,
            paste(sprintf("%.1f%%", 100 * vp[seq_len(K_FUN)]), collapse = " ")))
cat(sprintf("  acumulada con %d componentes: %.1f%%\n",
            K_FUN, 100 * sum(vp[seq_len(K_FUN)])))
if (sum(vp[seq_len(K_FUN)]) < 0.90)
  cat("  [!] no se alcanza el 90 %: ajusta K_FUN o revisa el suavizado\n")

cat("\n== Modelo A1 (clinicas) seleccionado en la seccion anterior ==\n")
cat("  incidencia:", paste(A1$covs_inc, collapse = " + "), "\n")
cat("  latencia  :", paste(A1$covs_lat, collapse = " + "), "\n")

# ---------------------------------------------------------------------------
# M2: solo la senal funcional
# ---------------------------------------------------------------------------
cat("\n\n############ M2: backward sobre FPC1..FPC", K_FUN, " ############\n", sep = "")
t0 <- Sys.time()
bw2 <- backward_curacion(dat, tt, dd, vars = VARS_FPC,
                         alpha = ALPHA_SALIDA, nboot = B_BUSQUEDA,
                         seed = SEED_BUSQUEDA)
cat(sprintf("\n[M2] busqueda: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
M2 <- ajustar_curacion2(bw2$covs_inc, bw2$covs_lat, dat, tt, dd,
                        nboot = B_FINAL, seed = SEED, Var = TRUE)
stopifnot(!fallo(M2))

# ---------------------------------------------------------------------------
# M3: clinicas seleccionadas + senal funcional
# ---------------------------------------------------------------------------
cat("\n\n############ M3: backward sobre clinicas seleccionadas + FPC1..FPC", K_FUN, " ############\n", sep = "")
ini_inc <- union(A1$covs_inc, VARS_FPC)
ini_lat <- union(A1$covs_lat, VARS_FPC)
cat("  arranque incidencia:", paste(ini_inc, collapse = " + "), "\n")
cat("  arranque latencia  :", paste(ini_lat, collapse = " + "), "\n")
# Cada componente arranca con las clinicas que sobrevivieron EN ESA componente
# en A1, mas las ocho funcionales. Una clinica descartada de la incidencia no
# reaparece en ella.
t0 <- Sys.time()
bw3 <- backward_curacion(dat, tt, dd,
                         vars_inc = ini_inc, vars_lat = ini_lat,
                         alpha = ALPHA_SALIDA, nboot = B_BUSQUEDA,
                         seed = SEED_BUSQUEDA)
cat(sprintf("\n[M3] busqueda: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
M3 <- ajustar_curacion2(bw3$covs_inc, bw3$covs_lat, dat, tt, dd,
                        nboot = B_FINAL, seed = SEED, Var = TRUE)
stopifnot(!fallo(M3))

# ---------------------------------------------------------------------------
# Resumen
# ---------------------------------------------------------------------------
frac <- function(fit, covs) 1 - mean(1 / (1 + exp(-drop(
  cbind(1, as.matrix(dat[, covs, drop = FALSE])) %*% fit$b))))

MOD <- list(M1 = list(fit = A1$fit, inc = A1$covs_inc, lat = A1$covs_lat),
            M2 = list(fit = M2,     inc = bw2$covs_inc, lat = bw2$covs_lat),
            M3 = list(fit = M3,     inc = bw3$covs_inc, lat = bw3$covs_lat))
NOM <- c(M1 = "M1: clínicas", M2 = "M2: funcional (ECG)", M3 = "M3: mixto")

cat("\n\n================ RESUMEN ================\n")
for (m in names(MOD)) {
  x <- MOD[[m]]
  cat(sprintf("\n%s\n  incidencia: %s\n  latencia  : %s\n  fraccion de curacion: %.3f\n",
              NOM[m], paste(x$inc, collapse = " + "),
              paste(x$lat, collapse = " + "), frac(x$fit, x$inc)))
  print(candidatas(x$fit, x$inc, x$lat), digits = 3, row.names = FALSE)
}

# ------------------------------ Tablas -------------------------------------
et <- function(v) ifelse(is.na(ETIQ[v]), v, ETIQ[v])

tabla_traza <- function(bw, fichero, label, cual) {
  filas <- c()
  for (h in bw$historia) {
    if (is.null(h$eliminada)) next
    filas <- c(filas, sprintf("%d & %d & %d & %s & %s & %s \\\\", h$paso,
      length(h$covs_inc) + ifelse(h$eliminada$comp == "incidencia", 1L, 0L),
      length(h$covs_lat) + ifelse(h$eliminada$comp == "latencia", 1L, 0L),
      et(h$eliminada$var),
      ifelse(h$eliminada$comp == "incidencia", "Incidencia", "Latencia"),
      if (isTRUE(h$converge)) fmt_p(h$eliminada$p) else "---"))
  }
  if (!length(filas)) { cat("  [tab]", fichero, "-> sin eliminaciones\n"); return(invisible()) }
  escribir_tabla_tex(fichero, "cccllc",
    paste("\\textbf{Paso} & \\textbf{$k_{\\text{inc}}$} & \\textbf{$k_{\\text{lat}}$}",
          "& \\textbf{Variable eliminada} & \\textbf{Componente} & \\textbf{$p$}"),
    filas,
    sprintf("Eliminación hacia atrás para el modelo %s ($n=%d$, %d eventos; $B=%d$ en la búsqueda, $B=%d$ en el ajuste final).",
            cual, length(tt), sum(dd), B_BUSQUEDA, B_FINAL),
    label)
}
tabla_traza(bw2, "tab_mod_backward_m2.tex", "tab:mod-backward-m2", "M2")
tabla_traza(bw3, "tab_mod_backward_m3.tex", "tab:mod-backward-m3", "M3")

filas <- c()
for (m in names(MOD)) {
  f <- MOD[[m]]$fit
  filas <- c(filas, sprintf("\\multicolumn{6}{l}{\\textbf{%s}} \\\\", NOM[m]))
  for (i in seq_along(f$b))
    filas <- c(filas, sprintf("\\quad %s & Incidencia & %s & %s & [%s; %s] & %s \\\\",
      et(f$bnm[i]), fmt_n(f$b[i], 3), fmt_n(exp(f$b[i]), 3),
      fmt_n(exp(f$b[i] - 1.96 * f$b_sd[i]), 3),
      fmt_n(exp(f$b[i] + 1.96 * f$b_sd[i]), 3), fmt_p(f$b_pvalue[i])))
  for (i in seq_along(f$beta))
    filas <- c(filas, sprintf("\\quad %s & Latencia & %s & %s & [%s; %s] & %s \\\\",
      et(f$betanm[i]), fmt_n(f$beta[i], 3), fmt_n(exp(f$beta[i]), 3),
      fmt_n(exp(f$beta[i] - 1.96 * f$beta_sd[i]), 3),
      fmt_n(exp(f$beta[i] + 1.96 * f$beta_sd[i]), 3), fmt_p(f$beta_pvalue[i])))
}
escribir_tabla_tex("tab_mod_curacion_v2.tex", "llcccc",
  "\\textbf{Covariable} & \\textbf{Componente} & \\textbf{Coef.} & \\textbf{OR / HR} & \\textbf{IC 95\\,\\%} & \\textbf{$p$}",
  filas,
  sprintf("Modelos de mixtura de curación ajustados con \\textit{smcure} ($n=%d$, %d eventos), con las covariables seleccionadas por eliminación hacia atrás sobre las nueve variables clínicas y las %d primeras componentes principales funcionales (%.1f\\,\\%% de la variabilidad de la señal). Errores estándar por \\textit{bootstrap} ($B=%d$). Covariables estandarizadas.",
          length(tt), sum(dd), K_FUN, 100 * sum(vp[seq_len(K_FUN)]), B_FINAL),
  "tab:mod-curacion")

saveRDS(list(M1 = MOD$M1, M2 = MOD$M2, M3 = MOD$M3,
             bw2 = bw2, bw3 = bw3, K_FUN = K_FUN,
             varprop = vp, fracciones = sapply(names(MOD), function(m)
               frac(MOD[[m]]$fit, MOD[[m]]$inc))),
        file.path(ROOT, "cache", "modelos_v2.rds"))

cat("\n== Figuras ==\n")
cat("  NO se regeneran aqui a proposito: figuras_modelo_m3() dibuja un panel\n")
cat("  por covariable y con", length(MOD$M3$inc), "covariables en M3 la figura\n")
cat("  mod_prob_curacion habria que rediseñarla. Lo vemos cuando esten los\n")
cat("  conjuntos finales.\n")
cat("\nTODO OK -> cache/modelos_v2.rds\n")
