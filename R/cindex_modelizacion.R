# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# cindex_modelizacion.R -- Rehace SOLO la validacion cruzada (seccion 5 de
# modelizacion_curacion.R), la tabla tab_mod_cindex.tex y la figura
# mod_cindex.png. No repite los ajustes finales ni el bootstrap: esos no se
# ven afectados.
#
# Motivos por los que ha habido que rehacerla:
#   1. (jul-2026) ajustar_curacion() fijaba la semilla en CADA llamada, tambien
#      sin bootstrap. Dentro del bucle eso reiniciaba el generador y 19 de las
#      20 repeticiones acababan usando la MISMA particion, de modo que la
#      desviacion tipica salia practicamente nula. Corregido en tfg_comun.R:
#      la semilla solo se fija si Var = TRUE.
#   2. (ago-2026) predecir_curacion() evaluaba la supervivencia basal con
#      approx(..., ties = "ordered") sobre un vector de tiempos SIN ordenar, y
#      ese argumento le dice a approx que no ordene. La S0(t0) que se usaba era
#      arbitraria, sin error ni aviso. Corregido en tfg_comun.R con s0_en().
#      >> Esta correccion CAMBIA el C-index de M1, M2 y M3. <<
#
# El grueso del trabajo esta en validacion_cruzada(), en tfg_comun.R, que es la
# misma funcion que llama el guion maestro: asi no puede volver a pasar que las
# dos copias diverjan.
#
# Uso:  Rscript R/cindex_modelizacion.R
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

RDS <- file.path(ROOT, "RESULTADOS_ECG", "modelos_curacion.rds")
if (!file.exists(RDS))
  stop("Falta ", RDS, ". Ejecuta antes modelizacion_curacion.R")
M <- readRDS(RDS)
COVS <- M$covs

cat("== Preparacion de datos ==\n")
D <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta

K <- 5; R <- 20
t0 <- median(tt[dd == 1])

cat("== C-index por validacion cruzada", K, "x", R, "==\n")
res  <- validacion_cruzada(dat, tt, dd, COVS, t0, K = K, R = R)
cind <- tabla_cindex(res, length(tt), sum(dd), t0, K = K, R = R)
cat("\n== RESULTADO ==\n"); print(cind, digits = 3)

figura_cindex(res)

M$cindex <- res
M$t0     <- t0
saveRDS(M, RDS)
cat("\nListo. Copiame el bloque RESULTADO y vuelve a subir a Overleaf\n",
    " tab_mod_cindex.tex y mod_cindex.png\n")
