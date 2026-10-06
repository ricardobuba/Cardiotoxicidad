# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# figuras_modelizacion.R -- Redibuja SOLO las figuras del capitulo de
# modelizacion, en PNG a 300 ppp. No reajusta ningun modelo: los lee de
# RESULTADOS_ECG/modelos_curacion.rds, que genera modelizacion_curacion.R.
# Tarda segundos (lo unico que recalcula es la FPCA, para tener las
# covariables de cada paciente).
#
# El codigo de dibujo esta en tfg_comun.R (figuras_modelo_m3, figura_cindex),
# que es el mismo que usa el guion maestro: antes estaba duplicado aqui y las
# dos copias podian divergir.
#
# La figura del modelo single-index (mod_sicure_beta.png) NO se rehace aqui:
# la genera modelizacion_sicure.R, que es quien tiene el ajuste cacheado.
#
# Uso:  Rscript R/figuras_modelizacion.R
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))

RDS <- file.path(ROOT, "RESULTADOS_ECG", "modelos_curacion.rds")
if (!file.exists(RDS))
  stop("Falta ", RDS, ". Ejecuta antes modelizacion_curacion.R")
M <- readRDS(RDS)

cat("== Preparacion de datos ==\n")
D   <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta

cat("== Figuras (PNG, 300 ppp) ==\n")
figuras_modelo_m3(M$curacion$M3, M$covs$M3, dat, tt, dd)
figura_cindex(M$cindex)

cat("\nListo. Sube a Overleaf (imaxes/modelizacion) los 3 PNG de\n  ",
    file.path(SAL, "figuras"), "\n")
