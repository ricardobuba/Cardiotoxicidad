# ---------------------------------------------------------------------------
#  diag_cribado.R                                              31/08/2026
#
#  Cierra el unico cabo suelto que queda: por que el p-valor del cribado
#  univariante de la frecuencia cardiaca se movio entre versiones intermedias
#  del codigo (0,018 -> 0,0108 -> 0,0544).
#
#  La idea del contraste: el coeficiente b lo determina el algoritmo EM y NO
#  depende del bootstrap. El error estandar, en cambio, es una estimacion de
#  Monte Carlo sobre B = 500 remuestreos, y cambia con la semilla. Ajustando el
#  mismo modelo con varias semillas se separa lo uno de lo otro:
#
#    - si b sale identico en las 8 semillas (debe salirlo) y el p-valor oscila
#      en un rango que contiene 0,018 y 0,0108, entonces aquel movimiento era
#      ruido de Monte Carlo, no un cambio de resultado;
#    - lo que ya NO puede ser ruido es el salto a 0,0544, porque ahi tambien se
#      movio el coeficiente (0,730 -> 0,782), y eso solo lo explica el cambio
#      de cohorte que trajo el suavizado por GCV (183/26 -> 181/25).
#
#  No modifica ningun resultado ni sobrescribe nada. Duracion: 10-15 minutos.
#
#  Ejecutar desde la carpeta del proyecto:
#     Rscript R/diagnostico/diag_cribado.R > logs/diag_cribado.txt 2>&1
# ---------------------------------------------------------------------------

# ROOT se deduce de la ubicacion del propio fichero, de modo que da igual desde
# que carpeta se lance. Para forzarlo a mano, definir ROOT antes de ejecutar.
if (!exists("ROOT")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..", "..")) else getwd()
}
setwd(ROOT); cat("ROOT =", ROOT, "\n")
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

D   <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta
cat(sprintf("\ncohorte actual: n = %d, eventos = %d\n", length(tt), sum(dd)))
cat(sprintf("B = %d remuestreos por ajuste, cota de descarte C = %d\n\n",
            NBOOT, MAX_COEF_BOOT))

SEMILLAS <- 1:8
res <- do.call(rbind, lapply(SEMILLAS, function(s) {
  cat(sprintf("  semilla %d ...\n", s)); flush.console()
  f <- ajustar_curacion("heart_rate", dat, tt, dd, seed = s)
  data.frame(semilla = s, b = f$b[2], EE = f$b_sd[2], p = f$b_pvalue[2],
             descartadas = f$nboot_descartadas, row.names = NULL)
}))

cat("\n--- cribado univariante: heart_rate, incidencia ---\n")
print(res, digits = 4, row.names = FALSE)

cat(sprintf("\n  b   : min %.4f  max %.4f  (rango %.2g)\n",
            min(res$b), max(res$b), diff(range(res$b))))
cat("        ^ si el rango es ~0, el coeficiente NO depende del bootstrap.\n")
cat(sprintf("  EE  : min %.4f  max %.4f  sd %.4f\n",
            min(res$EE), max(res$EE), sd(res$EE)))
cat(sprintf("  p   : min %.4f  max %.4f  sd %.4f\n",
            min(res$p), max(res$p), sd(res$p)))
cat(sprintf("  replicas descartadas por divergencia: %d en total sobre %d\n",
            sum(res$descartadas), length(SEMILLAS) * NBOOT))

amp <- diff(range(res$p))
cat("\n--- lectura ---\n")
cat(sprintf("  La oscilacion del p-valor entre semillas es de %.4f.\n", amp))
if (amp >= abs(0.018 - 0.0108)) {
  cat("  Es MAYOR que la diferencia 0,018 vs 0,0108, de modo que aquel\n")
  cat("  movimiento entre versiones queda explicado por el error de Monte\n")
  cat("  Carlo del bootstrap, sin cambio de resultado.\n")
} else {
  cat("  Es MENOR que la diferencia 0,018 vs 0,0108: ese movimiento NO se\n")
  cat("  explica solo por la semilla, y habria que buscar otra causa.\n")
}
cat("  El salto a 0,0544 es aparte: ahi cambio tambien el coeficiente\n")
cat("  (0,730 -> 0,782), lo que solo puede venir del cambio de cohorte.\n\n")
