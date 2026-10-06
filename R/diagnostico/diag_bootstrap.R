# ---------------------------------------------------------------------------
# diag_bootstrap.R -- elegir el criterio correcto para el bootstrap de smcure
#
# El bucle de smcure solo se queda con las replicas que CONVERGEN
#     bootfit <- em(...); if (bootfit$tau < eps) i <- i+1
# y descarta el resto. Mi primera version las aceptaba todas, y las replicas en
# las que la parte logistica de incidencia se va a separacion (coeficientes
# enormes pero finitos) inflaban el error estandar hasta valores absurdos
# (EE = 23,7 para un coeficiente de 0,99).
#
# Este guion genera las 500 replicas de M1, M2 y M3, guarda TODOS los
# coeficientes y compara el error estandar bajo varios criterios de descarte,
# incluyendo -- para M1 -- el bootstrap propio de smcure como patron oro.
#
# Uso:  Rscript R/diagnostico/diag_bootstrap.R      (unos 15 minutos)
# Salidas: logs/diag_bootstrap.txt (si se redirige) y boot_M1.csv/M2/M3
# ---------------------------------------------------------------------------
if (!exists("ROOT")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..", "..")) else getwd()
}
setwd(ROOT); cat("ROOT =", ROOT, "\n")
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

D <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta
COVS <- list(M1 = c("heart_rate", "weight"),
             M2 = c("FPC1", "FPC2"),
             M3 = c("heart_rate", "FPC1", "FPC2"))
B <- 500

ajustar1 <- function(covs, datos_) {
  f_lat  <- as.formula(paste("Surv(time, delta) ~", paste(covs, collapse = " + ")))
  f_cure <- as.formula(paste("~", paste(covs, collapse = " + ")))
  r <- NULL
  invisible(capture.output(
    r <- smcure::smcure(f_lat, cureform = f_cure, data = datos_,
                        model = "ph", Var = FALSE)))
  r
}

replicas <- function(covs, B = 500, seed = SEED) {
  d  <- cbind(time = tt, delta = dd, dat[, covs, drop = FALSE])
  set.seed(seed)
  i1 <- which(dd == 1); i0 <- which(dd == 0)
  Mb <- NULL; Mg <- NULL; fallos <- 0L
  for (k in seq_len(B)) {
    idx <- c(sample(i1, length(i1), replace = TRUE),
             sample(i0, length(i0), replace = TRUE))
    r <- tryCatch(suppressWarnings(ajustar1(covs, d[idx, , drop = FALSE])),
                  error = function(e) NULL)
    if (is.null(r)) { fallos <- fallos + 1L; next }
    Mb <- rbind(Mb, r$b); Mg <- rbind(Mg, r$beta)
  }
  list(b = Mb, beta = Mg, fallos = fallos)
}

comparar <- function(M, est, titulo) {
  cat("\n  --", titulo, "--\n")
  cat(sprintf("     estimacion: %s\n", paste(sprintf("%+8.4f", est), collapse = " ")))
  cat(sprintf("     max|coef| por replica: p50 = %.2f  p90 = %.2f  p99 = %.2f  max = %.2f\n",
      quantile(apply(abs(M), 1, max), .50), quantile(apply(abs(M), 1, max), .90),
      quantile(apply(abs(M), 1, max), .99), max(abs(M))))
  linea <- function(nm, S) {
    if (nrow(S) < 5) { cat(sprintf("     %-16s n=%3d  (insuficiente)\n", nm, nrow(S))); return(invisible()) }
    sd_ <- apply(S, 2, sd); p <- 2 * (1 - pnorm(abs(est / sd_)))
    cat(sprintf("     %-16s n=%3d  EE %s   p %s\n", nm, nrow(S),
        paste(sprintf("%7.3f", sd_), collapse = " "),
        paste(sprintf("%6.4f", p), collapse = " ")))
  }
  linea("todas",        M)
  for (C in c(20, 10, 5, 3))
    linea(sprintf("|coef| <= %d", C), M[apply(abs(M) <= C, 1, all), , drop = FALSE])
  sd_iqr <- apply(M, 2, function(x) IQR(x) / 1.349)
  sd_mad <- apply(M, 2, mad)
  cat(sprintf("     %-16s        EE %s   p %s\n", "IQR/1.349", "",
      paste(sprintf("%6.4f", 2*(1-pnorm(abs(est/sd_iqr)))), collapse = " ")))
  cat(sprintf("     %-16s        EE %s   p %s\n", "  (valores)",
      paste(sprintf("%7.3f", sd_iqr), collapse = " "), ""))
  cat(sprintf("     %-16s        EE %s   p %s\n", "MAD*1.4826",
      paste(sprintf("%7.3f", sd_mad), collapse = " "),
      paste(sprintf("%6.4f", 2*(1-pnorm(abs(est/sd_mad)))), collapse = " ")))
}

for (m in names(COVS)) {
  cat("\n===============================", m, ":",
      paste(COVS[[m]], collapse = " + "), "===============================\n")
  fit <- ajustar1(COVS[[m]], cbind(time = tt, delta = dd, dat[, COVS[[m]], drop = FALSE]))
  t0  <- Sys.time()
  R   <- replicas(COVS[[m]], B)
  cat(sprintf("  %d replicas en %.0f s (%d dieron error y se descartaron)\n",
              nrow(R$b), as.numeric(difftime(Sys.time(), t0, units = "secs")), R$fallos))
  write.csv(cbind(R$b, R$beta), file.path(ROOT, sprintf("boot_%s.csv", m)), row.names = FALSE)
  comparar(R$b,    fit$b,    "INCIDENCIA (parte logistica)")
  comparar(R$beta, fit$beta, "LATENCIA (parte Cox)")
}

cat("\n=== PATRON ORO: bootstrap propio de smcure para M1, probando semillas ===\n")
covs <- COVS$M1
d <- cbind(time = tt, delta = dd, dat[, covs, drop = FALSE])
f_lat  <- as.formula(paste("Surv(time, delta) ~", paste(covs, collapse = " + ")))
f_cure <- as.formula(paste("~", paste(covs, collapse = " + ")))
for (s in 1:3) {
  set.seed(s); t0 <- Sys.time()
  r <- tryCatch(suppressWarnings({
        out <- NULL
        invisible(capture.output(
          out <- smcure::smcure(f_lat, cureform = f_cure, data = d, model = "ph",
                                Var = TRUE, nboot = 500)))
        out}), error = function(e) e)
  seg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) {
    cat(sprintf("  semilla %d: FALLA (%.0f s) -- %s\n", s, seg, conditionMessage(r)))
  } else {
    cat(sprintf("  semilla %d: OK (%.0f s)\n", s, seg))
    cat(sprintf("     EE incidencia smcure: %s\n", paste(sprintf("%7.3f", r$b_sd), collapse = " ")))
    cat(sprintf("     p  incidencia smcure: %s\n", paste(sprintf("%6.4f", r$b_pvalue), collapse = " ")))
    cat(sprintf("     EE latencia   smcure: %s\n", paste(sprintf("%7.3f", r$beta_sd), collapse = " ")))
    cat(sprintf("     p  latencia   smcure: %s\n", paste(sprintf("%6.4f", r$beta_pvalue), collapse = " ")))
    break
  }
}
cat("\nFIN.\n")
