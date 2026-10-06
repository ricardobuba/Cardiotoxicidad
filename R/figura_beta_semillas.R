# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# figura_beta_semillas.R - Correcciones A7 y A9.
#
# Sustituye la figura mod_sicure_beta.png (un solo beta(s), el de la semilla 1)
# por una que superpone los TRES ajustes cacheados. No reajusta nada: lee
# cache/sicure.rds. Tarda segundos.
#
# Motivo: los tres arranques aleatorios dan pesos funcionales cualitativamente
# distintos (dos de ellos con IDENTICO valor de la funcion objetivo, 211), de
# modo que presentar solo uno de ellos, con su interpretacion por tramos del
# ciclo, no es sostenible. La figura honesta es la que muestra la dispersion.
#
# Los coeficientes se toman directamente de par (ver verificar_indice_sicure.R):
#   theta = c(1, par[1:16]);  beta(s) = sum_k (theta_{9+k} / sd_k) phi_k(s)
# donde sd_k es la desviacion tipica del score k (sicure estandariza los
# scores internamente) y phi_k los armonicos de la FPCA.
#
# Uso:  Rscript R/figura_beta_semillas.R
# Salida: salidas/figuras/mod_sicure_beta_semillas.png
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))

K_FUN <- 8
CACHE <- file.path(ROOT, "cache", "sicure.rds")
stopifnot(file.exists(CACHE))

cat("== Preparacion de datos ==\n")
D    <- preparar_datos()
fits <- readRDS(CACHE)

sc_s  <- scale(D$pca$scores[D$ok, seq_len(K_FUN), drop = FALSE])
sd_sc <- attr(sc_s, "scaled:scale")
s     <- seq(0, 1, length.out = 1001)
arm   <- fda::eval.fd(s, D$pca$harmonics)[, seq_len(K_FUN), drop = FALSE]

B <- sapply(fits, function(f) {
  th <- c(1, f$par[1:16])
  as.numeric(arm %*% (th[10:(9 + K_FUN)] / sd_sc))
})
val <- sapply(fits, function(f) f$value)
cat("\nvalores de la funcion objetivo:", paste(round(val, 1), collapse = ", "), "\n")
for (j in seq_len(ncol(B)))
  cat(sprintf("  semilla %d: beta en [%.4f, %.4f], media %.4f, maximo en s = %.3f\n",
              j, min(B[, j]), max(B[, j]), mean(B[, j]), s[which.max(B[, j])]))

# Curva ECG media, para situar los tramos del ciclo (panel inferior).
ecg_medio <- colMeans(as.matrix(D$E)[D$ok, ])
s_ecg     <- seq(0, 1, length.out = length(ecg_medio))

COL <- c(PINK, BLUE, GRAY)
figura("mod_sicure_beta_semillas.png", 8.3, 6.7, {
  op <- par(mfrow = c(2, 1), mar = c(4.2, 4.6, 2.2, 1.2), mgp = c(2.6, 0.8, 0))
  matplot(s, B, type = "l", lty = 1, lwd = 2.2, col = COL,
          xlab = "", ylab = expression(hat(beta)(s)),
          main = "Peso funcional estimado según el punto de arranque")
  abline(h = 0, col = "grey70", lty = 3)
  legend("topright", bty = "n", lwd = 2.2, col = COL, cex = 0.85,
         legend = sprintf("arranque %d  (objetivo = %.0f)", seq_along(val), val))
  plot(s_ecg, ecg_medio, type = "l", lwd = 1.8, col = "grey25",
       xlab = "s (ciclo cardíaco normalizado)", ylab = "ECG medio",
       main = "")
  par(op)
})

cat("\nSube a Overleaf (imaxes/modelizacion):\n  ",
    file.path(SAL, "figuras", "mod_sicure_beta_semillas.png"), "\n")
