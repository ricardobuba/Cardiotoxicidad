# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# verificar_indice_sicure.R - Correccion A7 de la tutora ("Esto no es
# correcto", sobre el R2 = 1,000 de la recuperacion del indice).
#
# QUE DICE EL CODIGO FUENTE DE sicure (CRAN, R/sicure.R):
#
#   * sicure.vf estandariza INTERNAMENTE las covariables vectoriales y los
#     scores funcionales:   x_cov <- cbind(scale(x_cov_v), scale(scores))
#   * la FPCA interna es fda::pca.fd sobre x_cov_f, primero con nharm = 20
#     para medir la variabilidad y luego con el minimo nharm que alcanza
#     propvar (0,90).
#   * la restriccion de identificacion fija el PRIMER coeficiente a 1, y el
#     indice es, literalmente:
#         si <- colSums(c(1, par[1:d]) * t(x_cov))
#   * par tiene longitud d + 4: los CUATRO ULTIMOS elementos son los
#     logaritmos de los anchos de banda. En el ajuste cacheado par tiene
#     longitud 20, luego d = 16 = 9 clinicas + K funcionales - 1  =>  K = 8.
#     (Confirma que K_FUN = 8 es correcto.)
#
# CONSECUENCIA: el vector de coeficientes del indice, en la escala de las
# covariables estandarizadas, es EXACTAMENTE  c(1, par[1:16]).  No hay nada
# que "recuperar": regresar si sobre esas mismas 17 variables da R2 = 1 por
# construccion, y por eso la tutora marca la frase como incorrecta. El R2 no
# verificaba nada; era una tautologia.
#
# Este guion (a) comprueba empiricamente que la identidad se cumple, que es
# lo unico que hay que verificar de verdad, porque depende de que TU FPCA y
# la interna de sicure sean la misma; y (b) calcula beta(s) directamente a
# partir de par, sin regresion.
#
# OJO AL SIGNO: cada armonico de la FPCA esta definido salvo signo. Si tu
# pca.fd y el interno de sicure devolvieran armonicos con signos distintos,
# el beta(s) actual estaria mal. La regresion por minimos cuadrados absorbia
# ese posible cambio de signo sin que se notara: de ahi que "funcionase".
# El contraste de abajo lo detecta.
#
# Uso:  Rscript R/diagnostico/verificar_indice_sicure.R
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))

K_FUN  <- 8
CACHE  <- file.path(ROOT, "cache", "sicure.rds")
stopifnot(file.exists(CACHE))

cat("== Datos ==\n")
D <- preparar_datos()

fits <- readRDS(CACHE)
cat("\nAjustes cacheados:", length(fits), "\n")

# Matriz de covariables tal y como la construye sicure.vf internamente:
# scale() de las 9 clinicas y scale() de los K scores funcionales.
clin_s   <- scale(as.matrix(D$clin_fun[D$ok, VARS_CLIN]))
sc_raw   <- D$pca$scores[D$ok, seq_len(K_FUN), drop = FALSE]
sc_s     <- scale(sc_raw)
sd_sc    <- attr(sc_s, "scaled:scale")
X        <- cbind(clin_s, sc_s)
stopifnot(ncol(X) == 9 + K_FUN)

for (j in seq_along(fits)) {
  f  <- fits[[j]]
  d  <- length(f$par) - 4L                 # los 4 ultimos son log-anchos
  cat(sprintf("\n---------- ajuste %d  (par: %d = %d coef + 4 log-h) ----------\n",
              j, length(f$par), d))
  if (d != ncol(X) - 1L) {
    cat("  [!] d =", d, "no coincide con ncol(X)-1 =", ncol(X) - 1L,
        "-> K_FUN esta mal; sicure retuvo", d + 1L - 9L, "componentes\n")
    next
  }
  theta <- c(1, f$par[1:d])                # coeficientes en escala estandarizada
  si_re <- as.numeric(X %*% theta)

  err   <- max(abs(si_re - as.numeric(f$si)))
  esc   <- diff(range(as.numeric(f$si)))
  r     <- cor(si_re, as.numeric(f$si))
  cat(sprintf("  max |si_reconstruido - si_sicure| = %.3e   (rango de si: %.3f)\n", err, esc))
  cat(sprintf("  correlacion                        = %.8f\n", r))

  if (err < 1e-6 * max(1, esc)) {
    cat("  => IDENTIDAD EXACTA. Tu FPCA y la interna de sicure coinciden,\n")
    cat("     armonico a armonico y con el mismo signo. beta(s) se puede\n")
    cat("     calcular directamente desde par: no hace falta la regresion.\n")
  } else if (r > 0.999) {
    cat("  => casi identico pero no exacto: revisa el orden de las columnas\n")
    cat("     o la version de scale(); la reconstruccion es utilizable.\n")
  } else {
    cat("  [!] NO COINCIDE. Muy probablemente algun armonico tiene el signo\n")
    cat("      cambiado respecto al de sicure, o la FPCA interna usa otra\n")
    cat("      base. Comprueba los signos uno a uno:\n")
    a <- theta[10:(9 + K_FUN)]
    for (k in seq_len(K_FUN)) {
      th2 <- theta; th2[9 + k] <- -th2[9 + k]
      cat(sprintf("      cambiando el signo de FPC%d -> corr = %.6f\n",
                  k, cor(as.numeric(X %*% th2), as.numeric(f$si))))
    }
  }

  # ---- beta(s) directamente desde par -------------------------------------
  # El indice usa scores ESTANDARIZADOS: z_k = (xi_k - m_k)/s_k. Como
  # xi_k = <X - Xbar, phi_k>, la parte funcional del indice es
  #   sum_k a_k z_k = < X - Xbar , sum_k (a_k/s_k) phi_k >,
  # de modo que beta(s) = sum_k (a_k / s_k) phi_k(s).
  # (El codigo actual, al regresar sobre los scores SIN estandarizar, obtiene
  # directamente a_k/s_k y por eso no dividia; con par SI hay que dividir.)
  a_std <- theta[10:(9 + K_FUN)]
  a_ori <- a_std / sd_sc
  s     <- seq(0, 1, length.out = 1001)
  arm   <- fda::eval.fd(s, D$pca$harmonics)[, seq_len(K_FUN), drop = FALSE]
  beta  <- as.numeric(arm %*% a_ori)
  cat(sprintf("  beta(s): min = %.5f  max = %.5f  media = %.5f  (max/min = %.2f)\n",
              min(beta), max(beta), mean(beta), max(beta) / min(beta)))
  cat(sprintf("  cruza su media en s = %.3f y s = %.3f ; maximo en s = %.3f\n",
              s[which(diff(sign(beta - mean(beta))) != 0)[1]],
              s[rev(which(diff(sign(beta - mean(beta))) != 0))[1]],
              s[which.max(beta)]))

  # Comparacion con el metodo antiguo (la regresion), solo para documentar
  # la diferencia en la memoria.
  aj    <- lm(f$si ~ cbind(clin_s, sc_raw))
  a_lm  <- coef(aj)[11:(10 + K_FUN)]
  beta0 <- as.numeric(arm %*% a_lm)
  cat(sprintf("  [metodo antiguo, regresion] R2 = %.6f ; correlacion entre ambos beta(s) = %.6f ; max|dif| = %.2e\n",
              summary(aj)$r.squared, cor(beta, beta0), max(abs(beta - beta0))))
  if (j == 1L) saveRDS(list(s = s, beta = beta, theta = theta,
                            clinicos = theta[1:9], funcionales = a_ori),
                       file.path(ROOT, "cache", "beta_sicure_directo.rds"))
}
cat("\nHecho.\n")
