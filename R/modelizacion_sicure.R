# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# modelizacion_sicure.R - Modelo single-index (sicure) sobre el subconjunto
# con ECG. Reproduce TODAS las cifras del apartado 5.4.2 de la memoria y
# genera la figura del peso funcional.
#
# Un ajuste tarda unos 9 minutos y se hacen tres, asi que el resultado se
# cachea en cache/sicure.rds. Si se toca el preprocesado hay que borrarlo a
# mano: el cache no detecta el cambio.
#
# Uso:  Rscript R/modelizacion_sicure.R
# Salidas: salidas/figuras/mod_sicure_beta.png
#          y por consola las cifras que cita el Capitulo 5.
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(sicure))

PROPVAR  <- 0.90                 # varianza retenida -> 8 componentes
SEMILLAS <- 1:3                  # puntos de arranque de la busqueda aleatoria
K_FUN    <- 8                    # componentes que retiene PROPVAR
# OJO: con el suavizado por GCV (30/08/2026) el 90 % de varianza ya no se
# alcanza con 6 componentes sino con 8. Verificar en la primera ejecucion:
# si K_FUN no coincide con lo que retiene sicure.vf, coeficientes_indice()
# avisa porque el R2 de la recuperacion del indice baja de 0,99.
CACHE    <- file.path(ROOT, "cache", "sicure.rds")

# --------------------------------- Datos -----------------------------------

# E y clin_fun contienen las 183 pacientes con curva ECG valida; ok recorta a
# las 181 sin valores clinicos ausentes, que son las que usan todos los
# modelos del capitulo.  (Antes del cambio de suavizado eran 185 y 183.)
datos_sicure <- function(D) {
  stopifnot(sum(D$ok) == 181, nrow(D$E) == 183, nrow(D$clin_fun) == 183)
  dat <- list(curvas = as.matrix(D$E)[D$ok, ],
              clin   = scale(as.matrix(D$clin_fun[D$ok, VARS_CLIN])),
              tiempo = D$time,
              delta  = D$delta)
  stopifnot(nrow(dat$curvas) == nrow(dat$clin),
            length(dat$tiempo) == nrow(dat$clin),
            !anyNA(dat$clin), !anyNA(dat$delta))
  dat
}

ajustar_sicure <- function(dat, semilla = 1) {
  set.seed(semilla)
  sicure.vf(x_cov_v = dat$clin, x_cov_f = dat$curvas, time = dat$tiempo,
            delta = dat$delta, propvar = PROPVAR, randomsearch = TRUE)
}

# --------------------- Coeficientes del indice y beta(s) -------------------

# `par` no da los coeficientes en la escala original: omite el que la
# restriccion de identificacion fija a 1 y expresa la parte funcional sobre los
# scores estandarizados. Como el indice es, por construccion, una combinacion
# lineal de las 9 clinicas y los K scores, se recuperan sin aproximacion por
# minimos cuadrados a partir de si; el R2 lo verifica (sale 0,9998).
coeficientes_indice <- function(fit, dat, D, K = K_FUN) {
  sc <- D$pca$scores[D$ok, seq_len(K)]
  aj <- lm(fit$si ~ cbind(dat$clin, sc))
  r2 <- summary(aj)$r.squared
  if (r2 < 0.99) warning("recuperacion del indice imprecisa: R2 = ", round(r2, 4))
  list(r2 = r2,
       clinicos    = coef(aj)[2:10],
       funcionales = coef(aj)[11:(10 + K)])
}

# Peso funcional beta(s): los coeficientes funcionales por los armonicos.
beta_funcional <- function(coefs, D, m = 1001) {
  s   <- seq(0, 1, length.out = m)
  arm <- fda::eval.fd(s, D$pca$harmonics)[, seq_along(coefs$funcionales)]
  list(s = s, beta = as.numeric(arm %*% coefs$funcionales))
}

# Se representa beta(s) centrado: el peso es positivo en todo el ciclo y su
# variacion es pequena frente a su nivel, de modo que sin centrar parece plano.
# El panel inferior anade el ECG medio para situar los tramos.
figura_beta <- function(bf, dat, fichero = "mod_sicure_beta.png") {
  openfig(fichero, 8.3, 6.7)
  op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 1))
  plot(bf$s, bf$beta - mean(bf$beta), type = "l", lwd = 2, xlab = "",
       ylab = expression(hat(beta)(s) - bar(beta)))
  abline(h = 0, lty = 3); grid()
  plot(bf$s, colMeans(dat$curvas), type = "l", lwd = 2, col = GRAY,
       xlab = "Fracción del ciclo cardíaco", ylab = "ECG medio")
  grid()
  par(op); done(fichero)
}

# ------------------- Resultados que cita el Capitulo 5 ---------------------

resultados_sicure <- function(D, usar_cache = TRUE) {
  dat <- datos_sicure(D)

  if (usar_cache && file.exists(CACHE)) {
    cat("  (usando", CACHE, "-- borralo si has tocado el preprocesado)\n")
    aj <- readRDS(CACHE)
  } else {
    aj <- lapply(SEMILLAS, function(s) {
      cat("  ajustando semilla", s, "...\n"); ajustar_sicure(dat, s)
    })
    names(aj) <- SEMILLAS
    dir.create(dirname(CACHE), recursive = TRUE, showWarnings = FALSE)
    saveRDS(aj, CACHE)
  }
  fit <- aj[["1"]]

  # El indice frente a las clinicas: a cual sigue, y cuanto explican en conjunto
  cor_clin <- setNames(round(as.numeric(cor(fit$si, dat$clin)), 3), VARS_CLIN)
  r2_clin  <- summary(lm(fit$si ~ dat$clin))$r.squared

  # Estabilidad entre puntos de arranque de la busqueda aleatoria
  cor_semillas <- min(cor(sapply(aj, function(f) f$si)))

  # Separacion por terciles del indice (en la muestra de ajuste: p optimista)
  g  <- cut(fit$si, quantile(fit$si, 0:3 / 3), include.lowest = TRUE,
            labels = 1:3)
  lr <- survdiff(Surv(dat$tiempo, dat$delta) ~ g)
  p_logrank <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

  coefs <- coeficientes_indice(fit, dat, D)
  bf    <- beta_funcional(coefs, D)
  figura_beta(bf, dat)

  list(n = nrow(dat$clin), eventos = sum(dat$delta),
       cor_clinicas = cor_clin, r2_clinicas = r2_clin,
       r2_indice = coefs$r2, cor_semillas = cor_semillas,
       p_logrank = p_logrank,
       beta_media = mean(bf$beta),
       beta_min = min(bf$beta), beta_max = max(bf$beta),
       s_min = bf$s[which.min(bf$beta)],
       cruce_beta = bf$s[which.min(abs(bf$beta - mean(bf$beta)))])
}

# ------------------------------- Ejecucion ---------------------------------

if (sys.nframe() == 0) {
  D <- preparar_datos()
  R <- resultados_sicure(D)

  cat("\n--- Cifras citadas en el Capitulo 5 ---\n")
  cat("n =", R$n, " eventos =", R$eventos, "\n")
  cat("correlacion del indice con las clinicas:\n"); print(R$cor_clinicas)
  cat("las 9 clinicas explican del indice:", round(100 * R$r2_clinicas, 1), "%\n")
  cat("resto, atribuible a la curva:       ", round(100 * (1 - R$r2_clinicas), 1), "%\n")
  cat("correlacion minima entre semillas:  ", round(R$cor_semillas, 3), "\n")
  cat("log-rank por terciles del indice:    p =", signif(R$p_logrank, 3), "\n")
  cat("beta(s): media", round(R$beta_media, 4),
      "| max", round(R$beta_max, 4),
      "| min", round(R$beta_min, 4), "en s =", round(R$s_min, 3), "\n")
  cat("beta(s) cruza su media en s =       ", round(R$cruce_beta, 3), "\n")
  cat("verificacion de la recuperacion: R2 =", round(R$r2_indice, 4), "\n")
}
