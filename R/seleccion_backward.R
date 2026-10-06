# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# seleccion_backward.R - Correccion A1 de la tutora (10/09/2026):
#
#   "Para ver que variables incluir en el modelo semiparametrico, puedes
#    incluir todas e ir eliminando, paso a paso, las menos significativas.
#    Yo votaria por la primera opcion: partes del modelo con todas las
#    variables clinicas y vas eliminando"
#
# Sustituye al cribado univariante (tab_mod_screening) como criterio de
# seleccion de M1.
#
# NO modifica tfg_comun.R ni modelizacion_curacion.R: solo los lee.
#
# Decisiones de diseno (acordadas 10/09/2026):
#   * Backward DESACOPLADO: smcure admite una formula para la latencia y otra
#     (cureform) para la incidencia, de modo que en cada paso se elimina la
#     variable con mayor p-valor de ENTRE las 18 candidatas (9 en incidencia
#     + 9 en latencia) y se elimina SOLO de esa componente. Asi, el resultado
#     "el efecto se concentra en la incidencia y no en la latencia" sale del
#     propio procedimiento en lugar de afirmarse a mano.
#   * Criterio: p-valor de Wald con error estandar bootstrap. Durante la
#     busqueda se usa un bootstrap reducido (B_BUSQUEDA), porque ahi el
#     bootstrap solo ORDENA candidatas; el modelo final se reajusta con
#     B_FINAL = NBOOT (500) y es de ese reajuste de donde salen todas las
#     cifras que van a la memoria.
#   * Salida cuando todos los p < ALPHA_SALIDA.
#
# Uso:   Rscript R/seleccion_backward.R
# Salidas:
#   salidas/tablas/tab_mod_backward.tex   (traza de la eliminacion)
#   salidas/tablas/tab_mod_m1b.tex        (modelo final)
#   cache/backward_clinicas.rds                    (ajuste final + historia)
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
stopifnot(file.exists(file.path(ROOT, "R", "tfg_comun.R")))
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(smcure))

B_BUSQUEDA    <- 100          # replicas bootstrap en cada paso de la busqueda
B_FINAL       <- NBOOT        # replicas del reajuste final (500)
ALPHA_SALIDA  <- 0.05
SEED_BUSQUEDA <- 20260910

dir.create(file.path(ROOT, "cache"), showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Ajuste con covariables DISTINTAS en incidencia y en latencia
# ---------------------------------------------------------------------------
# Copia adaptada de tfg_comun.R::ajustar_curacion. Unico cambio de fondo: dos
# conjuntos de covariables en lugar de uno. El bootstrap tolerante a replicas
# divergentes (criterio MAX_COEF_BOOT) es identico, para que las cifras sean
# comparables con las de modelizacion_curacion.R.
ajustar_curacion2 <- function(covs_inc, covs_lat, datos, time, delta,
                              nboot = B_FINAL, seed = SEED, Var = TRUE,
                              max_intentos = 4L) {
  stopifnot(length(covs_lat) >= 1L)   # smcure exige latencia no vacia
  usadas <- union(covs_inc, covs_lat)
  d <- cbind(time = time, delta = delta, datos[, usadas, drop = FALSE])
  f_lat  <- as.formula(paste("Surv(time, delta) ~", paste(covs_lat, collapse = " + ")))
  f_cure <- as.formula(paste("~", paste(covs_inc, collapse = " + ")))

  ajuste1 <- function(dd) {
    r <- NULL
    invisible(utils::capture.output(
      r <- smcure::smcure(f_lat, cureform = f_cure, data = dd,
                          model = "ph", Var = FALSE)))
    r
  }

  fit <- tryCatch(suppressWarnings(ajuste1(d)), error = function(e) e)
  if (inherits(fit, "error")) return(structure(list(error = conditionMessage(fit)),
                                               class = "curacion_fallida"))
  fit$covs_inc <- covs_inc; fit$covs_lat <- covs_lat
  if (!isTRUE(Var)) return(fit)

  set.seed(seed)
  i1 <- which(delta == 1); i0 <- which(delta == 0)
  Bb <- matrix(NA_real_, nboot, length(fit$b))
  Bg <- matrix(NA_real_, nboot, length(fit$beta))
  ok <- 0L; intentos <- 0L; tope <- as.integer(max_intentos) * nboot
  while (ok < nboot && intentos < tope) {
    intentos <- intentos + 1L
    idx <- c(sample(i1, length(i1), replace = TRUE),
             sample(i0, length(i0), replace = TRUE))
    r <- tryCatch(suppressWarnings(ajuste1(d[idx, , drop = FALSE])),
                  error = function(e) NULL)
    if (is.null(r) || !all(is.finite(r$b)) || !all(is.finite(r$beta))) next
    if (max(abs(r$b)) > MAX_COEF_BOOT || max(abs(r$beta)) > MAX_COEF_BOOT) next
    ok <- ok + 1L
    Bb[ok, ] <- r$b; Bg[ok, ] <- r$beta
  }
  if (ok < 2L) return(structure(list(error = sprintf(
      "bootstrap imposible: %d replicas validas en %d intentos", ok, intentos)),
      class = "curacion_fallida"))
  Bb <- Bb[seq_len(ok), , drop = FALSE]; Bg <- Bg[seq_len(ok), , drop = FALSE]

  fit$b_sd        <- sqrt(apply(Bb, 2, var))
  fit$b_zvalue    <- fit$b / fit$b_sd
  fit$b_pvalue    <- (1 - pnorm(abs(fit$b_zvalue))) * 2
  fit$beta_sd     <- sqrt(apply(Bg, 2, var))
  fit$beta_zvalue <- fit$beta / fit$beta_sd
  fit$beta_pvalue <- (1 - pnorm(abs(fit$beta_zvalue))) * 2
  names(fit$b_sd) <- names(fit$b_pvalue) <- names(fit$b)
  names(fit$beta_sd) <- names(fit$beta_pvalue) <- names(fit$beta)
  fit$nboot_validas     <- ok
  fit$nboot_descartadas <- intentos - ok
  fit
}

fallo <- function(x) inherits(x, "curacion_fallida")

# Candidatas de un ajuste: p-valores de incidencia (sin intercepto) y latencia.
candidatas <- function(fit, covs_inc, covs_lat) {
  ii <- which(fit$bnm %in% covs_inc)         # excluye "(Intercept)"
  rbind(
    data.frame(comp = "incidencia", var = fit$bnm[ii],
               coef = as.numeric(fit$b[ii]), p = as.numeric(fit$b_pvalue[ii]),
               stringsAsFactors = FALSE),
    data.frame(comp = "latencia", var = fit$betanm,
               coef = as.numeric(fit$beta), p = as.numeric(fit$beta_pvalue),
               stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------------------
# 2. Eliminacion hacia atras
# ---------------------------------------------------------------------------
backward_curacion <- function(datos, time, delta,
                              vars = VARS_CLIN,
                              alpha = ALPHA_SALIDA,
                              nboot = B_BUSQUEDA,
                              seed = SEED_BUSQUEDA,
                              min_lat = 1L) {
  covs_inc <- vars; covs_lat <- vars
  historia <- list(); paso <- 0L

  repeat {
    paso <- paso + 1L
    cat(sprintf("\n-- paso %d: incidencia (%d) = %s\n              latencia   (%d) = %s\n",
                paso, length(covs_inc), paste(covs_inc, collapse = ", "),
                length(covs_lat), paste(covs_lat, collapse = ", ")))
    t0 <- Sys.time()
    fit <- ajustar_curacion2(covs_inc, covs_lat, datos, time, delta,
                             nboot = nboot, seed = seed, Var = TRUE)
    cat(sprintf("   (%.1f s)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    if (fallo(fit)) {
      # El modelo de este paso no converge. Con 25 eventos esto puede pasar en
      # los primeros pasos, con muchas covariables en la incidencia. Se elimina
      # la variable con menor |t| del ultimo ajuste que SI convergio; si no hay
      # ninguno (falla el modelo completo), se elimina de la incidencia la peor
      # segun un Cox univariante, que es la unica informacion disponible.
      cat("   [!] no converge:", fit$error, "\n")
      if (length(historia) > 0L) {
        prev <- historia[[length(historia)]]$cand
      } else {
        z <- sapply(covs_inc, function(v)
          abs(summary(coxph(Surv(time, delta) ~ datos[[v]]))$coefficients[, "z"]))
        prev <- data.frame(comp = "incidencia", var = covs_inc,
                           coef = NA_real_, p = 1 - stats::pnorm(z),
                           stringsAsFactors = FALSE)
      }
      peor <- prev[which.max(prev$p), ]
      cat(sprintf("   -> se elimina %s de %s (rescate por no convergencia)\n",
                  peor$var, peor$comp))
      if (peor$comp == "incidencia") covs_inc <- setdiff(covs_inc, peor$var)
      else                            covs_lat <- setdiff(covs_lat, peor$var)
      historia[[paso]] <- list(paso = paso, covs_inc = covs_inc, covs_lat = covs_lat,
                               cand = prev, eliminada = peor, converge = FALSE)
      if (length(covs_inc) == 0L) stop("la incidencia se ha quedado vacia")
      next
    }

    cand <- candidatas(fit, covs_inc, covs_lat)
    print(cand[order(-cand$p), ], digits = 3, row.names = FALSE)

    # No se puede vaciar la latencia: smcure exige al menos una covariable.
    elegibles <- cand
    if (length(covs_lat) <= min_lat)
      elegibles <- elegibles[elegibles$comp != "latencia", ]

    peor <- elegibles[which.max(elegibles$p), ]
    if (nrow(elegibles) == 0L || peor$p < alpha) {
      cat("\n== parada: todos los p < ", alpha, " (o no quedan candidatas) ==\n", sep = "")
      historia[[paso]] <- list(paso = paso, covs_inc = covs_inc, covs_lat = covs_lat,
                               cand = cand, eliminada = NULL, converge = TRUE)
      break
    }
    cat(sprintf("   -> se elimina %s de %s (p = %.3f)\n", peor$var, peor$comp, peor$p))
    if (peor$comp == "incidencia") covs_inc <- setdiff(covs_inc, peor$var)
    else                            covs_lat <- setdiff(covs_lat, peor$var)
    historia[[paso]] <- list(paso = paso, covs_inc = covs_inc, covs_lat = covs_lat,
                             cand = cand, eliminada = peor, converge = TRUE)
    if (length(covs_inc) == 0L) { cat("\n== la incidencia se ha quedado vacia ==\n"); break }
  }
  list(historia = historia, covs_inc = covs_inc, covs_lat = covs_lat)
}

# ---------------------------------------------------------------------------
# 3. Ejecucion
# ---------------------------------------------------------------------------
cat("== Preparacion de datos ==\n")
D <- preparar_datos()
dat <- D$datos; tt <- D$time; dd <- D$delta

cat("\n== Backward desacoplado con smcure (B =", B_BUSQUEDA, "en la busqueda) ==\n")
t_ini <- Sys.time()
bw <- backward_curacion(dat, tt, dd)
cat(sprintf("\nBusqueda terminada en %.1f min\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))
cat("Modelo seleccionado:\n  incidencia:", paste(bw$covs_inc, collapse = " + "),
    "\n  latencia  :", paste(bw$covs_lat, collapse = " + "), "\n")

cat("\n== Reajuste final con B =", B_FINAL, "==\n")
M1b <- ajustar_curacion2(bw$covs_inc, bw$covs_lat, dat, tt, dd,
                         nboot = B_FINAL, seed = SEED, Var = TRUE)
stopifnot(!fallo(M1b))
cat(sprintf("  bootstrap: %d replicas validas, %d descartadas por divergencia\n",
            M1b$nboot_validas, M1b$nboot_descartadas))
frac <- 1 - mean(1 / (1 + exp(-drop(cbind(1, as.matrix(dat[, bw$covs_inc, drop = FALSE])) %*% M1b$b))))
cat(sprintf("  fraccion de curacion media estimada: %.3f\n", frac))
print(candidatas(M1b, bw$covs_inc, bw$covs_lat), digits = 3, row.names = FALSE)

# ------------------------- Tabla: traza de la eliminacion ------------------
et <- function(v) ifelse(is.na(ETIQ[v]), v, ETIQ[v])
filas <- c()
for (h in bw$historia) {
  if (is.null(h$eliminada)) next
  filas <- c(filas, sprintf("%d & %d & %d & %s & %s & %s \\\\",
    h$paso,
    length(h$covs_inc) + ifelse(h$eliminada$comp == "incidencia", 1L, 0L),
    length(h$covs_lat) + ifelse(h$eliminada$comp == "latencia", 1L, 0L),
    et(h$eliminada$var),
    ifelse(h$eliminada$comp == "incidencia", "Incidencia", "Latencia"),
    if (isTRUE(h$converge)) fmt_p(h$eliminada$p) else "---"))
}
escribir_tabla_tex("tab_mod_backward.tex", "cccllc",
  paste("\\textbf{Paso} & \\textbf{$k_{\\text{inc}}$} & \\textbf{$k_{\\text{lat}}$}",
        "& \\textbf{Variable eliminada} & \\textbf{Componente} & \\textbf{$p$}"),
  filas,
  sprintf("Eliminación hacia atrás sobre el modelo de mixtura de curación, partiendo de las nueve variables clínicas en las dos componentes ($n=%d$, %d eventos). En cada paso se elimina la variable con mayor $p$-valor de entre todas las candidatas, y sólo de la componente en la que aparece. Los $p$-valores de la búsqueda se obtienen con un \\textit{bootstrap} reducido de $B=%d$ réplicas; el modelo final se reajusta con $B=%d$.",
          length(tt), sum(dd), B_BUSQUEDA, B_FINAL),
  "tab:mod-backward",
  nota = "$k_{\\text{inc}}$ y $k_{\\text{lat}}$: número de covariables en incidencia y en latencia antes de la eliminación.")

# ------------------------- Tabla: modelo final -----------------------------
filas <- c()
for (i in seq_along(M1b$b))
  filas <- c(filas, sprintf("%s & Incidencia & %s & %s & [%s; %s] & %s \\\\",
    et(M1b$bnm[i]), fmt_n(M1b$b[i], 3), fmt_n(exp(M1b$b[i]), 3),
    fmt_n(exp(M1b$b[i] - 1.96 * M1b$b_sd[i]), 3),
    fmt_n(exp(M1b$b[i] + 1.96 * M1b$b_sd[i]), 3), fmt_p(M1b$b_pvalue[i])))
for (i in seq_along(M1b$beta))
  filas <- c(filas, sprintf("%s & Latencia & %s & %s & [%s; %s] & %s \\\\",
    et(M1b$betanm[i]), fmt_n(M1b$beta[i], 3), fmt_n(exp(M1b$beta[i]), 3),
    fmt_n(exp(M1b$beta[i] - 1.96 * M1b$beta_sd[i]), 3),
    fmt_n(exp(M1b$beta[i] + 1.96 * M1b$beta_sd[i]), 3), fmt_p(M1b$beta_pvalue[i])))
escribir_tabla_tex("tab_mod_m1b.tex", "llcccc",
  "\\textbf{Covariable} & \\textbf{Componente} & \\textbf{Coef.} & \\textbf{OR / HR} & \\textbf{IC 95\\,\\%} & \\textbf{$p$}",
  filas,
  sprintf("Modelo de mixtura de curación con variables clínicas seleccionadas por eliminación hacia atrás ($n=%d$, %d eventos, $B=%d$). Fracción de curación media estimada: $%.3f$.",
          length(tt), sum(dd), B_FINAL, frac),
  "tab:mod-m1b")

saveRDS(list(historia = bw$historia, covs_inc = bw$covs_inc, covs_lat = bw$covs_lat,
             fit = M1b, fraccion = frac,
             B_busqueda = B_BUSQUEDA, B_final = B_FINAL, alpha = ALPHA_SALIDA),
        file.path(ROOT, "cache", "backward_clinicas.rds"))
cat("\nTODO OK -> cache/backward_clinicas.rds\n")
