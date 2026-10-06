# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# curacion_lib.R - Funciones compartidas por los guiones de correccion de
# septiembre de 2026 (A1, A2, A3). Solo definiciones: no ejecuta nada.
#
# Contiene las funciones de seleccion_backward.R, con un unico anadido:
# backward_curacion() acepta ahora conjuntos de arranque distintos en
# incidencia y en latencia (vars_inc / vars_lat). Por defecto se comporta
# exactamente igual que la version de seleccion_backward.R.
#   ajustar_curacion2()  - smcure con covariables DISTINTAS en incidencia y
#                          latencia, con el bootstrap tolerante a replicas
#                          divergentes de tfg_comun.R::ajustar_curacion.
#   candidatas()         - p-valores de las dos componentes.
#   backward_curacion()  - eliminacion hacia atras desacoplada.
#
# Requiere que tfg_comun.R este ya cargado (NBOOT, MAX_COEF_BOOT, SEED,
# VARS_CLIN) y la libreria smcure.
#
# NOTA: seleccion_backward.R lleva hoy su propia copia de estas funciones,
# porque se escribio antes de que hubiera un segundo guion que las usara.
# Cuando termine la ejecucion en curso conviene sustituir alli el bloque 1-2
# por  source(file.path(ROOT, "R", "curacion_lib.R"))  para que no haya dos copias.
# ---------------------------------------------------------------------------

if (!exists("B_BUSQUEDA"))   B_BUSQUEDA   <- 100
if (!exists("B_FINAL"))      B_FINAL      <- NBOOT
if (!exists("ALPHA_SALIDA")) ALPHA_SALIDA <- 0.05

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
                              min_lat = 1L, min_inc = 1L,
                              vars_inc = vars, vars_lat = vars) {
  # vars_inc / vars_lat permiten arrancar con conjuntos DISTINTOS en cada
  # componente. Se usa en M3, que parte de las clinicas que sobrevivieron en
  # cada componente por separado en A1, y no de la union de ambas: una
  # variable descartada de la incidencia no debe reaparecer en ella.
  covs_inc <- vars_inc; covs_lat <- vars_lat
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
    # Ni la incidencia ni la latencia pueden quedarse vacias: smcure exige al
    # menos una covariable en cada componente. Cuando solo queda una, deja de
    # ser candidata a eliminacion y se conserva aunque no sea significativa;
    # eso hay que declararlo en la memoria (es lo que ocurre con la FEVI en la
    # latencia de M1 y con la ultima componente funcional en M2).
    elegibles <- cand
    if (length(covs_lat) <= min_lat)
      elegibles <- elegibles[elegibles$comp != "latencia", ]
    if (length(covs_inc) <= min_inc)
      elegibles <- elegibles[elegibles$comp != "incidencia", ]
    if (nrow(elegibles) == 0L) {
      cat("\n== parada: solo queda una covariable en cada componente ==\n")
      historia[[paso]] <- list(paso = paso, covs_inc = covs_inc, covs_lat = covs_lat,
                               cand = cand, eliminada = NULL, converge = TRUE)
      break
    }

    peor <- elegibles[which.max(elegibles$p), ]
    if (peor$p < alpha) {
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
