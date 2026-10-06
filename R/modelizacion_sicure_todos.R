# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# modelizacion_sicure_todos.R - Correccion A3 / reestructuracion del Cap. 5
# pedida por la tutora (TEAMS, 10/09/2026):
#
#   5.3 Modelos de curacion con variables clinicas   -> smcure Y sicure
#   5.4 Modelos con las componentes funcionales      -> smcure Y sicure
#   5.5 Modelos con todas las variables              -> smcure Y sicure
#
# La parte de smcure la cubren seleccion_backward.R (5.3) y
# modelizacion_funcional.R (5.4 y 5.5). La parte de sicure:
#
#   5.3 -> sicure.v(x_cov = 9 clinicas)          [ESTE GUION]
#   5.4 -> sicure.f(x_cov = curvas ECG)          [ESTE GUION]
#   5.5 -> sicure.vf(clinicas + curvas)          [YA ESTA: cache/sicure.rds]
#
# Es decir, de los tres ajustes de sicure que pide, uno ya lo tienes hecho.
#
# Cada ajuste se repite con SEMILLAS puntos de arranque distintos, porque
# verificar_indice_sicure.R demostro que la direccion del indice no queda
# identificada: dos arranques con el mismo valor objetivo dan pesos
# funcionales opuestos. La estabilidad hay que medirla, no suponerla.
#
# Uso:  Rscript R/modelizacion_sicure_todos.R
# Salidas: cache/sicure_v.rds, cache/sicure_f.rds
#          salidas/tablas/tab_mod_sicure.tex
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))
suppressPackageStartupMessages(library(sicure))

SEMILLAS <- 1:2        # subir a 1:3 si sobra tiempo
PROPVAR  <- 0.90
dir.create(file.path(ROOT, "cache"), showWarnings = FALSE)

cat("== Preparacion de datos ==\n")
D   <- preparar_datos()
ok  <- D$ok
tt  <- D$time; dd <- D$delta
clin   <- scale(as.matrix(D$clin_fun[ok, VARS_CLIN]))
curvas <- as.matrix(D$E)[ok, ]
stopifnot(nrow(clin) == length(tt), nrow(curvas) == length(tt))
cat(sprintf("  n = %d, eventos = %d, curvas = %d x %d\n",
            length(tt), sum(dd), nrow(curvas), ncol(curvas)))

ajustar <- function(etiqueta, fn, cache, ...) {
  f <- file.path(ROOT, "cache", cache)
  if (file.exists(f)) { cat("\n[", etiqueta, "] cacheado en ", cache, "\n", sep = "")
                        return(readRDS(f)) }
  cat("\n########", etiqueta, "########\n")
  res <- list()
  for (s in SEMILLAS) {
    cat("  semilla", s, "... ")
    t0 <- Sys.time()
    set.seed(s)
    res[[as.character(s)]] <- suppressWarnings(fn(...))
    cat(sprintf("%.1f min (objetivo = %.2f)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                res[[as.character(s)]]$value))
  }
  saveRDS(res, f); res
}

# --------------------------- 5.3: solo clinicas ----------------------------
FV <- ajustar("sicure.v  (9 variables clinicas)", function()
                sicure.v(x_cov = clin, time = tt, delta = dd, randomsearch = TRUE),
              "sicure_v.rds")

# ------------------------ 5.4: solo senal funcional ------------------------
FF <- ajustar("sicure.f  (curva ECG)", function()
                sicure.f(x_cov = curvas, time = tt, delta = dd,
                         propvar = PROPVAR, randomsearch = TRUE),
              "sicure_f.rds")

# --------------------- 5.5: todo (ya ajustado antes) -----------------------
FVF <- if (file.exists(file.path(ROOT, "cache", "sicure.rds")))
         readRDS(file.path(ROOT, "cache", "sicure.rds")) else NULL

# ------------------------------- Resumen -----------------------------------
resumen <- function(nombre, fits, k_vec) {
  if (is.null(fits)) return(NULL)
  val <- sapply(fits, function(f) f$value)
  si  <- sapply(fits, function(f) as.numeric(f$si))
  d   <- length(fits[[1]]$par) - 4L
  cat(sprintf("\n--- %s ---\n  ajustes: %d | d = %d coeficientes (par = %d)\n",
              nombre, length(fits), d, length(fits[[1]]$par)))
  cat("  valor objetivo:", paste(sprintf("%.2f", val), collapse = ", "), "\n")
  if (ncol(si) > 1) {
    cm <- cor(si); cat(sprintf("  correlacion entre indices: min = %.3f\n",
                               min(cm[upper.tri(cm)])))
    th <- sapply(fits, function(f) c(1, f$par[seq_len(d)]))
    cmt <- cor(th); cat(sprintf("  correlacion entre vectores de coeficientes: min = %.3f\n",
                                min(cmt[upper.tri(cmt)])))
  }
  # correlacion del indice con las clinicas: la lectura que usa la memoria
  if (!is.null(k_vec)) {
    r <- apply(clin, 2, function(x) cor(x, si[, 1]))
    o <- order(-abs(r))[1:3]
    cat("  indice (arranque 1) vs clinicas, 3 mayores:",
        paste(sprintf("%s = %.2f", VARS_CLIN[o], r[o]), collapse = "; "), "\n")
  }
  list(valor = val, d = d)
}
R1 <- resumen("5.3  sicure.v  (clinicas)",  FV,  TRUE)
R2 <- resumen("5.4  sicure.f  (funcional)", FF,  TRUE)
R3 <- resumen("5.5  sicure.vf (todo)",      FVF, TRUE)

filas <- c()
add <- function(sec, mod, fits, R) {
  if (is.null(R)) return(invisible())
  val <- sprintf("%.1f", range(R$valor))
  filas <<- c(filas, sprintf("%s & \\texttt{%s} & %d & %s & %d \\\\",
    sec, mod, R$d, if (identical(val[1], val[2])) val[1] else paste(val, collapse = "--"),
    length(fits)))
}
add("5.3", "sicure.v",  FV,  R1)
add("5.4", "sicure.f",  FF,  R2)
add("5.5", "sicure.vf", FVF, R3)
if (length(filas))
  escribir_tabla_tex("tab_mod_sicure.tex", "llccc",
    "\\textbf{Sección} & \\textbf{Función} & \\textbf{Coef. del índice} & \\textbf{Función objetivo} & \\textbf{Arranques}",
    filas,
    sprintf("Modelos de curación \\textit{single-index} ajustados con el paquete \\textit{sicure} ($n=%d$, %d eventos). La columna de la función objetivo recoge el rango alcanzado entre los distintos puntos de arranque.",
            length(tt), sum(dd)),
    "tab:mod-sicure")

cat("\nTODO OK -> cache/sicure_v.rds, cache/sicure_f.rds\n")
