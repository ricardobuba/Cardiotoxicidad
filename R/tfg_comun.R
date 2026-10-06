# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# tfg_comun.R - Utilidades compartidas por los scripts de R del TFG.
#
# Antes, el bloque de carga de datos, deteccion de atipicos funcionales y FPCA
# estaba DUPLICADO en eda_funcional_superv.R y en aplicacion_modelos_curacion.R.
# Aqui se centraliza (equivalente en R de lo que config.py hace en Python), de
# modo que ambos scripts parten exactamente del mismo subconjunto de curvas y
# de los mismos scores: si se cambia el factor del fbplot o NHARM, cambia en
# los dos a la vez y no se pueden desincronizar.
#
# Uso:  source(file.path(ROOT, "R", "tfg_comun.R"))  (desde la raiz del repo)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(survival)
  library(fda)
})

# ---------------------------- 1. Rutas y estilo ----------------------------
if (!exists("ROOT")) ROOT <- getwd()
SAL <- file.path(ROOT, "salidas")
DATOS <- file.path(ROOT, "datos")   # CSV originales del CHUAC
FIG <- file.path(SAL, "figuras")
TAB <- file.path(SAL, "tablas")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

PINK <- "#C5006E"   # CTRCD (evento)
BLUE <- "#4C9BD4"   # No CTRCD (censurado)
GRAY <- "#8C8C8C"
PAL  <- c("0" = BLUE, "1" = PINK)

# Constantes del analisis (unicas para todos los scripts de R)
NHARM     <- 8      # componentes principales funcionales
FACTOR_FB <- 1.5    # factor del functional boxplot (analogo al 1.5*IQR)
# 500, y no los 100 que smcure trae por defecto: con B = 100 el error de Monte
# Carlo del p-valor es de +-0,02, demasiado para los resultados de las FPC, que
# caen en p ~ 0,06. Solo afecta a los 12 ajustes con Var = TRUE.
NBOOT     <- 500    # remuestreos bootstrap de smcure
# Umbral para descartar replicas bootstrap divergentes (ver ajustar_curacion).
MAX_COEF_BOOT <- 10
SEED      <- 1      # semilla del bootstrap de smcure
ALPHA     <- 0.05

# Variables clinicas usadas en la modelizacion: las mismas que selecciona el
# script de la directora ("Codigo_aplicacion_smcure_sicure.R").
VARS_CLIN <- c("heart_rate", "age", "weight", "height",
               "LVEF", "PWT", "LAd", "LVDd", "LVSd")

# Etiquetas legibles para las tablas de la memoria
# OJO: heart_rate es la FRECUENCIA cardiaca; el ritmo cardiaco es heart_rhythm,
# otra variable del mismo fichero.
ETIQ <- c(heart_rate = "Frecuencia cardíaca", age = "Edad", weight = "Peso",
          height = "Estatura", LVEF = "FEVI", PWT = "PWT", LAd = "LAd",
          LVDd = "LVDd", LVSd = "LVSd", FPC1 = "FPC1", FPC2 = "FPC2",
          FPC3 = "FPC3", "(Intercept)" = "Intercepto")

# ------------------------- 2. Helpers de figuras ---------------------------
openpdf <- function(name, w, h) pdf(file.path(FIG, name), width = w, height = h)

# Abre un dispositivo grafico eligiendo por la extension del nombre: .png usa
# el dispositivo raster (a 300 ppp, con w y h en pulgadas, igual que pdf()) y
# .pdf el vectorial. Se usa PNG en las figuras del capitulo de modelizacion
# porque el recorte del PDF resultaba irregular al incrustarlo en LaTeX.
FIG_RES <- 300
openfig <- function(name, w, h) {
  if (grepl("\\.png$", name, ignore.case = TRUE)) {
    # type="cairo" da mejor antialiasing, pero no esta en todas las
    # instalaciones; si no lo hay, se usa el dispositivo por defecto.
    args <- list(file.path(FIG, name), width = w * FIG_RES, height = h * FIG_RES,
                 res = FIG_RES, bg = "white")
    if (isTRUE(capabilities("cairo"))) args$type <- "cairo"
    do.call(png, args)
  } else {
    pdf(file.path(FIG, name), width = w, height = h)
  }
}
done    <- function(name) { dev.off(); cat("  [fig]", name, "\n") }

# Dibuja una figura cerrando SIEMPRE el dispositivo, tambien si el codigo de
# dentro falla. Sin esto, un error a mitad de un plot deja el dispositivo
# abierto y la figura siguiente se escribe en el fichero equivocado.
figura <- function(name, w, h, expr) {
  openfig(name, w, h)
  on.exit(done(name), add = TRUE)
  force(expr)
  invisible(name)
}

# ------------------------- 3. Helper de tablas LaTeX -----------------------
# Mismo estilo (rowcolors udcgray/udcpink + resizebox) que tab_fpca.tex y que
# las tablas que genera eda_analysis.py, para que la memoria sea homogenea.
# Coma decimal en las tablas, para que concuerden con la prosa de la memoria.
# Se emplea "{,}" y no una coma suelta porque en modo matematico esta ultima se
# compone como puntuacion. El patron solo actua entre digitos.
coma_decimal <- function(x) gsub("(?<=[0-9])\\.(?=[0-9])", "{,}", x, perl = TRUE)

escribir_tabla_tex <- function(fichero, colspec, cabecera, filas, caption, label,
                               nota = NULL, coma = TRUE) {
  if (isTRUE(coma)) {
    cabecera <- coma_decimal(cabecera)
    filas    <- coma_decimal(filas)
    caption  <- coma_decimal(caption)
    if (!is.null(nota)) nota <- coma_decimal(nota)
  }
  tex <- paste0(
    "\\begin{table}[htbp]\n  \\centering\n",
    "  \\rowcolors{2}{white}{udcgray!25}\n",
    "  \\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%\n",
    "  \\begin{tabular}{", colspec, "}\n",
    "    \\rowcolor{udcpink!25}\n    ", cabecera, " \\\\ \\hline\n",
    paste0("    ", filas, collapse = "\n"), "\n",
    "  \\end{tabular}}\n",
    if (is.null(nota)) "" else paste0("  \\\\[2pt]\n  {\\footnotesize ", nota, "}\n"),
    "  \\caption{", caption, "}\n  \\label{", label, "}\n\\end{table}\n")
  writeLines(tex, file.path(TAB, fichero), useBytes = TRUE)
  cat("  [tab]", fichero, "\n")
}

# Formatea un p-valor al estilo de la memoria (coma decimal la pone LaTeX).
# Las dos ramas van en modo matematico: si una fuera texto plano, en el PDF se
# compondrian con fuentes distintas dentro de la misma columna.
fmt_p <- function(p) ifelse(p < 0.001, "$<0{,}001$",
                            sprintf("$%s$", sub("\\.", "{,}", sprintf("%.3f", p))))
fmt_n <- function(x, d = 3) formatC(x, format = "f", digits = d)

# --------------------------- 4. Carga de datos -----------------------------
cargar_clinicas <- function() {
  f <- file.path(DATOS, "BC_cardiotox_clinical_variables.csv")
  stopifnot(file.exists(f))
  d <- as.data.frame(read_delim(f, delim = ";",
                                locale = locale(decimal_mark = ","),
                                show_col_types = FALSE))
  stopifnot(nrow(d) == 531, all(c("time", "CTRCD") %in% names(d)))
  d
}

cargar_ecg <- function() {
  f <- file.path(ROOT, "RESULTADOS_ECG", "senales_ecg_full.csv")
  stopifnot(file.exists(f))
  d <- as.data.frame(read_csv(f, show_col_types = FALSE))
  stopifnot(all(c("num", "CTRCD", "time") %in% names(d)))
  d
}

matriz_curvas <- function(ecg) {
  pcols <- grep("^p[0-9]+$", names(ecg), value = TRUE)
  stopifnot(length(pcols) == 1001)
  as.matrix(ecg[, pcols])
}

# --------------- 5. Atipicos funcionales (functional boxplot) --------------
# Modified Band Depth (Sun & Genton, 2011): analogo funcional de la regla
# 1.5*IQR usada con las variables escalares.
atipicos_funcionales <- function(E, factor = FACTOR_FB) {
  n <- nrow(E)
  ranks   <- apply(E, 2, rank)
  depth_t <- (ranks - 1) * (n - ranks)
  MBD     <- rowMeans(depth_t) / (n * (n - 1) / 2)
  central <- order(MBD, decreasing = TRUE)[1:ceiling(0.5 * n)]
  env_min <- apply(E[central, , drop = FALSE], 2, min)
  env_max <- apply(E[central, , drop = FALSE], 2, max)
  w  <- env_max - env_min
  lo <- env_min - factor * w
  hi <- env_max + factor * w
  apply(E, 1, function(v) any(v < lo | v > hi))
}

# ------------------------------ 6. FPCA ------------------------------------
# Mismo enfoque que el script del Departamento de Estadistica de la UDC:
# objeto funcional con fd() y pca.fd(..., nharm = 8).
fpca_ecg <- function(E, nharm = NHARM) {
  stopifnot(requireNamespace("fda", quietly = TRUE))
  fda::pca.fd(fda::fd(t(E)), nharm = nharm)
}

# ------------- 7. Preparacion completa para la modelizacion ----------------
# Devuelve el conjunto emparejado clinicas + scores FPCA, ya estandarizado.
# "num" es la fila (1-indexada) del paciente en el CSV clinico, que guarda
# extraer_ecg.py; se comprueba que CTRCD y time coinciden en ambos ficheros.
preparar_datos <- function(verbose = TRUE) {
  clin <- cargar_clinicas()
  ecg  <- cargar_ecg()
  E_all <- matriz_curvas(ecg)

  out <- atipicos_funcionales(E_all)
  E      <- E_all[!out, , drop = FALSE]
  ecg_ok <- ecg[!out, ]
  if (verbose) cat(sprintf("  curvas: %d validas, %d atipicas -> %d\n",
                           nrow(E_all), sum(out), nrow(E)))

  clin_fun <- clin[ecg_ok$num, ]
  stopifnot(all(clin_fun$CTRCD == ecg_ok$CTRCD),
            all(clin_fun$time  == ecg_ok$time))

  pca <- fpca_ecg(E)
  scores <- as.data.frame(scale(pca$scores))
  names(scores) <- paste0("FPC", seq_len(ncol(scores)))

  ok <- complete.cases(clin_fun[, VARS_CLIN])
  if (verbose) cat(sprintf("  pacientes: %d con ECG (%d eventos); %d completos en clinicas (%d eventos)\n",
                           nrow(clin_fun), sum(clin_fun$CTRCD),
                           sum(ok), sum(clin_fun$CTRCD[ok])))

  list(clin = clin, ecg = ecg, clin_fun = clin_fun, ecg_ok = ecg_ok,
       E_all = E_all, E = E, atipicos = out, pca = pca,
       # conjunto de modelizacion (casos completos, estandarizado)
       datos = cbind(as.data.frame(scale(clin_fun[ok, VARS_CLIN])), scores[ok, ]),
       time = clin_fun$time[ok], delta = clin_fun$CTRCD[ok],
       ok = ok)
}

# ------------- 8. Contraste de Maller-Zhou (seguimiento suficiente) --------
# H0: tau_{F0} > tau_G (seguimiento INsuficiente) frente a H1: tau_{F0} <= tau_G.
# Un p-valor pequeno permite rechazar H0 y concluir que el seguimiento es
# suficiente. Se usa npcure::testmz si esta disponible; si no, la formula
# original de Maller y Zhou (1992, p. 736), que es la que implementa npcure.
maller_zhou <- function(t, d) {
  n <- length(t)
  tmax  <- max(t)
  tFmax <- max(t[d == 1])
  dn    <- tmax - tFmax
  lo    <- tFmax - dn
  # OJO: N cuenta TODAS las observaciones del intervalo (censuradas incluidas),
  # no solo los eventos. Es lo que hace npcure::testmz y lo que corresponde al
  # contraste: lo que informa sobre si el seguimiento es suficiente es cuanta
  # masa observada hay junto al mayor tiempo de evento, no cuantos eventos.
  N     <- sum(t > lo & t <= tFmax)
  list(statistic = N, delta = dn, interval = c(lo, tFmax),
       pvalue = (1 - N / n)^n)
}

# ---------------- 9. Ajuste del modelo de mixtura de curacion --------------
# Envoltorio unico sobre smcure::smcure, para que todos los modelos del
# capitulo se ajusten exactamente igual (misma semilla, mismo nboot, misma
# parametrizacion: incidencia y latencia con el mismo conjunto de covariables).
ajustar_curacion <- function(covs, datos, time, delta,
                             nboot = NBOOT, seed = SEED, Var = TRUE,
                             max_intentos = 4L) {
  d <- cbind(time = time, delta = delta, datos[, covs, drop = FALSE])
  f_lat  <- as.formula(paste("Surv(time, delta) ~", paste(covs, collapse = " + ")))
  f_cure <- as.formula(paste("~", paste(covs, collapse = " + ")))

  # smcure escribe "Program is running..be patient..." en CADA llamada; en el
  # bootstrap serian cientos de lineas, asi que se captura su salida.
  ajuste1 <- function(dd) {
    r <- NULL
    invisible(utils::capture.output(
      r <- smcure::smcure(f_lat, cureform = f_cure, data = dd,
                          model = "ph", Var = FALSE)))
    r
  }

  fit <- tryCatch(suppressWarnings(ajuste1(d)), error = function(e) e)
  if (inherits(fit, "error"))
    stop("el ajuste principal de smcure no converge con (",
         paste(covs, collapse = " + "), "): ", conditionMessage(fit))

  # OJO: la semilla solo se fija cuando hay bootstrap (Var = TRUE), que es el
  # unico paso aleatorio. Fijarla siempre rompe la validacion cruzada: cada
  # llamada reiniciaria el generador y todas las repeticiones acabarian usando
  # la MISMA particion en folds, con lo que la desviacion tipica entre
  # repeticiones saldria practicamente nula.
  if (!isTRUE(Var)) return(fit)

  # --------------------------------------------------------------------------
  # Bootstrap propio, equivalente al de smcure (remuestreo con reemplazamiento
  # estratificado por estado de evento; EE = desviacion tipica entre replicas;
  # p por aproximacion normal), pero TOLERANTE a replicas que no convergen.
  #
  # Motivo: cuando el EM de una replica devuelve NA, smcure aborta el ajuste
  # entero con "while (convergence > eps & i < emmax): valor ausente donde
  # TRUE/FALSE es necesario". Con 25 eventos eso ocurre. Aqui esa replica se
  # descarta y se sortea otra, hasta reunir `nboot` validas o agotar
  # max_intentos*nboot intentos. El recuento queda en fit$nboot_descartadas
  # para poder declararlo en la memoria.
  # --------------------------------------------------------------------------
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
    # Se descarta la replica si el EM ha divergido. Criterio: algun coeficiente
    # con |valor| > MAX_COEF_BOOT, senal inequivoca de separacion en la parte
    # logistica. El diagnostico de 30/08/2026 (diag_bootstrap.R) justifica el
    # umbral: en M1 la mayor replica sana tiene max|b| = 7,1 y la unica
    # divergente 639,9; en M3, 5,3 frente a 20,8 y 53,4. Entre 8 y 600 el
    # resultado es identico, de modo que el valor exacto del corte es
    # irrelevante. Umbrales mas estrictos (3 o 5) o los estimadores robustos
    # (MAD, IQR) empiezan a recortar la cola legitima del bootstrap y reducen
    # artificialmente los errores estandar: NO usarlos.
    if (is.null(r) || !all(is.finite(r$b)) || !all(is.finite(r$beta))) next
    if (max(abs(r$b)) > MAX_COEF_BOOT || max(abs(r$beta)) > MAX_COEF_BOOT) next
    ok <- ok + 1L
    Bb[ok, ] <- r$b; Bg[ok, ] <- r$beta
  }
  if (ok < 2L)
    stop("bootstrap imposible: solo ", ok, " replicas validas en ", intentos,
         " intentos para (", paste(covs, collapse = " + "), ")")
  if (ok < nboot)
    warning(sprintf("solo %d replicas bootstrap validas de las %d pedidas (%d intentos)",
                    ok, nboot, intentos))
  Bb <- Bb[seq_len(ok), , drop = FALSE]
  Bg <- Bg[seq_len(ok), , drop = FALSE]

  fit$b_var       <- apply(Bb, 2, var)
  fit$b_sd        <- sqrt(fit$b_var)
  fit$b_zvalue    <- fit$b / fit$b_sd
  fit$b_pvalue    <- (1 - pnorm(abs(fit$b_zvalue))) * 2
  fit$beta_var    <- apply(Bg, 2, var)
  fit$beta_sd     <- sqrt(fit$beta_var)
  fit$beta_zvalue <- fit$beta / fit$beta_sd
  fit$beta_pvalue <- (1 - pnorm(abs(fit$beta_zvalue))) * 2
  for (nm in c("b_var", "b_sd", "b_zvalue", "b_pvalue"))
    names(fit[[nm]]) <- names(fit$b)
  for (nm in c("beta_var", "beta_sd", "beta_zvalue", "beta_pvalue"))
    names(fit[[nm]]) <- names(fit$beta)
  fit$nboot_validas     <- ok
  fit$nboot_intentos    <- intentos
  fit$nboot_descartadas <- intentos - ok
  cat(sprintf("    [bootstrap] %d replicas validas de %d intentos (%d descartadas por divergencia)\n",
              ok, intentos, intentos - ok))
  fit
}


# Fraccion de curacion media estimada por un ajuste
fraccion_curacion <- function(fit, covs, datos) {
  Z  <- cbind(1, as.matrix(datos[, covs, drop = FALSE]))
  pi <- 1 / (1 + exp(-drop(Z %*% fit$b)))
  1 - mean(pi)
}

# Supervivencia basal de las susceptibles en t0, a partir del vector s que
# devuelve smcure (una estimacion por observacion, en el orden de los datos de
# entrenamiento).
#
# OJO: smcure NO devuelve s ordenado por tiempo. Hay que ordenar antes de
# evaluar. La version anterior llamaba a approx(..., ties = "ordered") sobre el
# vector sin ordenar, y ese argumento le dice justamente a approx que NO ordene
# (en regularize.values el bloque que ordena esta dentro de
# `if (!identical(ties, "ordered"))`). La busqueda binaria interna daba
# entonces un valor arbitrario, sin error ni aviso.
#
# S0 es una funcion escalonada, asi que se toma el ultimo valor observado con
# tiempo <= t0 en lugar de interpolar linealmente entre saltos.
s0_en <- function(time_train, s_train, t0) {
  stopifnot(length(time_train) == length(s_train))
  o <- order(time_train)
  tt <- time_train[o]; ss <- s_train[o]
  i <- sum(tt <= t0)
  if (i == 0L) 1 else ss[i]     # antes del primer tiempo observado, S0 = 1
}

# Supervivencia poblacional predicha  S(t|z) = (1-pi) + pi * S0(t)^exp(beta'x)
predecir_curacion <- function(fit, newdata, t0, time_train, s_train) {
  # El envoltorio ajustar_curacion() usa las mismas covariables, y en el mismo
  # orden, en incidencia y en latencia; b lleva ademas el intercepto delante.
  # Si eso dejara de cumplirse, Z y b se desalinearian en silencio.
  stopifnot(identical(fit$betanm, fit$bnm[-1]))
  s0 <- s0_en(time_train, s_train, t0)
  X  <- as.matrix(newdata[, fit$betanm, drop = FALSE])
  Z  <- cbind(1, X)
  pi <- 1 / (1 + exp(-drop(Z %*% fit$b)))
  (1 - pi) + pi * s0^exp(drop(X %*% fit$beta))
}

# --------- 10. Etiquetas de los modelos del capitulo de modelizacion -------
# Una sola definicion para las tres salidas (tabla, figura y consola). Antes
# estaban duplicadas en modelizacion_curacion.R, cindex_modelizacion.R y
# figuras_modelizacion.R, y las copias habian divergido: una con tildes y otra
# sin ellas, de modo que reejecutar el guion maestro degradaba la tabla de la
# memoria.
NOM_MODELOS <- c(Cox9  = "Cox (9 clínicas)",
                 CoxHR = "Cox (frecuencia cardíaca)",
                 M1    = "Curación M1 (clínicas)",
                 M2    = "Curación M2 (ECG)",
                 M3    = "Curación M3 (mixto)")
# Version sin tildes para los ejes de las figuras (el dispositivo grafico no
# siempre resuelve bien los acentos en Windows).
NOM_MODELOS_ASCII <- c(Cox9  = "Cox (9 clinicas)",
                       CoxHR = "Cox (frecuencia cardiaca)",
                       M1    = "Curacion M1 (clinicas)",
                       M2    = "Curacion M2 (ECG)",
                       M3    = "Curacion M3 (mixto)")

# ------------- 11. Validacion cruzada repetida (C-index) -------------------
# Devuelve una lista con los R valores del C-index de cada modelo. La FPCA se
# calcula UNA vez fuera (paso no supervisado); t0 se pasa desde fuera para que
# el horizonte sea el mismo en todas las repeticiones.
#
# Si smcure no converge en un fold concreto, ese fold se anota como fallido y
# la repeticion se descarta, en lugar de tumbar el bucle entero (que son
# R x K x 3 ajustes y tarda mucho).
validacion_cruzada <- function(datos, time, delta, covs, t0,
                               K = 5, R = 20, semilla = 20260810,
                               verbose = TRUE) {
  cindex <- function(score, t, d)
    as.numeric(concordance(Surv(t, d) ~ score)$concordance)
  nombres <- c("Cox9", "CoxHR", names(covs))
  res <- setNames(vector("list", length(nombres)), nombres)
  fallos <- 0L
  set.seed(semilla)
  for (r in seq_len(R)) {
    fold <- sample(rep(1:K, length.out = length(time)))
    sc <- setNames(lapply(nombres, function(z) numeric(length(time))), nombres)
    ok_rep <- TRUE
    for (k in 1:K) {
      tr <- fold != k; te <- !tr
      c9 <- coxph(as.formula(paste("Surv(time[tr], delta[tr]) ~",
                                   paste(VARS_CLIN, collapse = " + "))),
                  data = datos[tr, ])
      sc$Cox9[te]  <- -predict(c9, newdata = datos[te, ], type = "lp")
      c1 <- coxph(Surv(time[tr], delta[tr]) ~ heart_rate, data = datos[tr, ])
      sc$CoxHR[te] <- -predict(c1, newdata = datos[te, ], type = "lp")
      for (m in names(covs)) {
        f <- tryCatch(
          suppressMessages(ajustar_curacion(covs[[m]], datos[tr, ],
                                            time[tr], delta[tr], Var = FALSE)),
          error = function(e) NULL)
        if (is.null(f)) { fallos <- fallos + 1L; ok_rep <- FALSE; next }
        sc[[m]][te] <- predecir_curacion(f, datos[te, , drop = FALSE], t0,
                                         time[tr], f$s)
      }
    }
    if (!ok_rep) next
    for (nm in nombres) res[[nm]] <- c(res[[nm]], cindex(sc[[nm]], time, delta))
    if (verbose) cat("  repeticion", r, "de", R, "\n")
  }
  if (fallos > 0)
    cat("  [AVISO]", fallos, "ajustes no convergieron; sus repeticiones se descartan.\n")
  # Comprobacion de que la semilla no se esta fijando dentro del bucle: si se
  # fijara, todas las repeticiones usarian la MISMA particion y la desviacion
  # tipica saldria practicamente nula.
  nd <- sapply(res, function(x) length(unique(round(x, 8))))
  if (any(nd < length(res[[1]]) / 2))
    cat("  [AVISO] valores repetidos entre repeticiones: revisar la semilla.\n")
  res
}

# Tabla y figura del C-index, a partir de la salida de validacion_cruzada()
tabla_cindex <- function(res, n, eventos, t0, K = 5, R = 20) {
  cind <- data.frame(modelo = names(res),
                     C = sapply(res, mean), sd = sapply(res, sd))
  escribir_tabla_tex("tab_mod_cindex.tex", "lcc",
    "\\textbf{Modelo} & \\textbf{C-index} & \\textbf{DT}",
    sprintf("%s & %s & %s \\\\", NOM_MODELOS[cind$modelo],
            fmt_n(cind$C, 3), fmt_n(cind$sd, 3)),
    sprintf(paste("Índice de concordancia estimado por validación cruzada de %d particiones",
                  "repetida %d veces ($n=%d$, %d eventos). Para los modelos de curación la",
                  "puntuación de riesgo es la supervivencia poblacional predicha en",
                  "$t_0=%.0f$ días (la mediana del tiempo hasta el evento). La FPCA y el",
                  "horizonte $t_0$ se calculan una única vez sobre toda la muestra, por ser",
                  "pasos no supervisados. La desviación típica recoge la variabilidad entre",
                  "repeticiones, no la incertidumbre de la estimación."),
            K, R, n, eventos, t0),
    "tab:mod-cindex")
  cind
}

figura_cindex <- function(res) {
  figura("mod_cindex.png", 9, 5, {
    op <- par(mar = c(4.5, 11, 3, 1.5), mgp = c(2.6, 0.8, 0))
    boxplot(rev(res), names = rev(NOM_MODELOS_ASCII[names(res)]),
            horizontal = TRUE, las = 1,
            col = adjustcolor(BLUE, 0.6), border = GRAY,
            xlab = "C-index (validacion cruzada)",
            main = "Discriminacion de los modelos")
    abline(v = 0.5, col = PINK, lty = 2, lwd = 2)
    par(op)
  })
}

# ------------- 12. Figuras del modelo M3 (probabilidad y KM) ---------------
# Se generan en PNG (300 ppp) y no en PDF: el recorte del PDF salia irregular
# al incrustarlo con \includegraphics.
figuras_modelo_m3 <- function(fit, covs_m3, datos, time, delta) {
  pi3 <- 1 / (1 + exp(-drop(cbind(1, as.matrix(datos[, covs_m3])) %*% fit$b)))

  figura("mod_prob_curacion.png", 12, 4.5, {
    op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1.2), mgp = c(2.6, 0.8, 0))
    for (v in covs_m3) {
      plot(datos[[v]], 1 - pi3, pch = 19,
           col = adjustcolor(PAL[as.character(delta)], 0.7),
           xlab = paste(ETIQ[v], "(estandarizada)"),
           ylab = if (v == covs_m3[1]) "Prob. estimada de curacion" else "",
           main = ETIQ[v], ylim = c(0, 1))
      if (v == covs_m3[1])
        legend("topright", c("No CTRCD", "CTRCD"), pch = 19,
               col = c(BLUE, PINK), bty = "n")
    }
    par(op)
  })

  grp <- cut(pi3, quantile(pi3, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
             labels = c("Riesgo bajo", "Riesgo medio", "Riesgo alto"))
  lr   <- survdiff(Surv(time, delta) ~ grp)
  p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  figura("mod_superv_predicha.png", 9, 5.5, {
    op <- par(mar = c(4.5, 4.5, 3, 1.2), mgp = c(2.8, 0.8, 0))
    plot(survfit(Surv(time, delta) ~ grp), col = c(BLUE, GRAY, PINK), lwd = 2,
         ylim = c(0.4, 1), xlab = "Tiempo (dias)",
         ylab = "Prob. de no sufrir cardiotoxicidad",
         main = "Kaplan-Meier por tercil de riesgo predicho (modelo M3)")
    legend("bottomleft", levels(grp), col = c(BLUE, GRAY, PINK), lwd = 2, bty = "n")
    legend("topright", sprintf("log-rank p = %.4f", p_lr), bty = "n")
    par(op)
  })
  cat("  log-rank terciles de riesgo (M3): p =", signif(p_lr, 3), "\n")
  invisible(list(pi = pi3, p_logrank = p_lr))
}
