# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# representatividad.R - Compara las submuestras del estudio, para sustituir por
# evidencia la afirmacion del Capitulo 3 sobre la ausencia de sesgo:
#
#   (A) 194 curvas conservadas frente a 76 descartadas
#   (B) 270 pacientes con imagen frente a 261 sin imagen
#
# U de Mann-Whitney en las continuas y prueba exacta de Fisher en las binarias.
# Son unos treinta contrastes sin correccion por multiplicidad: lo informativo
# es el patron de conjunto, no cada p-valor por separado.
#
# Uso:  Rscript R/representatividad.R
# Salidas: salidas/tablas/tab_repr_{curvas,imagen}.tex
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))

# Variables binarias con prevalencia suficiente para que el contraste diga
# algo. Las que se registran en una o dos pacientes se omiten a proposito.
VARS_BIN <- c(AC = "Antraciclinas", antiHER2 = "Terapia anti-HER2",
              HTA = "Hipertensión", DL = "Dislipidemia",
              DM = "Diabetes mellitus", smoker = "Fumadora",
              exsmoker = "Exfumadora", ACprev = "Antraciclinas previas",
              RTprev = "Radioterapia previa")

# --------------------------- 1. Los tres grupos ----------------------------

cargar_malas <- function() {
  f <- file.path(ROOT, "RESULTADOS_ECG", "senales_malas_full.csv")
  stopifnot(file.exists(f))
  d <- as.data.frame(readr::read_csv(f, show_col_types = FALSE))
  stopifnot("num" %in% names(d))
  d
}

clin  <- cargar_clinicas()
buena <- cargar_ecg()
mala  <- cargar_malas()

n_ok  <- sort(unique(buena$num))          # filas de clin con curva valida
n_ko  <- sort(unique(mala$num))           # filas de clin con imagen descartada
n_img <- sort(union(n_ok, n_ko))          # filas de clin con imagen

cat(sprintf("curvas conservadas: %d | descartadas: %d | con imagen: %d | sin imagen: %d\n",
            length(n_ok), length(n_ko), length(n_img), nrow(clin) - length(n_img)))
stopifnot(length(intersect(n_ok, n_ko)) == 0,
          all(n_img %in% seq_len(nrow(clin))))

# ------------------------- 2. Motor de comparacion -------------------------

# Devuelve las tres celdas de una fila: resumen de cada grupo y p-valor. Las
# continuas se resumen con mediana [Q1-Q3], que es lo que corresponde a un
# contraste de rangos; las binarias, con la prevalencia.
comparar <- function(x, g, binaria = FALSE) {
  a <- x[g];  a <- a[!is.na(a)]
  b <- x[!g]; b <- b[!is.na(b)]
  if (length(a) < 3 || length(b) < 3) return(c("--", "--", "--"))
  if (binaria) {
    tab <- matrix(c(sum(a == 1), sum(a == 0), sum(b == 1), sum(b == 0)),
                  nrow = 2)
    p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    c(sprintf("%.1f\\%% (%d/%d)", 100 * mean(a == 1), sum(a == 1), length(a)),
      sprintf("%.1f\\%% (%d/%d)", 100 * mean(b == 1), sum(b == 1), length(b)),
      if (is.na(p)) "--" else fmt_p(p))
  } else {
    p <- tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                  error = function(e) NA_real_)
    q <- function(v) sprintf("%.2f [%.2f--%.2f]", median(v),
                             quantile(v, 0.25), quantile(v, 0.75))
    c(q(a), q(b), if (is.na(p)) "--" else fmt_p(p))
  }
}

# Construye y escribe la tabla completa para una particion de `d`.
#   d : data.frame clinico sobre el que se compara
#   g : vector logico de longitud nrow(d); TRUE = primer grupo
tabla_comparacion <- function(d, g, etq_a, etq_b, fichero, caption, label) {
  stopifnot(length(g) == nrow(d), !anyNA(g))
  fila <- function(etiqueta, x, binaria = FALSE)
    sprintf("%s & %s \\\\", etiqueta,
            paste(comparar(x, g, binaria = binaria), collapse = " & "))

  filas <- vapply(VARS_CLIN, function(v) fila(ETIQ[[v]], d[[v]]), character(1))

  # Respuesta y seguimiento: un desequilibrio aqui afectaria directamente a la
  # fraccion de curacion estimada.
  filas <- c(filas, "\\hline",
             fila("Cardiotoxicidad", d$CTRCD, binaria = TRUE),
             fila("Tiempo de seguimiento (días)", d$time),
             "\\hline")

  presentes <- intersect(names(VARS_BIN), names(d))
  filas <- c(filas,
             vapply(presentes, function(v) fila(VARS_BIN[[v]], d[[v]], TRUE),
                    character(1)))

  escribir_tabla_tex(
    fichero, "lccc",
    sprintf("\\textbf{Variable} & \\textbf{%s} & \\textbf{%s} & \\textbf{$p$}",
            etq_a, etq_b),
    unname(filas), caption, label,
    nota = paste("Variables continuas: mediana [Q1--Q3] y prueba U de",
                 "Mann--Whitney. Variables binarias: prevalencia (casos sobre",
                 "observaciones válidas) y prueba exacta de Fisher. Se",
                 "realizan del orden de treinta contrastes sin corrección por",
                 "multiplicidad, de modo que conviene atender al patrón de",
                 "conjunto antes que a cada $p$-valor por separado."))
}

# ------------------------------ 3. Las dos tablas --------------------------

# (A) Entre las 270 con imagen: curva conservada frente a imagen descartada.
d_img    <- clin[n_img, , drop = FALSE]
es_buena <- n_img %in% n_ok
tabla_comparacion(
  d_img, es_buena,
  sprintf("Conservadas ($n=%d$)", sum(es_buena)),
  sprintf("Descartadas ($n=%d$)", sum(!es_buena)),
  "tab_repr_curvas.tex",
  sprintf("Comparación de las pacientes cuya curva ECG se extrajo con éxito ($n=%d$) frente a aquellas cuya imagen hubo que descartar ($n=%d$), sobre el total de %d pacientes con imagen disponible.",
          sum(es_buena), sum(!es_buena), nrow(d_img)),
  "tab:repr-curvas")

# (B) Sobre las 531: con imagen frente a sin imagen.
tiene_img <- seq_len(nrow(clin)) %in% n_img
tabla_comparacion(
  clin, tiene_img,
  sprintf("Con imagen ($n=%d$)", sum(tiene_img)),
  sprintf("Sin imagen ($n=%d$)", sum(!tiene_img)),
  "tab_repr_imagen.tex",
  sprintf("Comparación de las pacientes con imagen de ecocardiograma disponible ($n=%d$) frente a las que no la tienen ($n=%d$), sobre el total de la cohorte ($n=%d$).",
          sum(tiene_img), sum(!tiene_img), nrow(clin)),
  "tab:repr-imagen")

cat("\nListo. Sube a Overleaf tab_repr_curvas.tex y tab_repr_imagen.tex, y\n",
    "sustituye en el Capítulo 3 la frase sobre los 'estudios preliminares'\n",
    "por lo que digan estas dos tablas.\n")
