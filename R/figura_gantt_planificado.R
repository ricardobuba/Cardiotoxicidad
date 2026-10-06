# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# figura_gantt_planificado.R - Correccion de la tutora:
#   "Incluye un diagrama planificado, con intervalos que sean razonables,
#    segun lo que habiamos hablado de pasos a seguir a inicio del proyecto."
#
# La planificacion NO se inventa: son las cuatro fases del ANTEPROYECTO
# (abril de 2026), mas la indicacion de la tutora de que la redaccion de la
# memoria avanzase en paralelo a los analisis y no al final.
#
# Uso:  Rscript R/figura_gantt_planificado.R
# Salida: salidas/figuras/gantt_planificado.png
# ---------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- getwd()
source(file.path(ROOT, "R", "tfg_comun.R"))

tareas <- data.frame(
  etiqueta = c(
    "F1. Revisión bibliográfica",
    "F2. Digitalización de la señal",
    "F3. Modelización y evaluación",
    "F4. Interpretación y redacción",
    "Redacción (en paralelo)"),
  inicio = as.Date(c("2026-03-02", "2026-04-06", "2026-05-18", "2026-07-13", "2026-03-16")),
  fin    = as.Date(c("2026-04-17", "2026-05-29", "2026-07-31", "2026-09-11", "2026-09-11")),
  grupo  = c(1, 1, 1, 1, 2),
  stringsAsFactors = FALSE)

# Hitos acordados con la direccion
hitos <- data.frame(
  etiqueta = c("Anteproyecto", "Entrega"),
  fecha    = as.Date(c("2026-04-21", "2026-09-11")),
  stringsAsFactors = FALSE)

n   <- nrow(tareas)
y   <- rev(seq_len(n))
COL <- ifelse(tareas$grupo == 1, PINK, BLUE)
xlim <- range(c(tareas$inicio, tareas$fin))
meses <- seq(as.Date("2026-03-01"), as.Date("2026-10-01"), by = "month")

figura("gantt_planificado.png", 10, 3.9, {
  # El margen izquierdo se calcula a partir del ancho real de las etiquetas,
  # para que no salte "figure margins too large" si el tamano de letra por
  # defecto difiere entre instalaciones de R.
  par(mgp = c(2, 0.6, 0), xaxs = "i")
  anchos <- strwidth(tareas$etiqueta, units = "inches", cex = 0.85)
  izq <- min(11, max(6, ceiling(max(anchos) / par("csi")) + 1.5))
  par(mar = c(3.0, izq, 1.4, 1.0))
  plot(NA, xlim = xlim, ylim = c(0.4, n + 0.6), axes = FALSE, xlab = "", ylab = "")
  abline(v = meses, col = "grey88", lty = 1)
  for (i in seq_len(n))
    rect(tareas$inicio[i], y[i] - 0.28, tareas$fin[i], y[i] + 0.28,
         col = COL[i], border = NA)
  for (h in seq_len(nrow(hitos))) {
    abline(v = hitos$fecha[h], col = "grey30", lty = 2)
    text(hitos$fecha[h], n + 0.55, hitos$etiqueta[h], cex = 0.72,
         col = "grey25", adj = c(1.02, 0.5))
  }
  axis.Date(1, at = meses, format = "%b", cex.axis = 0.8, tick = FALSE)
  axis(2, at = y, labels = tareas$etiqueta, las = 1, tick = FALSE,
       cex.axis = 0.85, line = -0.4)
  box(col = "grey75")
})

cat("\nSube a Overleaf (imaxes/gestion):\n  ",
    file.path(SAL, "figuras", "gantt_planificado.png"), "\n")
