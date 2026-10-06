# Cuenta cuantas replicas bootstrap se descartaron por divergencia en los
# modelos NUEVOS (M1, M2, M3). Sirve para actualizar el Anexo del bootstrap,
# cuyas cifras ("una en M1, ninguna en M2 y dos en M3") corresponden a los
# modelos anteriores. Tarda un segundo: solo lee los .rds.
if (!exists("ROOT")) ROOT <- getwd()
A1 <- readRDS(file.path(ROOT, "cache", "backward_clinicas.rds"))
V2 <- readRDS(file.path(ROOT, "cache", "modelos_v2.rds"))
f <- list(M1 = A1$fit, M2 = V2$M2$fit, M3 = V2$M3$fit)
for (m in names(f)) {
  x <- f[[m]]
  cat(sprintf("%s: %d replicas validas, %d descartadas por divergencia\n",
              m, ifelse(is.null(x$nboot_validas), NA, x$nboot_validas), ifelse(is.null(x$nboot_descartadas), NA, x$nboot_descartadas)))
}
