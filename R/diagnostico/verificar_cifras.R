# ---------------------------------------------------------------------------
# verificar_cifras.R  --  comprobacion rapida (segundos) de las cifras
# discutidas de la memoria. NO ajusta modelos: solo recalcula numeros.
#
# Uso:  Rscript R/diagnostico/verificar_cifras.R
# Requiere: readr, survival. Opcionalmente npcure (para el test de Maller-Zhou).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(readr); library(survival) })

if (!exists("ROOT")) ROOT <- getwd()

clin <- as.data.frame(read_delim(file.path(ROOT, "datos", "BC_cardiotox_clinical_variables.csv"),
        delim = ";", locale = locale(decimal_mark = ","), show_col_types = FALSE))
ecg  <- as.data.frame(read_csv(file.path(ROOT, "RESULTADOS_ECG", "senales_ecg_full.csv"),
        show_col_types = FALSE))

linea <- function(x) cat(strrep("-", 70), "\n", x, "\n", strrep("-", 70), "\n", sep="")

# ------------------------------------------------------------------ 1
linea(" 1. COHORTE COMPLETA (memoria: 531 pacientes, 54 eventos, 10,2%)")
cat("  n =", nrow(clin), "| eventos =", sum(clin$CTRCD),
    sprintf("(%.1f%%)\n", 100*mean(clin$CTRCD)))
cat("  mediana de time =", median(clin$time), "| maximo =", max(clin$time), "\n")

# ------------------------------------------------------------------ 2
linea(" 2. TEST DE MALLER-ZHOU  (memoria: p = 3,9e-14 -- CONFIRMADO)")
n  <- nrow(clin); Y <- clin$time; d <- clin$CTRCD
tF <- max(Y[d == 1]); tM <- max(Y)
lo <- 2*tF - tM
N  <- sum(Y > lo & Y <= tF)   # TODAS las obs del intervalo, no solo eventos
p_manual <- (1 - N/n)^n
cat("  mayor tiempo de EVENTO   T^F_max =", tF, "dias\n")
cat("  mayor tiempo OBSERVADO   T_max   =", tM, "dias\n")
cat("  intervalo del contraste  (", lo, ",", tF, "]\n")
cat("  observaciones en el intervalo  N =", N,
      " (de ellas eventos:", sum(Y > lo & Y <= tF & d == 1), ")\n")
cat(sprintf("  p (formula Maller-Zhou 1992, p.736) = %.4g\n", p_manual))

if (requireNamespace("npcure", quietly = TRUE)) {
  mz <- npcure::testmz(time, CTRCD, clin)
  cat(sprintf("  p (npcure::testmz)                  = %.6g\n", mz$pvalue))
  cat("  estadistico =", format(mz$aux$statistic),
      "| delta =", format(mz$aux$delta), "\n")
  cat("  intervalo npcure = [", format(mz$aux$interval[1]), ",",
      format(mz$aux$interval[2]), "]\n")
  if (abs(mz$pvalue - p_manual) > 1e-6)
    cat("  >> AVISO: npcure y la formula NO coinciden. Copiame esta seccion.\n")
  else
    cat("  >> Coinciden. La cifra de la memoria es correcta.\n")
} else {
  cat("  [npcure no instalado]  install.packages('npcure') para confirmarlo.\n")
}

# ------------------------------------------------------------------ 3
linea(" 3. EXTRACCION DE LA CURVA (memoria: 194 validas / 76 descartadas)")
cat("  curvas validas =", nrow(ecg),
    "| sin CTRCD =", sum(ecg$CTRCD == 0), "| con CTRCD =", sum(ecg$CTRCD == 1), "\n")

# ------------------------------------------------------------------ 4
linea(" 4. MESETA DE KAPLAN-MEIER (memoria, cap. modelizacion: 0,607)")
sub <- clin[ecg$num, ]
km  <- survfit(Surv(time, CTRCD) ~ 1, data = sub)
cat("  subconjunto con imagen: n =", nrow(sub), "| eventos =", sum(sub$CTRCD), "\n")
cat(sprintf("  meseta S(t_max) = %.4f\n", min(km$surv)))
ok <- complete.cases(sub[, c("heart_rate","age","weight","height",
                             "LVEF","PWT","LAd","LVDd","LVSd")])
cat("  tras exigir clinicas completas: n =", sum(ok), "| eventos =", sum(sub$CTRCD[ok]), "\n")
cat("  (la memoria usa n=181 y 25 eventos tras excluir tambien 11 curvas atipicas)\n")

# ------------------------------------------------------------------ 5
linea(" 5. BIVARIANTE: ritmo cardiaco (memoria: 80,3 vs 74,0 lpm; p=0,001)")
hr <- clin$heart_rate
cat(sprintf("  media CTRCD=1: %.1f | CTRCD=0: %.1f\n",
            mean(hr[clin$CTRCD==1], na.rm=TRUE), mean(hr[clin$CTRCD==0], na.rm=TRUE)))
cat(sprintf("  Mann-Whitney p = %.4f\n",
            wilcox.test(hr ~ clin$CTRCD)$p.value))

cat("\nFIN. Copiame toda esta salida.\n")
