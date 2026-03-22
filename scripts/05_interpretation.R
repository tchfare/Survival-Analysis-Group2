# ============================================================
# Script  : 05_interpretation.R
# Author  : Josiane Kazanenda
# Task    : Extract MLE parameters (Table 6), print full
#           model interpretation summary
# Covers  : Lines 378-427 of Group_2_SA1.R
# ============================================================

library(survival)
library(flexsurv)
library(knitr)
library(kableExtra)

colon_clean   <- read.csv("data/colon_clean.csv")
SurvObj_colon <- Surv(time=colon_clean$time, event=colon_clean$status)

base_dir <- getwd()
tbl_dir  <- file.path(base_dir, "outputs", "tables")
dir.create(tbl_dir, recursive=TRUE, showWarnings=FALSE)

# ── Refit best model ─────────────────────────────────────────
best_model      <- flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="gengamma")
best_model_name <- "GenGamma"

cat("\n=== Full Model Summary ===\n")
print(best_model)

# ── Table 6: MLE parameters ───────────────────────────────────
param_raw <- as.data.frame(best_model$res)

tbl6_data <- data.frame(
  Parameter      = rownames(param_raw),
  Estimate       = round(param_raw$est,    4),
  `Std. Error`   = round(param_raw$se,     4),
  `95% CI Lower` = round(param_raw$`L95%`, 4),
  `95% CI Upper` = round(param_raw$`U95%`, 4),
  check.names    = FALSE
)

write.csv(tbl6_data, file.path(tbl_dir,"Table6_MLE_Parameters.csv"), row.names=FALSE)

save_kable(
  kable(tbl6_data,
        caption=paste("Table 6: MLE Parameter Estimates —", best_model_name, "Model"),
        align="lrrrr") %>%
    kable_styling(bootstrap_options=c("striped","hover","condensed","bordered"),
                  full_width=FALSE, position="center") %>%
    column_spec(1, bold=TRUE),
  file=file.path(tbl_dir,"Table6_MLE_Parameters.html")
)

# ── Load mean/variance from Table 5 ──────────────────────────
tbl5 <- read.csv(file.path(tbl_dir,"Table5_Mean_Variance.csv"))
mean_T <- as.numeric(tbl5$Value[1])
var_T  <- as.numeric(tbl5$Value[3])

# ── Interpretation summary ────────────────────────────────────
n_total    <- nrow(colon_clean)
n_events   <- sum(colon_clean$status == 1)
n_censored <- sum(colon_clean$status == 0)

cat("\n========================================================\n")
cat("INTERPRETATION SUMMARY\n")
cat("========================================================\n")
cat("Model     :", best_model_name, "\n")
cat("Dataset   : Colon Cancer\n")
cat("Patients  :", n_total, " | Events:", n_events,
    " | Censored:", n_censored, "\n\n")
cat("Parameters (Generalised Gamma):\n")
cat("  mu    = 7.545  → log-scale location; median ≈ exp(7.545) ≈ 1,900 days\n")
cat("  sigma = 1.575  → spread of the survival distribution\n")
cat("  Q     = -0.481 → Q < 0: hazard rises then falls (non-monotone)\n\n")
cat("Life Functions:\n")
cat("  S(t) : Probability of surviving beyond time t\n")
cat("  h(t) : Instantaneous risk — peaks near day 300 then declines\n")
cat("  f(t) : Density — most deaths concentrated around day 300\n")
cat("  F(t) : 50% of patients die within ~2,500 days (~6.8 years)\n\n")
cat(sprintf("Mean survival : %.2f days (%.1f years)\n", mean_T, mean_T/365))
cat(sprintf("Std deviation : %.2f days\n", sqrt(var_T)))
cat("========================================================\n")
