# ============================================================
# Script  : 01_data_cleaning.R
# Author  : Tchandikou Ouadja Fare
# Task    : Load colon data, apply LOCF, export clean dataset
#           and create survival object
# Covers  : Lines 32-148 of Group_2_SA1.R
# ============================================================

library(survival)
library(dplyr)

# ── Load data ────────────────────────────────────────────────
colon <- survival::colon

colon_raw <- colon %>%
  filter(etype == 2) %>%
  distinct(id, .keep_all = TRUE)

cat("=== Dimensions (death records, one row per patient) ===\n")
print(dim(colon_raw))

cat("\n=== Missing Values Before Imputation ===\n")
mv_before <- colSums(is.na(colon_raw))
print(mv_before[mv_before > 0])

# ── LOCF imputation (nodes variable, <2% missing) ────────────
locf <- function(x) {
  for (i in seq_along(x))
    if (is.na(x[i]) && i > 1) x[i] <- x[i - 1]
  for (i in rev(seq_along(x)))
    if (is.na(x[i]) && i < length(x)) x[i] <- x[i + 1]
  return(x)
}

colon_clean <- colon_raw %>%
  arrange(id) %>%
  mutate(nodes = locf(nodes))

cat("\n=== Missing Values After LOCF ===\n")
print(colSums(is.na(colon_clean[, c("time","status","age","sex","nodes")])))

# ── Save clean dataset ────────────────────────────────────────
write.csv(colon_clean, "data/colon_clean.csv", row.names = FALSE)
cat("\nClean dataset saved to: data/colon_clean.csv\n")

# ── Survival object ───────────────────────────────────────────
SurvObj_colon <- Surv(time = colon_clean$time, event = colon_clean$status)
cat("Survival object created.\n")
