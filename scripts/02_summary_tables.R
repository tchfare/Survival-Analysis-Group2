# ============================================================
# Script  : 02_summary_tables.R
# Author  : Rabecca Kanini Kating'u
# Task    : Produce Table 1 (variable statistics) and
#           Table 2 (sample summary), save as CSV and HTML
# Covers  : Lines 49-143 of Group_2_SA1.R
# ============================================================

library(dplyr)
library(knitr)
library(kableExtra)

colon_raw   <- survival::colon %>% filter(etype == 2) %>% distinct(id, .keep_all = TRUE)
colon_clean <- read.csv("data/colon_clean.csv")

base_dir <- getwd()
tbl_dir  <- file.path(base_dir, "outputs", "tables")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ── Table 1: Variable summary ─────────────────────────────────
tbl1_data <- data.frame(
  Variable = c("Time (days)", "Status", "Age (years)", "Sex", "Nodes"),
  Min      = c(min(colon_raw$time,   na.rm=TRUE), min(colon_raw$status, na.rm=TRUE),
               min(colon_raw$age,    na.rm=TRUE), min(colon_raw$sex,    na.rm=TRUE),
               min(colon_raw$nodes,  na.rm=TRUE)),
  Max      = c(max(colon_raw$time,   na.rm=TRUE), max(colon_raw$status, na.rm=TRUE),
               max(colon_raw$age,    na.rm=TRUE), max(colon_raw$sex,    na.rm=TRUE),
               max(colon_raw$nodes,  na.rm=TRUE)),
  Mean     = round(c(mean(colon_raw$time,   na.rm=TRUE), mean(colon_raw$status, na.rm=TRUE),
                     mean(colon_raw$age,    na.rm=TRUE), mean(colon_raw$sex,    na.rm=TRUE),
                     mean(colon_raw$nodes,  na.rm=TRUE)), 2),
  Median   = round(c(median(colon_raw$time,   na.rm=TRUE), median(colon_raw$status, na.rm=TRUE),
                     median(colon_raw$age,    na.rm=TRUE), median(colon_raw$sex,    na.rm=TRUE),
                     median(colon_raw$nodes,  na.rm=TRUE)), 2),
  Missing  = c(sum(is.na(colon_raw$time)),   sum(is.na(colon_raw$status)),
               sum(is.na(colon_raw$age)),    sum(is.na(colon_raw$sex)),
               sum(is.na(colon_raw$nodes)))
)

# ── Table 2: Sample summary ───────────────────────────────────
n_total    <- nrow(colon_clean)
n_events   <- sum(colon_clean$status == 1)
n_censored <- sum(colon_clean$status == 0)

tbl2_data <- data.frame(
  Statistic = c("Total Patients","Events (Deaths)","Censored Observations",
                "Percentage Events (%)","Percentage Censored (%)",
                "Minimum Follow-up (days)","Maximum Follow-up (days)",
                "Mean Follow-up (days)","Median Follow-up (days)"),
  Value = c(n_total, n_events, n_censored,
            paste0(round(100*n_events/n_total,1),"%"),
            paste0(round(100*n_censored/n_total,1),"%"),
            min(colon_clean$time), max(colon_clean$time),
            round(mean(colon_clean$time),1), round(median(colon_clean$time),1))
)

# ── Save ──────────────────────────────────────────────────────
write.csv(tbl1_data, file.path(tbl_dir, "Table1_Variable_Summary.csv"), row.names=FALSE)
write.csv(tbl2_data, file.path(tbl_dir, "Table2_Sample_Summary.csv"),   row.names=FALSE)

save_kable(
  kable(tbl1_data, caption="Table 1: Summary Statistics of Key Variables (Pre-Imputation)",
        align="lrrrrr") %>%
    kable_styling(bootstrap_options=c("striped","hover","condensed","bordered"),
                  full_width=FALSE, position="center") %>%
    column_spec(1, bold=TRUE),
  file = file.path(tbl_dir, "Table1_Variable_Summary.html")
)

save_kable(
  kable(tbl2_data, caption="Table 2: Sample Summary — Colon Cancer Dataset",
        col.names=c("Statistic","Value"), align="lr") %>%
    kable_styling(bootstrap_options=c("striped","hover","condensed","bordered"),
                  full_width=FALSE, position="center") %>%
    column_spec(1, bold=TRUE),
  file = file.path(tbl_dir, "Table2_Sample_Summary.html")
)

cat("Tables 1 and 2 saved (CSV + HTML).\n")
