
# PARAMETRIC SURVIVAL ANALYSIS — GROUP 2
# Dataset : Colon Cancer (from survival package)
# Event   : Death (etype == 2)
# ============================================================

rm(list = ls())

# ============================================================
library(survival)
library(flexsurv)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(knitr)
library(kableExtra)

# ============================================================
# ── Saving outputs ─────────────

base_dir  <- getwd()
fig_dir   <- file.path(base_dir, "outputs", "figures")
tbl_dir   <- file.path(base_dir, "outputs", "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

cat("Figures will be saved to :", fig_dir, "\n")
cat("Tables  will be saved to :", tbl_dir, "\n\n")


# Data Exploration
# ============================================================

colon <- survival::colon

colon_raw <- colon %>%
  filter(etype == 2) %>%
  distinct(id, .keep_all = TRUE)

cat("=== Dimensions (death records, one row per patient) ===\n")
dim(colon_raw)

cat("\n=== Missing Values per Variable (before imputation) ===\n")
mv_before <- colSums(is.na(colon_raw))
print(mv_before[mv_before > 0])


# ── Table 1: Summary Statistics of Key Variables ────────────
tbl1_data <- data.frame(
  Variable  = c("Time (days)", "Status", "Age (years)", "Sex", "Nodes"),
  Min       = c(min(colon_raw$time,   na.rm = TRUE),
                min(colon_raw$status, na.rm = TRUE),
                min(colon_raw$age,    na.rm = TRUE),
                min(colon_raw$sex,    na.rm = TRUE),
                min(colon_raw$nodes,  na.rm = TRUE)),
  Max       = c(max(colon_raw$time,   na.rm = TRUE),
                max(colon_raw$status, na.rm = TRUE),
                max(colon_raw$age,    na.rm = TRUE),
                max(colon_raw$sex,    na.rm = TRUE),
                max(colon_raw$nodes,  na.rm = TRUE)),
  Mean      = round(c(mean(colon_raw$time,   na.rm = TRUE),
                      mean(colon_raw$status, na.rm = TRUE),
                      mean(colon_raw$age,    na.rm = TRUE),
                      mean(colon_raw$sex,    na.rm = TRUE),
                      mean(colon_raw$nodes,  na.rm = TRUE)), 2),
  Median    = round(c(median(colon_raw$time,   na.rm = TRUE),
                      median(colon_raw$status, na.rm = TRUE),
                      median(colon_raw$age,    na.rm = TRUE),
                      median(colon_raw$sex,    na.rm = TRUE),
                      median(colon_raw$nodes,  na.rm = TRUE)), 2),
  Missing   = c(sum(is.na(colon_raw$time)),
                sum(is.na(colon_raw$status)),
                sum(is.na(colon_raw$age)),
                sum(is.na(colon_raw$sex)),
                sum(is.na(colon_raw$nodes)))
)

kable(tbl1_data,
      caption = "Table 1: Summary Statistics of Key Variables (Pre-Imputation)",
      align   = "lrrrrr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") %>%
  column_spec(1, bold = TRUE)


# Data Preparation — LOCF and Survival Object Definition
# ============================================================

# ── Handling missings ──────────────────
locf <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1) x[i] <- x[i - 1]
  }

  for (i in rev(seq_along(x))) {
    if (is.na(x[i]) && i < length(x)) x[i] <- x[i + 1]
  }
  return(x)
}

colon_clean <- colon_raw %>%
  arrange(id) %>%
  mutate(nodes = locf(nodes))

cat("=== Missing Values After LOCF Imputation ===\n")
print(colSums(is.na(colon_clean[, c("time","status","age","sex","nodes")])))

# ── Table 2: Sample Summary ───────────────────────────────────
n_total    <- nrow(colon_clean)
n_events   <- sum(colon_clean$status == 1, na.rm = TRUE)
n_censored <- sum(colon_clean$status == 0, na.rm = TRUE)
pct_event  <- round(100 * n_events  / n_total, 1)
pct_cens   <- round(100 * n_censored / n_total, 1)

tbl2_data <- data.frame(
  Statistic = c("Total Patients",
                "Events (Deaths)",
                "Censored Observations",
                "Percentage Events (%)",
                "Percentage Censored (%)",
                "Minimum Follow-up (days)",
                "Maximum Follow-up (days)",
                "Mean Follow-up (days)",
                "Median Follow-up (days)"),
  Value     = c(n_total,
                n_events,
                n_censored,
                paste0(pct_event,  "%"),
                paste0(pct_cens,   "%"),
                min(colon_clean$time,    na.rm = TRUE),
                max(colon_clean$time,    na.rm = TRUE),
                round(mean(colon_clean$time,   na.rm = TRUE), 1),
                round(median(colon_clean$time, na.rm = TRUE), 1))
)

kable(tbl2_data,
      caption   = "Table 2: Sample Summary — Colon Cancer Dataset (Post-Cleaning)",
      col.names = c("Statistic", "Value"),
      align     = "lr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") %>%
  column_spec(1, bold = TRUE)


# ── Create Survival Object ────────────────────────────────────
SurvObj_colon <- Surv(time  = colon_clean$time,
                      event = colon_clean$status)


# Fit Multiple Parametric Models
# ============================================================

fits <- list(
  Exponential  = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "exp"),
  Weibull      = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "weibull"),
  LogNormal    = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "lnorm"),
  LogLogistic  = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "llogis"),
  Gamma        = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "gamma"),
  GenGamma     = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "gengamma"),
  Gompertz     = flexsurvreg(SurvObj_colon ~ 1, data = colon_clean, dist = "gompertz")
)

# ── Table 3: Model Comparison ─────────────────────────────────
extract_params <- function(model_name, fit) {
  params     <- fit$res[, "est"]
  param_str  <- paste(
    paste(rownames(fit$res), round(params, 4), sep = " = "),
    collapse = "; "
  )
  data.frame(
    Model      = model_name,
    Parameters = param_str,
    LogLik     = round(as.numeric(logLik(fit)), 3),
    AIC        = round(AIC(fit), 3),
    BIC        = round(BIC(fit), 3),
    stringsAsFactors = FALSE
  )
}

tbl3_data <- do.call(rbind, mapply(extract_params,
                                   names(fits), fits,
                                   SIMPLIFY = FALSE))
tbl3_data <- tbl3_data[order(tbl3_data$AIC), ]
rownames(tbl3_data) <- NULL

kable(tbl3_data,
      caption = "Table 3: Parametric Model Comparison — MLE Parameters, Log-Likelihood, AIC & BIC",
      align   = "llrrr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") #%>%
  # column_spec(1, bold = TRUE) %>%
  # row_spec(1, bold = TRUE, color = "white", background = "#2C6E49")


# Model Selection
# ============================================================

best_model_name <- tbl3_data$Model[1]
best_model      <- fits[[best_model_name]]

cat("\n=== Best Model by AIC:", best_model_name, "===\n")
print(best_model)

# ── Figure 1: AIC Bar Chart ───────────────────────────────────
aic_df        <- data.frame(Model = tbl3_data$Model, AIC = tbl3_data$AIC)
aic_df$Model  <- factor(aic_df$Model, levels = aic_df$Model[order(aic_df$AIC)])

p_aic <- ggplot(aic_df, aes(x = Model, y = AIC,
                             fill = Model == best_model_name)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = round(AIC, 1)), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#2C6E49", "FALSE" = "#A8DADC"),
                    guide = "none") +
  labs(#title    = "Figure 1: AIC Comparison Across Parametric Models",
       subtitle = paste("Best model:", best_model_name, "(lowest AIC, highlighted)"),
       x = "Model", y = "AIC") +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 20, hjust = 1))

print(p_aic)

# ── Figure 2: Log-Likelihood Dot Plot ────────────────────────
loglik_df        <- data.frame(Model = tbl3_data$Model, LogLik = tbl3_data$LogLik)
loglik_df$Model  <- factor(loglik_df$Model,
                            levels = loglik_df$Model[order(loglik_df$LogLik,
                                                           decreasing = TRUE)])

p_loglik <- ggplot(loglik_df, aes(x = LogLik, y = Model,
                                   color = Model == best_model_name)) +
  geom_point(size = 5) +
  geom_vline(xintercept = max(loglik_df$LogLik),
             linetype = "dashed", color = "#2C6E49") +
  geom_text(aes(label = round(LogLik, 1)), vjust = -1.2, size = 3.5) +
  scale_color_manual(values = c("TRUE" = "#2C6E49", "FALSE" = "black"),
                     guide = "none") +
  labs(#title    = "Figure 2: Log-Likelihood Comparison Across Parametric Models",
       subtitle = "Higher (less negative) = better fit",
       x = "Log-Likelihood", y = "Model") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_loglik)

# ── Figure 3: Survival Curves — All Models ───────────────────
time_seq_all <- seq(0, max(colon_clean$time, na.rm = TRUE), by = 10)

colors_models <- c("#E63946","#457B9D","#2C6E49","#F4A261",
                   "#9B2226","#6A0572","#3D405B")
names(colors_models) <- names(fits)

surv_curves <- do.call(rbind, lapply(names(fits), function(m) {
  s <- summary(fits[[m]], t = time_seq_all, type = "survival", ci = FALSE)[[1]]
  data.frame(Time = s$time, Survival = s$est, Model = m)
}))

p_surv_compare <- ggplot(surv_curves, aes(x = Time, y = Survival, color = Model)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_models) +
  labs(#title  = "Figure 3: Fitted Survival Curves — All Parametric Models",
       x      = "Time (days)",
       y      = "Survival Probability S(t)",
       color  = "Model") +
  theme_bw(base_size = 13) +
  theme(plot.title     = element_text(face = "bold"),
        legend.position = "right")

print(p_surv_compare)


# Life Functions of the Best Model
# ============================================================

time_seq <- seq(1, max(colon_clean$time, na.rm = TRUE), by = 10)

surv_best   <- summary(best_model, t = time_seq, type = "survival", ci = FALSE)[[1]]
haz_best    <- summary(best_model, t = time_seq, type = "hazard",   ci = FALSE)[[1]]

S_t <- surv_best$est
h_t <- haz_best$est
F_t <- 1 - S_t
f_t <- h_t * S_t


# ── Mean & Variance ─────────────────
dt     <- diff(time_seq)[1]
mean_T <- sum(S_t) * dt
ET2    <- 2 * sum(time_seq * S_t) * dt
var_T  <- ET2 - mean_T^2

cat(sprintf("\nMean survival time  : %.2f days (%.1f years)\n",
            mean_T, mean_T / 365))
cat(sprintf("Variance            : %.2f days^2\n", var_T))
cat(sprintf("Standard deviation  : %.2f days\n",  sqrt(var_T)))



# ── Table 4: Life Functions at Selected Time Points ───────────
selected_times <- c(100, 200, 365, 500, 730, 1000, 1500, 2000, 2500, 3000)
lf_idx         <- sapply(selected_times,
                         function(tp) which.min(abs(time_seq - tp)))

tbl4_data <- data.frame(
  `Time (days)` = time_seq[lf_idx],
  `S(t)`        = round(S_t[lf_idx], 4),
  `h(t)`        = round(h_t[lf_idx], 6),
  `f(t)`        = round(f_t[lf_idx], 6),
  `F(t)`        = round(F_t[lf_idx], 4),
  check.names   = FALSE
)

kable(tbl4_data,
      caption = paste("Table 4: Life Functions at Selected Time Points —",
                      best_model_name, "Model"),
      align   = "rrrrr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, "Survival" = 1, "Hazard" = 1,
                     "Density" = 1, "CDF" = 1)) %>%
  column_spec(1, bold = TRUE)

# ── Table 5: Mean & Variance ──────────────────────────────────
tbl5_data <- data.frame(
  Statistic = c("Mean Survival Time (days)",
                "Mean Survival Time (years)",
                "Variance (days\u00b2)",
                "Standard Deviation (days)"),
  Value     = c(round(mean_T,        2),
                round(mean_T / 365,  2),
                round(var_T,         2),
                round(sqrt(var_T),   2))
)

kable(tbl5_data,
      caption = paste("Table 5: Mean and Variance of Survival Time —",
                      best_model_name, "Model"),
      align   = "lr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") %>%
  column_spec(1, bold = TRUE)

# ── Figure 4: 2×2 Life Function Plots ────────────────────────

df_life <- data.frame(Time = time_seq, S = S_t, h = h_t, f = f_t, F = F_t)

p_St <- ggplot(df_life, aes(x = Time, y = S)) +
  geom_line(color = "#457B9D", linewidth = 1.2) +
  labs(title = paste("Survival Function —", best_model_name),
       x = "Time (days)", y = "S(t)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_ht <- ggplot(df_life, aes(x = Time, y = h)) +
  geom_line(color = "#E63946", linewidth = 1.2) +
  labs(title = paste("Hazard Function —", best_model_name),
       x = "Time (days)", y = "h(t)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_ft <- ggplot(df_life, aes(x = Time, y = f)) +
  geom_line(color = "#2C6E49", linewidth = 1.2) +
  labs(title = paste("Density Function —", best_model_name),
       x = "Time (days)", y = "f(t)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_Ft <- ggplot(df_life, aes(x = Time, y = F)) +
  geom_line(color = "#9B2226", linewidth = 1.2) +
  labs(title = paste("CDF —", best_model_name),
       x = "Time (days)", y = "F(t)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_life_grid <- gridExtra::grid.arrange(p_St, p_ht, p_ft, p_Ft, nrow = 2)


# MLE Parameter Table & Interpretation
# ============================================================

# ── Table 6: MLE Parameter Estimates ─────────────────────────

param_raw <- as.data.frame(best_model$res)

tbl6_data <- data.frame(
  Parameter       = rownames(param_raw),
  Estimate        = round(param_raw$est,    4),
  `Std. Error`    = round(param_raw$se,     4),
  `95% CI Lower`  = round(param_raw$`L95%`, 4),
  `95% CI Upper`  = round(param_raw$`U95%`, 4),
  check.names     = FALSE
)

kable(tbl6_data,
      caption = paste("Table 6: MLE Parameter Estimates —",
                      best_model_name, "Model"),
      align   = "lrrrr") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                full_width = FALSE, position = "center") %>%
  column_spec(1, bold = TRUE)

cat("\n=== Full Model Summary ===\n")
print(best_model)


# ── Printed Interpretation ───────────────────────────────────
cat("\n========================================================\n")
cat("INTERPRETATION SUMMARY\n")
cat("========================================================\n")
cat("Model     :", best_model_name, "\n")
cat("Dataset   : Colon Cancer\n")
cat("Patients  :", n_total, " | Events:", n_events,
    " | Censored:", n_censored, "\n\n")
cat("Shape parameter:\n")
cat("  > 1 → Hazard INCREASES over time (aging / deterioration)\n")
cat("  < 1 → Hazard DECREASES over time (early risk dominates)\n")
cat("  = 1 → Constant hazard (Exponential)\n\n")
cat("Scale parameter: sets the characteristic time scale of survival.\n\n")
cat("Life Functions:\n")
cat("  S(t) : Probability of surviving beyond time t\n")
cat("  h(t) : Instantaneous risk of the event at time t\n")
cat("  f(t) : Probability density of the event near time t\n")
cat("  F(t) : Probability of having experienced the event by time t\n\n")
cat(sprintf("Mean survival : %.2f days (%.1f years)\n", mean_T, mean_T / 365))
cat(sprintf("Std deviation : %.2f days\n", sqrt(var_T)))
cat("========================================================\n")



# Saving all figs & tables
# ============================================================

cat("\n\n=== SAVING OUTPUTS ===\n")

# ── Figures ──────────────────────────────────────────────────
ggsave(file.path(fig_dir, "Figure1_AIC_Comparison.png"),
       plot = p_aic, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: Figure1_AIC_Comparison.png\n")

ggsave(file.path(fig_dir, "Figure2_LogLikelihood_Comparison.png"),
       plot = p_loglik, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: Figure2_LogLikelihood_Comparison.png\n")

ggsave(file.path(fig_dir, "Figure3_Survival_Curves_AllModels.png"),
       plot = p_surv_compare, width = 9, height = 5, dpi = 300, bg = "white")
cat("Saved: Figure3_Survival_Curves_AllModels.png\n")

ggsave(file.path(fig_dir, "Figure4_LifeFunctions_BestModel.png"),
       plot = p_life_grid, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: Figure4_LifeFunctions_BestModel.png\n")

# ── Tables as CSV ─────────────────────────────────────────────
write.csv(tbl1_data, file.path(tbl_dir, "Table1_Variable_Summary.csv"),  row.names = FALSE)
write.csv(tbl2_data, file.path(tbl_dir, "Table2_Sample_Summary.csv"),    row.names = FALSE)
write.csv(tbl3_data, file.path(tbl_dir, "Table3_Model_Comparison.csv"),  row.names = FALSE)
write.csv(tbl4_data, file.path(tbl_dir, "Table4_Life_Functions.csv"),    row.names = FALSE)
write.csv(tbl5_data, file.path(tbl_dir, "Table5_Mean_Variance.csv"),     row.names = FALSE)
write.csv(tbl6_data, file.path(tbl_dir, "Table6_MLE_Parameters.csv"),    row.names = FALSE)
cat("All CSV tables saved.\n")


# ── Tables as HTML ───────────────────
save_kable(
  kable(tbl1_data,
        caption = "Table 1: Summary Statistics of Key Variables (Pre-Imputation)",
        align   = "lrrrrr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    column_spec(1, bold = TRUE),
  file = file.path(tbl_dir, "Table1_Variable_Summary.html")
)

save_kable(
  kable(tbl2_data,
        caption   = "Table 2: Sample Summary — Colon Cancer Dataset",
        col.names = c("Statistic", "Value"), align = "lr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    column_spec(1, bold = TRUE),
  file = file.path(tbl_dir, "Table2_Sample_Summary.html")
)

save_kable(
  kable(tbl3_data,
        caption = "Table 3: Parametric Model Comparison — MLE, Log-Likelihood, AIC & BIC",
        align   = "llrrr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(1, bold = TRUE, color = "white", background = "#2C6E49"),
  file = file.path(tbl_dir, "Table3_Model_Comparison.html")
)

save_kable(
  kable(tbl4_data,
        caption = paste("Table 4: Life Functions at Selected Time Points —",
                        best_model_name, "Model"),
        align   = "rrrrr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    add_header_above(c(" " = 1, "Survival" = 1, "Hazard" = 1,
                       "Density" = 1, "CDF" = 1)) %>%
    column_spec(1, bold = TRUE),
  file = file.path(tbl_dir, "Table4_Life_Functions.html")
)

save_kable(
  kable(tbl5_data,
        caption = paste("Table 5: Mean and Variance of Survival Time —",
                        best_model_name, "Model"),
        align   = "lr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    column_spec(1, bold = TRUE),
  file = file.path(tbl_dir, "Table5_Mean_Variance.html")
)

save_kable(
  kable(tbl6_data,
        caption = paste("Table 6: MLE Parameter Estimates —",
                        best_model_name, "Model"),
        align   = "lrrrr") %>%
    kable_styling(bootstrap_options = c("striped","hover","condensed","bordered"),
                  full_width = FALSE, position = "center") %>%
    column_spec(1, bold = TRUE),
  file = file.path(tbl_dir, "Table6_MLE_Parameters.html")
)

cat("All HTML tables saved.\n")


# Final confrmation
# ============================================================
cat("\n========================================================\n")
cat("ALL OUTPUTS SAVED SUCCESSFULLY\n")
cat("========================================================\n")
cat("Root folder :", file.path(base_dir, "outputs"), "\n\n")
cat("FIGURES (", fig_dir, ")\n", sep = "")
cat("  Figure1_AIC_Comparison.png\n")
cat("  Figure2_LogLikelihood_Comparison.png\n")
cat("  Figure3_Survival_Curves_AllModels.png\n")
cat("  Figure4_LifeFunctions_BestModel.png\n\n")
cat("TABLES  (", tbl_dir, ")\n", sep = "")
cat("  Table1_Variable_Summary     .csv  .html\n")
cat("  Table2_Sample_Summary       .csv  .html\n")
cat("  Table3_Model_Comparison     .csv  .html\n")
cat("  Table4_Life_Functions       .csv  .html\n")
cat("  Table5_Mean_Variance        .csv  .html\n")
cat("  Table6_MLE_Parameters       .csv  .html\n")
cat("========================================================\n")
