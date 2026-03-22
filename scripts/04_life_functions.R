# ============================================================
# Script  : 04_life_functions.R
# Author  : Djadida Uwituze
# Task    : Compute life functions for best model (GenGamma),
#           produce Tables 4 & 5 and Figure 4 (2x2 grid)
# Covers  : Lines 272-375 of Group_2_SA1.R
# ============================================================

library(survival)
library(flexsurv)
library(ggplot2)
library(gridExtra)

colon_clean   <- read.csv("data/colon_clean.csv")
SurvObj_colon <- Surv(time=colon_clean$time, event=colon_clean$status)

base_dir <- getwd()
fig_dir  <- file.path(base_dir, "outputs", "figures")
tbl_dir  <- file.path(base_dir, "outputs", "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tbl_dir, recursive=TRUE, showWarnings=FALSE)

# ── Fit best model (GenGamma) ─────────────────────────────────
best_model      <- flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="gengamma")
best_model_name <- "GenGamma"

# ── Life functions ────────────────────────────────────────────
time_seq <- seq(1, max(colon_clean$time, na.rm=TRUE), by=10)

S_t <- summary(best_model, t=time_seq, type="survival", ci=FALSE)[[1]]$est
h_t <- summary(best_model, t=time_seq, type="hazard",   ci=FALSE)[[1]]$est
F_t <- 1 - S_t
f_t <- h_t * S_t

# ── Mean & Variance (trapezoidal integration) ─────────────────
dt     <- diff(time_seq)[1]
mean_T <- sum(S_t) * dt
ET2    <- 2 * sum(time_seq * S_t) * dt
var_T  <- ET2 - mean_T^2

cat(sprintf("Mean survival : %.2f days (%.1f years)\n", mean_T, mean_T/365))
cat(sprintf("Variance      : %.2f days^2\n", var_T))
cat(sprintf("SD            : %.2f days\n",   sqrt(var_T)))

# ── Table 4 ───────────────────────────────────────────────────
selected_times <- c(100, 200, 365, 500, 730, 1000, 1500, 2000, 2500, 3000)
lf_idx         <- sapply(selected_times, function(tp) which.min(abs(time_seq - tp)))

tbl4_data <- data.frame(
  `Time (days)` = time_seq[lf_idx],
  `S(t)`        = round(S_t[lf_idx], 4),
  `h(t)`        = round(h_t[lf_idx], 6),
  `f(t)`        = round(f_t[lf_idx], 6),
  `F(t)`        = round(F_t[lf_idx], 4),
  check.names   = FALSE
)

# ── Table 5 ───────────────────────────────────────────────────
tbl5_data <- data.frame(
  Statistic = c("Mean Survival Time (days)", "Mean Survival Time (years)",
                "Variance (days\u00b2)",     "Standard Deviation (days)"),
  Value     = c(round(mean_T,2), round(mean_T/365,2),
                round(var_T,2),  round(sqrt(var_T),2))
)

write.csv(tbl4_data, file.path(tbl_dir,"Table4_Life_Functions.csv"), row.names=FALSE)
write.csv(tbl5_data, file.path(tbl_dir,"Table5_Mean_Variance.csv"),  row.names=FALSE)

# ── Figure 4: 2x2 life function plots ────────────────────────
df_life <- data.frame(Time=time_seq, S=S_t, h=h_t, f=f_t, F=F_t)

p_St <- ggplot(df_life, aes(x=Time, y=S)) +
  geom_line(color="#457B9D", linewidth=1.2) +
  labs(title=paste("Survival Function —", best_model_name), x="Time (days)", y="S(t)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold", size=11))

p_ht <- ggplot(df_life, aes(x=Time, y=h)) +
  geom_line(color="#E63946", linewidth=1.2) +
  labs(title=paste("Hazard Function —", best_model_name), x="Time (days)", y="h(t)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold", size=11))

p_ft <- ggplot(df_life, aes(x=Time, y=f)) +
  geom_line(color="#2C6E49", linewidth=1.2) +
  labs(title=paste("Density Function —", best_model_name), x="Time (days)", y="f(t)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold", size=11))

p_Ft <- ggplot(df_life, aes(x=Time, y=F)) +
  geom_line(color="#9B2226", linewidth=1.2) +
  labs(title=paste("CDF —", best_model_name), x="Time (days)", y="F(t)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold", size=11))

p_life_grid <- gridExtra::grid.arrange(p_St, p_ht, p_ft, p_Ft, nrow=2)

ggsave(file.path(fig_dir,"Figure4_LifeFunctions_BestModel.png"),
       p_life_grid, width=10, height=8, dpi=300, bg="white")

cat("Tables 4-5 and Figure 4 saved.\n")
