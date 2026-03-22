# ============================================================
# Script  : 03_model_fitting.R
# Author  : Joshua Pius Opio
# Task    : Fit 7 parametric models, produce Table 3,
#           Figure 1 (AIC), Figure 2 (LogLik), Figure 3 (curves)
# Covers  : Lines 151-269 of Group_2_SA1.R
# ============================================================

library(survival)
library(flexsurv)
library(ggplot2)

colon_clean   <- read.csv("data/colon_clean.csv")
SurvObj_colon <- Surv(time = colon_clean$time, event = colon_clean$status)

base_dir <- getwd()
fig_dir  <- file.path(base_dir, "outputs", "figures")
tbl_dir  <- file.path(base_dir, "outputs", "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tbl_dir, recursive=TRUE, showWarnings=FALSE)

# ── Fit all models ────────────────────────────────────────────
fits <- list(
  Exponential = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="exp"),
  Weibull     = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="weibull"),
  LogNormal   = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="lnorm"),
  LogLogistic = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="llogis"),
  Gamma       = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="gamma"),
  GenGamma    = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="gengamma"),
  Gompertz    = flexsurvreg(SurvObj_colon ~ 1, data=colon_clean, dist="gompertz")
)

# ── Table 3 ───────────────────────────────────────────────────
extract_params <- function(model_name, fit) {
  params    <- fit$res[, "est"]
  param_str <- paste(paste(rownames(fit$res), round(params,4), sep=" = "), collapse="; ")
  data.frame(Model=model_name, Parameters=param_str,
             LogLik=round(as.numeric(logLik(fit)),3),
             AIC=round(AIC(fit),3), BIC=round(BIC(fit),3),
             stringsAsFactors=FALSE)
}

tbl3_data <- do.call(rbind, mapply(extract_params, names(fits), fits, SIMPLIFY=FALSE))
tbl3_data <- tbl3_data[order(tbl3_data$AIC), ]
rownames(tbl3_data) <- NULL

best_model_name <- tbl3_data$Model[1]
cat("Best model:", best_model_name, "(AIC =", tbl3_data$AIC[1], ")\n")

write.csv(tbl3_data, file.path(tbl_dir, "Table3_Model_Comparison.csv"), row.names=FALSE)

# ── Figure 1: AIC bar chart ───────────────────────────────────
aic_df       <- data.frame(Model=tbl3_data$Model, AIC=tbl3_data$AIC)
aic_df$Model <- factor(aic_df$Model, levels=aic_df$Model[order(aic_df$AIC)])

p_aic <- ggplot(aic_df, aes(x=Model, y=AIC, fill=Model==best_model_name)) +
  geom_bar(stat="identity", width=0.6, color="black", linewidth=0.3) +
  geom_text(aes(label=round(AIC,1)), vjust=-0.4, size=3.5, fontface="bold") +
  scale_fill_manual(values=c("TRUE"="#2C6E49","FALSE"="#A8DADC"), guide="none") +
  labs(subtitle=paste("Best model:", best_model_name, "(lowest AIC, highlighted)"),
       x="Model", y="AIC") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=20,hjust=1))

ggsave(file.path(fig_dir,"Figure1_AIC_Comparison.png"),
       p_aic, width=8, height=5, dpi=300, bg="white")

# ── Figure 2: Log-likelihood dot plot ────────────────────────
loglik_df       <- data.frame(Model=tbl3_data$Model, LogLik=tbl3_data$LogLik)
loglik_df$Model <- factor(loglik_df$Model,
                           levels=loglik_df$Model[order(loglik_df$LogLik, decreasing=TRUE)])

p_loglik <- ggplot(loglik_df, aes(x=LogLik, y=Model, color=Model==best_model_name)) +
  geom_point(size=5) +
  geom_vline(xintercept=max(loglik_df$LogLik), linetype="dashed", color="#2C6E49") +
  geom_text(aes(label=round(LogLik,1)), vjust=-1.2, size=3.5) +
  scale_color_manual(values=c("TRUE"="#2C6E49","FALSE"="black"), guide="none") +
  labs(subtitle="Higher (less negative) = better fit", x="Log-Likelihood", y="Model") +
  theme_bw(base_size=13)

ggsave(file.path(fig_dir,"Figure2_LogLikelihood_Comparison.png"),
       p_loglik, width=8, height=5, dpi=300, bg="white")

# ── Figure 3: Survival curves overlay ────────────────────────
colors_models        <- c("#E63946","#457B9D","#2C6E49","#F4A261","#9B2226","#6A0572","#3D405B")
names(colors_models) <- names(fits)
time_seq_all         <- seq(0, max(colon_clean$time, na.rm=TRUE), by=10)

surv_curves <- do.call(rbind, lapply(names(fits), function(m) {
  s <- summary(fits[[m]], t=time_seq_all, type="survival", ci=FALSE)[[1]]
  data.frame(Time=s$time, Survival=s$est, Model=m)
}))

p_surv_compare <- ggplot(surv_curves, aes(x=Time, y=Survival, color=Model)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=colors_models) +
  labs(x="Time (days)", y="Survival Probability S(t)", color="Model") +
  theme_bw(base_size=13) +
  theme(legend.position="right")

ggsave(file.path(fig_dir,"Figure3_Survival_Curves_AllModels.png"),
       p_surv_compare, width=9, height=5, dpi=300, bg="white")

cat("Table 3 and Figures 1-3 saved.\n")
