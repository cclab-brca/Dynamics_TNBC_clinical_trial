# Load necessary libraries
library(dplyr)
library(lmerTest)
library(ggplot2)
library(forcats)
library(ggpubr)
library(DescTools)

# Data preprocessing
PDTX_TP1 <- PDTX_TP1 %>% filter(!Mouse %in% c("AN17CUK012247"))  # Exclude specific mice due to experimental issue



# Initialize dataframe to store fixed effects summary
fixed_effects_summary <- data.frame(Model = character(),
                                    Parameter = character(),
                                    Estimate = numeric(),
                                    StdError = numeric(),
                                    Pr = numeric(),
                                    stringsAsFactors = FALSE)

# Iterate over each unique model in the dataset
for (i in unique(P1MT_NZ$Model)) {
  tumrs <- subset(P1MT_NZ, Model == i)
  tumrs$Size <- log(tumrs$Size)
  tumrs$Size[which(!is.finite(tumrs$Size))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")  # Adjust reference group as needed
  
  # Fit the mixed-effects model
  m_t <- lmerTest::lmer(Size ~ Time * Group + (Time | Mouse), data = tumrs, 
                        control = lmerControl(optimizer = "Nelder_Mead"))
  
  # Extract fixed effects and store in summary dataframe
  summary_m_t <- summary(m_t)
  fixed_eff_df <- as.data.frame(summary_m_t$coefficients)
  fixed_eff_df$Parameter <- rownames(summary_m_t$coefficients)
  fixed_eff_df$Model <- i
  fixed_effects_summary <- rbind(fixed_effects_summary, fixed_eff_df[, c("Model", "Parameter", "Estimate", "Std. Error", "Pr(>|t|)")])
}

# Reset row names for the fixed effects summary
rownames(fixed_effects_summary) <- NULL

# Add confidence intervals for fixed effects
fixed_effects_summary <- fixed_effects_summary %>%
  mutate(CI_low = NA, CI_high = NA)

# Iterate over each unique model to calculate confidence intervals
for (i in unique(fixed_effects_summary$Model)) {
  tumrs <- subset(P1MT_NZ, Model == i)
  tumrs$Size <- log(tumrs$Size)
  tumrs$Size[which(!is.finite(tumrs$Size))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")
  
  m_t <- lmerTest::lmer(Size ~ Time * Group + (Time | Mouse), data = tumrs, 
                        control = lmerControl(optimizer = "bobyqa"))
  
  ci <- confint(m_t, method = "boot", level = 0.95)
  
  for (param in fixed_effects_summary$Parameter[fixed_effects_summary$Model == i]) {
    if (param %in% rownames(ci)) {
      ci_param <- ci[param, ]
      fixed_effects_summary$CI_low[fixed_effects_summary$Model == i & fixed_effects_summary$Parameter == param] <- ci_param[1]
      fixed_effects_summary$CI_high[fixed_effects_summary$Model == i & fixed_effects_summary$Parameter == param] <- ci_param[2]
    }
  }
}

# Assign pCR status to models
fixed_effects_summary$pCR <- "Yes"
fixed_effects_summary[fixed_effects_summary$Model %in% c("PAR1040", "PAR1008", "PAR1141"), ]$pCR <- "No"
fixed_effects_summary$pCR <- factor(fixed_effects_summary$pCR)

# Extract drug names from interaction terms in the 'Parameter' column
fixed_effects_summary <- fixed_effects_summary %>%
  mutate(Drug = case_when(
    str_detect(Parameter, "^Time:") ~ str_extract(Parameter, "(?<=Time:).*$"),
    TRUE ~ NA_character_
  )) %>%
  mutate(Drug = str_replace(Drug, "^Group", ""))

# Save the summary data to a CSV file
write.csv2(fixed_effects_summary, file = "Trial1_metric1_stats.csv", row.names = FALSE)

# Plot: Interaction term for all models in Trial 1
curve_colors <- c("CT" = "#C43826", "CTO" = "#499A87", "Olaparib" = "#858FB0", 
                  "AZD1775" = "#C95C2D", "Olap+AZD1775" = "#6EB8CE")

pdf("Trial1_metric1_all_drugs.pdf", width = 18, height = 7)
ggplot(fixed_effects_summary %>% filter(!is.na(Drug) & grepl("^Time:", Parameter)), 
       aes(x = fct_reorder(Drug, Estimate, .desc = FALSE), y = Estimate, color = Drug, group = Drug)) +
  geom_point(size = 3, position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "", y = "Change in growth rate T - UT", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.06, 0.01)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none") +
  facet_grid(~ Model, scales = "free_y", space = "free_y")
dev.off()

# Plot: Interaction term with pCR and non-pCR groups
filtered_df <- fixed_effects_summary %>%
  inner_join(data.frame(Model = c("PAR1022", "PAR1053", "PAR1006", "PAR1040", "PAR1008", "PAR1141"),
                        Drug = c("CT", "CT", "CTO", "CT", "CT", "CTO")), by = c("Model", "Drug"))

pdf("Trial1_metric1_pCR.pdf", width = 10, height = 5)
ggplot(filtered_df, aes(x = fct_reorder(Model, Estimate, .desc = FALSE), y = Estimate, color = Drug)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "Model", y = "Change in growth rate T - UT", color = "Drug", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.2, 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16)) +
  facet_wrap(~ pCR, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon test for pCR and non-pCR groups
wilcox_test <- wilcox.test(Estimate ~ pCR, data = filtered_df)

pdf("Trial1_metric1_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = pCR, y = Estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  scale_color_manual(values = curve_colors) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Comparison of Estimates between pCR and non-pCR groups",
       x = "pCR Status", y = "Change in growth rate T - UT") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = max(filtered_df$Estimate) + 0.005, 
           label = paste("Wilcoxon Test:",
                         "\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# RCB Score and Analysis for Metric 1
P1MT_NZ$RCB_Score <- case_when(
  Model == "PAR1006" ~ 0,
  Model == "PAR1040" ~ 2,
  Model == "PAR1053" ~ 0,
  Model == "PAR1022" ~ 0,
  Model == "PAR1008" ~ 1,
  Model == "PAR1141" ~ 2,
  TRUE ~ NA_real_
)

merged_df <- merge(filtered_df, unique(P1MT_NZ[, c("Model", "RCB_Score")]), by = "Model", all.x = TRUE)
merged_df$RCB_Category <- factor(case_when(
  merged_df$RCB_Score == 0 ~ "pCR",
  merged_df$RCB_Score == 1 ~ "RCB-I",
  merged_df$RCB_Score == 2 ~ "RCB-II"
), levels = c("pCR", "RCB-I", "RCB-II"))

pdf("Trial1_metric1_RCB.pdf", width = 11, height = 5)
ggplot(merged_df, aes(x = fct_reorder(Model, Estimate, .desc = FALSE), y = Estimate, color = Drug)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "Model", y = "Change in growth rate T - UT", color = "Drug", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.2, 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16),
        panel.spacing = unit(1, "lines")) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere-Terpstra Test for Metric 1
merged_df$RCB_Category <- factor(merged_df$RCB_Category, levels = c("pCR", "RCB-I", "RCB-II"), ordered = TRUE)
jt_result <- JonckheereTerpstraTest(Estimate ~ RCB_Category, data = merged_df, alternative = "increasing")

# Summary of mean Estimate values for each RCB category
summary_data <- merged_df %>%
  group_by(RCB_Category) %>%
  summarise(mean_Estimate = mean(Estimate, na.rm = TRUE), .groups = 'drop')

pdf("Trial1_metric1_RCB_JT.pdf", width = 8, height = 6)
ggplot(merged_df, aes(x = RCB_Category, y = Estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_Estimate, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_Estimate), size = 3, color = "darkblue") +
  labs(title = "Change in Growth Rate by RCB Category",
       x = "RCB Category", y = "Change in Growth Rate UT - T") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = curve_colors) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(merged_df$Estimate) + 0.005, 
           label = paste("Jonckheere-Terpstra Test:",
                         "\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 0.8)
dev.off()

# Load necessary libraries
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)
library(DescTools)

# Initialize empty lists and data frames to store models and summary statistics
model_list <- list()
emm_pairs_summary <- data.frame()

# Iterate over each unique model in the P1MT_NZ dataset
for (i in unique(P1MT_NZ$Model)) {
  tumrs <- subset(P1MT_NZ, Model == i)
  tumrs$Size <- log(tumrs$Size)
  tumrs$Size[which(!is.finite(tumrs$Size))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")
  
  # Fit the linear mixed-effects model for each model
  m_t <- lmer(Size ~ Time * Group + (Time | Mouse), data = tumrs, 
              control = lmerControl(optimizer = "bobyqa"))
  model_list[[i]] <- m_t
  
  # Compute estimated marginal means (EMM) for 'Group' at Time c(0, 77)
  emm <- emmeans(m_t, specs = ~ Group, at = list(Time = c(0, 77)))
  
  # Generate pairwise comparisons of EMMs
  emm_pairs <- pairs(emm, scale = 77)
  
  # Convert the EMM pairwise comparison summary to a data frame and append model identifier
  emm_pairs_df <- summary(emm_pairs, infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    mutate(Model = i)
  
  # Append the data frame to the overall summary
  emm_pairs_summary <- rbind(emm_pairs_summary, emm_pairs_df)
}

# Label pCR status for models
emm_pairs_summary$pCR <- "Yes"
emm_pairs_summary[emm_pairs_summary$Model %in% c("PAR1040", "PAR1008", "PAR1141"), ]$pCR <- "No"
emm_pairs_summary$pCR <- factor(emm_pairs_summary$pCR)

# Extract drug names from the 'contrast' column
emm_pairs_summary$Drugs <- sub("Untreated - ", "", emm_pairs_summary$contrast)

# Save the summary data frame to a CSV file
write.csv2(emm_pairs_summary, file = "Trial1_metric2_stats.csv", row.names = FALSE)

# Plot: Estimated Differences in AUC for Single Drug Treatments
curve_colors <- c("CT" = "#C43826", "CTO" = "#499A87", "Olaparib" = "#858FB0", 
                  "AZD1775" = "#C95C2D", "(Olap+AZD1775)" = "#6EB8CE")

single_drug_emm_summary <- emm_pairs_summary %>%
  filter(Drugs %in% c("AZD1775", "CTO", "CT", "(Olap+AZD1775)", "Olaparib")) %>%
  mutate(Model = fct_reorder(Model, estimate, .desc = TRUE))

pdf("Trial1_metric2_all_drugs.pdf", width = 35, height = 7)
ggplot(single_drug_emm_summary, aes(x = Drugs, y = estimate, fill = Drugs)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Drug", y = "Estimated Difference in AUC", fill = "Drug", title = "PDTX Models - Single Drug Treatments") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),  # Increase font size for strip text
    axis.text.x = element_text(size = 12),  # Increase font size for x-axis text
    axis.text.y = element_text(size = 12),  # Increase font size for y-axis text
    axis.title = element_text(face = "bold", size = 14),  # Increase font size for axis titles
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Increase font size for title
    legend.position = "none",  # Remove legend
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ Model, nrow = 1, scales = "free_x")
dev.off()

# Plot: AUC Estimates with pCR and non-pCR groups
single_drug_emm_summary$pCR <- factor(single_drug_emm_summary$pCR, levels = c("Yes", "No"), labels = c("pCR", "non-pCR"))

desired_combinations <- data.frame(
  Model = c("PAR1022", "PAR1053", "PAR1006", "PAR1040", "PAR1008", "PAR1141"),
  Drugs = c("CT", "CT", "CTO", "CT", "CT", "CTO")
)

filtered_df <- single_drug_emm_summary %>%
  inner_join(desired_combinations, by = c("Model", "Drugs"))

pdf("Trial1_metric2_pCR.pdf", width = 10, height = 5)
ggplot(filtered_df, aes(x = Model, y = estimate, fill = Drugs)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Estimated Difference in AUC", fill = "Drug", title = "PDTX Models by Drug") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),  # Increased font size for strip text
    axis.text.x = element_text(size = 14),  # Increased font size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased font size for y-axis text
    axis.title = element_text(face = "bold", size = 16),  # Increased font size for axis titles
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Increased font size for title
    legend.text = element_text(size = 14),  # Increased font size for legend text
    legend.title = element_text(face = "bold", size = 16)  # Increased font size for legend title
  ) +
  facet_wrap(~ pCR, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon test: Metric 2, pCR vs non-pCR
wilcox_test <- wilcox.test(estimate ~ pCR, data = filtered_df)

pdf("Trial1_metric2_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = pCR, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Comparison of Estimates between pCR and non-pCR groups",
       x = "pCR Status", y = "Estimated Difference in AUC") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),  # Increased font size for strip text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Increased font size for title
    axis.title = element_text(face = "bold", size = 16),  # Increased font size for axis titles
    axis.text.x = element_text(size = 14),  # Increased font size for x-axis text
    axis.text.y = element_text(size = 14)  # Increased font size for y-axis text
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = max(filtered_df$estimate) * 1.1, 
           label = paste("Wilcoxon Test:", "\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 0.8)
dev.off()

# RCB Score and Analysis for Metric 2
merged_df <- merge(filtered_df, unique(P1MT_NZ[, c("Model", "RCB_Score")]), by = "Model", all.x = TRUE)
merged_df$RCB_Category <- factor(ifelse(merged_df$RCB_Score == 0, "pCR",
                                        ifelse(merged_df$RCB_Score == 1, "RCB-I", "RCB-II")),
                                 levels = c("pCR", "RCB-I", "RCB-II"))

pdf("Trial1_metric2_RCB.pdf", width = 11, height = 5)
ggplot(merged_df, aes(x = fct_reorder(Model, estimate, .desc = TRUE), y = estimate, fill = Drugs)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Estimated Difference in AUC", fill = "Drug", title = "PDTX Models by Drug and RCB Category") +
  expand_limits(y = c(min(merged_df$estimate), max(merged_df$estimate))) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere-Terpstra Test for Metric 2
merged_df$RCB_Category <- factor(merged_df$RCB_Category, levels = c("pCR", "RCB-I", "RCB-II"), ordered = TRUE)
jt_result <- JonckheereTerpstraTest(estimate ~ RCB_Category, data = merged_df, alternative = "decreasing")

# Summarize mean Estimate values for each RCB category
summary_data <- merged_df %>%
  group_by(RCB_Category) %>%
  summarise(mean_estimate = mean(estimate, na.rm = TRUE), .groups = 'drop')

pdf("Trial1_metric2_RCB_JT.pdf", width = 8, height = 6)
ggplot(merged_df, aes(x = RCB_Category, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_estimate, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_estimate), size = 3, color = "darkblue") +
  labs(title = "Estimated Difference in AUC by RCB Category", x = "RCB Category", y = "Estimated Difference in AUC") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(merged_df$estimate) + 10, 
           label = paste("Jonckheere-Terpstra Test:", "\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# Load required libraries
library(lmerTest)
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)
library(DescTools)

# Initialize an empty data frame to store results
results_summary <- data.frame(Model = character(), 
                              Drug = character(),
                              Lb_hat = numeric(),
                              CI_low = numeric(), 
                              CI_high = numeric(),
                              stringsAsFactors = FALSE)

alpha = 0.05 # Confidence level

# Iterate over each unique model in the dataset
for (i in unique(P1MT_NZ$Model)) {
  # Subset data for the current model
  current_data <- filter(P1MT_NZ, Model == i)
  current_data$Size <- log(current_data$Size)
  current_data$Group <- relevel(current_data$Group, ref = "Untreated")
  
  # Fit the mixed model
  m_t <- lmer(Size ~ Time*Group + (Time | Mouse), 
              data = current_data, control = lmerControl(optimizer = "bobyqa"))
  
  # List of possible drugs (interaction terms) to consider
  drugs <- c("AZD1775", "CT", "CTO", "Olap+AZD1775", "Olaparib")
  
  # For each drug, set up K vector, calculate Lb.hat, and append results
  for (drug in drugs) {
    interaction_term <- paste0("Time:Group", drug)
    
    # Check if the model includes the interaction term for the current drug
    if (interaction_term %in% names(fixef(m_t))) {
      K <- setNames(rep(0, length(fixef(m_t))), names(fixef(m_t)))
      K["Time"] <- 1
      K[interaction_term] <- 1
      
      # Calculation
      Va <- vcov(m_t)
      b.hat <- fixef(m_t)
      Lb.hat <- sum(K * b.hat)
      SE <- sqrt(sum((K %*% Va) * K))
      ddf <- get_Lb_ddf(m_t, K) 
      confint <- Lb.hat + c(-1, 1) * SE * qt(1 - 0.5 * alpha, ddf)
      
      # Append the results
      results_summary <- rbind(results_summary, data.frame(Model = i, Drug = drug, Lb_hat = Lb.hat, CI_low = confint[1], CI_high = confint[2]))
    }
  }
}

# Add pCR information
results_summary$pCR <- "Yes"
results_summary[results_summary$Model %in% c("PAR1141", "PAR1040", "PAR1008"),]$pCR <- "No"
results_summary$pCR <- factor(results_summary$pCR)

# Save results to CSV
write.csv2(results_summary, file = "Trial1_metric3_stats.csv", row.names = FALSE)

# Plot Time and Time*Group Interaction for All Drugs
curve_colors <- c("CT" = "#C43826", "CTO" = "#499A87", "Olaparib" = "#858FB0", 
                  "AZD1775" = "#C95C2D", "Olap+AZD1775" = "#6EB8CE")

pdf("Trial1_metric3_all_drugs.pdf", width = 12, height = 7)
ggplot(results_summary, aes(x = fct_reorder(Drug, Lb_hat, .desc = FALSE), y = Lb_hat, color = Drug, group = Drug)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "", y = "Growth rate under treatment", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.2, max(results_summary$CI_high) * 1.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  ) +
  facet_grid(~ Model, scales = "free_y", space = "free_y") 
dev.off()

# Plot: Time+Time*GroupCT with pCR and non-pCR groups
results_summary$pCR <- factor(results_summary$pCR, levels = c("Yes", "No"), labels = c("pCR", "non-pCR"))

desired_combinations <- data.frame(
  Model = c("PAR1022", "PAR1053", "PAR1006", "PAR1040", "PAR1008", "PAR1141"),
  Drug = c("CT", "CT", "CTO", "CT", "CT", "CTO")
)

filtered_df <- results_summary %>%
  inner_join(desired_combinations, by = c("Model", "Drug"))

pdf("Trial1_metric3_pCR.pdf", width = 11, height = 5)
ggplot(filtered_df, aes(x = fct_reorder(Model, Lb_hat, .desc = FALSE), y = Lb_hat, color = Drug)) +
  geom_point(show.legend = TRUE, size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "Model", y = "Growth rate under treatment", color = "Drug", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.2, max(filtered_df$CI_high) * 1.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  ) +
  facet_wrap(~ pCR, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon Test: Metric 3 pCR vs non-pCR
wilcox_test <- wilcox.test(Lb_hat ~ pCR, data = filtered_df)

pdf("Trial1_metric3_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = pCR, y = Lb_hat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Growth Rate Under Treatment by pCR Status",
       x = "pCR Status",
       y = "Growth Rate Under Treatment") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = max(filtered_df$Lb_hat) * 1.1, 
           label = paste("Wilcoxon Test:", "\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# RCB Score and Analysis for Metric 3
filtered_df <- merge(filtered_df, unique(P1MT_NZ[, c("Model", "RCB_Score")]), by = "Model", all.x = TRUE)

filtered_df$RCB_Category <- ifelse(filtered_df$RCB_Score == 0, "pCR",
                                   ifelse(filtered_df$RCB_Score == 1, "RCB-I", "RCB-II"))
filtered_df$RCB_Category <- factor(filtered_df$RCB_Category, levels = c("pCR", "RCB-I", "RCB-II"))

pdf("Trial1_metric3_RCB.pdf", width = 11, height = 5)
ggplot(filtered_df, aes(x = fct_reorder(Model, Lb_hat, .desc = FALSE), y = Lb_hat, color = Drug)) +
  geom_point(show.legend = TRUE, size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "Model", y = "Growth rate under treatment", color = "Drug", 
       title = "PDTX Models - Growth Rate Under Treatment by RCB Category") +
  expand_limits(y = c(-0.2, max(filtered_df$CI_high) * 1.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  ) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere-Terpstra Test for Metric 3
filtered_df$RCB_Category <- factor(filtered_df$RCB_Category, 
                                   levels = c("pCR", "RCB-I", "RCB-II"), 
                                   ordered = TRUE)

jt_result <- JonckheereTerpstraTest(Lb_hat ~ RCB_Category, data = filtered_df, 
                                    alternative = "increasing")

summary_data <- filtered_df %>%
  group_by(RCB_Category) %>%
  summarise(mean_Lb_hat = mean(Lb_hat, na.rm = TRUE), .groups = 'drop')

pdf("Trial1_metric3_RCB_JT.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = RCB_Category, y = Lb_hat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_Lb_hat, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_Lb_hat), size = 3, color = "darkblue") +
  labs(title = "Growth Rate Under Treatment by RCB Category",
       x = "RCB Category",
       y = "Growth Rate Under Treatment") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(filtered_df$Lb_hat) * 1.1, 
           label = paste("Jonckheere-Terpstra Test Results:",
                         "\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# Load required libraries
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(DescTools)

# Initialize lists for storing models and combined data
model_list <- list()
combined_data_list <- list()

# Data processing and model fitting for each model
for (i in unique(P1MT_NZ$Model)) {
  tumrs <- subset(P1MT_NZ, Model == i)
  tumrs$lSize <- log(tumrs$Size)
  tumrs$lSize[which(!is.finite(tumrs$lSize))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")
  
  # Fit the mixed model
  m_t <- lmerTest::lmer(lSize ~ Time * Group + (Time | Mouse), data = tumrs, 
                        control = lmerControl(optimizer = "Nelder_Mead"))
  model_list[[i]] <- m_t 
  
  # Create prediction frame for Time 78
  dc <- distinct(dplyr:::select(tumrs, Mouse, Group))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 78) %>% 
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  # Combine actual and predicted data
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Model <- i
  
  combined_data_list[[i]] <- comb
}

# Combine all models' data
final_combined_data <- bind_rows(combined_data_list)

# Filter and rename data for final boxplot analysis
predicted_final_bxp <- final_combined_data %>%
  filter(Time == 78 & Group != "Untreated") %>%
  distinct(Model, Mouse, .keep_all = TRUE) %>%
  rename(Drugs = Group, Size = ModelSize)

# Add pCR information
predicted_final_bxp$pCR <- "Yes"
predicted_final_bxp[predicted_final_bxp$Model %in% c("PAR1141", "PAR1040", "PAR1008"),]$pCR <- "No"
predicted_final_bxp$pCR <- factor(predicted_final_bxp$pCR)

# Save results to CSV
write.csv(predicted_final_bxp, file = "Trial1_metric4_stats.csv", row.names = FALSE)

# Boxplot of predicted values for Trial 1 (ON)
curve_colors <- c("CT" = "#C43826", "CTO" = "#499A87", "Olaparib" = "#858FB0", 
                  "AZD1775" = "#C95C2D", "Olap+AZD1775" = "#6EB8CE")

pdf("Trial1_metric4_all_drugs.pdf", width = 17, height = 9)
ggplot(predicted_final_bxp, aes(x = fct_reorder(Model, Size, .fun = mean, .desc = TRUE), y = Size)) +
  geom_boxplot(aes(fill = Drugs), position = position_dodge(width = 1), show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Predicted volume at the end of the treatment", fill = "Drug", title = "PDTX Models") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.8, face = "bold", size = 16),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ Model, nrow = 1, scales = "free_x") +
  ylim(0, 800)
dev.off()

# Individual plots by model
all_combinations <- expand_grid(Model = unique(predicted_final_bxp$Model), Drugs = unique(predicted_final_bxp$Drugs))
predicted_final_bxp_complete <- full_join(all_combinations, predicted_final_bxp, by = c("Model", "Drugs"))

plot_list <- list()
models <- unique(predicted_final_bxp_complete$Model)

for (model in models) {
  model_data <- predicted_final_bxp_complete[predicted_final_bxp_complete$Model == model, ]
  p <- ggplot(model_data, aes(x = Drugs, y = Size)) +
    geom_boxplot(aes(fill = Drugs), position = position_dodge(width = 1), show.legend = FALSE, outlier.shape = NA) +
    scale_fill_manual(values = curve_colors, drop = FALSE) +
    labs(x = "Drug", y = "Predicted volume at \n the end of the treatment") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines")
    ) +
    coord_cartesian(ylim = c(0, 1050)) +
    facet_wrap(~ Model, scales = "fixed", ncol = 1)
  
  plot_list[[model]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(title = "PDTX Models by Drug - Metric 4", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))) &
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"), strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

pdf("Trial1_metric4_all_drugs_modf.pdf", width = 15, height = 12)
print(combined_plot)
dev.off()

# Boxplot for pCR and non-pCR for Metric 4
predicted_final_bxp$pCR <- factor(predicted_final_bxp$pCR, levels = c("Yes", "No"), labels = c("pCR", "non-pCR"))
desired_combinations <- data.frame(Model = c("PAR1022", "PAR1053", "PAR1006", "PAR1040", "PAR1008", "PAR1141"), Drugs = c("CT", "CT", "CTO", "CT", "CT", "CTO"))
filtered_df <- predicted_final_bxp %>% inner_join(desired_combinations, by = c("Model", "Drugs"))

pdf("Trial1_metric4_pCR.pdf", width = 11, height = 5)
ggplot(filtered_df, aes(x = fct_reorder(Model, lSize, .fun = mean, .desc = TRUE), y = Size, fill = Drugs)) +
  geom_boxplot(position = position_dodge(0.5), show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Predicted volume at the end \n of the treatment (Day 77)", fill = "Drug", title = "PDTX Models by Drug") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  ) +
  facet_wrap(~ pCR, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon Test for Metric 4
wilcox_test <- wilcox.test(Size ~ pCR, data = filtered_df)

pdf("Trial1_metric4_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = pCR, y = Size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  scale_color_manual(values = curve_colors) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Tumour Size by pCR Status", x = "pCR Status", y = "Predicted Tumour Size") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = max(filtered_df$Size) * 1.1, 
           label = paste("Wilcoxon Test:", "\nW =", round(wilcox_test$statistic, 2), "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# RCB_Score plots and analysis for Metric 4
filtered_df <- merge(filtered_df, unique(P1MT_NZ[, c("Model", "RCB_Score")]), by = "Model", all.x = TRUE)
filtered_df$RCB_Category <- factor(ifelse(filtered_df$RCB_Score == 0, "pCR", ifelse(filtered_df$RCB_Score == 1, "RCB-I", "RCB-II")), levels = c("pCR", "RCB-I", "RCB-II"))

pdf("Trial1_metric4_RCB.pdf", width = 11, height = 5)
ggplot(filtered_df, aes(x = fct_reorder(Model, Size, .fun = mean, .desc = TRUE), y = Size, fill = Drugs)) +
  geom_boxplot(position = position_dodge(0.5), show.legend = TRUE) +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Predicted volume at the end \n of the treatment (Day 77)", fill = "Drug", title = "PDTX Models Boxplots by RCB Category") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere's Test for Metric 4
filtered_df$RCB_Category <- factor(filtered_df$RCB_Category, levels = c("pCR", "RCB-I", "RCB-II"), ordered = TRUE)
jt_result <- JonckheereTerpstraTest(Size ~ RCB_Category, data = filtered_df, alternative = "increasing")

summary_data <- filtered_df %>%
  group_by(RCB_Category) %>%
  summarise(mean_Size = mean(Size, na.rm = TRUE), .groups = 'drop')

pdf("Trial1_metric4_RCB_JT.pdf", width = 8, height = 6)
ggplot(filtered_df, aes(x = RCB_Category, y = Size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_Size, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_Size), size = 3, color = "darkblue") +
  labs(title = "Predicted Tumour Size by RCB Category", x = "RCB Category", y = "Predicted Tumour Size") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(filtered_df$Size) * 1.1, 
           label = paste("Jonckheere-Terpstra Test Results:", "\nJT =", round(jt_result$statistic, 2), "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()


# Load necessary libraries
library(dplyr)
library(tidyr)
library(lmerTest)
library(gridExtra)
library(ggplot2)
library(forcats)
library(stringr)

# Step 1: Subset data for untreated group and groups with treatment "off"
untreated <- subset(PDTX_TP1, Group == "Untreated")
other_groups <- subset(PDTX_TP1, Group != "Untreated" & Treatment == "off")
P1MPT_UTNZ <- rbind(untreated, other_groups)

# Step 2: Initialize a data frame to store fixed effect parameters
fixed_effects_summary <- data.frame(Model = character(),
                                    Parameter = character(),
                                    Estimate = numeric(),
                                    StdError = numeric(),
                                    Pr = numeric(),
                                    stringsAsFactors = FALSE)

# Step 3: Iterate over each unique model to fit mixed models and extract fixed effects
for (i in unique(P1MPT_UTNZ$Model)) {
  tumrs <- subset(P1MPT_UTNZ, Model == i)
  tumrs$Size <- log(tumrs$Size)
  tumrs$Size[which(!is.finite(tumrs$Size))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")
  
  m_t <- lmerTest::lmer(Size ~ Time * Group + (Time | Mouse), data = tumrs, 
                        control = lmerControl(optimizer = "bobyqa"))
  
  summary_m_t <- summary(m_t)
  fixed_eff_df <- as.data.frame(summary_m_t$coefficients)
  fixed_eff_df$Parameter <- rownames(summary_m_t$coefficients)
  fixed_eff_df$Model <- i
  
  fixed_effects_summary <- rbind(fixed_effects_summary, fixed_eff_df[, c("Model", "Parameter", "Estimate", "Std. Error", "Pr(>|t|)")])
}

rownames(fixed_effects_summary) <- NULL

# Step 4: Add columns for confidence intervals
fixed_effects_summary <- fixed_effects_summary %>%
  mutate(CI_low = NA, CI_high = NA)

# Step 5: Iterate over each unique model to calculate confidence intervals using bootstrapping
for (i in unique(fixed_effects_summary$Model)) {
  tumrs <- subset(P1MPT_UTNZ, Model == i)
  tumrs$Size <- log(tumrs$Size)
  tumrs$Size[which(!is.finite(tumrs$Size))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "Untreated")
  
  m_t <- lmerTest::lmer(Size ~ Time * Group + (Time | Mouse), data = tumrs,
                        control = lmerControl(optimizer = "bobyqa"))
  
  ci <- confint(m_t, method = "boot", level = 0.95)
  
  for (param in fixed_effects_summary$Parameter[fixed_effects_summary$Model == i]) {
    if (param %in% rownames(ci)) {
      ci_param <- ci[param, ]
      fixed_effects_summary$CI_low[fixed_effects_summary$Model == i & fixed_effects_summary$Parameter == param] <- ci_param[1]
      fixed_effects_summary$CI_high[fixed_effects_summary$Model == i & fixed_effects_summary$Parameter == param] <- ci_param[2]
    }
  }
}

# Step 6: Label pCR status per model
fixed_effects_summary$pCR <- "Yes"
fixed_effects_summary[fixed_effects_summary$Model %in% c("PAR1040", "PAR1008", "PAR1141"),]$pCR <- "No"

# Step 7: Extract drug names from interaction terms
fixed_effects_summary <- fixed_effects_summary %>%
  mutate(Drug = ifelse(str_detect(Parameter, "^Time:"), str_replace(Parameter, "Time:", ""), NA_character_)) %>%
  mutate(Drug = case_when(
    str_detect(Parameter, "^Time:") ~ str_extract(Parameter, "(?<=Time:).*$"),
    TRUE ~ NA_character_
  )) %>%
  mutate(Drug = str_replace(Drug, "^Group", ""))

# Step 8: Save the results to a CSV file
write.csv2(fixed_effects_summary, file = "Trial1_metric1_PT_stats.csv", row.names = FALSE)

# Step 9: Create interaction plot for all models
# Filter data to only include interaction terms
fixefs <- fixed_effects_summary %>%
  filter(!is.na(Drug) & grepl("^Time:", Parameter))

# Define color palette for drugs
curve_colors <- c("CT" = "#C43826",
                  "CTO" = "#499A87",
                  "Olaparib" = "#858FB0",
                  "AZD1775" = "#C95C2D",
                  "Olap+AZD1775" = "#6EB8CE")

# Step 10: Create PDF output for the interaction plot
pdf("Trial1_PT_metric1.pdf", width = 12, height = 9)

ggplot(fixefs, aes(x = fct_reorder(Model, Estimate, .desc = FALSE), y = Estimate, color = Drug)) +
  geom_point(aes(shape = pCR), size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors, guide = "none") +
  scale_shape_manual(values = c("No" = 16, "Yes" = 17)) +
  labs(x = "Model", y = "Change in growth rate T - UT", 
       shape = "pCR", 
       title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.05, 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right") +
  facet_wrap(~ Drug, nrow = 1, scales = "free_x")

dev.off()

