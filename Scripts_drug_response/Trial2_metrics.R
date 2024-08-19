# Load necessary libraries
library(dplyr)
library(readr)

#
PDTX_TP2 <- PDTX_TP2 %>%
  mutate(Size = as.numeric(Size),
         Response = case_when(
           Model %in% c("PAR1040", "PAR1008", "PAR1141", "PAR1221") ~ "No",
           TRUE ~ "Yes"
         ),
         Response = factor(Response),
         Group = factor(Group))

# Filter and prepare the data
PDTX_TP2 <- PDTX_TP2 %>%
  drop_na(Size) %>%
  mutate(lSize = log(Size)) %>%
  filter(Size > 0 & lSize > 1)

# Specify conditions for filtering specific rows
conditions <- data.frame(
  Mouse = c("AN4", "AN4", "AN3", "AN3", "AN2", "AN5"),
  Time = c(56, 63, 56, 63, 63, 63),
  Model = rep("PAR1022", 6)
)

# Function to check if a row matches one of the deletion conditions
matches_condition <- function(row, conditions) {
  any(apply(conditions, 1, function(cond) {
    row$Mouse == cond["Mouse"] &&
      row$Time == cond["Time"] &&
      row$Model == cond["Model"]
  }))
}

# Use the filter function to exclude rows that match the deletion conditions
PDTX_TP2_filtered <- PDTX_TP2 %>%
  rowwise() %>%
  filter(!matches_condition(cur_data(), conditions)) %>%
  ungroup()

# Load necessary libraries
library(dplyr)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(forcats)
library(DescTools)

# Initialize the emmeans summary dataframe
emm_pairs_summary <- data.frame()

# Iterate over each unique model to fit the model and compute emmeans
model_list <- list()
for (model_name in unique(PDTX_TP2nz$Model)) {
  tumrs <- subset(PDTX_TP2nz, Model == model_name)
  tumrs$lSize <- log(tumrs$Size)
  tumrs$lSize[!is.finite(tumrs$lSize)] <- NA
  tumrs$Group <- relevel(factor(tumrs$Group), ref = "UT")
  
  m_t <- lmer(lSize ~ Time * Group + (Time | Mouse), data = tumrs, 
              control = lmerControl(optimizer = "Nelder_Mead"))
  model_list[[model_name]] <- m_t
  
  # Compute estimated marginal means and pairwise comparisons
  emm <- emmeans(m_t, "Group", at = list(Time = c(0, 84)))
  emm_pairs <- pairs(emm, scale = 84)
  
  emm_pairs_df <- as.data.frame(summary(emm_pairs, infer = c(TRUE, TRUE))) %>%
    mutate(Model = model_name)
  
  emm_pairs_summary <- rbind(emm_pairs_summary, emm_pairs_df)
}

# Add pCR and Drug classification based on Model
emm_pairs_summary <- emm_pairs_summary %>%
  mutate(
    pCR = factor(ifelse(Model %in% c("PAR1040", "PAR1008", "PAR1141", "PAR1221"), "No", "Yes")),
    Drug = factor(ifelse(Model %in% c("PAR1006", "PAR1177", "PAR1141", "PAR1221"), "CTO", "CT"))
  )

# Save emmeans summary to CSV
write.csv2(emm_pairs_summary, file = "Trial2_metric2_stats.csv", row.names = FALSE)

# Plot AUC differences by Model and Drug
pdf("Trial2_metric2.pdf", width = 12, height = 9)
ggplot(emm_pairs_summary, aes(x = fct_reorder(Model, estimate, .desc = TRUE), y = estimate, fill = pCR)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  scale_fill_manual(values = c("No" = "#0E417F", "Yes" = "#E5AA27")) +
  labs(x = "Model - Drug", y = "AUC T ~ UT", fill = "pCR", title = "PDTX Models") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(~ Drug, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon test for pCR and non-pCR groups (Metric 2)
wilcox_test <- wilcox.test(estimate ~ pCR, data = emm_pairs_summary)

pdf("Trial2_metric2_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(emm_pairs_summary, aes(x = pCR, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Comparison of Estimates between pCR and non-pCR groups",
       x = "pCR Status", y = "Estimate Differences in AUC") +
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
  annotate("text", x = 1.5, y = max(emm_pairs_summary$estimate) * 1.1, 
           label = paste("Wilcoxon Test:\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# RCB Score plots and analysis for Metric 2
emm_pairs_summary <- emm_pairs_summary %>%
  mutate(
    RCB_Score = case_when(
      Model == "PAR1006" ~ 0, Model == "PAR1008" ~ 1, Model == "PAR1022" ~ 0,
      Model == "PAR1040" ~ 2, Model == "PAR1053" ~ 0, Model == "PAR1141" ~ 2,
      Model == "PAR1177" ~ 0, Model == "PAR1221" ~ 2
    ),
    RCB_Category = factor(case_when(
      RCB_Score == 0 ~ "pCR", RCB_Score == 1 ~ "RCB-I", RCB_Score == 2 ~ "RCB-II"
    ), levels = c("pCR", "RCB-I", "RCB-II"))
  )

pdf("Trial2_metric2_RCB.pdf", width = 11, height = 5)
ggplot(emm_pairs_summary, aes(x = fct_reorder(Model, estimate, .desc = TRUE), y = estimate, fill = Drug)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  scale_fill_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  labs(x = "Model", y = "Estimated Difference", fill = "Drug", title = "PDTX Models by Drug and RCB Category") +
  expand_limits(y = c(min(emm_pairs_summary$estimate), max(emm_pairs_summary$estimate))) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere-Terpstra Test for Metric 2 (Trial 2)
emm_pairs_summary$RCB_Category <- factor(emm_pairs_summary$RCB_Category, 
                                         levels = c("pCR", "RCB-I", "RCB-II"), 
                                         ordered = TRUE)

jt_result <- JonckheereTerpstraTest(estimate ~ RCB_Category, data = emm_pairs_summary, 
                                    alternative = "decreasing")

summary_data <- emm_pairs_summary %>%
  group_by(RCB_Category) %>%
  summarise(mean_estimate = mean(estimate, na.rm = TRUE), .groups = 'drop')

pdf("Trial2_metric2_RCB_JT.pdf", width = 8, height = 6)
ggplot(emm_pairs_summary, aes(x = RCB_Category, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_estimate, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_estimate), size = 3, color = "darkblue") +
  labs(title = "Estimated Difference by RCB Category",
       x = "RCB Category", y = "Estimated Difference in AUC") +
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
  annotate("text", x = 2, y = max(emm_pairs_summary$estimate) * 1.1, 
           label = paste("Jonckheere-Terpstra Test Results:\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# Load necessary libraries
library(dplyr)
library(lmerTest)
library(ggplot2)
library(forcats)
library(DescTools)

# Initialize results_summary dataframe
results_summary <- data.frame(Model = character(), 
                              Lb_hat = numeric(),
                              CI_low = numeric(), 
                              CI_high = numeric(),
                              stringsAsFactors = FALSE)

alpha <- 0.05  # Confidence level

# Iterate over each unique model in the dataset
for (model_name in unique(PDTX_TP2nz$Model)) {
  # Subset data for the current model
  current_data <- filter(PDTX_TP2nz, Model == model_name)
  current_data$lSize <- log(current_data$Size)
  
  # Fit the mixed-effects model
  m_t <- lmer(lSize ~ Time * Group + (Time | Mouse), data = current_data, 
              control = lmerControl(optimizer = "bobyqa"))
  
  # Dynamically set the K vector based on the presence of GroupCT or GroupCTO
  fixed_effects <- names(fixef(m_t))
  K <- c("(Intercept)" = 0, "Time" = 1, "GroupCT" = 0, "GroupCTO" = 0) #
  interaction_term <- intersect(fixed_effects, c("Time:GroupCT", "Time:GroupCTO"))
  if (length(interaction_term) > 0) {
    K[interaction_term] <- 1
  }
  
  # Calculate the estimates and confidence intervals
  Va <- vcov(m_t)
  b.hat <- fixef(m_t)
  Lb.hat <- sum(K * b.hat)
  SE <- sqrt(sum((K %*% Va) * K))
  ddf <- get_Lb_ddf(m_t, K) 
  confint <- Lb.hat + c(-1, 1) * SE * qt(1 - 0.5 * alpha, ddf)
  
  results_summary <- rbind(results_summary, data.frame(Model = model_name, Lb_hat = Lb.hat, CI_low = confint[1], CI_high = confint[2]))
}

# Add pCR and Drug classification based on Model
results_summary <- results_summary %>%
  mutate(
    pCR = factor(ifelse(Model %in% c("PAR1040", "PAR1008", "PAR1141", "PAR1221"), "No", "Yes")),
    Drug = factor(ifelse(Model %in% c("PAR1006", "PAR1177", "PAR1141", "PAR1221"), "CTO", "CT"))
  )

# Save the results to CSV
write.csv2(results_summary, file = "Trial2_metric3_stats.csv", row.names = FALSE)

# Plot Growth Rate under Treatment by Model and Drug
pdf("Trial2_metric3.pdf", width = 12, height = 9)
ggplot(results_summary, aes(x = fct_reorder(Model, Lb_hat, .desc = FALSE), y = Lb_hat, color = pCR)) +
  geom_point(show.legend = FALSE, size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5), show.legend = TRUE) +
  scale_color_manual(values = c("No" = "#0E417F", "Yes" = "#E5AA27")) +
  labs(x = "Model", y = "Growth rate under treatment", color = "pCR", title = "PDTX Models") +
  expand_limits(y = c(-0.2, 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(~ Drug, nrow = 1, scales = "free_x")
dev.off()

# Plot Growth Rate by pCR status
results_summary$pCR <- factor(results_summary$pCR, levels = c("Yes", "No"), labels = c("pCR", "non-pCR"))

pdf("Trial2_metric3_pCR.pdf", width = 11, height = 5)
ggplot(results_summary, aes(x = fct_reorder(Model, Lb_hat, .desc = FALSE), y = Lb_hat, color = Drug)) +
  geom_point(show.legend = TRUE, size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  labs(x = "Model", y = "Growth rate under treatment", color = "Drug", title = "PDTX Models by Drug") +
  expand_limits(y = c(-0.2, 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  facet_wrap(~ pCR, nrow = 1, scales = "free_x")
dev.off()

# Wilcoxon test for pCR and non-pCR groups
wilcox_test <- wilcox.test(Lb_hat ~ pCR, data = results_summary)

pdf("Trial2_metric3_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(results_summary, aes(x = pCR, y = Lb_hat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(
    title = "Comparison between pCR and non-pCR groups",
    x = "pCR Status",
    y = "Growth Rate under treatment"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.5, y = max(results_summary$CI_high) * 1.1, 
           label = paste("Wilcoxon Test:", 
                         "\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 7, size = 6)
dev.off()

# RCB Score plots and analysis for Metric 3
results_summary <- results_summary %>%
  mutate(
    RCB_Score = case_when(
      Model == "PAR1006" ~ 0, Model == "PAR1008" ~ 1, Model == "PAR1022" ~ 0,
      Model == "PAR1040" ~ 2, Model == "PAR1053" ~ 0, Model == "PAR1141" ~ 2,
      Model == "PAR1177" ~ 0, Model == "PAR1221" ~ 2
    ),
    RCB_Category = factor(case_when(
      RCB_Score == 0 ~ "pCR", RCB_Score == 1 ~ "RCB-I", RCB_Score == 2 ~ "RCB-II"
    ), levels = c("pCR", "RCB-I", "RCB-II"))
  )

# Plot 1: Growth Rate Under Treatment by RCB Category
pdf("Trial2_metric3_RCB.pdf", width = 11, height = 5)
ggplot(results_summary, aes(x = fct_reorder(Model, Lb_hat, .desc = FALSE), y = Lb_hat, color = Drug)) +
  geom_point(show.legend = TRUE, size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = curve_colors) +
  labs(x = "Model", y = "Growth rate under treatment", color = "Drug", 
       title = "PDTX Models - Growth Rate Under Treatment by RCB Category") +
  expand_limits(y = c(-0.2, max(results_summary$CI_high) * 1.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Jonckheere-Terpstra Test for Metric 3
jt_result <- JonckheereTerpstraTest(Lb_hat ~ RCB_Category, data = results_summary, 
                                    alternative = "increasing")

# Plot with Jonckheere-Terpstra Test results
summary_data <- results_summary %>%
  group_by(RCB_Category) %>%
  summarise(mean_Lb_hat = mean(Lb_hat, na.rm = TRUE), .groups = 'drop')

pdf("Trial2_metric3_RCB_JT.pdf", width = 8, height = 6)
ggplot(results_summary, aes(x = RCB_Category, y = Lb_hat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drug, color = Drug), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_Lb_hat, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_Lb_hat), size = 3, color = "darkblue") +
  labs(title = "Growth Rate under Treatment by RCB Category",
       x = "RCB Category",
       y = "Growth rate under treatment") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_color_manual(values = curve_colors) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(results_summary$Lb_hat) * 1.1, 
           label = paste("Jonckheere-Terpstra Test Results:",
                         "\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = 7)
dev.off()

# Load necessary libraries
library(dplyr)
library(lmerTest)
library(ggplot2)
library(forcats)
library(DescTools)

# Initialize lists to store models and data
model_list <- list()
combined_data_list <- list()

# Iterate over each unique model in the dataset
for (model_name in unique(PDTX_TP2nz$Model)) {
  tumrs <- subset(PDTX_TP2nz, Model == model_name)
  tumrs$lSize <- log(tumrs$Size)
  tumrs$lSize[which(!is.finite(tumrs$lSize))] <- NA
  tumrs$Group <- factor(tumrs$Group)
  tumrs$Group <- relevel(tumrs$Group, ref = "UT")
  
  # Fit the mixed-effects model
  m_t <- lmerTest::lmer(lSize ~ Time * Group + (Time | Mouse), data = tumrs, 
                        control = lmerControl(optimizer = "bobyqa"))
  model_list[[model_name]] <- m_t 
  
  # Create a prediction frame
  dc <- distinct(dplyr::select(tumrs, Mouse, Group))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 84) %>% 
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  # Combine actual data with predicted data
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Model <- model_name
  
  combined_data_list[[model_name]] <- comb
}

# Combine all models' data
final_combined_data <- bind_rows(combined_data_list)

# Filter data for Time = 84 and Group CT or CTO
predicted_final_bxp <- final_combined_data %>%
  filter(Time == 84 & Group %in% c("CT", "CTO")) %>%
  distinct(Model, Mouse, .keep_all = TRUE) %>%
  rename(Drugs = Group, Size = ModelSize)

# Assign pCR status per model
predicted_final_bxp$pCR <- "Yes"
predicted_final_bxp[predicted_final_bxp$Model %in% c("PAR1040", "PAR1008", "PAR1141", "PAR1221"), ]$pCR <- "No"
predicted_final_bxp$pCR <- factor(predicted_final_bxp$pCR)

# Save the results to CSV
write.csv2(predicted_final_bxp, file = "Trial2_metric4_stats.csv", row.names = FALSE)

# Plot 1: Boxplots for predicted values at Time 84
pdf("Trial2_metric4.pdf", width = 12, height = 9)
ggplot(predicted_final_bxp, aes(x = fct_reorder(Model, Size, .fun = median, .desc = TRUE), y = Size)) +
  geom_boxplot(aes(fill = pCR), position = position_dodge(0.5), show.legend = TRUE) +
  scale_fill_manual(values = c("No" = "#0E417F", "Yes" = "#E5AA27")) +
  labs(x = "Model - Drug", y = "Predicted volume at the end of treatment", fill = "pCR", title = "PDTX Models Boxplots") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(~ Drugs, nrow = 1, scales = "free_x")
dev.off()

# Plot 2: Boxplots by pCR status
pdf("Trial2_metric4_pCR.pdf", width = 11, height = 5)
ggplot(predicted_final_bxp, aes(x = fct_reorder(Model, Size, .fun = median, .desc = TRUE), y = Size, fill = Drugs)) +
  geom_boxplot(position = position_dodge(0.5), show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  labs(x = "Model", y = "Predicted volume at the end of treatment (Day 84)", fill = "Drug", title = "PDTX Models by Drugs") +
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
wilcox_test <- wilcox.test(Size ~ pCR, data = predicted_final_bxp)

pdf("Trial2_metric4_pCR_wilcox.pdf", width = 8, height = 6)
ggplot(predicted_final_bxp, aes(x = pCR, y = Size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  labs(title = "Predicted Tumour Sizes between pCR and non-pCR groups",
       x = "pCR Status",
       y = "Predicted Tumour Sizes") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +
  annotate("text", x = 1.5, y = 500 * 0.9, 
           label = paste("Wilcoxon Test:",
                         "\nW =", round(wilcox_test$statistic, 2), 
                         "\np-value =", format.pval(wilcox_test$p.value, digits = 4)),
           hjust = 0.5, vjust = 1)
dev.off()

# Add RCB_Score and RCB_Category columns
predicted_final_bxp <- predicted_final_bxp %>%
  mutate(RCB_Score = case_when(
    Model == "PAR1006" ~ 0,
    Model == "PAR1008" ~ 1,
    Model == "PAR1022" ~ 0,
    Model == "PAR1040" ~ 2,
    Model == "PAR1053" ~ 0,
    Model == "PAR1141" ~ 2,
    Model == "PAR1177" ~ 0,
    Model == "PAR1221" ~ 2
  ),
  RCB_Category = factor(case_when(
    RCB_Score == 0 ~ "pCR",
    RCB_Score == 1 ~ "RCB-I",
    RCB_Score == 2 ~ "RCB-II"
  ), levels = c("pCR", "RCB-I", "RCB-II")))

# Plot 1: Boxplots by RCB Category
pdf("Trial2_metric4_RCB.pdf", width = 11, height = 5)
ggplot(predicted_final_bxp, aes(x = fct_reorder(Model, Size, .fun = median, .desc = TRUE), y = Size, fill = Drugs)) +
  geom_boxplot(position = position_dodge(0.5), show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = curve_colors) +
  labs(x = "Model", y = "Predicted volume at the end of treatment (Day 84)", 
       fill = "Drug", title = "PDTX Models by RCB Category") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) +
  facet_wrap(~ RCB_Category, nrow = 1, scales = "free_x")
dev.off()

# Plot 2: Boxplot by RCB Category for Metric 4
ggplot(predicted_final_bxp, aes(x = RCB_Category, y = Size)) +
  geom_boxplot() +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.2, size = 3) +
  scale_color_manual(values = curve_colors) +
  labs(x = "RCB Category", 
       y = "Predicted volume at the end of treatment (Day 84)", 
       title = "PDTX Models by RCB Category - Metric 4",
       shape = "Drug",
       color = "Drug") +
  expand_limits(y = c(0, max(predicted_final_bxp$Size) * 1.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold"))

# Perform Jonckheere-Terpstra test
jt_result <- JonckheereTerpstraTest(Size ~ RCB_Category, data = predicted_final_bxp, 
                                    alternative = "increasing")

summary_data <- predicted_final_bxp %>%
  group_by(RCB_Category) %>%
  summarise(mean_Size = mean(Size, na.rm = TRUE), .groups = 'drop')

# Plot Jonckheere-Terpstra test results
pdf("Trial2_metric4_RCB_JT.pdf", width = 8, height = 6)
ggplot(predicted_final_bxp, aes(x = RCB_Category, y = Size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Drugs, color = Drugs), width = 0.3, size = 3) +
  geom_line(data = summary_data, aes(y = mean_Size, group = 1), color = "darkblue") +
  geom_point(data = summary_data, aes(y = mean_Size), size = 3, color = "darkblue") +
  labs(title = "Predicted Tumour Size by RCB Category",
       x = "RCB Category",
       y = "Predicted Tumour Size") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = c("CT" = "#C43826", "CTO" = "#499A87")) +
  scale_shape_manual(values = c("CT" = 16, "CTO" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 2, y = max(predicted_final_bxp$Size, na.rm = TRUE) * 1.1, 
           label = paste("Jonckheere-Terpstra Test Results:",
                         "\nJT =", round(jt_result$statistic, 2), 
                         "\np-value =", format.pval(jt_result$p.value, digits = 4)),
           hjust = 0.5, vjust = -5.9)
dev.off()

