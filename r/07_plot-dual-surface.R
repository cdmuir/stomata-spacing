##### Plot Dual Surface Results #####
source("r/header.R")

dual_surface_results = readr::read_rds("objects/dual_surface_results.rds") |>
  # add light column
  dplyr::mutate(
    light = dplyr::case_when(
      str_detect(treatment, "HL") ~ "high",
      str_detect(treatment, "LL") ~ "low"
    ),
    light = tidyr::replace_na(light, "medium")) |>
  mutate(corr_diff = corr - mu_corr_sim) |>
  ungroup() |>
  mutate(light_reorder = fct_relevel(light, "low", "medium", "high"))

# anova
correlation_anova = aov(corr ~ light, data = dual_surface_results)
summary(correlation_anova)

stat.test.correlation <- aov(corr ~ light, data = dual_surface_results) |>
  tukey_hsd()
stat.test.correlation

plot = ggplot(dual_surface_results, aes(light_reorder, corr)) +
  stat_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Light treatment") +
  ylab(expression(paste("Correlation between NS", D^2))) +
  theme_pubr()


ggsave("dual-surface.pdf", plot = plot, path = "ms/figures",
       width = 3, height = 3)

