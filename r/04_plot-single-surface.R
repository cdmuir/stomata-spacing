source("r/header.R")

single_surface_results = readr::read_rds("objects/single_surface_results.rds") |>
  mutate(light_reorder = fct_relevel(light, "low", "medium", "high"))

# anova
single_surface_anova = aov(dispersion ~ light * surface, data = single_surface_results)
write_rds(single_surface_anova, "objects/single_surface_anova.rds")

stat.test <- aov(dispersion ~ light * surface, data = single_surface_results) |>
  tukey_hsd()
stat.test = stat.test[1:3,]

# anova
density_anova = aov(n_stomata ~ light * surface, data = single_surface_results)
write_rds(density_anova, "objects/density_anova.rds")

stat.test.density <- aov(n_stomata ~ light * surface, data = single_surface_results) |>
  tukey_hsd()
stat.test.density = stat.test.density[1:3,]

# 1. plot dispersion index ----
plot = ggplot(
    single_surface_results,
    aes(light_reorder, dispersion, fill = surface)
  ) +
  stat_boxplot() +
  xlab("Light treatment") +
  ylab("Dispersion Index\n0 = random; 1 = maximum") +
  scale_fill_grey(start = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif",
    y.position = c(.5, .575, .65),
    inherit.aes = FALSE
  ) +
  theme_pubr()

ggsave("single-surface.pdf", plot = plot, path = "ms/figures",
       width = 4, height = 4)

# 2. plot stomatal density ----
density_plot = ggplot(
    single_surface_results,
    aes(light_reorder, n_stomata, fill = surface)
  ) +
  stat_boxplot() +
  xlab("Light treatment") +
  ylab(expression(paste("Stomatal density [pores are", a^-2, "]"))) +
  scale_fill_grey(start = 0.3) +
  stat_pvalue_manual(
    stat.test.density, label = "p.adj.signif",
    y.position = c(100, 110, 50),
    inherit.aes = FALSE
  ) +
  theme_pubr()

ggsave("density.pdf", plot = density_plot, path = "ms/figures",
       width = 4, height = 4)
