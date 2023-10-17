source("r/header.R")

fit = read_rds("objects/fit.rds")

emtrends(fit, specs = c("light", "surface"), var = "length") |>
  write_rds("objects/fit_emtrends.rds")

df_new = fit$data |>
  dplyr::group_by(light, surface) |>
  summarise(length_min = min(length), length_max = max(length)) |>
  rowwise() |>
  pmap_dfr(~{
    tibble(light = ..1, surface = ..2, length = seq(..3, ..4, length.out = 1e2))
  })

df = posterior_epred(fit, newdata = df_new, re_formula = NA)

df_med = df |>
  apply(2, median) |>
  as_tibble(.name_repair = "minimal") |>
  rename(sqrt_area = value)

df_qi = df |>
  apply(2, ggdist::qi) |>
  t() |>
  as_tibble(.name_repair = "minimal") |>
  magrittr::set_colnames(c("lower", "upper"))

df_pred = df_new |>
  bind_cols(df_med, df_qi) |>
  arrange(light)

plot = ggplot(
  fit$data,
  aes(length, sqrt_area, color = light, fill = light, linetype = surface)
) +
  geom_ribbon(data = df_pred, mapping = aes(ymin = lower, ymax = upper),
              alpha = 0.35) +
  geom_line(data = df_pred) +
  geom_point(alpha = 0.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  xlab("Stomata length [UNITS]") +
  ylab(expression(paste(sqrt(Tessellation~zone), " [UNITS]"))) +
  theme_pubr() +
  theme(legend.position = "right", panel.background = element_rect(fill = "grey"))

ggsave("length-area.pdf", plot = plot, path = "ms/figures",
       width = 5, height = 4)
