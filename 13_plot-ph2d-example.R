source("r/header.R")
source("r/photo_2d_pm.R")

r = list(I_0 = 0.001, n_x = 101, n_z = 201, phi_pal = 0.2)
ph2d_offset = photo_2d_pm("offset", fit_rubisco = fit_rubisco, replace = r)
ph2d_aligned = photo_2d_pm("aligned", fit_rubisco = fit_rubisco, replace = r)

# Plot results
df_C = expand.grid(
  z = seq_len(ph2d_offset$parms[["n_z"]]),
  x = seq_len(ph2d_offset$parms[["n_x"]]),
  name = c("italic(C)[ias]", "italic(C)[liq]"),
  arrangement = c("offset", "aligned")
) |>
  mutate(value = c(ph2d_offset$fit$y, ph2d_aligned$fit$y))

ggplot(df_C, aes(x, z, z = value)) +
  facet_grid(arrangement ~ name, labeller = label_parsed) +
  # geom_contour_filled() +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  scale_y_continuous(breaks = c(0, 100, 200)) +
  geom_contour_filled(breaks = seq(0.004, 0.015, 0.0005)) +
  xlab(expression(distance~from~abaxial~stomate~bgroup("[", paste(mu, 'm'), "]"))) +
  ylab(expression(distance~from~abaxial~epidermis~bgroup("[", paste(mu, 'm'), "]"))) +
  coord_equal() +
  guides(fill = guide_legend(expression(group('[', CO[2], ']')~mol~m^-3))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ms/figures/ph2d-example.pdf", width = 6.5, height = 4.5)
write_rds(ph2d_offset, "objects/ph2d_offset.rds")
