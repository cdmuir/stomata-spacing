source("r/header.R")
source("r/photo_2d_pm.R")

ph2d = read_rds("objects/model_output.rds")

# STUFF for processing
vars = c("stomatal_arrangement", "n_z", "n_x", "I_0", "A_n", "phi_pal")

ph2d_vars = ph2d |>
  map("parms") |>
  map_dfr(\(.x) as_tibble(.x[vars]))

ph2d_vars |>
  pivot_wider(names_from = "stomatal_arrangement", values_from = "A_n") |>
  # coordination advantage
  mutate(ca = log(offset / aligned)) |>
  ggplot(aes(I_0, ca, color = phi_pal)) +
  facet_grid(n_z ~ n_x) +
  geom_point()

ph2d_vars |>
  filter(n_z == 101, n_x == 169)
