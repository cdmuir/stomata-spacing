# I think I did this before, but calculating g_liq from Evans et al. 2009
# Measured range of r_liq = 25-300 m^2 chloroplast s bar / mol
# 1 s / m  = 0.025 m^2 chloroplast s bar / mol
# SO, 1000 - 12000 s / m
# SO g_liq = 0.001 - 8.33e-5

source("r/header.R")
source("r/photo_2d_pm.R")

library(furrr)
plan(multisession, workers = 10)

library(progressr)
handlers(global = TRUE)
handlers("progress", "beepr")

input = crossing(
  n_x = c(50, 100),
  n_z = c(250, 500),
  stomatal_arrangement = c("aligned", "offset"),
  I_0 = round(seq(0.0, 0.001, 0.0001), 4)
) %>%
  split(~ seq_len(nrow(.)))

run_photo_2d_pm = function(input, fit_rubisco) {

  p = progressor(along = input)
  out = future_map(input, \(.x, fit_rubisco) {

    r = with(.x, list(I_0 = I_0, n_x = n_x, n_z = n_z))
    ph2d = photo_2d_pm(.x$stomatal_arrangement, fit_rubisco = fit_rubisco, replace = r)
    ph2d$parms$stomatal_arrangement = .x$stomatal_arrangement
    ph2d$parms$A_n = calc_2d_pm_An(ph2d)
    p()
    ph2d

  }, fit_rubisco = fit_rubisco)

  out

}

tmp = run_photo_2d_pm(input, fit_rubisco)

tmp |>
  map("parms") |>
  map_dfr(\(x) as_tibble(x[c("stomatal_arrangement", "n_z", "n_x", "I_0", "A_n")])) |>
  ggplot(aes(I_0, A_n, color = stomatal_arrangement)) +
  facet_grid(n_z ~ n_x) +
  geom_point()


r = list(I_0 = 0.0003, n_x = 50, n_z = 500)
ph2d = photo_2d_pm("aligned", replace = r)

ph2d = photo_2d_pm("aligned", replace = list(I_0 = 0.001))

# Plot results
df_C = expand.grid(
  z = seq_len(ph2d$parms[["n_z"]]),
  x = seq_len(ph2d$parms[["n_x"]]),
  name = c("C_ias", "C_liq")
) |>
  mutate(value = ph2d$fit$y)

ggplot(df_C, aes(x, z, z = value)) +
  facet_grid(~ name) +
  geom_contour_filled() +
  coord_equal()

# Calculate area-based A_n

df_A = expand.grid(
  z = seq_len(parms[["n_z"]]),
  x = seq_len(parms[["n_x"]])
) |>
  mutate(a = c(a_n))

ggplot(df_A, aes(x, z, z = a)) +
  geom_contour_filled() +
  coord_equal()




# SIMPLE VERSION I SENT TO TOM
photo_2d_pm = function(t, Y, parms) {

  # Arguments:
  # * t, time step used by steady.2d() to find steady-state solution
  # * Y, vector of length 2 * n_x * n_z. The first half are the elements
  # corresponding to C_ias[i,j]; the second half are correspoding elements for
  # C_liq[i,j].
  # * parms, list of model parameters. Use `get_2d_pm_default_parms()` for
  #  default parameter values.

  # Set up:
  # The leaf is divided to an area n_x elements wide, n_z units deep
  n = parms[["n_z"]] * parms[["n_x"]]

  # Create empty matrices for computation
  C_ias = matrix(nrow = parms[["n_z"]], ncol = parms[["n_x"]],
                 data = Y[1:n])
  C_liq = matrix(nrow = parms[["n_z"]], ncol = parms[["n_x"]],
                 data = Y[(n + 1):(2 * n)])
  dC_ias = dC_liq = numeric(length = n) # empty vectors

  # FLUX ----
  D_e = parms[["D_c"]] * parms[["phi"]] / parms[["tau"]]

  # Boundary conditions
  bound_bottom = C_ias[1,]
  bound_top = C_ias[parms[["n_z"]],]
  bound_left = C_ias[,1]
  bound_right = C_ias[, parms[["n_x"]]]
  bound_top[1] = parms[["C_stom"]]
  bound_bottom[parms[["n_x"]]] = parms[["C_stom"]]
  # bound_bottom[1] = parms[["C_stom"]]

  # diffusion in Z-direction; boundaries=imposed concentration
  Flux = -D_e / parms[["t_elem"]] * rbind(
    C_ias[1, ] - bound_bottom,
    (C_ias[2:parms[["n_z"]], ] - C_ias[1:(parms[["n_z"]] - 1), ]),
    bound_top - C_ias[parms[["n_z"]], ]
  )
  dC_ias = dC_ias -
    (Flux[2:(parms[["n_z"]] + 1), ] - Flux[1:parms[["n_z"]], ]) / parms[["t_elem"]]

  # diffusion in X-direction
  Flux = -D_e / parms[["t_elem"]] * cbind(
    C_ias[, 1] - bound_left,
    (C_ias[, 2:parms[["n_x"]]] - C_ias[, 1:(parms[["n_x"]] - 1)]),
    bound_right - C_ias[, parms[["n_x"]]]
  )
  dC_ias = dC_ias -
    (Flux[, 2:(parms[["n_x"]] + 1)] - Flux[, 1:parms[["n_x"]]]) /
    parms[["t_elem"]]

  dC_ias = dC_ias + parms[["g_liq"]] * (C_liq - C_ias) /
    (parms[["T_leaf"]] / parms[["S_m"]])

  # CARBOXYLATION AND RESPIRATION ----
  # Calculate volumetric j_max from area-based J_max. Here, I assume a single
  # j_max for every part of the leaf, but in the final model I will have a
  # gradient of j_max following Earles et al. (2017).
  #
  j_max = parms[["J_max"]] / (parms[["S_m"]] * parms[["V_strom"]])

  # Carboxylation (n.b. Rubisco concentration is assumed constant throughout
  # the leaf, but in the final model I will have a garudent of X_c following
  # Earles et al. (2017))
  w_c = (parms[["k_c"]] * parms[["X_c"]] * C_liq) / (parms[["K_m"]] + C_liq)
  w_j = C_liq * j_max / (4 * C_liq + 8 * parms[["gamma_star"]])
  r_c = pmin(w_c, w_j)
  r_d = parms[["r_d"]]
  r_p = r_c * parms[["gamma_star"]] / C_liq

  # I'M NOT SURE THIS IS THE RIGHT WAY TO SCALE TO STROMA VOLUME
  # r_c, r_p, and r_d are mol CO2 per stroma volume. Need to convert to per leaf volume
  # (mol CO2 / s / m^3 stroma) * (m^3 stroma / m^2 mesophyll) * (m^2 mesophyll * m^2 leaf) * (m^2 leaf / m^3 leaf)
  # 2e4 assumes 200 um thick leaf because 1 m^2 has volume of 2e-4 m^3
  # S_m [m^2 meso / m^2 leaf] * (1 m^2 / (1 m x 1 m x t_leaf m) [m^2 leaf / m^3 leaf]
  dC_liq = dC_liq +
    parms[["g_liq"]] * (C_ias - C_liq) / (parms[["T_leaf"]] / parms[["S_m"]]) +
    (-r_c + r_p + r_d) * (parms[["S_m"]] / parms[["T_leaf"]]) * parms[["V_strom"]]

  return(list(c(dC_ias, dC_liq)))

}

parms = get_2d_pm_default_parms() |>
  derive_2d_pm_parms()

# Initial values
C_ias_mat = matrix(nrow = parms[["n_z"]], ncol = parms[["n_x"]], parms[["C_stom"]])
C_liq_mat = matrix(nrow = parms[["n_z"]], ncol = parms[["n_x"]], parms[["C_stom"]])

# Solve for C_ias and C_liq
soln = steady.2D(c(C_ias_mat, C_liq_mat), func = photo_2d_pm, parms = parms,
                 pos = FALSE, dimens = c(parms[["n_z"]], parms[["n_x"]]),
                 nspec = 2, lrw = 1e8, atol = 1e-10, rtol = 1e-10, ctol = 1e-10)

# Plot results
df_C = expand.grid(
  z = seq_len(parms[["n_z"]]),
  x = seq_len(parms[["n_x"]]),
  name = c("C_ias", "C_liq")
) |>
  mutate(value = soln$y)

ggplot(df_C, aes(x, z, z = value)) +
  facet_grid(~ name) +
  geom_contour_filled() +
  coord_equal()

# Calculate area-based net photosynthesis
C_liq = df_C |>
  filter(name == "C_liq") |>
  pull(value)

j_max = parms[["J_max"]] / (parms[["S_m"]] * parms[["V_strom"]])

w_c = (parms[["k_c"]] * parms[["X_c"]] * C_liq) / (parms[["K_m"]] + C_liq)
w_j = C_liq * j_max / (4 * C_liq + 8 * parms[["gamma_star"]])
r_c = pmin(w_c, w_j)
r_d = parms[["r_d"]]
r_p = r_c * parms[["gamma_star"]] / C_liq

T_leaf = parms[["n_z"]] * parms[["t_elem"]]
a_n = (r_c - r_p - r_d) * (parms[["S_m"]] / T_leaf) * parms[["V_strom"]]

# 1 m^2 of 200 um thick leaf is 2e-04 m^3
mean(a_n) * T_leaf * 1e6
