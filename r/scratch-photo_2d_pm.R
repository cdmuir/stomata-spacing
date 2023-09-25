# I think I did this before, but calculating g_liq from Evans et al. 2009
# Measured range of r_liq = 25-300 m^2 chloroplast s bar / mol
# 1 s / m  = 0.025 m^2 chloroplast s bar / mol
# SO, 1000 - 12000 s / m
# SO g_liq = 0.001 - 8.33e-5
source("r/header.R")
source("r/photo_2d_pm.R")

parms = get_2d_pm_default_parms() |>
  derive_2d_pm_parms(fit = fit_rubisco)

parms$X_c_z

diffusion2D <- function(t, Y, par)   {

  n = n_z * n_x
  C_ias = matrix(nrow = n_z, ncol = n_x, data = Y[1:n])
  C_liq = matrix(nrow = n_z, ncol = n_x, data = Y[(n + 1):(2 * n)])
  dC_ias = dC_liq = numeric(length = n) # empty vectors

  # hard-coding parameters for now
  D_c = 1.54e-5
  phi = 0.2
  tau = 1.55
  S_m = 20
  V_strom = 1.74e-6
  k_c = 3
  X_c = 2.5
  K_m = 18.7e-3
  Gamma = 1.75e-3
  g_liq = 0.25e-3
  J_max = 0.000275
  j_max = J_max / (S_m * V_strom)

  w_c = (k_c * X_c * C_liq) / (K_m + C_liq)        # carboxylation
  w_j = C_liq * j_max / (4 * C_liq + 8 * Gamma)
  r_c = pmin(w_c, w_j)
  r_d = 0.066
  r_p = r_c * Gamma / C_liq

  D_e = D_c * phi / tau

  # FLUX ----
  bound_bottom <- C_ias[1,] # boundary concentration
  bound_top <- C_ias[n_z,]   # boundary concentration
  bound_left = C_ias[,1]
  bound_right = C_ias[, n_x]
  bound_top[1] = 0.015
  # bound_bottom[n_x] = 0.015
  bound_bottom[1] = 0.015

  # diffusion in Z-direction; boundaries=imposed concentration
  Flux <- -D_e * rbind(
    C_ias[1, ] - bound_bottom,
    (C_ias[2:n_z, ] - C_ias[1:(n_z - 1), ]),
    bound_top - C_ias[n_z,]
  ) / dz
  dC_ias = dC_ias - (Flux[2:(n_z + 1), ] - Flux[1:n_z, ]) / dz

  # diffusion in X-direction
  Flux <- -D_e * cbind(
    C_ias[, 1] - bound_left,
    (C_ias[, 2:n_x] - C_ias[, 1:(n_x - 1)]),
    bound_right - C_ias[, n_x]
  ) / dx
  dC_ias = dC_ias - (Flux[, 2:(n_x + 1)] - Flux[, 1:n_x]) / dx

  dC_ias = dC_ias + g_liq * (C_liq - C_ias) / V_strom

  # CARBOXYLATION ----
  # I *THINK* V_strom is the right thing to divide here because g_liq is
  #  conductance per m^2 chloroplast area and that needs to be converted to
  #  volume. Specifically, how diffusion from IAS to stroma changes
  #  concentration
  #
  # I'M NOT SURE THIS IS THE RIGHT WAY TO SCALE TO STROMA VOLUME
  # r_c, r_p, and r_d are mol CO2 per stroma volume. Need to convert to per leaf volume.
  # (mol CO2 / s / m^3 stroma) * (m^3 stroma / m^2 mesophyll) * (m^2 mesophyll * m^2 leaf) * (m^2 leaf / m^3 leaf)
  # 2e4 assumes 200 um thick leaf because 1 m^2 has volume of 2e-4 m^3
  # S_m [m^2 meso / m^2 leaf] * (1 m^2 / (1 m x 1 m x t_leaf m) [m^2 leaf / m^3 leaf]
  t_leaf = n_z * dz # [m]
  dC_liq = dC_liq + g_liq * (C_ias - C_liq) / V_strom + (-r_c + r_p + r_d) * (S_m * 1 / t_leaf) * V_strom

  return(list(c(dC_ias, dC_liq)))

}

# parameters
n_x = 200
n_z = 400
dz    <- dx <- 0.5e-6   # grid size

C_ias_mat = matrix(nrow = n_z, ncol = n_x, 0.015)
C_liq_mat = matrix(nrow = n_z, ncol = n_x, 0.015)

# stodes is used, so we should specify the size of the work array (lrw)
# We take a rather large value
# By specifying pos = TRUE, only positive values are allowed.

system.time(
  ST2 <- steady.2D(c(C_ias_mat, C_liq_mat), func = diffusion2D, parms = NULL,
                   pos = FALSE, dimens = c(n_z, n_x), nspec = 2, lrw = 1e8,
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
)

df_C = expand.grid(z = seq_len(n_z), x = seq_len(n_x), name = c("C_ias", "C_liq")) |>
  mutate(value = ST2$y) |>
  pivot_wider()

ggplot(df_C, aes(x, z, z = C_ias)) +
  geom_contour_filled()

ggplot(df_C, aes(C_ias, C_liq)) +
  geom_point()

S_m = 20
V_strom = 1.74e-6

k_c = 3
X_c = 2.5
K_m = 18.7e-3
Gamma = 1.75e-3

J_max = 0.000275
j_max = J_max / (S_m * V_strom)

w_c = (k_c * X_c * df_C$C_liq) / (K_m + df_C$C_liq)        # carboxylation
w_j = df_C$C_liq * j_max / (4 * df_C$C_liq + 8 * Gamma)
r_c = pmin(w_c, w_j)
r_d = 0.066
r_p = r_c * Gamma / df_C$C_liq

t_leaf = n_z * dz # [m]
a_n = (r_c - r_p - r_d) * (S_m * 1 / t_leaf) * V_strom

mean(a_n) # 0.194078

# 1 m^2 of 200 um thick leaf is 2e-04 m^3
mean(a_n) * t_leaf * 1e6 # OK...maybe this is right
# 38.8156 (offset) vs. 38.80798 (overtop)
