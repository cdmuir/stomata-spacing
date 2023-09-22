source("r/header.R")
library(tidyr)

diffusion2D <- function(t, Y, par)   {

  n = n_z * n_x
  C_ias = matrix(nrow = n_z, ncol = n_x, data = Y[1:n])
  C_liq = matrix(nrow = n_z, ncol = n_x, data = Y[(n + 1):(2 * n)])
  dC_ias = dC_liq = numeric(length = n) # empty vectors

  # hard-coding parameters for now
  D_c = 1.54e-5
  por = 0.2
  tort = 1.55
  S_m = 20
  V_strom = 1.74e-6
  k_c = 3
  X_c = 2.5
  K_m = 18.7e-3
  Gamma = 1.75e-3

  g_liq = 0.25e-3
  r_c = (k_c * X_c * C_liq) / (K_m + C_liq)        # carboxylation
  r_d = 0.066
  r_p = r_c * Gamma / C_liq

  D_e = D_c * por / tort

  # FLUX ----
  bound_bottom <- C_ias[1,] # boundary concentration
  bound_top <- C_ias[n_z,]   # boundary concentration
  bound_left = C_ias[,1]
  bound_right = C_ias[, n_x]
  bound_top[1] = 0.015
  bound_bottom[n_x] = 0.015

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

  dC_ias = dC_ias + g_liq * (C_liq - C_ias) / dx

  # CARBOXYLATION ----
  # I'M NOT SURE THIS IS THE RIGHT WAY TO SCALE TO STROMA VOLUME
  # 2e4 assumes 200 um thick leaf because 1 m^2 has volume of 2e-4 m^3
  # S_m [m^2 meso / m^2 leaf] * 2e4 [m^2 leaf / m^3 leaf]
  dC_liq = dC_liq + g_liq * (C_ias - C_liq) / dx + (-r_c + r_p + r_d) * (S_m * 2e4) * V_strom

  return(list(c(dC_ias, dC_liq)))

}

# parameters
n_x = 100
n_z = 200
dz    <- dx <- 1e-6   # grid size

C_ias_mat = matrix(nrow = n_z, ncol = n_x, 0.015)
C_liq_mat = matrix(nrow = n_z, ncol = n_x, 0.015)

# stodes is used, so we should specify the size of the work array (lrw)
# We take a rather large value
# By specifying pos = TRUE, only positive values are allowed.

system.time(
  ST2 <- steady.2D(c(C_ias_mat, C_liq_mat), func = diffusion2D, parms = NULL,
                   pos = TRUE, dimens = c(n_z, n_x), nspec = 2, lrw = 100000000,
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
)

df_C = expand.grid(z = seq_len(n_z), x = seq_len(n_x), name = c("C_ias", "C_liq")) |>
  mutate(value = ST2$y) |>
  pivot_wider()

ggplot(df_C, aes(x, z, z = C_liq)) +
  geom_contour_filled()

ggplot(df_C, aes(C_ias, C_liq)) +
  geom_point()

S_m = 20
V_strom = 1.74e-6

k_c = 3
X_c = 2.5
K_m = 18.7e-3
Gamma = 1.75e-3

r_c = (k_c * X_c * df_C$C_liq) / (K_m + df_C$C_liq)        # carboxylation
r_d = 0.066
r_p = r_c * Gamma / df_C$C_liq

a_n = (r_c - r_p - r_d) * (S_m * 2e4) * V_strom

mean(a_n)

# 1 m^2 of 200 um thick leaf is 2e-04 m^3
#
# 2.0801 mol / m^3 / s * (1 m^2 / 2e-4 m^3) * (1e6 umol / mol)
mean(a_n) * 2e-4 * 1e6
