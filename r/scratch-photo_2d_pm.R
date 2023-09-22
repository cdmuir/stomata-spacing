source("r/header.R")

diffusion2D <- function(t, Y, par)   {

  y    <- matrix(nrow = n_z, ncol = n_x, data = Y)  # vector to 2-D matrix

  # hard-coding parameters for now
  S_m = 20
  Vstrom = 1.74e-6
  k_c = 3
  X_c = 2.5
  K_m = 18.7e-3
  Gamma = 1.75e-3

  r_c = (k_c * X_c * y) / (K_m + y)        # carboxylation
  r_d = 0.066
  r_p = r_c * Gamma / y

  # I'M NOT SURE THIS IS THE RIGHT WAY TO SCALE TO STROMA VOLUME
  # 2e4 assumes 200 um thick leaf because 1 m^2 has volume of 2e-4 m^3
  # S_m [m^2 meso / m^2 leaf] * 2e4 [m^2 leaf / m^3 leaf]
  dY   <- -(r_c + r_p + r_d) * (S_m * 2e4) * Vstrom

  bound_bottom <- y[1,] # boundary concentration
  bound_top <- y[n_z,]   # boundary concentration
  bound_left = y[,1]
  bound_right = y[, n_x]
  bound_top[1] = 0.015
  bound_bottom[n_x] = 0.015

  #diffusion in Z-direction; boundaries=imposed concentration
  Flux <- -Dz * rbind(
    y[1, ] - bound_bottom,
    (y[2:n_z, ] - y[1:(n_z - 1), ]),
    bound_top - y[n_z,]
  ) / dz
  dY   <- dY - (Flux[2:(n_z+1),]-Flux[1:n_z,])/dz

  #diffusion in X-direction
  Flux <- -Dx * cbind(
    y[, 1] - bound_left,
    (y[, 2:n_x] - y[, 1:(n_x - 1)]),
    bound_right - y[, n_x]
  ) / dx
  dY    <- dY - (Flux[,2:(n_x+1)]-Flux[,1:n_x])/dx

  return(list(as.vector(dY)))
}

# parameters
n_x = 100
n_z = 200
dz    <- dx <- 1e-6   # grid size
Dz    <- Dx <- 1.54e-5   # diffusion coeff, X- and Y-direction

C_ias_mat  <- matrix(nrow = n_z, ncol = n_x, 0.1)

# stodes is used, so we should specify the size of the work array (lrw)
# We take a rather large value
# By specifying pos = TRUE, only positive values are allowed.

system.time(
  ST2 <- steady.2D(C_ias_mat, func = diffusion2D, parms = NULL, pos = TRUE,
                   dimens = c(n_z, n_x), lrw = 100000000,
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
)

mean(ST2$y) # 0.009950982

expand.grid(z = seq_len(n_z), x = seq_len(n_x)) |>
  mutate(C_ias = ST2$y) |>
  ggplot(aes(x, z, z = C_ias)) +
  geom_contour_filled()



S_m = 20
Vstrom = 1.74e-6

k_c = 3
X_c = 2.5
K_m = 18.7e-3
Gamma = 1.75e-3

r_c = (k_c * X_c * ST2$y) / (K_m + ST2$y)        # carboxylation
r_d = 0.066
r_p = r_c * Gamma / ST2$y

a_n = (r_c - r_p - r_d) * (S_m * 2e4) * Vstrom

mean(a_n)

# 1 m^2 of 200 um thick leaf is 2e-04 m^3
#
# 2.0801 mol / m^3 / s * (1 m^2 / 2e-4 m^3) * (1e6 umol / mol)
mean(a_n) * 2e-4 * 1e6
