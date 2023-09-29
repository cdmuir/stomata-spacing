# Model photosynthesis using a 2D porous medium approximation
#
# The concept is that this file will contain all functions to run 2D porous
# medium model and can be transitioned to a package script once it's ready

# Note that n_z = 0 is bottom (abaxial) side of leaf

# 2D photosynthesis model using a porous medium approximation
photo_2d_pm = function(stomatal_arrangment, replace = NULL, ...) {

  # Set parameters
  parms = set_2d_pm_parms(replace) |>
    derive_2d_pm_parms(fit_rubisco = fit_rubisco) |>
    make_2d_pm_mat() |>
    add_stomatal_position(stomatal_arrangment)

  # Find solution
  C_init = rep(parms[["C_stom"]], 2 * parms[["n_z"]] * parms[["n_x"]])
  fit = steady.2D(C_init, func = rd_2p_pm, parms = parms,
                  pos = FALSE, dimens = c(parms[["n_z"]], parms[["n_x"]]),
                  nspec = 2, lrw = 1e8, atol = 1e-10, rtol = 1e-10, ctol = 1e-10)

  # Structure output
  list(fit = fit, parms = parms)

}

# Function to add stomatal position to parms.
# For this project, it makes an amphistomatous leaf with stomata aligned or
# offset, but could be expanded for other purposes
add_stomatal_position = function(parms, stomatal_arrangment) {

  parms[["x_abaxial_stomata"]] = 1
  parms[["x_adaxial_stomata"]] = switch(
    stomatal_arrangment,
    aligned = 1,
    offset = parms[["n_x"]]
  )
  parms

}

# Set parameters for photo_2d_pm()
set_2d_pm_parms = function(replace, ...) {

  parms = get_2d_pm_default_parms()
  if (!is.null(replace)) {
    parms[names(replace)] = replace[names(replace)]
  }
  parms

}

# Get set of parameters for photo_2d_pm()
get_2d_pm_default_parms = function(...) {

  read_csv("raw-data/2d-pm-parameters.csv", col_types = "cccdccl") |>
    # read_csv("../raw-data/2d-pm-parameters.csv", col_types = "cccdccl") |>
    filter(!calculated) |>
    select(r, default_value) |>
    with(split(default_value, r))

}

# Calculate derived parameters for photo_2d_pm()
derive_2d_pm_parms = function(parms, fit_rubisco, ...) {

  # derived quantities
  parms$T_leaf = parms[["n_z"]] * parms[["t_elem"]] # leaf thickness [m]
  parms$U = parms[["n_x"]] * parms[["t_elem"]] # interstomatal distance [m]
  parms$f_spg = 1 - parms[["f_pal"]] # Fraction spongy mesophyll
  parms$n_z_pal = round(parms[["n_z"]] * parms[["f_pal"]])
  parms$n_z_spg = parms[["n_z"]] - parms[["n_z_pal"]]
  parms$S_m_z = rep(c(parms[["S_m_spg"]], parms[["S_m_pal"]]),
                 c(parms[["n_z_spg"]], parms[["n_z_pal"]]))

  parms$phi_z = rep(c(parms[["phi_spg"]], parms[["phi_pal"]]),
                  c(parms[["n_z_spg"]], parms[["n_z_pal"]]))

  parms %<>%
    derive_2d_pm_rubisco(fit = fit_rubisco) %>%
    derive_2d_pm_chl() %>%
    derive_2d_pm_j()

  parms

}

# Derive intraleaf Rubisco gradient for photo_2d_pm()
#
# fit: fitted model object that will predict relative Rubisco concentration as
# a function of relative depth (the explanatory variable needs to be called
# `rel_depth` in the model formula). It must work with the predict() function.
#
# integral01: the integral of fit from 0 to 1. Default is NULL, in which case
# the integral will be calculated. The integral is used to scale the Rubisco
# concentration in each element by the bulk leaf Rubisco concentration.
#
# ... : ignored, for extensibility
#
derive_2d_pm_rubisco = function(parms, fit, integral01 = NULL, ...) {

  rel_depth = seq(0 + 1 / parms[["n_z"]] / 2, 1 - 1 / parms[["n_z"]] / 2,
                  length.out = parms[["n_z"]])
  rel_rubisco = predict(fit, newdata = tibble(rel_depth = rel_depth))
  if (is.null(integral01)) {
    integral01 = integrate(\(x, fit) predict(fit, newdata = tibble(rel_depth = x)),
                           lower = 0, upper = 1, fit = fit)[["value"]]
  }
  parms$X_c_z = rel_rubisco * parms[["X_c"]] / integral01

  parms

}

# Derive intraleaf chlorophyll gradient for photo_2d_pm()
derive_2d_pm_chl = function(parms, ...) {

  # Using quadratic equation proposed by Johnson et al. (2005), see also
  # Borsuk and Brodersen (2019)
  rel_depth = seq(0 + 1 / parms[["n_z"]] / 2, 1 - 1 / parms[["n_z"]] / 2,
                  length.out = parms[["n_z"]])

  # analytical solution for integral of quadratic from 0 to 1
  integral01 = parms[["b0_chl"]] + parms[["b1_chl"]] / 2 +
    parms[["b2_chl"]] / 3

  # **NEED TO REVERSE BECAUSE BORSUK AND BRODERSEN HAVE 0 AT ADAXIAL SURFACE**
  parms$F_chl_z = rev((parms[["b0_chl"]] + parms[["b1_chl"]] * rel_depth +
    parms[["b2_chl"]] * rel_depth ^ 2) / integral01)

  parms

}

# Derive intraleaf maximum electron transport rate (j_max) gradient for photo_2d_pm()
# Assume volumetric electron transport rate is proportional Rubisco concentration
derive_2d_pm_jmax = function(parms, ...) {

  j_max_z1 = parms[["n_z"]] * parms[["X_c_z"]] / sum(parms[["X_c_z"]])
  mu_j_max_z1 = mean(j_max_z1 * parms[["S_m_z"]] * parms[["V_strom"]])
  jscale = parms[["J_max"]] / mu_j_max_z1
  parms$j_max_z = jscale * j_max_z1
  # Check
  checkmate::assert_true(abs(parms[["J_max"]] - mean(parms$j_max_z * parms[["S_m_z"]] * parms[["V_strom"]])) < 1e-6)

  parms

}

# Derive intraleaf potential electron transport rate (j_inf) gradient for photo_2d_pm()
# Assume potential electron transport rate is proportional Chl concentration
derive_2d_pm_jinf = function(parms, ...) {

  J_inf = parms[["I_0"]] * parms[["beta"]] * parms[["alpha"]] * parms[["phi_PSII"]]
  j_inf_z1 = parms[["n_z"]] * parms[["F_chl_z"]] / sum(parms[["F_chl_z"]])
  mu_j_inf_z1 = mean(j_inf_z1 * parms[["S_m_z"]] * parms[["V_strom"]])
  jscale = J_inf / mu_j_inf_z1
  parms$j_inf_z = jscale * j_inf_z1
  # Check
  checkmate::assert_true(abs(J_inf - mean(parms$j_inf_z * parms[["S_m_z"]] * parms[["V_strom"]])) < 1e-6)

  parms

}

# Derive intraleaf electron transport rate (j) gradient for photo_2d_pm()
derive_2d_pm_j = function(parms, ...) {

  parms %<>%
    derive_2d_pm_jmax() %>%
    derive_2d_pm_jinf()

  parms$j_e_z = pmin(parms[["j_inf_z"]], parms[["j_max_z"]])

  parms

}

# Make matrices of porosity, Sm, X_c, j_e for photo_2d_pm()
make_2d_pm_mat = function(parms, ...) {

  parms$phi_mat = matrix(parms[["phi_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$S_m_mat = matrix(parms[["S_m_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$X_c_mat = matrix(parms[["X_c_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$j_e_mat = matrix(parms[["j_e_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])

  parms

}

calc_2d_pm_rc = function(C_liq_mat, parms, ...) {

  w_c = (parms[["k_c"]] * parms[["X_c_mat"]] * C_liq_mat) / (parms[["K_m"]] + C_liq_mat)
  w_j = C_liq_mat * parms[["j_e_mat"]] / (4 * C_liq_mat + 8 * parms[["gamma_star"]])

  r_c = pmin(w_c, w_j)

  r_c

}

# Flux calculation for photo_2d_pm()
calc_2d_pm_flux = function(C_ias_mat, parms, ...) {

  dC_ias = numeric(length = length(C_ias_mat)) # empty vector

  # Boundary conditions
  ## No flux from top or bottom except where stomata present
  bnd_z_bottom = C_ias_mat[1, ]
  bnd_z_top = C_ias_mat[parms[["n_z"]], ]
  bnd_z_bottom[parms[["x_abaxial_stomata"]]] = parms[["C_stom"]]
  bnd_z_top[parms[["x_adaxial_stomata"]]] = parms[["C_stom"]]

  ## No flux from left or right because of symmetry condition
  bnd_x_left = C_ias_mat[, 1]
  bnd_x_right = C_ias_mat[, parms[["n_x"]]]

  D_e = parms[["D_c"]] * parms[["phi_mat"]] / parms[["tau"]]

  # diffusion in Z-direction (without D_e)
  flux_z = -rbind(
    C_ias_mat[1, ] - bnd_z_bottom,
    (C_ias_mat[2:parms[["n_z"]], ] - C_ias_mat[1:(parms[["n_z"]] - 1), ]),
    bnd_z_top - C_ias_mat[parms[["n_z"]], ]
  ) / parms[["t_elem"]]
  dC_ias = dC_ias -
    (flux_z[2:(parms[["n_z"]] + 1), ] - flux_z[1:parms[["n_z"]], ]) /
    parms[["t_elem"]]

  # diffusion in X-direction
  flux_x = -cbind(
    C_ias_mat[, 1] - bnd_x_left,
    (C_ias_mat[, 2:parms[["n_x"]]] - C_ias_mat[, 1:(parms[["n_x"]] - 1)]),
    bnd_x_right - C_ias_mat[, parms[["n_x"]]]
  ) / parms[["t_elem"]]
  dC_ias = dC_ias -
    (flux_x[, 2:(parms[["n_x"]] + 1)] - flux_x[, 1:parms[["n_x"]]]) /
    parms[["t_elem"]]

  # return
  D_e * dC_ias

}

# Reaction-diffusion function for steady.2D()
# Arguments:
# * t, time step used by steady.2d() to find steady-state solution
# * y, vector of length 2 * n_x * n_z. The first half are the elements
# corresponding to C_ias[i,j]; the second half are correspoding elements for
# C_liq[i,j].
# * parms, list of model parameters. Use `get_2d_pm_default_parms()` for
#  default parameter values.
# * ... Ignored. For extensibility

rd_2p_pm = function(t, y, parms, ...) {

  # Set up:
  # The leaf is divided to an area n_x elements wide, n_z units deep
  n = parms[["n_z"]] * parms[["n_x"]]

  # Create empty matrices for computation
  C_ias_mat = matrix(nrow = parms[["n_z"]],
                     ncol = parms[["n_x"]],
                     data = y[1:n])
  C_liq_mat = matrix(nrow = parms[["n_z"]],
                     ncol = parms[["n_x"]],
                     data = y[(n + 1):(2 * n)])

  # Calculate flux:
  dC_ias = calc_2d_pm_flux(C_ias_mat, parms)

  # Calculate carboxylation and respiration:
  r_c = calc_2d_pm_rc(C_liq_mat, parms)
  r_d = parms[["r_d"]]
  r_p = r_c * parms[["gamma_star"]] / C_liq_mat

  # Mass balance from IAS to stroma
  dC_ias = dC_ias + parms[["g_liq"]] * (C_liq_mat - C_ias_mat) /
    (parms[["T_leaf"]] / parms[["S_m_mat"]])

  dC_liq = parms[["g_liq"]] * (C_ias_mat - C_liq_mat) / (parms[["T_leaf"]] / parms[["S_m_mat"]]) +
    (-r_c + r_p + r_d) * (parms[["S_m_mat"]] / parms[["T_leaf"]]) * parms[["V_strom"]]

  return(list(c(dC_ias, dC_liq))) # return output from function

}

# Calculate A_n per area from output
calc_2d_pm_An = function(ph2d) {

  # ph2d is object returned from photo_2d_pm()
  df_C = expand.grid(
    z = seq_len(ph2d$parms[["n_z"]]),
    x = seq_len(ph2d$parms[["n_x"]]),
    name = c("C_ias", "C_liq")
  ) |>
    mutate(value = ph2d$fit$y)

  C_liq_mat = matrix(
    ph2d$fit$y[(ph2d$parms$n_x * ph2d$parms$n_z + 1):(2 * ph2d$parms$n_x * ph2d$parms$n_z)],
    nrow = ph2d$parms$n_z, ncol = ph2d$parms$n_x
  )
  r_c = calc_2d_pm_rc(C_liq_mat, ph2d$parms)
  r_d = ph2d$parms[["r_d"]]
  r_p = r_c * ph2d$parms[["gamma_star"]] / C_liq_mat

  a_n = (r_c - r_p - r_d) * (ph2d$parms[["S_m_mat"]] / ph2d$parms[["T_leaf"]]) * ph2d$parms[["V_strom"]]

  mean(a_n) * ph2d$parms[["T_leaf"]] * 1e6

}
