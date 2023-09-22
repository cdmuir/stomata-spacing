# Model photosynthesis using a 2D porous medium approximation
#
# The concept is that this file will contain all functions to run 2D porous
# medium model and can be transitioned to a package script once it's ready

# Note that n_z = 0 is bottom (abaxial) side of leaf

# 2D photosynthesis model using a porous medium approximation
photo_2d_pm = function(...) {

  # Set parmseters
  parmss = set_2d_pm_parms(...)

  # Find solution
  # Structure output

}

# Set parmseters for photo_2d_pm()
set_2d_pm_parms = function(...) {

}

# Get set of parmseters for photo_2d_pm()
get_2d_pm_default_parms = function(...) {

  # Stomta spaced 100 um apart, 200 um thick leaf
  # With hexagonal grid, this is stomatal density of:
  # 2 / sqrt(3) / 100 ^ 2 = 0.0001154701 um^-2 = 115.4701 mm^-2
  list(
    n_x = 100,
    n_z = 200,
    t_elem = 1e-6,    # thickness of element [m]
    D_c = 1.54e-5,     # Diffusivity of CO2 in air [m^2 / s]
    por_int = 0.2,    # Porosity at interface [m^3 air / m^3 leaf]
    por_spg = 0.3,    # Porosity of the spongy mesophyll [m^3 air / m^3 leaf]
    por_pal = 0.1,    # Porosity of the palisade mesophyll [m^3 air / m^3 leaf]
    tort = 1.55,      # Tortuosity of the palisade and spongy mesophyll [m / m]

    # Note that 400 ppm = 0.01764706 is in mol / m^3, according to Earles et al. (2017)
    # Implies 44.11765 mol air / m^3 - find eq for this based on temp and P
    #
    # THERE IS SOME ISSUE WITH THE WAY I AM CODING THE BOUNDARY LAYER. I SET C_S = 0.21
    # AS A HACK TO GET RESULTS SIMILAR TO EARLES ET AL. THEY ASSUMED 0.015 (PG. 1093)
    C_s = 0.015,      # CO2 concentration at stomate [mol / m^3]

    k_c = 3,          # Rubisco turnover rate [1 / s]
    X_c = 2.5,        # Rubisco concentration [mol / m^3]
    K_m = 18.7e-3,    # Rubisco effective Michaelis-Menten constant [mol / m^3]
    Gamma = 1.75e-3,  # CO2 compensation point [mol / m^3]
    Theta = 1,        # Curvature factor in FvCB model [1]
    r_d = 0.066,       # Dark respiratory rate [mol / m^3 / s]
    Jmax = 275e-6,    # Maximum electron transport rate [mol / m^2 / s]
    Beta = 0.44,      # Fraction of light absorbed by PSII [mol / mol]
    Alpha = 0.72,     # Leaf level absorption [mol / mol]
    phiPSII = 0.85,   # quantum efficiency of PSII [1]
    I_0 = 2e-3,       # Incident irradiance [mol / m^2 / s]

    # I don't think I'm using this any more...
    # k_i so 90% absorption for 300 um
    # exp()
    # k_i = 7675.284,   # light extinction coefficient [1/m]

    # Quadratic parameters describing relative chlorophyll content as a function
    # of relative leaf depth (Johnson et al. [2005]; Borsuk and Brodersen [2019])
    b0_chl = 67.52,
    b1_chl = 100 * 0.4149, # converts parmseters to 0-1 rather than 0-100 scale
    b2_chl = 100 ^ 2 * -0.0029, # converts parmseters to 0-1 rather than 0-100 scale

    frac_pal = 0.6,   # Fraction palisade mesophyll [1]
    Sm_spg = 6.5,     # Sm spongy mesophyll [m^2 / m^2]
    Sm_pal = 40,      # Sm palisade mesophyll [m^2 / m^2]
    Vstrom = 1.74e-6, # Stroma volume per mesophyll surface area [m^3 / m^2]
    Vmito = 0.27e-7,  # Mitochondrial volume per mesophyll surface area [m^3 / m^2]
    Sm_std = 30,      # Sm at which assumed J_max occurs
    g_liq = 0.25e-3,   # Cell wall + liquid conductivity into stroma [m s-1]

    # position of stomata in terms of element number along x-direction
    x_abaxial_stomata = 1,
    x_adaxial_stomata = 100

  )

}

# Calculate derived parmseters for photo_2d_pm()
derive_2d_pm_parms = function(parms, fit_rubisco, ...) {

  # parms = get_2d_pm_default_parms()
  # derived quantities
  parms$t_leaf = parms[["n_z"]] * parms[["t_elem"]] # leaf thickness [m]
  parms$frac_spg = 1 - parms[["frac_pal"]] # Fraction spongy mesophyll
  parms$n_z_pal = round(parms[["n_z"]] * parms[["frac_pal"]])
  parms$n_z_spg = parms[["n_z"]] - parms[["n_z_pal"]]
  parms$Sm_z = rep(c(parms[["Sm_spg"]], parms[["Sm_pal"]]),
                 c(parms[["n_z_spg"]], parms[["n_z_pal"]]))

  parms$por_z = rep(c(parms[["por_spg"]], parms[["por_pal"]]),
                  c(parms[["n_z_spg"]], parms[["n_z_pal"]]))
  parms$por_z[parms[["n_z_spg"]]] = parms$por_int

  parms %<>%
    derive_2d_pm_rubisco(fit = fit_rubisco, ...)# |>
    # derive_2d_pm_j(...) #|>
    # derive_2d_pm_chl(...)

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
  parms$Fchl_z = rev((parms[["b0_chl"]] + parms[["b1_chl"]] * rel_depth +
    parms[["b2_chl"]] * rel_depth ^ 2) / integral01)

  parms

}

# Derive intraleaf maximum electron transport rate (j_max) gradient for photo_2d_pm()
# Assume volumetric electron transport rate is proportional Rubisco concentration
derive_2d_pm_jmax = function(parms, ...) {

  j_max_z1 = parms[["n_z"]] * parms[["X_c_z"]] / sum(parms[["X_c_z"]])
  mu_j_max_z1 = mean(j_max_z1 * parms[["Sm_z"]] * parms[["Vstrom"]])
  jscale = parms[["Jmax"]] / mu_j_max_z1
  parms$j_max_z = jscale * j_max_z1
  # Check
  checkmate::assert_true(parms[["Jmax"]] == mean(parms$j_max_z * parms[["Sm_z"]] * parms[["Vstrom"]]))

  parms

}

# Derive intraleaf potential electron transport rate (j_inf) gradient for photo_2d_pm()
# Assume potential electron transport rate is proportional Chl concentration
derive_2d_pm_jinf = function(parms, ...) {

  Jinf = parms[["I_0"]] * parms[["Beta"]] * parms[["Alpha"]] * parms[["phiPSII"]]
  j_inf_z1 = parms[["n_z"]] * parms[["Fchl_z"]] / sum(parms[["Fchl_z"]])
  mu_j_inf_z1 = mean(j_inf_z1 * parms[["Sm_z"]] * parms[["Vstrom"]])
  jscale = Jinf / mu_j_inf_z1
  parms$j_inf_z = jscale * j_inf_z1
  # Check
  checkmate::assert_true(Jinf == mean(parms$j_inf_z * parms[["Sm_z"]] * parms[["Vstrom"]]))

  parms

}

# Derive intraleaf electron transport rate (j) gradient for photo_2d_pm()
derive_2d_pm_j = function(parms, ...) {

  parms %<>%
    derive_2d_pm_jmax(...) %>%
    derive_2d_pm_jinf(...)

  parms$j_e_z = pmin(parms[["j_inf_z"]], parms[["j_max_z"]])

  parms

}

# Make matrices of porosity, Sm, X_c, j_e for photo_2d_pm()
make_2d_pm_mat = function(parms, ...) {

  parms$por_mat = matrix(parms[["por_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$Sm_mat = matrix(parms[["Sm_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$X_c_mat = matrix(parms[["X_c_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])
  parms$j_e_mat = matrix(parms[["j_e_z"]], nrow = parms[["n_z"]], ncol = parms[["n_x"]])

  parms

}

calc_2d_pm_rc = function(C_liq_mat, parms, ...) {

  w_c = (parms[["k_c"]] * parms[["X_c_mat"]] * C_liq_mat) / (parms[["K_m"]] + C_liq_mat)
  w_j = C_liq_mat * parms[["j_e_mat"]] / (4 * C_liq_mat + 8 * parms[["Gamma"]])

  r_c = pmin(w_c, w_j)

  r_c

}

# Flux calculation for photo_2d_pm()
calc_2d_pm_flux = function(C_ias_mat, parms, ...) {

  # Boundary conditions
  ## No flux from top or bottom except where stomata present
  bnd_z_bottom = C_ias_mat[1, ]
  bnd_z_top = C_ias_mat[parms[["n_z"]], ]
  bnd_z_bottom[parms[["x_abaxial_stomata"]]] = parms[["C_s"]]
  bnd_z_top[parms[["x_adaxial_stomata"]]] = parms[["C_s"]]

  ## No flux from left or right because of symmetry condition
  bnd_x_left = C_ias_mat[, 1]
  bnd_x_right = C_ias_mat[, parms[["n_x"]]]

  D_e = parms[["D_c"]] * parms[["por_mat"]] / parms[["tort"]]

  # diffusion in Z-direction; boundaries = imposed concentration
  ## flux from bottom to top
  flux_z = -rbind(
    C_ias_mat[1, ] - bnd_z_bottom,
    (C_ias_mat[2:parms[["n_z"]], ] - C_ias_mat[1:(parms[["n_z"]] - 1), ]),
    bnd_z_top - C_ias_mat[parms[["n_z"]],]
  ) / parms[["t_elem"]]

  ## flux from top to bottom
  flux_z = -(flux_z[2:(parms[["n_z"]] + 1), ] - flux_z[1:parms[["n_z"]], ]) /
    parms[["t_elem"]]

  # diffusion in X-direction; boundaries = imposed concentration
  ## flux from left to right
  flux_x = -cbind(
    C_ias_mat[, 1] - bnd_x_left,
    (C_ias_mat[, 2:parms[["n_x"]]] - C_ias_mat[, 1:(parms[["n_x"]] - 1)]),
    bnd_x_right - C_ias_mat[, parms[["n_x"]]]
  ) / parms[["t_elem"]]

  ## flux from right to left
  flux_x = -(flux_x[, 2:(parms[["n_x"]] + 1)] - flux_x[, 1:parms[["n_x"]]]) /
    parms[["t_elem"]]

  flux = D_e * (flux_z + flux_x)

  flux

}

# Reaction-diffusion function for steady.2D()
rd_2p_pm = function(t, y, parms, ...) {

  C_ias_mat = matrix(y[1:(parms[["n_x"]] * parms[["n_z"]])],
                       ncol = parms[["n_x"]], nrow = parms[["n_z"]])
  C_liq_mat = matrix(y[(parms[["n_x"]] * parms[["n_z"]] + 1):(2 * parms[["n_x"]] * parms[["n_z"]])], ncol = parms[["n_x"]], nrow = parms[["n_z"]])

  flux = calc_2d_pm_flux(C_ias_mat, parms, ...)

  r_c = calc_2d_pm_rc(C_liq_mat, parms, ...)

  r_p = r_c * parms[["Gamma"]] / C_liq_mat

  dC_ias = dC_liq = flux
  # dC_ias = (flux + parms[["g_liq"]] * (C_liq_mat - C_ias_mat) / parms[["t_elem"]]) /
  #   parms[["por_mat"]]

  # Based on Earles et al. (2017). I don't understand why it's multiple by t_leaf
  # dC_liq = (parms[["g_liq"]] * (C_ias_mat - C_liq_mat) / parms[["t_elem"]] -
  #             r_c + r_p + parms[["r_d"]]) * parms[["t_leaf"]] /
  #   sum(parms[["Sm_z"]]) / parms[["Vstrom"]]

  return(list(c(dC_ias, dC_liq))) # return output from function

}
