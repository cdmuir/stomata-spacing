# Model photosynthesis using a 2D porous medium approximation
#
# The concept is that this file will contain all functions to run 2D porous
# medium model and can be transitioned to a package script once it's ready

# Note that n_z = 0 is bottom (abaxial) side of leaf

# 2D photosynthesis model using a porous medium approximation
photo_2d_pm = function(...) {

  # Set parameters
  params = set_2d_pm_param(...)

  # Find solution
  # Structure output

}

# Set parameters for photo_2d_pm()
set_2d_pm_param = function(...) {

}

# Get set of parameters for photo_2d_pm()
get_2d_pm_default_param = function(...) {

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
    R_d = 0.066,       # Dark respiratory rate [mol / m^3 / s]
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
    b1_chl = 100 * 0.4149, # converts parameters to 0-1 rather than 0-100 scale
    b2_chl = 100 ^ 2 * -0.0029, # converts parameters to 0-1 rather than 0-100 scale

    frac_pal = 0.6,   # Fraction palisade mesophyll [1]
    Sm_spg = 6.5,     # Sm spongy mesophyll [m^2 / m^2]
    Sm_pal = 40,      # Sm palisade mesophyll [m^2 / m^2]
    Vstrom = 1.74e-6, # Stroma volume per mesophyll surface area [m^3 / m^2]
    Vmito = 0.27e-7,  # Mitochondrial volume per mesophyll surface area [m^3 / m^2]
    Sm_std = 30,      # Sm at which assumed J_max occurs
    g_liq = 0.25e-3   # Cell wall + liquid conductivity into stroma [m s-1]

  )

}

# Calculate derived parameters for photo_2d_pm()
derive_2d_pm_param = function(param, fit_rubisco, ...) {

  # param = get_2d_pm_default_param()
  # derived quantities
  param$t_leaf = param[["n_z"]] * param[["t_elem"]] # leaf thickness [m]
  param$frac_spg = 1 - param[["frac_pal"]] # Fraction spongy mesophyll
  param$n_z_pal = round(param[["n_z"]] * param[["frac_pal"]])
  param$n_z_spg = param[["n_z"]] - param[["n_z_pal"]]
  param$Sm_z = rep(c(param[["Sm_spg"]], param[["Sm_pal"]]),
                 c(param[["n_z_spg"]], param[["n_z_pal"]]))

  param$por_z = rep(c(param[["por_spg"]], param[["por_pal"]]),
                  c(param[["n_z_spg"]], param[["n_z_pal"]]))
  param$por_z[param[["n_z_spg"]]] = param$por_int

  param %<>%
    derive_2d_pm_rubisco(fit = fit_rubisco, ...)# |>
    # derive_2d_pm_j(...) #|>
    # derive_2d_pm_chl(...)

  param

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
derive_2d_pm_rubisco = function(param, fit, integral01 = NULL, ...) {

  rel_depth = seq(0 + 1 / param[["n_z"]] / 2, 1 - 1 / param[["n_z"]] / 2,
                  length.out = param[["n_z"]])
  rel_rubisco = predict(fit, newdata = tibble(rel_depth = rel_depth))
  if (is.null(integral01)) {
    integral01 = integrate(\(x, fit) predict(fit, newdata = tibble(rel_depth = x)),
                           lower = 0, upper = 1, fit = fit)[["value"]]
  }
  param$X_c_z = rel_rubisco * param[["X_c"]] / integral01

  param

}

# Derive intraleaf chlorophyll gradient for photo_2d_pm()
derive_2d_pm_chl = function(param, ...) {

  # Using quadratic equation proposed by Johnson et al. (2005), see also
  # Borsuk and Brodersen (2019)
  rel_depth = seq(0 + 1 / param[["n_z"]] / 2, 1 - 1 / param[["n_z"]] / 2,
                  length.out = param[["n_z"]])

  # analytical solution for integral of quadratic from 0 to 1
  integral01 = param[["b0_chl"]] + param[["b1_chl"]] / 2 +
    param[["b2_chl"]] / 3

  # **NEED TO REVERSE BECAUSE BORSUK AND BRODERSEN HAVE 0 AT ADAXIAL SURFACE**
  param$Fchl_z = rev((param[["b0_chl"]] + param[["b1_chl"]] * rel_depth +
    param[["b2_chl"]] * rel_depth ^ 2) / integral01)

  param

}

# Derive intraleaf maximum electron transport rate (j_max) gradient for photo_2d_pm()
# Assume volumetric electron transport rate is proportional Rubisco concentration
derive_2d_pm_jmax = function(param, ...) {

  j_max_z1 = param[["n_z"]] * param[["X_c_z"]] / sum(param[["X_c_z"]])
  mu_j_max_z1 = mean(j_max_z1 * param[["Sm_z"]] * param[["Vstrom"]])
  jscale = param[["Jmax"]] / mu_j_max_z1
  param$j_max_z = jscale * j_max_z1
  # Check
  checkmate::assert_true(param[["Jmax"]] == mean(param$j_max_z * param[["Sm_z"]] * param[["Vstrom"]]))

  param

}

# Derive intraleaf potential electron transport rate (j_inf) gradient for photo_2d_pm()
# Assume potential electron transport rate is proportional Chl concentration
derive_2d_pm_jinf = function(param, ...) {

  Jinf = param[["I_0"]] * param[["Beta"]] * param[["Alpha"]] * param[["phiPSII"]]
  j_inf_z1 = param[["n_z"]] * param[["Fchl_z"]] / sum(param[["Fchl_z"]])
  mu_j_inf_z1 = mean(j_inf_z1 * param[["Sm_z"]] * param[["Vstrom"]])
  jscale = Jinf / mu_j_inf_z1
  param$j_inf_z = jscale * j_inf_z1
  # Check
  checkmate::assert_true(Jinf == mean(param$j_inf_z * param[["Sm_z"]] * param[["Vstrom"]]))

  param

}

# Derive intraleaf electron transport rate (j) gradient for photo_2d_pm()
derive_2d_pm_j = function(param, ...) {

  param %<>%
    derive_2d_pm_jmax(...) %>%
    derive_2d_pm_jinf(...)

  param$j_e_z = pmin(param[["j_inf_z"]], param[["j_max_z"]])

  param

}

# Calculation A for photo_2d_pm()
calc_2d_pm_An = function(C_liq, param) {
  # NEED TO CONVERT X_c_z to n_z x n_x matrix
  pmin((param[["k_c"]] * param[["X_c_z"]] * C_liq) / (param[["K_m"]] + C_liq),
       parms[["n_y"]] * (C_liq * param[["j_e"]] / (4 * C_liq + 8 * param[["Gamma"]])))
}

# Make matrices of porosity, Sm, X_c, j_e for photo_2d_pm()
make_2d_pm_mat = function(param, ...) {

  param$por_mat = matrix(param[["por_z"]], nrow = param[["n_z"]], ncol = param[["n_x"]])
  param$Sm_mat = matrix(param[["Sm_z"]], nrow = param[["n_z"]], ncol = param[["n_x"]])
  param$X_c_mat = matrix(param[["X_c_z"]], nrow = param[["n_z"]], ncol = param[["n_x"]])
  param$j_e_mat = matrix(param[["j_e_z"]], nrow = param[["n_z"]], ncol = param[["n_x"]])

  param

}

#
# WORKING HERE
param = get_2d_pm_default_param() |>
  derive_2d_pm_param(fit_rubisco = fit_rubisco) |>
  derive_2d_pm_chl() |>
  derive_2d_pm_j() |>
  make_2d_pm_mat()

# ADD THIS TO set_2dm_pm_param()
param[["x_abaxial_stomata"]] = 1
param[["x_adaxial_stomata"]] = 100


# state
C_ias = rep(0.015 / 2, param[["n_z"]] * param[["n_x"]])
C_liq = rep(0.015 / 4, param[["n_z"]] * param[["n_x"]])
C_ias_mat = matrix(C_ias, nrow = param[["n_z"]], ncol = param[["n_x"]])
C_liq_mat = matrix(C_liq, nrow = param[["n_z"]], ncol = param[["n_x"]])

# Boundary conditions
## No flux from top or bottom except where stomata present
bnd_z_bottom = C_ias_mat[1, ]
bnd_z_top = C_ias_mat[param[["n_z"]], ]
bnd_z_bottom[param[["x_abaxial_stomata"]]] = param[["C_s"]]
bnd_z_top[param[["x_adaxial_stomata"]]] = param[["C_s"]]

## No flux from left or right because of symmetry condition
bnd_x_left = C_ias_mat[, 1]
bnd_x_right = C_ias_mat[, param[["n_x"]]]

D_e = param[["D_c"]] * param[["por_mat"]] / param[["tort"]]

# diffusion in Z-direction; boundaries = imposed concentration
## flux from bottom to top
flux_z = -rbind(
  C_ias_mat[1, ] - bnd_z_bottom,
  (C_ias_mat[2:param[["n_z"]], ] - C_ias_mat[1:(param[["n_z"]] - 1), ]),
  bnd_z_top - C_ias_mat[param[["n_z"]],]
) / param[["t_elem"]]

## flux from top to bottom
flux_z = -(flux_z[2:(param[["n_z"]] + 1), ] - flux_z[1:param[["n_z"]], ]) /
  param[["t_elem"]]

# diffusion in X-direction; boundaries = imposed concentration
## flux from left to right
flux_x = -cbind(
  C_ias_mat[, 1] - bnd_x_left,
  (C_ias_mat[, 2:param[["n_x"]]] - C_ias_mat[, 1:(param[["n_x"]] - 1)]),
  bnd_x_right - C_ias_mat[, param[["n_x"]]]
) / param[["t_elem"]]

## flux from right to left
flux_x = -(flux_x[, 2:(param[["n_x"]] + 1)] - flux_x[, 1:param[["n_x"]]]) /
  param[["t_elem"]]

flux = D_e * (flux_z + flux_x)

#CHECK
flux[1:2, 1:2] # bottom-left
flux[199:200, 99:100] # top-right

W_c = (param[["k_c"]] * param[["X_c_mat"]] * C_liq_mat) / (param[["K_m"]] + C_liq_mat)
W_j = C_liq_mat * param[["j_e_mat"]] / (4 * C_liq_mat + 8 * param[["Gamma"]])

A_n = pmin(W_c, W_j)

R_p = A_n * param[["Gamma"]] / C_liq_mat

dC_ias = (flux + param[["g_liq"]] * (C_liq_mat - C_ias_mat) / param[["t_elem"]]) /
  param[["por_mat"]]

dC_liq = (param[["g_liq"]] * (C_liq_mat - C_ias_mat) / param[["t_elem"]] -
  A_n + R_p + param[["R_d"]])

#CHECK
C_ias_mat[1:2, 1:2] # bottom-left
dC_ias[1:2, 1:2] # bottom-left
dC_ias[199:200, 99:100] # top-right
dC_liq[1:2, 1:2] # bottom-left
dC_liq[199:200, 99:100] # top-right


add_z_boundaries = function(C_ias_mat, ...) {

  C_ias_mat1 = matrix(NA, nrow = nrow(C_ias_mat) + 2L, ncol = ncol(C_ias_mat))
  C_ias_mat1[1,] = C_ias_mat[1,]
  C_ias_mat1[seq_len(nrow(C_ias_mat)) + 1L,] = C_ias_mat
  C_ias_mat1[nrow(C_ias_mat) + 2L,] = C_ias_mat[nrow(C_ias_mat),]

  param$x_abaxial_stomata = 1
  param$x_adaxial_stomata = 4
  C_ias_mat1[1, param[["x_abaxial_stomata"]]] = param[["C_s"]]
  C_ias_mat1[nrow(C_ias_mat) + 2L, param[["x_adaxial_stomata"]]] = param[["C_s"]]

  C_ias_mat1

}

# Calculate C flux for photo_2d_pm()
calc_2d_pm_flux = function(C_ias, param, ...) {

  # Ignoring stomatal distribution argument for now.
  # Only using amphi model.
  # sr = match.arg(stomatal_distribution, c("amphi", "hypo"))
  sr = "amphi"

  # assume C_ias is vector with length n_x * n_z
  C_ias_mat = matrix(C_ias, nrow = param[["n_z"]], ncol = param[["n_x"]])
  dC = numeric(length(C_ias))


  C_ias_mat1 = add_z_boundaries(C_ias_mat)
  C_ias_mat1[2:nrow(C_ias_mat1), ] - C_ias_mat1[seq_len(nrow(C_ias_mat1) - 1)]

  #diffusion in Y-direction
  Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
  dY    <- dY - (Flux[,2:(n+1)]-Flux[,1:n])/dy

  tmp = matrix(1:12, ncol = 3, nrow = 4)
  tmp[2:4,] - tmp[1:3,]

  if (sr == "amphi") {
    # need to change boundary condition for offset
    boundary = matrix(c(param[["C_s"]], C_ias[1, 2:param[["n_x"]]]),
                      nrow = 1, ncol = param[["n_x"]])

    # CHECK THE POROSITY HERE - use intermediate porosity for boundaries
    D_e = matrix(param[["Dz"]] * (param[["por_int"]] / param[["tort"]]),
                 nrow = param[["n_z"]] + 1, ncol = param[["n_x"]])
    flux1 = -D_e * rbind(
      C_ias[1, ] - boundary,
      C_ias[2:param[["n_z"]], ] - C_ias[1:(param[["n_z"]] - 1), ],
      rev(boundary) - C_ias[param[["n_z"]], ]
    ) / param[["t_elem"]]

    dC = dC - (flux1[2:(param[["n_z"]] + 1), ] - flux1[1:param[["n_z"]], ]) / param[["t_elem"]]

    D_e = matrix(param[["Dz"]] * (param[["por_z"]] / param[["tort"]]),
                 nrow = param[["n_z"]], ncol = param[["n_x"]] + 1)
    flux2 = -D_e * cbind(
      rep(0, param[["n_z"]]),
      C_ias[, 2:param[["n_x"]]] - C_ias[, 1:(param[["n_x"]] - 1)],
      rep(0, param[["n_z"]])
    ) / param[["t_elem"]]

    dC = dC - (flux2[, 2:(param[["n_x"]] + 1)] - flux2[, 1:param[["n_x"]]]) / param[["t_elem"]]

  }

  dC

}
# This should be "generic" since it would presumably be the same calculation
# with different state/parameter inputs
calc_2d_flux

n = param[["n_x"]] * param[["n_z"]]
C_ias = state[seq_len(n)] # Define vapor [CO2]
C_liq = state[n + seq_len(n)] # Define liquid [CO2]

flux_cias = calc_2d_pm_flux(C_ias, param)

dC_ias = (flux_cias + param[["g_liq"]] * (C_liq - C_ias) / param[["t_elem"]]) / param[["por_z"]] # original eqn divided by porosity, but I don't understand why

# Calculate carboxylation rate
An = calc_2d_pm_An(C_liq, param)

Rp = (An * param[["Gamma"]]) / C_liq # Calculate oxygenation rate

dC_liq = (param[["g_liq"]] * (C_ias - C_liq) / param[["t_elem"]] - An + Rp + param[["R_d"]]) * param[["t_leaf"]] / param[["Sm_z"]] / param[["Vstrom"]]


state = rep(c(0.5, 0.25), each = param[["n_x"]] * param[["n_z"]]) * 1.72e-2 # CO2 concentration [mol / m^3]
