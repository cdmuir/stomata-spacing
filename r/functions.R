# Custom functions

scientize <- function(x, threshold = -1L, digits = 2L) {

  purrr::map_chr(x, function(.x, threshold, digits) {

    if (is.na(.x)) {return(NA)}
    oom <- log10(abs(.x))
    if (oom < threshold) {
      x1 <- .x |>
        magrittr::multiply_by(10 ^ -floor(oom)) |>
        round(digits)

      x2 <- sprintf(glue::glue("%.{digits}f"), x1) |>
        stringr::str_c(" \\times 10^{", floor(oom), "}")

      return(x2)
    } else {
      x1 <- .x |>
        signif(digits + 1L) |>
        as.character()
      return(x1)
    }

  }, threshold = threshold, digits = digits)

}
## Milestone 1: Simulate and plot stomatal grids ----

### Functions for converting from D to U ----
# D = stomatal density in 1 / pixel^2
# U = interstomatal distance
D2U <- function(D) {

  checkmate::assert_numeric(D, lower = 0, finite = TRUE, any.missing = FALSE,
                            min.len = 1L)
  # D must have units of 1 / pixel ^ 2
  # U = 2 * Rhex where Rhex is major radius of a regular hexagon
  # Ahex = 3 * sqrt(3) / 2 * Rhex ^ 2

  # Ahex is the area of the hexagon that each stomata would occupy in a
  # perfectly dispersed array of stomata with density (D).
  Ahex <- 1 / D

  # Rhex is the outcircle radius Ahex calculated using equations which describe
  # the radius of a regular hexagon
  Rhex <- sqrt(2 * Ahex / (3 * sqrt(3)))

  # a is the apothem, or incircle radius of a regular hexagon.
  a <- (sqrt(3) / 2) * Rhex

  # U is the distance between stomata assuming they are separated by two a.
  U <- 2 * a # in pixels

  # return U
  return(U)

}

# a function that converts distance to density

U2D <- function(U) {

  # checks
  checkmate::assert_numeric(U, lower = 0, finite = TRUE, any.missing = FALSE,
                            min.len = 1L)
  # U must have units of pixel
  # a is the apothem or incircle radius of a regular hexagon
  a <- U / 2

  # Rhex is the outcircle radius of a regular hexagon
  Rhex <- (2 * a) / sqrt(3)

  # Ahex is the area of a regular hexagon
  Ahex <- 3 * sqrt(3) / 2 * Rhex ^ 2

  # D is the density of stomata assuming distance U and perfect dispersion
  D <- 1 / Ahex # 1 / pixel ^ 2

  # return D
  return(D)

}

# Rotate a set of XY coordinates by an angle (in radians)
# From spdep::Rotation
rotate_xy <- function (xy, angle) {
  xy <- as.matrix(xy)
  cos.angle <- cos(angle)
  sin.angle <- sin(angle)
  xy.rot <- xy %*% t(matrix(c(cos.angle, sin.angle, -sin.angle,
                              cos.angle), 2, 2))
  return(xy.rot)
}

make_equal_stomata_grid <- function(
  U, # distance between stomata assuming perfect dispersion
  pixels_x, # number of pixels on the horizontal picture axis
  pixels_y, # number of pixels on the vertical picture axis
  jitter, # Logical, if TRUE, jitter is added to the stomatal grid
  rotate # Logical, if TRUE, rotation is added to the stomatal grid
) {

  # Checks
  checkmate::assert_numeric(U, lower = 0L, len = 1L)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_flag(jitter)
  checkmate::assert_flag(rotate)

  # Calculation x and y positions of stomata separated by U distance
  # U equals distance between stomata in pixels
  # Equilateral triangular grid ---- optimally dispersed stomata

  # where to start and end - scaled for any grid size and stomatal density
  start_end = ceiling((((pixels_x * pixels_y) ^ 0.5) / 100) * U * -1)

  # sequence of x values
  x_seq = c(
    seq(from = start_end, to = (pixels_x - start_end), by = U),
    seq(
      from = start_end + (U / 2),
      to = ((pixels_x - start_end + (U / 2))),
      by = U
    )
  )

  # number of reps to rep the sequence
  x_rep = ceiling(((pixels_x * pixels_y) ^ 0.5) / 150) + ceiling((sqrt(2) * (pixels_y) / (2 * (sqrt(3)) / 2 * U)))
  # vector of x values
  x = rep(x_seq, x_rep)

  # sequence of y values
  y_seq = seq(
    from = start_end,
    by = ((sqrt(3)) / 2) * U,
    length.out = (2 * x_rep)
  )
  # number of reps to rep the sequence
  y_rep = length(x) / (2 * x_rep)
  # sort the sequence to match the x values
  y = sort(rep(y_seq, y_rep))

  # building a dataframe of x, y coordinates
  # if jitter = TRUE, at a uniform random number between 0 and U / 2 (x)
  # or 0 and sqrt(3)/4 * U (y) across each vector
  # if else statement for rotation
  if (!rotate) {

    df_stomata = tibble(
      x = x + jitter * runif(1, 0, U / 2),
      y = y + jitter * runif(1, 0, sqrt(3) / 2 * U)
    ) |>
      dplyr::filter(x <= pixels_x, y <= pixels_y, x >= 0, y >= 0) |>
      dplyr::mutate(n_stomata = n())

    # Return tibble with x and y positions of stomata
    return(df_stomata)

  } else {

    # center at origin
    df_stomata = tibble(x = x, y = y) |>
      dplyr::mutate(
        x_centered = x - pixels_x / 2,
        y_centered = y - pixels_y / 2
      )

    df_rotated = df_stomata |>
      dplyr::select(x_centered, y_centered) |>
      rotate_xy(runif(1, 0, 45) * pi / 180) |>
      as_tibble(.name_repair = ~ c("x", "y")) |>
      dplyr::mutate(
        x = x + (pixels_x / 2) + jitter * runif(1, 0, U / 2),
        y = y + (pixels_y / 2) + jitter * runif(1, 0, sqrt(3) / 2 * U)
      ) |>
      dplyr::filter(x <= pixels_x, y <= pixels_y, x >= 0, y >= 0) |>
      dplyr::mutate(n_stomata = n())

    # return
    return(df_rotated)

  }
}

# a function that makes a random stomata grid

make_random_stomata_grid = function(
  n, # number of stomata
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
) {

  # Checks
  checkmate::assert_integerish(n, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # randomly sample with uniform distribution using runif() and put into a
  # tibble
  df_stomata = tibble(
    x = runif(n = n, min = 0, max = pixels_x),
    y = runif(n = n, min = 0, max = pixels_y),
    n_stomata = n
  )

  # Return tibble with x and y positions of stomata
  return(df_stomata)

}

# function for plotting stomatal grid

plot_stomatal_grid = function(
  df_stomata, # a tibble of x, y coordinates
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
) {

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # plot
  stom_plot = ggplot(df_stomata, mapping = aes(x, y),
                     lims = c(pixels_x, pixels_y)) +
    geom_point() +
    geom_rect(aes(xmin = 0, xmax = pixels_x, ymin = 0, ymax = pixels_y),
              fill = NA, color = "black") +
    ggtitle("Stomatal Grid") +
    coord_equal() +
    theme_void()

  # return the plot
  return(stom_plot)

}

# make ideal amphi grid

make_amphi_grid = function(
  U, #interstomatal distance
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  jitter, # logical; jitter point locations
  rotate # logical; rotate points
) {

  # checks
  checkmate::assert_numeric(U, lower = 0, len = 1L)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_flag(jitter)
  checkmate::assert_flag(rotate)

  # make abaxial surface
  abaxial = make_equal_stomata_grid(U, pixels_x, pixels_y, jitter, rotate)

  # make adaxial surface; add to optimally offset
  adaxial = abaxial
  adaxial$y = adaxial$y - (U * sqrt(3)/3)
  # filter out points outside leaf surface
  adaxial = dplyr::filter(adaxial, x <= pixels_x, y <= pixels_y, x >= 0, y >= 0)

  # define surfaces
  abaxial$surface = rep("abaxial", nrow(abaxial))
  adaxial$surface = rep("adaxial", nrow(adaxial))

  # combine datasets
  df_stomata = rbind(abaxial, adaxial)

  # return
  return(df_stomata)

}

# plot amphi grid

plot_amphi_grid = function(
  df_stomata, # a tibble of x, y coordinates for abaxial and adaxial stomata
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
) {

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 3L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y", "surface"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # plot
  stom_plot = ggplot(df_stomata, mapping = aes(x, y, color = surface), lims = c(pixels_x, pixels_y)) +
    coord_equal() +
    geom_ellipse(aes(x0 = x, y0 = y,
                     a = 15, b = 7.5, angle = pi/2)) +
    geom_ellipse(aes(x0 = x, y0 = y,
                     a = 1.5, b = 0.75, angle = pi/2)) +
    theme_void() +
    scale_color_manual(values = c("black", "grey")) +
    geom_rect(aes(xmin = 0, xmax = pixels_x, ymin = 0, ymax = pixels_y),
              fill = NA, color = "black") +
    ggtitle("    Ideal Stomatal Patterning")

  # return the plot
  return(stom_plot)

}

## Milestone 2: Simulate stomatal grids based on data ----

simulate_grids_from_data <- function(
  data, # a tibble of x, y positions of stomata
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  n_grid, # number of grids wanted
  progress = FALSE # progress bar logical
) {

  # Checks
  checkmate::assert_tibble(data, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(data))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_integerish(n_grid, lower =  1L, len = 1L)
  checkmate::assert_flag(progress)

  n_stomata <- nrow(data)

  # Draw n_stomata for synthetic data - removed this because I realized even though we don't know true D, we should condition on observed n_stomata
  # n_stomata_sim <- rpois(n_grid, lambdas)

  # Simulate uniform random grids
  random_grids <- rep(n_stomata, n_grid) |>
    purrr::map(make_random_stomata_grid, pixels_x = pixels_x,
               pixels_y = pixels_y) |>
    purrr::imap_dfr(~ {dplyr::mutate(.x, sim = .y)}) |>
    dplyr::mutate(
      sim = stringr::str_pad(sim, pad = "0",
                             width = floor(log10(n_grid) + 1)),
      grid = "random"
    )

  # Simulate equal distance grids with jitter and rotation
  # Use Gamma sampling distribution to simulate from plausible range of density
  # values (lambda)
  # I think this is appropriate sampling distribution from the mean of Poisson
  # with flat prior (alpha -> 0, beta -> 0)
  # see https://en.wikipedia.org/wiki/Poisson_distribution#Bayesian_inference
  # lambdas = rgamma(n_grid, shape = n_stomata, scale = 1)
  # Calculate U from density
  # U = D2U(lambdas / (pixels_x * pixels_y))

  # OLD WAY - much faster
  # equal_grids = U |>
  #   purrr::map(make_equal_stomata_grid, pixels_x = pixels_x,
  #              pixels_y = pixels_y, jitter = TRUE, rotate = TRUE) |>
  #   purrr::imap_dfr(~ {dplyr::mutate(.x, sim = .y)}) |>
  #   dplyr::mutate(
  #     sim = stringr::str_pad(sim, pad = "0",
  #                            width = floor(log10(n_grid) + 1)),
  #     grid = "equal"
  #   )

  # NEW WAY - much slower
  equal_grids = vector(mode = "list", length = n_grid)
  i <- 0

  while (i < n_grid) {
    lambda = rgamma(1, shape = n_stomata, scale = 1)
    U = D2U(lambda / (pixels_x * pixels_y))
    g = make_equal_stomata_grid(U, pixels_x = pixels_x, pixels_y = pixels_y,
                                 jitter = TRUE, rotate = TRUE)
    if (g$n_stomata[1] == n_stomata) {
      i = i + 1
      equal_grids[[i]] = g
    }
  }

  equal_grids <- equal_grids |>
    purrr::imap_dfr(~ {dplyr::mutate(.x, sim = .y)}) |>
    dplyr::mutate(
      sim = stringr::str_pad(sim, pad = "0",
                             width = floor(log10(n_grid) + 1)),
      grid = "equal"
    )

  # Return tibble
  simulations = bind_rows(random_grids, equal_grids)
  if (progress) pb$tick()
  return(simulations)

}

# get nearest neighbor index (nni)

get_nni = function(
  df_stomata, # tibble of x, y coordinates of stomata
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  progress = FALSE # progress bar logical
){

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_flag(progress)
  checkmate::assert_numeric(df_stomata$x, lower = 0, upper = pixels_x)
  checkmate::assert_numeric(df_stomata$y, lower = 0, upper = pixels_y)

  df_stomata_sp <- df_stomata |>
    dplyr::select(x, y) |>
    sp::SpatialPoints()

  df_stomata_sp@bbox <- matrix(c(0, 0, pixels_x, pixels_y), ncol = 2)

  # perform average nearest neighbor test
  nni = nni(df_stomata_sp, win = "extent")
  nni$n_stomata <- nrow(df_stomata)

  # return result
  if (progress) pb$tick()
  return(as_tibble(nni))

}

# a function that gets the variance in area occupied by each stomate

#UNFINISHED
get_variance = function(
  df_stomata,
  pixels_x,
  pixels_y,
  progress = FALSE
){

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_flag(progress)
  checkmate::assert_numeric(df_stomata$x, lower = 0, upper = pixels_x)
  checkmate::assert_numeric(df_stomata$y, lower = 0, upper = pixels_y)

  #
}


## Functions for stomata on both leaf surfaces ----

# starting from raw data

# a function that uses the sp and raster packages to make a kde
kerneldensity = function(
  data, # x, y stomatal coordinates
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
) {

  # checks
  checkmate::assert_tibble(data, min.cols = 2L, any.missing = F,
                               col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(data))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)


  # change to spatial points data frame using the sp package
  xy = data[, c("x", "y")]
  stomata_sp = SpatialPointsDataFrame(xy, data)

  # kernel density function
  # define extent of analysis
  e = extent(0, pixels_x, 0, pixels_y)
  # run sp.kde which builds a kernel density estimate
  kernelden = sp.kde(x = stomata_sp, bw = 120,
                      nr = pixels_y, nc = pixels_x,
                      newdata = e, standardize = TRUE,
                      mask = FALSE)

  # return result
  return(kernelden)

}

# a function that builds a base raster for dist_sq function

create_raster = function(
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
){

  # checks
  checkmate::assert_integerish(pixels_x)
  checkmate::assert_integerish(pixels_y)

  # build a raster
  leaf = raster(ncol = pixels_x, nrow = pixels_y,
                xmx = pixels_x, xmn = 0,
                ymx = pixels_y, ymn = 0)

  # return raster
  return(leaf)

}

# a function that calculates the distance squared from each stomate

dist_sq = function(
  data, # a tibble of x, y coordinates of stomata
  leaf # null raster of dimensions pixels_x, pixels_y
){

  # checks
  checkmate::assert_tibble(data, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(data))

  # prepare data
  # change to spatial points data frame using the sp package
  xy = data[, c("x", "y")]
  stomata_sp = SpatialPointsDataFrame(xy, data)

  # calculate distance from nearest stomata using function from raster package
  dist_raster = raster::distanceFromPoints(leaf, stomata_sp)

  # square that distance
  dist_sq_raster = dist_raster^2

  # return raster
  return(dist_sq_raster)

}

# function that returns correlation value for paired dataframe

amphi_corr = function(
  data, # grouped dataframe
  leaf, # baseline raster of dimensions pixels_x, pixels_y
  progress = FALSE # progress bar logical
){

  checkmate::assert_data_frame(data, min.cols = 3L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("surface"), colnames(data)) # x,y will be checked in dist_sq
  checkmate::assert_class(leaf, "RasterLayer")
  checkmate::assert_flag(progress)

  abaxial = data |>
    dplyr::filter(surface == "abaxial")
  adaxial = data |>
    dplyr::filter(surface == "adaxial")

  abaxial_raster = dist_sq(abaxial, leaf)
  adaxial_raster = dist_sq(adaxial, leaf)

  # stack rasters using raster package stack function
  stack = raster::stack(abaxial_raster, adaxial_raster)

  # pearson correlation using raster's layerStats function
  coefficients = layerStats(stack, stat = 'pearson')

  corr_coef = coefficients$`pearson correlation coefficient`["layer.1", "layer.2"]

  # return coefficients
  if (progress) pb$tick()
  return(corr_coef)

}

corr_pair = function(
  df_pair, # list of stomatal grids
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
){

  # checks
  checkmate::assert_list(df_pair)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # make kernel density plots
  kernels = kerneldensities(df_pair, pixels_x, pixels_y)

  # stack rasters using raster package stack function
  stack = raster::stack(kernels[[1]], kernels[[2]])

  # pearson correlation using raster's layerStats function
  coefficients = layerStats(stack, stat = 'pearson')

  # return coefficients
  return(coefficients)

}

# a function that perfroms corr_pairs for every item in a list

corr_pairs = function(
  list, # list of lists: paired stomatal grids
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
){

  # checks
  checkmate::assert_list(list)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # perform tests
  corr_pairs = lapply(list, corr_pair,
                        pixels_x, pixels_y)

  # return result
  return(corr_pairs)

}

# a function that runs a correlation test on paired dist_sq rasters

dist_sq_corr = function(
  df_pair,
  leaf
){

  # checks
  checkmate::assert_list(df_pair)

  # build dist_sq rasters
  dist_sq = dist_sqs(df_pair, leaf)

  # stack rasters using raster package stack function
  stack = raster::stack(dist_sq[[1]], dist_sq[[2]])

  # pearson correlation using raster's layerStats function
  coefficients = layerStats(stack, stat = 'pearson')

  # return coefficients
  return(coefficients)

}

# a function that applies dist_sq_corr over each element in a list

dist_sq_corrs = function(
  list,
  leaf
){

  # checks
  checkmate::assert_list(list)

  # lapply
  dist_sq_corrs = lapply(list, dist_sq_corr, leaf)

  # return
  return(dist_sq_corrs)

}

# function that converts corr_pairs result to dataframe

corr_to_frame = function(
  list # output of corr_pairs() or dist_sq_corrs function; list of correlation coefficients and raster averages
){

  #checks
  checkmate::assert_list(list)

  # create a dataframe
  result = data.frame(matrix(unlist(list), nrow = (sum(lengths(list) / 2)),
                             byrow = TRUE), stringsAsFactors = FALSE)
  result = result[, -c(1,2,4)] |>
    rename(cor_coef = X3,
           abax_mean = X5,
           adax_mean = X6)

  # make dataframe into tibble
  result = tibble(result)

  # return the result
  return(result)
}

# a function that compares correlation values to simulated data

real_corr_test = function(
  df_pair, # list of paired stomatal grids; x, y coordinates
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  n_grid, # number of grids to simulate,
  dist_or_kern # logical; TRUE = dist_sq; FALSE = kernel density
){

  # checks
  checkmate::assert_list(df_pair)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_integerish(n_grid, lower =  1L, len = 1L)
  checkmate::assert_flag(dist_or_kern)


  # define
  abaxial = df_pair[[1]]
  adaxial = df_pair[[2]]

  # perform corr test on real data
  # if else statement to determine to dist_sq or kernel
  if (dist_or_kern == TRUE) {
    leaf = create_raster(pixels_x, pixels_y)
    corr = dist_sq_corr(df_pair, leaf)
    corr_coef = corr$`pearson correlation coefficient`[2,1]
  } else {
    corr = corr_pair(df_pair, pixels_x, pixels_y)
    corr_coef = corr$`pearson correlation coefficient`[2,1]
  }

  # run simulation
  simulation = simulate_grids_from_pairs(df_pair, pixels_x, pixels_y, n_grid = n_grid)
  if (dist_or_kern == TRUE) {
    leaf2 = create_raster(pixels_x, pixels_y)
    simulation = dist_sq_corrs(simulation, leaf2)
    simulation = corr_to_frame(simulation)
    } else {
    simulation = corr_pairs(simulation, pixels_x, pixels_y)
    simulation = corr_to_frame(simulation)
  }

  # compare real to simulations
  # p-value is the proportion of random simulated data greater than result for real
  # calculate p value
  p = sum(simulation$cor_coef < corr_coef) / nrow(simulation)

  # number of stomata
  n_stom_ab = as.character(nrow(abaxial))
  n_stom_ad = as.character(nrow(adaxial))

  # metadata
  treatment = as.character(abaxial$treatment[1])
  leaf_number = as.character(abaxial$leaf_number[1])
  image_number = as.character(gsub("^.{0,3}", "", abaxial$image_number[1]))

  result = list(p, corr_coef, n_stom_ab, n_stom_ad, treatment, leaf_number, image_number)

  # return the result
  return(result)

}

# a function that performs real_corr_test for every item in a list

real_corr_tests = function(
  list, # list of lists of associated stomatal grids
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  n_grid, # number of simulations desired
  dist_or_kern, # logical; TRUE = dist_sq; FALSE = kernel density
  num_cores # number of cores to run code on
) {

  # checks
  checkmate::assert_list(list)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_integerish(n_grid, lower =  1L, len = 1L)
  checkmate::assert_flag(dist_or_kern)
  checkmate::assert_integerish(num_cores, lower =  1L, len = 1L)

  #perform tests
  tests = mclapply(list, real_corr_test,
                 pixels_x, pixels_y, n_grid, dist_or_kern, mc.cores = num_cores)

  # return result
  return(tests)

}


# list to frame function for final result
corr_to_frame_final = function(
  list # output from real_corr_tests(); list of lists of test results
){

  # checks
  checkmate::assert_list(list)

  # create a dataframe
  result <- data.frame(matrix(unlist(list), nrow = sum(lengths(list))/7,
                              byrow = TRUE), stringsAsFactors = FALSE) |>
    rename(
      p = X1,
      corr_coef = X2,
      n_stom_ab = X3,
      n_stom_ad = X4,
      treatment = X5,
      leaf_number = X6,
      image_number = X7
    )
  # add column for light treatment
  result$light = ifelse(grepl("HL", result$treatment), "High Light",
                        ifelse(grepl("LL", result$treatment), "Low Light",
                               "Medium Light"))
  result$stom_ratio = as.numeric(result$n_stom_ad) / as.numeric(result$n_stom_ab)

  # make dataframe into tibble
  result = tibble(result)

  # return the result
  return(result)

}

# function that plots a histogram for amphistomatic data
amphi_stom_hist = function(
  result, # output of corr_to_frame_final; real data and p-value
  simulation # output of corr_to_frame; simulated data
) {

  # checks
  checkmate::assert_tibble(result, min.cols = 9L, any.missing = FALSE,
                               col.names = "named")
  checkmate::assert_subset(c("corr_coef", "n_stom_ab","n_stom_ad", "p", "treatment",
                             "leaf_number", "light", "stom_ratio", "image_number"), colnames(result))
  checkmate::assert_tibble(simulation, min.cols = 3L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("cor_coef", "abax_mean","adax_mean"), colnames(simulation))


  # build histogram
  histogram <- simulation |>
    ggplot(aes(x=cor_coef)) +
    geom_histogram(color = "#e9ecef", alpha=0.6, position = 'identity',
                   binwidth = 0.05) +
    geom_vline(xintercept = as.numeric(result$corr_coef)) +
    annotate("text", x = 0, y = 20, label = paste(result$p)) +
    annotate("text", x = 0, y = 25, label = paste(result$corr_coef)) +
    theme_pubr() +
    labs(fill="")
  histogram

  # return histogram
  return(histogram)

}

## Tessellation Functions ----

tessellate_stomata = function(
  df_stomata, # tibble of stomatal x, y positions
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
){

  #checks
  checkmate::assert_tibble(df_stomata, min.cols = 2L, any.missing = FALSE,
                          col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # create a voronoi object
  voro_stom = voronoi.mosaic(df_stomata$x, df_stomata$y, duplicate = "error")

  # calculate cell areas
  voro_cell = tripack::cells(voro_stom)

  # convert to a dataframe
  voro_cell_df = as.data.frame(do.call(rbind, voro_cell))

  # identify rejects
  reject = voronoi.findrejectsites(voro_stom, 0, pixels_x, 0, pixels_y)

  # combine with df_stomata, select only locations and area
  tess_area = cbind(voro_cell_df, df_stomata) |>
    dplyr::select(x, y, area)

  # make area a numeric
  tess_area$area = as.numeric(tess_area$area)

  # add reject column
  tess_area$reject = reject

  # return result
  return(tess_area)

}

# a function that performs tessellation for every leaf surface in raw data

tessellate_all = function(
  df_stomata,
  pixels_x,
  pixels_y,
  length_cutoff
){

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 5L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y", "trt_leaf_number", "image_number")
                           , colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_integerish(length_cutoff, lower =  1L, len = 1L)

  # group by surface, leaf number, image number
  stom_group = df_stomata |> group_by(trt_leaf_number, image_number)

  # run function over every group of the grouped dataset
  stom_tess = stom_group |>
    dplyr::group_modify(~ tessellate_stomata(.x, pixels_x = pixels_x, pixels_y = pixels_y))

  # merge dataframes
  df_merge = merge(stom_tess, df_stomata)

  # define NAs based on reject logical
  tessout1 = df_merge

  rejcol = grep("reject", names(df_merge))
  areacol = grep("area", names(df_merge))

  tessout1[areacol][df_merge[rejcol]] <- NA

  # define NAs based on cut off value of 5 (from distribution of length values)
  is.na(tessout1[["length"]]) = tessout1[["length"]] < length_cutoff

  # return final dataframe with all values
  return(tessout1)

}

# a function that performs tessellate_stomata() over every element of a list

tessellate_stomata_list = function(
  list, # list of x, y stomatal grids
  pixels_x, # number of horizontal pixels
  pixels_y # number of vertical pixels
){

  # checks
  checkmate::assert_list(list)
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)

  # lapply
  random_tess = lapply(list[[1]], tessellate_stomata, pixels_x, pixels_y)

  equal_tess = lapply(list[[2]], tessellate_stomata, pixels_x, pixels_y)

  # make list
  list = list(random_tess, equal_tess)

  # return
  return(list)

}

# a function that finds the mean, standard deviation, median absolute variation,
# standard error of tessellation area values

tessellation_summary = function(
  tessellate_out # output of tessellation function
){

  # checks
  checkmate::assert_data_frame(tessellate_out, min.cols = 2L, any.missing = TRUE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y", "area", "reject"), colnames(tessellate_out))

  # filter out NAs and TRUEs
  filtered_tess = tessellate_out |> dplyr::filter(reject == FALSE)
  tess_final = na.omit(filtered_tess)

  # calculate mean, standard deviation, and standard error and sample number
  number = nrow(tess_final)
  mean = mean(tess_final$area)
  median = median(tess_final$area)
  sd = sd(tess_final$area)
  se = sd(tess_final$area)/sqrt(number)
  mad = mad(tess_final$area)

  # make list
  list = list(mean, median, mad, sd, se, number)

  # return list
  return(list)

}

# a function that performs tessellation_summary over every element of list

tessellation_summaries = function(
  list # list of tessellate_stomata outputs
){

  # checks
  checkmate::assert_list(list)

  # lapply
  random_summ = lapply(list[[1]], tessellation_summary)

  equal_summ = lapply(list[[2]], tessellation_summary)

  # make a list
  list = list(random_summ, equal_summ)

  # return
  return(list)

}

# a function that build a dataframe from a list of tessellation_summaries

tess_list_to_frame = function(
  tessellation_summaries_out # output from tessellation summaries
){

  # checks
  checkmate::assert_list(tessellation_summaries_out)

  # build a dataframe
  frame = data.frame(matrix(unlist(tessellation_summaries_out),
                                     nrow = sum(lengths(tessellation_summaries_out)),
                                     byrow = TRUE), stringsAsFactors = FALSE) %>%
    rename(
      mean = X1,
      median = X2,
      mad = X3,
      sd = X4,
      se = X5,
      number = X6
    ) %>%
    dplyr::mutate(type = c(rep("random", length(tessellation_summaries_out[[1]])),
                   rep("equal", length(tessellation_summaries_out[[2]]))))

  # remove NAs caused by a single stomate
  frame = na.omit(frame)

  # return dataframe
  return(frame)

}

# a function that plots tessellation areas

plot_tessellation = function(
  df_stomata, # a dataframe of stomatal positions (x, y)
  pixels_x, # number of horizontal pixels
  pixels_y, # number of vertical pixels
  real, # logical; if FALSE, image used, if TRUE,
       # ellipse used
  image_path = "objects/open-stomate.png"
){

  # checks
  checkmate::assert_tibble(df_stomata, min.cols = 2L, any.missing = FALSE,
                           col.names = "named")
  checkmate::assert_subset(c("x", "y"), colnames(df_stomata))
  checkmate::assert_integerish(pixels_x, lower =  1L, len = 1L)
  checkmate::assert_integerish(pixels_y, lower =  1L, len = 1L)
  checkmate::assert_flag(real)
  checkmate::assertTRUE(file.exists(image_path))

  # tessellate to calculate area
  tess = tessellate_stomata(df_stomata, pixels_x, pixels_y) |>
    mutate(image = image_path) |>
    mutate(area = ifelse(test = reject, NA, area))

  # plot
  # define bounding box
  box = data.frame(x = c(0, pixels_x, pixels_x, 0),
                   y = c(0, 0, pixels_y, pixels_y))
  if (!real) {

    tess = tess |>
      mutate(angle = runif(n(), 0, 2 * pi), length = 2 * pixels_x / 35)

  } else {

    tess = full_join(tess, df_stomata, by = c("x", "y")) |>
      mutate(angle = pi * angle / 180)
  }

  plot = ggplot(tess, aes(x, y, fill = area)) +
    geom_voronoi(color = "black", outline = box) +
    geom_ellipse(aes(x0 = x, y0 = y, a = length / 2, b = length / 4,
                     angle = angle)) +
    geom_ellipse(aes(x0 = x, y0 = y, a = length / 10, b = length / 20,
                     angle = angle)) +
    # geom_image(mapping = aes(image = image), size = 0.06) +
    scale_fill_gradient(low = "#FFFF00", high = "#FF0000", na.value = "white",
                        name= "zone area") +
    geom_rect(aes(xmin = 0, xmax = pixels_x, ymin = 0, ymax = pixels_y),
              fill = NA, color = "black") +
    #ggtitle("Zone defined by each stomate") +
    #coord_cartesian(xlim = c(-1, 513)) +
    coord_equal() +
    theme_void()

  # return plot
  return(plot)

}


## Functions for diffusion model ----

### Helper functions ----
calc_flux = function(C_ias, parms, stomatal_distribution) {

  sr = match.arg(stomatal_distribution, c("amphi", "hypo"))

  # assume C_ias is vector with length n_x * n_y
  C_ias = matrix(C_ias, nrow = parms[["n_y"]], ncol = parms[["n_x"]])
  dC = numeric(length(C_ias))

  if (sr == "amphi") {
    # need to change boundary condition for offset
    boundary = matrix(c(parms[["C_s"]], C_ias[1, 2:parms[["n_x"]]]),
                      nrow = 1, ncol = parms[["n_x"]])

    D_e = matrix(parms[["Dz"]] * (parms[["Intpor"]] / parms[["tort"]]),
                 nrow = parms[["n_y"]] + 1, ncol = parms[["n_x"]])
    flux1 = -D_e * rbind(
      C_ias[1, ] - boundary,
      C_ias[2:parms[["n_y"]], ] - C_ias[1:(parms[["n_y"]] - 1), ],
      rev(boundary) - C_ias[parms[["n_y"]], ]
    ) / parms[["t_elem"]]

    dC = dC - (flux1[2:(parms[["n_y"]] + 1), ] - flux1[1:parms[["n_y"]], ]) / parms[["t_elem"]]

    D_e = matrix(parms[["Dz"]] * (parms[["por"]] / parms[["tort"]]),
                 nrow = parms[["n_y"]], ncol = parms[["n_x"]] + 1)
    flux2 = -D_e * cbind(
      rep(0, parms[["n_y"]]),
      C_ias[, 2:parms[["n_x"]]] - C_ias[, 1:(parms[["n_x"]] - 1)],
      rep(0, parms[["n_y"]])
    ) / parms[["t_elem"]]

    dC = dC - (flux2[, 2:(parms[["n_x"]] + 1)] - flux2[, 1:parms[["n_x"]]]) / parms[["t_elem"]]

  }

  if (sr == "hypo") {

    # need to change boundary condition for offset
    boundary = matrix(c(parms[["C_s"]], C_ias[1, 2:(parms[["n_x"]] - 1)],
                      parms[["C_s"]]), nrow = 1, ncol = parms[["n_x"]])

    D_e = matrix(parms[["Dz"]] * (parms[["Intpor"]] / parms[["tort"]]),
                 nrow = parms[["n_y"]] + 1, ncol = parms[["n_x"]])

    flux1 = -D_e * rbind(
      C_ias[1, ] - boundary,
      C_ias[2:parms[["n_y"]], ] - C_ias[1:(parms[["n_y"]] - 1), ],
      rep(0, parms[["n_x"]])
    ) / parms[["t_elem"]]

    dC = dC - (flux1[2:(parms[["n_y"]] + 1), ] - flux1[1:parms[["n_y"]], ]) / parms[["t_elem"]]

    D_e = matrix(parms[["Dz"]] * (parms[["por"]] / parms[["tort"]]),
                 nrow = parms[["n_y"]], ncol = parms[["n_x"]] + 1)
    flux2 = -D_e * cbind(
      rep(0, parms[["n_y"]]),
      C_ias[, 2:parms[["n_x"]]] - C_ias[, 1:(parms[["n_x"]] - 1)],
      rep(0, parms[["n_y"]])
    ) / parms[["t_elem"]]

    dC = dC - (flux2[, 2:(parms[["n_x"]] + 1)] - flux2[, 1:parms[["n_x"]]]) / parms[["t_elem"]]

  }

  dC
}

calc_An = function(C_liq, parms) {
  pmin((parms[["k_c"]] * parms[["X_c"]] * C_liq) / (parms[["K_m"]] + C_liq),
       parms[["n_y"]] * (C_liq * parms[["j_e"]] / (4 * C_liq + 8 * parms[["Gamma"]])))
}

calc_light <- function(parms) {

  # depth: depth, 0 is abaxial surface [m]
  # t_leaf: leaf thickness [m]

  n = parms[["n_y"]]
  depth = seq(0, parms[["t_leaf"]], by = parms[["t_elem"]])
  cum_irradiance = pexp(depth, parms[["k_i"]])
  I_prop = cum_irradiance[2:(n + 1)] - cum_irradiance[1:n]
  rev(I_prop) # order so that first element for abaxial surface

}

# Reaction-Diffusion function ----
diffusion2D <- function(t, state, parms) {

  n = parms[["n_x"]] * parms[["n_y"]]
  C_ias <- state[seq_len(n)] # Define vapor [CO2]
  C_liq <- state[n + seq_len(n)] # Define liquid [CO2]

  # Calculate change in vapor [CO2] over time;
  # set constant [CO2] = C_s at lower/upper boundary of leaf
  flux_cias <- calc_flux(C_ias, parms)

  # Calculate the divergence of the C_ias flux
  # original eqn divided by porosity, but I don't understand why. We think we figured it out.
  dC_ias <- (flux_cias + parms[["g_liq"]] * (C_liq - C_ias) / parms[["t_elem"]]) / parms[["por"]]

  # Calculate carboxylation rate
  An <- calc_An(C_liq, parms)

  Rp <- (An * parms[["Gamma"]]) / C_liq # Calculate oxygenation rate

  # Calculate the corresponding change in C_liq flux
  dC_liq <- (parms[["g_liq"]] * (C_ias - C_liq) / parms[["t_elem"]] - An + Rp + parms[["Rd"]]) * parms[["t_leaf"]] / sum(parms[["Sm"]]) / parms[["Vstrom"]]

  return(list(c(dC_ias, dC_liq))) # return output from function

}

# function to run 2D CO2 diffusion model

co2d = function(
  t_stomate, # distance between stomata [m]
  t_leaf, # leaf thickness excluding cuticle [m]
  t_elem, # length of element [m],
  g_sc, # stomatal conductance to CO2 [m/s]
  stomata_offset, # are stomata offset?
  co2d_pars # list of parameters
) {

  # For testing
  # t_stomate = 100 / 1e6 # [m]
  # t_leaf = 200 / 1e6 # [m]
  # t_elem = 1 / 1e6 # [m]
  # g_sc = 0.005 # [m/s] # convert to mol/m^2/s by multiplying P / Temp / R_gas. Yields g_sc = 0.2 mol/m^2/s at 101 kPa and 25 C.
  # offset = TRUE
  # co2d_pars = make_co2d_pars()

  # make variables into list ----
  co2d_vars = list(
    t_stomate = t_stomate,
    t_leaf = t_leaf,
    t_elem = t_elem,
    g_sc = g_sc
  )

  # checks ----
  # should there be a warning if variables do not include units? currently, assume units are um
  # checkmate::assert_flag(stomata_offset)
  # co2d_vars = check_variables(co2d_vars)
  # co2d_pars = check_parameters(co2d_pars)

  # additional calculations ----
  co2d_all = make_co2d_all(co2d_vars, co2d_pars)

  # run ----
  co2d_out = steady.2D(
    y = initialize_grid(co2d_all),
    func = diffusion2D,
    parms = co2d_all,
    dimens = c(co2d_all[["n_y"]], co2d_all[["n_x"]]),
    nspec = 2,
    # stomata_offset = stomata_offset, # ignored right now
    method = "runsteady",
    lrw = 2e3 * co2d_all[["n_y"]] * co2d_all[["n_x"]]
  )

  ## check results
  # need to do
  ## summarize
  co2d_out$summary = summarize_co2d(co2d_out, co2d_all)

  co2d_out

}

### summarize output ----
summarize_co2d = function(co2d_out, co2d_all) {

  # Leaf-level photosynthesis ----
  n = prod(attr(co2d_out, "dimens"))
  C_ias = co2d_out$y[seq_len(n)]
  C_liq = co2d_out$y[n + seq_len(n)]
  An = calc_An(C_liq, co2d_all)
  Rp = (An * co2d_all[["Gamma"]]) / C_liq

  # Multiply by 1e6 to convert to [umol m-2 s-1]
  An_area = 1e6 / co2d_all[["n_x"]] * sum((An - Rp - co2d_all[["Rd"]]) * co2d_all[["Vstrom"]] * co2d_all[["Sm"]])

  list(C_ias = C_ias, C_liq = C_liq, An_area = An_area)

}

### making functions ----
make_co2d_pars  = function() {

  # We want to change this to allow customization
  parms = list(

    P = 101.325, # air pressure at sea level [kPa]
    temp = 298.15, # assume constant temperature [K]
    R_gas = 8.314, # ideal gas constant [J/K/mol]

    # FYI, how to calculate for C_a and O for different P or temp
    # 21% O2
    # set_units(0.21 * P / (R_gas * temp), mol/m^3)
    # 415 ppm
    # set_units((415/1e6) * P / (R_gas * temp), mmol/m^3)

    # Environmental
    C_a = 1.72e-2, # CO2 concentration in air from Earles et al. [mol/m^3]
    # O = set_units(8.584027, mol/m^3),

    Dz = 1.54e-5, # Diffusivity of CO2 in air [m2 s-1]
    Intpor = 0.2, # Porosity at interface [m3 m-3]
    por_spg = 0.3, # Porosity of the spongy mesophyll [m3 air m3 leaf]
    por_pal = 0.1, # Porosity of the palisade mesophyll [m3 air m3 leaf]
    tort = 2, # Tortuosity of the palisade and spongy mesophyll [m m-1]
    C_s = 0.85 * 1.72e-2, # CO2 concentration at stomate [mol / m^3]
    k_c = 3, # Rubisco turnover rate [s-1]
    X_c = 2.5, # Rubisco concentration [mol m-3]
    K_m = 18.7e-3, # Rubisco effective Michaelis-Menten constant [mol m-3]
    Gamma = 1.75e-3, # CO2 compensation point [mol m-3]
    Theta = 1, # Curvature factor in FvCB model
    Rd = 0.066, # Dark respiratory rate [mol m-3 s-1]
    Jmax = 275e-6, # [mol/m^2/s]
    Beta = 0.44, # Fraction of light absorbed by PSII [mol mol-1]
    Alpha = 0.72, # Leaf level absorption [mol mol-1]
    phiPSII = 0.85, # quantum efficiency of PSII
    I_0 = 1e-3, # Incident irradiance [mol m-2 s-1]
    # k_i so 90% absorption for 300 um
    k_i = 7675.284, # light extinction coefficient [units]
    frac_pal = 0.6, # Fraction palisade mesophyll
    Sm_spg = 6.5, # Sm spongy mesophyll [m2 m-2]
    Sm_pal = 40,  # Sm palisade mesophyll [m2 m-2]
    Vstrom = 1.74e-6, # Stroma volume per mesophyll surface area [m3 m-2]
    Vmito = 0.27e-7, # Mitochondrial volume per mesophyll surface area [m3 m-2]
    Sm_std = 30, # Sm at which assumed J_max occurs
    g_liq = 0.25e-3 # Cell wall + liquid conductivity into stroma [m s-1]
  )

  parms$frac_spg = 1 - parms[["frac_pal"]] # Fraction spongy mesophyll

  parms

}

# Use variable and parameter information to complete calculations
make_co2d_all = function(co2d_vars, co2d_pars) {

  parms = c(co2d_vars, co2d_pars)

  # dimensions
  parms$n_x = round(parms[["t_stomate"]] / 2 / parms[["t_elem"]])
  parms$n_y = round(parms[["t_leaf"]] / parms[["t_elem"]])
  parms$n_y_pal = round(parms[["n_y"]] * parms[["frac_pal"]])
  parms$n_y_spg = parms[["n_y"]] - parms[["n_y_pal"]]

  # element-wise porosity and Sm
  parms$Sm = rep(c(parms[["Sm_spg"]], parms[["Sm_pal"]]),
                 c(parms[["n_y_spg"]], parms[["n_y_pal"]])) / parms[["n_y"]]

  parms$por = rep(c(parms[["por_spg"]], parms[["por_pal"]]),
                  c(parms[["n_y_spg"]], parms[["n_y_pal"]]))
  parms$Intpor = c(parms[["por"]], dplyr::last(parms[["por"]]))

  # Light
  parms$I_i = calc_light(parms)

  # Define electron transport rate function
  # Proportionality factor between mesophyll surface area and maximum e- transport rate
  parms$Jprop = parms[["Jmax"]] / parms[["Sm_std"]]
  # Maximum e- transfer rate on a leaf area basis [mol m-2 s-1]
  parms$J_max = parms[["Jprop"]] * sum(parms[["Sm"]])
  # Electron transport rate on a stroma volume basis at position z [mol m-3 s-1]
  parms$j_max_z = parms[["J_max"]] / sum(parms[["Sm"]]) / parms[["Vstrom"]]
  # Scenario 2: Jmax proportional to chloroplast distribution
  parms$j_max = parms[["j_max_z"]] * parms[["I_i"]] / sum(parms[["I_i"]])

  parms$j_inf =
    parms[["I_0"]] * parms[["Beta"]] * parms[["Alpha"]] * parms[["phiPSII"]] *
    (parms[["I_i"]] / sum(parms[["I_i"]])) /
    sum(parms[["Sm"]]) / parms[["Vstrom"]]

  parms$j_e = pmin(parms[["j_inf"]], parms[["j_max"]])

  parms

}

initialize_grid = function(co2d_all) {

  # CO2 concentration [mol / m^3]
  # For initial condition, assume C_ias / C_a = 0.5; C_liq / C_ias = 0.5
  rep(
    c(0.5, 0.25) * co2d_all[["C_a"]],
    each = co2d_all[["n_x"]] * co2d_all[["n_y"]]
  )

}

### checking functions ----
check_variables = function(co2d_vars) {
  checkmate::assert_list(co2d_vars)
  co2d_vars = within(co2d_vars, {
    interstomatal_distance = set_units(interstomatal_distance, um) |>
      drop_units()
    leaf_thickness = set_units(leaf_thickness, um) |>
      drop_units()
    node_length = set_units(node_length, um) |>
      drop_units()
    g_sc = set_units(g_sc, m/s) |>
      drop_units()
  })
  checkmate::assert_true(all(purrr::map_lgl(co2d_vars, is.numeric)))
  checkmate::assert_true(all(purrr::map_lgl(co2d_vars, ~ {length(.x) == 1L})))

  checkmate::assert_true(
    (co2d_vars$interstomatal_distance / 2) %% co2d_vars$node_length == 0,
    .var.name = "interstomatal_distance / 2 is divisible by node_length"
  )
  checkmate::assert_true(
    co2d_vars$leaf_thickness %% co2d_vars$node_length == 0,
    .var.name = "leaf_thickness is divisible by node_length"
  )

  within(co2d_vars, {
    n_row = leaf_thickness / node_length
    n_col = (interstomatal_distance / 2) / node_length
  })

}

check_parameters = function(co2d_pars) {
  checkmate::assert_list(co2d_pars)
  # no checks yet
  co2d_pars
}



# TESTING OUT FUNCTION TO DRAW LEAF FOR MODEL EXPLANATION

plot_photo_2d_pm = function(parms, soln = NULL, ...) {

  # scaling ratio of U / 2 to T_leaf
  # assumes that dx == dz
  xz_ratio = parms[["n_x"]] / parms[["n_z"]]

  df_stomata = crossing(
    nesting(x = c(0, 1), y = c(0, 1) / xz_ratio),
    r = 0.04,
    x_offset = c(-1.5, 1.5)
  ) |>
    mutate(x0 = x + x_offset * r)

  df_epidermis = df_stomata |>
    summarise(stomata_min = min(x0 - r), stomata_max = max(x0 + r), .by = "y") |>
    mutate(
      xmin = ifelse(stomata_min < 0, stomata_max, 0),
      xmax = ifelse(stomata_max > 1, stomata_min, 1)
    )

  ggplot() +
    geom_circle(data = df_stomata, mapping = aes(x0 = x0, y0 = y, r = r)) +
    geom_segment(data = df_epidermis,
                 mapping = aes(
                   x = xmin,
                   xend = xmax,
                   y = y,
                   yend = y
                 ), linewidth = 1.5) +
    geom_text(
      data = df_stomata |>
        select(x, y) |>
        distinct(),
      mapping = aes(x, y, label = "stomate"),
      position = position_nudge(y = c(-0.05, 0.05) / xz_ratio)
    ) +
    geom_text(
      mapping = aes(x = 0.5, y = c(-0.1, 0.5, 1.1) / xz_ratio, label = c("abaxial", "mesophyll", "adaxial"))
    ) +
    geom_segment(mapping = aes(x = 0, y = -0.15 / xz_ratio, xend = 1, yend = -0.15 / xz_ratio),
                 arrow = arrow(angle = 90, length = unit(0.025 / xz_ratio, "npc"), ends = "both") ) +
    geom_segment(mapping = aes(x = -0.15, y = 0 / xz_ratio, xend = -0.15, yend = 1 / xz_ratio),
                 arrow = arrow(angle = 90, length = unit(0.025 / xz_ratio, "npc"), ends = "both") ) +
    geom_text(mapping = aes(x = 0.5, y = -0.2 / xz_ratio, label = "paste(italic(U) / 2, ', half the interstomatal distance [', mu, 'm]')"),
              parse = TRUE) +
    geom_text(mapping = aes(x = -0.2, y = 0.5 / xz_ratio, label = "paste(italic(T)[leaf], ', leaf thickness [', mu, 'm]')"),
              parse = TRUE, angle = 90) +
    coord_equal() +
    theme_void()

}

# Convenience function for reporting F-statistics
report_fstat = function(m, x) {

  assert_class(m, c("aov", "lm"))
  s = summary(m)
  r = str_detect(rownames(s[[1]]), glue("^{x}\\s*$"))
  assert_true(sum(r) == 1)

  glue(
    "F_{{{df1},{df2}}} = {Fstat}, P = {pval}",
    df1 = s[[1]][r, "Df"],
    df2 = m$df.residual,
    Fstat = scientize(s[[1]][r, "F value"]),
    pval = scientize(s[[1]][r, "Pr(>F)"])
  )

}
