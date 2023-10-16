# testing code to make table of 2d fem model parameters
source("r/header.R")

model_output = read_rds("objects/model_output.rds") |>
  map_dfr(\(.x) {
    l = map(.x$parms, length)
    .x$parms[l == 1] |>
      as_tibble()
  })

pm_pars = readr::read_csv("raw-data/2d-pm-parameters.csv", show_col_types = FALSE)

# Select variables
## Find all elements that vary. It should include desired variables and
## and additional response variables and redundant variables (e.g. T_leaf is
## derived from n_z). The last step selects variables for table.
focal_vars = c("I_0", "phi_pal", "T_leaf", "U")
focal_vars_units = c(
  I_0 = "$\\mu$mol m$^{-2}$ s$^{-1}$",
  phi_pal = "m$^3$ airspace m$^{-3}$ leaf",
  T_leaf = "$\\mu$m",
  U = "$\\mu$m"
)


vars = model_output |>
  as.list() |>
  sapply(\(.x) length(unique(.x))) |>
  sapply(\(.x) .x > 1) %>%
  and(names(.) %in% focal_vars)

# Variable transformations
var_transformation = list(
  # convert from mol to umol
  I_0 = function(.x) {.x * 1e6},
  # no conversion
  phi_pal = function(.x) .x,
  # convert from m to um
  T_leaf = function(.x) {.x * 1e6},
  # convert from m to um
  U = function(.x) {.x * 1e6}
)

vars |>
  which() |>
  names() %>%
  glue("model_output[, '{var}'] %<>% var_transformation${var}()", var = .) %>%
  parse(text = .) |>
  eval()

model_output[, vars] |>
  as.list() |>
  map_dfr(\(.x) {
    glue("${min}-{max}$", min = min(.x), max = max(.x))
  }) |>
  pivot_longer(
    everything(),
    names_to = "r",
    values_to = "Parameter range"
  ) |>
  left_join(select(pm_pars, r, Variable = symbol), by = join_by(r)) |>
  full_join(
    tibble(r = focal_vars, Units = focal_vars_units[focal_vars]),
    by = join_by(r)
  ) |>
  select(Variable, `Parameter range`, Units) |>
  write_rds("objects/model_var_table.rds")


