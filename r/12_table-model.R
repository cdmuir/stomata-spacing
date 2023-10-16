# testing code to make table of 2d fem model parameters
source("r/header.R")

model_output = read_rds("objects/model_output.rds") |>
  map_dfr(\(.x) {
    l = map(.x$parms, length)
    .x$parms[l == 1] |>
      as_tibble()
  })

# Select variables
vars = model_output |>
  as.list() |>
  sapply(\(.x) length(unique(.x))) |>
  sapply(\(.x) .x > 1)

# Variable transformations
var_transformation = list(
  # convert from mol to umol
  I_0 = function(.x) {.x * 1e6},

)
model_output[, vars] |>
  as.list() |>
  map(\(.x) {
    glue::glue("${min}--{max}$", min = min(.x), max = max(.x))
  })

