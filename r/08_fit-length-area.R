source("r/header.R")

stomata = readr::read_rds("objects/stomata_position_length_area.rds") |>
  dplyr::ungroup() |>
  filter(length > 5,!reject,!is.na(area)) |>
  mutate(
    sqrt_area = sqrt(area),
    image = str_c(trt_leaf_number, image_number),
    sep = "-"
  )

fit = brm(
  bf(
    sqrt_area ~ surface * light * length + (surface * length |
                                              trt_leaf_number) +
      (surface * length | image),
    sigma ~ light
  ),
  data = stomata,
  backend = "cmdstanr",
  chains = 4L,
  cores = 4L,
  seed = 226140189,
  iter = 4e3,
  thin = 2,
  silent = 0
)

write_rds(fit, "objects/fit.rds")
