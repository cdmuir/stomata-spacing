source("r/header.R")

# load processed data
stomata = readr::read_rds("objects/stomata_position.rds") |>
  dplyr::group_by(trt_leaf_number, image_number1)

# load synthetic data and group by
synthetic_data = readr::read_rds("objects/synthetic_data.rds") |>
  # add image_number1 to synthetic data for abaxial, adaxial pairing
  dplyr::mutate(image_number1 = stringr::str_extract(image_number, "[0-9]{1}$")) |>
  # group synthetic data
  dplyr::group_by(treatment, leaf_number, image_number1, sim) |>
  # exclude equal grids - only concerned with random grids
  dplyr::filter(grid == "random")

# run correlation on real data
leaf = create_raster(pixels_x = 512, pixels_y = 512)
corr_observed = stomata |>
  purrrlyr::by_slice(amphi_corr, leaf = leaf, .collate = "rows",
                     .to = "corr_coefficient")

# run correlation on synthetic data
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
  total = nrow(attr(synthetic_data, "groups"))
)

corr_synthetic = synthetic_data |>
  purrrlyr::by_slice(amphi_corr, leaf = leaf, progress = TRUE,
                     .collate = "rows", .to = "corr_coefficient")

# write rds files
readr::write_rds(corr_observed, "objects/corr_observed.rds")
readr::write_rds(corr_synthetic, "objects/corr_synthetic.rds")
