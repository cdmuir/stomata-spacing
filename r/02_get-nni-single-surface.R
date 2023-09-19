source("r/header.R")

stomata = readr::read_csv("raw-data/stomata_position.csv") |>
  dplyr::group_by(treatment, leaf_number, image_number, surface)

synthetic_data = readr::read_rds("objects/synthetic_data.rds") |>
  dplyr::group_by(treatment, leaf_number, image_number, surface, sim, grid)

nni_observed <- stomata |>
  purrrlyr::by_slice(get_nni, pixels_x = 512, pixels_y = 512, .collate = "rows")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
  total = nrow(attr(synthetic_data, "groups"))
)

nni_synthetic <- synthetic_data |>
  purrrlyr::by_slice(get_nni, pixels_x = 512, pixels_y = 512, progress = TRUE,
                     .to = "NNI_synthetic")

readr::write_rds(nni_observed, "objects/nni_observed.rds")
readr::write_rds(nni_synthetic, "objects/nni_synthetic.rds")
