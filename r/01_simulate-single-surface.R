source("r/header.R")

stomata = readr::read_csv("raw-data/stomata_position.csv",
                          show_col_types = FALSE) |>
  group_by(treatment, leaf_number, image_number, surface)

pb = progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
  total = nrow(attr(stomata, "groups"))
)

set.seed(197789263)
synthetic_data = stomata |>
  purrrlyr::by_slice(
    simulate_grids_from_data,
    pixels_x = 512, pixels_y = 512, n_grid = 1e3, progress = TRUE,
    .collate = "rows"
  )

readr::write_rds(synthetic_data, "objects/synthetic_data.rds")
