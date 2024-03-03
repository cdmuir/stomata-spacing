source("r/header.R")

stomata_position = readr::read_csv("raw-data/stomata_position.csv",
                                   show_col_types = FALSE) |>

  # add light column
  dplyr::mutate(
    light = dplyr::case_when(
      str_detect(treatment, "HL") ~ "high",
      str_detect(treatment, "LL") ~ "low"
    ),
    light = tidyr::replace_na(light, "medium")
  ) |>

  # modify leaf_number
  tidyr::unite(trt_leaf_number, treatment, leaf_number) |>
  dplyr::group_by(trt_leaf_number, image_number) |>
  dplyr::mutate(n_stomata = dplyr::n())

# add image_number1
stomata_position$image_number1 = as.character(gsub("^.{0,3}", "", stomata_position$image_number))

# Check for exact duplicates
stomata_position |>
  dplyr::group_by(trt_leaf_number, image_number, surface) |>
  dplyr::reframe(is_dup = duplicated(tibble(x, y))) |>
  dplyr::filter(is_dup)

# identify near duplicates
stomata_position |>
  # group by
  dplyr::group_by(trt_leaf_number, image_number, surface) |>
  #first, arrange by theta
  arrange(x) |>
  #mark rows with theta <=5 of next row
  mutate(x_mark = if_else(abs(x - lag(x)) <= 7, 1, 0)) |>
  #mark rows with rho <=15 of next row
  mutate(y_mark = if_else(abs(y - lag(y)) <= 7, 1, 0)) |>
  #replace NA's with 0
  tidyr::replace_na(list(x_mark = 0, y_mark = 0)) |>
  #filter all unwanted rows
  filter(x_mark == 1 & y_mark == 1) |>
  #drop unwanted columns
  select(-x_mark, -y_mark)

# add tessellation info
stomata_tess = stomata_position |>
  dplyr::group_modify(~ tessellate_stomata(.x, pixels_x = 512, pixels_y = 512)) |>
  dplyr::full_join(stomata_position,
                   by = c("trt_leaf_number", "image_number", "x", "y")) |>
  dplyr::mutate(area = ifelse(test = reject, NA, area))

readr::write_rds(stomata_position, "objects/stomata_position.rds")
readr::write_rds(stomata_tess, "objects/stomata_position_length_area.rds")
