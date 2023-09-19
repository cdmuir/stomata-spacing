source("r/header.R")

nni_observed = readr::read_rds("objects/nni_observed.rds") |>
  dplyr::select(treatment, leaf_number, image_number, surface, n_stomata,
                NNI_observed = NNI)

nni_synthetic = readr::read_rds("objects/nni_synthetic.rds") |>
  dplyr::select(treatment, leaf_number, image_number, surface, #n_stomata,
                sim, grid, NNI_synthetic)

stomata_area_summ = readr::read_rds("objects/area_summ.rds")


single_surface_results <- nni_observed |>
  dplyr::rowwise() |>
  purrr::pmap_dfr(~ {

    sim_ran = nni_synthetic |>
      dplyr::filter(
        treatment == ..1,
        leaf_number == ..2,
        image_number == ..3,
        surface == ..4,
        grid == "random"
      ) |>
      dplyr::pull(NNI_synthetic) |>
      purrr::map(magrittr::extract, "NNI") |>
      purrr::map_dbl(unlist)

    sim_equ = nni_synthetic |>
      dplyr::filter(
        treatment == ..1,
        leaf_number == ..2,
        image_number == ..3,
        surface == ..4,
        grid == "equal"
      ) |>
      dplyr::pull(NNI_synthetic) |>
      purrr::map(magrittr::extract, "NNI") |>
      purrr::map_dbl(unlist)

    ret = tibble(
      treatment = ..1,
      leaf_number = ..2,
      image_number = ..3,
      surface = ..4,
      n_stomata = ..5,
      NNI_observed = ..6,
      p_value = length(which(sim_ran >= ..6)) / length(sim_ran),
      dispersion = (..6 - median(sim_ran)) / (median(sim_equ - median(sim_ran)))
    ) #|>
      #full_join(as.tibble(stomata_area_summ), by = c("treatment", "leaf_number", "image_number", "surface"))
      #cbind(stomata_area_summ)

    # I'm trying to merge these datasets using cbind and full_join. but I really
    # don't know why it's not working.

    return(ret)

  }) |>
  # add light column
  dplyr::mutate(
    light = dplyr::case_when(
      str_detect(treatment, "HL") ~ "high",
      str_detect(treatment, "LL") ~ "low"
    ),
    light = tidyr::replace_na(light, "medium"))

readr::write_rds(single_surface_results, "objects/single_surface_results.rds")



