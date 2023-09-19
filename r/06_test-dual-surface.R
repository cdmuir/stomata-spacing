source("r/header.R")

# load data
corr_observed = readr::read_rds("objects/corr_observed.rds") |>
  # separate columns
  tidyr::separate(trt_leaf_number, c("treatment", "leaf_number"), sep = "_",
                  remove = TRUE) |>
  dplyr::select(treatment, leaf_number, image_number1, corr = corr_coefficient)

# synthetic data
corr_synthetic = readr::read_rds("objects/corr_synthetic.rds") |>
  dplyr::select(treatment, leaf_number, image_number1, sim,
                corr_sim = corr_coefficient) |>
  dplyr::mutate(leaf_number = as.character(leaf_number))

# make dataframe of results

dual_surface_results <- corr_observed |>
  dplyr::full_join(
    corr_synthetic,
    by = c("treatment", "leaf_number", "image_number1")
  ) |>
  dplyr::group_by(treatment, leaf_number, image_number1) |>
  # two-sided test
  dplyr::mutate(corrsim_ge_corr = abs(corr_sim) >= abs(corr)) |>
  dplyr::summarise(
    corr = dplyr::first(corr),
    mu_corr_sim = mean(corr_sim),
    n_corr_ge_corrsim = sum(corrsim_ge_corr),
    n = n(),
    p_value = n_corr_ge_corrsim / n
  )

readr::write_rds(dual_surface_results, "objects/dual_surface_results.rds")

