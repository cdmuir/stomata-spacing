# woorking on stuff to automate results for paper
source("r/header.R")

dual_surface_results = readr::read_rds("objects/dual_surface_results.rds")

nsd1 = dual_surface_results |>
  ungroup() |>
  dplyr::transmute(sig = p_value < 0.05) |>
  dplyr::summarize(n1 = sum(sig), n2 = sum(!sig))

adjp = mt.rawp2adjp(dual_surface_results$p_value, proc = "BH")

nsd2 = bind_cols(dual_surface_results[adjp$index,], adjp$adjp) |>
  ungroup() |>
  dplyr::transmute(sig = BH < 0.05) |>
  dplyr::summarize(n1 = sum(sig), n2 = sum(!sig))
