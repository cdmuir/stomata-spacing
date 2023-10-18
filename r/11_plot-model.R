source("r/header.R")
source("r/photo_2d_pm.R")

ph2d = read_rds("objects/model_output.rds")

vars = c("stomatal_arrangement", "n_z", "n_x", "I_0", "A_n", "phi_pal")

ph2d_vars = ph2d |>
  map("parms") |>
  map_dfr(\(.x) as_tibble(.x[vars]))

tikz("ms/figures/model_summary.tex", standAlone = TRUE, width = 5, height = 5)
options(tikzLatexPackages = c(getOption("tikzLatexPackages")))

ph2d_vars |>
  pivot_wider(names_from = "stomatal_arrangement", values_from = "A_n") |>
  # coordination advantage
  mutate(
    ca = log(offset / aligned),
    U = glue("$U = {U}~\\mu \\mathrm{{m}}$", U = n_x * 2),
    T_leaf = glue("$T_\\mathrm{{leaf}} = {T_leaf}~\\mu \\mathrm{{m}}$", T_leaf = n_z - 1),
    U = fct_reorder(U, n_x),
    T_leaf = fct_reorder(T_leaf, n_z),
    `$\\varphi_\\mathrm{pal}$` = as.factor(phi_pal)
  ) |>
  ggplot(aes(1e6 * I_0, ca, color = `$\\varphi_\\mathrm{pal}$`)) +
  facet_grid(T_leaf ~ U) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1000), breaks = c(0, 500, 1000)) +
  xlab("$I_0~[\\mu \\mathrm{mol}~\\mathrm{m}^{-2}~\\mathrm{s}^{-1}]$") +
  ylab("coordination advantage") +
  theme(legend.position = "top")

dev.off()
system("pdflatex -output-directory=ms/figures/ ms/figures/model_summary.tex")

file.remove(paste0("ms/figures/model_summary.", c("aux", "log")))
