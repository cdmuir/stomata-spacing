# woorking on stuff to automate results for paper
source("r/header.R")
single_surface_results = read_rds("objects/single_surface_results.rds")

hist(tmp$p_value)
plot(tmp$NNI_observed, log(tmp$p_value))

nni = tmp |>
  dplyr::transmute(sig = p_value < 0.05) |>
  dplyr::summarize(n1 = sum(sig), n2 = sum(!sig))

adjp = mt.rawp2adjp(single_surface_results$p_value, proc = "BH")

bind_cols(single_surface_results[adjp$index,], adjp$adjp)

single_surface_anova = read_rds("objects/single_surface_anova.rds")


single_surface_anova

summary(single_surface_anova)["light", ]

str(single_surface_anova)

single_surface_anova$df.residual

m = single_surface_anova
x = "light"

glue(
  "F_{{{df1},{df2}}} = {Fstat}, P = {pval}",
  df1 = 1,
  df2 = m$df.residual
)
