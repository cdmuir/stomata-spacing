# Clear workspace
rm(list = ls())
graphics.off()

# hello world
# Source custom functions
source("r/functions.R")

# Load libraries

library(arulesViz)
library(brms)
library(caret)
library(cmdstanr)
library(checkmate)
library(corrplot)
library(cowplot)
library(dplyr)
library(forcats)
library(furrr)
library(ggdist)
library(ggforce)
library(ggimage)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggvoronoi)
library(glue)
library(Hmisc)
library(magrittr)
library(maptools)
library(metR)
library(mgcv)
library(parallel)
library(progress)
library(progressr)
library(purrr)
library(purrrlyr)
library(RANN)
library(raster)
library(readr)
library(rgdal)
library(rlist)
library(rootSolve)
library(rstatix)
library(sp)
library(spatialEco)
library(spatstat)
library(stringr)
library(stats)
library(terra)
library(tibble)
library(tidyr)
library(tripack)
library(units)

theme_set(theme_cowplot())

# Fit Rubisco model
x_depth = c(145, 830.5)
y_depth = c(0, 700)
x_rubisco = c(558, 56.5)
y_rubisco = c(0, 100)
nishio_carbon_1993_fig5 = read_csv("raw-data/nishio_carbon_1993_fig5.csv",
                                   show_col_types = FALSE) |>
  mutate(
    depth = raw_depth * (diff(y_depth) / diff(x_depth)) - x_depth[1],
    rubisco = raw_rubisco * (diff(y_rubisco) / diff(x_rubisco)) -
      x_rubisco[1] * (diff(y_rubisco) / diff(x_rubisco)),
    # Reversing so that 0 is at abaxial, 1 is at adaxial, to match Earles et al.
    rel_depth = rev(depth / (min(depth) + max(depth)))
  )

fit_rubisco = gam(rubisco ~ s(rel_depth), data = nishio_carbon_1993_fig5)

# Resolve conflicts
# library(conflicted)
# conflict_prefer("abbreviate", "arules")
# conflict_prefer("write", "arules")
# conflict_prefer("intersect", "dplyr")
# conflict_prefer("recode", "dplyr")
# conflict_prefer("setdiff", "dplyr")
# conflict_prefer("setequal", "dplyr")
# conflict_prefer("union", "dplyr")
# conflict_prefer("filter", "dplyr")
# conflict_prefer("lag", "dplyr")
# conflict_prefer("src", "dplyr")
# conflict_prefer("summarize", "dplyr")
# conflict_prefer("select", "dplyr")
# conflict_prefer("combine", "dplyr")
# conflict_prefer("format.pval", "base")
# conflict_prefer("units", "base")
# conflict_prefer("label", "Hmisc")
# conflict_prefer("lift", "purrr")
# conflict_prefer("mask", "raster")
# conflict_prefer("zoom", "raster")
# conflict_prefer("shift", "raster")
# conflict_prefer("rotate", "raster")
# conflict_prefer("cluster", "survival")
# conflict_prefer("concordance", "survival")
# conflict_prefer("getData", "raster")
# conflict_prefer("collapse", "dplyr")
# conflict_prefer("panel.histogram", "lattice")
# conflict_prefer("cells", "terra")
# conflict_prefer("rescale", "terra")
# conflict_prefer("project", "rgdal")
# conflict_prefer("describe", "Hmisc")
# conflict_prefer("size", "arules")

