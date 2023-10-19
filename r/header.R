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
library(emmeans)
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
library(multtest)
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
library(tikzDevice)
library(tidyr)
library(tripack)
library(units)

theme_set(theme_cowplot())

# Scale for conversion from pixel to um according to image metadata
pixels_per_um = 621.21 / 512

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
