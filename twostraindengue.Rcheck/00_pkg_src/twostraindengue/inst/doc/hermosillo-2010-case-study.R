## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(twostraindengue)

analysis <- run_hermosillo_2010_analysis(
  observed_rows = 1:3,
  full_time_grid = seq(0, 5, by = 0.5),
  calibration_grid = seq(0, 2, by = 0.5),
  fit_time_grid = seq(0, 5, by = 0.5),
  beta_h_values = seq(1, 1.4, length.out = 3),
  beta_m_values = seq(0.001, 0.003, length.out = 3),
  r01_values = seq(1.7, 1.9, length.out = 3),
  r02_values = seq(0.003, 0.0034, length.out = 3),
  r0_values = seq(1.3, 1.4, length.out = 3),
  write_outputs = FALSE
)

analysis$reproduction_numbers

