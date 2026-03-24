## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(twostraindengue)

parameters <- default_parameters()
state <- default_state(parameters)
time_grid <- seq(0, 2, by = 0.1)

model_output <- run_model(state, time_grid, parameters)
head(extract_reported_cases(model_output))

