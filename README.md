# twostraindengue

`twostraindengue` is a research-focused R package for simulating and calibrating
a two-strain dengue transmission model. This repository includes reusable model
code, raw and cleaned data assets, reproducible analysis scripts, tests, and a
pkgdown-ready documentation layout.

## Core Workflow

```r
library(twostraindengue)

parameters <- default_parameters()
initial_state <- default_state(parameters)
time_grid <- seq(0, 6, by = 0.1)

model_output <- run_model(initial_state, time_grid, parameters)
reported_cases <- extract_reported_cases(model_output)

study_window <- load_study_data()
analysis <- run_hermosillo_2010_analysis(
  config = list(outputs = list(write_outputs = FALSE, show_progress = FALSE))
)
```

## Repository Layout

- `R/`: package functions
- `data/`: packaged example datasets
- `inst/extdata/raw/`: raw study inputs
- `inst/scripts/`: reproducible analysis scripts
- `outputs/`: generated figures and calibration tables
- `vignettes/`: long-form analysis and package documentation
- `tests/testthat/`: automated tests
