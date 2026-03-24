source("R/model_utils.R", local = FALSE)

calibration_example <- list(
  parameters = default_parameters(),
  initial_state = default_state(),
  observation_window = 0:5
)

save(
  calibration_example,
  file = "data/calibration_example.rda",
  compress = "bzip2"
)
