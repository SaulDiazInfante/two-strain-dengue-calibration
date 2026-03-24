testthat::test_that("parameter vectors preserve the expected order", {
  fixed_parameters <- default_parameters()[fixed_parameter_names]
  parameter_vector <- build_parameter_vector(fixed_parameters, beta_h = 2, beta_m = 0.01)

  testthat::expect_equal(names(parameter_vector), full_parameter_names)
  testthat::expect_equal(unname(parameter_vector[["betaH"]]), 2)
  testthat::expect_equal(unname(parameter_vector[["betaM"]]), 0.01)
})

testthat::test_that("the packaged study data can be loaded with the expected columns", {
  study_data <- load_study_data(rows = 35:40)

  testthat::expect_s3_class(study_data, "data.frame")
  testthat::expect_named(study_data, c("week", "classical", "hemorrhagic"))
  testthat::expect_equal(study_data$week, 35:40)
})

testthat::test_that("the packaged example dataset matches the default study window", {
  utils::data(
    "dengue_hermosillo_2010",
    package = "twostraindengue",
    envir = environment()
  )

  testthat::expect_equal(load_study_data(), dengue_hermosillo_2010)
})

testthat::test_that("observation times map to exact rows in the solver output", {
  model_output <- cbind(
    time = seq(0, 1, by = 0.1),
    matrix(0, nrow = 11, ncol = 2)
  )
  colnames(model_output)[2:3] <- c("z", "Y1h")

  matched_rows <- match_observation_times(model_output, c(0, 0.5, 1))

  testthat::expect_equal(matched_rows, c(1, 6, 11))
})

testthat::test_that("reparameterised likelihoods are consistent with direct betaM evaluation", {
  testthat::skip_if_not_installed("deSolve")

  fixed_parameters <- default_parameters()[fixed_parameter_names]
  initial_state <- default_state(default_parameters())
  observation_times <- 0:5
  time_grid <- seq(0, 5, by = 0.1)
  observed_cases <- as.matrix(load_study_data()[1:6, c("classical", "hemorrhagic")])

  beta_h <- 2.05
  beta_m <- 0.0061
  direct_value <- negative_loglikelihood_original_scale(
    par = c(beta_h, beta_m),
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = time_grid,
    fixed_parameters = fixed_parameters
  )

  components <- compute_reproduction_components(initial_state, fixed_parameters)
  r01 <- beta_h * beta_m * components[["C1"]] * components[["C2"]]
  r02 <- beta_h * beta_m * components[["C1"]] * components[["C3"]]
  r0 <- sqrt(r01 + r02)

  testthat::expect_equal(
    negative_loglikelihood_r01(
      par = c(beta_h, r01),
      observation_times = observation_times,
      observed_cases = observed_cases,
      initial_state = initial_state,
      time_grid = time_grid,
      fixed_parameters = fixed_parameters
    ),
    direct_value,
    tolerance = 1e-8
  )

  testthat::expect_equal(
    negative_loglikelihood_r02(
      par = c(beta_h, r02),
      observation_times = observation_times,
      observed_cases = observed_cases,
      initial_state = initial_state,
      time_grid = time_grid,
      fixed_parameters = fixed_parameters
    ),
    direct_value,
    tolerance = 1e-8
  )

  testthat::expect_equal(
    negative_loglikelihood_r0(
      par = c(beta_h, r0),
      observation_times = observation_times,
      observed_cases = observed_cases,
      initial_state = initial_state,
      time_grid = time_grid,
      fixed_parameters = fixed_parameters
    ),
    direct_value,
    tolerance = 1e-8
  )
})

testthat::test_that("the Hermosillo analysis runner returns a structured result", {
  testthat::skip_if_not_installed("deSolve")
  testthat::skip_if_not_installed("fields")

  result <- run_hermosillo_2010_analysis(
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

  testthat::expect_true(is.list(result))
  testthat::expect_named(result$reproduction_numbers, c("R01", "R02", "R0"))
  testthat::expect_true(all(vapply(result$output_files, length, integer(1)) == 0L))
})

testthat::test_that("the Hermosillo analysis runner writes deterministic timestamped outputs", {
  testthat::skip_if_not_installed("deSolve")
  testthat::skip_if_not_installed("fields")

  output_dir <- file.path(tempdir(), "hermosillo-analysis-outputs")
  timestamp <- "2026_03_23__14_30"

  result <- run_hermosillo_2010_analysis(
    observed_rows = 1:3,
    full_time_grid = seq(0, 5, by = 0.5),
    calibration_grid = seq(0, 2, by = 0.5),
    fit_time_grid = seq(0, 5, by = 0.5),
    beta_h_values = seq(1, 1.4, length.out = 3),
    beta_m_values = seq(0.001, 0.003, length.out = 3),
    r01_values = seq(1.7, 1.9, length.out = 3),
    r02_values = seq(0.003, 0.0034, length.out = 3),
    r0_values = seq(1.3, 1.4, length.out = 3),
    output_dir = output_dir,
    timestamp = timestamp,
    write_outputs = TRUE
  )

  all_files <- unlist(result$output_files, recursive = TRUE, use.names = FALSE)

  testthat::expect_length(all_files, 22L)
  testthat::expect_true(all(file.exists(all_files)))
  testthat::expect_true(all(grepl(
    paste0("^", timestamp, "__"),
    basename(all_files)
  )))
})
