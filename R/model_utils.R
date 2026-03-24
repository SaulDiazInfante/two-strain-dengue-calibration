fixed_parameter_names <- c(
  "LambdaM",
  "LambdaS",
  "Lambda1",
  "alphaC",
  "alphaH",
  "b",
  "muH",
  "muM",
  "vsigma",
  "vtheta",
  "p"
)

full_parameter_names <- c(
  fixed_parameter_names[1:6],
  "betaH",
  "betaM",
  fixed_parameter_names[7:11]
)

package_path <- function(..., must_work = TRUE) {
  installed_path <- system.file(..., package = "twostraindengue", mustWork = FALSE)

  if (nzchar(installed_path)) {
    return(installed_path)
  }

  source_path <- file.path("inst", ...)

  if (!must_work || file.exists(source_path)) {
    return(source_path)
  }

  stop("Cannot locate packaged file: ", file.path(...), call. = FALSE)
}

#' Model Defaults and Simulation Helpers
#'
#' These functions expose the default parameterisation, initial conditions, raw
#' Hermosillo study data loader, and the ODE simulation workflow used
#' throughout the package.
#'
#' @name model_setup
NULL

#' @rdname model_setup
#'
#' @return A named numeric vector containing the default model parameters.
#'
#' @examples
#' parameters <- default_parameters()
#' state <- default_state(parameters)
#' time_grid <- seq(0, 1, by = 0.1)
#'
#' model_output <- run_model(state, time_grid, parameters)
#' head(extract_reported_cases(model_output))
#' @export
default_parameters <- function() {
  c(
    LambdaM = 30702.6139006,
    LambdaS = 76.89246,
    Lambda1 = 2834.92,
    alphaC = 1.1655,
    alphaH = 1.1655,
    b = 21.875,
    betaH = 0.95,
    betaM = 0.01503,
    muH = 0.000273973,
    muM = 0.5075,
    vsigma = 2.5,
    vtheta = 0.05,
    p = 0.025
  )
}

#' @rdname model_setup
#'
#' @param parameters A named numeric vector, typically returned by
#'   [default_parameters()].
#'
#' @return A named numeric vector containing the default initial state for the
#'   ODE system, including the cumulative reported-case state `z`.
#' @export
default_state <- function(parameters = default_parameters()) {
  state <- c(
    Ms = 120000,
    M1 = 10,
    M2 = 10,
    S = 278931,
    I1 = 120,
    I2 = 40,
    S1 = 4400,
    Y1c = 0,
    Y1h = 1,
    R = 0
  )
  c(state, z = parameters[["p"]] * sum(state[c("I1", "I2", "Y1c")]))
}

#' @rdname model_setup
#'
#' @param path Path to the raw study CSV. By default, this resolves to the file
#'   shipped under `inst/extdata/raw/`.
#' @param rows Integer row indices to select from the raw study data.
#'
#' @return A data frame with columns `week`, `classical`, and `hemorrhagic`.
#' @export
load_study_data <- function(
                            path = system.file(
                              "extdata",
                              "raw",
                              "study_data_dengue_hemorrhagic_fever.csv",
                              package = "twostraindengue"
                            ),
                            rows = 35:50) {
  if (!nzchar(path)) {
    path <- package_path(
      "extdata",
      "raw",
      "study_data_dengue_hemorrhagic_fever.csv"
    )
  }

  raw_data <- utils::read.csv(path)
  selected <- raw_data[rows, 2:4]
  names(selected) <- c("week", "classical", "hemorrhagic")
  selected
}

normalise_fixed_parameters <- function(parameters) {
  if (is.null(names(parameters))) {
    names(parameters) <- fixed_parameter_names[seq_along(parameters)]
  }

  missing_names <- setdiff(fixed_parameter_names, names(parameters))
  if (length(missing_names) > 0) {
    stop("Missing fixed parameters: ", paste(missing_names, collapse = ", "))
  }

  parameters[fixed_parameter_names]
}

build_parameter_vector <- function(fixed_parameters, beta_h, beta_m) {
  fixed_parameters <- normalise_fixed_parameters(fixed_parameters)
  c(
    fixed_parameters[1:6],
    betaH = beta_h,
    betaM = beta_m,
    fixed_parameters[7:11]
  )[full_parameter_names]
}

#' @rdname model_setup
#'
#' @param initial_state Named numeric state vector passed to [deSolve::ode()].
#' @param time_grid Numeric vector of solver times.
#' @param parameters Named numeric parameter vector.
#' @param ode_function ODE callback used by [deSolve::ode()]. Defaults to the
#'   package's internal two-strain dengue system.
#'
#' @return A matrix of ODE output as returned by [deSolve::ode()].
#' @export
run_model <- function(initial_state,
                      time_grid,
                      parameters,
                      ode_function = dengue_model_ode) {
  deSolve::ode(
    y = initial_state,
    times = time_grid,
    func = ode_function,
    parms = parameters
  )
}

#' @rdname model_setup
#'
#' @param model_output Matrix or data frame returned by [run_model()].
#'
#' @return A two-column matrix with `classical` and `hemorrhagic` reported-case
#'   trajectories extracted from the ODE output.
#' @export
extract_reported_cases <- function(model_output) {
  output <- as.data.frame(model_output)
  cbind(classical = output[["z"]], hemorrhagic = output[["Y1h"]])
}

match_observation_times <- function(model_output, observation_times) {
  model_times <- round(model_output[, "time"], digits = 8)
  matched_rows <- match(round(observation_times, digits = 8), model_times)

  if (anyNA(matched_rows)) {
    stop("Observation times are not available in the ODE output grid.")
  }

  matched_rows
}

poisson_negative_loglikelihood <- function(expected_cases, observed_cases) {
  safe_expected <- pmax(expected_cases, .Machine$double.eps)
  -(sum(observed_cases * log(safe_expected)) - sum(safe_expected))
}

#' Calibration Utilities
#'
#' These helpers support likelihood evaluation, reparameterisation by
#' reproduction numbers, likelihood surfaces, and profile-likelihood summaries.
#'
#' @name calibration_tools
NULL

#' @rdname calibration_tools
#'
#' @param beta_h Human-to-mosquito transmission coefficient.
#' @param beta_m Mosquito-to-human transmission coefficient.
#' @param observation_times Observation times on the simulation grid.
#' @param observed_cases Two-column matrix of observed classical and
#'   hemorrhagic case counts.
#' @param initial_state Named numeric state vector used by [run_model()].
#' @param time_grid Numeric vector of solver times.
#' @param fixed_parameters Named numeric vector containing the fixed model
#'   parameters.
#'
#' @return A scalar negative Poisson log-likelihood.
#' @export
evaluate_negative_loglikelihood <- function(beta_h,
                                            beta_m,
                                            observation_times,
                                            observed_cases,
                                            initial_state,
                                            time_grid,
                                            fixed_parameters) {
  parameters <- build_parameter_vector(fixed_parameters, beta_h, beta_m)
  model_output <- run_model(initial_state, time_grid, parameters)
  matched_rows <- match_observation_times(model_output, observation_times)
  fitted_cases <- extract_reported_cases(model_output)[matched_rows, , drop = FALSE]

  poisson_negative_loglikelihood(
    fitted_cases[, "classical"],
    observed_cases[, 1]
  ) +
    poisson_negative_loglikelihood(
      fitted_cases[, "hemorrhagic"],
      observed_cases[, 2]
    )
}

#' @rdname calibration_tools
#'
#' @param initial_state Named numeric state vector used by [run_model()].
#' @param fixed_parameters Named numeric vector containing the fixed model
#'   parameters.
#'
#' @return A named numeric vector with the reproduction-number components
#'   `C1`, `C2`, and `C3`.
#' @export
compute_reproduction_components <- function(initial_state, fixed_parameters) {
  fixed_parameters <- normalise_fixed_parameters(fixed_parameters)
  human_population <- sum(initial_state[4:10])
  secondary_population <- sum(initial_state[7:10])

  c(
    C1 = fixed_parameters[["LambdaM"]] *
      (fixed_parameters[["b"]] / (fixed_parameters[["muM"]] * human_population)) ^ 2,
    C2 = (human_population - secondary_population +
      (1 - fixed_parameters[["vtheta"]]) *
      fixed_parameters[["vsigma"]] *
      secondary_population) / (fixed_parameters[["alphaC"]] + fixed_parameters[["muH"]]),
    C3 = (fixed_parameters[["vsigma"]] *
      fixed_parameters[["vtheta"]] *
      secondary_population) / (fixed_parameters[["alphaH"]] + fixed_parameters[["muH"]])
  )
}

#' @rdname calibration_tools
#'
#' @param reproduction_value Value of `R01`, `R02`, or `R0`, depending on
#'   `mode`.
#' @param mode Character scalar selecting `"R01"`, `"R02"`, or `"R0"`.
#'
#' @return The implied mosquito-to-human transmission coefficient `betaM`.
#' @export
beta_m_from_reproduction <- function(beta_h,
                                     reproduction_value,
                                     mode,
                                     initial_state,
                                     fixed_parameters) {
  components <- compute_reproduction_components(initial_state, fixed_parameters)

  if (mode == "R01") {
    return(reproduction_value / (beta_h * components[["C1"]] * components[["C2"]]))
  }

  if (mode == "R02") {
    return(reproduction_value / (beta_h * components[["C1"]] * components[["C3"]]))
  }

  if (mode == "R0") {
    return((reproduction_value ^ 2) /
      (beta_h * components[["C1"]] * (components[["C2"]] + components[["C3"]])))
  }

  stop("Unsupported reproduction mode: ", mode)
}

#' @rdname calibration_tools
#'
#' @param reproduction_value Value of `R01`, `R02`, or `R0`, depending on
#'   `mode`.
#' @param mode Character scalar selecting `"R01"`, `"R02"`, or `"R0"`.
#'
#' @return A scalar negative Poisson log-likelihood after converting the
#'   reproduction-number parameterisation back to `betaM`.
#' @export
evaluate_reparameterised_loglikelihood <- function(beta_h,
                                                   reproduction_value,
                                                   mode,
                                                   observation_times,
                                                   observed_cases,
                                                   initial_state,
                                                   time_grid,
                                                   fixed_parameters) {
  beta_m <- beta_m_from_reproduction(
    beta_h = beta_h,
    reproduction_value = reproduction_value,
    mode = mode,
    initial_state = initial_state,
    fixed_parameters = fixed_parameters
  )

  evaluate_negative_loglikelihood(
    beta_h = beta_h,
    beta_m = beta_m,
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = time_grid,
    fixed_parameters = fixed_parameters
  )
}

#' @rdname calibration_tools
#'
#' @param x_values Numeric vector for the first parameter axis.
#' @param y_values Numeric vector for the second parameter axis.
#' @param evaluator Function taking a two-element numeric vector and returning a
#'   negative log-likelihood.
#'
#' @return A matrix of log-likelihood values on the supplied grid.
#' @export
evaluate_loglikelihood_grid <- function(x_values, y_values, evaluator) {
  if (length(x_values) == 0L || length(y_values) == 0L) {
    return(matrix(numeric(0), nrow = length(x_values), ncol = length(y_values)))
  }

  t(vapply(
    x_values,
    function(x_value) {
      vapply(
        y_values,
        function(y_value) -evaluator(c(x_value, y_value)),
        numeric(1)
      )
    },
    numeric(length(y_values))
  ))
}

#' @rdname calibration_tools
#'
#' @param loglikelihood Numeric vector or matrix of log-likelihood values.
#'
#' @return Relative likelihood values normalised by the maximum entry.
#' @export
relative_likelihood <- function(loglikelihood) {
  exp(loglikelihood - max(loglikelihood))
}

#' @rdname calibration_tools
#'
#' @param values Numeric vector defining the profile grid.
#' @param relative_profile Numeric vector of relative profile likelihood values.
#' @param cutoff Relative-likelihood cutoff used to define the interval.
#'
#' @return A named numeric vector with entries `lower` and `upper`.
#' @export
profile_interval <- function(values, relative_profile, cutoff = 0.146) {
  inside <- which(relative_profile >= cutoff)

  if (length(inside) == 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  c(lower = values[min(inside)], upper = values[max(inside)])
}
