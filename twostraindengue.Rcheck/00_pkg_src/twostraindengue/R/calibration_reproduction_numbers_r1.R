#' @rdname calibration_objectives
#'
#' @param par Numeric vector with `betaH` in the first position and the
#'   wrapper-specific second parameter in the second position: `betaM` for
#'   [negative_loglikelihood_original_scale()], `R01` for
#'   [negative_loglikelihood_r01()], `R02` for [negative_loglikelihood_r02()],
#'   or `R0` for [negative_loglikelihood_r0()].
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
negative_loglikelihood_r01 <- function(par,
                                       observation_times,
                                       observed_cases,
                                       initial_state,
                                       time_grid,
                                       fixed_parameters) {
  evaluate_reparameterised_loglikelihood(
    beta_h = par[1],
    reproduction_value = par[2],
    mode = "R01",
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = time_grid,
    fixed_parameters = fixed_parameters
  )
}

FuncionMenosLogVerobetaHR01 <- function(VecPar,
                                        VecTiemposObs,
                                        VecCasosObs,
                                        Yinicial,
                                        RejillaTiempo,
                                        VecParFijo) {
  negative_loglikelihood_r01(
    par = VecPar,
    observation_times = VecTiemposObs,
    observed_cases = VecCasosObs,
    initial_state = Yinicial,
    time_grid = RejillaTiempo,
    fixed_parameters = VecParFijo
  )
}
