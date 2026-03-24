#' Hermosillo 2010 Analysis Workflow
#'
#' These functions expose plotting helpers and a reusable analysis runner for
#' the Hermosillo 2010 case study bundled with the package.
#'
#' @name hermosillo_analysis
NULL

build_output_path <- function(output_dir, subdir, stem, extension, timestamp) {
  directory <- file.path(output_dir, subdir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, sprintf("%s__%s.%s", timestamp, stem, extension))
}

save_plot <- function(path, plot_call) {
  grDevices::jpeg(path, width = 1800, height = 1200, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_call()
  invisible(path)
}

write_likelihood_outputs <- function(loglikelihood, log_file, relative_file) {
  relative <- relative_likelihood(loglikelihood)
  utils::write.csv(loglikelihood, log_file)
  utils::write.csv(relative, relative_file)
  relative
}

format_interval_labels <- function(interval) {
  if (anyNA(interval)) {
    return(NULL)
  }

  c(
    sprintf("LI=%.4f", interval[["lower"]]),
    sprintf("LS=%.4f", interval[["upper"]])
  )
}

evaluate_surface_for_mode <- function(beta_values,
                                      reproduction_values,
                                      mode,
                                      observation_times,
                                      observed_cases,
                                      initial_state,
                                      time_grid,
                                      fixed_parameters,
                                      log_file = NULL,
                                      relative_file = NULL) {
  likelihood_fn <- switch(
    mode,
    R01 = negative_loglikelihood_r01,
    R02 = negative_loglikelihood_r02,
    R0 = negative_loglikelihood_r0,
    stop("Unsupported reproduction mode: ", mode, call. = FALSE)
  )

  loglikelihood <- evaluate_loglikelihood_grid(
    beta_values,
    reproduction_values,
    function(par) {
      likelihood_fn(
        par = par,
        observation_times = observation_times,
        observed_cases = observed_cases,
        initial_state = initial_state,
        time_grid = time_grid,
        fixed_parameters = fixed_parameters
      )
    }
  )

  relative <- if (is.null(log_file) || is.null(relative_file)) {
    relative_likelihood(loglikelihood)
  } else {
    write_likelihood_outputs(loglikelihood, log_file, relative_file)
  }

  list(
    loglikelihood = loglikelihood,
    relative = relative,
    files = list(loglikelihood = log_file, relative = relative_file)
  )
}

#' @rdname hermosillo_analysis
#'
#' @param data A data frame with columns `week`, `classical`, and
#'   `hemorrhagic`.
#' @param rows Integer row indices to display.
#' @param title Plot title.
#' @param ylim Optional y-axis limits.
#'
#' @return Invisibly returns the plotted subset of `data`.
#' @export
plot_observed_cases <- function(data,
                                rows = seq_len(nrow(data)),
                                title = "Observed Cases",
                                ylim = NULL) {
  selected <- data[rows, , drop = FALSE]

  if (is.null(ylim)) {
    ylim <- c(0, max(selected$classical, selected$hemorrhagic))
  }

  graphics::matplot(
    selected$week,
    cbind(selected$classical, selected$hemorrhagic),
    type = "p",
    pch = c(19, 19),
    col = c(1, 2),
    cex.lab = 1.35,
    cex.axis = 1.35,
    main = title,
    xlab = "Week",
    ylab = "Infecteds",
    cex.main = 1.30,
    ylim = ylim
  )
  graphics::legend(
    "topleft",
    c("Classical", "Hemorrhagic"),
    cex = 1.2,
    pch = c(19, 19),
    col = c(1, 2),
    bty = "n",
    y.intersp = 0.35,
    inset = -0.025
  )

  invisible(selected)
}

#' @rdname hermosillo_analysis
#'
#' @param model_output Matrix or data frame returned by [run_model()].
#' @param week_origin Week number corresponding to simulation time zero.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return Invisibly returns the fitted case matrix used in the overlay.
#' @export
plot_model_fit <- function(model_output,
                           data,
                           rows = seq_len(nrow(data)),
                           title = "Model Fit",
                           week_origin = min(data$week[rows]),
                           xlim = NULL,
                           ylim = NULL) {
  fitted_cases <- extract_reported_cases(model_output)
  selected <- data[rows, , drop = FALSE]

  if (is.null(ylim)) {
    ylim <- c(0, max(fitted_cases, selected$classical, selected$hemorrhagic))
  }

  graphics::matplot(
    model_output[, "time"] + week_origin,
    fitted_cases,
    type = "l",
    lty = 1,
    lwd = 2,
    col = c(1, 2),
    cex.lab = 1.25,
    cex.axis = 1.25,
    main = title,
    xlab = "Week",
    ylab = "Infecteds",
    cex.main = 1.30,
    xlim = xlim,
    ylim = ylim
  )
  graphics::legend(
    "topleft",
    c("Classical", "Hemorrhagic"),
    lty = c(1, 1),
    col = c(1, 2),
    bty = "n",
    y.intersp = 0.35,
    inset = -0.025
  )
  graphics::points(selected$week, selected$classical, pch = 19, col = 1, cex = 0.8)
  graphics::points(selected$week, selected$hemorrhagic, pch = 19, col = 2, cex = 0.8)

  invisible(fitted_cases)
}

#' @rdname hermosillo_analysis
#'
#' @param x_values Numeric vector for the x-axis.
#' @param y_values Numeric vector for the y-axis.
#' @param surface Numeric matrix of likelihood values.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point Optional numeric vector of length 2 to highlight on the plot.
#' @param label Optional legend label for `point`.
#' @param main Plot title.
#'
#' @return Invisibly returns `surface`.
#' @export
plot_likelihood_surface <- function(x_values,
                                    y_values,
                                    surface,
                                    x_label,
                                    y_label,
                                    point = NULL,
                                    label = NULL,
                                    main = "Likelihood Contours") {
  fields::image.plot(
    x_values,
    y_values,
    surface,
    xlab = x_label,
    ylab = y_label,
    cex.lab = 1.25,
    main = main
  )
  graphics::contour(
    x_values,
    y_values,
    surface,
    levels = c(0.05),
    add = TRUE,
    drawlabels = TRUE,
    labcex = 1
  )

  if (!is.null(point)) {
    graphics::points(point[1], point[2], pch = 19, col = 2, cex = 1.4)
  }

  if (!is.null(label)) {
    graphics::legend(
      "topright",
      legend = label,
      pch = 19,
      col = 2,
      cex = 1.2,
      bty = "n",
      inset = -0.04
    )
  }

  invisible(surface)
}

#' @rdname hermosillo_analysis
#'
#' @param values Numeric vector defining the profile grid.
#' @param relative_profile Numeric vector of relative profile likelihood values.
#' @param x_label X-axis label.
#' @param interval Named numeric vector with `lower` and `upper`.
#' @param interval_label Optional character vector for the interval legend.
#' @param main Plot title.
#'
#' @return Invisibly returns `relative_profile`.
#' @export
plot_profile_likelihood <- function(values,
                                    relative_profile,
                                    x_label,
                                    interval = c(lower = NA_real_, upper = NA_real_),
                                    interval_label = NULL,
                                    main = "Relative Profile Likelihood") {
  if (is.null(interval_label)) {
    interval_label <- format_interval_labels(interval)
  }

  graphics::plot(
    values,
    relative_profile,
    type = "l",
    lty = 1,
    lwd = 2,
    col = 1,
    cex.lab = 1.25,
    cex.axis = 1.25,
    main = main,
    xlab = x_label,
    ylim = c(0, 1),
    ylab = "Relative Profile Likelihood",
    cex.main = 1.30
  )

  if (!anyNA(interval)) {
    graphics::points(interval[["lower"]], 0.015, pch = -9658, col = 1, cex = 1)
    graphics::points(interval[["upper"]], 0.015, pch = -9668, col = 1, cex = 1)
  }

  if (!is.null(interval_label)) {
    graphics::legend(
      "topright",
      interval_label,
      pch = c(-9658, -9668)[seq_along(interval_label)],
      col = rep(1, length(interval_label)),
      bty = "n",
      title.adj = 0.35,
      y.intersp = 0.5,
      inset = -0.025,
      title = "95%CI"
    )
  }

  invisible(relative_profile)
}

#' @rdname hermosillo_analysis
#'
#' @param data Study data used for calibration and plotting.
#' @param observed_rows Integer row indices defining the calibration window.
#' @param full_time_grid Numeric solver grid used for the baseline simulation.
#' @param calibration_grid Numeric solver grid used during likelihood
#'   evaluation.
#' @param fit_time_grid Optional solver grid used for the fitted trajectory. If
#'   `NULL`, a grid covering the observed study window is built automatically.
#' @param beta_h_values Grid of `betaH` values for the original-scale surface.
#' @param beta_m_values Grid of `betaM` values for the original-scale surface.
#' @param r01_values Grid of `R01` values for the reparameterised surface.
#' @param r02_values Grid of `R02` values for the reparameterised surface.
#' @param r0_values Grid of `R0` values for the reparameterised surface.
#' @param output_dir Directory where figures and calibration tables are written.
#' @param timestamp Timestamp prefix used in output file names.
#' @param write_outputs Logical; if `TRUE`, writes figures and CSV surfaces to
#'   `output_dir`.
#' @param optim_start Starting values for the constrained optimisation of
#'   `betaH` and `betaM`.
#'
#' @return A named list containing the fitted parameters, likelihood surfaces,
#'   profile summaries, reproduction-number estimates, and any generated output
#'   file paths.
#'
#' @examples
#' analysis <- run_hermosillo_2010_analysis(
#'   observed_rows = 1:3,
#'   full_time_grid = seq(0, 5, by = 0.5),
#'   calibration_grid = seq(0, 2, by = 0.5),
#'   fit_time_grid = seq(0, 5, by = 0.5),
#'   beta_h_values = seq(1, 1.4, length.out = 3),
#'   beta_m_values = seq(0.001, 0.003, length.out = 3),
#'   r01_values = seq(1.7, 1.9, length.out = 3),
#'   r02_values = seq(0.003, 0.0034, length.out = 3),
#'   r0_values = seq(1.3, 1.4, length.out = 3),
#'   write_outputs = FALSE
#' )
#'
#' analysis$reproduction_numbers
#' @export
run_hermosillo_2010_analysis <- function(
    data = NULL,
    observed_rows = 1:6,
    full_time_grid = seq(0, 18, by = 0.1),
    calibration_grid = seq(0, 6, by = 0.1),
    fit_time_grid = NULL,
    beta_h_values = seq(1, 3.5, length.out = 250),
    beta_m_values = seq(0.001, 0.0175, length.out = 250),
    r01_values = seq(1.65, 3, length.out = 250),
    r02_values = seq(0.003, 0.006, length.out = 250),
    r0_values = seq(1.25, 1.75, length.out = 250),
    output_dir = "outputs",
    timestamp = NULL,
    write_outputs = TRUE,
    optim_start = c(0.9, 0.004)) {
  if (is.null(data)) {
    data <- load_study_data()
  }

  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y_%m_%d__%H_%M")
  }

  all_parameters <- default_parameters()
  fixed_parameters <- all_parameters[fixed_parameter_names]
  initial_state <- default_state(all_parameters)
  full_rows <- seq_len(nrow(data))
  week_origin <- min(data$week[observed_rows])

  if (is.null(fit_time_grid)) {
    fit_time_grid <- seq(0, max(data$week) - week_origin, by = 0.1)
  }

  observation_times <- data$week[observed_rows] - week_origin
  observed_cases <- as.matrix(data[observed_rows, c("classical", "hemorrhagic")])

  baseline_output <- run_model(initial_state, full_time_grid, all_parameters)

  global_loglikelihood <- evaluate_loglikelihood_grid(
    beta_h_values,
    beta_m_values,
    function(par) {
      negative_loglikelihood_original_scale(
        par = par,
        observation_times = observation_times,
        observed_cases = observed_cases,
        initial_state = initial_state,
        time_grid = calibration_grid,
        fixed_parameters = fixed_parameters
      )
    }
  )

  output_files <- list(figures = list(), calibration = list())

  if (write_outputs) {
    beta_h_beta_m_log <- build_output_path(
      output_dir,
      "calibration",
      "loglik_surface_beta_h_beta_m_poisson",
      "csv",
      timestamp
    )
    beta_h_beta_m_relative <- build_output_path(
      output_dir,
      "calibration",
      "relative_loglik_surface_beta_h_beta_m_poisson",
      "csv",
      timestamp
    )
    relative_global_loglikelihood <- write_likelihood_outputs(
      global_loglikelihood,
      beta_h_beta_m_log,
      beta_h_beta_m_relative
    )
    output_files$calibration$beta_h_beta_m <- list(
      loglikelihood = beta_h_beta_m_log,
      relative = beta_h_beta_m_relative
    )
  } else {
    relative_global_loglikelihood <- relative_likelihood(global_loglikelihood)
  }

  mle_fit <- stats::constrOptim(
    theta = optim_start,
    f = negative_loglikelihood_original_scale,
    grad = NULL,
    ui = diag(2),
    ci = c(0, 0),
    mu = 1e-04,
    method = "Nelder-Mead",
    outer.iterations = 100,
    outer.eps = 1e-05,
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = calibration_grid,
    fixed_parameters = fixed_parameters
  )
  mle_parameters <- stats::setNames(mle_fit$par, c("betaH", "betaM"))

  fit_output <- run_model(
    initial_state,
    fit_time_grid,
    build_parameter_vector(fixed_parameters, mle_parameters[["betaH"]], mle_parameters[["betaM"]])
  )

  beta_h_profile <- relative_likelihood(apply(global_loglikelihood, 1, max))
  beta_h_interval <- profile_interval(beta_h_values, beta_h_profile)
  beta_m_profile <- relative_likelihood(apply(global_loglikelihood, 2, max))
  beta_m_interval <- profile_interval(beta_m_values, beta_m_profile)

  reproduction_components <- compute_reproduction_components(initial_state, fixed_parameters)
  pi_r <- mle_parameters[["betaH"]] * mle_parameters[["betaM"]] * reproduction_components[["C1"]]
  reproduction_numbers <- c(
    R01 = pi_r * reproduction_components[["C2"]],
    R02 = pi_r * reproduction_components[["C3"]],
    R0 = sqrt((pi_r * reproduction_components[["C2"]]) + (pi_r * reproduction_components[["C3"]]))
  )

  r01_surface <- evaluate_surface_for_mode(
    beta_values = beta_h_values,
    reproduction_values = r01_values,
    mode = "R01",
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = calibration_grid,
    fixed_parameters = fixed_parameters,
    log_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "loglik_surface_beta_h_r1_poisson", "csv", timestamp)
    } else {
      NULL
    },
    relative_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "relative_loglik_surface_beta_h_r1_poisson", "csv", timestamp)
    } else {
      NULL
    }
  )
  r02_surface <- evaluate_surface_for_mode(
    beta_values = beta_h_values,
    reproduction_values = r02_values,
    mode = "R02",
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = calibration_grid,
    fixed_parameters = fixed_parameters,
    log_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "loglik_surface_beta_h_r2_poisson", "csv", timestamp)
    } else {
      NULL
    },
    relative_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "relative_loglik_surface_beta_h_r2_poisson", "csv", timestamp)
    } else {
      NULL
    }
  )
  r0_surface <- evaluate_surface_for_mode(
    beta_values = beta_h_values,
    reproduction_values = r0_values,
    mode = "R0",
    observation_times = observation_times,
    observed_cases = observed_cases,
    initial_state = initial_state,
    time_grid = calibration_grid,
    fixed_parameters = fixed_parameters,
    log_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "loglik_surface_beta_h_r0_poisson", "csv", timestamp)
    } else {
      NULL
    },
    relative_file = if (write_outputs) {
      build_output_path(output_dir, "calibration", "relative_loglik_surface_beta_h_r0_poisson", "csv", timestamp)
    } else {
      NULL
    }
  )

  if (write_outputs) {
    output_files$calibration$r01 <- r01_surface$files
    output_files$calibration$r02 <- r02_surface$files
    output_files$calibration$r0 <- r0_surface$files
  }

  r01_profile <- relative_likelihood(apply(r01_surface$loglikelihood, 2, max))
  r01_interval <- profile_interval(r01_values, r01_profile)
  r02_profile <- relative_likelihood(apply(r02_surface$loglikelihood, 2, max))
  r02_interval <- profile_interval(r02_values, r02_profile)
  r0_profile <- relative_likelihood(apply(r0_surface$loglikelihood, 2, max))
  r0_interval <- profile_interval(r0_values, r0_profile)

  if (write_outputs) {
    output_files$figures$observed_cases_full <- save_plot(
      build_output_path(output_dir, "figures", "observed_cases_full", "jpeg", timestamp),
      function() {
        plot_observed_cases(
          data = data,
          rows = full_rows,
          title = "Observed Cases: Hermosillo 2010",
          ylim = c(0, 200)
        )
      }
    )
    output_files$figures$model_fit_full_window <- save_plot(
      build_output_path(output_dir, "figures", "model_fit_full_window", "jpeg", timestamp),
      function() {
        plot_model_fit(
          model_output = baseline_output,
          data = data,
          rows = full_rows,
          title = "Baseline Model Fit: Hermosillo 2010",
          week_origin = week_origin
        )
      }
    )
    output_files$figures$observed_cases_calibration_window <- save_plot(
      build_output_path(output_dir, "figures", "observed_cases_calibration_window", "jpeg", timestamp),
      function() {
        plot_observed_cases(
          data = data,
          rows = observed_rows,
          title = "Observed Cases: Calibration Window",
          ylim = c(0, 150)
        )
      }
    )
    output_files$figures$likelihood_surface_beta_h_beta_m <- save_plot(
      build_output_path(output_dir, "figures", "likelihood_surface_beta_h_beta_m", "jpeg", timestamp),
      function() {
        plot_likelihood_surface(
          beta_h_values,
          beta_m_values,
          relative_global_loglikelihood,
          x_label = expression(beta[H]),
          y_label = expression(beta[M]),
          point = unname(mle_parameters),
          label = sprintf(
            "(betaH, betaM) = (%.4f, %.4f)",
            mle_parameters[["betaH"]],
            mle_parameters[["betaM"]]
          )
        )
      }
    )
    output_files$figures$model_fit_calibration_window <- save_plot(
      build_output_path(output_dir, "figures", "model_fit_calibration_window", "jpeg", timestamp),
      function() {
        plot_model_fit(
          model_output = fit_output,
          data = data,
          rows = observed_rows,
          title = "Poisson-Strain Model Fit: Calibration Window",
          week_origin = week_origin,
          xlim = range(data$week[observed_rows])
        )
      }
    )
    output_files$figures$model_fit_full_data <- save_plot(
      build_output_path(output_dir, "figures", "model_fit_full_data", "jpeg", timestamp),
      function() {
        plot_model_fit(
          model_output = fit_output,
          data = data,
          rows = full_rows,
          title = "Poisson-Strain Model Fit: Full Study Window",
          week_origin = week_origin,
          xlim = range(data$week[full_rows])
        )
      }
    )
    output_files$figures$profile_likelihood_beta_h <- save_plot(
      build_output_path(output_dir, "figures", "profile_likelihood_beta_h", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          beta_h_values,
          beta_h_profile,
          x_label = expression(beta[H]),
          interval = beta_h_interval,
          main = expression(paste("Profile Likelihood for ", beta[H]))
        )
      }
    )
    output_files$figures$profile_likelihood_beta_m <- save_plot(
      build_output_path(output_dir, "figures", "profile_likelihood_beta_m", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          beta_m_values,
          beta_m_profile,
          x_label = expression(beta[M]),
          interval = beta_m_interval,
          main = expression(paste("Profile Likelihood for ", beta[M]))
        )
      }
    )
    output_files$figures$likelihood_surface_beta_h_r1 <- save_plot(
      build_output_path(output_dir, "figures", "likelihood_surface_beta_h_r1", "jpeg", timestamp),
      function() {
        plot_likelihood_surface(
          beta_h_values,
          r01_values,
          r01_surface$relative,
          x_label = expression(beta[H]),
          y_label = expression(R[1]),
          point = c(mle_parameters[["betaH"]], reproduction_numbers[["R01"]]),
          label = sprintf(
            "(betaH, R1) = (%.4f, %.4f)",
            mle_parameters[["betaH"]],
            reproduction_numbers[["R01"]]
          )
        )
      }
    )
    output_files$figures$likelihood_surface_beta_h_r2 <- save_plot(
      build_output_path(output_dir, "figures", "likelihood_surface_beta_h_r2", "jpeg", timestamp),
      function() {
        plot_likelihood_surface(
          beta_h_values,
          r02_values,
          r02_surface$relative,
          x_label = expression(beta[H]),
          y_label = expression(R[2]),
          point = c(mle_parameters[["betaH"]], reproduction_numbers[["R02"]]),
          label = sprintf(
            "(betaH, R2) = (%.4f, %.4f)",
            mle_parameters[["betaH"]],
            reproduction_numbers[["R02"]]
          )
        )
      }
    )
    output_files$figures$likelihood_surface_beta_h_r0 <- save_plot(
      build_output_path(output_dir, "figures", "likelihood_surface_beta_h_r0", "jpeg", timestamp),
      function() {
        plot_likelihood_surface(
          beta_h_values,
          r0_values,
          r0_surface$relative,
          x_label = expression(beta[H]),
          y_label = expression(R[0]),
          point = c(mle_parameters[["betaH"]], reproduction_numbers[["R0"]]),
          label = sprintf(
            "(betaH, R0) = (%.4f, %.4f)",
            mle_parameters[["betaH"]],
            reproduction_numbers[["R0"]]
          )
        )
      }
    )
    output_files$figures$profile_likelihood_r1 <- save_plot(
      build_output_path(output_dir, "figures", "profile_likelihood_r1", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          r01_values,
          r01_profile,
          x_label = expression(R[1]),
          interval = r01_interval,
          main = expression(paste("Profile Likelihood for ", R[1]))
        )
      }
    )
    output_files$figures$profile_likelihood_r2 <- save_plot(
      build_output_path(output_dir, "figures", "profile_likelihood_r2", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          r02_values,
          r02_profile,
          x_label = expression(R[2]),
          interval = r02_interval,
          main = expression(paste("Profile Likelihood for ", R[2]))
        )
      }
    )
    output_files$figures$profile_likelihood_r0 <- save_plot(
      build_output_path(output_dir, "figures", "profile_likelihood_r0", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          r0_values,
          r0_profile,
          x_label = expression(R[0]),
          interval = r0_interval,
          main = expression(paste("Profile Likelihood for ", R[0]))
        )
      }
    )
  }

  list(
    timestamp = timestamp,
    week_origin = week_origin,
    data = data,
    observed_rows = observed_rows,
    initial_state = initial_state,
    fixed_parameters = fixed_parameters,
    baseline_output = baseline_output,
    mle = list(
      optimisation = mle_fit,
      parameters = mle_parameters,
      fit_output = fit_output
    ),
    reproduction_numbers = reproduction_numbers,
    surfaces = list(
      beta_h_beta_m = list(
        loglikelihood = global_loglikelihood,
        relative = relative_global_loglikelihood
      ),
      r01 = r01_surface,
      r02 = r02_surface,
      r0 = r0_surface
    ),
    profiles = list(
      beta_h = list(values = beta_h_values, relative = beta_h_profile, interval = beta_h_interval),
      beta_m = list(values = beta_m_values, relative = beta_m_profile, interval = beta_m_interval),
      r01 = list(values = r01_values, relative = r01_profile, interval = r01_interval),
      r02 = list(values = r02_values, relative = r02_profile, interval = r02_interval),
      r0 = list(values = r0_values, relative = r0_profile, interval = r0_interval)
    ),
    output_files = output_files
  )
}
