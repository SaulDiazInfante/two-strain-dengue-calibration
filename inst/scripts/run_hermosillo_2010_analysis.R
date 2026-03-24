load_analysis_function <- function() {
  if (dir.exists("R")) {
    invisible(lapply(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)), source))
    return(run_hermosillo_2010_analysis)
  }

  if (requireNamespace("twostraindengue", quietly = TRUE)) {
    return(twostraindengue::run_hermosillo_2010_analysis)
  }

  stop(
    "Cannot locate the package sources or an installed twostraindengue package.",
    call. = FALSE
  )
}

analysis_function <- load_analysis_function()
analysis_args <- list()
script_args <- commandArgs(trailingOnly = TRUE)
formal_names <- names(formals(analysis_function))

if ("config" %in% formal_names) {
  analysis_args$config <- list(outputs = list(show_progress = TRUE))

  if (length(script_args) > 0L) {
    analysis_args$config_path <- script_args[[1]]
  }
} else if ("show_progress" %in% formal_names) {
  analysis_args$show_progress <- TRUE
}

analysis_results <- do.call(analysis_function, analysis_args)

print(analysis_results$reproduction_numbers)
invisible(analysis_results)
