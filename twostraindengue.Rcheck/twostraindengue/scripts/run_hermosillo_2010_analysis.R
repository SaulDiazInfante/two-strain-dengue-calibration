load_analysis_function <- function() {
  if (requireNamespace("twostraindengue", quietly = TRUE)) {
    return(twostraindengue::run_hermosillo_2010_analysis)
  }

  invisible(lapply(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)), source))
  run_hermosillo_2010_analysis
}

analysis_function <- load_analysis_function()
analysis_results <- analysis_function()

print(analysis_results$reproduction_numbers)
invisible(analysis_results)
