#' Load Hermosillo 2010 Analysis Configuration
#'
#' @description
#' Loads the packaged JSON configuration for the Hermosillo 2010 analysis
#' workflow and optionally merges overrides from another JSON file.
#'
#' @param path Optional path to a JSON file containing configuration overrides.
#'
#' @return A normalized configuration list ready to pass to
#'   [run_hermosillo_2010_analysis()].
#'
#' @examples
#' config <- load_hermosillo_2010_analysis_config()
#' config$outputs$write_outputs
#'
#' @family hermosillo_analysis
#' @export
load_hermosillo_2010_analysis_config <- function(path = NULL) {
  config <- read_hermosillo_2010_analysis_config(
    default_hermosillo_2010_analysis_config_path()
  )

  if (!is.null(path)) {
    config <- merge_nested_lists(
      config,
      read_hermosillo_2010_analysis_config(path)
    )
  }

  normalise_hermosillo_2010_analysis_config(config)
}
