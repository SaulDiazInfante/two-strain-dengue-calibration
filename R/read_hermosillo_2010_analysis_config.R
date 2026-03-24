#' Read a Hermosillo 2010 Analysis Config JSON File
#'
#' @param path Path to a JSON configuration file.
#'
#' @return A raw nested list parsed from JSON.
#' @keywords internal
#' @noRd
read_hermosillo_2010_analysis_config <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot locate analysis config file: ", path, call. = FALSE)
  }

  jsonlite::fromJSON(path, simplifyVector = FALSE)
}
