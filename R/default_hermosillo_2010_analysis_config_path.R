#' Locate the Packaged Hermosillo 2010 Analysis Config
#'
#' @return A character scalar path to the packaged JSON configuration file.
#' @keywords internal
#' @noRd
default_hermosillo_2010_analysis_config_path <- function() {
  installed_path <- system.file(
    "extdata",
    "hermosillo_2010_analysis_config.json",
    package = "twostraindengue"
  )

  if (nzchar(installed_path)) {
    return(installed_path)
  }

  package_path("extdata", "hermosillo_2010_analysis_config.json")
}
