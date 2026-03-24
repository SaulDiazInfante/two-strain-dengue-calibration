#' Hermosillo 2010 Dengue Observations
#'
#' Weekly counts of classical and hemorrhagic dengue cases used in the package
#' examples, tests, and case-study workflow.
#'
#' @format A data frame with 16 rows and 3 variables:
#' \describe{
#'   \item{week}{Epidemiological week number.}
#'   \item{classical}{Observed classical dengue cases.}
#'   \item{hemorrhagic}{Observed hemorrhagic dengue cases.}
#' }
#' @source Raw study file `inst/extdata/raw/study_data_dengue_hemorrhagic_fever.csv`.
"dengue_hermosillo_2010"

#' Calibration Example Inputs
#'
#' A compact list of package defaults used by examples and tests.
#'
#' @format A named list with 3 entries:
#' \describe{
#'   \item{parameters}{Named numeric vector returned by [default_parameters()].}
#'   \item{initial_state}{Named numeric vector returned by [default_state()].}
#'   \item{observation_window}{Integer vector defining a short calibration window.}
#' }
"calibration_example"
