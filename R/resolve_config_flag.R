#' Resolve a Logical Config Flag
#'
#' @param value A logical scalar or a character scalar such as `"interactive"`,
#'   `"true"`, or `"false"`.
#'
#' @return A logical scalar.
#' @keywords internal
#' @noRd
resolve_config_flag <- function(value) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    return(value)
  }

  if (is.character(value) && length(value) == 1L) {
    lower_value <- tolower(value)

    if (lower_value == "interactive") {
      return(interactive())
    }

    if (lower_value %in% c("true", "false")) {
      return(lower_value == "true")
    }
  }

  stop(
    "Config flags must be logical values or one of 'interactive', 'true', or 'false'.",
    call. = FALSE
  )
}
