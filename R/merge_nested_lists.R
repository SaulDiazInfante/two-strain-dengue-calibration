#' Recursively Merge Nested Lists
#'
#' @param base Base list.
#' @param override Override list.
#'
#' @return A recursively merged list.
#' @keywords internal
#' @noRd
merge_nested_lists <- function(base, override) {
  if (is.null(override)) {
    return(base)
  }

  if (!is.list(base) || !is.list(override)) {
    return(override)
  }

  for (name in names(override)) {
    if (name %in% names(base)) {
      base[[name]] <- merge_nested_lists(base[[name]], override[[name]])
    } else {
      base[[name]] <- override[[name]]
    }
  }

  base
}
