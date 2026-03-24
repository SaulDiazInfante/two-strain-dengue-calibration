#' Drop NULL Entries from Nested Lists
#'
#' @param x A list-like object.
#'
#' @return `x` with `NULL` entries removed recursively.
#' @keywords internal
#' @noRd
drop_null_entries <- function(x) {
  if (!is.list(x)) {
    return(x)
  }

  null_entries <- vapply(x, is.null, logical(1))
  x <- x[!null_entries]

  for (name in names(x)) {
    if (is.list(x[[name]])) {
      x[[name]] <- drop_null_entries(x[[name]])
    }
  }

  x
}
