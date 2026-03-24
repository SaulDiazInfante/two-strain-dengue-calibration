#' Resolve a Sequence Specification
#'
#' @param spec A numeric vector or a list containing either `values`, or
#'   `start`/`end` with `by` or `length.out`.
#' @param allow_null Logical scalar indicating whether `NULL` is allowed.
#'
#' @return A numeric vector, or `NULL` when `allow_null` is `TRUE`.
#' @keywords internal
#' @noRd
resolve_sequence_spec <- function(spec, allow_null = FALSE) {
  if (is.null(spec)) {
    if (isTRUE(allow_null)) {
      return(NULL)
    }

    stop("A non-null sequence specification is required.", call. = FALSE)
  }

  if (is.atomic(spec)) {
    return(as.numeric(spec))
  }

  if (!is.list(spec)) {
    stop("Sequence specifications must be numeric vectors or lists.", call. = FALSE)
  }

  if (!is.null(spec[["values"]])) {
    return(as.numeric(unlist(spec[["values"]], use.names = FALSE)))
  }

  start <- spec[["start"]]
  end <- spec[["end"]]
  length_out <- spec[["length.out"]]

  if (is.null(length_out)) {
    length_out <- spec[["length_out"]]
  }

  if (is.null(start) || is.null(end)) {
    stop("Sequence specs must include `start` and `end`.", call. = FALSE)
  }

  if (!is.null(length_out)) {
    return(seq(start, end, length.out = length_out))
  }

  if (!is.null(spec[["by"]])) {
    return(seq(start, end, by = spec[["by"]]))
  }

  stop(
    "Sequence specs must include either `values`, `by`, or `length.out`.",
    call. = FALSE
  )
}
