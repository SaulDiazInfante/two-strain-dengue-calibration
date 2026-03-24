#' Build a Timestamped Output File Path
#'
#' @param output_dir Base directory for generated outputs.
#' @param subdir Child directory created under `output_dir`.
#' @param stem File stem appended after the timestamp.
#' @param extension File extension without a leading dot.
#' @param timestamp Timestamp prefix used in generated file names.
#'
#' @return A character scalar containing the full output path.
#' @keywords internal
#' @noRd
build_output_path <- function(output_dir, subdir, stem, extension, timestamp) {
  directory <- file.path(output_dir, subdir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, sprintf("%s__%s.%s", timestamp, stem, extension))
}
