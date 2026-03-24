if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install the 'testthat' package to run the test suite.")
}

testthat::test_check("twostraindengue")
