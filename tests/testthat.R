Sys.setenv("R_TESTS" = "")
library(testthat)
library(GA)

test_check("GA")
