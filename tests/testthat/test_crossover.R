library(GA)
context("crossover")

P_combn_test = combn(x = 1:100,m = 2,simplify = FALSE)
selection_ix_test = c(45,27,80,86,85,20,67,54,13,13,54,74,67,67,13,3,42,57,48,17,31,59,3,43,12,100,76,37,38,45,49,80,45,56,55,40,18,67,57,38,38,41,97,24,2,46,40,13,65,38,40,70,49,80,31,51,38,37,45,88,83,52,95,62,49,31,91,54,64,2,70,29,96,86,32,86,58,43,33,40,43,78,83,46,83,74,73,7,93,90,7,50,100,16,24,96,24,46,97,58)
pop_test = replicate(n = 100, expr = {sample(x = c(0,1), size = 100, replace = TRUE)}, simplify = FALSE)
pop_test = pop_test[selection_ix_test]

test_that('crossover works', {
  expect_that(crossover(pop_test, P_combn_test, 100, 100), is_a("list"))
  expect_length(crossover(pop_test, P_combn_test, 100, 100), 100)
  expect_that(crossover(pop_test, P_combn_test, 100, 100)[[1]], is_a("numeric"))
  expect_length(crossover(pop_test, P_combn_test, 100, 100)[[1]], 100)
})

test_that('crossover input error', {
  expect_error(crossover(pop = "character", P_combn_test, 100, 100))
})
