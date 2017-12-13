context("select")

y <- as.matrix(mtcars$mpg)
x <- as.matrix(mtcars[2:11])

test_that('selection works',
          expect_that(select(y = y, x = x, k = 2, family = "gaussian", P=10, maxIter = 100)$count, is_less_than(100)))

test_that('select fails',
          expect_that(select(y = y, x = x, k = 2, family = "gaussian", P=10, maxIter = 2)$warning, equals("Reached Maximum Iterations prior to convergence")))

test_that('select input error', {
          expect_error(select(y = y,x = "something", k = 2, family = "gaussian", maxIter = 10))
          expect_error(select(y = rep(-1,10), x = x, k = 2, family = "gamma", maxIter = 10))
          expect_error(select(y = y,x = x, P = c(10,20), k = 2, family = "gaussian", maxIter = 10))
          expect_error(select(y = y,x = x, k = 2, fitness_function= "something", family = "gaussian", maxIter = 10))
          expect_error(select(y = y,x = x, k = 2, fitness= "other that weight or rank", family = "gaussian", maxIter = 10))
})


