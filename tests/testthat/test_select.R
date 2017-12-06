context("select")

library(simrel)
set.seed(42L)
N = 500 # number of observations
p = 100 # number of covariates
q = floor(p/4) # number of relevant predictors
m = q # number of relevant components

ix = sample(x = 1:p,size = m)
data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)
y = data$Y
x = data$X

test_that('select works',
          expect_true(select(y = y, x = x,family = "gaussian", maxIter = 10)))

test_that('select fails',
          expect_false(select(y = y, x = x, family = "gaussian", maxIter = 10)))

test_that('select input warning',
          expect_warning(select(y = y,x = "something",family = "gaussian", maxIter = 10)),
          expect_warning(select(y = rep(-1,10), x = x,family = "gamma", maxIter = 10)),
          expect_warning(select(y = y,x = x, P = c(10,20), family = "gamma", maxIter = 10)),
          expect_warning(select(y = y,x = x, fitness_function= "something", family = "gamma", maxIter = 10))
)