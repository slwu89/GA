library(GA)
context("select")

y <- as.matrix(mtcars$mpg)
x <- as.matrix(mtcars[2:ncol(mtcars)])
C <- ncol(x) #number of chromosomes (models considered)
C_ix <- 1:C #chromosome sequence
pop <- replicate(n = 10, expr = {
  chromosome = sample(x = c(0,1), size = C, replace = TRUE)
  while(all(chromosome==0)){ # check that we dont get an all-0 chromosome in small problems
    chromosome[sample(x = C_ix,size = C,replace = FALSE)] = 1L
  }
  return(chromosome)
}, simplify = FALSE)

test_that('serial fitness works',
          expect_is(fitness_serial(pop = pop,y = y,x = x,family = "gaussian",fitness_function = stats::AIC), "numeric")
          )

test_that('parallel fitness works',
          expect_is(fitness_parallel(pop = pop,y = y,x = x,family = "gaussian",fitness_function = stats::AIC,ncores = 2), "numeric")
)


