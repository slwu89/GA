###############################################################################
#       _____ __        __ ___  __ __ _____      _________
#      / ___// /_____ _/ /|__ \/ // /|__  /_    / ____/   |
#     \__ \/ __/ __ `/ __/_/ / // /_ /_ <(_)  / / __/ /| |
#    ___/ / /_/ /_/ / /_/ __/__  __/__/ /    / /_/ / ___ |
#   /____/\__/\__,_/\__/____/ /_/ /____(_)   \____/_/  |_|
#
#   Main GA Function
#   Luna Luan, Xinyu Liu, Nina Magnuson, Sean Wu
#   Statistics 243 GA Package
#
###############################################################################

#' Genetic Algorithm based Variable Selection
#'
#' Select utilizes the genetic algoritim to optimize model selection for generalized linear models.
#' Relying on functions: \code{\link{fitness_serial}}, \code{\link{selection}}, \code{\link{crossover}} and \code{\link{mutation}}, new populations of chromosomes correspondng to models are generated, selecting based on AIC, until the algorithm is assumed to converge based on user-specified tolerance or maximum iterations are reached.
#'
#'  1. Burnham, K. P. and D. R. Anderson. 2002. Model Selection and Multimodel Inference. Springer-Verlag, New York
#'  2. Geof H. Givens, Jennifer A. Hoeting (2013) Combinatorial Optimization (italicize). Chapter 3 of Computational Statistics (italicize).
#'
#' @param y a vector or matrix of responses
#' @param x a matrix of covariates
#' @param family a description of the error distribution to be used in the glm fitting. Default is to use gaussian.
#' @param k size of disjoin subset (must be smaller than P and whose quotient with P is 0, eg; P mod k = 0)
#' @param P population size (number of chromosomes to evaluate)
#' @param mutation an optional value specifying the rate at which mutation occurs in the new population.
#' @param ncores an optional value indicating the number of cores to use in parallelization. Should be numeric.
#' @param fitness_function optional function to evaluate model fitness. Must be of type closure. Default is to use AIC. Options: "rank" uses relative rank of models based on AIC. Option: "weight" uses absolute value of AIC to determine probability of reproduction in the preceding generation. This option should be used with caution because it can become stuck at a local minimum, as a single model with very low AIC will have large probability of reproduction.
#' @param fitness character in 'rank' or 'weight' used in \code{\link{selection}} or \code{\link{selection_tournament}}; select between using AIC weight or ranks to evaluate fitness. For general problems (eg; if fitness_function is \code{\link[stats]{BIC}}) using 'rank' is preferred
#' @param selection character in 'fitness' or 'tournament' to select either standard selection \code{\link{selection}} or tournament selection \code{\link{selection_tournament}}
#' @param tol Optional value indicating relative convergence tolerance. Should be of class numeric.
#' @param maxIter Optional value indicating the maximum number of iterations. Default is 100.
#'
#' @examples
#' # GA test with mtcars
#' rm(list=ls());gc()
#' set.seed(42L)
#'
#' # generate data from mtcars
#' y <- as.matrix(mtcars$mpg)
#' x <- as.matrix(mtcars[2:11])
#'
#' \dontrun{select(y = y, x = x, k = 2, family = "gaussian", P=10, maxIter = 10)}
#'
#' # GA test with simrel data
#' library(simrel)
#'
#' N = 500 # number of observations
#' p = 100 # number of covariates
#' q = floor(p/4) # number of relevant predictors
#' m = q # number of relevant components
#' ix = sample(x = 1:p,size = m) # location of relevant components
#' data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)
#'
#' fitness = "rank" # character in 'value' or 'rank'
#' family = "gaussian"
#' y = data$Y
#' x = data$X
#'
#' \dontrun{GAoptim = GA::select(y = y,x = x,family = family,k = 20,P = 100,ncores = 0,fitness = "rank",selection = "tournament")}
#'
#' @export
select <- function(y,x,family,k,P,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="rank",selection="fitness",tol=0.0005,maxIter=100L){

  # sanity checks
  # check x
  if (!is.matrix(x)) stop("x should be a matrix of numbers")
  if (!is.numeric(x)) stop("x should be a matrix of numbers")

  # check y
  if (!(is.vector(y) | is.matrix(y))) stop("y should be a matrix or a vector of numbers")
  if (!is.numeric(y)) stop("y should be a vector of numbers")

  # check type of regression
  if(family=="gaussian" & all(y %% 1 == 0)){cat("Warning: vector of integer responses but family 'Gaussian' error distribution selected\n")}
  if(family=="gamma" & any(y < 0)){stop("Error: family 'Gamma' error distribution selected with negative responses")}

  # check mutation rate
  if (length(mutation) != 1) stop("Please provide only one mutation rate")
  if (mutation < 0 | mutation > 1) stop("The mutation rate should be between 0 and 1")

  # check population size
  if (length(P) != 1) stop("Please provide only one population size")
  if (!is.numeric(P)) stop("Population size should be a number")
  if (round(P) != P) stop("Population size should be an integer")

  # check maximum iteration
  if (length(maxIter) != 1) stop("Please provide only one maximum iteration")
  if (!is.numeric(maxIter)) stop("Maximum iteration should be a number")
  if (round(maxIter) != maxIter) stop("Maximum iteration should be an integer")

  # check tolerance rate
  if (length(tol) != 1) stop("Please provide only one convergence rate")

  # check fitness criterion
  if(!fitness %in% c("rank","weight")){stop("'fitness' must be either 'rank' or 'weight'")}
  if(typeof(fitness_function)!="closure"){stop("fitness_function must be a function that returns an objective value to minimize")}

  # check ncores
  if(!is.numeric(ncores)){stop("the nubmer of cores used in a parallelization should be a number")}
  if(ncores<0){stop("the number of cores used in a parallelization should be non negative")}
  if(ncores>0){
    ncores = min(ncores,parallel::detectCores()) # ncores should not be greater than the number of cores available
  }

  # tournament selection; check that choice of k will produce valid disjoint subsets of the P chromosomes
  if(!selection %in% c("tournament","fitness")){stop("selection must be either 'fitness' or 'tournament'")}
  if(selection=="tournament"){
    check = (P%%k != 0 | k >= P)
    if(length(check)==0 | check){
      stop("k must be chosen such that it is smaller than P and whose quotient with P is 0, eg; P mod k = 0")
    }
  }

  # set up parameters once rather than on each iteration
  C = ncol(x) #number of chromosomes (models considered)
  P_ix = 1:P #population sequence
  C_ix = 1:C #chromosome sequence
  P_combn = combn(x = 1:P,m = 2,simplify = FALSE) # used during crossover
  stop_condition = FALSE
  hist_fit = vector(mode="numeric")

  # check C: just warn the user, don't stop them from doing this
  if(P <= C | P >= 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C, ",C," < P <= ",2*C,"\n")
  }

  # make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
  new_pop = pop = replicate(n = P, expr = {
    chromosome = sample(x = c(0,1), size = C, replace = TRUE)
    while(all(chromosome==0)){ # check that we dont get an all-0 chromosome in small problems
      chromosome[sample(x = C_ix,size = C,replace = FALSE)] = 1L
    }
    return(chromosome)
  }, simplify = FALSE)
  pop_fitness = vector(mode = "numeric", length = P)

  # main loop over generations
  old_fitness = 1
  i = 0
  while(TRUE){
    i = i + 1

    # evaluate fitness
    if(ncores > 0){
      pop_fitness = fitness_parallel(pop,y,x,family,fitness_function,ncores)
    } else {
      pop_fitness = fitness_serial(pop,y,x,family,fitness_function)
    }

    # fittest chromosome
    best_ix = which.min(pop_fitness)
    best_fitness = min(pop_fitness)
    best_chromosome = pop[best_ix]

    # stopping condition
    cat("iteration",i,"\n",'best fitness score ',best_fitness,"\n")
    if(i>1){
      if((abs(min(hist_fit)-best_fitness)/abs(best_fitness))< tol){
        stop_condition=TRUE
      }
    }

    # if stopping criteria met record the best chromosome and break loop
    if(i >= maxIter | stop_condition){
      index = unlist(best_chromosome)
      break()
    # otherwise record best fitness value and continue
    } else {
      hist_fit = c(hist_fit,best_fitness)
    }

    # selection
    if(selection=="fitness"){
      selection_ix = selection(pop_fitness,fitness,P,P_ix)
    } else {
      selection_ix = selection_tournament(pop_fitness,fitness,P,P_ix,k)
    }
    pop = pop[selection_ix]

    # crossover
    new_pop = crossover(pop,P_combn,P,C)

    # mutation
    pop = mutation(new_pop,C,C_ix,mutation)

  } # end loop

  # give the user a warning if returning simply due to reaching maxIter
  if(i == maxIter){
    out = list(
      best_chromosome = best_chromosome[[1]],
      best_fitness = best_fitness,
      best_model = stats::glm(y~x[,as.logical(index)],family=family),
      count = i,
      warning = "Reached Maximum Iterations prior to convergence"
    )
  } else {
    out = list(
      best_chromosome = best_chromosome[[1]],
      best_fitness = best_fitness,
      best_model = stats::glm(y~x[,as.logical(index)],family=family),
      count = i
    )
  }

  return(out)
}
