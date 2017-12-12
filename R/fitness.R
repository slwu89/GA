###############################################################################
#       _____ __        __ ___  __ __ _____      _________
#      / ___// /_____ _/ /|__ \/ // /|__  /_    / ____/   |
#     \__ \/ __/ __ `/ __/_/ / // /_ /_ <(_)  / / __/ /| |
#    ___/ / /_/ /_/ / /_/ __/__  __/__/ /    / /_/ / ___ |
#   /____/\__/\__,_/\__/____/ /_/ /____(_)   \____/_/  |_|
#
#   Fitness Evaluation
#   Luna Luan, Xinyu Liu, Nina Magnuson, Sean Wu
#   Statistics 243 GA Package
#
###############################################################################


#' Serial Fitness Evaluation
#'
#' Evaluate fitness in serial using \code{\link{vapply}} and return a numeric vector of results. Parallel fitness evaluation is provided through \code{\link{fitness_parallel}}
#'
#' @param pop list of chromosomes to evaluate fitness of
#' @param y vector of responses
#' @param x matrix of covariates
#' @param family error distribution family for \code{\link[stats]{glm}}
#' @param fitness_function fitness function passed from \code{\link{select}}
#'
#' @export
fitness_serial <- function(pop,y,x,family,fitness_function){
  # evaluate fitness of each chromosome
  pop_fitness = vapply(X = pop,FUN = function(xx, y, x, family){
    ix_mod = as.logical(xx)
    mod = stats::glm(y ~ x[,ix_mod], family=family)
    fitness_function(mod)
  },FUN.VALUE = numeric(1), y=y, x=x, family=family)
  return(pop_fitness)
}

#' Parallel Fitness Evaluation
#'
#' Evaluate fitness in parallel using \code{\link[parallel]{parLapply}} and return a numeric vector of results. Parallel evaluation
#' of chromosomes may be slower than calling serial \code{\link{fitness_serial}} unless P is very large.
#'
#' @param pop list of chromosomes to evaluate fitness of
#' @param y vector of responses
#' @param x matrix of covariates
#' @param family error distribution family for \code{\link[stats]{glm}}
#' @param fitness_function fitness function passed from \code{\link{select}}
#' @param ncores number of cores
#'
#' @export
fitness_parallel <- function(pop,y,x,family,fitness_function,ncores){
  cl = parallel::makePSOCKcluster(ncores)
  pop_fitness <- unlist(parallel::parLapply(cl, X=pop, fun = function(xx, y, x, family){
    ix_mod = as.logical(xx)
    mod = stats::glm(y ~ x[,ix_mod], family)
    fitness_function(mod)
  }, y, x, family))
  parallel::stopCluster(cl)
  return(pop_fitness)
}
