
#' select is used to optimize model selection. It can be used to carry out regression utilizing the genetic algorithm to return the optimal model based on some criteria, default: AIC.
#' select(y, x, family, mutation = .01, ncores = 0, fitness_function = stats::AIC, P = 100, tol = .0005, maxIter = 100L)
#'
#' @param y a vector or matrix of responses
#' @param x a matrix of covariates
#' @param family a description of the error distribution to be used in the glm fitting. Default is to use gaussian.
#' @param mutation an optional value specifying the rate at which mutation occurs in the new population.
#' @param ncores an optional value indicating the number of cores to use in parallelization. Should be numeric.
#' @param fitness_function optional function to evaluate model fitness. Must be of type closure. Default is to use AIC.
#' @param P Population size.
#' @param tol Optional value indicating relative convergence tolerance. Should be numeric.
#' @param maxIter Optional value indicating the maximum number of iterations. Defaults to 100.
#' @export
