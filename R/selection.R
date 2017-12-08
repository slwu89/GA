###############################################################################
#       _____ __        __ ___  __ __ _____      _________
#      / ___// /_____ _/ /|__ \/ // /|__  /_    / ____/   |
#     \__ \/ __/ __ `/ __/_/ / // /_ /_ <(_)  / / __/ /| |
#    ___/ / /_/ /_/ / /_/ __/__  __/__/ /    / /_/ / ___ |
#   /____/\__/\__,_/\__/____/ /_/ /____(_)   \____/_/  |_|
#
#   Selection
#   Luna Luan, Xinyu Liu, Nina Magnuson, Sean Wu
#   Statistics 243 GA Package
#
###############################################################################

#' Selection Of Models for Recombination
#'
#' Standard selection is used to determine parent chromosomes fitness for recombination in the genetic algorithm for model optimization.
#' Fitness is assessed based on AIC or an inputted fitness function (must be function to minimize).
#' AIC is used in two ways: one through rank \eqn{\frac{2\ r_{i}}{P(P+1)}} and another through AIC weights.
#' Weight takes the absolute fitness to determine probability of recombination while rank uses relative ranking of fitness of models.
#'  * Selection of chromosomes that contribute their genes to generation t+1 done by resampling the population with probability proportional to weight or rank, with replacement.
#'
#' @param pop_fitness a vector indicating fitness of each chrosomosome in the current generation based on the evaluated fitness function.
#' @param fitness a character indicating whether to use rank or weight when comparing model fitness.
#' @param P Population size, corresponding to the number of genes or covariates on each chromosome or model considered for recombination.
#' @param P_ix a sequence of integers, ending value P.
#'
#' References:
#' Khalid Jebari, Mohammed Madiafi (2013) Selection Methods for Genetic Algorithms. https://www.researchgate.net/publication/259461147_Selection_Methods_for_Genetic_Algorithms.
#' @export
selection <- function(pop_fitness,fitness,P,P_ix){
  if(fitness == "weight"){
    # sort by objective function (minimum of AIC)
    minimum=min(pop_fitness)
    relative = vapply(X = pop_fitness,FUN = function(x){exp(-0.01*(x-minimum))},FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    weights = vapply(X = relative,FUN = function(x){x/sum(relative)},FUN.VALUE = numeric(1))
    # selection
    selection_ix = sample(x = P_ix, size = P, replace = TRUE, prob = weights)
  }
  else {
    # sort by rank (given by formula: 2*ri / P(P+1))
    pop_rank = rank(-pop_fitness)
    pop_rank_final = (2*pop_rank) / (P*(P+1))
    # selection
    selection_ix = sample(x = P_ix, size = P, replace = TRUE, prob = pop_rank_final)
  }

  return(selection_ix)
}


#' Tournament Selection Of Models for Recombination
#'
#' Standard selection is used to determine parent chromosomes fitness for recombination in the genetic algorithm for model optimization.
#' Fitness is assessed based on AIC or an inputted fitness function (must be function to minimize).
#' AIC is used in two ways: one through rank \eqn{\frac{2\ r_{i}}{P(P+1)}} and another through AIC weights.
#'  1. in tournament selection the population is first divded into a random partition of subsets of equal size k
#'  2. then within each partition we select the most fit chromosome for propagation in the new population
#'  3. iterate until we have a new population equal in size to the old population
#'
#' @param pop_fitness a vector indicating fitness of each chrosomosome in the current generation based on the evaluated fitness function.
#' @param fitness a character indicating whether to use rank or weight when comparing model fitness.
#' @param P Population size, corresponding to the number of genes or covariates on each chromosome or model considered for recombination.
#' @param P_ix a sequence of integers, ending value P.
#' @param k size of disjoin subset (must be smaller than P and whose quotient with P is 0, eg; P mod k = 0)
#'
#' Author(s):
#' References:
#' Khalid Jebari, Mohammed Madiafi (2013) Selection Methods for Genetic Algorithms. https://www.researchgate.net/publication/259461147_Selection_Methods_for_Genetic_Algorithms.
#' @export
#'
selection_tournament <- function(pop_fitness,fitness,P,P_ix,k){
  if(P %% k != 0 | k >= P){stop("tournament selection parameter k must be less than P and have P mod k = 0")}
  if(fitness == "weight"){
    # sort by objective function (minimum of AIC)
    minimum=min(pop_fitness)
    relative = vapply(X = pop_fitness,FUN = function(x){exp(-0.01*(x-minimum))},FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    weights = vapply(X = relative,FUN = function(x){x/sum(relative)},FUN.VALUE = numeric(1))
  } else {
    # sort by rank (given by formula: 2*ri / P(P+1))
    pop_rank = rank(-pop_fitness)
    weights = (2*pop_rank) / (P*(P+1))
  }
  # implement tournament
  selection_ix = rep(NaN,P)
  while(any(is.nan(selection_ix))){
    random_partition = split(x = P_ix,f = sample(x = P_ix,size = P/k,replace = FALSE))
    champions = vapply(X = random_partition,FUN = function(x){
      x[order(weights[x],decreasing = TRUE)][1]
    },FUN.VALUE = integer(1),USE.NAMES = FALSE)
    i = which(is.nan(selection_ix))[1]
    selection_ix[i:(i+k-1)] = champions
  }
  return(selection_ix)
}
