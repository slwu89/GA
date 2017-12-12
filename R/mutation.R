###############################################################################
#       _____ __        __ ___  __ __ _____      _________
#      / ___// /_____ _/ /|__ \/ // /|__  /_    / ____/   |
#     \__ \/ __/ __ `/ __/_/ / // /_ /_ <(_)  / / __/ /| |
#    ___/ / /_/ /_/ / /_/ __/__  __/__/ /    / /_/ / ___ |
#   /____/\__/\__,_/\__/____/ /_/ /____(_)   \____/_/  |_|
#
#   Mutation
#   Luna Luan, Xinyu Liu, Nina Magnuson, Sean Wu
#   Statistics 243 GA Package
#
###############################################################################

#' Mutation of Chromosomes in Genetic Algoritm
#'
#' Mutation randomly determines if and where on the chromosome, or model, a mutation of one of the alleles, or covariates occurs.
#' This is used in the genetic algorithm to introduce randomness not present in the genetic makeup of the population.
#' Mutation is according to a Poisson process because in the usual case the number of alleles on the chromosome, C,
#' is large and probability of mutation is small. The number of mutations is first sampled by \code{\link[stats]{rpois}} and then
#' the mutations are scattered randomly across the chromosome.
#'
#' @param new_pop a nested list of chromosomes representing the child population. The the result of selection and crossover.
#' @param C a value representing the population size. Should be of class integer, corresponding to the number of columns of input matrix x for initial regression.
#' @param C_ix a sequence of integers from 1 to number of models considered for mutation (C).
#' @param mutation a value indicating the rate at which mutation should occur. Input as rate parameter lambda, or rpois for generating indices of chromosome list where mutation occurs.
#'
#' @export
mutation <- function(new_pop,C,C_ix,mutation){
  for(k in 1:length(new_pop)){
    # sample how many mutations on this chromosome
    mutation_N = rpois(n = 1, lambda = mutation*C)
    if(mutation_N > C){ # if more mutations than alleles, resample
      mutation_N = rpois(n = 1, lambda = mutation*C)
    }
    if(mutation_N<1) { # no mutations this chromosome
      next()
    }
    # scatter the mutations across this chromosome
    mutation_ix = sample(x = C_ix,size = mutation_N,replace = FALSE) # we don't assume same site can mutate twice; has biological plausibility (probability on order of 1/N^2 too small to consider)
    new_pop[[k]][mutation_ix] = 1 - new_pop[[k]][mutation_ix] #changes value 0 or 1
  }

  # check if any all-0 chromosome (possible on small genomes)
  c_sums = sapply(X = new_pop,FUN = sum)
  if(any(c_sums==0)){
    ix = which(c_sums==0)
    for(i in ix){
      new_pop[[i]][sample(x = C_ix,size = 1,replace = FALSE)] = 1L
    }
  }
  return(new_pop)
}
