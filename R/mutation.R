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
#' mutation randomly determines if and where on the chromosome, or model, a mutation of one of the genes, or covariates occurs. This is used in the genetic algorithm to introduce randomness not present in the genetic makeup of the population.
#'
#' mutation(new_pop,C,C_ix,mutation)
#'
#' @param new_pop a nested list of chromosomes representing the child population. The the result of selection and crossover.
#' @param C a value representing the population size. Should be of class integer, corresponding to the number of columns of input matrix x for initial regression.
#' @param  C_ix a sequence of integers from 1 to number of models considered for mutation (C).
#' @param  mutation a value indicating the rate at which mutation should occur. Input as rate parameter lambda, or rpois for generating indices of chromosome list where mutation occurs.
#' Author(s):
#' Examples (?)
#' @export
#'
mutation<-function(new_pop,C,C_ix,mutation){
  for(k in 1:length(new_pop)){
    # sample how many mutations on this chromosome
    mutation_N = rpois(n = 1, lambda = mutation*C) #rpois for rarity
    if(mutation_N<1) {
      next()
    } #no mutations this chromosome
    # scatter the mutations across this chromosome
    mutation_ix = sample(x = C_ix,size = mutation_N,replace = FALSE) # we don't assume          same site can mutate twice
    new_pop[[k]][mutation_ix] = 1 - new_pop[[k]][mutation_ix] #changes value 0 or 1
  }

  # check if any all-0 chromosome
  c_sums = sapply(X = new_pop,FUN = sum)
  if(any(c_sums==0)){
    ix = which(c_sums==0)
    for(i in ix){
      cat("DEBUG: putting a random one in the chromosome!\n")
      new_pop[[i]][sample(x = C_ix,size = 1,replace = FALSE)] = 1L
    }
  }
  return(new_pop)
}
