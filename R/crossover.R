###############################################################################
#       _____ __        __ ___  __ __ _____      _________
#      / ___// /_____ _/ /|__ \/ // /|__  /_    / ____/   |
#     \__ \/ __/ __ `/ __/_/ / // /_ /_ <(_)  / / __/ /| |
#    ___/ / /_/ /_/ / /_/ __/__  __/__/ /    / /_/ / ___ |
#   /____/\__/\__,_/\__/____/ /_/ /____(_)   \____/_/  |_|
#
#   Crossover
#   Luna Luan, Xinyu Liu, Nina Magnuson, Sean Wu
#   Statistics 243 GA Package
#
###############################################################################

#' Crossover of Chromosomes in Genetic Algorithm
#'
#' Crossover is used in genetic algorithms for optimization when two parent chromosomes recombine to produce a child chromosome for the following generation.
#' Crossover randomly picks a location on the chrosome to split the genetic alleles of both parents and then recombines the segments of genes to produce a new chromosome.
#'
#' @param pop a nested list represented the parent population of models or chromosomes
#' @param P_combn a list of all unique pairwise combinations of chromosome index (see \code{\link[utils]{combn}})
#' @param C a value representing chromosome size (number of 'alleles'). Should be of class integer, corresponding to the number of columns of input matrix x.
#'
#' @export
crossover <- function(pop,P_combn,P,C){
  new_pop = pop
  pairs = P_combn[sample(x = 1:length(P_combn),size = P,replace = FALSE)]

  j = 1
  for(i in pairs){
    split = sample(x = 1:(C-1), size = 1) # split point in chromosome
    new_pop[[j]] = c(pop[[i[1]]][1:split], pop[[i[2]]][(split+1):C])
    j = j + 1
  }
  return(new_pop)
}
