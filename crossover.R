#' Crossover of Chromosomes in Genetic Algorithm 
#' Crossover is used in genetic algorithm for optimization when two parent chromosomes recombine to produce a child chromosome for the following generation. Crossover randomly picks a location on the chrosome to split the genetic alleles of both parents and then recombines the segments of genes to produce a new chromosome. 
#' 
#' crossover<-function(pop,odd_seq,even_seq,C)
#' 
#' @param pop a nested list represented the parent population of models or chromosomes
#' @param odd_seq a sequence of odd integers. Up to P, corresponding to the the number of genes per chromosome (or covariates per model). Used for splitting of the chromosome for recombination
#' @param even_seq a sequence of even integers. Used for splitting the chromosome for recombination.
#' @param C a value representing the population size (number of models). Should be of class integer, corresponding to the number of columns of input matrix x for initial regression. 
#' Author(s)
#' Examples (?)
#' @export
#' 
crossover<-function(pop,odd_seq,even_seq,C){
  new_pop = pop
  j = 1
  for(k in 1:length(odd_seq)){
    split = sample(x = 1:(C-1), size = 1) # split point in chromosome
    #recombination
    new_pop[[j]] = c(pop[[odd_seq[k]]][1:split], pop[[even_seq[k]]][(split+1):C])
    j = j + 1
    new_pop[[j]] = c(pop[[even_seq[k]]][1:split], pop[[odd_seq[k]]][(split+1):C])
    j = j + 1
  }
  return(new_pop)
}