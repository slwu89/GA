#' Crossover of Chromosomes in Genetic Algorithm 
#' Crossover is used in genetic algorithm for optimization when two parent chromosomes recombine to produce a child chromosome for the following generation. Crossover randomly picks a location on the chrosome to split the genetic alleles of both parents and then recombines the segments of genes to produce a new chromosome. 
#' 
#' crossover(pop,odd_seq,even_seq,C)
#' 
#' @param pop a nested list represented the parent population of models or chromosomes
#' @param P_combn all combinations of 
#' @param C a value representing the population size (number of models). Should be of class integer, corresponding to the number of columns of input matrix x for initial regression. 
#' Author(s)
#' Examples (?)
#' @export
#' 
crossover<-function(pop,P_combn,P,C){
  new_pop = pop
  pairs = P_combn[sample(x = 1:length(P_combn),size = P,replace = FALSE)]
  
  j = 1
  for(i in pairs){
    split = sample(x = 1:(C-1), size = 1) # split point in chromosome
    new_pop[[j]] = c(pop[[i[1]]][1:split], pop[[i[2]]][(split+1):C])
    j = j + 1
  }
  
  
  # j = 1
  # for(k in 1:length(odd_seq)){
  #   split = sample(x = 1:(C-1), size = 1) # split point in chromosome
  #   #recombination
  #   new_pop[[j]] = c(pop[[odd_seq[k]]][1:split], pop[[even_seq[k]]][(split+1):C])
  #   j = j + 1
  #   new_pop[[j]] = c(pop[[even_seq[k]]][1:split], pop[[odd_seq[k]]][(split+1):C])
  #   j = j + 1
  # }
  return(new_pop)
}