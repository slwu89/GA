#' Selection Of Models for Recombination
#' Selection is used to determine parent chromosomes fitness for recombination in the genetic algorithm for model optimization. Fitness is assessed based on AIC or an inputted fitness function. AIC is used in two ways: one through rank and another through weight. Weight takes the absolute fitness to determine probability of recombination while rank uses relative ranking of fitness of models. This fitness assessment determines likelihood to 'reproduce' or carry on genetic makeup into the t+1 generation. 
#' 
#' selection(pop_fitness,fitness,P,P_ix)
#'
#' 
#' @param pop_fitness a vector indicating fitness of each chrosomosome in the current generation based on the evaluated fitness function. 
#' @param fitness a character indicating whether to use rank or weight when comparing model fitness.
#' @param P Population size, corresponding to the number of genes or covariates on each chromosome or model considered for recombination. 
#' @param P_ix a sequence of integers, ending value P. 
#' 
#' Author(s):
#' References: 
#' Khalid Jebari, Mohammed Madiafi (2013) Selection Methods for Genetic Algorithms. https://www.researchgate.net/publication/259461147_Selection_Methods_for_Genetic_Algorithms. 
#' @export
#' 


selection <- function(pop_fitness,fitness,P,P_ix){
  if(fitness == "weight"){
    # sort by objective function (minimum of AIC)
    minimum=min(pop_fitness)
    relative = lapply(pop_fitness, function(x) exp(-0.01*(x-minimum))) #minimum?
    weights  = lapply(relative, function(x) x / sum(unlist(relative)))
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
