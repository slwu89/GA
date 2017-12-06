# selection
#' 
#'
#' here is my function
#' @param 
#'
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
