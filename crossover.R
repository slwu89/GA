# Crossover
#' 
#'
#' here is my function
#' @param 
#'
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
