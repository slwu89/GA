# mutation
#' 
#'
#' here is my function
#' @param 
#'
#' @export
#' 
mutation<-function(new_pop,C,C_ix,mutation){
  #temp<-new_pop
  #C_ix = 1:C
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
  return(new_pop)
}
