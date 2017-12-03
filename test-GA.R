# test function


GA_test <- function(y,x,family,mutation=0.01,fitness_function=stats::AIC,fitness="rank",P=100,tol=0.01,maxIter=100L){
  
  # set up parameters (stuff that is made each iteration just make once and save)
  C = ncol(x)
  P_ix = 1:P 
  C_ix = 1:C
  odd_seq = seq(from=1,to=P,by=2) 
  even_seq = seq(from=2,to=P,by=2)
  stop_condition = FALSE
  
  # sanity checks
  if(family=="gassian" & all(y %% 1 == 0)){cat("vector of integer responses but family 'Gaussian' error distribution selected\n")}
  if(family=="Gamma" & any(y < 0)){stop("family 'Gamma' error distribution selected but have negative responses")}
  if(P < C | P > 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  }
  if(!fitness %in% c("rank","weight")){stop("'fitness' must be either 'rank' or 'weight'")}
  if(typeof(fitness_function)!="closure"){stop("fitness_function must be a function that returns an objective value to minimize")}
  
  # make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
  new_pop = pop = replicate(n = P,expr = {sample(x = c(0,1),size = C,replace = TRUE)},simplify = FALSE)
  pop_fitness = vector(mode = "numeric",length = P)
  
  i = 0
  while(TRUE){
    i = i + 1
    
    # evaluate fitness of each chromosome
    pop_fitness = vapply(X = pop,FUN = function(x,data,family){
      ix_mod = as.logical(x)
      mod = stats::glm(data$Y~data$X[,ix_mod],family=family)
      fitness_function(mod)
    },FUN.VALUE = numeric(1),data=data,family=family)
    
    # the best fitness
    best_ix = which.min(pop_fitness)
    best_fitness = min(pop_fitness)
    best_chromosome = pop[best_ix]
    
    # selection
    
    if(fitness=="rank"){
      # # sort by objective function (minimum of AIC)
      # new_minimum=min(pop_fitness)
      # ### here I use 0.01 as weights instead of 0.5 from the paper
      # relative=lapply(pop_fitness,function(x)exp(-0.01*(x-minimum)))
      # weights=lapply(relative,function(x)x/sum(unlist(relative)))
      # #selection
      # selection_ix = sample(x = P_ix,size = P,replace = TRUE,prob = weights)
      # pop = pop[selection_ix]
      
      
      # if(something){
      #  stop_condition = TRUE
      # }
      } else {
      # sort by rank (given by formula: 2*ri / P(P+1))
      pop_rank = rank(-pop_fitness)
      pop_rank_final = (2*pop_rank) / (P*(P+1))
      # selection
      selection_ix = sample(x = P_ix,size = P,replace = TRUE,prob = pop_rank_final)
      pop = pop[selection_ix]
    }
    
    # crossover
    new_pop = pop
    j = 1
    for(i in 1:length(odd_seq)){
      split = sample(x = C_ix,size = 1) # split point in chromosome
      new_pop[[j]] = c(pop[[odd_seq[i]]][1:split],pop[[even_seq[i]]][(split+1):C])
      j = j + 1
      new_pop[[j]] = c(pop[[even_seq[i]]][1:split],pop[[odd_seq[i]]][(split+1):C])
      j = j + 1
    }
    
    # mutation
    for(i in 1:length(new_pop)){
      # sample how many mutations on this chromosome
      mutation_N = rpois(n = 1,lambda = mutation*C)
      # scatter the mutations across this chromosome
      mutation_ix = sample(x = C_ix,size = mutation_N,replace = FALSE) # we don't assume same site can mutate twice
      new_pop[[i]][mutation_ix] = 1-new_pop[[1]][mutation_ix]
    }
    
    pop = new_pop

    cat("iteration ",i,"\n")
    if(i >= maxIter | stop_condition){
      break()
    }
  }
  
  # return(NULL) # give back something
  return(list(
    best_chromosome=best_chromosome,
    best_fitness = best_fitness
  ))
}


#  make up some data
library(simrel)
N = 500 # number of observations
p = 100 # number of covariates
q = floor(p/4) # number of relevant predictors
m = q # number of relevant components

ix = sample(x = 1:p,size = m)
data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)
y = data$Y
x = data$X

out = GA_test(y = y,x = x,family = )
