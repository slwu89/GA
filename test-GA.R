# test function


GA_test <- function(y,x,family,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="rank",P=100,tol=0.0005,maxIter=100L){
  # browser()
  # set up parameters (stuff that is made each iteration just make once and save)
  C = ncol(x) #number of chromosomes (models considered)
  P_ix = 1:P #population sequence 
  C_ix = 1:C #chromosome sequence
  odd_seq = seq(from=1,to=P,by=2) 
  even_seq = seq(from=2,to=P,by=2)
  stop_condition = FALSE
  
  # sanity checks
  # check x
  if (!is.matrix(x)) stop("x should be a matrix of numbers")
  if (!is.numeric(x)) stop("x should be a matrix of numbers")
  
  # check y
  if (!(is.matrix(y) | is.vector(y))) stop("y should be a matrix or a vector of numbers")
  if (!is.numeric(y)) stop("y should be a vector of numbers")
  
  # check type of regression
  if(family=="gassian" & all(y %% 1 == 0)){cat("Warning: vector of integer responses but family 'Gaussian' error distribution selected\n")}
  if(family=="Gamma" & any(y < 0)){stop("Error: family 'Gamma' error distribution selected with negative responses")}
  
  # check mutation rate
  if (length(mutation) != 1) stop("Please provide only one mutation rate")
  if (mutation < 0 | mutation > 1) stop("The mutation rate should be between 0 and 1")
  
  # check population size
  if (length(P) != 1) stop("Please provide only one population size")
  if (!is.numeric(P)) stop("Population size should be a number")
  if (round(P) != P) stop("Population size should be an integer")
  if(P < C | P > 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  }
  
  # check maximum iteration
  if (length(maxIter) != 1) stop("Please provide only one maximum iteration")
  if (!is.numeric(maxIter)) stop("Maximum iteration should be a number")
  if (!is.integer(maxIter)) stop("Maximum iteration should be an integer")
  
  # check tolerance rate
  if (length(tol) != 1) stop("Please provide only one convergence rate")
  
  # check fitness criterion
  if(!fitness %in% c("rank","weight")){stop("'fitness' must be either 'rank' or 'weight'")}
  if(typeof(fitness_function)!="closure"){stop("fitness_function must be a function that returns an objective value to minimize")}
  
  # check ncores
  if(!is.numeric(ncores)){stop("the nubmer of cores used in a parallelization should be a number")}
  if(ncores<0){stop("the number of cores used in a parallelization should be non negative")}
  
  
  # make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
  new_pop = pop = replicate(n = P, expr = {sample(x = c(0,1), size = C, replace = TRUE)}, simplify = FALSE)
  pop_fitness = vector(mode = "numeric", length = P)
  
  #initialize stopping condition
  old_fitness = 1
  i = 0
  while(TRUE){
    i = i + 1
    
    if(ncores == 0){
      # evaluate fitness of each chromosome
      pop_fitness = vapply(X = pop,FUN = function(xx, y, x, family){
        ix_mod = as.logical(xx)
        mod = stats::glm(y ~ x[,ix_mod], family=family)
        fitness_function(mod)
      },FUN.VALUE = numeric(1), y=y, x=x, family=family)
    } else if(ncores>0){ #evaluate fitness in parallel 
      require(parallel)
      nCores <- ncores
      cl <- makeCluster(nCores) # by default this uses sockets
      pop_fitness <- unlist(parLapply(cl, X=pop, fun = function(xx, y, x, family){
        #browser()
        ix_mod = as.logical(xx)
        mod = stats::glm(y ~ x[,ix_mod], family)
        stats::AIC(mod)
      }, y, x, family))
      stopCluster(cl)
    }
    # chromosome/model with the best fitness
    best_ix = which.min(pop_fitness)
    best_fitness = min(pop_fitness)
    best_chromosome = pop[best_ix]
    
    #stopping condition
    cat("iteration ",i,"\n",best_fitness,"\n")
    if( abs((best_fitness-old_fitness) / old_fitness) < tol){ #relative change in fitness
      stop_condition=TRUE
    }
    
    if(i >= maxIter | stop_condition){
      best_fitness=old_fitness
      best_chromosome=old_chromosome
      index=unlist(best_chromosome)
      break()
    }
    
    #reiterate if no stopping
    old_fitness=best_fitness
    old_chromosome=best_chromosome
    
    if(fitness == "weight"){
      # sort by objective function (minimum of AIC)
      relative = lapply(pop_fitness, function(x) exp(-0.01*(x-minimum))) #minimum?
      weights  = lapply(relative, function(x) x / sum(unlist(relative)))
      # selection
      selection_ix = sample(x = P_ix, size = P, replace = TRUE, prob = weights)
      pop = pop[selection_ix]
    } else {
      # sort by rank (given by formula: 2*ri / P(P+1))
      pop_rank = rank(-pop_fitness)
      pop_rank_final = (2*pop_rank) / (P*(P+1))
      # selection
      selection_ix = sample(x = P_ix, size = P, replace = TRUE, prob = pop_rank_final)
      #likelihood for selection based on rank
      pop = pop[selection_ix] 
    }
    
    # crossover
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
    
    # mutation
    for(k in 1:length(new_pop)){
      # sample how many mutations on this chromosome
      mutation_N = rpois(n = 1, lambda = mutation*C) #rpois for rarity 
      if(mutation_N<1) {next()} #no mutations this chromosome
      # scatter the mutations across this chromosome
      mutation_ix = sample(x = C_ix,size = mutation_N,replace = FALSE) # we don't assume          same site can mutate twice
      new_pop[[k]][mutation_ix] = 1 - new_pop[[k]][mutation_ix] #changes value 0 or 1
    }
    
    pop = new_pop 
  }
  
  # return(NULL) # give back something
  #goes all the way to max iteration, failing message?
  return(list(
    best_chromosome = best_chromosome,
    best_fitness = best_fitness,
    best_model = stats::glm(y~x[,as.logical(index)],family=family),
    count = i #number of iterations
  ))
}


#  make up some data
library(simrel)
set.seed(42L)
N = 500 # number of observations
p = 100 # number of covariates
q = floor(p/4) # number of relevant predictors
m = q # number of relevant components

ix = sample(x = 1:p,size = m)
data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)
y = data$Y
x = data$X

out = GA_test(y = y,x = x,family = "gaussian")

best_mod = stats::glm(y~x[,as.logical(out$best_chromosome[[1]])],family = "gaussian")
out<-GA_test(y = y,x = x,family = "gaussian",ncores=4)

mod = lm(y~x)
mod1 = lm(y~x[,ix])
AIC(mod)
AIC(mod1)
AIC(best_mod)
AIC(out$best_model)