#' Select utilizes the genetic algoritim to optimize model selection. Relying on functions: selection, crossover and mutation, new populations of chromosomes correspondng to models are generated, selecting based on AIC, until the optimal model is achieved. 
#'
#' select(y,x,family,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="rank",P=100,tol=0.0005,maxIter=100L)
#' 
#' @param y a vector or matrix of responses
#' @param x a matrix of covariates
#' @param family a description of the error distribution to be used in the glm fitting. Default is to use gaussian.
#' @param mutation an optional value specifying the rate at which mutation occurs in the new population.
#' @param ncores an optional value indicating the number of cores to use in parallelization. Should be numeric.
#' @param fitness_function optional function to evaluate model fitness. Must be of type closure. Default is to use AIC. Options: "rank" uses relative rank of models based on AIC. Option: "weight" uses absolute value of AIC to determine probability of reproduction in the preceding generation. This option should be used with caution because it can become stuck at a local minimum, as a single model with very low AIC will have large probability of reproduction. 
#' @param P Population size, corresponding to the number of genes or covariates on each chromosome or model. 
#' @param tol Optional value indicating relative convergence tolerance. Should be of class numeric.
#' @param maxIter Optional value indicating the maximum number of iterations. Default is 100.
#' Author(s):
#' References
#' Burnham, K. P. and D. R. Anderson. 2002. Model Selection and Multimodel Inference. Springer-Verlag, New York
#' 
#' Geof H. Givens, Jennifer A. Hoeting (2013) Combinatorial Optimization (italicize). Chapter 3 of Computational Statistics (italicize). 
#' 
#' 
#' Examples:
#' @export
#' 

select <- function(y,x,family,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="rank",P=100,tol=0.0005,maxIter=100L){
  
  # sanity checks
  # check x
  if (!is.matrix(x)) stop("x should be a matrix of numbers")
  if (!is.numeric(x)) stop("x should be a matrix of numbers")
  
  # check y
  if (!is.matrix(y)) stop("y should be a matrix of numbers")
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
  
  # check maximum iteration
  if (length(maxIter) != 1) stop("Please provide only one maximum iteration")
  if (!is.numeric(maxIter)) stop("Maximum iteration should be a number")
  if (round(maxIter) != maxIter) stop("Maximum iteration should be an integer")
  
  # check tolerance rate
  if (length(tol) != 1) stop("Please provide only one convergence rate")
  
  # check fitness criterion
  if(!fitness %in% c("rank","weight")){stop("'fitness' must be either 'rank' or 'weight'")}
  if(typeof(fitness_function)!="closure"){stop("fitness_function must be a function that returns an objective value to minimize")}
  
  # check ncores
  if(!is.numeric(ncores)){stop("the nubmer of cores used in a parallelization should be a number")}
  if(ncores<0){stop("the number of cores used in a parallelization should be non negative")}
  
  
  # set up parameters (stuff that is made each iteration just make once and save)
  C = ncol(x) #number of chromosomes (models considered)
  P_ix = 1:P #population sequence 
  C_ix = 1:C #chromosome sequence
  odd_seq = seq(from=1,to=P,by=2) 
  even_seq = seq(from=2,to=P,by=2)
  stop_condition = FALSE
  
  # check C
  if(P < C | P > 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  }
  
  # make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
  new_pop = pop = replicate(n = P, expr = {sample(x = c(0,1), size = C, replace = TRUE)}, simplify = FALSE)
  pop_fitness = vector(mode = "numeric", length = P)
  
  #initialize stopping condition
  #old_fitness = 1
  i = 0
  hist_fit=vector(mode="numeric")
  while(TRUE){
    i = i + 1
    if(ncores == 0){
      # evaluate fitness of each chromosome
      pop_fitness = vapply(X = pop,FUN = function(xx, y, x, family){
        ix_mod = as.logical(xx)
        mod = stats::glm(y ~ x[,ix_mod], family=family)
        fitness_function(mod)
      },FUN.VALUE = numeric(1), y=y, x=x, family=family)
    } 
    else if(ncores>0){ #evaluate fitness in parallel 
      require(parallel)
      nCores <- ncores
      cl <- makeCluster(nCores) # by default this uses sockets
      pop_fitness <- unlist(parLapply(cl, X=pop, fun = function(xx, y, x, family){
        #browser()
        ix_mod = as.logical(xx)
        mod = stats::glm(y ~ x[,ix_mod], family)
        fitness_function(mod)
      }, y, x, family))
      stopCluster(cl)
      #cat('para',"\n",pop_fitness,"\n",i)
    }
    
    # chromosome/model with the best fitness
    best_ix = which.min(pop_fitness)
    best_fitness = min(pop_fitness)
    best_chromosome = pop[best_ix]
    
    #stopping condition
    #cat("iteration ",i,"\n",best_fitness,"\n")
    #if( abs((best_fitness-old_fitness) / old_fitness) < tol){ #relative change in fitness
      #stop_condition=TRUE
    #}
    #old_fitness=min(hist_fit)
    #hist_fit=c(hist_fit,best_fitness)
    
    #stopping condition
    cat("iteration ",i,"\n",best_fitness,"\n")
    if(i>1){
    if ((abs(min(hist_fit)-best_fitness)/abs(best_fitness))< tol){ 
      #print(i)  
      stop_condition=TRUE
      }
    }
    if(i>=maxIter | stop_condition){
      index=unlist(best_chromosome)
      break()
    }
    hist_fit<-c(hist_fit,best_fitness)
    #if(i >= maxIter | stop_condition){
    #  best_fitness=old_fitness
    #  best_chromosome=old_chromosome
    #  index=unlist(best_chromosome)
    #  break()
    #}
    
    #reiterate if not stopping
    #old_fitness=best_fitness
    #old_chromosome=best_chromosome
    
    selection_ix<-selection(pop_fitness,fitness,P,P_ix)
    pop = pop[selection_ix]
    # crossover
    new_pop<-crossover(pop,odd_seq,even_seq,C)
    
    # mutation
    pop<-mutation(new_pop,C,C_ix,mutation)
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
