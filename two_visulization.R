plot_ga<-function(y,x,family,k,P,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="rank",selection="fitness",tol=0.0005,maxIter=100L){
  # set up parameters (stuff that is made each iteration just make once and save)
  C = ncol(x) #number of chromosomes (models considered)
  P_ix = 1:P #population sequence
  C_ix = 1:C #chromosome sequence
  # odd_seq = seq(from=1,to=P,by=2)
  # even_seq = seq(from=2,to=P,by=2)
  P_combn = combn(x = 1:P,m = 2,simplify = FALSE)
  stop_condition = FALSE
  hist_fit=vector(mode="numeric")

  # check C
  #if(P <= C | P >= 2*C){
  #  cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  #}

  # make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
  new_pop = pop = replicate(n = P, expr = {
    chromosome = sample(x = c(0,1), size = C, replace = TRUE)
    if(all(chromosome==0)){
      chromosome[sample(x = C_ix,size = 1,replace = FALSE)] = 1L
    }
    return(chromosome)
  }, simplify = FALSE)
  pop_fitness = vector(mode = "numeric", length = P)

  #initialize stopping condition
  old_fitness = 1
  i = 0
  out = vector(mode="list",length=maxIter)

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
    }
    out[[i]]$AIC = pop_fitness
    out[[i]]$chromosomes = pop
    # chromosome/model with the best fitness
    best_ix = which.min(pop_fitness)
    best_fitness = min(pop_fitness)
    best_chromosome = pop[best_ix]

    #stopping condition
    #cat("iteration ",i,"\n",best_fitness,"\n")
    #if( abs((best_fitness-old_fitness) / old_fitness) < tol){ #relative change in fitness
    #  stop_condition=TRUE
    #}

    #stopping condition
    cat("iteration",i,"\n",'Fitness Scores',best_fitness,"\n")
    #if(i>1){
    #  if ((abs(min(hist_fit)-best_fitness)/abs(best_fitness))< tol){
    #    #print(i)
    #    stop_condition=TRUE
    #  }
    #}

    if(i >= maxIter | stop_condition){
      #best_fitness=old_fitness
      #best_chromosome=old_chromosome
      #index=unlist(best_chromosome)
      break()
    }

    hist_fit<-c(hist_fit,best_fitness)

    #reiterate if not stopping
    #old_fitness=best_fitness
    #old_chromosome=best_chromosome

    # selection
    switch(selection,
           fitness = {selection_ix = selection(pop_fitness,fitness,P,P_ix)},
           tournament = {selection_ix = selection_tournament(pop_fitness,fitness,P,P_ix,k)},
           {stop("invalid 'selecion' argument: must be character in 'fitness' or 'tournament'")}
    )

    pop = pop[selection_ix]
    # crossover
    new_pop<-crossover(pop,P_combn,P,C)

    # mutation
    pop<-mutation(new_pop,C,C_ix,mutation)


  }

  return(out)

}

x=mtcars[,c(2:11)]
x=as.matrix(x)
y=as.matrix(mtcars$mpg)
out = plot_ga(y,x,family='gaussian',k=20,P=5,mutation=0.05,ncores=0,fitness_function=stats::AIC,fitness="weight",selection="fitness",tol=0.0005,maxIter=5L)

library(ggplot2)
library(reshape2)

datMelt = lapply(out,function(x){x$AIC})

datMelt = melt(datMelt)
colnames(datMelt) = c("AIC","Generation")
historical<-lapply(out,function(x){min(x$AIC)})
historical<-melt(historical)
colnames(historical)=c('minimum','Generation')
data<-merge(historical,datMelt,by='Generation')
data$Generation<-factor(data$Generation)
p1<-ggplot(data, aes(x=Generation, y=AIC, col=Generation)) +
    geom_point()+
    geom_point(aes(y=minimum), color="blue")+
  #scale_color_grey()
    scale_color_brewer(palette='RdGy')+
  ggtitle('Fitness evolution with 5 iterations,\n Data Source: MTCars')
p1
library(simrel)
set.seed(2L)
N = 500 # number of observations
p = 100 # number of covariates
q = floor(p/4) # number of relevant predictors
m = q # number of relevant components

ix = sample(x = 1:p,size = m)
data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)
y = data$Y
x = data$X
out = plot_ga(y,x,family='gaussian',k=20,P=50,mutation=0.01,ncores=0,fitness_function=stats::AIC,fitness="weight",selection="fitness",tol=0.0005,maxIter=10L)

datMelt = lapply(out,function(x){x$AIC})

datMelt = melt(datMelt)
colnames(datMelt) = c("AIC","Generation")
historical<-lapply(out,function(x){min(x$AIC)})
historical<-melt(historical)
colnames(historical)=c('minimum','Generation')
data<-merge(historical,datMelt,by='Generation')
data$Generation<-factor(data$Generation)
p2<-ggplot(data, aes(x=Generation, y=AIC, col=Generation)) +
  geom_point()+
  geom_point(aes(y=minimum), color="blue")+
  #scale_color_grey()
  scale_color_brewer(palette='RdGy')+
  ggtitle('Fitness evolution with 10 iterations,\n Data Source: use simrel package to generate')
p2
