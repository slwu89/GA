# GA test
rm(list=ls());gc()
set.seed(42L)

# generate true regression model
library(simrel)
N = 500 # number of observations
p = 100 # number of covariates
q = floor(p/4) # number of relevant predictors
m = q # number of relevant components

ix = sample(x = 1:p,size = m)
data = simrel(n=N, p=p, m=m, q=q, relpos=ix, gamma=0.2, R2=0.75)

mod = lm(data$Y~data$X)
mod1 = lm(data$Y~data$X[,ix])
AIC(mod)
AIC(mod1)

# set up parameters
P = 50 # number of chromosomes
P_ix = 1:P # so we don't have to waste time making the vector over and over
C = ncol(data$X) # chromosome length (number of sites in genome)
C_ix = 1:C # so we don't have to waste time making the vector over and over
mutation = 0.01 # per-site, per-generation probability of mutation
fitness = "rank" # character in 'value' or 'rank'
family = "gaussian"

# argument sanity checks here maybe?

# make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
pop = replicate(n = P,expr = {sample(x = c(0,1),size = C,replace = TRUE)},simplify = FALSE)
pop_fitness = vector(mode = "numeric",length = P)

# fitnesss
# this could be parallelized (but if X and Y too large, sending out data to cores might be too inefficient)
pop_fitness = vapply(X = pop,FUN = function(x,data,family){
  ix_mod = as.logical(x)
  mod = stats::glm(data$Y~data$X[,ix_mod],family)
  stats::AIC(mod)
},FUN.VALUE = numeric(1),data=data,family=family)

# sort by objective function (minimum of AIC)

# sort by rank (given by formula: 2*ri / P(P+1))
pop_rank = rank(-pop_fitness)
pop_rank_final = (2*pop_rank) / (P*(P+1))

# selection
selection_ix = sample(x = P_ix,size = 50,replace = TRUE,prob = pop_rank_final)
pop = pop[selection_ix]

# crossover
odd_seq = seq(from=1,to=P,by=2)
even_seq = seq(from=2,to=P,by=2)
# pop_new = pop # this is inefficient!

# # this could be parallelized (might not be any faster though)
# mapply(i=odd_seq,j=even_seq,FUN = function(i,j){
#   split = sample(x = C_ix,size = 1) # split point in chromosome
#   new_i = c(pop[[i]][1:split],pop[[j]][(split+1):C])
#   new_j = c(pop[[j]][1:split],pop[[i]][(split+1):C])
#   # c(list(new_i),list(new_j))
#   list(new_i,new_j)
# },SIMPLIFY = FALSE,USE.NAMES = FALSE)

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

pop = new_pop # keep on iterating
