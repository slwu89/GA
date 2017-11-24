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
C = ncol(data$X) # chromosome length (number of sites in genome)
mutation = 0.01 # per-site, per-generation probability of mutation
fitness = "rank" # character in 'value' or 'rank'
family = "gaussian"

# argument sanity checks here maybe?

# make a population of candidate solutions (use a list because we can apply over it quickly with vapply,lapply,sapply...unlike a matrix)
pop = replicate(n = P,expr = {sample(x = c(0,1),size = C,replace = TRUE)},simplify = FALSE)
pop_fitness = vector(mode = "numeric",length = P)

# fitness
pop_fitness = vapply(X = pop,FUN = function(x,data,family){
  ix_mod = as.logical(x)
  mod = stats::glm(data$Y~data$X[,ix_mod],family)
  stats::AIC(mod)
},FUN.VALUE = numeric(1),data=data,family=family)

if(fitness=="rank"){
  pop_fitness = 2*rank(pop_fitness) / P*(P+1)
}

# selection
pop = pop

# mutation