# make synthetic data to validate our method

# beta2, beta5, beta6 are truly null
beta = c(1.52,12.5,0,9.12,21.34,0,0,88.65)

N = 20 # number of observations

x = matrix(NaN,nrow=N,ncol=length(beta)-1)
for(i in 1:ncol(x)){
  x[,i] = rnorm(n = N,mean = i,sd = 1.5)
}

y = (cbind(1,x) %*% beta) + rnorm(n = N)

mod = lm(y~x)
AIC(mod)

mod1 = lm(y~x[,which(beta[-1]!=0)])
AIC(mod1)
