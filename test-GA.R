# test function


GA_test <- function(y,x,family,mutation=0.01,fitness="rank",P=100,tol=0.01,maxIter=100L){
  
  # set up parameters
  C = ncol(x)
  P_ix = 1:P
  C_ix = 1:C
  
  # sanity checks
  if(family=="gassian"){
    
  }
  if(P < C | P > 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  }
  
  while(TRUE){
    
    
    
    
    # if(something){
      # break()
    # }
  }
}

# set up parameters
P = 50 # number of chromosomes
P_ix = 1:P # so we don't have to waste time making the vector over and over
C = ncol(data$X) # chromosome length (number of sites in genome)
C_ix = 1:C # so we don't have to waste time making the vector over and over
mutation = 0.01 # per-site, per-generation probability of mutation
fitness = "rank" # character in 'value' or 'rank'
family = "gaussian"