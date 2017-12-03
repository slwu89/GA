# test function


GA_test <- function(y,x,family,P=100,tol=0.01,maxIter=100L){
  
  C = ncol(x)
  
  # sanity checks
  if(P < C | P > 2*C){
    cat("P ",P," not within suggested population size range C <= P <= 2C\n")
  }
  
  while(TRUE){
    
    
    
    
    # if(something){
      # break()
    # }
  }
}