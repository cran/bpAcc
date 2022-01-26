##########################################################################
# These functions are 
# Copyright (C) 1998-2021 T. Chandel, V. Miranda,  A. Lowe, . Lee, 
# Auckland University of Technology.
# All rights reserved.


## Below is the CDF of sigmaami  P(sigaami0.85 < any value)
#psigmaami(sigmaami = 5, mu = 4 , std.dev = 6, lower.tail = TRUE)


toint <- function(xbar, mean, sd, n) {
  dnorm(xbar, mean = mean, sd = sd/sqrt(n))
}

psigmaami <- function(sigmaami, mu, std.dev, n, ptolerror = 0.85, 
                      lower.tail = TRUE) {
  
  mu <- abs(mu)
  
  if (lower.tail) {
    myret <- integrate(toint, lower = rootinv(sigmaami, ptolerror), 
                       upper = Inf,
                       mean = mu, sd = std.dev, n = n)$value 
    
    return(myret)
    
  } else {
    myret <- integrate(toint, lower = -Inf, 
                       upper = rootinv(sigmaami, ptolerror),
                       mean = mu, sd = std.dev, n = n)$value 
    
    return(myret)
    
  }
  cat("Stop. Something went wrong.")
}

