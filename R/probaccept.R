##########################################################################
# These functions are 
# Copyright (C) 1998-2021 T. Chandel, V. Miranda,  A. Lowe, . Lee, 
# Auckland University of Technology.
# All rights reserved.

if (FALSE)
ProbAccept <- function(n, mu, sd, ptolerror = 0.85) {
  mu <- abs(mu)
  
  if (mu >= 9.999)
    stop("Unable to find 'sd' with mu = ", mu,
         ", and ptolerror = ", ptolerror, ".")
  
  if (sd <= 0)
    stop("Wrong input for 'sd'.")
  
  if ( (ptolerror <=0 ) || (ptolerror >=1))
    stop("Wrong input for 'ptolerror'.")
  
  integrate(innerIntegral3, lower = 0.0000, upper = root(0),
            mu = mu, std.dev= sd, n = n, ptolerror = ptolerror)$value
  
}


innerIntegral3 <- function(mysd, mu, std.dev, n, ptolerror) {
  
  dnorm(mysd, mean = std.dev * sqrt((2*n - 3)/(2*(n - 1))),
        sd = sqrt(std.dev^2 / (2 * (n - 1)  ))) *
    sapply(mysd, function(mysd) {
      psigmaami(mysd, mu = mu, std.dev = std.dev, n = n, 
                ptolerror = ptolerror, lower.tail = FALSE)
    })
}


ProbTolError <- function(distribution = "normal", 
                         mu, std.dev, delta) {
  
  check_in <- match.arg(distribution, c("normal"))
  pnorm(q = delta, mean = mu, sd = std.dev, lower.tail = TRUE) -
    pnorm(q = -delta, mean = mu, sd = std.dev, lower.tail  = TRUE)   
  
}





#### normal distribution #####
fsup1 <- function(mu, sd.aami, ptolerror){
  delta <- 10
  dbar <- mu
  E1 <- pnorm(-delta,dbar,sd.aami)
  E2 <- pnorm(delta,dbar,sd.aami)
  return( (E2 - E1) - ptolerror)
}

root <- function(mu, ptolerror = 0.85) {
  

#    stop("Unable to find 'sd' with mu = ", mu,
#         ", and ptolerror = ", ptolerror, ".")
  
  out <- tryCatch(
    {uniroot(function(sd.aami, mu) fsup1(mu = mu, sd.aami, ptolerror = ptolerror),
             interval = c(0.0001, 10), mu = mu)$root
    },
    error=function(cond) {
      return(NA) # cat ("Unable to find roots for simulation with mu", mu, "\n", sep = " ")
    })    
  
  if (is.na(out)) {
    #  if (mu >= 9.999)
    stop("Unable to find 'sd' with mu = ", mu,
         ", and ptolerror = ", ptolerror, ".")
    
  } else{
    return(out)
  }
  
}

##############################
fsup2 <- function(mu.aami, sdev, ptolerror) {
  delta <- 10
  E1 <- pnorm(-delta, mu.aami, sdev)
  E2 <- pnorm(delta, mu.aami, sdev)
  
  E2 - E1 - ptolerror
}

rootinv  <- function(sdev, ptolerror = 0.85) {
  


  out <- tryCatch(
    {uniroot(function(mu.aami) fsup2(mu.aami, sdev, ptolerror = ptolerror),
             interval = c(1e-08, 10))$root
    },
    error=function(cond) {
      return(NA) # cat ("Unable to find roots for simulation with mu", mu, "\n", sep = " ")
    })    
  
  if(is.na(out)) {
    #  if (sdev > 6.9467046755)  # approx root(0)
    stop("Unable to find 'mu' with sdev = ", sdev,
         ", and ptolerror = ", ptolerror, ".")
  } else{
    return(out)
  }
  
}
