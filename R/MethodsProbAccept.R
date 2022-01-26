##########################################################################
# These functions are 
# Copyright (C) 1998-2021 T. Chandel, V. Miranda,  A. Lowe, . Lee, 
# Auckland University of Technology.
# All rights reserved.



#setClass("test", slots = list(a = "numeric", b = "list"))
#not <- new("test", a = rnorm(20), b = list(x = 5))
#test <- function(n) new("test", a = rnorm(n), b = list(x = 9))
#setMethod("plot", signature(x = "test", y = "ANY"), function(x, y, ...) plot(x@a))



setClass("ProbAcceptance", slots = 
           list(probabilities = "numeric",
                extra = "list"))

ProbAccept <- function(n, mu, sd, ptolerror = 0.85,
                       distribution = "normal",
                       criteria = "SP10:2006",
                       simulate = FALSE, sim.count = 1e4,
                       noshow = FALSE) {
  
  if (!is.logical(simulate))
    stop("Wrong input for 'simulate'.")
  
  
  check_crit <- match.arg(criteria, c("SP10:2006"))
  check_in <- match.arg(distribution, c("normal"))
  
  if ((sim.count < 0) || !(sim.count = sim.count))
    stop("Bad input for 'sim.count'.")
  
 # if (sd > root(0, ptolerror = ptolerror))
#    stop("Unable to find 'mu' with sd = ", sd, ", and ptolerror =",
#         ptolerror, ".")
  
  if (mu >= rootinv(0.0001, ptolerror = ptolerror))
    stop("'sigmaAAMI' not defined for mu > ",
         round(rootinv(0.0001, ptolerror = ptolerror), 3),
         ", subject to ptolerror = ", ptolerror, ".")
  
  if (sd <= 0)
    stop("Wrong input for 'sd'.")
  
  if ( (ptolerror <=0 ) || (ptolerror >=1))
    stop("Wrong input for 'ptolerror'.")
  
  exProb <- integrate(innerIntegral3, lower = 0.0000, upper = root(0),
            mu = mu, std.dev= sd, n = n, ptolerror = ptolerror)$value
  
  
  if (simulate) {
    std.dev <- sd
    pass.rates <- sim.mu <- sim.sd <- NULL
    lower.limit <- -10
    upper.limit <- 10

    for(i in 1:length(mu)){
      for (j in 1:sim.count){
        data <- rnorm(n = n, mean = mu[i], sd = std.dev[i])
        pass.rate <- sum(data >= lower.limit && data <= upper.limit) / n
        pass.rates <- append(pass.rates, pass.rate)
        sim.mu <- append(sim.mu, mean(data))
        sim.sd <- append(sim.sd, sd(data))
      }
    }
    
    atest <- TRUE
    if (atest) {
      sim.sd.aami = mapply(root,sim.mu, ptolerror)
    } else {
      sim.sd.aami = root(mu = mu) 
    }
    
    
    audit.df <- data.frame(Mean=sim.mu, StdDev=sim.sd, StdDevaami=sim.sd.aami)
    audit.df['Pass'] = ifelse (audit.df$StdDevaami >= audit.df$StdDev, 1, 0)
    audit.check = round(sum(audit.df['Pass'], na.rm = TRUE) / 
                          sum(!is.na(audit.df['Pass'])), 3)
    
    out <- new("ProbAcceptance", 
               probabilities = audit.check,
               extra = list(sim.sd = sim.sd,
                            sim.sd.aami = sim.sd.aami,
                            simulate = simulate,
                            n = n,
                            noshow = noshow,
                            ptolerror = ptolerror,
                            sim.count = sim.count,
                            mu = mu, std.dev = sd))
    

    
  } else {
    
    
    out <- new("ProbAcceptance", 
               probabilities = exProb,
               extra = list(sim.sd = NULL,
                            sim.sd.aami = NULL,
                            simulate = simulate,
                            n = n, 
                            noshow = noshow,
                            sim.count = sim.count,
                            mu = mu, std.dev = sd))
    
  }
  
  
  out@extra[["exactProb"]] <- exProb
  
  if (!noshow) {
    message("\n 
                     --- Exact results ---
            \n
            The exact probability of accepting the device is ", exProb, 
            ".
            \n")
  #  cat("           --- Exact results ---")
  #  cat("\n")
  #  cat("The exact probability of accepting the device is ", exProb, ".")
  #  cat("\n")
    
  }
 
  
 
  out
  
}


setMethod("show", signature = signature(object = "ProbAcceptance"),
          function(object) {
            
            
            if (object@extra$simulate) {
              cat("\n")
              cat("\n")
              cat("      --- Simulation results ---")
              cat("\n")
              cat("Based on n =", object@extra$sim.count,
                  "samples from a normal distribution, the probability"
                  ,"of accepting the device is", 
                  round(object@probabilities, 5),".")
              
              cat("\n")  
              cat("Input: sample size (i.e., 'n') = ", object@extra$n, ";",
                  "'mu' = ", object@extra$mu, "  ; 'sd' = ",
                  object@extra$std.dev,
                  " ; ptolerror =", object@extra$ptolerror, ".")
              
              cat("\n")
              cat("\n")
              
              cat("** Run plot() for this object to get the scatter of the simulated samples.")
              cat("\n")
            } else {
              cat("\n")
              cat(">>> No simulations to run.")
              cat("\n")
            }
           
            
          })

setMethod("plot", signature = signature(x ="ProbAcceptance", y = "ANY"), 
      function(x, y, ...) {
        
        if (x@extra$simulate) {
          
          sim.sd <- x@extra$sim.sd
          sim.sd.aami <- x@extra$sim.sd.aami
          #plot(sim.sd, sim.sd.aami, ylim = c(0, 10), xlim = c(0, 10),
          #     main = "Scatter SD-AAMI SP10 vs. Sample SDs.",
          #     ylab = "SD-AAMI SP10)", xlab = "Sample standard deviations")
          
          xx <- c(0, 2*(1:5), 2*(5:1), 0)
          yy <- c(0, 2*(1:5), 0, 0, 0, 0, 0 ,0)
          plot(xx, yy, type = "l",
               ylab = "SD-AAMI SP10", xlab = "Sample standard deviations",
               main = "Scatter SD-AAMI SP10 vs. Sample SDs.", ...)
          polygon(xx, yy, col = "#65BFFF")
          points(sim.sd.aami, sim.sd, ylim = c(0, 10), xlim = c(0, 10),
                 pch = 19, cex = 0.65)
          abline(a = 0, b =1, lwd = 2, col = "Red")
          text(x = 7, y = 2, paste("ACCEPTANCE REGION ( mu = ",
                                   x@extra$mu, ")."))
          text(x = 7, y = 1, paste("Prob. accepting device (simulated) =",
                                   x@probabilities))
          text(x = 2, y = 8, paste("Prob. tolerable error = ", 
                                   x@extra$ptolerror))
        } else {
          cat("\n")
          cat("No plot to show. Set 'simulate = TRUE' to get simulated results.")
          cat("\n")
        }
           
            
      })

