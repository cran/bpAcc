
PAccept = function(xbar,sd,N,distribution = "normal", criteria = "SP10:2006"){
  
  E1 = pnorm(10,xbar,sd)
  E2 = pnorm(-10,xbar,sd)
  phat = E1 - E2
  
  P_accept = 1-pnorm(0.78,phat,sqrt((phat*(1-phat))/N))
  cat("-------------------------------------------------------\n\n")
  cat("The probability of acceptance as per SP10 is",P_accept)
  
  if(P_accept < 0.95){
    cat("\n\n The device is not meeting the SP10 criteria.\n")
    cat("------------------------------------------------------")
  }
  else if(P_accept >= 0.95){
    cat("\n\n The device is meeting the SP10 criteria.\n")
    cat("------------------------------------------------------")
  }
}