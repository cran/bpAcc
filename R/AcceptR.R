phat_rev = function(N){
  rev_phat = uniroot(function(phat) {
    E1 <- phat - 1.645*(sqrt((phat*(1-phat))/N))
    return(E1-0.78)
  },
  interval = c(0, 1))$root
}





AcceptR = function(n,distribution = "normal",criteria = "SP10:2006"){
  
  N <- n
  cat("-----------------------------------")
  cat("\n")
  cat("\n For",N,"samples,",phat_rev(N)*100, "% of errors must be within -10 mmHg to 10 mmHg \n")
  cat("\n")
  cat("-----------------------------------\n")
  roots2 <- function(xbar){
    uniroot(
      func2 <- function(s){
        delta <- 10 
        dbar <- xbar
       
        
        E1 <- pnorm(-delta,dbar,s)
        E2 <- pnorm(delta,dbar,s)
        return((E2-E1)-phat_rev(N))
      },interval = c(0.0001, 10))$root
   
  }
  
  
  xbar = seq(0, 5, by=0.5)
  sd = mapply(roots2,xbar)
  df <- data.frame(xbar,sd)
  #print(phat_rev(N))
 
  print(df)
}

