\name{bpAcc-methods}
%\name{undocumented-methods}
\docType{methods}
%\alias{ccoef,ANY-method}
%\alias{ccoef-method}
%
%
% 201607, 201608:
%\alias{psint,pvgam-method}
%\alias{summary,pvgam-method}
%\alias{show,summary.pvgam-method}
%\alias{df.residual,pvgam-method}
\alias{show,ProbAcceptance-method}
\alias{plot,ProbAcceptance-method}
%
%
% This does not work (need one line for each one):
  %\alias{trplot,qrrvglm,uqo-method}
%
%
%
%\title{ Undocumented Methods Functions in \pkg{bpAcc} }
\title{Methods functions in \pkg{bpAcc}}

\description{
  %Undocumented methods functions are aliased here.
  %In the \pkg{bpAcc} package there are currently some
  %objects/methods/classes which are currently internal and/or
  %undocumented. The help file suppresses the warnings when the package is
  %'CHECK'ed.
  The function \code{\link[bpAcc:ProbAccept]{ProbAccept}}
  returns an object of (S4) class '\code{ProbAcceptance}', a list
  with elements required by methods defined 
  for \code{ProbAcceptance} objects. See below.
  
}
%\usage{
  %  \S4method{ccoef}{rrvgam,Coef.rrvgam,rrvglm,qrrvglm,
    %                   Coef.rrvglm,Coef.qrrvglm}(object, ...)
  %}

\section{Methods}{

  Methods for objects of class '\code{ProbAcceptance}'.
  \describe{
    
    \item{\code{show}}{
      To display the probability of tolerable error in a friendly version.
      }
      \item{\code{plot}}{
      To plot the estimated probabiliy of tolerable error based on simulated data,
      as described in \code{\link[bpAcc:ProbAccept]{ProbAccept}}
      
    }
  }
}
\section{Examples}{

# Creates an object of class 'ProbAcceptance'
fit1 <- ProbAccept(n = 85, mu = 5, sd = 5, ptolerror = 0.85)

# Methods 'show' and 'plot'
# plot(fit1)
# fit1
}
\keyword{methods}
\keyword{classes}
%\keyword{ ~~ other possible keyword(s)}
\keyword{models}
\keyword{regression}
\keyword{internal}
