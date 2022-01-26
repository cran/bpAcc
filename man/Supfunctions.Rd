%\name{f}
\name{Supplementary functions}
\alias{fsup1}
\alias{fsup2}
\alias{toint}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Supplementary fuctions of the package \pkg{bpAcc}.
}
\description{
Supplementary function required to compute 
the acceptance probability of a device for BP measurement.
}
\usage{

  fsup1(mu, sd.aami, ptolerror)
  fsup2(mu.aami, sdev, ptolerror)
  toint(xbar, mean, sd, n)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
 True mean error (bias). This is 'mu' of the errors distribution (normal).
}
\item{sd.aami}{
 Standard deviation of a normal distribution meeting the ANSI/AAMI-SP10 standards.
}
\item{ptolerror}{
  Probability of tolerable error. Default is 0.85.
}
\item{mu.aami}{
  Mean of a normal distribution meeting the ANSI/AAMI-SP10 standards.
}
\item{sdev}{
  True standard deviation of the errors (normal distribution).
}
\item{xbar, mean, sd, n}{
Auxiliary arguments.
}


}
\details{
Supplementary functions to compute 
the \emph{probability of accepting the device D}.
The errors (average of the differences from
three \emph{device} measurements and the three corresponding reference
readings) are by default normally distributed with mean
\code{mean} and \code{sd}. Function \code{fsup1} is called 
by \code{\link[bpAcc:root]{root}}, \code{fsup2} is
needed by \code{\link[bpAcc:rootinv]{rootinv}}, and 
\code{toint} is required by \code{\link[bpAcc:psigmaami]{psigmaami}}.

See \code{\link[bpAcc:psigmaami]{psigmaami}} for further
details on \code{ptolerror}.

}
\value{

No value returned.
These are funtcion internally called by 
\code{\link[bpAcc:ProbAccept]{ProbAccept}}.

}
\references{
Chandel, T. and Lee, TC. and Lowe, A. and Miranda, V. (2022)
Blood Pressure Device Accuracy Evaluation: Statistical
Considerations with an Implementarion in R,
\emph{Statistical Methods in Medical Research} (under review).
}
\author{
Tanvi Chandel, Tet-Chuan Lee, Andrew Lowe, Victor Miranda.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[bpAcc:ProbAccept]{ProbAccept}}
  \code{\link[bpAcc:psigmaami]{psigmaami}}
}
\examples{
#### 'sd' that meets the AAMI-sp10 standars
root(mu = 5, ptolerror = 0.9)
# root(mu = 11, ptolerror = 0.85) # Error

#### 'mu' that meets the AAMI-sp10 standars
rootinv(sdev = 3, ptolerror = 0.9)
# rootinv(sdev = 2, ptolerror = 0.85) # Error

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
