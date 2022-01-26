%\name{f}
\name{psigmaami}
\alias{psigmaami}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution function of the sigma-aami transformation (not vectorized).
}
\description{
Distribution function of the sigma-aami transformation (not vectorized).
}
\usage{

  psigmaami(sigmaami, mu, std.dev, n, ptolerror = 0.85, lower.tail = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sigmaami}{

A single, positive, quantile. See below for further details.
 
}
\item{mu, std.dev}{
mean and standar deviation of the average errors distribution (normal).
}
\item{n}{
sample size
}
\item{ptolerror}{
Probability of tolerable error, Default is 0.85
}
\item{lower.tail}{
logical; it \emph{TRUE} (default), probabilities are
\eqn{P[sigmaami \le x]} otherwise, \eqn{P[sigmaami > x]}.
}
}
\details{
This is the distribution function of the \code{sigmaami} transformation.

A size-n sample of blood pressure average errors (average of
the difference between three \emph{device} measurements
and the three corresponding reference readings) 
is drawn from a normal distribution 
with mean \code{mu} and standard deviation \code{std.dev}.
The r.v. \code{sigmaami} results from transforming \eqn{xbar}, the 
sample mean, assuming the proportion of drawn errors lie in 
(-\code{delta}, \code{delta})
with probability \code{ptolerror}, that is 
\code{sigmaami}\eqn{=sigmaami(xbar)}.

%A \emph{tolerable error}, according the AAMI-SP10 is and of 10mmHg or
%less on a single person, using the average of that person's readings.


}
\value{
\code{psigmaami} gives the distribution function.
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
%  \code{\link[AAMIsp1]{psigmaami}}
}
\examples{

%The probability of tolerable error with mu = 0, std.dev= 6, as per AAMI-SP10
%ProbTolError(mu =0, std.dev = 6, delta = 10)

## Sample of n = 85 average errors from a normal distribution with mean 3 and st.dev = 2.
psigmaami(sigmaami = 4, mu = 3, std.dev = 2, n = 85)



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
