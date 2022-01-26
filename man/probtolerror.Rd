%\name{f}
\name{ProbTolError}
\alias{ProbTolError}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of tolerable error - AAMI-SP10 standard.
}
\description{
Probability of tolerable error - AAMI-SP10 standard.
}
\usage{

  ProbTolError(distribution = "normal", mu, std.dev, delta)
  

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{distribution}{
The errors distirbution. Default and only option is
"normal". 
}
\item{mu, std.dev}{

Mean and standar deviation of the average errors distribution (normal).

}
\item{delta}{
Maximum average error allowed for the device D.
}
}
\details{

Computes the \emph{probability of tolerable error} for a device D for
blood pressure measuring under the ANSI/AAMI-SP10 standard
for a size-\code{n} sample of average errors from a normal
distribution with mean \code{mu} and stadard deviation \code{sd}.
The maximum error accepted is 0.85,
A \emph{tolerable error}, according to the ANSI/AAMI-SP10,
is a an error of 10mmHg or
less on a single person, using the average of that person's readings.

Currently, only normally distirbuted errors are handled.
Further choices will be implemented over time.

Fuller details in Chandel, et al. (2022).



}
\value{
It returns the probability of tolerable error based on 
a normal distribution.
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


## Probability of tolerable error, mu = 4, sd = 5, delta = 10 (ANSI/AAMI-SP10)
ProbTolError(mu = 4, std.dev = 5, delta = 10)



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
