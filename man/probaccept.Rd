%\name{f}
\name{ProbAccept}
\alias{ProbAccept}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of device acceptance - AANSI/AAMI-SP10 standard.
}
\description{
Probability of device acceptance - AANSI/AAMI-SP10 standard.
}
\usage{

  ProbAccept(n, mu, sd, ptolerror = 0.85,
             distribution = "normal",
             criteria = "SP10:2006",
             simulate = FALSE, sim.count = 1e4,
             noshow = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{

Sample size (evaluated people, as defined by Chandel et al. (2022)).
 
}
\item{mu, sd}{
mean and standar deviation of the average errors distribution (normal).
}
\item{ptolerror}{
Probability of tolerable error. Default is 0.85
}
\item{distribution}{
Underlying distribution the errors \eqn{\delta} are pulled from.
The default is \code{normal}, i.e. a normal distribution.
}
\item{criteria}{
SP10 criteria used.
}
\item{simulate}{
Logical. If \code{TRUE}, the acceptance probability is simulated from
\code{sim.count} samples form a normal distribution. Else,
the exact probability is returned.
Default is \code{FALSE}.
}
\item{sim.count}{
Integer, positive. Number of samples taken from normal distribution to
estimate the probability of accepting the device.
}
\item{noshow}{
Logical. If \code{FALSE} then results are prompted on terminal.
}
}
\details{

Computes the \emph{acceptance probability} of a device D for
blood pressure measuring under the ANSI/AAMI-SP10 standards
for a size-\code{n} sample of average errors from a normal
distribution with mean \code{mu} and stadard deviation \code{sd}.
The probability of tolerable error is set to 0.85, by default.
A \emph{tolerable error}, according to the ANSI/AAMI-SP10,
is a an error of 10mmHg or
less on a single person, using the average of that person's readings.

Fuller details in Chandel, et al. (2022).
The paper outlines the mathematical and
statistical aspects behind \code{ProbAccept}. Two random variables
are involved: the sample standard deviation and a transformation 
of the sample mean, resulting in a double integral
over a two-dimensional region.


}
\value{
It returns the probability of accepting the device
(either simulated or exact).
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


## Probability of accepting a device with bias (mu) = 5, and true standard
## deviation = 5.
ProbAccept(n = 85, mu = 5, sd = 5, ptolerror = 0.85)


## Probability of accepting a device with bias (mu) = 5, and true standard
## deviation = 7.
ProbAccept(n = 85, mu = 5, sd = 7, ptolerror = 0.85)



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
