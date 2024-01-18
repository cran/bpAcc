%\name{f}
\name{PAccept}
\alias{PAccept}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of acceptance - ANSI/AAMI-SP10 standard.
}
\description{
Probability of acceptance - ANSI/AAMI-SP10 standard.
}
\usage{

  PAccept(xbar,sd,N,distribution = "normal", criteria = "SP10:2006")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{xbar, sd}{
mean and standar deviation of the average errors distribution (normal).
}
\item{N}{

Sample size (number of  participants).
 
}
\item{distribution}{
Underlying distribution the errors \eqn{\delta} are pulled from.
The default is \code{normal}, i.e. a normal distribution.
}
\item{criteria}{
SP10 criteria used.
}


}
\details{

Computes the \emph{acceptance probability} of a device D for
blood pressure measuring under the ANSI/AAMI-SP10 standards
for a size-\code{n} sample of average errors from an asymptotically normal
distribution with mean \code{xbar} and stadard deviation \code{sd}.

\deqn{1-\frac{1}{2}[1+erf(\frac{0.78-\hat{p}}{\sqrt(2)(\frac{\hat{p}(1-\hat{p})}{N})})]}

The distribution is of the true probability of tolerable error \code{p} where the 
tolerable error according to the ANSI/AAMI-SP10, is an error between -10 mmHg 
to 10 mmHg on a single person, using average of that person's readings. Using 
the sampling distribution of sampling proportion, the probabilty of 
\code{p>=0.78} is evaluated, which is called as the probabilty of acceptance or 
probability that for a given sample size \code{n}, sample mean \code{xbar} and 
sample standard deviation \code{sd}, the device is meeting the SP10 criteria.

Complete details in Chandel, et al. (2023).
The paper outlines the mathematical and
statistical aspects behind \code{PAccept}. The threshold probability for acceptance according to ANSI/AAMI-SP10 is 95\% i.e.,\code{Prob(p>=0.78) >= 0.95}


}
\value{
It returns the probability of a device meeting the SP10 criteria

}
\references{
Chandel, T. and Lee, TC. and Lowe, A. and Miranda, V. (2023)
Blood Pressure Device Accuracy Evaluation: Statistical
Considerations with an Implementarion in R,
\emph{Statistics in Medicine} (under review).
}
\author{
Tanvi Chandel, Tet-Chuan Lee, Andrew Lowe, Victor Miranda.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{


## Probability of acceptance of a device for a sample size (N) = 30 with sample 
## mean (xbar) = 5, standard deviation = 5.
PAccept(xbar=5,sd=5,N=30,distribution = "normal", criteria = "SP10:2006")


## Probability of acceptance of a device for a sample size (N) = 60 with sample 
## mean (xbar) = 2, standard deviation = 7.
PAccept(xbar=2,sd=7,N=60,distribution = "normal", criteria = "SP10:2006")



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
