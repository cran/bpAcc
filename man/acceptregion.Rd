%\name{f}
\name{AcceptR}
\alias{AcceptR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Acceptance Region - ANSI/AAMI-SP10 standard.
}
\description{
Acceptance Region for a given sample size to comply with ANSI/AAMI-SP10 standard.
}
\usage{

    AcceptR(n,distribution = "normal",criteria = "SP10:2006")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{

Sample size .
 
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

Computes the revised probability of tolerable error i.e. the minimum probability of errors \eqn{\delta} within tolerable error for a given sample size-\code{n}. The \emph{tolerable error}, according to the ANSI/AAMI-SP10, is an error between -10 mmHg to 10 mmHg on a single person, using average of that person's readings. The revised probability of tolerable error varies for different sample sizes. Thus to meet the SP10 criteria, an acceptance region based on the sample size is provided by the package.

\deqn{\int_{-\infty}^{\delta} f_N(x;\mu_0,\sigma)-\int_{-\infty}^{-\delta} f_N(x;\mu_0,\sigma) \; dx = \hat{p}}
Complete details in Chandel, et al. (2023).
The paper outlines he mathematical and
statistical aspects behind \code{AcceptR}. 

}
\value{
It returns the acceptance region for a sample size i.e the upper limit of \code{sd} for a sample mean for a given \code{n}. This gives clinicians a flexible way to comprehend how the upper limit of the standard deviation -\code{sd} fluctuates depending on the sample size.
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

%\seealso{
%  \code{\link[bpAcc:ProbAccept]{ProbAccept}}
%  \code{\link[bpAcc:psigmaami]{psigmaami}}
%}

\examples{


## SP10-Acceptance Region for a sample size (n) = 85

AcceptR(n=85,distribution = "normal",criteria = "SP10:2006")


## SP10-Acceptance Region for a sample size (n) = 50 

AcceptR(n=50,distribution = "normal",criteria = "SP10:2006")

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
