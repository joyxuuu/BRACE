\name{RuleIII}
\alias{RuleIII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BRACE RAR Rule III 
}
\description{
Calculate the number of patients allocated to each arm at the interim stage following BRACE RAR Rule III.
}
\usage{
RuleIII(nj, pvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nj}{
Total sample size at jth interim stage
}
  \item{pvec}{
Vector of allocation probability for all treatment arms
}
}
\value{
A vector of number of patients allocated to each arm.
}
\examples{
nj <- 200
pvec <- c(0.8,0.1,0.1)
n_assign <- RuleIII(nj=nj,pvec=pvec)	
}
