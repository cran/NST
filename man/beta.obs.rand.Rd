\name{beta.obs.rand}
\alias{beta.obs.rand}
\docType{data}
\title{
Test data B observed and null beta diversity
}
\description{
A simple dataset of observed and null beta diversity values, with sample grouping information.
}
\usage{data("beta.obs.rand")}
\format{
  A list object with 3 elements.
  \describe{
    \item{\code{obs}}{matrix, pairwise values of beta diversity (dissimilarity).}
    \item{\code{rand}}{list, each element shows the beta diversity of randomized communities from a null model algorithm.}
    \item{\code{group}}{data.frame, only one column showing which samples are controls and which are under treatment.}
  }
}
\examples{
data(beta.obs.rand)
beta.obs=beta.obs.rand$obs
beta.rand.list=beta.obs.rand$rand
group=beta.obs.rand$group
}
\keyword{datasets}
