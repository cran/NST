\name{dist.3col}
\alias{dist.3col}
\title{
Transform distance matrix to 3-column matrix
}
\description{
Transform a distance matrix to a 3-column matrix in which the first 2 columns indicate the pairwised samples/species names.
}
\usage{
dist.3col(dist)
}
\arguments{
  \item{dist}{a square matrix or distance object with column names and row names.}
}
\details{
In many cases, a 3-column matrix is easier to use than a distance matrix.
}
\value{
  \item{name1}{1st column, the first item of pairwised two items}
  \item{name2}{2nd column, the second item of pairwised two items}
  \item{dis}{3rd column, distance value of the pairwised two itmes}
}
\author{
Daliang Ning
}
\note{
Version 1: 2015.5.17
}
\examples{
data(tda)
comm=tda$comm
bray=beta.g(comm,dist.method="bray")
bray.3col=dist.3col(bray)
}

\keyword{ Tool }
