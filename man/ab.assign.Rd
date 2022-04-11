\name{ab.assign}
\alias{ab.assign}
\title{
Randomly draw individuals into species according to specified probabilities
}
\description{
This funciton is to assign abundances to species when randomizing communities based on null models considering abundances. Individuals are randomly drawn into species according to the specified probabilities.
}
\usage{
ab.assign(comm.b, samp.ab=NULL, prob.ab)
}
\arguments{
  \item{comm.b}{numeric matrix, binary (present/absent) community data, rownames are sample/site names, colnames are species names.}
  \item{samp.ab}{numeric vector, total abundances (total individual numbers) in each sample. If samp.ab=NULL, Dirichlet distribution will be used to generate randomized community matrix with relative abundance (proportion) of each taxon in each sample.}
  \item{prob.ab}{numeric matrix, probability of each species into which the individuals in a certain sample are drawn.}
}
\details{
This function is called by the function \code{\link{taxo.null}} to generate randomized communities.
}
\value{
A matrix of community data with abundances (or relative abundance) is returned. rownames are sample/site names, and colnames are species names.
}
\references{
Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, Rockhold ML, and Konopka A. Quantifying community assembly processes and identifying features that impose them. Isme Journal 7, 2069-2079 (2013).
}
\author{
Daliang Ning
}
\note{
Version 3: 2021.7.27, debug, if samp.ab is lower than samp.rich, no need to assign abundance.
Version 2: 2021.4.16, add new algorithm based on Dirichlet distribution.
Version 1: 2015.10.22.
}

\seealso{
\code{\link{taxo.null}}
}
\examples{
data(tda)
comm=tda$comm
comm.b=comm
comm.b[comm.b>0]=1
samp.ab=rowSums(comm)
prob.ab=matrix(colSums(comm),nrow=nrow(comm),ncol=ncol(comm),byrow=TRUE)
comm.rand=ab.assign(comm.b,samp.ab,prob.ab)
}

\keyword{ Null model }
