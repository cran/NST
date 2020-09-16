\name{tda}
\alias{tda}
\docType{data}
\title{
Test dataset A
}
\description{
A simple test data with a community matrix and treatment information
}
\usage{data("tda")}
\format{
  A list object with 2 elements.
  \describe{
    \item{\code{comm}}{matrix, community table; each row is a sample, thus rownames are sample IDs; each column is a taxon, thus colnames are OTU IDs.}
    \item{\code{group}}{matrix with only one column. treatment information; rownames are sample IDs; the only column shows treatment IDs.}
    \item{\code{tree}}{phylogenetic tree.}
  }
}
\examples{
data(tda)
comm=tda$comm
group=tda$group
}

\keyword{ Dataset }
