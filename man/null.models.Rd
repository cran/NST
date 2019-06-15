\name{null.models}
\alias{null.models}
\docType{data}
\title{
Options of null model algorithms
}
\description{
The parameters passing to function \code{\link{taxo.null}} for each null model algorithm
}
\usage{data("null.models")}
\format{
  A data frame with 13 rows on the following 3 variables. Rownames are null model algorithm IDs.
  \describe{
    \item{\code{sp.freq}}{character, how the species occurrence frequency will be constrainted in the null model.}
    \item{\code{samp.rich}}{character, how the species richness in each sample will be constrainted in the null model.}
    \item{\code{swap.method}}{character, method for fixed sp.freq and fixed samp.rich.}
  }
}
\references{
Gotelli NJ. Null model analysis of species co-occurrence patterns. Ecology 81, 2606-2621 (2000) doi:10.1890/0012-9658(2000)081[2606:nmaosc]2.0.co;2.
}
\examples{
data(null.models)
}
\keyword{null model}
\keyword{datasets}
