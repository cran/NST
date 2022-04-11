\name{beta.g}
\alias{beta.g}
\alias{chaosorensen}
\alias{chaojaccard}
\title{
Various taxonomic beta diversity indexes
}
\description{
This function can simultaneously calculate various taxonomic dissimilarity indexes, mainly based on vegdist from package vegan.
}
\usage{
beta.g(comm, dist.method="bray", abundance.weighted=TRUE,
       as.3col=FALSE,out.list=TRUE, transform.method=NULL, logbase=2)
chaosorensen(comm, dissimilarity=TRUE, to.dist=TRUE)
chaojaccard(comm, dissimilarity=TRUE, to.dist=TRUE)
}
\arguments{
  \item{comm}{Community data matrix. rownames are sample names. colnames are species names.}
  \item{dist.method}{A character or vector indicating one or more index(es). match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis", "mGower", "mEuclidean", "mManhattan", "chao.jaccard", "chao.sorensen". default is "bray"}
  \item{abundance.weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{as.3col}{Logic, output a 3-column matrix (TRUE) or a square matrix (FALSE) for each index. default is FALSE.}
  \item{out.list}{Logic, if using multiple indexes, output their results as a list (TRUE) or a matrix combining all 3-column matrixes (FALSE). if out.list=FALSE, as.3col will be forced to be TRUE. default is TRUE.}
  \item{dissimilarity}{Logic, calculate dissimilarity or similarity. default is TRUE, means to return dissimilarity.}
  \item{to.dist}{Logic, return distance object or squared matrix. default is TRUE, means to return distance object.}
  \item{transform.method}{character or a defined function, to specify how to transform community matrix before calculating dissimilarity. if it is a characher, it should be a method name as in the function 'decostand' in package 'vegan', including 'total','max','freq','normalize','range','standardize','pa','chi.square','cmdscale','hellinger','log'.}
  \item{logbase}{numeric, the logarithm base used when transform.method='log'.}
}
\details{
All the taxonomic beta diversity indexes are mainly calculated by vegdist in package vegan, except following methods:

mGower, mEuclidean, and mManhattan are modified from Gower, Euclidean, and Manhattan, respectively, according to the method reported previously (Anderson et al 2006).

chao.jaccard and chao.sorensen are calculated as described previously (Chao et al 2005), using open-source code from R package "fossil" (Vavrek 2011), but output as dissimilarity for each pairwise comparison.
}
\value{
beta.g will return a square matrix of each index if as.3col=FALSE, and combined as a list if out.list=TRUE (default). A 3-column matrix with first 2 columns indicating the pairwised samples will be output for each index if as.3col=TRUE, and combined as a list if out.list=TRUE  or integrated into one matrix if out.list=FALSE.

chaosorensen and chaojaccard will return a distance object (if to.dist=TRUE) or a squared matrix (if to.dist=FALSE).
}
\references{
Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2019). vegan: Community Ecology Package. R package version 2.5-4.

Anderson MJ, Ellingsen KE, & McArdle BH (2006) Multivariate dispersion as a measure of beta diversity. Ecol Lett 9(6):683-693.

Chao, A., R. L. Chazdon, et al. (2005) A new statistical approach for assessing similarity of speciescomposition with incidence and abundance data. Ecology Letters 8: 148-159

Vavrek, Matthew J. 2011. fossil: palaeoecological and palaeogeographical analysis tools. Palaeontologia Electronica, 14:1T.

Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful transformations for ordination of species data. Oecologia 129, 271--280.

Others cited in the help document of vegdist in R package vegan.
}
\author{
Daliang Ning
}
\note{
Version 3: 2021.4.16, add option to transform community matrix.
Version 2: 2019.5.10.
Version 1: 2015.9.25.
}
\seealso{
\code{\link{tNST}}
}
\examples{
data(tda)
comm=tda$comm
# calculate one index
beta.bray=beta.g(comm=comm,as.3col=TRUE)

# calculate multiple indexes
beta.td=beta.g(comm=comm,dist.method=c("bray","jaccard","euclidean",
              "manhattan","binomial","chao","cao"),
              abundance.weighted = TRUE,out.list=FALSE)
}

\keyword{ Indexes }
