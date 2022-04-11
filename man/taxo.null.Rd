\name{taxo.null}
\alias{taxo.null}
\title{
Null models of taxonomic beta diversity
}
\description{
to randomize the taxonomic structures based on one of various null model algorithms.
}
\usage{
taxo.null(comm,sp.freq=c("not","equip","prop","prop.ab","fix"),
          samp.rich=c("not","equip","prop","fix"),
          swap.method=c("not","swap","tswap","quasiswap",
                        "backtrack"),burnin=0,
          abundance=c("not","shuffle","local","region"),
          region.meta=NULL,region.freq=NULL,dirichlet=FALSE)
}
\arguments{
  \item{comm}{matrix, community data, rownames are sample/site names, colnames are species names}
  \item{sp.freq}{character, the constraint of species occurrence frequency when randomizing taxonomic structures, see details.}
  \item{samp.rich}{character, the constraint of sample richness when randomizing taxonomic structures, see details.}
  \item{swap.method}{character, the swap method for fixed sp.freq and fixed samp.rich, see \code{\link[vegan]{commsim}} for details.}
  \item{burnin}{Nonnegative integer, specifying the number of steps discarded before starting simulation. Active only for sequential null model algorithms. Ignored for non-sequential null model algorithms. also see \code{\link[vegan]{nullmodel}}.}
  \item{abundance}{character, the method to draw individuals (abundance) into present species when randomizing taxonomic structures, see details.}
  \item{region.meta}{a numeric vector, to define the (relative) abundance of each species in metacommunity/regional pool. The names should be species IDs. If no name, it should be in exact the same order as columns of comm. Default is NULL, the relative abundance in metacommunity will be calculated from comm.}
  \item{region.freq}{a numeric vector, to define the occurrence frequency of each species in metacommunity/regional pool. The names should be species IDs. If no name, it should be in exact the same order as columns of comm. Default is NULL, the occurrence frequency in metacommunity will be calculated from comm. If sp.freq='fix', the input region.freq must be integers. If sp.freq='fix' and samp.rich='fix', since no applicable algorithm now, region.freq will be ignored.}
  \item{dirichlet}{Logic. If TRUE, the taxonomic null model will use Dirichlet distribution to generate relative abundances in randomized community matrix. If the input community matrix has all row sums no more than 1, the function will automatically set dirichlet=TRUE. default is FALSE.}
}
\details{
This function returns a randomized community dataset (one time randomization), used by the function \code{\link{tNST}}. The null models differentiated by how to deal with species occurrence frequency (sp.freq), species richness in each sample (samp.rich), relative abundances (abundance), and which swap method used if both sp.freq and samp.rich are fixed.

Options of sp.freq and samp.rich (Gotelli 2000):
not: the whole co-occurrence pattern (present/absent) is not randomized;
equip: all the species or samples have equal probability when randomizing;
prop: randomization according to probability proportional to observed species occurrence frequency or sample richness;
prop.ab: randomization according to probability proportional to observed regional abundance sum of each species, only for sp.freq;
fix: randomization maintains the species  occurrence frequency or sample richness exactly the same as observed.

Options of abundance:
not: not abundance weighted;
shuffle: randomly assign observed abundance values of observed species in a sample to species in this sample after the present/absent pattern has been randomized, thus shuffle can only be used if the richness is fixed. Similar to "richness" algorithm in R package picante (Kembel et al 2010);
local: randomly draw individuals into randomized species in a sample on the probablities proportional to observed species-abundance-rank curve in this sample. If randomized species number in this sample is more than observed, the probabilities of exceeding species will be proportional to observed minimum abundance. If randomized species number (rN) in this sample is less than observed, the probabilities will be proportional to the observed abundances of top rN observed species. The rank of randomized species in a sample is randomly assigned.
region: randomly draw individuals into each ranodmized species in each sample on the probabilities proportional to observed relative abundances of each species in the whole region, as described previously (Stegen et al 2013).
}
\value{
a matrix of community data, e.g. an randomized OTU table, is returned. Rownames are sample/site names, and colnames are species names.
}
\references{
Gotelli NJ. Null model analysis of species co-occurrence patterns. Ecology 81, 2606-2621 (2000) doi:10.1890/0012-9658(2000)081[2606:nmaosc]2.0.co;2.

Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD, Blomberg SP, and Webb CO. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26, 1463-1464 (2010) doi:10.1093/bioinformatics/btq166.

Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, Rockhold ML, and Konopka A. Quantifying community assembly processes and identifying features that impose them. Isme Journal 7, 2069-2079 (2013).

Others cited in \code{\link[vegan]{commsim}}.
}
\author{
Daliang Ning
}
\note{
Version 3: 2021.5.11, add option region.freq to specify occurrence frequency in regional pool.
Version 2: 2021.4.16, add option dirichlet to handle community matrix with relative abundance values rather than counts.
Version 1: 2015.10.22
}

\seealso{
\code{\link{tNST}}, \code{\link{ab.assign}}, \code{\link{null.models}}
}
\examples{
data(tda)
comm=tda$comm
comm.rand=taxo.null(comm,sp.freq="prop",samp.rich="fix",abundance="region")
}

\keyword{ Null model }
