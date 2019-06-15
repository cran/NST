\name{tNST}
\alias{tNST}
\title{
Taxonomic Normalized Stochasticity Ratio (tNST)
}
\description{
Calculate normalized stochasticity ratio (NST)  based on specified taxonomic dissimilarity index and null model algorithm.
}
\usage{
tNST(comm, group, meta.group=NULL, meta.com=NULL,
     dist.method="jaccard", abundance.weighted=TRUE,
     rand=1000, output.rand=FALSE, nworker=4,
     LB=FALSE, null.model="PF", between.group=FALSE,
     SES=FALSE, RC=FALSE)
}

\arguments{
  \item{comm}{matrix or data.frame, local community data, each row is a sample or site, each colname is a species or OTU or gene, thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{group}{matrix or data.frame, a one-column (n x 1) matrix indicating the group or treatment of each sample, rownames are sample IDs. if input a n x m matrix, only the first column is used.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity IDs. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL.}
  \item{meta.com}{matrix or data.frame, metacommunity data. Each row can be a samlple or a metacommunity, thus rownames are sample IDs or metacommunity IDs. Such that the relative abundance of each taxa in metacommunity can be different from the average relative abundance in the observed samples. This can be useful for uneven sampling design. NULL means the relative abundance of each taxa in the metacommunity can be directly calculated from the local community data (comm). Default is NULL.}
  \item{dist.method}{A character indicating dissimilarity index, including "manhattan","mManhattan", "euclidean","mEuclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "mGower", "morisita", "horn", "binomial", "chao", "cao". default is "jaccard"}
  \item{abundance.weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{rand}{integer, randomization times. default is 1000}
  \item{output.rand}{Logic, whether to output dissimilarity results of each randomization. Default is FALSE.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{LB}{Logic, whether to use a load balancing version for parallel computation. If TRUE, this can result in better cluster utilization, but increased communication can reduce performance. default is FALSE.}
  \item{null.model}{Character, indicates null model algorithm, including "EE", "EP", "EF", "PE", "PP", "PF", "FE", "FP", "FF", etc. The first letter indicate how to constraint species occurrence frequency, the second letter indicate how to constraint richness in each sample. see \code{\link{null.models}} for details. default is "PF".}
  \item{between.group}{Logic, whether to calculate stochasticity for between-group turnovers. default is FALSE.}
  \item{SES}{Logic, whether to calculate standardized effect size, which is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.}
  \item{RC}{Logic, whether to calculate modified Raup-Crick metric, which is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.}
}
\details{
NST is a metric to estimate ecological stochasticity based on null model analysis of dissimilarity. It is improved from previous index ST (Zhou et al 2014). Detailed description will be published soon (Ning et al 2019). Modified stochasticity ratio (MST) is also calculated (Liang et al 2019; Guo et al 2018), which can be regarded as a spcial transformation of NST under assumption that observed similarity can be equal to mean of null similarity under pure stochastic assembly.}
\value{
Output is a list.
  \item{index.pair}{indexes for each pairwise comparison. D.ij, observed dissimilarity, not standardized; G.ij, average null expectation of dissimilarity, not standardized; Ds.ij, observed dissimilarity, standardized to range from 0 to 1; Gs.ij,  average null expectation of dissimilarity, standardized; C.ij and E.ij are similarity and average null expectation of simmilarity, standardized if the dissimilarity has no fixed upper limit; ST.ij, stochasticity ratio calculated by previous method (Zhou et al 2014); MST.ij, modified stochasticity ratio calculated by a modified method (Liang et al 2019; Guo et al 2018); SES.ij, standard effect size of difference between observed and null dissimilarity (Kraft et al 2011); RC.ij, modified Roup-Crick metrics (Chase et al 2011, Stegen et al 2013).}
  \item{index.grp}{mean value of each index in each group. group, group name; size, number of pairwise comparisons in this group; ST.i, group mean of stochasticity ratio, not normalized; NST.i, group mean of normalized stochasticity ratio; MST.i, group mean of modified stochasticity ratio.}
  \item{index.pair.grp}{pairwise values of each index in each group. group, group name; C.ij, E.ij, ST.ij, and MST.ij have the same meaning as in index.pair; NST.ij, the pairwise values of NST. Since NST is normalized ST calculated from ST.ij, NST pairwise values NST.ij have no ecological meaning, just help to roughly estimate statistical variation of NST estimation. Variation of NST from bootstrapping test is preferred, see \code{\link{nst.boot}}.}
  \item{index.between}{mean value of each index between each two groups. Similar to index.grp, but calcualted from comparisons between each two groups.}
  \item{index.pair.between}{pairwise values of each index between each two groups. Similar to index.pair.grp, but calcualted from comparisons between each two groups.}
  \item{Dmax}{The maximum or upper limit of dissimilarity before standardized, which is used to standardize the dissimilarity with upper limit not equal to one. See \code{\link{beta.limit}} for details.}
  \item{dist.method}{dissimilarity index name.}
  \item{details}{detailed results. rand.mean, mean of null dissimilarity for each pairwise comparison, not standardized; Dmax, the maximum or upper limit of dissimilarity before standardized; obs3, observed dissimilarity, not standardized; dist.ran, alll null dissimilarity values, each row is a pairwise comparison, each column is results from one randomization; group, input group informaiton; meta.group, input metacommunity information.}
}
\references{
Ning D, Deng Y, Tiedjee JM, and Zhou J. (2019) A General Framework for Quantitatively Assessing Ecological Stochasticity. in review.

Zhou J, Deng Y, Zhang P, Xue K, Liang Y, Van Nostrand JD, Yang Y, He Z, Wu L, Stahl DA, Hazen TC, Tiedje JM, and Arkin AP. (2014) Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America 111, E836-E845. doi:10.1073/pnas.1324044111.

Liang Y, Ning D, Lu Z, Zhang N, Hale L, Wu L, Clark IM, McGrath SP, Storkey J, Hirsch PR, Sun B, and Zhou J. (2019) Century long fertilization reduces stochasticity controlling grassland microbial community succession. bioRxiv, doi:10.1101/638908.

Guo X, Feng J, Shi Z, Zhou X, Yuan M, Tao X, Hale L, Yuan T, Wang J, Qin Y, Zhou A, Fu Y, Wu L, He Z, Van Nostrand JD, Ning D, Liu X, Luo Y, Tiedje JM, Yang Y, and Zhou J. (2018) Climate warming leads to divergent succession of grassland microbial communities. Nature Climate Change 8, 813-818. doi:10.1038/s41558-018-0254-2.

Kraft NJB, Comita LS, Chase JM, Sanders NJ, Swenson NG, Crist TO, Stegen JC, Vellend M, Boyle B, Anderson MJ, Cornell HV, Davies KF, Freestone AL, Inouye BD, Harrison SP, and Myers JA. (2011) Disentangling the drivers of beta diversity along latitudinal and elevational gradients. Science 333, 1755-1758. doi:10.1126/science.1208584.

Chase JM, Kraft NJB, Smith KG, Vellend M, and Inouye BD. (2011) Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere 2, art24. doi:10.1890/es10-00117.1.

Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, Rockhold ML, and Konopka A. (2013) Quantifying community assembly processes and identifying features that impose them. The Isme Journal 7, 2069. doi:10.1038/ismej.2013.93.
}
\author{
Daliang Ning
}
\note{
version 1: 2019.5.10
}
\seealso{
\code{\link{nst.boot}}, \code{\link{nst.panova}}, \code{\link{taxo.null}}, \code{\link{beta.limit}}
}
\examples{
data(tda)
comm=tda$comm
group=tda$group
tnst=tNST(comm=comm, group=group, meta.group=NULL, meta.com=NULL,
          dist.method="jaccard", abundance.weighted=TRUE, rand=20,
          output.rand=FALSE, nworker=1, LB=FALSE, null.model="PF",
          between.group=TRUE, SES=TRUE, RC=TRUE)
# rand is usually set as 1000, here set rand=20 to save test time.
tnst.sum=tnst$NSTi
}

\keyword{ Indexes }
