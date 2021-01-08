\name{pNST}
\alias{pNST}
\title{
Normalized Stochasticity Ratio based on phylogenetic beta diversity
}
\description{
Calculate normalized stochasticity ratio according to method improved from Zhou et al (2014, PNAS), based on phylogenetic beta diversity index.
}
\usage{
pNST(comm, tree=NULL, pd=NULL,pd.desc=NULL,pd.wd=NULL,pd.spname=NULL,
     group, meta.group=NULL, abundance.weighted=TRUE, rand=1000,
     output.rand=FALSE, taxo.null.model=NULL, phylo.shuffle=TRUE,
     exclude.conspecifics=FALSE, nworker=4, LB=FALSE,
     between.group=FALSE, SES=FALSE, RC=FALSE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, rows are samples/sites, colnames are taxa (species/OTUs/ASVs)}
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{pd}{matrix, phylogenetic distance matrix.}
  \item{pd.desc}{character, the name of the file to hold the backingfile description of the phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function. If it is NULL and 'pd' is not given either, the fucntion pd.big will be used to calculate the phylogenetic distance matrix from tree, and save it in pd.wd as a big.memory file..}
  \item{pd.wd}{folder path, where the bigmemmory file of the phylogenetic distance matrix are saved.}
  \item{pd.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances.}
  \item{group}{a n x 1 matrix indicating the group or treatment of each sample, rownames are sample names. if input a n x m matrix, only the first column is used.}
  \item{meta.group}{a n x 1 matrix, to specify the metacommunity ID that each sample belongs to. NULL means the samples are from the same metacommunity.}
  \item{abundance.weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.rand}{Logic, whether to output dissimilarity results of each randomization. Default is FALSE.}
  \item{taxo.null.model}{Character, indicates null model algorithm to randomize the community matrix 'comm', including "EE", "EP", "EF", "PE", "PP", "PF", "FE", "FP", "FF", etc. The first letter indicate how to constraint species occurrence frequency, the second letter indicate how to constraint richness in each sample. see \code{\link{null.models}} for details. default is NULL, means not to randomze community but just randomize the tips, i.e. phylogeny shuffle, also named taxa shuffle.}
  \item{phylo.shuffle}{Logic, if TRUE, the null model algorithm "taxa shuffle" (Kembel 2009) is used, i.e. shuffling taxa labels across the tips of the phylogenetic tree to randomize phylogenetic relationships among species.}
  \item{exclude.conspecifics}{Logic, should conspecific taxa in different communities be exclude from MNTD calculations? default is FALSE.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{LB}{logic, whether to use a load balancing version of parallel computing code.}
  \item{between.group}{Logic, whether to calculate stochasticity for between-group turnovers. default is FALSE.}
  \item{SES}{Logic, whether to calculate standardized effect size, which is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.}
  \item{RC}{Logic, whether to calculate modified Raup-Crick metric, which is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.}
}
\details{
NST is a metric to estimate ecological stochasticity based on null model analysis of dissimilarity (Ning et al 2019). NST is improved from previous index ST (Zhou et al 2014). Modified stochasticity ratio (MST) is also calculated (Liang et al 2020; Guo et al 2018), which can be regarded as a spcial transformation of NST under assumption that observed similarity can be equal to mean of null similarity under pure stochastic assembly.

pNST is NST based on phylogenetic beta diversity (Ning et al 2019, Guo et al 2018), here, beta mean nearest taxon distance (bMNTD). pNST showed better performance in stochasticity estimation than tNST in some cases (Ning et al 2020).
}
\value{
Output is a list. Please DO NOT use NST.ij values in index.pair.grp and index.between.grp which can be out of [0,1] without ecologcial meanning. Please use \code{\link{nst.boot}} to get variation of NST.

  \item{index.pair}{indexes for each pairwise comparison. D.ij, observed dissimilarity, not standardized; G.ij, average null expectation of dissimilarity, not standardized; Ds.ij, observed dissimilarity, standardized to range from 0 to 1; Gs.ij,  average null expectation of dissimilarity, standardized; C.ij and E.ij are similarity and average null expectation of simmilarity, standardized if the dissimilarity has no fixed upper limit; ST.ij, stochasticity ratio calculated by previous method (Zhou et al 2014); MST.ij, modified stochasticity ratio calculated by a modified method (Liang et al 2020; Guo et al 2018); bNTI, beta nearest taxon index, i.e. standard effect size of difference between observed and null betaMNTD (Webb et al 2008); RC.bMNTD, modified Roup-Crick metrics (Chase et al 2011) but based on betaMNTD.}
  \item{index.grp}{mean value of each index in each group. group, group name; size, number of pairwise comparisons in this group; ST.i, group mean of stochasticity ratio, not normalized; NST.i, group mean of normalized stochasticity ratio; MST.i, group mean of modified stochasticity ratio.}
  \item{index.pair.grp}{pairwise values of each index in each group. group, group name; C.ij, E.ij, ST.ij, and MST.ij have the same meaning as in index.pair; NST.ij, the pairwise values of NST, for reference only, DO NOT use. Since NST is normalized ST calculated from ST.ij, NST pairwise values NST.ij have no ecological meaning. Variation of NST from bootstrapping test is preferred, see \code{\link{nst.boot}}.}
  \item{index.between}{mean value of each index between each two groups. Similar to index.grp, but calcualted from comparisons between each two groups.}
  \item{index.pair.between}{pairwise values of each index between each two groups. Similar to index.pair.grp, but calcualted from comparisons between each two groups.}
  \item{Dmax}{The maximum or upper limit of dissimilarity before standardized, which is used to standardize the dissimilarity with upper limit not equal to one.}
  \item{dist.method}{dissimilarity index name.}
  \item{details}{detailed results. rand.mean, mean of null dissimilarity for each pairwise comparison, not standardized; Dmax, the maximum or upper limit of dissimilarity before standardized; obs3, observed dissimilarity, not standardized; dist.ran, alll null dissimilarity values, each row is a pairwise comparison, each column is results from one randomization; group, input group informaiton; meta.group, input metacommunity information.}
}
\references{
Ning D., Deng Y., Tiedje J.M. & Zhou J. (2019) A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892-16898. doi:10.1073/pnas.1904623116.

Zhou J, Deng Y, Zhang P, Xue K, Liang Y, Van Nostrand JD, Yang Y, He Z, Wu L, Stahl DA, Hazen TC, Tiedje JM, and Arkin AP. (2014) Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America 111, E836-E845. doi:10.1073/pnas.1324044111.

Liang Y, Ning D, Lu Z, Zhang N, Hale L, Wu L, Clark IM, McGrath SP, Storkey J, Hirsch PR, Sun B, and Zhou J. (2020) Century long fertilization reduces stochasticity controlling grassland microbial community succession. Soil Biology and Biochemistry 151, 108023. doi:10.1016/j.soilbio.2020.108023.

Guo X, Feng J, Shi Z, Zhou X, Yuan M, Tao X, Hale L, Yuan T, Wang J, Qin Y, Zhou A, Fu Y, Wu L, He Z, Van Nostrand JD, Ning D, Liu X, Luo Y, Tiedje JM, Yang Y, and Zhou J. (2018) Climate warming leads to divergent succession of grassland microbial communities. Nature Climate Change 8, 813-818. doi:10.1038/s41558-018-0254-2.

Webb, C.O., Ackerly, D.D. & Kembel, S.W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics, 24, 2098-2100.

Chase JM, Kraft NJB, Smith KG, Vellend M, and Inouye BD. (2011) Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere 2, art24. doi:10.1890/es10-00117.1.

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 3: 2020.9.9, remove setwd; change dontrun to donttest and revise save.wd in help doc.
Version 2: 2020.8.22, add to NST package, update help document.
Version 1: 2018.1.9
}
\seealso{
\code{\link{tNST}}, \code{\link{nst.boot}}, \code{\link{nst.panova}}
}
\examples{
data("tda")
comm=tda$comm
group=tda$group
tree=tda$tree

# since it needs to save some file to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  save.wd=tempdir() # please change to the folder you want to use.
  nworker=2 # parallel computing thread number
  rand.time=20 # usually use 1000 for real data.
  pnst=pNST(comm=comm, tree=tree, group=group,
            pd.wd=save.wd, rand=rand.time, nworker=nworker)
}
}
\keyword{ Indexes }
