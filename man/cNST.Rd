\name{cNST}
\alias{cNST}
\title{
Normalized Stochasticity Ratio based on customized metrics and null results
}
\description{
Calculate normalized stochasticity ratio (NST) based on given values of observed and null beta diveresity metrics.
}
\usage{
cNST(beta.obs, beta.rand.list, group,
     Dmax = 1, between.group = FALSE,
     SES = FALSE, RC = FALSE, output.detail = FALSE)
}

\arguments{
  \item{beta.obs}{square matrix or distance object, to provide the observed pairwise values of  beta diversity (dissimilarity).}
  \item{beta.rand.list}{a list object. Each element of the list is a square matrix or a distance object, to provide null values of beta diversity from a null model.}
  \item{group}{matrix or data.frame, a one-column (n x 1) matrix indicating the group or treatment of each sample, rownames are sample IDs. if input a n x m matrix, only the first column is used. Attention: different group setting will change NST values.}
  \item{Dmax}{The maximum or upper limit of dissimilarity before standardized, which is used to standardize the dissimilarity with upper limit not equal to one.}
  \item{between.group}{Logic, whether to calculate stochasticity for between-group turnovers. default is FALSE.}
  \item{SES}{Logic, whether to calculate standardized effect size, which is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.}
  \item{RC}{Logic, whether to calculate modified Raup-Crick metric, which is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.}
  \item{output.detail}{Logic, whether to output some details, including dissimilarity results of each randomization. Default is FALSE.}
}
\details{
NST is a metric to estimate ecological stochasticity based on null model analysis of dissimilarity. When using the function tNST or pNST, you can only select the metrics or null model from the given options. This function gives more flexibility if your beta diversity metric or null model algorithm is not included in tNST or pNST's options.
}
\value{
Output is a list. Please DO NOT use NST.ij values in index.pair.grp and index.between.grp which can be out of [0,1] without ecologcial meanning. Please use \code{\link{nst.boot}} to get variation of NST.

  \item{index.pair}{indexes for each pairwise comparison. D.ij, observed dissimilarity, not standardized; G.ij, average null expectation of dissimilarity, not standardized; Ds.ij, observed dissimilarity, standardized to range from 0 to 1; Gs.ij,  average null expectation of dissimilarity, standardized; C.ij and E.ij are similarity and average null expectation of simmilarity, standardized if the dissimilarity has no fixed upper limit; ST.ij, stochasticity ratio calculated by previous method (Zhou et al 2014); MST.ij, modified stochasticity ratio calculated by a modified method (Liang et al 2020; Guo et al 2018); SES.ij, standard effect size of difference between observed and null dissimilarity (Kraft et al 2011); RC.ij, modified Roup-Crick metrics (Chase et al 2011, Stegen et al 2013).}
  \item{index.grp}{mean value of each index in each group. group, group name; size, number of pairwise comparisons in this group; ST.i, group mean of stochasticity ratio, not normalized; NST.i, group mean of normalized stochasticity ratio; MST.i, group mean of modified stochasticity ratio; SES.i, group mean of standard effect size; RC.i, group mean of modified Roup-Crick metric.}
  \item{index.pair.grp}{pairwise values of each index in each group. group, group name; C.ij, E.ij, ST.ij, MST.ij, SES.ij, and RC.ij have the same meaning as in index.pair; NST.ij, the pairwise values of NST, for reference only, DO NOT use. Since NST is normalized ST calculated from ST.ij, NST pairwise values NST.ij have no ecological meaning. Variation of NST from bootstrapping test is preferred, see \code{\link{nst.boot}}.}
  \item{index.between}{mean value of each index between each two groups. Similar to index.grp, but calcualted from comparisons between each two groups.}
  \item{index.pair.between}{pairwise values of each index between each two groups. Similar to index.pair.grp, but calcualted from comparisons between each two groups.}
  \item{Dmax}{The maximum or upper limit of dissimilarity before standardized, which is used to standardize the dissimilarity with upper limit not equal to one. See \code{\link{beta.limit}} for details.}
  \item{details}{detailed results. rand.mean, mean of null dissimilarity for each pairwise comparison, not standardized; Dmax, the maximum or upper limit of dissimilarity before standardized; obs3, observed dissimilarity, not standardized; dist.ran, alll null dissimilarity values, each row is a pairwise comparison, each column is results from one randomization; group, input group informaiton; meta.group, input metacommunity information.}
}
\references{
Ning D., Deng Y., Tiedje J.M. & Zhou J. (2019) A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892-16898. doi:10.1073/pnas.1904623116.

Zhou J, Deng Y, Zhang P, Xue K, Liang Y, Van Nostrand JD, Yang Y, He Z, Wu L, Stahl DA, Hazen TC, Tiedje JM, and Arkin AP. (2014) Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America 111, E836-E845. doi:10.1073/pnas.1324044111.

Liang Y, Ning D, Lu Z, Zhang N, Hale L, Wu L, Clark IM, McGrath SP, Storkey J, Hirsch PR, Sun B, and Zhou J. (2020) Century long fertilization reduces stochasticity controlling grassland microbial community succession. Soil Biology and Biochemistry 151, 108023. doi:10.1016/j.soilbio.2020.108023.

Guo X, Feng J, Shi Z, Zhou X, Yuan M, Tao X, Hale L, Yuan T, Wang J, Qin Y, Zhou A, Fu Y, Wu L, He Z, Van Nostrand JD, Ning D, Liu X, Luo Y, Tiedje JM, Yang Y, and Zhou J. (2018) Climate warming leads to divergent succession of grassland microbial communities. Nature Climate Change 8, 813-818. doi:10.1038/s41558-018-0254-2.

Kraft NJB, Comita LS, Chase JM, Sanders NJ, Swenson NG, Crist TO, Stegen JC, Vellend M, Boyle B, Anderson MJ, Cornell HV, Davies KF, Freestone AL, Inouye BD, Harrison SP, and Myers JA. (2011) Disentangling the drivers of beta diversity along latitudinal and elevational gradients. Science 333, 1755-1758. doi:10.1126/science.1208584.

Chase JM, Kraft NJB, Smith KG, Vellend M, and Inouye BD. (2011) Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere 2, art24. doi:10.1890/es10-00117.1.

Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, Rockhold ML, and Konopka A. (2013) Quantifying community assembly processes and identifying features that impose them. The Isme Journal 7, 2069. doi:10.1038/ismej.2013.93.
}
\author{
Daliang Ning
}
\note{
Version 3: 2021.10.29, add summary of SES and RC.
Version 2: 2021.8.25, revised to avoid error for special cases in MST calculation.
Version 1: 2021.7.29
}

\seealso{
\code{\link{nst.boot}}, \code{\link{nst.panova}}, \code{\link{taxo.null}}, \code{\link{tNST}}, \code{\link{pNST}}
}
\examples{
data(beta.obs.rand)
beta.obs=beta.obs.rand$obs
beta.rand.list=beta.obs.rand$rand
group=beta.obs.rand$group
nst=cNST(beta.obs=beta.obs, beta.rand.list=beta.rand.list,
         group=group,Dmax = 1,between.group = TRUE,
         SES = TRUE, RC = TRUE, output.detail = FALSE)
}
\keyword{ Indexes }
