\name{nst.panova}
\alias{nst.panova}
\title{
Permutational multivariate ANOVA test for ST and NST
}
\description{
Permutational multivariate ANOVA test for stochasticity ratio and normalized stochasticity ratio between treatments
}
\usage{
nst.panova(nst.result, group=NULL, rand=999, trace=TRUE, SES=FALSE, RC=FALSE)
}
\arguments{
  \item{nst.result}{list object, the output of nsto, must have "details"}
  \item{group}{nx1 matrix, if the grouping is different from the nst.result. default is NULL, means to use the grouping in nst.result.}
  \item{rand}{integer, randomization times for permuational test}
  \item{trace}{logic, whether to show message when randomizing.}
  \item{SES}{Logic, whether to perform the test for standardized effect size (SES). SES is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.}
  \item{RC}{Logic, whether to perform the test for modified Raup-Crick metric (RC). RC is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.}
}
\details{
PERMANOVA for stochasticity ratio (ST or NST or MST) or SES or RC is based on the comparison of F values between observed pattern and the permutated patterns where samples are randomly shuffled regardless of treatments.  However, it is a bit different from PERMANOVA for dissimilarity. The PERMANOVA of stochasticity ratio here is to ask whether the ST values within a group is higher than those within another group. But the PERMANOVA of dissimilarity is to ask whether the between-group dissimilarity is higher than within-group dissimilarity.
}
\value{
Output is a data.frame object.
  \item{index}{name of index}
  \item{group1}{treatment/group name}
  \item{group2}{treatment/group name}
  \item{Index.group1}{index value in group1}
  \item{Index.group2}{index value in group2}
  \item{Difference}{index.group1 - index.group2}
  \item{F.obs}{F value}
  \item{P.anova}{P value of parametric ANOVA test}
  \item{P.panova}{P value of permutational ANOVA test}
  \item{P.perm}{P value of permutational test of the difference}
}
\references{
Ning D., Deng Y., Tiedje J.M. & Zhou J. (2019) A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892-16898. doi:10.1073/pnas.1904623116. 
}
\author{
Daliang Ning
}
\note{
Version 7: 2021.10.29, add PERMANOVA test for SES and RC.
Version 6: 2021.9.28, avoid error for special cases in permutation.
Version 5: 2021.8.25, revised to avoid error for special cases in MST calculation.
Version 4: 2020.10.14, debug some error when replecate number is low and edit details in help.
Version 3: 2019.10.8, Update reference.
Version 2: 2019.5.10
Version 1: 2017.12.30
}
\seealso{
\code{\link{tNST}}, \code{\link{nst.boot}}
}
\examples{
data(tda)
comm=tda$comm
group=tda$group
tnst=tNST(comm=comm, group=group, rand=20,
          output.rand=TRUE, nworker=1)
# rand is usually set as 1000, here set rand=20 to save test time.

nst.pova=nst.panova(nst.result=tnst, rand=99)
# rand is usually set as 999, here set rand=99 to save test time.
}

\keyword{ Significance }
