\name{nst.boot}
\alias{nst.boot}
\title{
Bootstrapping test for ST and NST
}
\description{
To test the distribution of ST and NST in each group, and the significance of ST and NST difference between each pair of groups.
}
\usage{
nst.boot(nst.result, group=NULL, rand=999, trace=TRUE,
         two.tail=FALSE, out.detail=FALSE, between.group=FALSE,
         nworker=1)
}
\arguments{
  \item{nst.result}{list object, the output of tNST, must have "details", thus the output.rand must be TRUE in tNST function.}
  \item{group}{n x 1 matrix, if the grouping is different from the nst.result. default is NULL, means to use the grouping in nst.result.}
  \item{rand}{integer, random draw times for bootstrapping test.}
  \item{trace}{logic, whether to show message when randomizing.}
  \item{two.tail}{logic, the p value is two-tail or one-tail.}
  \item{out.detail}{logic, whether to output details rather than just summarized results.}
  \item{between.group}{Logic, whether to calculate for between-group turnovers. default is FALSE.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 1, means not to use parallel computing.}
}
\details{
Normalized stochasticity ratio (NST, Ning et al 2019) is a index to estimate average stochasticity within a group of samples. Bootstrapping is an excellent method to evaluate the statistical variation. Since the observed/null dissimilarity values are not independent (pairwise comparisons), bootstrapping should be random draw of samples rather than the pairwise values. Bootstrapping for stochasticity ratio (ST, Zhou et al 2014) is also performed.
}
\value{
Output is a list object, includes
  \item{summary}{Index, based on ST, NST, or MST; Group, group/treatment; obs, the index value of observed samples; mean, mean value of bootstrapping results; stdev, standard deviation; Min, minimal value; Quantile25, quantile of 25 percent; Median, median value; Quantile75, 75 percent quantile; Max, maximum value; LowerWhisker, LowerHinge, Median, HigherHinge, HigherWhisker, values for box-and-whisker plot; Outliers, outliers in bootstrapping values which out of the range of 1.5 fold of IQR.}
  \item{compare}{Comparison between each pair of groups, and p values. p.wtest, p value from wilcox.test; w.value, w value from wilcom.test; p.count, p value calculated by directly comparing all values in two groups; ..noOut, means outliers were not included for significance test. In principle, p.count or p.count.noOutis preferred, and others have defects.}
  \item{detail}{a list object. ST.boot, a list of bootstrapping detail results of ST for each group, each element in the list means the result of one random draw; NST.boot and MST.boot, bootstrapping results of NST and MST; ST.boot.rmout, bootstrapping results of ST without outliers; NST.boot.rmout and MST.boot.rmout, bootstrapping results of NST and MST without outliers; STb.boot, NSTb.boot, MSTb.boot, STb.boot.rmout, NSTb.boot.rmout, and MSTb.boot.rmout have the same meanning but for between-group comparisons.}
}
\references{
Ning D., Deng Y., Tiedje J.M. & Zhou J. (2019) A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892-16898. doi:10.1073/pnas.1904623116. 
}
\author{
Daliang Ning
}
\note{
Version 4: 2020.9.19, Add MST results into output.
Version 3: 2019.10.8, Update reference.
Version 2: 2019.5.10
Version 1: 2018.1.9
}
\seealso{
\code{\link{tNST}}, \code{\link{nst.panova}}
}
\examples{
data(tda)
comm=tda$comm
group=tda$group
tnst=tNST(comm=comm, group=group, rand=20,
          output.rand=TRUE, nworker=1)
# rand is usually set as 1000, here set rand=20 to save test time.

nst.bt=nst.boot(nst.result=tnst, group=NULL, rand=99,
                trace=TRUE, two.tail=FALSE, out.detail=FALSE,
                between.group=FALSE, nworker=1)
# rand is usually set as 999, here set rand=99 to save test time.

}

\keyword{ Significance }
