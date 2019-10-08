\name{NST-package}
\alias{NST-package}
\alias{NST}
\docType{package}
\title{
Normalized Stochasticity Ratio
}
\description{
This package is to estimate ecological stochasticity in community assembly based on beta diversity. Various indexes can be calculated, including Stochasticity Ratio (ST), Normalized Stochasticity Ratio (NST), Modified Stochasticity Ratio (MST), Standard Effect Size (SES), and modified Raup-Crick metrics (RC), based on various dissimilarity metrics and different null model algorithms.
Version 2.0.4: Update citation and references. Emphasize that NST variation should be calculated from nst.boot rather than pairwise NST.ij from tNST. Emphasize that different group setting in tNST may lead to different NST results.
}
\details{
\tabular{ll}{
Package: \tab NST\cr
Type: \tab Package\cr
Version: \tab 2.0.4\cr
Date: \tab 2019-10-08\cr
License: \tab GPL-2\cr
}
}

\author{
Daliang Ning <ningdaliang@ou.edu>
}
\references{
Ning D., Deng Y., Tiedje J.M. & Zhou J. (2019) A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences 116, 16892-16898. doi:10.1073/pnas.1904623116. 

Zhou J, Deng Y, Zhang P, Xue K, Liang Y, Van Nostrand JD, Yang Y, He Z, Wu L, Stahl DA, Hazen TC, Tiedje JM, and Arkin AP. (2014) Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America 111, E836-E845. doi:10.1073/pnas.1324044111.
}
\keyword{package}
\examples{
data(tda)
comm=tda$comm
group=tda$group
tnst=tNST(comm=comm, group=group, dist.method="jaccard",
          abundance.weighted=TRUE, rand=100,
          nworker=1, null.model="PF", between.group=TRUE,
          SES=TRUE, RC=TRUE)
}
