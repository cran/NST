\name{NST-package}
\alias{NST-package}
\alias{NST}
\docType{package}
\title{
Normalized Stochasticity Ratio
}
\description{
This package is to estimate ecological stochasticity in community assembly based on beta diversity. Various indexes can be calculated, including Stochasticity Ratio (ST), Normalized Stochasticity Ratio (NST), Modified Stochasticity Ratio (MST), Standard Effect Size (SES), and modified Raup-Crick metrics (RC), based on various taxonomic and phylogenetic dissimilarity metrics and different null model algorithms. All versions and examples are available from GitHub. URL: https://github.com/DaliangNing/NST

Version 2.0.4: Update citation and references. Emphasize that NST variation should be calculated from nst.boot rather than pairwise NST.ij from tNST. Emphasize that different group setting in tNST may lead to different NST results.
Version 3.0.1: Add NST based on phylogenetic beta diversity (pNST).
Version 3.0.2: debug pNST.
Version 3.0.3: remove setwd in functions; change dontrun to donttest and revise save.wd in help doc.
Version 3.0.4: update github link of NST; update nst.boot and nst.panova to include MST results.
Version 3.0.5: debug nst.panova.
Version 3.0.6: update references.
Version 3.1.1: add options to allow input propotional data (rather than counts) as community matrix, as well as community data transformation before dissimilarity calculation.
Version 3.1.2: provide temporary solution for the failure of makeCluster in some OS.
Version 3.1.3: add options to specify occurrence frequency in regional pool.
Version 3.1.4: debug ab.assign.
Version 3.1.5: add function cNST to calculate NST using user customized beta diversity and the null results.
Version 3.1.6: revise functions tNST, pNST, cNST, nst.boot, and nst.panova to avoid error for special cases in MST calculation.
Version 3.1.7(20210928): revise function nst.panova to avoid error for special cases in permutation.
Version 3.1.8(20211029): add summary and test for SES and RC in functions tNST, pNST, cNST, nst.boot, and nst.panova.
Version 3.1.9(20220410): address notes from package check.
Version 3.1.10(20220603): tested with the latest version of package iCAMP.
}
\details{
\tabular{ll}{
Package: \tab NST\cr
Type: \tab Package\cr
Version: \tab 3.1.10\cr
Date: \tab 2022-6-3\cr
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
