\name{findTh}

\alias{findTh}

\title{Find calibration thresholds}

\description{
The purpose of this function is to automatically find calibration thresholds
for a numerical causal condition, to be split into separate groups.
}

\usage{
findTh(x, n = 1, hclustm = "complete", distm = "euclidean", ...)
}

\arguments{
  \item{x}{A numerical causal condition.}
  \item{n}{The number of thresholds to find.}
  \item{hclustm}{The agglomeration (clustering) method to be used.}
  \item{distm}{The distance measure to be used.}
  \item{...}{Other arguments (mainly for backwards compatibility).}
}

\details{
The process of calibration into crisp sets assumes expert knowledge about the
best threshold(s) that separate the raw data into the most meaningful
groups.

In the absence of such knowledge, an automatic procedure might help grouping
the raw data according to statistical clustering techniques.

The number of groups to split depends on the number of thresholds: one
thresholds splits into two groups, two thresholds splits into three groups etc.

For more details about how many groups can be formed with how many thresholds,
see \bold{\code{\link[stats]{cutree}()}}.

More details about the clustering techniques used in this function are found
using \bold{\code{\link[stats]{hclust}()}}, and also more details about different
distance measures can be found with \bold{\code{\link[stats]{dist}()}}. This
function uses their default values.

}

\value{
A numeric vector of length \bold{\code{n}}.
}

\author{
Adrian Dusa
}


\seealso{\code{\link[stats]{cutree}}, \code{\link[stats]{hclust}}, \code{\link[stats]{dist}}}

\examples{

# hypothetical list of country GDPs
gdp <- c(460, 500, 900, 2000, 2100, 2400, 15000, 16000, 20000)


# find one threshold to separate into two groups
findTh(gdp)
# 8700


# find two thresholds to separate into two groups
findTh(gdp, n = 2)
# 8700 18000


# using different clustering methods
findTh(gdp, n = 2, hclustm = "ward.D2", distm = "canberra")
# 1450 8700

}

\keyword{functions}

