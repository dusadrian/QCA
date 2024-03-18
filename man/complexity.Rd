\name{complexity}

\alias{complexity}

\title{Number of combinations at a given complexity layer}

\description{
This function calculates the number of all possible combinations of conditions
(including all levels for each condition), at a given complexity layer.
}

\usage{
complexity(n, layers, noflevels, ...)
}

\arguments{
    \item{n}{Numeric scalar, the number of input conditions.}
    \item{layers}{Numeric vector, the complexity layer(s) with values from 1 to n.}
    \item{noflevels}{Numeric vector containing the number of levels for each
            of the n conditions.}
    \item{...}{Other arguments, mainly for internal use.}
}

\details{
These are the number of combinations which the CCubes algorithm (Dusa, 2018) checks
to determine the prime implicants from a minimization process.

In the bottom-up approach, CCubes first checks for single conditions (combinations of 
both presence and absence, or more levels if multi-value), then all possible combinations
of levels for two conditions etc.

The precise equation that partitions the search space into complexity layers is:

\deqn{
    \sum_{c = 1}^{k} {k\choose c} \prod_{s = 1}^{c} l_s
}

where \eqn{l} stands for the number of levels for each combination of \eqn{c}
conditions out of \eqn{k}.
}

\value{
A numeric vector.
}

\author{
Adrian Dusa
}


\references{
Dusa, A. (2018) \dQuote{Consistency Cubes: A Fast, Efficient Method for Boolean
Minimization}, R Journal, \doi{10.32614/RJ-2018-080}
}

\examples{
complexity(3) # all layers from 1 to 3

complexity(5, layers = 2)
}

\keyword{functions}
