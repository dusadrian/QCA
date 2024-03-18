\name{generate}
\alias{generate}

\title{Generate a custom data structure}

\description{
This function acts as a DGS - Data Generating Structure for a certain SOP expression.
}

\usage{
generate(expression = "", snames = "", noflevels, ...)
}

\arguments{
  \item{expression}{String: a SOP - sum of products expression.}
  \item{snames}{A string containing the sets' names, separated by commas.}
  \item{noflevels}{Numerical vector containing the number of levels for each set.}
  \item{...}{Other arguments, mainly for internal use.}
}

\details{
Using the power of SOP expressions, this function can generate the data for any type
of expressions, either Boolean or multi-value.

Causal conditions should always be separated by a product sign "*", unless:
- they are single letters, or
- the set names are provided, or
- the expression is multi-value

All conditions are considered binary crisp, unless the number of levels are provided
in conjunction with the set names, in the order of their specification from the
\bold{\code{snames}} argument.

This is an extension of the function \bold{\code{expand()}} from package admisc,
the process of data generating process being essentially a Quine expansion to a
Disjunctive Normal Form.
}

\value{
A data frame.
}

\author{
Adrian Dusa
}

\examples{
generate(D + ~AB + B~C -> Z)

# the positive configurations in their complete DNF expansion:
expanded <- expand(D + ~AB + B~C, snames = c(A, B, C, D))
# ~A~B~CD + ~A~BCD + ~AB~CD + ~ABCD + A~B~CD + A~BCD +
# AB~CD + ABCD + ~AB~C~D + ~ABC~D + AB~C~D

# which has the equivalent simpler, initial expression:
simplify(expanded)
# D + ~AB + B~C

# same structure with different set names
# (note the mandatory use of the product sign *)
generate(Alpha + ~Beta*Gamma + Gamma*~Delta -> Omicron)

# introducing an additional, irrelevant condition
# (note the product sign is not mandatory if providing the set names)
setnames <- "Alpha, Beta, Gamma, Delta, Epsilon"
dat <- generate(Alpha + ~BetaGamma + Gamma~Delta -> Omicron, snames = setnames)

head(dat)

#   Alpha Beta Gamma Delta Epsilon Omicron
# 1     0    0     0     0       0       0
# 2     0    0     0     0       1       0
# 3     0    0     0     1       0       0
# 4     0    0     0     1       1       0
# 5     0    0     1     0       0       1
# 6     0    0     1     0       1       1

minimize(dat, outcome = Omicron)

# M1: Alpha + ~Beta*Gamma + Gamma*~Delta <-> Omicron
}

\keyword{functions}
