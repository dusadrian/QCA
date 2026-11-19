\name{PI chart functions: makeChart, findmin, solveChart}

\alias{makeChart}
\alias{findmin}
\alias{solveChart}

\title{Create and solve a prime implicants chart}

\description{
\code{makeChart()} shows which prime implicants cover each positive truth table
configuration. \code{findmin()} gives the smallest number of prime implicants
needed to cover them all. \code{solveChart()} shows which combinations do so.
}

\usage{
makeChart(primes = "", configs = "", snames = "", mv = FALSE, collapse = "*", ...)

findmin(chart, type = c("hybrid", "ilp"), ...)

solveChart(chart, row.dom = FALSE, all.sol = FALSE, depth = NULL, max.comb = 0,
           first.min = FALSE, ...)
}

\arguments{
  \item{primes}{A string containing prime implicants, separated by commas, or a
        matrix of implicants.}
  \item{configs}{A string containing  causal configurations, separated by commas,
        or a matrix of causal configurations in the implicants space.}
  \item{snames}{A string containing the sets' names, separated by commas.}
  \item{mv}{Logical, row and column names in multi-value notation.}
  \item{collapse}{Scalar character, how to collapse different parts of the row or
        column names.}
  \item{chart}{An object of class \code{"QCA_pic"} or a logical matrix.}
  \item{type}{How \code{findmin()} finds the minimum. The default,
        \code{"hybrid"}, lets QCA choose how to solve the chart. Use \code{"ilp"}
        to send it directly to an integer programming solver. Both give the
        same minimum number of prime implicants.}
  \item{row.dom}{Set to \code{TRUE} to drop a prime implicant if another
        covers all the same configurations and at least one more. This can
        leave out some equally small solutions.}
  \item{all.sol}{Set to \code{TRUE} to include larger solutions in which
        every prime implicant is needed to cover the chart.}
  \item{depth}{Accepted for compatibility with older calls.
        Changing it has no effect in \code{solveChart()}.}
  \item{max.comb}{Numeric, to stop searching for solutions (see Details).}
  \item{first.min}{Logical, to return only the very first minimal solution (see
        Details).}
  \item{...}{Other arguments (mainly for backwards compatibility).}
}

\details{
Each row of a PI chart is a prime implicant and each column is a positive
truth table configuration. An \code{x} means that the prime implicant covers
that configuration. If a column contains just one \code{x}, the corresponding
prime implicant is essential.

When entering configurations as expressions, use \code{snames} if the set
names are unclear, or separate the names with \code{*}, as in \code{EF*~GH}.
When entering matrices, use the same condition order in both. In a prime
implicant, \code{0} means that the condition has been eliminated.

For example, in the first chart below, \code{a} is the only prime implicant
that covers \code{a~bc}, so every solution must include \code{a}. Either
\code{b} or \code{~c} covers the remaining configuration. Thus
\code{findmin(chart)} returns \code{2}, and \code{solveChart(chart)} returns
both solutions, \code{a + b} and \code{a + ~c}.

The default search returns all the smallest solutions. Use
\code{first.min = TRUE} when one is enough. Use \code{all.sol = TRUE} to
also find larger solutions that would leave a configuration uncovered if
any of their terms were removed. This uses the whole chart, even when
\code{row.dom = TRUE}, and can take longer on large charts.

When searching for the smallest solutions, \code{max.comb} limits how many
combinations to check, in billions. For example, \code{max.comb = 0.001}
sets a limit of one million combinations. If the limit is reached, some
solutions may be missing. The default of \code{0} sets no limit.

With \code{type = "ilp"}, QCA uses Gurobi when available or the
\code{lp_solve} solver included with the package.
}

\value{
For \bold{\code{makeChart}}: a logical matrix of class \code{"QCA_pic"}.

For \bold{\code{findmin}}: the minimum number of prime implicants needed
to cover the chart. With \code{solind = TRUE}, a value of \code{1} marks
each row used in one smallest solution and \code{0} marks the other rows.

For \bold{\code{solveChart}}: a matrix with one solution per column.
The numbers identify the rows of the PI chart used in that solution.
}

\author{
Adrian Dusa
}

\references{
Quine, W.V. (1952) \emph{The Problem of Simplifying Truth Functions},
The American Mathematical Monthly, vol.59, no.8. (Oct., 1952), pp.521-531.

Ragin, Charles C. (1987) \emph{The Comparative Method. Moving beyond qualitative and
quantitative strategies}, Berkeley: University of California Press
}

\examples{
# non-standard products, it needs the set names
chart <- makeChart("a, b, ~c", "abc, a~b~c, a~bc, ~ab~c")

# same with unquoted expressions
chart <- makeChart(c(a, b, ~c), c(abc, a~b~c, a~bc, ~ab~c))

chart
#      abc   a~b~c a~bc  ~ab~c
# a      x     x     x     -
# b      x     -     -     x
# ~c     -     x     -     x

findmin(chart)
# 2

solveChart(chart)
# first and second rows (a + b)
# and first and third rows (a + ~c)
# a is an essential prime implicant
#      a + b  a + ~c
#      [,1]   [,2]
# [1,]    1      1
# [2,]    2      3

# using SOP standard product sign
rows <- "EF, ~GH, IJ"
cols <- "~EF*~GH*IJ, EF*GH*~IJ, ~EF*GH*IJ, EF*~GH*~IJ"
chart <- makeChart(rows, cols)
chart
#     ~EF*~GH*IJ EF*GH*~IJ ~EF*GH*IJ EF*~GH*~IJ
# EF      -          x         -         x
# ~GH     x          -         -         x
# IJ      x          -         x         -

solveChart(chart)
# ~GH is redundant
#     EF + IJ
#      [,1]
# [1,]    1
# [2,]    3


# using implicant matrices
primes <- matrix(c(2,2,1,0,2,2,0,2,2,2), nrow = 2)
configs <- matrix(c(2,2,2,1,1,2,2,2,2,1,2,2,2,2,2), nrow = 3)
colnames(primes) <- colnames(configs) <- letters[1:5]

# the prime implicants: a~bce and acde
primes
#      a b c d e
# [1,] 2 1 2 0 2
# [2,] 2 0 2 2 2

# the initial causal combinations: a~bc~de, a~bcde and abcde
configs
#      a b c d e
# [1,] 2 1 2 1 2
# [2,] 2 1 2 2 2
# [3,] 2 2 2 2 2

chartLC <- makeChart(primes, configs, collapse = "")
chartLC
#         a~bc~de a~bcde   abcde
# a~bce      x       x       -
# acde       -       x       x

}

\keyword{functions}
