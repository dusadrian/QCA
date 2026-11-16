\name{PI chart functions: makeChart, findmin, solveChart}

\alias{makeChart}
\alias{findmin}
\alias{solveChart}

\title{Create and solve a prime implicants chart}

\description{
\code{makeChart()} shows which prime implicants cover each positive truth
table configuration. \code{findmin()} tells you the fewest prime implicants
needed to cover them all. \code{solveChart()} shows which combinations of
prime implicants do that.
}

\usage{
makeChart(primes = "", configs = "", snames = "", mv = FALSE, collapse = "*", ...)

findmin(chart, type = c("hybrid", "ilp"), ...)

solveChart(chart, row.dom = FALSE, all.sol = FALSE, depth = NULL, max.comb = 0,
           first.min = FALSE, ...)
}

\arguments{
  \item{primes}{The prime implicants, given as expressions or as a matrix with
        one prime implicant per row.}
  \item{configs}{The positive truth table configurations to cover, given as
        expressions or as a matrix with one configuration per row.}
  \item{snames}{Condition names. Give them as a vector or as a single string
        with commas between the names.}
  \item{mv}{Set to \code{TRUE} to show condition values in multi-value QCA
        notation in the chart labels.}
  \item{collapse}{The separator between conditions in a chart label. The
        default is \code{"*"}. Use \code{""} to join them without a separator.}
  \item{chart}{A PI chart like the one returned by \code{makeChart()}, with
        \code{TRUE} and \code{FALSE} entries.}
  \item{type}{How \code{findmin()} finds the minimum number of prime implicants.
        Use \code{"hybrid"} (the default) or \code{"ilp"}.}
  \item{row.dom}{Set to \code{TRUE} to drop a prime implicant if another covers
        all the same configurations and at least one more.}
  \item{all.sol}{Set to \code{TRUE} to include larger solutions in which every
        prime implicant is needed for that solution to cover the chart.}
  \item{depth}{Accepted for compatibility. It currently has no effect in
        \code{solveChart()}.}
  \item{max.comb}{The maximum number of combinations to check when searching
        for smallest solutions, in billions. The default of \code{0} sets
        no limit.}
  \item{first.min}{Set to \code{TRUE} to stop at the first smallest solution
        when \code{all.sol = FALSE}.}
  \item{...}{Additional options. For example, \code{solveChart()} passes
        \code{type} along to \code{findmin()}.}
}

\details{
\subsection{Reading and building the chart}{
Each row is a prime implicant and each column is a positive truth table
configuration. A \code{TRUE} entry, printed as \code{x}, means the prime
implicant covers that configuration. If a column has just one \code{x},
its prime implicant is essential and must appear in every solution.

When writing conditions together, as in \code{A~BC}, supply \code{snames}
if the condition names cannot be inferred from the expressions.
Using \code{*} between conditions, as in \code{A*~B*C}, makes the
condition names explicit.

If you supply matrices, use one column per condition in the same order in
both matrices. In a prime implicant, \code{0} means the condition has been
eliminated. For a binary condition, use \code{1} for absence and \code{2}
for presence.
}

\subsection{Finding the fewest prime implicants}{
Both methods in \code{findmin()} aim to give the same minimum number of
prime implicants needed to cover every configuration in the chart. Both
can simplify the chart by setting aside prime implicants whose coverage
can be replaced by another row without changing that number.

With \code{"hybrid"}, QCA searches the chart itself first. If it reaches
its search limit before proving the minimum, the bundled
\code{lp_solve} solver finishes the search on the reduced chart. This can
save time on large charts or charts in which most rows cover only a few
configurations.

With \code{"ilp"} (integer linear programming), QCA sends the simplified
chart directly to a solver. It uses Gurobi if available and the bundled
\code{lp_solve} otherwise. Set \code{gurobi = FALSE} to use
\code{lp_solve} from the start.
}

\subsection{Finding the solutions}{
By default, \code{solveChart()} uses \code{findmin()} to learn how many
prime implicants it needs. It then looks for every combination of that
size that covers all the columns. It rules out groups of rows as soon as
it can tell they cannot cover the remaining configurations.

Different prime implicants can appear in different solutions even when
they cover the same configurations. Setting \code{row.dom = TRUE} can
speed up the search by dropping rows whose coverage is contained in a
row that covers more. The minimum number of prime implicants stays the
same, but some alternative solutions may disappear.

Use \code{all.sol = TRUE} to look for larger solutions too. A larger
solution qualifies if dropping any of its prime implicants would leave
a configuration uncovered. This search uses the whole chart even if
\code{row.dom = TRUE} and can take much longer.

When looking for smallest solutions, use \code{first.min = TRUE} if one
solution is enough. A positive \code{max.comb} limits how many combinations
of rows QCA checks. For example, \code{max.comb = 0.001} limits the search
to one million combinations. If the limit is reached, the result can
contain only the solutions found so far.
}
}

\value{
\code{makeChart()} returns a matrix showing whether each prime implicant
covers each configuration.

\code{findmin()} returns the minimum number of prime implicants, or
\code{0} if a configuration has no covering prime implicant. With
\code{solind = TRUE}, it marks the rows in one smallest solution with
\code{1} and the other rows with \code{0}.

\code{solveChart()} returns a matrix with one solution per column. The
numbers in each column refer to rows of the original PI chart. If
solutions have different sizes, shorter columns are padded with \code{0}.
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
# Write conditions together when the names are clear from the configurations
chart <- makeChart("a, b, ~c", "abc, a~b~c, a~bc, ~ab~c")

# The same chart using unquoted expressions
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

# Use * to separate conditions with longer names
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


# Enter prime implicants and configurations as matrices
primes <- matrix(c(2,2,1,0,2,2,0,2,2,2), nrow = 2)
configs <- matrix(c(2,2,2,1,1,2,2,2,2,1,2,2,2,2,2), nrow = 3)
colnames(primes) <- colnames(configs) <- letters[1:5]

# The prime implicants are a~bce and acde
primes
#      a b c d e
# [1,] 2 1 2 0 2
# [2,] 2 0 2 2 2

# The configurations to cover are a~bc~de, a~bcde and abcde
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
