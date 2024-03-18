\name{Implicant matrix functions: allExpressions, createMatrix, getRow}

\alias{allExpressions}

\alias{createMatrix}

\alias{getRow}

\title{Functions Related to the Implicant Matrix}

\description{
This is a set of functions dedicated to the implicant matrix, a space where all
causal configurations and their minimized solutions are found.

They can produce all possible implicants and prime implicants, or all possible
combinations for a specific number of causal conditions and their number
of values (either binary or multi-value).
}

\usage{
allExpressions(noflevels = NULL, arrange = FALSE, depth, raw = FALSE, ...)

createMatrix(noflevels = NULL, ...)

getRow(row.no = NULL, noflevels = NULL, zerobased = FALSE, ...)
}

\arguments{
  \item{noflevels}{The number of levels (values) for each causal condition.}
  \item{arrange}{Logical, if \code{TRUE} the result matrix is arranged for visual inspection.}
  \item{depth}{Integer, an upper number of causal conditions to form expressions with.}
  \item{raw}{Logical, if \code{TRUE} it returns the matrix indicating which conditions have
             been minimized, using \bold{\code{-1}}.}
  \item{row.no}{A vector, the desired row numbers.}
  \item{zerobased}{Logical, if \bold{\code{TRUE}} the first row number is zero.}
  \item{...}{Other arguments.}
}

\details{
A truth table for binary crisp conditions is a matrix with \eqn{2^k} rows, where
\eqn{k} is the number of causal conditions.

For multi-value causal conditions, the same equation can be generalised to:

\eqn{v_{1} \cdot v_{2} \cdot \dots \cdot v_{k}}{v1 * v2 * ... * vk}

where \eqn{v} is the number of values (levels) for every causal condition from \eqn{1}
to \eqn{k}.

Implicant matrices contain all rows from the truth table, plus all of their supersets,
(all implicants and prime implicants), including the empty set (Dusa 2007, 2010).

For a binary crisp set procedure, there are \eqn{3^k - 1} possible expressions (groupings),
see Ragin (2010). Including the empty set (the situation when all causal conditions have
been minimized), the implicant matrix consists of exactly \eqn{3^k} rows, including the
truth table configurations.

In fact, \eqn{3^k} is also obtained by the product:

\eqn{(2 + 1) \cdot (2 + 1) \cdot \dots \cdot (2 + 1)}{(2 + 1) * (2 + 1) * ... * (2 + 1)}

For multi-value causal conditions, the same equation can be generalised to:

\eqn{(v_{1} + 1) \cdot (v_{2} + 1) \cdot \dots \cdot (v_{k} + 1)}{(v1 + 1) * (v2 + 1) * ... * (vk + 1)}

where every number of levels in each causal conditions is incremented with 1, to allow
coding the minimization of literals in each (prime) implicant (see examples).

The function \bold{\code{allExpressions()}} creates a matrix which contains all possible implicants
and prime implicants, displayed in the original values form using the code \bold{\code{-1}} to
point the minimized literals, while the other functions use the code \bold{\code{0}}, all other
values being incremented with 1.

Specifying a smaller \bold{\code{depth}} automatically activates the argument \bold{\code{arrange}}.

When the argument \bold{\code{arrange}} is activated, the output is arranged in the increasing
order of the number of conditions which form conjunctions, up to the maximum number specified
by the argument \bold{\code{depth}} (which if \bold{\code{NULL}}, it is considered equal to the
number of columns in the matrix).

The function \bold{\code{createMatrix()}} creates a base matrix for truth tables and implicant
matrices.

The function \bold{\code{getRow()}} takes the number of a row in the truth table or implicant matrix
(in its decimal form), and transforms it into its binary (or multi-base) representation,
as a configuration of binary or multi-values for each causal condition.

Note that \eqn{\textsf{R}} is a 1-based language (all numbers start from 1), and similarly positions
in vectors and matrices start with 1. For this reason, although (mathematicall) the binary
representation of the decimal number 0 (for example, at three causal conditions) is 0 0 0,
in \eqn{\textsf{R}} that would be the \dQuote{first} line in the implicant matrix, therefore 0 0 0 is
translated into the number 1, unless the argument \bold{\code{zerobased}} is activated.
}

\value{
A matrix with \eqn{k} columns and:

\eqn{v_{1} \cdot v_{2} \cdot \dots \cdot v_{k}}{v1 * v2 * ... * vk} rows if a truth table;

\eqn{(v_{1} + 1) \cdot (v_{2} + 1) \cdot \dots \cdot (v_{k} + 1)}{(v1 + 1) * (v2 + 1) * ... * (vk + 1)} rows if an implicant matrix;

\eqn{x} rows, equal to the length of \bold{\code{row.no}}.
}

\author{
Adrian Dusa
}

\references{
Dusa, A. (2007b) \emph{Enhancing Quine-McCluskey}.
WP 2007-49, \href{https://compasss.org/working-papers-series/}{COMPASSS Working Papers series}.

Dusa, Adrian. 2010. \dQuote{A Mathematical Approach to the Boolean Minimization 
Problem.} \emph{Quality & Quantity} vol.44, no.1, pp.99-113.

Ragin, Charles C. (2000) \emph{Fuzzy-Set Social Science}. Chicago: University of 
Chicago Press.
}

\seealso{
\code{\link[base]{expand.grid}}
}

\examples{
# three binary causal conditions, having two levels each: 0 and 1
noflevels <- c(2, 2, 2)


# for three binary causal conditions
allExpressions(noflevels)

# the same matrix, this time arranged better
# (last rows represent the truth table)
allExpressions(noflevels, arrange = TRUE)

# show only the implicants (excluding the truth table)
allExpressions(noflevels, arrange = TRUE, depth = 2)

# using the raw form
allExpressions(noflevels, raw = TRUE)


# create a base truth table for 3 binary conditions
createMatrix(noflevels)

# its implicant matrix
createMatrix(noflevels + 1)

# create a base truth table where the second condition has three levels
createMatrix(c(2, 3, 2))



# deriving rows
rows <- c(2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17)
mat <- getRow(rows, noflevels + 1) # note the +1
rownames(mat) <- rows
colnames(mat) <- c("A", "B", "C")
mat

# implicant matrix     normal values
# 
#       A  B  C    |       A  B  C       
#    2  0  0  1    |    2  -  -  0       ~C      
#    4  0  1  0    |    4  -  0  -       ~B
#    5  0  1  1    |    5  -  0  0       ~B~C
#    7  0  2  0    |    7  -  1  -       B
#    8  0  2  1    |    8  -  1  0       B~C
#   10  1  0  0    |   10  0  -  -       ~A  
#   11  1  0  1    |   11  0  -  0       ~A~C                 
#   13  1  1  0    |   13  0  0  -       ~A~B   
#   14  1  1  1    |   14  0  0  0       ~A~B~C
#   16  1  2  0    |   16  0  1  -       ~AB
#   17  1  2  1    |   17  0  1  0       ~AB~C                   
}

\keyword{functions}
