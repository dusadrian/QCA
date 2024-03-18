\name{superSubset, findSubsets, findSupersets}

\alias{superSubset}
\alias{findSubsets}
\alias{findSupersets}

\title{Functions to find subsets or supersets}

\description{
Functions to find a list of implicants that satisfy some restrictions (see details),
or to find  the corresponding row numbers in the implicant matrix, for all subsets,
or supersets, of a (prime) implicant or an initial causal configuration.
}

\usage{
superSubset(data, outcome = "", conditions = "", relation = "necessity",
    incl.cut = 1, cov.cut = 0, ron.cut = 0, pri.cut = 0, depth = NULL,
    use.letters = FALSE, use.labels = FALSE, add = NULL, ...)

findSubsets(input, noflevels = NULL, stop = NULL, ...)

findSupersets(input, noflevels = NULL, ...)
}

\arguments{
  \item{data}{A data frame with crisp (binary and multi-value) or fuzzy
        causal conditions}
  \item{outcome}{The name of the outcome.}
  \item{conditions}{A string containing the conditions' names,
        separated by commas.}
  \item{relation}{The set relation to \bold{\code{outcome}}, either
        \bold{\code{"necessity"}}, \bold{\code{"sufficiency"}},
        \bold{\code{"necsuf"}} or \bold{\code{"sufnec"}}. Partial
        words like \bold{\code{"suf"}} are accepted.}
  \item{incl.cut}{The minimal inclusion score of the set relation.}
  \item{cov.cut}{The minimal coverage score of the set relation.}
  \item{ron.cut}{The minimal score for the \bold{\code{RoN}} - relevance of necessity.}
  \item{pri.cut}{The minimal score for the \bold{\code{PRI}} - proportional reduction in
        inconsistency.}
  \item{use.letters}{Logical, use simple letters instead of original conditions'
        names.}
  \item{use.labels}{Logical, use category labels if present.}
  \item{noflevels}{A vector containing the number of levels for each causal
        condition plus 1 (all subsets are located in the higher dimension, implicant
        matrix)}
  \item{input}{A vector of row numbers where the (prime) implicants are located, or a 
        matrix of configurations (only for supersets).}
  \item{stop}{The maximum line number (subset) to stop at, and return}
  \item{depth}{Integer, an upper number of causal conditions to form expressions with.}
  \item{add}{A function, or a list containing functions, to add more parameters of fit.}
  \item{...}{Other arguments, mainly for backward compatibility.}
}

\details{
The function \bold{\code{\link{superSubset}()}} finds a list of implicants that satisfy
some restrictions referring to the inclusion and coverage with respect to the outcome,
under given assumptions of necessity and/or sufficiency.

Ragin (2000) posits that under the necessity relation, instances of the outcome
constitute a subset of the instances of the cause(s). Conversely, under the
sufficiency relation, instances of the outcome constitute a superset of the instances
of the cause(s).

When \bold{\code{relation = "necessity"}} the function finds all implicants which are
supersets of the  outcome, then eliminates the redundant ones and returns the
surviving (minimal) supersets, provided they pass the inclusion and coverage
thresholds. If none of the surviving supersets pass these thresholds, the function
will find disjunctions of causal conditions, instead of conjunctions.

When \bold{\code{relation = "sufficiency"}} it finds all implicants which are subsets
of the outcome, and similarly eliminates the redundant ones and return the surviving
(minimal) subsets.

When \bold{\code{relation = "necsuf"}}, the relation is interpreted as necessity, and
\bold{\code{cov.cut}} is automatically set equal to the inclusion cutoff
\bold{\code{incl.cut}}. The same automatic equality is made for
\bold{\code{relation = "sufnec"}}, when relation is interpreted as sufficiency.

The argument \bold{\code{outcome}} specifies the name of the outcome, and if
multi-value the argument can also specify the level to explain, using square brackets
notation.

Outcomes can be negated using a tilde operator \bold{\code{~X}}. The logical argument
\bold{\code{neg.out}} is now deprecated, but still backwards compatible. Replaced by
the tilde in front of the outcome name, it controls whether \bold{\code{outcome}} is
to be explained or its negation. If \bold{\code{outcome}} is from a multivalent
variable, it has the effect that the disjunction of all remaining values becomes the
new outcome to be explained. \bold{\code{neg.out = TRUE}} and a tilde \bold{\code{~}}
in the outcome name don't cancel each other out, either one (or even both) signaling
if the \bold{\code{outcome}} should be negated.

If the argument \bold{\code{conditions}} is not specified, all other columns in
\bold{\code{data}} are used.

Along with the standard measures of inclusion and coverage, the function also returns
\bold{\code{PRI}} for sufficiency and \bold{\code{RoN}} (relevance of necessity, see
Schneider & Wagemann, 2012) for the necessity relation.

A subset is a conjunction (an intersection) of causal conditions, with respect to a
larger (super)set, which is another (but more parsimonious) conjunction of causal
conditions.

All subsets of a given set can be found in the so called \dQuote{implicant matrix},
which is a \eqn{n^k} space, understood as all possible combinations of values in any
combination of bases \eqn{n}, each causal condition having three or more levels
(Dusa, 2007, 2010).

For every two levels of a binary causal conditions (values 0 and 1), there are three
levels in the implicants matrix:

\tabular{rl}{
     0 \tab to mark a minimized literal\cr
     1 \tab to replace the value of 0 in the original binary condition\cr
    -1 \tab to replace the value of 1 in the original binary condition
}

A prime implicant is a superset of an initial combination of causal conditions, and 
the reverse is also true: the initial combination is a subset of a prime implicant.

Any normal implicant (not prime) is a subset of a prime implicant, and in the
same time a superset of some initial causal combinations.

Functions \bold{\code{\link{findSubsets}()}} and \bold{\code{\link{findSupersets}()}} find:

\tabular{rl}{
  \tab - all possible such subsets for a given (prime) implicant, or\cr
  \tab - all possible supersets of an implicant or initial causal combination
}

in the implicant matrix.

The argument \bold{\code{depth}} can be used to impose an upper number of causal
conditions to form expressions with, it is the complexity level where the search is
stopped. Depth is set to a maximum by default, and the algorithm will always stop at
the maximum complexity level where no new, non-redundant prime implicants are found.
Reducing the depth below that maximum will also reduce computation time.

For examples on how to add more parameters of fit via argument \bold{\code{add}}, see
the function \bold{\code{\link{pof}()}}.
}

\author{
Adrian Dusa
}

\value{
The result of the \bold{\code{\link{superSubset}()}} function is an object of class "ss",
which is a list with the following components:

\item{incl.cov}{A data frame with the parameters of fit.}
\item{coms}{A data frame with the (m)embersip (s)cores of the resulting
      (co)mbinations.}

For \bold{\code{\link{findSubsets}()}} and \bold{\code{\link{findSupersets}()}}, a vector with the
row numbers corresponding to all possible subsets, or supersets, of a (prime)
implicant.
}

\references{
Cebotari, V.; Vink, M.P. (2013) \dQuote{A Configurational Analysis of Ethnic
Protest in Europe}. \emph{International Journal of Comparative Sociology}
vol.54, no.4, pp.298-324, \doi{10.1177/0020715213508567}.

Cebotari, Victor; Vink, Maarten Peter (2015) \emph{Replication Data for: A
configurational analysis of ethnic protest in Europe}, Harvard Dataverse, V2,
\doi{10.7910/DVN/PT2IB9}.

Dusa, A. (2007b) \emph{Enhancing Quine-McCluskey}.
WP 2007-49, \href{https://compasss.org/working-papers-series/}{COMPASSS Working Papers series}.

Dusa, Adrian (2010) \dQuote{A Mathematical Approach to the Boolean Minimization 
Problem.} \emph{Quality & Quantity} vol.44, no.1, pp.99-113, \doi{10.1007/s11135-008-9183-x}.

Lipset, S. M. (1959) \dQuote{Some Social Requisites of Democracy: Economic Development
and Political Legitimacy}, \emph{American Political Science Review} vol.53, pp.69-105.

Schneider, Carsten Q.; Wagemann, Claudius (2012) \emph{Set-Theoretic Methods 
for the Social Sciences: A Guide to Qualitative Comparative Analysis (QCA)}. 
Cambridge: Cambridge University Press.
}

\seealso{\code{\link{createMatrix}}, \code{\link{getRow}}}

\examples{
    
# Lipset binary crisp sets
ssLC <- superSubset(LC, "SURV")

library(venn)
x = list("SURV" = which(LC$SURV == 1),
         "STB" = which(ssLC$coms[, 1] == 1),
         "LIT" = which(ssLC$coms[, 2] == 1))
venn(x, cexil = 0.7)

# Lipset multi-value sets
superSubset(LM, "SURV")

# Cebotari & Vink (2013) fuzzy data
# all necessary combinations with at least 0.9 inclusion and 0.6 coverage cut-offs
ssCVF <- superSubset(CVF, outcome = "PROTEST", incl.cut = 0.90, cov.cut = 0.6)
ssCVF

# the membership scores for the first minimal combination (GEOCON)
ssCVF$coms$GEOCON

# same restrictions, for the negation of the outcome
superSubset(CVF, outcome = "~PROTEST", incl.cut = 0.90, cov.cut = 0.6)

# to find supersets or supersets, a hypothetical example using
# three binary causal conditions, having two levels each: 0 and 1
noflevels <- c(2, 2, 2)

# second row of the implicant matrix: 0 0 1
# which in the "normal" base is:      - - 0
# the prime implicant being: ~C
(sub <- findSubsets(input = 2, noflevels + 1))
#  5  8 11 14 17 20 23 26 


getRow(sub, noflevels + 1)

# implicant matrix   normal values
#      a  b  c    |       a  b  c       
#   5  0  1  1    |    5  -  0  0      ~b~c    
#   8  0  2  1    |    8  -  1  0      b~c
#  11  1  0  1    |   11  0  -  0      ~a~c
#  14  1  1  1    |   14  0  0  0      ~a~b~c
#  17  1  2  1    |   17  0  1  0      ~ab~c
#  20  2  0  1    |   20  1  -  0      a~c
#  23  2  1  1    |   23  1  0  0      a~b~c               
#  26  2  2  1    |   26  1  1  0      ab~c 


# stopping at maximum row number 20
findSubsets(input = 2, noflevels + 1, stop = 20)
#  5  8 11 14 17 20


# -----
# for supersets
findSupersets(input = 14, noflevels + 1)
#  2  4  5 10 11 13 14

findSupersets(input = 17, noflevels + 1)
#  2  7  8 10 11 16 17

# input as a matrix
(im <- getRow(c(14, 17), noflevels + 1))

# implicant matrix   normal values
#  14  1  1  1    |   14  0  0  0       ~a~b~c
#  17  1  2  1    |   17  0  1  0       ~ab~c


sup <- findSupersets(input = im, noflevels + 1)
sup
#  2  4  5  7  8 10 11 13 14 16 17


getRow(sup, noflevels + 1)

# implicant matrix   normal values
#      a  b  c    |       a  b  c       
#   2  0  0  1    |    2  -  -  0       ~c      
#   4  0  1  0    |    4  -  0  -       ~b
#   5  0  1  1    |    5  -  0  0       ~b~c
#   7  0  2  0    |    7  -  1  -       b
#   8  0  2  1    |    8  -  1  0       b~c
#  10  1  0  0    |   10  0  -  -       ~a  
#  11  1  0  1    |   11  0  -  0       ~a~c                 
#  13  1  1  0    |   13  0  0  -       ~a~b   
#  14  1  1  1    |   14  0  0  0       ~a~b~c
#  16  1  2  0    |   16  0  1  -       ~ab
#  17  1  2  1    |   17  0  1  0       ~ab~c
                             
}

\keyword{functions}
