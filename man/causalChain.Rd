\name{causalChain}

\alias{causalChain}

\title{Perform CNA - coincidence analysis using QCA}

\description{
This function mimics the functionality in the package \pkg{cna}, finding all possible
necessary and sufficient solutions for all possible outcomes in a specific dataset.
}

\usage{
      
causalChain(data, ordering = NULL, strict = FALSE, pi.cons = 0, pi.depth = 0,
      sol.cons = 0, sol.cov = 1, sol.depth = 0, ...)
}

\arguments{
  \item{data}{A data frame containing calibrated causal conditions.}
  \item{ordering}{A character string, or a list of character vectors specifying
        the causal ordering of the causal conditions.}
  \item{strict}{Logical, prevents causal conditions on the same temporal level to
        act as outcomes for each other.}
  \item{pi.cons}{Numerical fuzzy value between 0 and 1, minimal consistency
        threshold for a prime implicant to be declared as sufficient.}
  \item{pi.depth}{Integer, the maximum number of causal conditions in each
        prime implicant in the CCubes search. Positive values impose an explicit
        limit regardless of pi.cons. The default 0 leaves the depth unrestricted
        by the user and retains automatic stopping.}
  \item{sol.cons}{Numerical fuzzy value between 0 and 1, minimal consistency 
        threshold for a model to be declared as sufficient.}
  \item{sol.cov}{Numerical fuzzy value between 0 and 1, minimal coverage 
        threshold for a model to be declared as necessary.}
  \item{sol.depth}{Integer, a maximum number of prime implicants in a
        consistency-based disjunctive model. The default 0 uses automatic
        enumeration without a user-imposed term-count limit.}
  \item{...}{Other arguments to be passed to functions \code{\link{minimize}()}
        and \code{\link{truthTable}()}.}
}

\details{
Although claiming to be a novel technique, coincidence analysis is yet another form
of Boolean minimization. What it does is very similar and results in the same set of
solutions as performing separate QCA analyses where every causal condition from the
\bold{\code{data}} is considered an outcome.

This function aims to demonstrate this affirmation and show that results from package
\pkg{cna} can be obtained with package \pkg{QCA}. It is not intended to offer a
complete replacement for the function \bold{\code{cna()}}, but only to
replicate its so called \dQuote{asf} - atomic solution formulas.

The three most important arguments from function \bold{\code{cna()}} have direct
correspondents in function \bold{\code{\link{minimize}()}}:

\tabular{rl}{
\bold{\code{con}} \tab corresponds to \bold{\code{sol.cons}}.\cr
\bold{\code{con.msc}} \tab corresponds to \bold{\code{pi.cons}}.\cr
\bold{\code{cov}} \tab corresponds to \bold{\code{sol.cov}}.
}

Two other arguments from function \bold{\code{cna()}} have been directly
imported in this function, to complete the list of arguments that generate the same
results.

The argument \bold{\code{ordering}} splits the causal conditions in different temporal
levels, where prior arguments can act as causal conditions, but not as outcomes for the
subsequent temporal conditions. One simple way to split conditions is to use a list
object, where different components act as different temporal levels, in the order of
their index in the list: conditions from the first component act as the oldest causal
factors, while those from the and the last component are part of the most recent temporal
level.

Another, perhaps simpler way to express the same thing is to use a single character,
where factors on the same level are separated with a comma, and temporal levels are
separated by the sign \bold{\code{<}}.

A possible example is: \bold{\code{"A, B, C < D, E < F"}}.

Here, there are three temporal levels and conditions A, B and C can act as causal factors
for the conditions D, E and F, while the reverse is not possible. Given that D, E and F
happen in a subsequent temporal levels, they cannot act as causal conditions for A, B or C.
The same thing is valid with D and E, which can act as causal conditions for F, whereas
F cannot act as a causal condition for D or E, and certainly not for A, B or C.

The argument \bold{\code{strict}} controls whether causal conditions from the same temporal
level may be outcomes for each other. If activated, none of A, B and C can act as causal
conditions for the other two, and the same thing happens in the next temporal level where
neither D nor E can be causally related to each other.

Although the two functions reach the same results, they follow different methods. The input
for the minimization behind the function \bold{\code{cna()}} is a coincidence list,
while in package \pkg{QCA} the input for the minimization procedure is a truth table. The
difference is subtle but important, with the most important difference that package
\pkg{cna} is not exhaustive.

To find a set of solutions in a reasonable time, the formal choice in package \pkg{cna} is
to deliberately stop the search at certain (default) depths of complexity. Users are free
to experiment with these depths from the argument \bold{\code{maxstep}}, but there is no
guarantee the results will be exhaustive.

On the other hand, the function \bold{\code{causalChain()}} and generally all related
functions from package \pkg{QCA} are spending more time to make sure the search is
exhaustive. Depths can be set via the arguments \bold{\code{pi.depth}} and
\bold{\code{sol.depth}}, but unlike package \pkg{cna} these are not mandatory.

By default, the package \pkg{QCA} employes a different search algorithm based on
Consistency Cubes (Dusa, 2018), analysing all possible combinations of causal conditions
and all possible combinations of their respective levels. The structure of the input
dataset (number of causal conditions, number of levels, number of unique rows in the
truth table) has a direct implication on the search time, as all of those characteristics
become entry parameters when calculating all possible combinations.

Consequently, two kinds of depth arguments are provided: 

\tabular{rl}{
\bold{\code{pi.depth}} \tab the maximum number of causal conditions needed to conjunctively\cr
 \tab construct a prime implicant; a positive value is a hard limit, even\cr
 \tab if no model can be found within that limit.\cr
\bold{\code{sol.depth}} \tab the maximum number of prime implicants needed to disjunctively\cr
 \tab build a solution model that covers all initial positive output configurations.
}

A positive \bold{\code{pi.depth}} limits the number of literals in each term
in the CCubes search, independently of \bold{\code{pi.cons}} and
\bold{\code{all.sol}}. For example, \bold{\code{pi.depth = 3}} excludes a term
such as \bold{\code{A*B*C*D}}, even if that term would reduce the number of
terms in a model. Any minimum found is relative to the permitted terms; the
unrestricted minimum can be smaller. The limit is not exceeded when the retained
prime implicants cannot cover the outcome.

With \bold{\code{pi.depth = 0}}, no user-imposed conjunction limit applies and
automatic stopping remains active. These depth controls concern the CCubes
search with remainders, which \bold{\code{causalChain()}} uses by default.
When calling \bold{\code{minimize()}} directly, include remainders with
\bold{\code{include = "?"}}.

For ordinary truth-table minimization, the minimum bound is the smallest number
of prime implicants needed to cover the positive configurations. It is searched
independently of \bold{\code{all.sol}}. The default
\bold{\code{all.sol = TRUE}} includes larger irredundant models as well;
\bold{\code{all.sol = FALSE}} returns only minimum-size models.

A complexity level producing no new prime implicants does not by itself prove
that the smallest bound has been reached: deeper terms can lower that bound.
CCubes uses the same automatic certification rules as \code{\link{minimize}}.
For binary tables, minimum-size enumeration can stop after proving the minimum
bound and that all models of that size are represented. Enumerating all
irredundant models instead requires a complete PI chart, established by a
completeness check or by finishing systematic generation. Proving the smallest
bound alone does not establish completeness of the larger models.

The certification checks have fixed work and memory budgets. Inconclusive
checks lead to continued systematic generation, not an assumption of optimality.
These rules do not impose a dataset-specific conjunction depth or guarantee a
short runtime on every input. Consistency-based searches do not use these
ordinary truth-table certificates. See \code{\link{minimize}} for resource
limits and the effects of explicit depth restrictions.

The argument \bold{\code{sol.cons}} introduces another method of solving the PI chart.
Normally, once the solution models are found among all possible combinations of
\bold{\code{k}} prime implicants, consistencies and coverages are subsequently
calculated. When \bold{\code{sol.cons}} is lower than 1, then models are searched
based on their consistencies, which should be at least equal to this threshold.

With \bold{\code{sol.cons > 0}}, models are searched in increasing numbers of
terms. The first successful size is the smallest bound for the generated prime
implicants, because all smaller sizes have already been checked or ruled out.
With \bold{\code{all.sol = TRUE}}, larger irredundant models are also retained.
An empty size level does not end the search: larger models can still exist.

When \bold{\code{max.comb = 0}}, automatic enumeration prunes a branch if it
contains an already accepted model, or if upper bounds on attainable consistency
and coverage rule it out. After completing a model size, enumeration stops only
if every remaining branch is redundant or cannot meet the thresholds at any
larger permitted size. Otherwise it continues. In particular, an unspecified
\bold{\code{sol.depth}} is no longer replaced by an arbitrary seven-term limit.

There is also a finite data-derived upper limit: a disjunction can be represented
on the observed data by selecting at most one membership-maximizing term per
case. A model with more terms therefore contains a proper subset with the same
memberships and fit. The search need not exceed the smaller of the number of
cases and the number of generated prime implicants. Positive truth-table row
count alone is not a valid model-size limit for fuzzy data.

A positive \bold{\code{sol.depth}} explicitly restricts this search. Likewise,
a positive \bold{\code{pi.depth}} restricts the available terms, and
\bold{\code{max.comb}} can truncate enumeration. Completeness and minimum-bound
claims are relative to those restrictions. Automatic enumeration stores the
models found rather than allocating an entry for every possible combination,
and checks for interrupts between work batches. It imposes no automatic time or
enumeration budget: difficult searches and large result sets can still require
substantial resources. The ordinary truth-table/Petrick path is unchanged.

When minimizing a dataset instead of a truth table, unless otherwise specified, the 
argument \bold{\code{incl.cut}} is automatically set to 0.5 whenever \bold{\code{pi.cons}}
and/or \bold{\code{sol.cons}} are used, then passed to the function
\bold{\code{truthTable()}}.
}

\value{
A list of length equal to the number of columns in the \bold{\code{data}}. Each component
contains the result of the QCA minimization for that specific column acting as an outcome.
}

\author{
Adrian Dusa
}


\seealso{\code{\link{minimize}}, \code{\link{truthTable}}}

\examples{
\dontrun{
# The following examples assume the package cna is installed
library(cna)
cna(d.educate, what = "a")

# same results with
cc <- causalChain(d.educate)
cc

# inclusion and coverage scores can be inspected for each outcome
cc$E$IC


# another example, function cna() requires specific complexity depths
cna(d.women, maxstep = c(3, 4, 9), what = "a")

# same results with, no specific depths are required
causalChain(d.women)


# multivalue data require a different function in package cna
mvcna(d.pban, ordering = list(c("C", "F", "T", "V"), "PB"),
      cov = 0.95, maxstep = c(6, 6, 10), what = "a")

# same results again, simpler command
causalChain(d.pban, ordering = "C, F, T, V < PB", sol.cov = 0.95)


# specifying a lower consistency threshold for the solutions
mvcna(d.pban, ordering = list(c("C", "F", "T", "V"), "PB"), con = .93,
      maxstep = c(6, 6, 10), what = "a")

# same thing with
causalChain(d.pban, ordering = "C, F, T, V < PB", pi.cons = 0.93,
            sol.cons = 0.95)


# setting consistency thresholds for the PIs, solutions and also
# a coverage threshold for the solution (note that an yet another
# function for fuzzy sets is needed in package cna)

dat2 <- d.autonomy[15:30, c("AU","RE", "CN", "DE")]
fscna(dat2, ordering = list("AU"), con = .9, con.msc = .85, cov = .85,
      what = "a")

# again, the same results using the same function:
causalChain(dat2, ordering = "AU", sol.cons = 0.9, pi.cons = 0.85,
            sol.cov = 0.85)
}
}

\keyword{functions}
