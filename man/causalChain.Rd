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
  \item{pi.depth}{Integer, a maximum number of causal conditions to be used
        when searching for conjunctive prime implicants.}
  \item{sol.cons}{Numerical fuzzy value between 0 and 1, minimal consistency 
        threshold for a model to be declared as sufficient.}
  \item{sol.cov}{Numerical fuzzy value between 0 and 1, minimal coverage 
        threshold for a model to be declared as necessary.}
  \item{sol.depth}{Integer, a maximum number of prime implicants to be used when
        searching for disjunctive models.}
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
 \tab construct a prime implicant; it is the complexity level where the search\cr
 \tab can be stopped, as long as the PI chart can be solved.\cr
\bold{\code{sol.depth}} \tab the maximum number of prime implicants needed to disjunctively\cr
 \tab build a solution model that covers all initial positive output configurations.
}

These arguments introduce a possible new way of deriving prime implicants and solution models,
that can lead to different results (i.e. even more parsimonious) compared to the classical
Quine-McCluskey. When either of them is modified from the default value of 0, the minimization
method is automatically set to \bold{\code{"CCubes"}} and the remainders are automatically
included in the minimization.

The higher these depths, the higher the search time. Connversely, the search time can be
significantly shorter if these depths are smaller. Irrespective of how large 
\bold{\code{pi.depth}} is, the algorithm will always stop at a maximum complexity level
where no new, non-redundant prime implicants are found. The argument \bold{\code{sol.depth}}
is relevant only when activating the argument \bold{\code{all.sol}} to solve the PI chart.

The argument \bold{\code{sol.cons}} introduces another method of solving the PI chart.
Normally, once the solution models are found among all possible combinations of
\bold{\code{k}} prime implicants, consistencies and coverages are subsequently
calculated. When \bold{\code{sol.cons}} is lower than 1, then models are searched
based on their consistencies, which should be at least equal to this threshold.

Exhaustiveness is guaranteed in package \pkg{QCA} precisely because it uses a truth table
as an input for the minimization procedure. The only exception is the option of finding
solutions based on their consistency, with the argument \bold{\code{sol.cons}}: for large
PI charts, time can quickly increase to infinity, to identify all possible irredundant
(disjunctions that are not subsets of previously found) disjunctive models. In such a
situation, the number of combinations of all possible numbers of prime implicants is
potentially too large to be solved in a polynomial time and if not otherwise specified in
the argument \bold{\code{sol.depth}} the function \bold{\code{causalChain()}} silently
sets a complexity level of 7 prime implicants per model.

When minimizing a dataset instead of a truth table, unless otherwise specified, the 
argument \bold{\code{incl.cut}} is automatically set to the minimum value between
\bold{\code{pi.cons}} and \bold{\code{sol.cons}}, then passed to the function
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
