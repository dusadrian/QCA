\name{findRows}

\alias{findRows}

\title{Find untenable configurations}

\description{
This function finds various types of untenable assumptions that are used when
excluding certain configurations from the minimization process.
}

\usage{
findRows(obj = NULL, expression = "", observed = FALSE, type = 1, ...)
}

\arguments{
  \item{obj}{A truth table (an object of class \code{"QCA_tt"}) or an equivalent
        numerical matrix.}
  \item{expression}{String: a QCA expression written in sum of products form.}
  \item{observed}{Logical: also return subset relations for observed
        configurations, when \code{obj} is a truth table.}
  \item{type}{Numeric vector, specifying the type(s) of untenable configurations.}
  \item{...}{Additional arguments to be passed to function
        \code{\link{truthTable}()}, for the negation of the outcome.}
}

\details{
The primary purpose is to find untenable assumptions to be excluded from the
Boolean minimization process. For this reason, the input is most of the times a
truth table, but for demonstration purposes it can also be a simple matrix
having column names.

It started as a function to find rows that are subsets of a given SOP expression,
and it developed to cover even more untenable asumptions.

Subset rows can be anything, from remainders to the observed configurations:
positive output, negative output and contradictions). By default, the function
returns only the subset configurations for the remaiders, but activating the
argument \bold{\code{observed}} adds the corresponding observed configurations
to the output.

It might ocasionally find negative output configurations or contradictions, but
that doesn't have any side effect because they are going to be excluded from the
minimization anyways, unless contradictions are included in the minimization.
The only category that really matters if they are identified or not, are the
positive output configurations.

The contradictory simplifying assumptions (CSAs) are those which are used for
both the presence and the absence of the outcome, while simultaneous subset
relations (SSRs) when observed configurations are sufficient for both the
presence and the absence of the outcome. CSAs and SSRs are incoherent
conterfactuals, part of a category called Untenable Assumptions.

This function takes does what is normally done with a series of commands, in a
more integrated and systematic way.

Providing a truth table is sufficient to perform all these tasks, because a
truth table already contains all necessary information of how it was produced,
most importantly the inclusion cut-off(s). By default, it uses the same options
to produce a truth table for the negation of the outcome (if the input truth
table was created for its presence, or the other way round), and minimizes both
to inspect their simplifying assumptions to detect which are contradictory.

Identical simplifying assumptions that found in both parsimonious solutions are
declared as contradictory. Observed configurations that are sufficient for both
the presence and the absence of the outcome are incoherent because of the
simultaneous subset relations problem.

The following types of untenable assumptions can be searched for:

\tabular{rl}{
    \code{0} \tab all of them\cr
    \code{1} \tab subsets of a given \bold{\code{expression}} (default)\cr
    \code{2} \tab contradictory simplifying assumptions\cr
    \code{3} \tab simultaneous subset relations\cr
}

To find contradictory simplifying assumptions, a truth table for the negated
outcome is constructed, using the \bold{\code{incl.cut}} argument from the
\bold{\code{obj}} input object. If the inclusion cut-off has a single value, the
same is used for the negated outcome, and if it has two values the second is used.

If very specific cutoff values are needed for the negation of the outcome, these
can be provided via the \bold{\code{...}} argument, that will be passed to
function \bold{\code{\link{truthTable}()}}.
}


\value{
A numeric vector of row numbers from the truth table.
}

\author{
Adrian Dusa
}

\seealso{\code{\link{truthTable}}, \code{\link{minimize}}}

\examples{
# Lipset's binary crisp version
ttLC <- truthTable(LC, "SURV", show.cases = TRUE)

findRows(obj = ttLC, "DEV*~IND*STB")

## all subset rows from the truth table, also for observed configurations
findRows(obj = ttLC, "DEV*~IND*STB", observed = TRUE)


# Lipset's fuzzy version
ttLF <- truthTable(LF, outcome = "SURV", incl.cut = 0.8)

findRows(obj = ttLF, type = 2) # contradictory simplifying assumptions


# Contradictory simplifying assumptions using different cutoff values
# for the _negation_ of the outcome

findRows(obj = ttLF, type = 2, incl.cut = 0.9, pri.cut = 0.7)
}

\keyword{functions}
