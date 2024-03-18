\name{Parameters of fit}

\alias{pof}
\alias{pofind}

\title{Calculate parameters of fit}

\description{
These functions returns inclusion (consistency) and coverage, plus PRI for sufficiency
and RoN for necessity. The function \bold{\code{pofind()}} is a stripped down version
of the \bold{\code{pof()}} function, to calculate parameters of fit for single conditions.
}

\usage{
pof(setms = NULL, outcome = NULL, data = NULL, relation = "necessity",
    use.labels = FALSE, inf.test = "", incl.cut = c(0.75, 0.5), add = NULL, ...)
    
pofind(data = NULL, outcome = "", conditions = "", relation = "necessity",
    use.labels = FALSE, ...)
}

\arguments{
  \item{setms}{A data frame or a single vector of (calibrated) set memberships, or an
        expression written in sum of products form.}
  \item{outcome}{The name of the outcome column from a calibrated data frame, or
        the actual numerical column from the data frame, representing the outcome.}
  \item{data}{A calibrated data frame.}
  \item{conditions}{A single string containing the conditions' (columns) names
        separated by commas, or a character vector of conditions' names.}
  \item{relation}{The set relation to \bold{\code{outcome}}, either \bold{\code{"necessity"}}
        or \bold{\code{"sufficiency"}}, partial words like \bold{\code{"suf"}} being
        accepted (see examples).}
  \item{use.labels}{Logical, use category labels if present.}
  \item{inf.test}{Specifies the statistical inference test to be performed (currently 
       only \bold{\code{"binom"}}) and the critical significance level. It can be either
       a vector of length 2, or a single string containing both, separated by a comma.}
  \item{incl.cut}{The inclusion cutoff(s): either a single value for the presence of the
        output, or a vector of length 2, the second for the absence of the output. Used
        only in conjunction with the argument \bold{\code{inf.test}}}
  \item{add}{A function, or a list containing functions, to add more parameters of fit.}
  \item{...}{Other arguments to be passed to the main function.}
}

\details{
The function \bold{\code{pof()}} is one of the most flexible functions in the
\pkg{QCA} package. Depending on particular situations, its arguments can be provided
in various formats which are automatically recognized and treated accordingly.

When specified as a data frame, the argument \bold{\code{setms}} contains any kind
of set membership scores:

- calibrated causal conditions from the original data,

- membership scores from the resulting combinations (component \bold{\code{coms}})
of function \bold{\code{\link{superSubset}()}},

- prime implicant membership scores (component \bold{\code{pims}}) from function
\bold{\code{\link{minimize}()}},

- any other, custom created combinations of set memberships.

When specified as a matrix, \bold{\code{setms}} contains the crisp causal combinations
similar to those found in the truth table. The number of columns in the matrix should be
equal to the number of causal conditions in the original \bold{\code{data}}. If some of
them are minimized, they can be replaced by the numerical value \bold{\code{-1}} (see
examples section). 

More generally, \bold{\code{setms}} can be a numerical vector of line numbers from the
implicant matrix (see function \bold{\code{\link{createMatrix}()}}), which are automatically
transformed into their corresponding set membership scores. 

The argument \bold{\code{setms}} can also be a string expression, written in SOP - sum of
products form.

For all other situations when \bold{\code{setms}} is something else than a data frame, it
requires the original \bold{\code{data}} to generate the set memberships.

If character, the argument \bold{\code{outcome}} is the name of the column from the
original \bold{\code{data}}, to be explained (it is a good practice advice to specify
it using upper case letters, although it will nevertheless be converted to upper case,
by default).

If the outcome column is multi-value, the argument \bold{\code{outcome}} should use
the standard curly-bracket notation \bold{\code{X{value}}}. Multiple values are
allowed, separated by a comma (for example \bold{\code{X{1,2}}}). Negation of the
outcome can also be performed using the tilde \bold{\code{~}} operator, for example
\bold{\code{~X{1,2}}}, which is interpreted as: "all values in X except 1 and 2" and
it becomes the new outcome to be explained.

The argument \bold{\code{outcome}} can also be a numerical vector of set membership
values, either directly from the original data frame, or a recoded version (if
originally multi-value).

The argument \bold{\code{inf.test}} provides the possibility to perform statistical
inference tests, comparing the calculated inclusion score with a pair of thresholds
(\bold{\code{ic1}} and \bold{\code{ic0}}) specified in the argument \bold{\code{incl.cut}}.
Currently, it can only perform binomial tests (\bold{\code{"binom"}}), which means that
data should only be provided as binary crisp (not multivalue, not fuzzy).

If the critical significance level is not provided, the default level of \bold{\code{0.05}}
is taken.

The resulting object will contain the calculated p-values (\code{pval1} and \code{pval0}) from two
separate, one-tailed tests with the alternative hypothesis that the true inclusion
score is:\cr
- greater than \bold{\code{ic1}} (the inclusion cut-off for an output value of 1)\cr
- greater than \bold{\code{ic0}} (the inclusion cut-off for an output value of 0)

It should be noted that statistical tests are performing well only when the number
of cases is large, otherwise they are usually not significant.

For the necessity relation, the standard measures of inclusion and coverage are
supplemented with the \bold{\code{RoN}} (Relevance of Necessity) measure, as suggested by
Schneider & Wagemann's (2012).

The negation of both \bold{\code{setms}} and \bold{\code{outcome}} is accepted and
recognized using the Boolean subtraction from 1. If the names of the conditions are
provided via an optional (undocumented) argument \bold{\code{conditions}}, the
column names of the \bold{\code{setms}} object are negated using the function
\bold{\code{\link[admisc]{negate}()}}.

The logical argument \bold{\code{neg.out}} is deprecated, but backwards compatible.
\bold{\code{neg.out = TRUE}} and a tilde \bold{\code{~}} in the outcome name don't
cancel each other out, either one (or even both) signaling if the \bold{\code{outcome}}
should be negated.

The arguments from function \bold{\code{pofind()}} are passed to the main function
\bold{\code{pof()}} to calculate parameters of fit.

}

\author{
Adrian Dusa
}

\references{
Cebotari, V.; Vink, M.P. (2013) \dQuote{A Configurational Analysis of Ethnic
Protest in Europe}. \emph{International Journal of Comparative Sociology}
vol.54, no.4, pp.298-324, \doi{10.1177/0020715213508567}

Schneider, C. and Wagemann, C. (2012) \emph{Set-Theoretic Metods for the Social Sciences.
A Guide to Qualitative Comparative Analysis}. Cambridge: Cambridge University Press.
}

\seealso{\code{\link{minimize}}, \code{\link{superSubset}}}

\examples{
\dontrun{
# -----
# Cebotari & Vink (2013) fuzzy data

conds <- CVF[, 1:5]
PROTEST <- CVF$PROTEST

# parameters of fit (default is necessity)
pof(conds, PROTEST)

# parameters of fit negating the conditions
pof(1 - conds, PROTEST)

# negating the outcome
pof(conds, 1 - PROTEST)

# parameters of fit for sufficiency
pof(conds, PROTEST, relation = "suf")

# also negating the outcome
pof(conds, 1 - PROTEST, relation = "suf")


# -----
# standard analysis of necessity
# using the "coms" component from superSubset()
nCVF <- superSubset(CVF, outcome = PROTEST, incl.cut = 0.90, cov.cut = 0.6)

# also checking their necessity inclusion score in the negated outcome
pof(nCVF$coms, 1 - PROTEST)


# -----
# standard analysis of sufficiency
# using the "pims" component from minimize()

# conservative solution
cCVF <- minimize(truthTable(CVF, outcome = PROTEST, incl.cut = 0.8))

# verify if their negations are also sufficient for the outcome
pof(1 - cCVF$pims, PROTEST, relation = "sufficiency")


# -----
# using a SOP expression, translated using the function translate()

pof(~NATPRIDE + GEOCON -> PROTEST, data = CVF)

# same for the negation of the outcome
pof(~NATPRIDE + GEOCON -> ~PROTEST, data = CVF)

# necessity is indicated by the reverse arrow
pof(~NATPRIDE + GEOCON <- PROTEST, data = CVF)


# -----
# more parameters of fit, for instance Haesebrouck' consistency

inclH <- function(x, y) {
    sum(fuzzyand(x, y)) /
    sum(fuzzyand(x, y) + sqrt(fuzzyor(x - y, 0) * x))
}

pof(~NATPRIDE + GEOCON -> ~PROTEST, data = CVF, add = inclH)
}
}

\keyword{functions}
