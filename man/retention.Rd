\name{retention}

\alias{retention}

\title{Compute the retention probability of a csQCA solution}

\description{
This function computes the retention probability for a csQCA solution, under various
perturbation scenarios. It only works with  bivalent crisp-set data, containing the
binary values 0 or 1.
}

\usage{
retention(data, outcome = "", conditions = "", incl.cut = 1, n.cut = 1,
         type = "corruption", dependent = TRUE, p.pert = 0.5, n.pert = 1, ...)
}

\arguments{
  \item{data}{A dataset of bivalent crisp-set factors.}
  \item{outcome}{The name of the outcome.}
  \item{conditions}{A string containing the condition variables' names,
        separated by commas.}
  \item{incl.cut}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{n.cut}{The minimum number of cases for a causal combination with a set
        membership score above 0.5, for an output function value of "0" or "1".}
  \item{type}{Simulate corruptions of values in the conditions ("corruption"),
        or cases deleted entirely ("deletion").}
  \item{dependent}{Logical, if \bold{\code{TRUE}} indicating DPA - Dependent
        Perturbations Assumption and if \bold{\code{FALSE}} indicating
        IPA - Independent Perturbations Assumption.}
  \item{p.pert}{Probability of perturbation under independent (IPA) assumption.}
  \item{n.pert}{Number of perturbations under dependent (DPA) assumption.}
  \item{...}{Other arguments, mainly for internal use.}
}

\details{

The argument \bold{\code{data}} requires a suitable data set, in the form of a data frame.
with the following structure: values of 0 and 1 for bivalent crisp-set variables.

The argument \bold{\code{outcome}} specifies the outcome to be explained, in upper-case
notation (e.g. \bold{\code{X}}). 

The argument \bold{\code{conditions}} specifies the names of the condition variables.
If omitted, all variables in \bold{\code{data}} are used except \bold{\code{outcome}}.

The argument \bold{\code{type}} controls which type of perturbations should be simulated
to calculate the retention probability.
When \bold{\code{type = "corruption"}}, it simulates changes of values in the conditions
(values of 0 become 1, and values of 1 become 0). When \bold{\code{type = "deletion"}},
it calculates the probability of retaining the same solution if a number of cases are
deleted from the original data.

The argument \bold{\code{dependent}} is a logical which choses between two categories of assumptions. 
If \bold{\code{dependent = TRUE}} (the default) it indicates DPA - Dependent Perturbations Assumption,
when perturbations depend on each other and are tied to a fixed number of cases, ex-ante
(see Thiem, Spohel and Dusa, 2016).
If \bold{\code{dependent = FALSE}}, it indicates IPA - Independent Perturbations Assumption, when
perturbations are assumed to occur independently of each other.

The argument \bold{\code{n.cut}} is one of the factors that decide which configurations
are coded as logical remainders or not, in conjunction with argument \bold{\code{incl.cut}}.
Those configurations that contain fewer than \bold{\code{n.cut}} cases with membership scores
above 0.5 are coded as logical remainders (\bold{\code{OUT = "?"}}). If the number of such
cases is at least \bold{\code{n.cut}}, configurations with an inclusion score of at least
\bold{\code{incl.cut}} are coded positive (\bold{\code{OUT = "1"}}), while configurations with an
inclusion score below \bold{\code{incl.cut}} are coded negative (\bold{\code{OUT = "0"}}).

The argument \bold{\code{p.pert}} specifies the probability of perturbation under the
IPA - independent perturbations assumption (when \bold{\code{dependent = FALSE}}).

The argument \bold{\code{n.pert}} specifies the number of perturbations under the
DPA - dependent perturbations assumption (when \bold{\code{dependent = TRUE}}). At least
one perturbation is needed to possibly change a csQCA solution, otherwise the solution will
remain the same (retention equal to 100\%) if zero perturbations occur under this argument.
}

\author{
Adrian Dusa
}


\references{
    
Thiem, A.; Spoehel, R.; Dusa, A. (2015) \dQuote{Replication Package for:
Enhancing Sensitivity Diagnostics for Qualitative Comparative Analysis: A Combinatorial
Approach}, Harvard Dataverse, V1. \doi{10.7910/DVN/QE27H9}

Thiem, A.; Spoehel, R.; Dusa, A. (2016) \dQuote{Enhancing Sensitivity
Diagnostics for Qualitative Comparative Analysis: A Combinatorial Approach.}
\emph{Political Analysis} vol.24, no.1, pp.104-120.
}

\examples{
# the replication data, see Thiem, Spohel and Dusa (2015)
dat <- data.frame(matrix(c(
    rep(1, 25), rep(0, 20), rep(c(0, 0, 1, 0, 0), 3),
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, rep(1, 7), 0, 1),
    nrow = 16, byrow = TRUE, dimnames = list(
    c("AT", "DK", "FI", "NO", "SE", "AU", "CA", "FR",
      "US", "DE", "NL", "CH", "JP", "NZ", "IE", "BE"),
    c("P", "U", "C", "S", "W"))
))


# calculate the retention probability, for 2.5% probability of data corruption
# under the IPA - independent perturbation assuption
retention(dat, outcome = "W", incl.cut = 1, type = "corruption",
       dependent = FALSE, p.pert = 0.025)

# the probability that a csQCA solution will change
1 - retention(dat, outcome = "W", incl.cut = 1, type = "corruption",
       dependent = FALSE, p.pert = 0.025)

}

\keyword{functions}
