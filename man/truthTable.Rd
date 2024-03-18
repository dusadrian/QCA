\name{truthTable}
\alias{truthTable}

\title{Create a truth table}

\description{
Function to create a truth table from all types of calibrated data (binary
crisp, multi-value crisp and fuzzy).
For fuzzy data, an improved verson of Ragin's (2008) procedure is applied to
assign cases to the vector space corners (the truth table rows).
}

\usage{
truthTable(data, outcome = "", conditions = "", incl.cut = 1, n.cut = 1, pri.cut = 0,
           exclude = NULL, complete = FALSE, use.letters = FALSE, use.labels = FALSE,
           show.cases = FALSE, dcc = FALSE, sort.by = "", inf.test = "", ...)
}

\arguments{
  \item{data}{A data frame containing calibrated causal conditions and an outcome.}
  \item{outcome}{String, the name of the outcome.}
  \item{conditions}{A single string containing the conditions' (columns) names
        separated by commas, or a character vector of conditions' names.}
  \item{incl.cut}{The inclusion cut-off(s): either a single value for the presence of the
        output, or a vector of length 2, the second for the absence of the output.}
  \item{n.cut}{The minimum number of cases under which a truth table row is
        declared as a remainder.}
  \item{pri.cut}{The minimal score for the \bold{\code{PRI}} - proportional reduction in
        inconsistency, under which a truth table row is declared as negative.}
  \item{exclude}{A vector of (remainder) row numbers from the truth table, to code as
        negative output configurations.}
  \item{complete}{Logical, print complete truth table.}
  \item{use.letters}{Logical, use letters instead of causal conditions' names.}
  \item{use.labels}{Logical, use category labels if present.}
  \item{show.cases}{Logical, print case names.}
  \item{dcc}{Logical, if \bold{\code{show.cases = TRUE}}, the cases being
        displayed are the deviant cases consistency in kind.}
  \item{sort.by}{Sort the truth table according to various columns.}
  \item{inf.test}{Specifies the statistical inference test to be performed
       (currently only \bold{\code{"binom"}}) and the critical significance level.
       It can be either a vector of length 2, or a single string containing both,
       separated by a comma.}
  \item{...}{Other arguments (mainly for backward compatibility).}
}

\details{
The \bold{\code{data}} should always be provided as a data frame, with calibrated columns.

Calibration can be either crisp, with 2 or more values starting from 0, or fuzzy with
continous scores from 0 to 1. Raw data containing relative frequencies can
also be continous between 0 and 1, but these are not calibrated, fuzzy data.

Some columns can contain the placeholder \bold{\code{"-"}} indicating a \dQuote{don't care},
which is used to indicate the temporal order between other columns in tQCA. These
special columns are not causal conditions, hence no parameters of fit will be
calculated for them.

The argument \bold{\code{outcome}} specifies the column name to be explained.
If the outcome is a multivalue column, it can be specified in curly bracket notation,
indicating the value to be explained (the others being automatically converted to
zero).

The outcome can be negated using a tilde operator \bold{\code{~X}}. The logical
argument \bold{\code{neg.out}} is now deprecated, but still backwards compatible.
Replaced by the tilde in front of the outcome name, it controls whether \bold{\code{outcome}}
is to be explained or its negation. Note that using both \bold{\code{neg.out = TRUE}} and a
tilde \bold{\code{~}} in the outcome name cancel each other out.

If the outcome column is multi-value, the argument \bold{\code{outcome}} should use
the standard curly-bracket notation \bold{\code{X{value}}}. Multiple values are
allowed, separated by a comma (for example \bold{\code{X{1,2}}}). Negation of the
outcome can also be performed using the tilde \bold{\code{~}} operator, for example
\bold{\code{~X{1,2}}}, which is interpreted as: "all values in X except 1 and 2"
and it becomes the new outcome to be explained.

The argument \bold{\code{conditions}} specifies the causal conditions' names among the other
columns in the data. When this argument is not specified, all other columns except for
the outcome are taken as causal conditions.

A good practice advice is to specify both \bold{\code{outcome}} and
\bold{\code{conditions}} as upper case letters. It is possible, in a next version,
to negate outcomes using lower case letters, a situation where it really does
matter how the outcome and/or conditions are specified.

The argument \bold{\code{incl.cut}} replaces both (deprecated, but still backwards
compatible) former arguments \bold{\code{incl.cut1}} and \bold{\code{incl.cut0}}.
Most of the analyses use the inclusion cutoff for the presence of the output
(code \bold{\code{"1"}}). When users need both inclusion cutoffs (see below),
\bold{\code{incl.cut}} can be specified as a vector of length 2, in the form:
\bold{\code{c(ic1, ic0)}} where:

\tabular{rl}{
\bold{\code{ic1}} \tab is the inclusion cutoff for the presence of the output,\cr
                  \tab a minimum sufficiency inclusion score above which the output
                       value is coded with \code{"1"}.\cr \cr
\bold{\code{ic0}} \tab is the inclusion cutoff for the absence of the output,\cr
                  \tab a maximum sufficiency inclusion score below which the output
                       value is coded with \code{"0"}.\cr
}

If not specifically declared, the argument \bold{\code{ic0}} is automatically set
equal to \bold{\code{ic1}}, but otherwise \bold{\code{ic0}} should always be lower
than \bold{\code{ic1}}.

Using these two cutoffs, as well as \bold{\code{pri.cut}} the observed combinations are
coded with:

\tabular{rl}{
\bold{\code{"1"}} \tab if they have an inclusion score of at least \bold{\code{ic1}}\cr
                  \tab and a PRI score of at least \bold{\code{pri.cut}}\cr \cr
\bold{\code{"C"}} \tab if they have an inclusion score below \bold{\code{ic1}} and
                       at least \bold{\code{ic0}} (contradiction)\cr \cr
\bold{\code{"0"}} \tab if they have an inclusion score below \bold{\code{ic0}} or \cr
                  \tab a PRI score below \bold{\code{pri.cut}}
}

The argument \bold{\code{n.cut}} specifies the frequency threshold under which a truth
table row is coded as a remainder, irrespective of its inclusion score.

When argument \bold{\code{show.cases}} is set to \bold{\code{TRUE}}, the case names will be
printed at their corresponding row in the truth table. The resulting object always contains
the cases for each causal combination, even if not printed on the screen (the print function
can later be used to print them).

The \bold{\code{sort.by}} argument orders all configurations by any of the columns
present in the truth table. Typically, sorting occurs by their outcome value,
and/or by their inclusion score, and/or by their frequency, in any order.

Sorting decreasingly (the default) or increasingly can be specified adding the signs
\bold{\code{-}} or \bold{\code{+}}, next after the column name in argument
\bold{\code{sort.by}} (see examples). Note that \bold{\code{-}} is redundant because
it is the default anyways.

The order specified in this vector is the order in which the configurations will
be sorted. When sorting based on the OUTput column, remainders will always be sorted last.

The argument \bold{\code{use.letters}} controls using the original names of the causal
conditions, or replace them by single letters in alphabetical order. If the
causal conditions are already named with single letters, the original letters
will be used.

The argument \bold{\code{inf.test}} combines the inclusion score with a statistical
inference test, in order to assign values in the output column OUT. For the moment, it
is only the binomial test, which needs crisp data (it doesn't work with fuzzy sets).
Following a similar logic as above, for a given (specified) critical significance level,
the output for a truth table row will be coded as:

\tabular{rl}{
\bold{\code{"1"}} \tab if the true inclusion score is significanly higher than
                       \bold{\code{ic1}},\cr \cr
\bold{\code{"C"}} \tab contradiction, if the true inclusion score is not significantly
                       higher than \bold{\code{ic1}}\cr
                  \tab but significantly higher than \bold{\code{ic0}},\cr \cr
\bold{\code{"0"}} \tab if the true inclusion score is not significantly higher than
                       \bold{\code{ic0}}.\cr
}

It should be noted that statistical tests perform well only when the number of cases is
large, otherwise they are usually not significant. For a low number of cases, depending
on the inclusion cutoff value(s), it will be harder to code a value of \bold{\code{"1"}}
in the output, and also harder to obtain contradictions if the true inclusion is not
signficantly higher than \bold{\code{ic0}}.

The argument \bold{\code{exclude}} is used to exclude truth table rows from the
minimization process, from the positive configurations and/or from the remainders. This
is achieved by coding those configurations with a value of 0 in the OUTput column (thus
treating them as if they were observed as negative output configurations).

The argument \bold{\code{complete}} controls how to print the table on the screen, either
complete (when set to \bold{\code{TRUE}}), or just the observed combinations (default).
For up to 7 causal conditions, the resulting object will always contain the complete
truth table, even if it's not printed on the screen. This is useful for multiple reasons:
researchers like to manually change output values in the truth table (sometimes including
in this way a remainder, for example), and it is also useful to plot Venn diagrams, each
truth table row having a correspondent intersection in the diagram.


}

\value{
An object of class \code{"tt"}, a list containing the following components:

\tabular{rl}{
\bold{tt} \tab The truth table itself.\cr
\bold{indexes} \tab The line numbers for the observed causal configurations.\cr
\bold{noflevels} \tab A vector with the number of values for each causal condition.\cr
\bold{initial.data} \tab The initial data.\cr
\bold{recoded.data} \tab The crisp version of the \bold{\code{initial.data}}, if fuzzy.\cr
\bold{cases} \tab The cases for each observed causal configuration.\cr
\bold{DCC} \tab Deviant cases for consistency.\cr
\bold{minmat} \tab The membership scores matrix of cases in the observed truth table combinations.\cr
\bold{categories} \tab Category labels, if present in the data.\cr
\bold{multivalue} \tab Logical flag, if either conditions or the outcome are multivalue.\cr
\bold{options} \tab The command options used.\cr
}

}

\author{
Adrian Dusa
}

\references{
Cronqvist, L.; Berg-Schlosser, D. (2009) \dQuote{Multi-Value QCA (mvQCA)}, in
Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods. Qualitative
Comparative Analysis (QCA) and Related Techniques}, SAGE.

Dusa, A. (2019) \emph{QCA with R. A Comprehensive Resource}.
Springer International Publishing, \doi{10.1007/978-3-319-75668-4}.

Lipset, S.M. (1959) \dQuote{Some Social Requisites of Democracy: Economic Development
and Political Legitimacy}, \emph{American Political Science Review} vol.53, pp.69-105.

Ragin, C.C. (1987) \emph{The Comparative Method: Moving beyond Qualitative and
Quantitative Strategies}. Berkeley: University of California Press.

Ragin, C.C. (2008) \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}.
Chicago: University of Chicago Press.

Ragin, C.C.; Strand, S.I. (2008) \dQuote{Using Qualitative Comparative
Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).}
\emph{Sociological Methods & Research} vol.36, no.4, pp.431-441.

Schneider, C.Q.; Wagemann, C. (2012) \emph{Set-Theoretic Methods
for the Social Sciences: A Guide to Qualitative Comparative Analysis (QCA)}.
Cambridge: Cambridge University Press.
}

\seealso{\code{\link{minimize}}}

\examples{
# -----
# Lipset binary crisp data
ttLC <- truthTable(LC, "SURV")

# inspect the truth table
ttLC

# print the cases too, even if not specifically asked for
print(ttLC, show.cases = TRUE)

# the printing function also supports the complete version
print(ttLC, show.cases = TRUE, complete = TRUE)

# formally asking the complete version
truthTable(LC, "SURV", complete = TRUE)

# sorting by multiple columns, decreasing by default
truthTable(LC, "SURV", complete = TRUE, sort.by = "incl, n")

# sort the truth table decreasing for inclusion, and increasing for n
# note that "-" is redundant, sorting is decreasing by default
truthTable(LC, "SURV", complete = TRUE, sort.by = "incl-, n+")


# -----
# Lipset multi-value crisp data (Cronqvist & Berg-Schlosser 2009, p.80)
truthTable(LM, "SURV", sort.by = "incl")

# using a frequency cutoff equal to 2 cases
ttLM <- truthTable(LM, "SURV", n.cut = 2, sort.by = "incl")
ttLM

# the observed combinations coded as remainders
ttLM$removed


# -----
# Cebotari & Vink fuzzy data
ttCVF <- truthTable(CVF, "PROTEST", incl.cut = 0.8, sort.by = "incl")

# view the Venn diagram for this truth table
library(venn)
venn(ttCVF)

# each intersection transparent by its inclusion score
venn(ttCVF, transparency = ttCVF$tt$incl)

# the truth table negating the outcome
truthTable(CVF, "~PROTEST", incl.cut = 0.8, sort.by = "incl")

# allow contradictions
truthTable(CVF, "PROTEST", incl.cut = c(0.8, 0.75), sort.by = "incl")


# -----
# Ragin and Strand data with temporal QCA
# truth table containing the "-" placeholder as a "don't care"
truthTable(RS, "REC")

}

\keyword{functions}
