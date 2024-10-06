\name{minimize}

\alias{minimize}
\alias{eqmcc}

\title{Minimize a truth table}

\description{
This function performs the QCA minimization of an input truth table, or if the
input is a dataset the minimization it minimizes a set of causal conditions with
respect to an outcome. Three minimization methods are available: the classical
Quine-McCluskey, the enhanced Quine-McCluskey and the latest Consistency Cubes
algorithm that is built for performance.

All algorithms return the same, exact solutions, see
\href{https://journal.r-project.org/archive/2018/RJ-2018-080/}{Dusa (2018)}
and \href{https://www.tandfonline.com/doi/full/10.1080/0022250X.2014.897949}{Dusa and Thiem (2015)}.
}

\usage{
minimize(input, include = "", dir.exp = NULL, details = FALSE, pi.cons = 0,
         sol.cons = 0, all.sol = FALSE, row.dom = FALSE, first.min = FALSE,
         max.comb = 0, use.labels = FALSE, method = "CCubes", ...)
}

\arguments{
  \item{input}{A truth table object (preferred) or a data frame containing
        calibrated causal conditions and an outcome.}
  \item{include}{A vector of other output values to include in the minimization
        process.}
  \item{dir.exp}{Character, a vector of directional expectations to derive the
        intermediate solution.}
  \item{details}{Logical, print more details about the solution.}
  \item{pi.cons}{Numerical fuzzy value between 0 and 1, minimal consistency
        threshold for a prime implicant to be declared as sufficient.}
  \item{sol.cons}{Numerical fuzzy value between 0 and 1, minimal consistency
        threshold for a model to be declared as sufficient.}
  \item{all.sol}{Logical, search for all possible models, including the non-minimal.}
  \item{row.dom}{Logical, perform row dominance in the prime implicants' chart to
        eliminate redundant prime implicants.}
  \item{first.min}{Logical, to return only the very first minimal solution (see
        Details).}
  \item{max.comb}{Numeric real, to limit the size of the PI chart (see Details).}
  \item{use.labels}{Logical, use category labels if present (see Examples).}
  \item{method}{Minimization method, one of "CCubes" (default), or "QMC" the
        classical Quine-McCluskey, or "eQMC" the enhanced Quine-McCluskey.}
  \item{...}{Other arguments, passed to other functions.}
}

\details{
Most of the times, this function takes a truth table object as the \bold{\code{input}}
for the minimization procedure, but the same argument can refer to a data frame
containing calibrated columns.

For the later case, the function \bold{\code{minimize()}} originally had some additional
formal arguments which were sent to the function \bold{\code{\link[QCA]{truthTable}()}}:
\bold{\code{outcome}}, \bold{\code{conditions}}, \bold{\code{n.cut}},
\bold{\code{incl.cut}}, \bold{\code{show.cases}}, \bold{\code{use.letters}} and
\bold{\code{inf.test}}.

All of these parameters are still possible with function \bold{\code{minimize()}}, but
since they are sent to the \bold{\code{\link[QCA]{truthTable}()}} function anyway, it is
unnecessary to duplicate their explanation here. The only situation which does need
an additional description relates to the argument \bold{\code{outcome}}, where
unlike \bold{\code{\link[QCA]{truthTable}()}} which accepts a single one, the function
\bold{\code{minimize()}} accepts multiple outcomes and performs a minimization for
each of them (a situation when all columns are considered causal conditions).

The argument \bold{\code{include}} specifies which other truth table rows are
included in the minimization process. Most often, the remainders are included but
any value accepted in the argument \bold{\code{explain}} is also accepted in the
argument \bold{\code{include}}.

The argument \bold{\code{dir.exp}} is used to specify directional expectations, as
described by Ragin (2003). They can be specified using SOP (sum of products)
expressions, which opens up the possibility to experiment with conjunctural
directional expectations. "Don't care" conditions are simply left unspecified.

If at least one of the conditions included in the analysis is multi-value, the
entire \bold{\code{dir.exp}} expression should be specified in multi-value
notation using squared brackets. If a condition X is crisp or fuzzy, the
multi-value notation X[0] is interpreted as its absence, as in the ~X notation.

Activating the \bold{\code{details}} argument has the effect of printing parameters
of fit for each prime implicant and each overall model, the essential prime
implicants being listed in the top part of the table. It also prints the truth
table, in case the argument \bold{\code{input}} has been provided as a data frame
instead of a truth table object.

The default method (when \bold{\code{all.sol = FALSE}}), is to find the minimal number
(\bold{\code{k}}) of prime implicants needed to cover all initial positive output
configurations (minterms), then exhaustively search through all possible disjunctions of
\bold{\code{k}} prime implicants which do cover those configurations.

Once the PI chart is constructed using the prime implicants found in the previous
stages, the argument \bold{\code{row.dom}} can be used to further eliminate
irrelevant prime implicants when solving the PI chart, applying the principle of row
dominance: if a prime implicant A covers the same (intial) positive output
configurations as another prime implicant B and in the same time covers
other configurations which B does not cover, then B is irrelevant and eliminated.

A large number of causal conditions (i.e. over 15), combined with a large number of cases
(i.e. hundreds) usually produce a very large number of prime implicants, resulting in a huge
and extremely complex PI chart with sometimes thousands of rows and hundreds of columns.

For such a complex PI chart, even finding a minimum is a formidable task, and exhaustively
solving it is very likely impossible in polynomial time. For this reason, after each level
of complexity the CCubes algorithm determines if the PI chart is too difficult, by
calculating the total number of combinations of minimum \bold{\code{k}} PIs necessary to
cover all columns.

The argument \bold{\code{max.comb}} controls this maximum number of combinations. It is a
rational number counted in (fractions of) billions, defaulted at zero to signal searching
to the maximum possible extent. If the total number of combinations exceeds a positive
value of \bold{\code{max.comb}}, the PI chart is determined as too complex, the search is
stopped and CCubes attempts to return all possible models using the PIs from the previous
levels of complexity, when the PI chart was still not too complex.

In the extreme situation even this is not feasible, the argument \bold{\code{first.min}}
controls returning only one (the very first found) minimal model, if at all possible.
}

\value{
An object of class \code{"qca"} when using a single outcome, or class \code{"mqca"}
when using multiple outcomes. These objects are lists having the following components:

\item{tt}{The truth table object.}
\item{options}{Values for the various options used in the function (including defaults).}
\item{negatives}{The line number(s) of the negative configuration(s).}
\item{initials}{The initial positive configuration(s).}
\item{PIchart}{A list containing the PI chart(s).}
\item{primes}{The prime implicant(s).}
\item{solution}{A list of solution model(s).}
\item{essential}{A list of essential PI(s).}
\item{pims}{A list of PI membership scores.}
\item{IC}{The matrix containing the inclusion and coverage scores for the model(s).}
\item{SA}{A list of simplifying assumptions.}
\item{i.sol}{A list of components specific to intermediate model(s), each having
      a PI chart, prime implicant membership scores, (non-simplifying) easy
      counterfactuals and difficult counterfactuals.}
\item{complex}{Flag solutions from a too complex PI chart}
\item{call}{The user's command which produced all these objects and result(s).}
}

\author{
Adrian Dusa
}

\references{
Cebotari, V.; Vink, M.P. (2013) \dQuote{A Configurational Analysis of Ethnic
Protest in Europe}. \emph{International Journal of Comparative Sociology}
vol.54, no.4, pp.298-324, \doi{10.1177/0020715213508567}.

Cebotari, V.; Vink, M.P. (2015) \dQuote{Replication Data for: A configurational
analysis of ethnic protest in Europe}, Harvard Dataverse, V2, \doi{10.7910/DVN/PT2IB9}.

Cronqvist, L.; Berg-Schlosser, D. (2009) \dQuote{Multi-Value QCA (mvQCA)}, in
Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods. Qualitative
Comparative Analysis (QCA) and Related Techniques}, SAGE.

Dusa, A.; Thiem, A. (2015) \dQuote{Enhancing the Minimization of Boolean and
Multivalue Output Functions With eQMC} \emph{Journal of Mathematical Sociology}
vol.39, no.2, pp.92-108,\cr
\doi{10.1080/0022250X.2014.897949}.

Dusa, A. (2018) \dQuote{Consistency Cubes: A Fast, Efficient Method for Boolean
Minimization}, R Journal vol.10, issue 2, pp. 357-370, \doi{10.32614/RJ-2018-080}

Dusa, A. (2019) \emph{QCA with R. A Comprehensive Resource}.
Springer International Publishing, \doi{10.1007/978-3-319-75668-4}.

Ragin, C. (2003) \emph{Recent Advances in Fuzzy-Set Methods and Their Application to
Policy Questions}. WP 2003-9,
\href{https://compasss.org/working-papers-series/}{COMPASSS Working Papers series}.

Ragin, C. (2009) \dQuote{Qualitative Comparative Analysis Using Fuzzy-Sets (fsQCA)},
in Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods.
Qualitative Comparative Analysis (QCA) and Related Techniques}, SAGE.

Ragin, C.C.; Strand, S.I. (2008) \dQuote{Using Qualitative Comparative
Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).}
\emph{Sociological Methods & Research} vol.36, no.4, pp.431-441,
\doi{10.1177/0049124107313903}.

Rihoux, B.; De Meur, G. (2009) \dQuote{Crisp Sets Qualitative Comparative Analysis
(mvQCA)}, in Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods.
Qualitative Comparative Analysis (QCA) and Related Techniques}, SAGE.

}

\seealso{\code{\link{truthTable}}, \code{\link[admisc]{factorize}}}

\examples{
\dontrun{
# -----
# Lipset binary crisp data

# the associated truth table
ttLC <- truthTable(LC, SURV, sort.by = "incl, n", show.cases = TRUE)
ttLC

# conservative solution (Rihoux & De Meur 2009, p.57)
cLC <- minimize(ttLC)
cLC

# view the Venn diagram for the associated truth table
library(venn)
venn(cLC)

# add details and case names
minimize(ttLC, details = TRUE)

# negating the outcome
ttLCn <- truthTable(LC, ~SURV, sort.by = "incl, n", show.cases = TRUE)
minimize(ttLCn)

# parsimonious solution, positive output
pLC <- minimize(ttLC, include = "?", details = TRUE)
pLC

# the associated simplifying assumptions
pLC$SA

# parsimonious solution, negative output
pLCn <- minimize(ttLCn, include = "?", details = TRUE)
pLCn


# -----
# Lipset multi-value crisp data (Cronqvist & Berg-Schlosser 2009, p.80)

# truth table, conditions all columns from DEV to IND
# note the sequence operator ":"
ttLM <- truthTable(LM, SURV, conditions = DEV:IND,
        sort.by = "incl", show.cases = TRUE)

# conservative solution, positive output
minimize(ttLM, details = TRUE)

# parsimonious solution, positive output
minimize(ttLM, include = "?", details = TRUE)

# negate the outcome
ttLMn <- truthTable(LM, ~SURV, conditions = DEV:IND,
         sort.by = "incl", show.cases = TRUE)

# conservative solution, negative output
minimize(ttLMn, details = TRUE)

# parsimonious solution, positive output
minimize(ttLMn, include = "?", details = TRUE)


# -----
# Lipset fuzzy sets data (Ragin 2009, p.112)

ttLF <- truthTable(LF, SURV, incl.cut = 0.8, sort.by = "incl", show.cases = TRUE)

# conservative solution
minimize(ttLF, details = TRUE)

# parsimonious solution
minimize(ttLF, include = "?", details = TRUE)

# intermediate solution
minimize(ttLF, include = "?", details = TRUE,
         dir.exp = c(DEV, URB, LIT, IND, STB))

# directional expectations can also be specified using a sequence
minimize(ttLF, include = "?", details = TRUE, dir.exp = DEV:STB)

# URB as a don't care condition (left unspecified) and
# conjunctural directional expectations
minimize(ttLF, include = "?", details = TRUE,
         dir.exp = c(DEV, STB, ~LIT*IND))


# -----
# Cebotari & Vink (2013, 2015)

ttCVF <- truthTable(CVF, outcome = PROTEST, incl.cut = 0.8,
                    sort.by = "incl, n", show.cases = TRUE)

pCVF <- minimize(ttCVF, include = "?", details = TRUE)
pCVF

# inspect the PI chart
pCVF$PIchart

# DEMOC*ETHFRACT*~POLDIS is dominated by DEMOC*ETHFRACT*GEOCON
# using row dominance to solve the PI chart
pCVFrd <- minimize(ttCVF, include = "?", row.dom = TRUE, details = TRUE)

# plot the prime implicants on the outcome
pims <- pCVFrd$pims

par(mfrow = c(2, 2))
for(i in 1:4) {
    XYplot(pims[, i], CVF$PROTEST, cex.axis = 0.6)
}


# -----
# temporal QCA (Ragin & Strand 2008) serving the input as a dataset,
# which will automatically be passed to truthTable() as an intermediary
# step before the minimization

minimize(RS, outcome = REC, details = TRUE)


# -----
# employ category labels

ttLF <- truthTable(LF, SURV, incl.cut = 0.8, sort.by = "incl", show.cases = TRUE)

minimize(ttLF, include = "?", use.labels = TRUE)

# or

ttLF <- truthTable(
    LF, SURV, incl.cut = 0.8, sort.by = "incl", show.cases = TRUE,
    use.labels = TRUE
)

minimize(ttLF, include = "?")

}
}

\keyword{functions}
