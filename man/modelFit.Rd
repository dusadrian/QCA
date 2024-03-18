\name{modelFit}
\alias{modelFit}

\title{Theory evaluation}

\description{
Function to enable theory evaluation, as introduced by Ragin (1987, p.118)
and extended Schneider & Wageman (2012, p.295), by producing parameters of fit
for all possible intersections between a given theoretical statement (a SOP
expresison) and the solutions found by function \bold{\code{\link[QCA]{minimize}()}}.

}

\usage{
modelFit(model, theory = "", select = NULL, ...)
}

\arguments{
  \item{model}{A minimization object of class \code{"QCA_min"}.}
  \item{theory}{Character, a SOP expression.}
  \item{select}{Character or numerical vector to select one or more models.}
  \item{...}{Other arguments, mainly for internal use.}
}

\details{
Following Ragin's (1987) original work, theory evaluation amounts to intersecting
a theoretical expectation with a model resulting from a minimization process.

There are in fact four intersections: presence - presence, presence - absence,
absence - presence and absence - absence, where by \dQuote{absence} is actually
meant a negation of an expression using the function \bold{\code{\link[admisc]{negate}()}}.

When multiple models exist, all of them are automatically detected, negated and
intersection with the theory. Intersections and parameters of fit are going to
be produced using a single theoretical expression.

In case of high model diversity, it is sometimes useful to select only a subset
to fit against theory. In such situations, the argument \code{select} can be
provided with either the name(s) of the model (for instance "C1P5" for intermediate
solutions) or simply the number of the model(s) of interest (for conservative
and parsimonious solutions).
}

\value{
A list containing objects of class \code{"QCA_pof"} with the parameters of fit. For a
single theoretical expression and a single model, the object is a simple
\code{"QCA_pof"} object.
}


\author{
Adrian Dusa
}

\references{

Ragin, C.C. (1987) \emph{The Comparative Method: Moving beyond Qualitative and
Quantitative Strategies}. Berkeley: University of California Press.

Schneider, C.Q.; Wagemann, C. (2012) \emph{Set-Theoretic Methods
for the Social Sciences: A Guide to Qualitative Comparative Analysis (QCA)}.
Cambridge: Cambridge University Press.
}

\seealso{\code{\link{intersection}}, \code{\link{negate}}, \code{\link{pof}}}

\examples{
# Lipset fuzzy version data

ttLF <- truthTable(LF, outcome = SURV, incl.cut = 0.8)


# parsimonious solution
pLF <- minimize(ttLF, include = "?")

# modelFit(model = pLF, theory = "DEV*STB")


# hypothetical exploration of intermediate solutions
iLF <- minimize(ttLF, include = "?", dir.exp = "1,0,0,0,0")

# modelFit(iLF, "DEV*STB")
}

\keyword{functions}
