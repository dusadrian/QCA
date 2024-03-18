\name{XYplot}

\alias{XYplot}

\title{Create an XY plot}

\description{
This function creates an XY plot from the first two columns of a dataframe/matrix, or from
two separate vectors of numeric values.
}

\usage{
XYplot(x, y, data, relation = "sufficiency", mguides = TRUE,
       jitter = FALSE, clabels, enhance = FALSE, model = FALSE, ...)
}

\arguments{
  \item{x}{Character, the name of the column from the \bold{\code{data}} for the X axis,
        or the coordinates of points in the plot (either a matrix/dataframe with at least
        two columns, or a vector of numerical values for the X axis), or a valid SOP
        expression.}
  \item{y}{Character, the name of the column from the \bold{\code{data}} for the Y axis,
        or the Y coordinates of points in the plot, optional if \bold{\code{x}}
        is a matrix/dataframe.}
  \item{data}{A calibrated dataset, only if \bold{\code{x}} and \bold{\code{y}} are names.}
  \item{relation}{The set relation to Y, either \bold{\code{"sufficiency"}} (default) or
        \bold{\code{"necessity"}}.}
  \item{mguides}{Logical, print the middle guides.}
  \item{jitter}{Logical, jitter the points.}
  \item{clabels}{A vector of case labels with the same length as \bold{\code{x}} and
        \bold{\code{y}}, or a logical vector of the same length as the number of rows
        in the \code{data} (if provided).}
  \item{enhance}{Logical, if \bold{\code{TRUE}} print the points using different characters for each
        of the five significant regions for process tracing.}
  \item{model}{Logical, for an enhanced plot specify if the SOP expression in argument
       \bold{\code{x}} is a solution model.}
  \item{...}{Other graphical parameters from \bold{\code{?par}}}
}

\details{
If \bold{\code{x}} is a dataframe or a matrix, the axes labels will be taken
from the column names of \bold{\code{x}}, otherwise they will be inferred from
the names of the \bold{\code{x}} and \bold{\code{y}} objects that are passed to
this function.

\bold{\code{x}} can also be a string containing either the name of the column
for the X axis, or two column names separated by a comma, referring to the X
and Y axis respectively. When \bold{\code{x}} contains both X and Y column names,
the next argument will be considered as the \bold{\code{data}}.

If \bold{\code{data}} is provided, and the names of the X and Y columns are valid
R statements, quoting them is not even necessary and they can be negated using
either a tilde \bold{\code{"~"}} or \bold{\code{"1 - "}}.

The numeric values should be restricted between 0 and 1, otherwise an error is
generated.

The XY plot will also provide inclusion and coverage scores for a sufficiency
(along with \bold{\code{PRI}}) or a necessity relation (along with \bold{\code{RoN}}).

The argument \bold{\code{x}} can also be a SOP - sum of products expression, in which
case the relation is determined by the usual forward arrow \bold{\code{"=>"}} for
sufficiency and backward arrow \bold{\code{"<="}} for necessity.

The argument \bold{\code{...}} is used to pass arguments to the various graphical
parameters from \code{?par}, and also to the settings from \code{?jitter}.

The points have a default \bold{\code{cex}} (character expansion) value of 0.8, and
a default \bold{\code{pch}} value of 21 (filled points), which can be modified
accordingly (for example with value 1 of empty points). When \bold{\code{pch = 21}},
the color for the margins of the points can be specified via the argument
\bold{\code{col}}, while the argument \bold{\code{bg}} will determine the fill color
of the points.

The axes' labels have a default \bold{\code{cex.axis}} value of 0.8, which affects
both the tickmarks labels and the axis labels.

When jittering the points, default values of 0.01 are used for the parameters
\bold{\code{factor}} and \bold{\code{amount}}, on both horizontal and vertical axes.

The argument \bold{\code{enhance}} does all the work for the shape of the points and
their colors, acording to the five regions specified by Schneider & Rohlfing (2016),
who augmented the classical XY plot with process tracing.

The default enhanced XY plot has even more settings when the input SOP expression is
a minimization model (different colors, different regions where to place the labels etc.),
available by activating the argument \bold{\code{model}}. The model is automatically
detected if the input for \bold{\code{x}} is a minimization object.

}

\value{
A list of x and y values, especially useful when the points are jittered. 
}

\author{
Adrian Dusa
}

\seealso{\code{\link{par}}, \code{\link{text}}, \code{\link{jitter}}}


\references{
Schneider, C.; Wagemann, C. (2012) \emph{Set-Theoretic Metods for the Social Sciences.
A Guide to Qualitative Comparative Analysis}. Cambridge: Cambridge University Press.

Cebotari, V.; Vink, M.P. (2013) \dQuote{A Configurational Analysis of Ethnic
Protest in Europe}. \emph{International Journal of Comparative Sociology}
vol.54, no.4, pp.298-324.

Schneider, C.; Rohlfing, I. (2016) \dQuote{Case Studies Nested in Fuzzy-set QCA on
Sufficiency. Formalizing Case Selection and Causal Inference}. \emph{Sociological
Methods and Research} vol.45, no.3, pp.536-568, \doi{10.1177/0049124114532446}
}

\examples{

# Cebotari & Vink (2013)
# necessity relation between NATPRIDE and PROTEST
XYplot(CVF[, 5:6])

# same using two numeric vectors
XYplot(CVF$NATPRIDE, CVF$PROTEST)

# same using two column names
XYplot(NATPRIDE, PROTEST, data = CVF)

# since they are valid R statements, it works even without quotes
# (this only works in normal R console, not in the GUI version)
XYplot(NATPRIDE, PROTEST, data = CVF)

# negating the X axis, using numeric vectors
XYplot(1 - CVF$NATPRIDE, CVF$PROTEST)

# same thing using quotes
XYplot(1 - NATPRIDE, PROTEST, data = CVF)

# using tilde for negation
XYplot(~NATPRIDE, PROTEST, data = CVF)

# different color for the points
XYplot(~NATPRIDE, PROTEST, data = CVF, col = "blue")

# using a different character expansion for the axes
XYplot(~NATPRIDE, PROTEST , data = CVF, cex.axis = 0.9)

# custom axis labels
XYplot(~NATPRIDE, PROTEST, data = CVF, xlab = "Negation of NATPRIDE",
       ylab = "Outcome: PROTEST")

# necessity relation
XYplot(~NATPRIDE, PROTEST, data = CVF, relation = "necessity")

# jitter the points
XYplot(~NATPRIDE, PROTEST, data = CVF, jitter = TRUE)

# jitter with more amount
XYplot(~NATPRIDE, PROTEST, data = CVF, jitter = TRUE, amount = 0.02)

# adding labels to points
XYplot(~NATPRIDE, PROTEST, data = CVF, jitter = TRUE, cex = 0.8,
       clabels = rownames(CVF))

# or just the row numbers, since the row names are too long
XYplot(~NATPRIDE, PROTEST, data = CVF, jitter = TRUE, cex = 0.8,
       clabels = seq(nrow(CVF)))


# using a SOP expression (necessity relation)
XYplot(NATPRIDE <- ~PROTEST, data = CVF, jitter = TRUE, cex = 0.8,
       clabels = seq(nrow(CVF)))


#-----
# enhanced XY plot for process tracing
XYplot(~NATPRIDE, PROTEST, data = CVF, enhance = TRUE, jitter = TRUE)


# enhanced XY plot for a solution model
ttCVF <- truthTable(CVF, outcome = PROTEST, incl.cut = 0.85)
pCVF <- minimize(ttCVF, include = "?")
XYplot(pCVF$solution[[1]], PROTEST, data = CVF, enhance = TRUE)


# same plot, using the solution as a SOP expression
XYplot(~NATPRIDE + DEMOC*GEOCON*POLDIS + DEMOC*ETHFRACT*GEOCON,
      PROTEST, data = CVF, enhance = TRUE, model = TRUE)

}


\keyword{functions}
