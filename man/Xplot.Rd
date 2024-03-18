\name{Xplot}

\alias{Xplot}

\title{Display the distribution of points for a single condition}

\description{
This function creates a plot for a single vector of numerical values, arranging them
horizontally on the X axis from minimum to maximum.
}

\usage{
Xplot(x, jitter = FALSE, at = pretty(x), ...)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{jitter}{Logical, vertically jitter the points.}
  \item{at}{The points at which tick-marks are to be drawn. Non-finite (infinite, NaN or NA)
            values are omitted. By default, tickmark locations are automatically
            computed, see the help file for \bold{\code{?pretty}}.}
  \item{...}{Other graphical parameters from \bold{\code{?par}}}
}

\details{
This is a special type of (scatter)plot, with points being arranged only on the
horizontal axis (it has no vertical axis). Useful when inspecting if points
are grouped into naturally occuring clusters, mainly for crisp calibration purposes.

The argument \bold{\code{...}} is used to pass arguments to the various graphical
parameters from \code{?par}, and also to the settings from \code{?jitter}.

The points have a default \bold{\code{cex}} (character expansion) value of 1, and
a default \bold{\code{pch}} value of 1 (empty points), which can be modified
accordingly (for instance value 21 for filled points). When \bold{\code{pch = 21}},
the color for the margins of the points can be specified via the argument
\bold{\code{col}}, while the argument \bold{\code{bg}} will determine the fill color
of the points.

The axis labels have a default \bold{\code{cex.axis}} value of 0.8, which affects
both the tickmarks labels and the axis labels.

When jittering the points, default values of 0.5 are used for the parameters
\bold{\code{factor}} and \bold{\code{amount}}, on the horizontal axis. More details
can be found in the base function \code{\link[base]{jitter}()}.

Although the points are displayed in a single dimension, on the horizontal axis,
the R graphical window will still have the default squared shape, with a lot of 
empty space on the vertical axis. Users are free to create their custom code to 
determine the size of the graphics window, or simply resize it to a suitable
height.
}

\author{
Adrian Dusa
}

\seealso{\code{\link[graphics]{par}}, \code{\link[graphics]{text}}, \code{\link[base]{jitter}}}

\examples{

# Lipset's raw data
# plot the DEV (level of developent) causal condition
Xplot(LR$DEV)

# jitter the points vertically
Xplot(LR$DEV, jitter = TRUE)

# clip plotting between the range of min and max
Xplot(LR$DEV, jitter = TRUE, at = range(LR$DEV))
}


\keyword{functions}
