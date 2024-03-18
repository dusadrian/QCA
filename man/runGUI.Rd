\name{runGUI}

\alias{runGUI}

\title{run the GUI shiny app for the QCA package}

\description{
Runs the graphical user interface app based on the \pkg{shiny} package.
}

\usage{
runGUI(x)
}

\arguments{
  \item{x}{Path to the shiny app.}
}

\details{
This function is a wrapper to the \bold{\code{\link[shiny]{runApp}()}} function in package
\pkg{shiny}. If \bold{\code{x}} is not provided, it automatically locates the \bold{\code{gui}}
directory in the path where the QCA package has been installed, and runs it.

The user interface has an interactive R console in the webpage. Commands are
parsed and evaluated into a dedicated environment, with efforts to capture errors and
warnings.
}

\author{
Adrian Dusa
}

\keyword{functions}
