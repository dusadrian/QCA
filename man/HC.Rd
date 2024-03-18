\name{_Hino}

\docType{data}

\alias{HC}

\title{Time-Difference}

\description{
This data set was used by Hino (2009), to demonstrate the Time-Difference QCA.
}

\usage{
data(HC)
}

\format{                  
A data frame containing 15 cases (countries) and the following 5 columns:

\tabular{rl}{
    
\bold{\code{FOREIGN}}
    \tab Percentage of foreign population. \cr

\bold{\code{UNEMP}}
    \tab Percentage of unemployed population. \cr
    
\bold{\code{CONV}}
    \tab Party system convergence. \cr

\bold{\code{PRES80}}
    \tab Presence of extreme-right parties in 1980s. \cr
    
\bold{\code{VOTE}}
    \tab Outcome, vote share of extreme-right parties.\cr
}
}

\details{
For all columns in the data, a value of 1 means a positive difference between 1990
and 1980, and a value of 0 means negative or zero difference, except for the condition
CONV, which is the inverse of the condition DIVERT in the raw data. The condition
PRES80 does not have a time difference, it represents a simple presence / absence
of extreme-right parties in the 1980s.
}

\references{
Hiro, A. (2009) \dQuote{Time-Series QCA. Studying Temporal Change through Boolean Analysis}.
\emph{Sociological Theory and Methods}, vol.24, no.2, pp.247-265.
}

\keyword{datasets}
