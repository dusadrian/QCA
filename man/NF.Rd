\name{_Nieuwbeerta}

\docType{data}

\alias{NF}

\title{Class voting in post-World War era}

\description{
This fuzzy dataset is an adaptation from Ragin (2005, 2008), the data itself being attributed
to Nieuwbeerta (1995).
}

\usage{
data(NF)
}

\format{                  
A data frame containing 12 cases (countries) and the following 5 columns:

\tabular{rl}{
    
\bold{\code{A}}
    \tab degree of membership in the set of highly \emph{affluent} countries \cr

\bold{\code{I}}
    \tab degree of membership in the set of countries with substantial levels of income \emph{inequality}\cr
    
\bold{\code{M}}
    \tab degree of membership in the set of countries with a high percentage of workers employed \cr
    \tab in \emph{manufacturing} \cr

\bold{\code{U}}
    \tab degree of membership in the set of countries with strong \emph{unions} \cr
    
\bold{\code{W}}
    \tab outcome: degree of membership in the set of countries with \emph{weak} class voting \cr
  }
}

\details{
All fuzzy sets in this data are constructed on a six-values scale, for demonstrative purposes.

In the original dataset, the outcome \code{W} is presented as the first column.
}

\references{
Nieuwbeerta, P. (1995) \emph{The Democratic Class Struggle in Twenty Countries: 1945:1990}. 
Amsterdam: Thesis Publishers.

Ragin, C.C. (2005) \dQuote{From fuzzy sets to crisp truth tables}.
WP 2004-28, \href{https://compasss.org/working-papers-series/}{COMPASSS Working Papers series}.

Ragin, C.C. (2008) \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}. 
Chicago: University of Chicago Press.


}

\keyword{datasets}
