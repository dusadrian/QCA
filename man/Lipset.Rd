\name{_Lipset}

\docType{data}

\alias{LR}
\alias{LC}
\alias{LM}
\alias{LF}

\title{Lipset's indicators for the survival of democracy during the inter-war period.}

\description{
This dataset is taken from Lipset (1959), as used by Rihoux and De Meur (2009),
Cronqvist and Berg-Schlosser (2009) and Ragin (2009).
}

\usage{
data(LR)
data(LC)
data(LM)
data(LF)
}

\format{
A data frame containing 18 rows and the following 6 columns:

\tabular{rl}{
\bold{\code{DEV}}
    \tab Level of development: it is the GDP per capita (USD) in the raw data, calibrated in the \cr
    \tab binary crisp version to 0 if below 550 USD and 1 otherwise. For the multi-value crisp \cr
    \tab version, two thresholds were used: 550 and 850 USD. \cr

\bold{\code{URB}}
    \tab Level of urbanization: percent of the population in towns with 20000 or more \cr
    \tab inhabitants, calibrated in the crisp versions to 0 if below 50\% and 1 if above. \cr

\bold{\code{LIT}}
    \tab Level of literacy: percent of the literate population, calibrated in the crisp versions \cr
    \tab to 0 if below 75\% and 1 if above. \cr

\bold{\code{IND}}
    \tab Level of industrialization: percent of the industrial labor force, calibrated in the \cr
    \tab crisp versions to 0 if below 30\% and 1 if above. \cr

\bold{\code{STB}}
    \tab Government stability: a \dQuote{political-institutional} condition added to the previous \cr
    \tab four \dQuote{socioeconomic} ones. The raw data has the number of cabinets which governed \cr
    \tab in the period under study, calibrated in the crisp versions to 0 if 10 or above and to 1 \cr
    \tab if below 10. \cr

\bold{\code{SURV}}
    \tab Outcome: survival of democracy during the inter-war period: calibrated to 0 if negative, \cr
    \tab and 1 if positive raw data.\cr
}
}

\details{
There are four different versions of the Lipset data:

\tabular{rl}{
\bold{\code{LR}}  \tab contains the raw data \cr
\bold{\code{LC}}  \tab is the same data calibrated to binary crisp sets \cr
\bold{\code{LM}}  \tab is calibrated to multi-value sets \cr
\bold{\code{LF}}  \tab is calibrated to fuzzy-sets \cr
}

}

\references{
Lipset, S. M. (1959) \dQuote{Some Social Requisites of Democracy: Economic
Development and Political Legitimacy}, \emph{American Political Science Review}
vol.53, pp.69-105.

Cronqvist, L.; Berg-Schlosser, D. (2009) \dQuote{Multi-Value QCA (mvQCA)}, in
Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods. Qualitative
Comparative Analysis (QCA) and Related Techniques}, SAGE.

Rihoux, B.; De Meur, G. (2009) \dQuote{Crisp Sets Qualitative Comparative Analysis
(mvQCA)}, in Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods.
Qualitative Comparative Analysis (QCA) and Related Techniques}, SAGE.

Ragin, C. (2009) \dQuote{Qualitative Comparative Analysis Using Fuzzy-Sets (fsQCA)},
in Rihoux, B.; Ragin, C. (eds.) \emph{Configurational Comparative Methods.
Qualitative Comparative Analysis (QCA) and Related Techniques}, SAGE.
}

\keyword{datasets}
