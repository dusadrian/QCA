\name{_Ragin and Strand}

\docType{data}

\alias{RS}

\title{University recognition of a graduate student union}

\description{
Original data used by Caren and Panofsky (2005), and reanalysed by Ragin 
and Strand (2008).
}

\usage{data(RS)}

\format{
A data frame containing 17 cases and the following 6 columns:

\tabular{rl}{
\bold{\code{P}}   \tab Public university               \cr
\bold{\code{E}}   \tab Support of elite allies         \cr
\bold{\code{A}}   \tab National union affiliation      \cr
\bold{\code{S}}   \tab Strike or a strike threat       \cr
\bold{\code{EBA}} \tab \bold{E} happens before \bold{A}\cr
\bold{\code{REC}} \tab Union recognition               \cr
}                                                                    
}

\details{
The causal conditions are \bold{\code{P}}, \bold{\code{E}}, \bold{\code{A}} and \bold{\code{S}}.
All of them are binary crisp with two values: \bold{\code{0 = No}} and \bold{\code{1 = Yes}}.

The column \bold{\code{EBA}} is not a causal condition, specifying in which case the
causal condition \bold{\code{E}} happens before the causal condition \bold{\code{A}}. It has
two values (\bold{\code{0 = No}} and \bold{\code{1 = Yes}}) plus the placeholder \bold{\code{"-"}}
to signal a \dQuote{don't care}.

The outcome is the union recognition \bold{\code{EBA}}, also binary crisp with two values:
\bold{\code{0 = No}} and \bold{\code{1 = Yes}}.
}

\source{
Caren, N; Panofsky, A. (2005) \dQuote{TQCA: A Technique for Adding Temporality 
to Qualitative Comparative Analysis.} \emph{Sociological Methods & Research}
vol. 34, no.2, pp.147-172, \cr
\doi{10.1177/0049124105277197}.

Ragin, C.C.; Strand, S.I. (2008) \dQuote{Using Qualitative Comparative 
Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).} 
\emph{Sociological Methods & Research} vol.36, no.4, pp.431-441,
\doi{10.1177/0049124107313903}.
}

\keyword{datasets}
