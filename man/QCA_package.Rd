\name{About the QCA package}

\alias{QCA-package}

\docType{package}

\title{
QCA: A Package for Qualitative Comparative Analysis
}

\description{
The package \pkg{QCA} contains functions to perform Qualitative Comparative
Analysis, complemented with a graphical user interface. It implements the
comparative method as first described by Ragin (1987), and extended by Cronqvist
and Berg-Schlosser (2009) and Ragin (2000, 2008). QCA is a bridge between the
qualitative and quantitative research methodologies, making use of the qualitative
procedures in a systematic, algorithmic way (therefore increasing the
\dQuote{confidence} in the results, as understood by quantitative researchers).

The Quine-McCluskey minimization algorithms implemented in this package are
mathematically exact, as described by Dusa (2007b), Dusa (2010), Dusa and Thiem
(2015) and Dusa (2018). They all return the same, relevant set of prime implicants
for \dfn{csQCA} (binary crisp sets QCA), \dfn{mvQCA} (multi-value QCA) and
\dfn{fsQCA} (fuzzy-sets QCA).

The package also showcases functionality for other types of QCA like \dfn{tsQCA}
(temporal QCA), see Caren and Panofsky (2005), Ragin and Strand (2008) and more
recently also causal chains similar to those from the package \bold{\pkg{cna}}
(see Ambuehl et al 2015).

The results of the \bold{\pkg{QCA}} package are consistent with (and sometimes
better than) the results of the other software packages for QCA, most notably
\bold{\pkg{fs/QCA}} by Ragin and Davey (2014) and \bold{\pkg{Tosmana}} by
Cronqvist and Berg-Schlosser (2009). A comparison of several such software is
provided by Thiem and Dusa (2013).

More recent versions bring major improvements and additions, most notably: - a new
minimization algorithm called CCubes (Consistency Cubes), that is hundreds of times
faster than the previous eQMC; - conjunctural directional expectations; - extension
to categorical data.
}

\details{
    \tabular{ll}{
    Package: \tab QCA\cr
    Type:    \tab Package\cr
    Version: \tab 3.23\cr
    Date:    \tab 2024-12-03\cr
    License: \tab GPL (>= 3)\cr
    }
}

\author{
Adrian Dusa\cr
Department of Sociology\cr
University of Bucharest\cr
\email{dusa.adrian@unibuc.ro}
}

\references{
Ambuehl, M. et al (2015) \emph{A Package for Coincidence Analysis (CNA),
R Package Version 2.0 [Computer Program]},
\href{https://cran.r-project.org/package=cna/}{CRAN}.

Caren, N.; Panofsky, A. (2005) \dQuote{TQCA: A Technique for Adding
Temporality to Qualitative Comparative Analysis.}
\emph{Sociological Methods & Research} vol.34, no.2, pp.147-172.

Cronqvist, L. (2016) \emph{Tosmana: Tool for Small-N Analysis, Version 1.522
[Computer Program]}. Trier: University of Trier. url: \url{https://www.tosmana.net/}

Dusa, A. (2007a) \dQuote{User manual for the QCA(GUI) package in R}.
\emph{Journal of Business Research} vol.60, no.5, pp.576-586,
\doi{10.1016/j.jbusres.2007.01.002}

Dusa, A. (2007b) \emph{Enhancing Quine-McCluskey}.
WP 2007-49, \href{https://compasss.org/working-papers-series/}{COMPASSS Working Papers series}.

Dusa, A. (2010) \dQuote{A Mathematical Approach to the Boolean Minimization
Problem.} \emph{Quality & Quantity} vol.44, no.1, pp.99-113,
\doi{10.1007/s11135-008-9183-x}

Dusa, A.; Thiem, A. (2015) \dQuote{Enhancing the Minimization of Boolean and
Multivalue Output Functions With eQMC} \emph{Journal of Mathematical Sociology}
vol.39, no.2, pp.92-108,\cr
\doi{10.1080/0022250X.2014.897949}

Dusa, A. (2018) \dQuote{Consistency Cubes: A Fast, Efficient Method for Boolean
Minimization}, R Journal, \doi{10.32614/RJ-2018-080}

Ragin, C.C. (1987) \emph{The Comparative Method: Moving beyond Qualitative and
Quantitative Strategies}. Berkeley: University of California Press.

Ragin, C.C. (2000) \emph{Fuzzy-Set Social Science}. Chicago: University of
Chicago Press.

Ragin, C.C. (2008) \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}.
Chicago: University of Chicago Press.

Ragin, C.C.; Strand, S.I. (2008) \dQuote{Using Qualitative Comparative
Analysis to Study Causal Order: Comment on Caren and Panofsky (2005)}.
\emph{Sociological Methods & Research} vol.36, no.4, pp.431-441.

Ragin, C.C.; Davey, S. (2014) \emph{fs/QCA: Fuzzy-Set/Qualitative
Comparative Analysis, Version 2.5 [Computer Program]}.
Irvine: Department of Sociology, University of California.

Thiem, A.; Dusa, A. (2013) \dQuote{Boolean Minimization in Social
Science Research: A Review of Current Software for Qualitative Comparative
Analysis (QCA).} \emph{Social Science Computer Review} vol.31, no.4, pp.505-521.
}


\keyword{package}
