\name{QCA internal functions}

\alias{.onAttach}

\alias{as.panel}
\alias{combint}
\alias{minimizeLoop}
\alias{print.aE}
\alias{print.chain}
\alias{print.fuzzy}
\alias{print.modelFit}
\alias{print.mqca}
\alias{print.panel}
\alias{print.pic}
\alias{print.pof}
\alias{print.qca}
\alias{print.sS}
\alias{print.tt}
\alias{rebuild}
\alias{removeRedundants}
\alias{replaceCategories}
\alias{verify.data}
\alias{verify.dir.exp}
\alias{verify.minimize}
\alias{verify.inf.test}
\alias{verify.qca}
\alias{verify.mqca}
\alias{verify.multivalue}
\alias{verify.tt}

\alias{d.AS}
\alias{d.autonomy}
\alias{d.Bas}
\alias{d.biodiversity}
\alias{d.BWB}
\alias{d.CS}
\alias{d.CZH}
\alias{d.educate}
\alias{d.education}
\alias{d.Emm}
\alias{d.graduate}
\alias{d.health}
\alias{d.HK}
\alias{d.HMN}
\alias{d.homeless}
\alias{d.jobsecurity}
\alias{d.Kil}
\alias{d.Kro}
\alias{d.napoleon}
\alias{d.pban}
\alias{d.partybans}
\alias{d.represent}
\alias{d.RS}
\alias{d.SA}
\alias{d.socialsecurity}
\alias{d.SS}
\alias{d.stakeholder}
\alias{d.transport}
\alias{d.urban}
\alias{d.women}
\alias{Emme}
\alias{HarKem}
\alias{Krook}
\alias{Quine}
\alias{RagStr}
\alias{Rokkan}


\alias{nec}
\alias{suf}


\alias{getSolution}
\alias{writeSolution}
\alias{prettyString}
\alias{rowDominance}
\alias{sortMatrix}
\alias{sortVector}
\alias{mvregexp}


\title{Internal Functions}

\description{
  The above functions are internal in the QCA package which are not designed to
  be called directly by the user. All of them are used by the \code{minimize} function,
  except \code{sortMatrix} which is used by \code{allExpressions}. The verification
  and error messages have been moved to the internal functions \code{verify.data}
  and \code{verify.tt}.
}

\keyword{internal}
