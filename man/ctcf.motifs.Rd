\name{ctcf.motifs}
\alias{ctcf.motifs}
\title{CTCF motifs predicted by motifRG}
\description{
  The output produced by \link{findMotif}.
}

\details{
  A list with following elements:
  motifs:a list motif descriptions of class \link{Motif-class}.
  category:input binary specification of foreground/background.
  mask.motifs:if \code{mask=T}, then mask.motifs contain the
  description of motif is based on motif matches after the input
  sequences being masked by previous motifs. In this case, "motifs"
  contained the unmasked motif descriptions.
}
\references{
  Unpublished
}
\keyword{datasets}
