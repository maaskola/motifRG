% File motifRG/man/refinePWMMotifExtend.Rd
\name{refinePWMMotifExtend}
\alias{refinePWMMotifExtend}
\title{create an extended PWM (Position Weight Model) model given a initial set of motif matches and
  input sequences}
\description{
  Create an extended PWM model given a initial set of motif matches and
  input sequences
}
\usage{
  refinePWMMotifExtend(motifs=NULL, seqs, pwm.ld=NULL, flank=3,  extend.tol=10^-3, trim.rel.entropy=0.2,  null=rep(0.25, 4), max.width=20, ...)
}
\arguments{
    \item{motifs}{The initial set of motif matches. character vector or
    DNAStringSet object}
  \item{seqs}{Input sequences. character vector or DNAStringSet
    object}
  \item{pwm.ld}{The initial PWM matrixes in logodds
    transformation. Either "motifs" or "pwm.ld" is not NULL}
  \item{flank}{The number of bases for extension on both sides of the
    motif. The extension will be iterated if the there is sufficient
    signal in the flanking region.}  
  \item{extend.tol}{Convergence criteria for extension.}  
  \item{trim.rel.entropy}{cutoff to be used to trim the uninformative flanking of a PWM
    model based on relative entropy against a null distribution.}
  \item{null}{NULL background distribution}
  \item{max.width}{The maximum width of PWM}
  \item{...}{other arguments passed to function refinePWMMotif}
}
\value{Same type of object returned by \link{refinePWMMotif}}
\details{Flanking regions with length equal to flank is still included in output for reference}

\seealso{
  \code{\link{findMotif}}
  \code{\link{refinePWMMotif}}
}
