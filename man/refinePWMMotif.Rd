% File motifRG/man/refinePWMMotif.Rd
\name{refinePWMMotif}
\alias{refinePWMMotif}
\title{create a PWM (Position Weight Model) model given a initial set of motif matches and
  input sequences}
\description{
  Create a PWM model given a initial set of motif matches and
  input sequences
}
\usage{
  refinePWMMotif(motifs=NULL, seqs,  pwm.ld= NULL, max.iter=50,
tol=10^-4, mod="oops", null=rep(0.25, 4),pseudo=1, weights=rep(1,
length(seqs)), motif.weights=NULL)
}
\arguments{
  \item{motifs}{The initial set of motif matches. character vector or
    DNAStringSet object}
  \item{seqs}{Input sequences. character vector or DNAStringSet
    object}
  \item{pwm.ld}{The initial PWM matrixes in logodds
    transformation. Either "motifs" or "pwm.ld" is not NULL}
  \item{max.iter}{Maximum number of iterations for refinement}
  \item{tol}{Convergence criteria. The percentage of total PWM scores
    improvement required for convergence.}
  \item{mod}{Motif ocurrence model.
    If \code{mod=="oops"}, assume one motif match per sequence.
    If \code{mod=="zoops"}, assume zero or one motif match per sequence.    
  }
  \item{null}{A numeric vector specifying the background model}
  \item{pseudo}{Pseudo counts for PWM construction}
  \item{weights}{a numeric vector specifying the weights for all
    sequences. Default: 1 for all sequences}
  \item{motif.weights}{a numeric vector specifying the weights for
    initial sets of motifs. Default: NULL} 
}
\value{
  Return a list with two elements:
  \item{model}{a list with two elements. "prob": PWM model, sum of columns
    add to 1. "logodd": PWM model in logodds form, log2 of original
    matrix substract the background model}.
  \item{match}{a data.frame specifying the motif matches in each
    sequence. Columns are:
    "match": the sequence of the match, "score": PWM score, "strand",
    the strand of the match in the input sequence, "pos": start position of the
    motif match. If multiple matches are allowed, then "seq.id"
    specifies the index of the input sequence for the motif match.}
  \item{score}{Total PWM score of the motif matches}
}

\seealso{
  \code{\link{findMotif}}
  \code{\link{refinePWMMotifExtend}}
}

\examples{
ctcf.seq <- readDNAStringSet(system.file("extdata","ctcf.fa", package="motifRG"))
data(ctcf.motifs)
### refine PWM model based on motif matches
pwm.match <- refinePWMMotif(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
### plot traditional motif logo
library("seqLogo")
seqLogo(pwm.match$model$prob)
### plot dinucleotide motif logo
plotMotif(pwm.match$match$pattern)
### automatically extend PWM model 
pwm.match.extend <-refinePWMMotifExtend(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
### plot the new motif matches
plotMotif(pwm.match.extend$match$pattern)
}


