% File motifRG/man/summaryMotif.Rd
\name{summaryMotif}
\alias{summaryMotif}
\title{summarize a list of motifs}
\description{
  Create a summary table of a list of motifs found by findMotif}
\usage{
summaryMotif(motifs, category)
}
\arguments{
  \item{motifs}{a list of motifs of class \link{Motif-class}}  
  \item{category}{a vector of 0 or 1, specifying which sequences are
    foreground and background. Input for findMotif}
}
\value{
  A data.frame with following columns:
  \item{scores}{scores for each Motif. All values are negative.
    The absolute scales of the scores reflect the discriminative power
    of the motif for separating the foreground and
    background. Statistically, they correspond to the Z-values of the
    predictor(counts of the motifs in this case) in the logistic
    regression model}
  \item{signs}{sign of the motifs. TRUE for enriched motifs, FALSE for
    depleted motifs}
  \item{fg.hits, bg.hits}{Total number of hits in the foreground, and
    background sequences. If the motif is scanned on both strands of the
    input sequences, the counts on both strands are added.}
  \item{fg.seq, bg.seq}{The number of sequences that contain at least one motif
    match in the foreground, and the background}.
  \item{ratio}{The enrichment/depleted ratio of motifs}
  \item{fg.frac,bg.frac}{The fraction of foreground/background sequences
    that contain at least one motif match}
}
\seealso{
  \code{\link{findMotif}}
  \code{\link{motifLatexTable}}
  \code{\link{motifHtmlTable}}
}

\examples{
data(ctcf.motifs)
###plot the summary statics of motif matches after masking previous motif ocurrences###
summaryMotif(ctcf.motifs$mask.motifs, ctcf.motifs$category)

###plot the summary statics of motif matches in the original sequences###
summaryMotif(ctcf.motifs$motifs, ctcf.motifs$category)
}


