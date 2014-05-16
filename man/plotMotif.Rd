% File motifRG/man/plotMotif.Rd
\name{plotMotif}
\alias{plotMotif}
\title{plot motif sequence matches}
\description{
  plot aligned sequences, revealing the independent position specificity and
  dependency among adjacent positions.}
\usage{
plotMotif(match, logodds=F, entropy=F, bg.ld=NULL,  alphabet=c("A", "C", "G", "T"),has.box=T,...)
}
\arguments{
  \item{match}{motif match to be plotted. character or \link{DNAStringSet} object}.
  \item{logodds}{logical; if true, plot the enrichment/depletion of a
    adjacent pair relative to the independent model.}
  \item{entropy}{logical; if true, areas outside the core region of the
    motif are dimmed}
  \item{bg.ld}{Experimental features: background dinucleotide logodds against independent model. if
    \code{logodds=T}, then the background logodds will be substracted}    
  \item{alphabet}{the alphabets used in the sequence. Do not change its value}
  \item{has.box}{logical; if true, plot the boundaries of the motif}
  \item{...}{other arguments passed to the lower level plot function}
}

\details{
  X-axis refers to the positions of the motifs.  
  
  Y-axis correspond to the alphabets.
  
  Letter sizes define the frequencies of the nucleotides at a given position.
  
  Edges between the letters specifiy the dinucleotide
  relationship. The depth of the color correspond to the dinucleotide
  frequency. If \code{logodds=T}, thinner edges will be plotted between
  dependent pairs. The edge is colored red if the pair is depleted
  (relative to the expected frequency if the pair is independent), and
  green if the pair is enriched. The gradient of color red/green
  correspond to the level of dependency.         
}

\examples{
data(ctcf.motifs)
### plot the dinucleotide representation of the first motif
plotMotif(ctcf.motifs$motifs[[1]]@match$pattern)
plotMotif(ctcf.motifs$motifs[[1]]@match$pattern, logodds=TRUE)
plotMotif(ctcf.motifs$motifs[[1]]@match$pattern, logodds=TRUE, entropy=TRUE)
}


