\name{findMotifFgBg}
\alias{findMotifFgBg}
\title{Wrapper function for findMotif using two sequence datasets}
\description{
  Perform motif search on two sequence datasets. This is a wrapper
  function for findMotif.
}
\usage{
findMotifFgBg(fg.seq, bg.seq, ...)
}
\arguments{
  \item{fg.seq}{DNAStringSet; foreground sequences}
  \item{bg.seq}{DNAStringSet; background sequences}
  \item{...}{Other parameters passed to findMotif}
}