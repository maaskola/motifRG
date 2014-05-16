\name{Motif-class}
\docType{class}
%Classes
\alias{class:Motif}
\alias{Motif-class}

\title{Motif objects}
\description{
  A Motif object contains general motif characteristics and and details of motif match
}

\details{
  A motif object has the following slots:
  score:absolute z-value based on the logistic regression model for the
    motif.
  sign:the sign of the motif:plus for enriched motif in the
    foreground sequences, and negative for depleted motif
  count:a numeric vector holding the number of matches in each
    sequence
  match:a data.frame with the following columns:
    match.strand:the strand on which the match is found;
    pattern:the motif match pattern;
    seq.id: on which sequence the match is found;
    pos:the position relative to sequence start of the match. 
  pattern:the motif pattern
  consensus:the motif consensus pattern determined by the
  majority votes at each position using the following rule: the most
  dominate single nucleotide if its frequency is greater than 0.6, or
  the two most dominate nucleotide if combined frequency is greater than 0.8, or 
  the three most dominate nucleotide if combined frequency is greater
  than 0.95
}

\seealso{
  \code{\link{findMotif}}
  \code{\link{summaryMotif}}
  \code{\link{plotMotif}}
  \code{\link{motifLatexTable}}
}
