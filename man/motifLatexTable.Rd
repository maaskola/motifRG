% File motifRG/man/plotMotif.Rd
\name{motifLatexTable}
\alias{motifLatexTable}
\alias{motifHtmlTable}
\title{create of table of motifs}
\description{create a latex table to be embedded in a latex document}
\usage{
 motifLatexTable(motifs,  main="", prefix="motif", dir=".", height=1,
 width=3,enriched.only=F, plot.pwm= TRUE,
summary.cols=c(1,7,8,9),use.mask=TRUE)
 motifHtmlTable(motifs, dir="html", prefix="motif", enriched.only=F,
plot.pwm= TRUE, summary.cols=c(1,7,8,9),use.mask=TRUE)
}
\arguments{
  \item{motifs}{result of findMotif}
  \item{main}{The title of table}
  \item{prefix}{The prefix for the filenames of motif logos}
  \item{dir}{The directory for storing motif logo files}
  \item{height, width}{size of the sequence logo}
  \item{enriched.only}{If true, list only enriched motifs}
  \item{plot.pwm}{If true, plot PWM logo instead of di-nucleotide logo}  
  \item{summary.cols}{The selected columns of summary table created by
    summaryMotif included in the table}
  \item{use.mask}{If true, use masked motif match summary statistics}
}
\value{
  motifLatexTable outputs a latex table to the stdout console.
  motifHtmlTable outputs a html file named as <preix>.html in "dir"
  directory. 
}

\seealso{
  \code{\link{findMotif}}
}

\examples{
  data(ctcf.motifs)
  ### Create table of motifs in Latex 
  motifLatexTable(ctcf.motifs, main="CTCF motifs", dir="motif")

  ### Create table of motifs in Html
  motifHtmlTable(ctcf.motifs, dir="Html")
}
