\name{findMotif}
\alias{findMotif}
\title{De-novo discovery of distriminative motifs}
\description{
  The function searches motifs that discriminate the given foregound and
  background sequences. 
}
\usage{
findMotif(all.seq,  category, weights = rep(1, length(all.seq)),
start.width=6,min.cutoff=5, min.ratio=1.3,
min.frac=0.01, both.strand=TRUE, flank=2, max.motif=5,
mask=TRUE,other.data=NULL, start.nmer=NULL,
enriched.only=F,n.bootstrap = 5, bootstrap.pvalue=0.1,is.parallel =
TRUE,mc.cores = 4,min.info=10,max.width=15,discretize=TRUE)
}

\arguments{
  \item{all.seq}{ DNAStringSet; foreground and background sequences.}
  \item{category}{numeric vector; specify which sequences are
    foreground (with value 1), and background (value 0).}
  \item{weights}{numeric vector: the weights for all sequences. Default: 1}
  \item{start.width}{logical; the width for enumerating seed patterns}
  \item{min.cutoff}{numeric; the score cutoff required for seed selection. All
    scores are negative, the lower the better.}
  \item{min.ratio}{numeric; the minimum fold change of motif occurences
    in foreground vs background.}
  \item{min.frac}{numeric; the minimum fraction of fg/bg sequences
    containing the candidate motifs}
  \item{both.strand}{logical; if true, search both strands}
  \item{flank}{integer; the length for step-wise pattern extension at
    both ends on candidate motifs}
  \item{max.motif}{integer; the maximum number of output motifs}
  \item{mask}{logical; if true, mask previous motifs when searching for
    the next motif}
  \item{other.data}{if not NULL, a matrix with additional terms for the regression
    model for bias adjustment}
  \item{start.nmer}{if not NULL, a matrix with counts for user specified seed pattern in
    each sequence}
  \item{enriched.only}{logical; if true, only predict enriched motif}
  \item{n.bootstrap}{integer; the number of bootstrapping tests to
    estimate score variance}
  \item{bootstrap.pvalue}{numeric: the bootstrap t.test pvalues to
    determine the significance of improvement}
  \item{is.parallel}{logical;if true, runs in parallel mode, and requires
    "parallel" library }
  \item{mc.cores}{integer; the number of CPUs for paralel run}
  \item{min.info}{minimal information content for the motif to prevent
    it from being too degenerate}
  \item{max.width}{maximum width of the motif for extension}
  \item{discretize}{\code{logical} default TRUE}
}

\value{
  return a list with following elements:
  \item{motifs}{a list motif descriptions of class \link{Motif-class}}. 
  \item{category}{input binary specification of foreground/background}
  \item{mask.motifs}{if \code{mask=T}, then mask.motifs contain the
  description of motif is based on motif matches after the input
  sequences being masked by previous motifs. In this case, "motifs"
  contained the unmasked motif descriptions.}
}


\examples{
MD.peak.seq <- readDNAStringSet(system.file("extdata","MD.peak.fa", package="motifRG"))
MD.control.seq <- readDNAStringSet(system.file("extdata","MD.control.fa", package="motifRG"))
category <- c(rep(1, length(MD.peak.seq)), rep(0, length(MD.control.seq)))
MD.motifs <- findMotif(append(MD.peak.seq, MD.control.seq),category, max.motif=3,enriched=TRUE)

### Get summary of motifs
summaryMotif(MD.motifs$motifs, MD.motifs$category)

### plot the dinucleotide representation of the first motif
plotMotif(MD.motifs$motifs[[1]]@match$pattern)

### Create table of motifs in Latex 
motifLatexTable(MD.motifs, main="MD motifs")

### Create table of motifs in Html
motifHtmlTable(MD.motifs)
}


