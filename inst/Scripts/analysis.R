library(motifRG)
##Use two fasta files as input
MD.motifs <- findMotifFasta(system.file("extdata", "MD.peak.fa"), system.file("extdata", "MD.control.fa"), max.motif=5,enriched=T)
###print the summary of motifs
summaryMotif(MD.motifs$motifs, MD.motifs$category)

##Load peak and control coordinates. 
data(YY1.peak)
data(YY1.control)
##Load human genome
library(BSgenome.Hsapiens.UCSC.hg19)
##Get foreground and background sequences.
YY1.peak.seq <- getSequence(YY1.peak, genome=Hsapiens)
YY1.control.seq <- getSequence(YY1.control, genome=Hsapiens)
##Search motifs.
YY1.motif.1 <- findMotifFgBg(YY1.peak.seq, YY1.control.seq, enriched=T)


##Narrow the YY1 peaks.
YY1.narrow.seq <- subseq(YY1.seq, pmax(round((width(YY1.seq) - 200)/2), 1), width=pmin(200, width(YY1.seq)))
YY1.control.narrow.seq <- subseq(YY1.control.seq, pmax(round((width(YY1.control.seq) - 200)/2), 1), width=pmin(200, width(YY1.control.seq)))

##concatenate sequences. 
all.seq <- append(YY1.narrow.seq, YY1.control.narrow.seq)
##calcualte GC content, and discretize it. 
gc <- as.integer(cut(letterFrequency(all.seq, "CG", as.prob=T),c(-1, 0.4, 0.45, 0.5, 0.55, 0.6, 2)))
##sequences category: 1 for foreground, 0 for background. 
category=c(rep(1, length(YY1.narrow.seq)), rep(0, length(YY1.control.narrow.seq)))
##Search motifs using the GC content as covariant for the regression model to adjust GC bias.
YY1.motif.2 <- findMotif(all.seq,category, other.data=gc, max.motif=5,enriched=T)

###Weight foreground sequences according to peak intensity and perform weighted regression.
all.weights = c(peak$weight, rep(1, length(YY1.control.seq)))
YY1.motif.3 <- findMotif(all.seq,category, other.data=gc, max.motif=5,enriched=T, weights=all.weights)

###load CTCF motifs predicted by motifRG###
data(ctcf.motifs)
summaryMotif(ctcf.motifs$motifs, ctcf.motifs$category)

###print sequence logo of the first motif
plotMotif(ctcf.motifs$motifs[[1]]@match$pattern)

###Find a refined PWM model given the motif matches as seed
ctcf.seq <- readDNAStringSet(system.file("inst/extdata", "ctcf.fa"))
pwm.match <- refinePWMMotif(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
library(seqLogo)
seqLogo(pwm.match$model$prob)


### Motifs found by findMotif tend to be relatively short, as longer and more specific
### motif models do not necessarily provide better discrimination of foreground background 
### vs background if they are already well separated.
### However, one can refine and extend a PWM model given the motif matches by findMotif as seed
### for more specific model. 
pwm.match.extend <- refinePWMMotifExtend(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
seqLogo(pwm.match.extend$model$prob)
###plot the dinucleotide logo fo PWM match
plotMotif(pwm.match.extend$match$pattern)

###Show dependency of adjacent positions, green for enriched pair, red for depleted pair
plotMotif(pwm.match.extend$match$pattern, logodds=T)







