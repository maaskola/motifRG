library(IRanges)
library(Biostrings)
logPWM <- function(pwm,  null=rep(0.25, 4))
{
  return(log2(pwm/null))
}

seqToMatrix <- function(seqs, encode= NULL)
  {
    if(!is.null(encode)){    
      width <- nchar(seqs[1])
      seq <- paste(seqs, collapse="")      
      seq.num <- as.integer(toComplex(DNAString(seq), encode))
      matrix(seq.num, nrow=length(seqs), ncol=width, byrow=T)        
    }
    else{
      as.matrix(seqs)
    }
  }

getPWM <- function(match, weights= NULL, pseudo=10, alphabet=DNA_BASES, null=rep(1/length(alphabet), length(alphabet)))
  {
    if(is.null(weights)){
      pwm <-  (consensusMatrix(match) + pseudo)[alphabet,]
      pwm <- t(t(pwm)/colSums(pwm))
    }
    else{
      mat <- seqToMatrix(DNAStringSet(match))
      pwm <- apply(mat, 2, function(x){
        tmp<- tapply(weights, x, sum)
        counts <- rep(0, length(alphabet))
        names(counts) <- alphabet
        tmp <- tmp[names(tmp)%in% alphabet]
        counts[names(tmp)] <- tmp
        counts
      })
      pwm <- t(t(pwm+ pseudo) / colSums(pwm+pseudo))
    }
    
    log.pwm <- logPWM(pwm, null)
    return(list(prob=pwm, logodd=log.pwm, bg=null))
  }

bestPWMMatch <- function(pwm, seqs, both.strand=T, flank=0) #pwm in logodd
  {
    seqs.view <- pasteSeq(seqs)
    seq <- subject(seqs.view)
    width <- ncol(pwm)
    len <- length(seq)
    tmp <- PWMscoreStartingAt(pwm, seq, 1:(len-width+1))
    strand <- rep(T, length(seqs))
    if(both.strand){
      seq.rev <- reverseComplement(seq)
      tmp2 <- rev(PWMscoreStartingAt(pwm, seq.rev, 1:(len-width+1)))
      tmp <- pmax(tmp, tmp2)    
      strand <- tmp > tmp2
    }
    tmp.view <- Views(Rle(tmp), start(seqs.view), end(seqs.view) - width+1)
    match.pos <- viewWhichMaxs(tmp.view)
    match <- DNAStringSet(Views(seq, match.pos - flank, match.pos+width-1 + flank))
    if(sum(!strand) > 0){     
      match.rev <- reverseComplement(match[!strand[match.pos]])
      match <- as.character(match)
      match[!strand[match.pos]] <- as.character(match.rev)
    }
    else{
      match <- as.character(match)
    }
    return(data.frame(pattern=match, score= tmp[match.pos],strand=strand[match.pos],pos = match.pos - start(seqs.view), stringsAsFactors=F))           
  }                      


scanPWMMatch <- function(pwm, seqs,  both.strand=T, flank=0,...) ##logscale
  {
    seqs.view <- pasteSeq(seqs)
    seq <- subject(seqs.view)
    width <- ncol(pwm)
    len <- length(seq)   
    match <- matchPWM(pwm, seq, ...)
    match <- IRanges(start(match), end(match))
    strand <- rep(T, length(match))
    if(both.strand){
      match.rev <- matchPWM(reverseComplement(pwm), seq, ...)
      match.rev <- IRanges(start(match.rev), end(match.rev))
      strand <- c(strand, rep(F, length(match.rev)))
      match <- c(match, match.rev)
    }    
    ## dup <- duplicated(start(match))
    ## match <- match[!dup]
    ## strand <- strand[!dup]    
    tmp <- Views(seq, start(match) - flank, end(match) + flank)   
    pattern <- as.character(tmp)
    pattern[!strand] <- as.character(reverseComplement(tmp[!strand]))
    seq.id <- findOverlaps(match, seqs.view, select="first")
    match.score <- bestPWMMatch(pwm, pattern, both.strand=F)$score
    df <- data.frame(pattern=pattern, score=match.score, seq.id=seq.id, strand=strand, pos=start(match) - start(seqs.view)[seq.id], stringsAsFactors=F)
    row.names(df) <- NULL
    select <- tapply(1:nrow(df), start(match), function(x){x[which.max(df$score[x])]})
    df <- df[select,]
    ord <- order(df$seq.id, df$pos)
    df[ord,]
  }

              

split.pwm <- function(matches,null=rep(0.25,4), min.split=0.5, min.size=100)
  {
    all.splits <- to.split <- list(all=1:length(matches))
    all.pwms <- list(all=getPWM(matches)$prob)
    while(length(to.split) > 0){      
      select <- to.split[[1]]
      name <- names(to.split)[[1]]
      pwm <- all.pwms[[name]]
      e <- sum(apply(pwm, 2, relativeEntropy, null=null))
      to.split <- to.split[-1]      
      to.split.scores <- sapply(1:ncol(pwm), function(i){
        sapply(names(IUPAC), function(a){
          tmp <- substring(matches[select], i,i) %in% IUPAC[[a]]
          total <- sum(tmp)
          if(total > min.size & total < length(select)- min.size){
            pwm1 <- getPWM(matches[select][tmp])$prob
            pwm2 <- getPWM(matches[select][!tmp])$prob
            w <- sum(tmp)/length(tmp)
            (sum(w*apply(pwm1, 2, relativeEntropy, null=null) +
                (1-w)*apply(pwm2, 2, relativeEntropy, null=null)) - e) * length(select)/length(matches)
          }
          else{
            0
          }
        })})
      if(max(to.split.scores) > min.split){
        a <- which.max(apply(to.split.scores, 1, max))
        i <- which.max(apply(to.split.scores, 2, max))
        desc <- paste(name, names(a), i)
        cat(desc, "\n")
        tmp <- substring(matches[select], i,i) %in% IUPAC[[a]]
        to.split <- c(list(select[tmp], select[!tmp]), to.split)
        names(to.split)[1:2] <- paste(desc, c(T,F))
        all.pwms <- c(list(getPWM(matches[select][tmp])$prob, getPWM(matches[select][!tmp])$prob), all.pwms)
        names(all.pwms)[1:2] <- paste(desc, c(T,F))
        all.splits <- c(to.split[1:2], all.splits)
      }
    }
    return(list(splits=all.splits, pwms = all.pwms))
  }



matchPWMCoor <- function(pwm, seqs, coor,  both.strand=T,...)
{
  seqs <- pasteSeq(seqs)
  seq <- subject(seqs)    
  match <- matchPWM(pwm, seq, ...)
  match.df <- data.frame(strand=rep("+", length(match)), pattern = as.character(match),stringsAsFactors=F)
  if(both.strand){    
    pwm.rev <- reverseComplement(pwm)
    match.rev <- matchPWM(pwm.rev, seq,...)
    match.df <- rbind(match.df, data.frame(strand=rep("-", length(match.rev)),
                                           pattern = as.character(reverseComplement(match.rev)), stringsAsFactors=F))
    match <- c(as(match,"IRanges"), as(match.rev, "IRanges"))                     
  }
  ord <- order(start(match))
  dup <- duplicated(start(match))
  ord <- ord[!dup[ord]]
  match<- match[ord]
  match.df <- match.df[ord,]  
  l <- findOverlaps(match, seqs, select="first")
  r <- IRanges(start(match) - start(seqs)[l] + start(coor)[l],end(match) - start(seqs)[l] + start(coor)[l])
  seqnames <- seqnames(coor)[l]
  score <- bestPWMMatch(pwm, match.df$pattern, both.strand=F)$score
  GRanges(r, seqnames=seqnames, strand=match.df$strand, pattern=match.df$pattern, score=score, seq.id=l)
}


bestPWMMatchCoor <- function(pwm, seqs, coor, both.strand=T, flank=0)
  {
    w <- ncol(pwm)
    match <- bestPWMMatch(pwm, seqs, both.strand=both.strand, flank=flank)
    seq.strand <- rep("+", length(seqs))
    strand <- rep("+", length(seqs))
    if("strand" %in% colnames(coor)){
      seq.strand <- coor$strand
    }
    start <- start(coor) + match$pos
    end <- start + w  - 1
    strand <- rep("+", nrow(match))
    
    select <- seq.strand== "-"
    end[select] <- (end(coor) - match$pos)[select]
    start[select] <- end[select] - w + 1
    
    select <- xor(seq.strand == "+",  match$strand)
    strand[select] <- "-"
    rd <- RangedData(IRanges(start, end), strand=strand, seq.strand=seq.strand, match.strand=match$strand, pattern=match$pattern,
                     space=space(coor), score=match$score, pos=match$pos)
    rd[names(coor)]
  }




refinePWMMotif <- function(motifs=NULL, seqs, pwm.ld=NULL,  max.iter=50, tol=10^-4, mod="oops", null=rep(0.25, 4),pseudo=1, weights=rep(1, length(seqs)), motif.weights=NULL)
  {
    
    total.score <- -Inf
    diff <- Inf
    iter = 1    
    if(mod=="zoops"){
      gamma0 <- gamma <- 0.5
    }
    if(is.null(pwm.ld)){
      patterns <- motifs
      if(is.null(motif.weights)){
        motif.weights=rep(1, length(patterns))
      }
      pwm.ld <- (getPWM(patterns, weights=motif.weights, null=null,pseudo=pseudo))$logodd
    }
    repeat{
      if(mod=="oops"){
        match<- bestPWMMatch(pwm.ld, seqs)
        match$weights = weights
      }
      else{
        match <- scanPWMMatch(pwm.ld, seqs, min.score=1)
        ratio <- 2^match$score
        tmp <- tapply(ratio, match$seq.id, sum)
        seq.sum <- rep(0, nrow(match))
        seq.sum[as.numeric(names(tmp))] <- tmp
        if(mod=="zoops"){
          lambda <- gamma/ width(seqs)[match$seq.id]          
          match$weights <- ratio * lambda/ (1 - gamma + seq.sum[match$seq.id]*lambda) * weights[match$seq.id]
          gamma <- (sum(match$weights) + gamma0)/(sum(weights)+1)
          #cat("gamma", gamma, "\n")          
        }
      }      
      patterns <- match$pattern
      new.score <- sum(match$score * match$weights)
      pwm <- (getPWM(patterns,weights=match$weights, null=null,pseudo=pseudo))
      pwm.ld <- pwm$logodd     
      cat(iter, new.score, "\n")
      iter = iter +1
      if((new.score - total.score)/abs(new.score) < tol || iter > max.iter ){
        break
      }
      total.score = new.score
    }
    return(list(model=pwm, match=match, score=total.score))
  }



refinePWMMotifExtend <- function(motifs=NULL, seqs, pwm.ld=NULL, flank=3, extend.tol=10^-3, trim.rel.entropy=0.2,  null=rep(0.25, 4), max.width=20, ...)
  {
    prev.score <- 0    
    repeat{
      result <- refinePWMMotif(motifs=motifs, seqs=seqs, pwm.ld=pwm.ld, null=null, ...)
      pwm <- result$model
      width <- ncol(pwm$prob)
      ###trim on entropy###
      e <- apply(pwm$prob, 2, relativeEntropy, null=null)
      trim.left <- 1
      trim.right <- width      
      if(!any(e > trim.rel.entropy)) {break}      
      tmp <- which(e > trim.rel.entropy)
      trim.left <- min(tmp)
      trim.right <- max(tmp)        
      if(trim.left > flank & width - trim.right > flank){ ## no information in the flanking region, no need to extend
        break
      }
      if(trim.right - trim.left > max.width){
        break
      }
      match <- bestPWMMatch(pwm$logodd[,trim.left:trim.right], seqs, flank=flank)      
      motifs <- match$pattern
      new.score <- result$score
      cat("Extend score", new.score, prev.score, "\n")
      if((new.score - prev.score)/abs(new.score) < extend.tol){break}
      prev.score = new.score
      pwm.ld <- NULL
    }
    return(result)
  }

