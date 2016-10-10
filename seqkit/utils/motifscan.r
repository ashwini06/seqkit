######################### MOTIF SCAN ###################################

#"~/Documents/Sequence_Analysis/Foxa2_Chip/Foxa2motifanalysis/wt.d6.overlap/Foxa2_PWM_toplogo1.txt
#"/media/ashwini/Ashwini/PROCESSED_DATA/ChIP-Seq/lizzy_andersson_ChIP_ESC/Chip_Pip_files/Foxa2_PWM_toplogo1.txt"

PWMfile = args[1]
Foxa2.pwm <- read.table(paste0(PWMfile))
Foxa2.pwm1 <- t(Foxa2.pwm)
rownames(Foxa2.pwm1) <- c("A","C","G","T")

seqLogo(Foxa2.pwm1,ic.scale=T)



### obtain the sequences nder the peaks  and prepare two heper function that will be called on each peak set. Those functions will return a matrix containing the chr
### and start position for each motif hit,as well as the sequence of each motif instance.


elementMetadata(macspeaks)$seqs <- getSeq(Mmusculus, macspeaks, as.character=TRUE)

# function to scan a string for PWM matches at the specified threshold by calling matchPWM()
# keeps only the closest of multiple hits
mymatchPWM <- function (pwm, myseq, threshold, summit) {
  # get all matches of PWM
  mymatch <- matchPWM(pwm, myseq, min.score=threshold)
  # collect starts/seqs into matrix (if any)
  if (length(mymatch)==0) {
    found <- cbind(NA,NA,0)
  } else {
    found <- cbind(start(mymatch), as.character(mymatch), length(start(mymatch)))
  }
  colnames(found) <- c("start","seq","nr")
  # keep only the match that is closest to the summit
  found <- found[order(abs(as.integer(found[,1])-summit),decreasing=FALSE)[1],]
  # return matrix
  return(found)
}


# function to call mymatchPWM() on each peak in a set
# returns a matrix with motif information per peak
ScanPeaks <- function(peak.GR, pwm, threshold) {
  # get all peak sequences
  myseqs <- elementMetadata(peak.GR)$seqs
  # get summit positions (relative to peak coordinates)
  summits <- elementMetadata(peak.GR)$maxpos - start(peak.GR)
  # apply mymatchPWM() to all peaks in the set
  motifmatrix <- sapply(1:length(myseqs), function(x) mymatchPWM(pwm, myseqs[x], threshold, summits[x]))
  motifmatrix <- t(motifmatrix)
  # set peak IDs as rownames
  rownames(motifmatrix) <- names(peak.GR)
  # return motif matrix
  return(motifmatrix)
}



# scan all peak sequences for the motif (using a cutoff of 80%)
macs.motifs <- ScanPeaks(macspeaks, as.matrix(Foxa2.pwm1), "70%")


elementMetadata(macspeaks)$motifs <- ScanPeaks(macspeaks, as.matrix(Foxa2.pwm1), "70%")

# add motif information (number of motifs per peaks) to the peak GRanges object
elementMetadata(macspeaks)$motif.no <- as.integer(macs.motifs[,3])


# # plot the proportion of peaks with 0, 1 or more motifs in the three peak sets
# 
# pie( c(sum(as.numeric(elementMetadata(macspeaks)$motif.no==0)),
#        sum(as.numeric(elementMetadata(macspeaks)$motif.no==1)),
#        sum(as.numeric(elementMetadata(macspeaks)$motif.no>1))),
#      labels=c("0","1","2+"), main=list("No of motifs:Foxa2_PWM1-70% threshold",cex=1.5),
#      col=c("white","lightgrey","darkgrey")
# )

# set the extension distance
mydist <- 10
# get absolute start positions of the motif hits
abs.starts <- as.integer(macs.motifs[,1]) + start(macspeaks)
abs.starts <- abs.starts[!is.na(abs.starts)]
# make GRanges
motifs.GR <- GRanges(
  seqnames=seqnames(macspeaks)[!is.na(macs.motifs[,1])],
  range=IRanges(start=abs.starts - mydist,end=abs.starts + mydist + 9),
  strand="+"
)



# get sequences
motif.seqs <- getSeq(Mmusculus, motifs.GR, as.character=TRUE)
# get letter counts per position
bs.matrix <- consensusMatrix(motif.seqs)
# transform counts to frequencies
bs.matrix <- apply(bs.matrix, 2, function(x){ x/sum(x) })
# plot the newly obtained sequence logo
seqLogo(bs.matrix, ic.scale=TRUE)


## Motif Localisation

# set window size
mydist <- 100

# function to determine the motif profile around peak summits
getProfile <- function(peaks.GR, pwm, window.size){
  #   get regions around summit
  summits.GR <- GRanges(
    seqnames=seqnames(peaks.GR),
    range=IRanges( start=elementMetadata(peaks.GR)$maxpos - window.size,
                   end=elementMetadata(peaks.GR)$maxpos + window.size),
    strand="+"
  )
  # create unique names for all peaks
  names(summits.GR) <- paste( seqnames(summits.GR), start(summits.GR), end(summits.GR), sep=":")
  # get sequence
  elementMetadata(summits.GR)$seqs <- getSeq(Mmusculus, summits.GR, as.character=TRUE)
  # scan sequences with the PWM
  summit.motifs <- ScanPeaks(summits.GR, pwm, "80%")
  summit.motifs <- summit.motifs[!is.na(summit.motifs[,1]),]
  # get all covered positions
  motif.pos <- sapply(as.integer(summit.motifs[,1]), function(x) seq(x, x+ncol(pwm)-1))
  motif.pos <- table(unlist(motif.pos))
  # convert to data.frame
  motif.pos <- data.frame(motif.pos)
  names(motif.pos) <- c("position","frequency")
  # shift positions relative to summit
  motif.pos$position <- as.integer(as.character(motif.pos$position)) - window.size
  # ensure all positions are present in the output
  profile <- data.frame(
    position=-window.size:window.size,
    frequency=motif.pos$frequency[match(-window.size:window.size, motif.pos$position)]
  )
  profile[is.na(profile)] <- 0
  return(profile)
}
# get the profiles
macs.profile <- getProfile(macspeaks, as.matrix(Foxa2.pwm1), mydist)

# generate plots
pl1 <- xyplot(frequency~position, data=macs.profile,
              type="l", main="motifs around MACS v2 peak summits",
              aspect=0.8
)

print(pl1, split=c(1,1,2,2), more=TRUE)

print(pl1)

outfile = args[2]
write.table(as.data.frame(macspeaks),paste(outfile),sep="\t",quote=F,row.names=F)



# # To obtain all the match of motifs in the target sequences
# 
# hits =matchPWM(Foxa2.pwm1,macspeaks$seqs[[2]],min.score="70%",with.score=T)

## Awk statement to add the binding site position to the last column in the dataframe
#awk -F\\t '{if ($15!~/^NA/){adjstart=$2+$14-1; end=adjstart+length($15);newpos=$1":"adjstart"-"end;print $0"\t"newpos} else {print $0"\t""NA"}}' TargetPeaks_Foxa2PWM1.txt > TargetPeaks_Foxa2PWM_v1.txt


