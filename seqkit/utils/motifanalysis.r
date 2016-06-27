## Program to identify the Motifs in the whole set of peaks; checkout the most common motifs in top,bottom and random peaks
## Once PWM is defined, checkout the locations of Binding sites in the given (Target) peaks

## Run as Rscript motifanalysis.r ip_file op_dir op_nm

.libPaths("/pica/h1/ashwini/anaconda2/lib/R/library")
library("BSgenome")
library("BSgenome.Mmusculus.UCSC.mm10.masked")
library("lattice")
library("seqLogo")

# function to convert MACS output table into a GRanges object
macs2GRanges <-function(peaks) {
  # generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
    length=peaks$length,
    pileup_summit=peaks$pileup,
    score=peaks$X.log10.pvalue,
    FE=peaks$fold_enrichment,
    fdr=peaks$X.log10.qvalue,
    maxpos=peaks$abs_summit,
    name=peaks$name  
  )
  return(myrange)
}
args=(commandArgs(TRUE))
macspeaks <-read.table(paste0(args[1]), header=TRUE, sep="\t", stringsAsFactors=F)#[,-c(11,12)]
macspeaks <- macs2GRanges(macspeaks)



###################### De-novo motif analysis ##############################


## De-novo motif discorvery -MEME analysis
## use a distance of 100 bp on either side of the summits
mydist <- 100

# define summit regions
macs.summits <- GRanges(
  seqnames=seqnames(macspeaks),
  range=IRanges( start=elementMetadata(macspeaks)$maxpos - mydist,
                 end=elementMetadata(macspeaks)$maxpos + mydist),
  strand="+"
)

# order the summit regions based on their score
macs.summits <- macs.summits[order(-elementMetadata(macspeaks)$score)]
# create unique names for all peaks
names(macs.summits) <- paste( seqnames(macs.summits), start(macs.summits), end(macs.summits), sep=":")
# take a look at the resulting GRanges centered around the MACS peak summit positions
macs.summits
out.folder <- args[2]
fl_nm <- args[3]

# ## Collect all peak sequences
seqs <- getSeq(Mmusculus, macs.summits, as.character=FALSE)
names(seqs) <- names(macs.summits)
writeXStringSet(seqs,file=paste(out.folder,fl_nm,sep="/"), format="fasta", append=FALSE)


## MOtif Analysis 
## generate the meme command and run this command on uppmax command line
# get list of fasta files
fa.files <- list.files(paste(out.folder,sep=""), pattern="*fa$", full.names=TRUE)
# define MEME parameters
#parameters <- " -nmotifs 3 -minsites 100 -minw 8 -maxw 35 -revcomp -maxsize 5000000 -dna -oc "
params <- " -db /sw/apps/bioinfo/MEMEsuite/4.11.1/milou/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme -meme-maxsize 500000 -meme-minsites 100 -oc "
for (fa.file in fa.files) {
meme.out <- substr(fa.file,1,nchar(fa.file)-3)
#construct full MEME command
mycommand <- paste("meme-chip",params,meme.out,fa.file,sep=" ")
system(mycommand)
print(mycommand)
}

