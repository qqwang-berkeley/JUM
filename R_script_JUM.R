args = commandArgs(trailingOnly=TRUE)

# test if experiment design input file is provided; if not, return an error
if (length(args)==0) {
  stop("An experiment design file must be supplied (input file).n", call.=FALSE)
}

source("http://bioconductor.org/biocLite.R")
biocLite("DEXSeq", dependencies=TRUE)
library(DEXSeq)

annotationfile = file.path("combined_AS_JUM.gff")

sampledesign = read.table(args[1], header=TRUE)
sampledesign
fullFilenames <- dir(pattern ="combined_count.txt")
fullFilenames

JUM = DEXSeqDataSetFromHTSeq(
fullFilenames,
sampleData=sampledesign,
design= ~ sample + exon + condition:exon,
flattenedfile=annotationfile)

JUM = estimateSizeFactors(JUM)
JUM = estimateDispersions(JUM)
JUM = testForDEU(JUM)
JUM = estimateExonFoldChanges(JUM, fitExpToVar="condition")
dxr1 = DEXSeqResults(JUM)
write.table(dxr1, "AS_differential.txt", sep="\t", quote=F)
q()

