#Rscript Simulate_Data.R

library("polyester")
library(Biostrings)

setwd("/home/steffen12/Penn_Internship/YanagiSuite")

args = commandArgs(trailingOnly=TRUE)
print(args)

# FASTA annotation
inputFile = args[1]
outDir = args[2]
fasta = readDNAStringSet(inputFile)
txsNum = length(fasta)
print(width(fasta))
readLen = 101.0

# ~40x coverage ----> reads per transcript = transcriptlength/readlength * 40
# here all transcripts will have ~equal FPKM
readspertx = round(40 * (width(fasta) - readLen + 1) / readLen)
print(readspertx)

simulate_experiment(inputFile, reads_per_transcript=readspertx,
                    num_reps=c(1,1), fold_changes=matrix(c(c(1,1,1),c(1,1,1)), nrow=txsNum), 
                    outdir=outDir,
                    readlen=readLen, paired=FALSE)
