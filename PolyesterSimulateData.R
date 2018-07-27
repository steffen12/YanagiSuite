#Rscript Simulate_Data.R

library("polyester")
library(Biostrings)

setwd("/home/steffen12/Penn_Internship/YanagiSuite")

args = commandArgs(trailingOnly=TRUE)

# FASTA annotation
inputFile = args[1]
outDir = args[2]
paired = (args[3] == "paired")
print(paste("Paired:", paired))
fasta = readDNAStringSet(inputFile)
txsNum = length(fasta)
readLen = 151.0

# ~40x coverage ----> reads per transcript = transcriptlength/readlength * 40
# here all transcripts will have ~equal FPKM
readspertx = round(40 * (width(fasta)) / readLen)
fold_changes_rep1 = runif(txsNum, min=0, max=2)
fold_changes = matrix(c(fold_changes_rep1, fold_changes_rep1), nrow=txsNum)
print(fold_changes)

simulate_experiment(inputFile, reads_per_transcript=readspertx,
                    num_reps=c(1, 0), fold_changes=fold_changes, 
                    outdir=outDir,
                    readlen=readLen, paired=paired, fraglen=250, fragsd=25)
