import sys
from Bio import SeqIO 
import numpy as np

#python YanagiReadStats.py ~/Penn_Internship/PennSeq/sim_reads/sample_01.fasta

if __name__ == "__main__":
	readsFilename = sys.argv[1]
	totalReadLengths = []
	for record in SeqIO.parse(readsFilename, "fasta"):
		totalReadLengths.append(len(record.seq))
	totalReadLengthArray = np.array(totalReadLengths)

	print("Mean: " + str(np.mean(totalReadLengthArray)))
	print("Std: " + str(np.std(totalReadLengthArray)))