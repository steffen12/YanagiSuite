import sys
from Bio import SeqIO 
import numpy as np

#python YanagiReadStats.py ~/Penn_Internship/PennSeq/sim_reads/sample_01.fasta

#Output Mean 

if __name__ == "__main__":
	readsFilename = sys.argv[1]
	inputType = sys.argv[2]
	totalReadLengths = []
	for record in SeqIO.parse(readsFilename, inputType):
		totalReadLengths.append(len(record.seq))
	totalReadLengthArray = np.array(totalReadLengths)

	print(str(int(np.mean(totalReadLengthArray))))