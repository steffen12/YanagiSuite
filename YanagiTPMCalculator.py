import sys
from Bio import SeqIO 

#python YanagiTPMCalculator.py sim_reads/sample_01.fasta txome_ENSG00000100842.fa 101 sim_reads/trueTPM.tsv

def calculateTPMs(readsType, readsFilename, transcriptomeFilename, readLength, outputFilename):
	outputFile = open(outputFilename, 'w')
	outputFile.write("transcriptId" + "\t" + "length" + "\t" + "effLength" + "\t" + 
			"count" + "\t" + "TPM" + "\n")
	transcriptCountDict = {}

	for record in SeqIO.parse(readsFilename, "fastq"):
		if readsType == "polyester":
			readId = record.id.split(";")[0]
			transcriptId = readId.split("/")[1]
		elif readsType == "rlsim":
			readId = record.id.split(" ")[0]
			transcriptId = readId.split("_")[2]

		if transcriptId in transcriptCountDict:
			transcriptCountDict[transcriptId] += 1
		else:
			transcriptCountDict[transcriptId] = 1

	transcriptLengthDict = {}

	for record in SeqIO.parse(transcriptomeFilename, "fasta"):
		transcriptId = record.id
		transcriptLength = len(record.seq)
		transcriptEffectiveLength = transcriptLength - readLength + 1
		transcriptLengthDict[transcriptId] = (transcriptLength, transcriptEffectiveLength)

	transcriptIds = sorted(transcriptLengthDict.keys())
	TPMTotal = 0
	for transcriptId in transcriptIds:
		if transcriptId in transcriptCountDict:
			transcriptCount = transcriptCountDict[transcriptId]
		else:
			transcriptCount = 0
		transcriptLength, transcriptEffectiveLength = transcriptLengthDict[transcriptId]
		TPMTotal += float(transcriptCount)/transcriptEffectiveLength

	for transcriptId in transcriptIds:
		if transcriptId in transcriptCountDict:
			transcriptCount = transcriptCountDict[transcriptId]
		else:
			transcriptCount = 0
		transcriptLength, transcriptEffectiveLength = transcriptLengthDict[transcriptId]
		transcriptTPM = (float(transcriptCount)/transcriptEffectiveLength)/TPMTotal*1e6
		outputFile.write(transcriptId + "\t" + str(transcriptLength) + "\t" + str(transcriptEffectiveLength) + "\t" + 
			str(transcriptCount) + "\t" + str(transcriptTPM) + "\n")

if __name__ == "__main__":
	readsType = sys.argv[1]
	readsFilename = sys.argv[2]
	transcriptomeFilename = sys.argv[3]
	readLength = float(sys.argv[4])
	outputFilename = sys.argv[5]
	calculateTPMs(readsType, readsFilename, transcriptomeFilename, readLength, outputFilename)