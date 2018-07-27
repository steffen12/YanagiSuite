#python YanagiInputSlicer.py YanagiSegmenter/homo_sapiens/homo_sapiens_segs_151.fa homo_sapiens_txome.fa 1 

#Watch out for transcriptname mismatching! In txome have to remove . 

import sys
from Bio import SeqIO

def getGeneId(record):
	ID, name, txs, segtype = record.description.split(" ")
	tokens = name.split(":")
	geneID = tokens[1]
	return geneID

def getTranscriptId(recordId):
	return recordId.split(".")[0]

def getTranscripts(record):
	ID, name, txs, segtype = record.description.split(" ")
	return set(txs.split(":")[1].split(","))

if __name__ == "__main__":
	segmentsFilename = sys.argv[1]
	transcriptomeFilename = sys.argv[2]
	numSegments = int(sys.argv[3])

	transcriptsSet = set()

	segmentsKeepList = []
	segmentIndex = 0

	currentGene = ""

	for record in SeqIO.parse(segmentsFilename, "fasta"):
		nextGene = getGeneId(record)
		if segmentIndex < numSegments or nextGene == currentGene:
			currentGene = nextGene
			segmentsKeepList.append(record)
			segmentIndex += 1
		else:
			break

	for segment in segmentsKeepList:
		transcriptsSet = transcriptsSet.union(getTranscripts(segment))

	transcriptomeKeepList = []

	for record in SeqIO.parse(transcriptomeFilename, "fasta"):
		record.id = getTranscriptId(record.id)
		record.description = ""
		if record.id in transcriptsSet:
			transcriptomeKeepList.append(record)
			transcriptsSet = transcriptsSet - set(record.id)

		if len(transcriptsSet) == 0:
			break

	newSegmentsFilename = segmentsFilename + "_" + str(segmentIndex) + ".fa"
	newTranscriptomeFilename = transcriptomeFilename + "_" + str(segmentIndex) + ".fa"

	SeqIO.write(segmentsKeepList, newSegmentsFilename, "fasta")
	SeqIO.write(transcriptomeKeepList, newTranscriptomeFilename, "fasta")