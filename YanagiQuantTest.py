import matplotlib.pyplot as plt
import numpy as np
import math, sys

#Run:
#python YanagiQuantTest.py sim_reads/trueTPM.tsv YanagiTPM.tsv KallistoOutput/abundance.tsv

def computeMSE(trueArray, predictedArray):
	MSE = 0.0
	for i in range(len(trueArray)):
		MSE += (trueArray[i] - predictedArray[i])**2
	MSE /= len(trueArray)
	return MSE

def computeAbsDiff(trueArray, predictedArray):
	absDiffArray = np.abs(np.subtract(trueArray, predictedArray))
	absDiff = np.mean(absDiffArray)
	return absDiff

if __name__ == "__main__":
	trueTPMFilename = sys.argv[1]
	YanagiOutputFilename = sys.argv[2]
	KallistoOutputFilename = sys.argv[3]

	trueTPMFile = open(trueTPMFilename, "rU")
	YanagiOutputFile = open(YanagiOutputFilename, "rU")
	KallistoOutputFile = open(KallistoOutputFilename, "rU")

	#Read in true TPM
	trueTranscriptTPMsDict = dict()
	line = trueTPMFile.readline() #Get Header
	line = trueTPMFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		TPM = splitLine[4]
		trueTranscriptTPMsDict[transcriptID] = float(TPM)
		line = trueTPMFile.readline()

	#Read in Yanagi TPM
	yanagiTranscriptTPMsDict = dict()
	yanagiTranscriptSegmentLengthsDict = dict()
	yanagiGeneTranscriptsDict = dict()
	line = YanagiOutputFile.readline() #Get Header
	line = YanagiOutputFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		geneID = splitLine[1]
		if (geneID not in yanagiGeneTranscriptsDict):
			yanagiGeneTranscriptsDict[geneID] = [transcriptID]
		else:
			yanagiGeneTranscriptsDict[geneID].append(transcriptID)
		TPM = splitLine[2]
		effectiveSegmentLengths = splitLine[4].split(",")
		for i in range(len(effectiveSegmentLengths)):
			effectiveSegmentLengths[i] = int(effectiveSegmentLengths[i])
		yanagiTranscriptSegmentLengthsDict[transcriptID] = effectiveSegmentLengths
		yanagiTranscriptTPMsDict[transcriptID] = float(TPM)
		line = YanagiOutputFile.readline()

	#Read in Kallisto TPM
	kallistoTranscriptTPMsDict = dict()
	line = KallistoOutputFile.readline() #Get Header
	line = KallistoOutputFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		TPM = splitLine[4]
		kallistoTranscriptTPMsDict[transcriptID] = float(TPM)
		line = KallistoOutputFile.readline()

	trueTranscriptIDs = set(trueTranscriptTPMsDict.keys())
	yanagiTranscriptIDs = set(yanagiTranscriptTPMsDict.keys())
	sharedYanagiTranscriptIDs = list(trueTranscriptIDs.intersection(yanagiTranscriptIDs))
	kallistoTranscriptIDs = set(kallistoTranscriptTPMsDict.keys())
	sharedKallistoTranscriptIDs = list(set(sharedYanagiTranscriptIDs).intersection(kallistoTranscriptIDs))

	trueYanagiTPMs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	yanagiTPMs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	trueKallistoTPMs = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))
	kallistoTPMs = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))

	numTPM = 0
	for transcriptID in sharedYanagiTranscriptIDs:
		trueYanagiTPMs[numTPM] = trueTranscriptTPMsDict[transcriptID]
		yanagiTPMs[numTPM] = yanagiTranscriptTPMsDict[transcriptID]
		numTPM += 1

	numTPM = 0
	for transcriptID in sharedKallistoTranscriptIDs:
		trueKallistoTPMs[numTPM] = trueTranscriptTPMsDict[transcriptID]
		kallistoTPMs[numTPM] = kallistoTranscriptTPMsDict[transcriptID]
		numTPM += 1

	#Plot #2
	print("R Yanagi: ", np.corrcoef(trueYanagiTPMs, yanagiTPMs))
	print("MSE Yanagi: ", computeMSE(trueYanagiTPMs, yanagiTPMs))
	print("Abs Diff Yanagi: ", computeAbsDiff(trueYanagiTPMs, yanagiTPMs))
	plt.scatter(np.log(trueYanagiTPMs+1), np.log(yanagiTPMs+1))
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1)$")
	plt.title("Yanagi")
	plt.savefig("figure_2_1.png")
	plt.clf()

	print("R Kallisto: ", np.corrcoef(trueKallistoTPMs, kallistoTPMs))
	print("MSE Kallisto: ", computeMSE(trueKallistoTPMs, kallistoTPMs))
	print("Abs Diff Kallisto: ", computeAbsDiff(trueKallistoTPMs, kallistoTPMs))
	plt.scatter(np.log(trueKallistoTPMs+1), np.log(kallistoTPMs+1))
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1)$")
	plt.title("Kallisto")
	plt.savefig("figure_2_2.png")
	plt.clf()

	yanagiPredictedDifference = np.subtract(np.log(yanagiTPMs+1), np.log(trueYanagiTPMs+1))

	#Plot #3
	transcriptNumSegs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptNumSegs[i] = len(yanagiTranscriptSegmentLengthsDict[transcriptID])

	plt.scatter(transcriptNumSegs, yanagiPredictedDifference)
	plt.xlabel(r"Num Segs")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("figure_3_1.png")
	plt.clf()

	transcriptMinSegLength = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptMinSegLength[i] = min(yanagiTranscriptSegmentLengthsDict[transcriptID])

	plt.scatter(np.log(transcriptMinSegLength), yanagiPredictedDifference)
	plt.xlabel(r"log(Min Length of Segment)")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("figure_3_2.png")
	plt.clf()

	transcriptMedianSegLength = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptMedianSegLength[i] = np.median(np.array(yanagiTranscriptSegmentLengthsDict[transcriptID]))

	plt.scatter(np.log(transcriptMedianSegLength), yanagiPredictedDifference)
	plt.xlabel(r"log(Median Length of Segment)")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("figure_3_3.png")
	plt.clf()

	geneNumTranscripts = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	yanagiPredictedDifference2 = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	numTranscript = 0
	for geneID in yanagiGeneTranscriptsDict:
		transcriptIDs = yanagiGeneTranscriptsDict[geneID]
		for transcriptID in transcriptIDs:
			if transcriptID in trueTranscriptTPMsDict:
				geneNumTranscripts[numTranscript] = len(transcriptIDs)
				yanagiPredictedDifference2[numTranscript] = math.log(yanagiTranscriptTPMsDict[transcriptID]+1) - math.log(trueTranscriptTPMsDict[transcriptID]+1)
				numTranscript += 1
	
	plt.scatter(geneNumTranscripts, yanagiPredictedDifference2)
	plt.xlabel(r"Num Transcripts")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("figure_3_4.png")
	plt.clf()

	#PLot #4
	plt.scatter(np.log(trueYanagiTPMs+1), yanagiPredictedDifference)
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.title("Yanagi")
	plt.savefig("figure_4_1.png")
	plt.clf()
