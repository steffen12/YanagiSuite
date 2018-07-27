from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import math, sys

#Run:
#python YanagiQuantTest.py sim_reads/trueTPM.tsv YanagiTPM.tsv KallistoOutput/abundance.tsv

#1) Resolve difference in TPM between Salmon and Kallisto
#2) Resolve Transcripts that YanagiPredicts that aren't in ref transcriptome

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
	SalmonOutputFilename = sys.argv[4]

	trueTPMFile = open(trueTPMFilename, "rU")
	YanagiOutputFile = open(YanagiOutputFilename, "rU")
	KallistoOutputFile = open(KallistoOutputFilename, "rU")
	SalmonOutputFile = open(SalmonOutputFilename, "rU")

	#Read in true TPM
	trueTranscriptTPMsDict = dict()
	trueTranscriptCountsDict = dict()
	line = trueTPMFile.readline() #Get Header
	line = trueTPMFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		count = splitLine[3]
		TPM = splitLine[4]
		trueTranscriptCountsDict[transcriptID] = float(count)
		trueTranscriptTPMsDict[transcriptID] = float(TPM)
		line = trueTPMFile.readline()

	#Read in Yanagi TPM
	yanagiTranscriptTPMsDict = dict()
	yanagiTranscriptCountsDict = dict()
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
		count = splitLine[2]
		TPM = splitLine[3]
		effectiveSegmentLengths = splitLine[5].split(",")
		for i in range(len(effectiveSegmentLengths)):
			effectiveSegmentLengths[i] = round(float(effectiveSegmentLengths[i]))
		yanagiTranscriptSegmentLengthsDict[transcriptID] = effectiveSegmentLengths
		yanagiTranscriptCountsDict[transcriptID] = float(count)
		yanagiTranscriptTPMsDict[transcriptID] = float(TPM)
		line = YanagiOutputFile.readline()

	#Read in Kallisto TPM
	kallistoTranscriptTPMsDict = dict()
	kallistoTranscriptCountsDict = dict()
	line = KallistoOutputFile.readline() #Get Header
	line = KallistoOutputFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		count = splitLine[3]
		TPM = splitLine[4]
		kallistoTranscriptCountsDict[transcriptID] = float(count)
		kallistoTranscriptTPMsDict[transcriptID] = float(TPM)
		line = KallistoOutputFile.readline()

	#Read in Kallisto TPM
	salmonTranscriptTPMsDict = dict()
	salmonTranscriptCountsDict = dict()
	line = SalmonOutputFile.readline() #Get Header
	line = SalmonOutputFile.readline()
	while line:
		splitLine = line.strip().split("\t")
		transcriptID = splitLine[0]
		TPM = splitLine[3]
		count = splitLine[4]
		salmonTranscriptTPMsDict[transcriptID] = float(TPM)
		salmonTranscriptCountsDict[transcriptID] = float(count)
		line = SalmonOutputFile.readline()

	trueTranscriptIDs = set(trueTranscriptTPMsDict.keys())
	yanagiTranscriptIDs = set(yanagiTranscriptTPMsDict.keys())
	kallistoTranscriptIDs = set(kallistoTranscriptTPMsDict.keys())
	salmonTranscriptIDs = set(salmonTranscriptTPMsDict.keys())

	#print(list(yanagiTranscriptIDs.difference(trueTranscriptIDs)))

	sharedYanagiTranscriptIDs = list(trueTranscriptIDs.intersection(yanagiTranscriptIDs))
	sharedKallistoTranscriptIDs = list(trueTranscriptIDs.intersection(kallistoTranscriptIDs))
	sharedSalmonTranscriptIDs = list(trueTranscriptIDs.intersection(salmonTranscriptIDs))

	trueYanagiTPMs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	yanagiTPMs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	trueYanagiCounts = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	yanagiCounts = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	
	trueKallistoTPMs = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))
	kallistoTPMs = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))
	trueKallistoCounts = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))
	kallistoCounts = np.zeros(shape=(len(sharedKallistoTranscriptIDs)))

	trueSalmonTPMs = np.zeros(shape=(len(sharedSalmonTranscriptIDs)))
	salmonTPMs = np.zeros(shape=(len(sharedSalmonTranscriptIDs)))
	trueSalmonCounts = np.zeros(shape=(len(sharedSalmonTranscriptIDs)))
	salmonCounts = np.zeros(shape=(len(sharedSalmonTranscriptIDs)))

	transcriptNum = 0
	for transcriptID in sharedYanagiTranscriptIDs:
		trueYanagiTPMs[transcriptNum] = trueTranscriptTPMsDict[transcriptID]
		yanagiTPMs[transcriptNum] = yanagiTranscriptTPMsDict[transcriptID]
		trueYanagiCounts[transcriptNum] = trueTranscriptCountsDict[transcriptID]
		yanagiCounts[transcriptNum] = yanagiTranscriptCountsDict[transcriptID]
		transcriptNum += 1

	transcriptNum = 0
	for transcriptID in sharedKallistoTranscriptIDs:
		trueKallistoTPMs[transcriptNum] = trueTranscriptTPMsDict[transcriptID]
		kallistoTPMs[transcriptNum] = kallistoTranscriptTPMsDict[transcriptID]
		trueKallistoCounts[transcriptNum] = trueTranscriptCountsDict[transcriptID]
		kallistoCounts[transcriptNum] = kallistoTranscriptCountsDict[transcriptID]
		transcriptNum += 1

	transcriptNum = 0
	for transcriptID in sharedSalmonTranscriptIDs:
		trueSalmonTPMs[transcriptNum] = trueTranscriptTPMsDict[transcriptID]
		salmonTPMs[transcriptNum] = salmonTranscriptTPMsDict[transcriptID]
		trueSalmonCounts[transcriptNum] = trueTranscriptCountsDict[transcriptID]
		salmonCounts[transcriptNum] = salmonTranscriptCountsDict[transcriptID]
		transcriptNum += 1

	#Plot #2
	print("R Yanagi: ", np.corrcoef(math.log(1+trueYanagiCounts), math.log(1+yanagiCounts))[1, 0])
	print("MSE Yanagi: ", computeMSE(trueYanagiCounts, yanagiCounts))
	print("Abs Diff Yanagi: ", computeAbsDiff(trueYanagiCounts, yanagiCounts))
	print()
	plt.scatter(np.log(trueYanagiCounts+1), np.log(yanagiCounts+1))
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1)$")
	plt.title("Yanagi")
	plt.savefig("YanagiQuantFigures/figure_2_1.png")
	plt.clf()

	print("R Kallisto: ", np.corrcoef(math.log(1+trueKallistoCounts), math.log(1+kallistoCounts))[1, 0])
	print("MSE Kallisto: ", computeMSE(trueKallistoCounts, kallistoCounts))
	print("Abs Diff Kallisto: ", computeAbsDiff(trueKallistoCounts, kallistoCounts))
	plt.scatter(np.log(trueKallistoCounts+1), np.log(kallistoCounts+1))
	print()
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1)$")
	plt.title("Kallisto")
	plt.savefig("YanagiQuantFigures/figure_2_2.png")
	plt.clf()

	print("R Salmon: ", np.corrcoef(math.log(1+trueSalmonCounts), math.log(1+salmonCounts))[1, 0])
	print("MSE Salmon: ", computeMSE(trueSalmonCounts, salmonCounts))
	print("Abs Diff Salmon: ", computeAbsDiff(trueSalmonCounts, salmonCounts))
	print()
	plt.scatter(np.log(trueSalmonCounts+1), np.log(salmonCounts+1))
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1)$")
	plt.title("Salmon")
	plt.savefig("YanagiQuantFigures/figure_2_3.png")
	plt.clf()

	yanagiPredictedDifference = np.log(np.subtract(yanagiTPMs, trueYanagiTPMs)+1)

	#Plot #3
	transcriptNumSegs = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptNumSegs[i] = len(yanagiTranscriptSegmentLengthsDict[transcriptID])

	plt.scatter(transcriptNumSegs, yanagiPredictedDifference)
	plt.xlabel(r"Num Segs")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("YanagiQuantFigures/figure_3_1.png")
	plt.clf()

	transcriptMinSegLength = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptMinSegLength[i] = min(yanagiTranscriptSegmentLengthsDict[transcriptID])

	plt.scatter(np.log(transcriptMinSegLength), yanagiPredictedDifference)
	plt.xlabel(r"log(Min Length of Segment)")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("YanagiQuantFigures/figure_3_2.png")
	plt.clf()

	transcriptMedianSegLength = np.zeros(shape=(len(sharedYanagiTranscriptIDs)))
	for i in range(len(sharedYanagiTranscriptIDs)):
		transcriptID = sharedYanagiTranscriptIDs[i]
		transcriptMedianSegLength[i] = np.median(np.array(yanagiTranscriptSegmentLengthsDict[transcriptID]))

	plt.scatter(np.log(transcriptMedianSegLength), yanagiPredictedDifference)
	plt.xlabel(r"log(Median Length of Segment)")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.savefig("YanagiQuantFigures/figure_3_3.png")
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
	plt.savefig("YanagiQuantFigures/figure_3_4.png")
	plt.clf()

	#PLot #4
	plt.scatter(np.log(trueYanagiTPMs+1), yanagiPredictedDifference)
	plt.xlabel(r"$log(\Theta+1)$")
	plt.ylabel(r"$log({\Theta}_{pred}+1) - log(\Theta+1)$")
	plt.title("Yanagi")
	plt.savefig("YanagiQuantFigures/figure_4_1.png")
	plt.clf()
