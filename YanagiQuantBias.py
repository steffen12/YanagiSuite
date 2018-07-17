#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
from Bio.SeqUtils import GC
from Bio import SeqIO 
import math, sys, os, time
import numpy as np
import matplotlib.pyplot as plt
#from multiprocessing import Process, Queue
from Queue import Queue

###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
def getArgs():
    validArgList = ["-counts", "-segments", "-readlen", "-precision", "-out"]
    for argIndex in range(1,len(sys.argv)):
        if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
            print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
            sys.exit()
            
    countFileExists = False
    segmentFileExists = False
    readLengthExists = False
    precisionExists = False
    outFileExists = False
    argIndex = 1

    while argIndex < len(sys.argv):
        if sys.argv[argIndex] == "-counts":  ## load in counts file
            argIndex += 1
            countFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
            countTmp = sys.argv[argIndex].split("/")
            countFile = countFileAbsPath + "/" + countTmp[len(countTmp)-1]
            countFileExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-segments":  ## load in counts file
            argIndex += 1
            segmentFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
            segmentTmp = sys.argv[argIndex].split("/")
            segmentFile = segmentFileAbsPath + "/" + segmentTmp[len(segmentTmp)-1]
            segmentFileExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-readlen":  ## load in read len
            argIndex += 1
            readLength = int(sys.argv[argIndex])
            readLengthExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-precision":  ## load in precision
            argIndex += 1
            diffMax = float(sys.argv[argIndex])
            precisionExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-out":  ## load in output file
            argIndex += 1
            outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
            outTmp = sys.argv[argIndex].split("/")
            outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
            outFileExists = True
            argIndex += 1                  

    if (not countFileExists) or (not segmentFileExists) or (not readLengthExists) or (not precisionExists) or (not outFileExists): ## lack enough arguments
        print("Please provide arguments:")
        print("-counts\tSegment count file")
        print("-segments\tSegment .fa file")
        print("-readlen\tRead Length")
        print("-precision\tPrecision")
        print("-out\tOutput file")
        sys.exit()

    return countFile, segmentFile, readLength, diffMax, outFile

def findBiasTrend(pairedEndMode, countFileLines, segmentFile):

    segmentGCDict = dict()
    for record in SeqIO.parse(segmentFile, "fasta"):
        segmentGCDict[record.id] = GC(record.seq)

    segmentIDs = []
    segmentCountsDict = dict()
    segmentCounts = []
    segmentLengths = []
    segmentIsoforms = dict()
    geneIsoforms = []
    segmentPositionsDict = dict()

    countFileLineIndex = 0
    countFileLine = countFileLines[countFileLineIndex]

    splitLine = countFileLine.strip().split("\t")

    if pairedEndMode:
        currGene = splitLine[4] #Get Current Gene
        while countFileLineIndex < countFileLineEnd and splitLine[4] == currGene:
            segmentID1 = splitLine[0]
            segmentID2 = splitLine[1]
            if segmentID1 == segmentID2:
                segmentIDs.append(segmentID1)
                segmentLengths.append(int(splitLine[5]))
                segmentIsoforms[segmentID1] = splitLine[9].split(",")
                for isoform in segmentIsoforms[segmentID1]:
                    if isoform not in geneIsoforms:
                        geneIsoforms.append(isoform)
            if segmentID1 not in segmentCountsDict:
                segmentCountsDict[segmentID1] = int(splitLine[2])
            else:
                segmentCountsDict[segmentID1] += int(splitLine[2])

            if segmentID2 not in segmentCountsDict:
                segmentCountsDict[segmentID2] = int(splitLine[2])
            else:
                segmentCountsDict[segmentID2] += int(splitLine[2])

            countFileLineIndex += 1
            if countFileLineIndex < countFileLineEnd:
                countFileLine = countFileLines[countFileLineIndex]
                splitLine = countFileLine.strip().split("\t")
    else:
        currGene = splitLine[3] #Get Current Gene
        while countFileLineIndex < countFileLineEnd and splitLine[3] == currGene:
            segmentID = splitLine[0]
            segmentIDs.append(segmentID)
            segmentCountsDict[segmentID] = int(splitLine[1])
            segmentLengths.append(int(splitLine[4]))

            segmentIsoforms[segmentID] = splitLine[6].split(",")
            for isoform in segmentIsoforms[segmentID]:
                if isoform not in geneIsoforms:
                    geneIsoforms.append(isoform)

            countFileLineIndex += 1
            if countFileLineIndex < countFileLineEnd:
                countFileLine = countFileLines[countFileLineIndex]
                splitLine = countFileLine.strip().split("\t")

    for segmentID in segmentIDs: #Only get segment counts of certain segments
        segmentCounts.append(segmentCountsDict[segmentID])

    segmentGC = [0.0 for x in range(len(segmentIDs))]
    segmentCount = [0.0 for x in range(len(segmentIDs))]

    for x in range(len(segmentIDs)):
        segmentID = segmentIDs[x]
        segmentGC[x] = segmentGCDict[segmentID]
        segmentCount[x] = math.log(segmentCountsDict[segmentID])

    plt.scatter(segmentGC, segmentCount)
    plt.show()

###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

def computeIsoformAbundances(countFileLines, readLength, diffMax):
    geneCount = 0

    countFileLineIndexStart = 0
    countFileLineEnd = len(countFileLines)
    print(countFileLineEnd)

    countFileLineIndex = countFileLineIndexStart
    countFileLine = countFileLines[countFileLineIndex]

    while countFileLineIndex < countFileLineEnd:
        segmentIDs = []
        segmentCountsDict = dict()
        segmentCounts = []
        segmentLengths = []
        segmentIsoforms = dict()
        geneIsoforms = []

        splitLine = countFileLine.strip().split("\t")

        if pairedEndMode:
            currGene = splitLine[4] #Get Current Gene
            while countFileLineIndex < countFileLineEnd and splitLine[4] == currGene:
                segmentID1 = splitLine[0]
                segmentID2 = splitLine[1]
                if segmentID1 == segmentID2:
                    segmentIDs.append(segmentID1)
                    segmentLengths.append(int(splitLine[5]))
                    segmentIsoforms[segmentID1] = splitLine[9].split(",")
                    for isoform in segmentIsoforms[segmentID1]:
                        if isoform not in geneIsoforms:
                            geneIsoforms.append(isoform)
                if segmentID1 not in segmentCountsDict:
                    segmentCountsDict[segmentID1] = int(splitLine[2])
                else:
                    segmentCountsDict[segmentID1] += int(splitLine[2])

                if segmentID2 not in segmentCountsDict:
                    segmentCountsDict[segmentID2] = int(splitLine[2])
                else:
                    segmentCountsDict[segmentID2] += int(splitLine[2])

                countFileLineIndex += 1
                if countFileLineIndex < countFileLineEnd:
                    countFileLine = countFileLines[countFileLineIndex]
                    splitLine = countFileLine.strip().split("\t")
        else:
            currGene = splitLine[3] #Get Current Gene
            while countFileLineIndex < countFileLineEnd and splitLine[3] == currGene:
                segmentID = splitLine[0]
                segmentIDs.append(segmentID)
                segmentCountsDict[segmentID] = int(splitLine[1])
                segmentLengths.append(int(splitLine[4]))

                segmentIsoforms[segmentID] = splitLine[6].split(",")
                for isoform in segmentIsoforms[segmentID]:
                    if isoform not in geneIsoforms:
                        geneIsoforms.append(isoform)

                countFileLineIndex += 1
                if countFileLineIndex < countFileLineEnd:
                    countFileLine = countFileLines[countFileLineIndex]
                    splitLine = countFileLine.strip().split("\t")

        for segmentID in segmentIDs: #Only get segment counts of certain segments
            segmentCounts.append(segmentCountsDict[segmentID])
        readCount = sum(segmentCounts)

        if readCount == 0:
            continue

        segmentIsoformIndicatorMatrix = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)), dtype="bool_")
        for i in range(len(segmentIDs)):
            segmentID = segmentIDs[i]
            segmentIsoformList = segmentIsoforms[segmentID]
            for isoform in segmentIsoformList:
                j = geneIsoforms.index(isoform)
                segmentIsoformIndicatorMatrix[i, j] = True

        effectiveSegmentLengths = np.zeros(shape=(len(segmentIDs)), dtype="int_")
        for i in range(len(segmentIDs)):
            effectiveSegmentLengths[i] = segmentLengths[i] - readLength + 1

        isoformCounts = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
        isoformLengths = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
        for j in range(len(geneIsoforms)):
            for i in range(len(segmentIDs)):
                if(segmentIsoformIndicatorMatrix[i,j]):
                    isoformCounts[j] += segmentCounts[i]
                    isoformLengths[j] += effectiveSegmentLengths[i]
        
        ############################################################################################################################################
        ## Find H for each segment
        ############################################################################################################################################
        
        segmentH = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
        for i in range(len(segmentIDs)):
            for j in range(len(geneIsoforms)):
                if(segmentIsoformIndicatorMatrix[i,j]):
                    segmentH[i, j] = float(segmentCounts[i])/(isoformCounts[j] * effectiveSegmentLengths[i])

        print("Time to Analyze Distribution", time.time() - startTime)

    #     #####################################################################################################
    #     ## EM algorithm
    #     #####################################################################################################

        Thetas = np.full(shape=len(geneIsoforms), fill_value=1.0/len(geneIsoforms))
        Alpha = np.zeros(shape=len(geneIsoforms)) #This is Theta with ~ over it
        for j in range(len(geneIsoforms)):
            Alpha[j] = Thetas[j] * isoformLengths[j]
        sumAlpha = np.sum(Alpha)
        for j in range(len(geneIsoforms)):
            Alpha[j] /= sumAlpha
        oldAlpha = np.zeros(shape=len(geneIsoforms))

         #########################################################################################################
         ## iteration begins        
                
        diff = 1.0
        iterCount = 0
        Z = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
        while diff > diffMax:
            #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+str(diff)+"\t"+str(tmpTime))

            #Expectation Step
            for i in range(len(segmentIDs)):
                probNumerators = np.multiply(Alpha, segmentH[i])
                probDenominator = np.sum(probNumerators)
                Z[i] = np.multiply(probNumerators, segmentCounts[i] / probDenominator)   

            #Maximization Step
            oldAlpha = Alpha
            Alpha = np.multiply(np.sum(Z, axis=0), 1/float(readCount))
            
            diffArray = np.absolute(np.subtract(Alpha, oldAlpha))
            diff = diffArray.max()

            iterCount += 1

                    
        sumAlpha = np.sum(Alpha)
        if sumAlpha == 0: 
            continue

        Alpha = np.multiply(Alpha, 1.0 / sumAlpha)

        Thetas = np.zeros(shape=(len(geneIsoforms),))
        for j in range(len(Alpha)):
            Thetas[j] = Alpha[j] / isoformLengths[j]
        sumTheta = np.sum(Thetas)

        print(str(geneCount) + ": " + str(currGene)+"\t"+str(iterCount)+" iterations\tDone!")

        outputStrings = ["" for i in range(len(Alpha))]
        for i in range(len(Alpha)):
            Thetas[i] /= sumTheta
            #print(currGene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+isoformNames[i]+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(tmpTime))
            outputString = currGene+"\t"+geneIsoforms[i]+"\t"+str(isoformCounts[i])+"\t"+str(Thetas[i])+"\n"
            outputStrings[i] = outputString

        outputFileLines.put(outputStrings)

        print("Time for EM", time.time() - startTime)
        geneCount += 1

if __name__ == "__main__":
    countFile, segmentFile, readLength, diffMax, outFile = getArgs()

    countFileHandle = open(countFile, 'r')
    OUT = open(outFile, 'w')
    outputFileLines = Queue()

    startTime = time.time()

    countFileLine = countFileHandle.readline() #Analyze header
    splitLine = countFileLine.strip().split("\t")

    if len(splitLine) == 7:
        pairedEndMode = False
    elif len(splitLine) == 10:
        pairedEndMode = True
    else:
        print("Invalid Counts Input File")
        sys.exit()

    OUT.write("GeneName\tIsoformName\tNumberOfReads\tRelativeAbundance\n") ## Header of Results

    countFileLines = countFileHandle.readlines()

    findGCBiasTrend(countFileLines, segmentFile)

    countFileLineIndexStart = 0
    countFileLineIndexEnd = len(countFileLines)
    computeIsoformAbundances(countFileLines[countFileLineIndexStart: countFileLineIndexEnd], readLength, diffMax)
    #p = Process(target=computeIsoformAbundances, args=(countFileLines[countFileLineIndexStart: countFileLineIndexEnd], readLength, diffMax))
    #p.start()
    #p.join()

    while not outputFileLines.empty():
        outputStrings = outputFileLines.get()
        for outputFileLine in outputStrings:
            #print(outputFileLine)
            OUT.write(outputFileLine)

    countFileHandle.close()
    OUT.close()
    print(time.time() - startTime)