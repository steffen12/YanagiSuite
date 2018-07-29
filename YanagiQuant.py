#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
from Bio.SeqUtils import GC
from Bio import SeqIO 
import math, sys, os, time
from sklearn import linear_model
from math import log
import numpy as np
from multiprocessing import Pool
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
#from Queue import Queue

#Todo:
#1) Small Transcripts, find length for TPM - ask about
#2) Paired End Reads
#3) YanagiTPM.tsv format - length and effective length
#4) Base on variance negative bionomial
#5) Test to see if multiple starting Thetas helps by doing min theta

#Run:
#python YanagiQuant.py -counts Hs_Counts.tsv -readlen 101 -fraglen 100 -out Output_Result_FILE2
#python YanagiQuant.py -counts unit_tests/Test_Input.txt -readlen 101 -out unit_tests/Test_Output.tsv
#python YanagiQuant.py -counts YanagiCountOutput/seg_counts.tsv -readlen 101 -fraglen 101 -out YanagiTPM.tsv

###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
def getArgs():
    validArgList = ["-counts", "-readlen", "-fraglen", "-precision", "-numprocs", "-out"]
    for argIndex in range(1,len(sys.argv)):
        if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
            print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
            sys.exit()
            
    countFileExists = False
    readLengthExists = False
    outFileExists = False
    argIndex = 1

    diffMax = 0.001 #Default
    numProcesses = 1 #Default
    fragmentLength = -1 #Default

    while argIndex < len(sys.argv):
        if sys.argv[argIndex] == "-counts":  ## load in counts file
            argIndex += 1
            countFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
            countTmp = sys.argv[argIndex].split("/")
            countFile = countFileAbsPath + "/" + countTmp[len(countTmp)-1]
            countFileExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-readlen":  ## load in read len
            argIndex += 1
            readLength = float(sys.argv[argIndex])
            readLengthExists = True
            argIndex += 1
        elif sys.argv[argIndex] == "-fraglen":  ## load in fragment len
            argIndex += 1
            fragmentLength = float(sys.argv[argIndex])
            argIndex += 1
        elif sys.argv[argIndex] == "-precision":  ## load in precision
            argIndex += 1
            diffMax = float(sys.argv[argIndex])
            argIndex += 1
        elif sys.argv[argIndex] == "-numprocs":  ## load in num processes
            argIndex += 1
            numProcesses = int(sys.argv[argIndex])
            argIndex += 1
        elif sys.argv[argIndex] == "-out":  ## load in output file
            argIndex += 1
            outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
            outTmp = sys.argv[argIndex].split("/")
            outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
            outFileExists = True
            argIndex += 1                  

    if (not countFileExists) or (not readLengthExists) or (not outFileExists): ## lack enough arguments
        print("Please provide arguments:")
        print("-counts\tSegment count file")
        print("-readlen\tRead Length")
        print("-out\tOutput file")
        sys.exit()

    return countFile, readLength, fragmentLength, diffMax, numProcesses, outFile

def detectPairedEnd(countFileHandle, fragmentLength):
    countFileLine = countFileHandle.readline() #Analyze header
    splitLine = countFileLine.strip().split("\t")

    pairedEndMode = False

    if len(splitLine) == 7:
        pairedEndMode = False
    elif len(splitLine) == 10:
        pairedEndMode = True
    else:
        print("Invalid Counts Input File")
        sys.exit()

    if not pairedEndMode and fragmentLength == -1:
        print("Please provide fragment length for single end mode:")
        print("-fraglen\tfragment length")
        sys.exit()

    return pairedEndMode

def calculatePairedFragmentLength(countFilename, readLength):
    countFile = open(countFilename, 'rU')
    countFile.readline()

    segmentCounts = dict()
    segmentLengths = dict()
    segmentJunctionCounts = dict()

    fragmentLengthArray = []

    for line in countFile.readlines():
        splitLine = line.strip().split("\t")
        segmentID1 = splitLine[0]
        segmentID2 = splitLine[1]
        segmentPairCount = int(splitLine[2])
        segment1Length = int(splitLine[5])
        segmentLengths[segmentID1] = segment1Length
        if segmentID1 == segmentID2:
            segmentCounts[segmentID1] = segmentPairCount
        else:
            if segmentID1 not in segmentJunctionCounts:
                segmentJunctionCounts[segmentID1] = segmentPairCount
            else:
                segmentJunctionCounts[segmentID1] += segmentPairCount

    for segmentID in segmentCounts:
        if segmentID in segmentJunctionCounts:
            fragmentLengthNumerator = segmentJunctionCounts[segmentID]*(segmentLengths[segmentID] + 1) + segmentCounts[segmentID]*(readLength - 1)
            fragmentLengthDenominator = float(segmentCounts[segmentID] + segmentJunctionCounts[segmentID])
            fragmentLength = fragmentLengthNumerator/fragmentLengthDenominator
            fragmentLengthArray.append(fragmentLength)

    #plt.hist(fragmentLengthArray)
    #plt.show()
    if len(fragmentLengthArray) != 0:
        fragmentDensity = gaussian_kde(fragmentLengthArray)
        x = np.linspace(0, max(fragmentLengthArray), num=5*max(fragmentLengthArray))
        #plt.plot(x, fragmentDensity(x))
        #plt.show()
        maxFragmentIndex = np.argmax(fragmentDensity(x))
        fragmentLength = x[maxFragmentIndex]
    else:
        fragmentLength = readLength
    print("Mean Fragment Length: ", fragmentLength)
    return fragmentLength

def partitionCountFile(countFileLines, countFileLength, pairedEndMode, numProcesses):
    processFileLength = int(countFileLength / numProcesses)
    countFileLineStarts = [0]
    countFileLineEnds = []
    for i in range(1, numProcesses):
        index = i*processFileLength
        line = countFileLines[index]
        splitLine = line.strip().split("\t")
        if pairedEndMode:
            currGene = splitLine[4] #Get Current Gene
            while index < countFileLength and splitLine[4] == currGene:
                index += 1
                if index < countFileLength:
                    line = countFileLines[index]
                    splitLine = line.strip().split("\t")
        else:
            currGene = splitLine[3] #Get Current Gene
            while index < countFileLength and splitLine[4] == currGene:
                index += 1
                if index < countFileLength:
                    line = countFileLines[index]
                    splitLine = line.strip().split("\t")
        countFileLineStarts.append(index)
        countFileLineEnds.append(index)
    countFileLineEnds.append(countFileLength)
    return countFileLineStarts, countFileLineEnds

def calculateIsoformTPMs(segmentIDs, geneIsoforms, segmentIsoformIndicatorMatrix, isoformLengths, segmentCounts, Thetas):
    isoformTPMs = np.zeros(shape=(len(geneIsoforms)))
    isoformCounts = np.zeros(shape=(len(geneIsoforms)))
    for i in range(len(segmentIDs)):
        totalProbSegment = 0
        for j in range(len(geneIsoforms)):
            if(segmentIsoformIndicatorMatrix[i,j]):
                totalProbSegment += Thetas[j]
        for j in range(len(geneIsoforms)):
            if(segmentIsoformIndicatorMatrix[i,j] and totalProbSegment != 0):
                isoformTPMs[j] += segmentCounts[i] * (Thetas[j] / totalProbSegment)
                isoformCounts[j] += segmentCounts[i] * (Thetas[j] / totalProbSegment)

    isoformTPMs = np.divide(isoformTPMs, isoformLengths)
    return isoformTPMs, isoformCounts
###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

def computeAlphaCost(Alpha, segmentCounts, readCount, segmentIsoformIndicatorMatrix, effectiveSegmentLengths, isoformLengths):
    cost = 0
    for i in range(len(segmentCounts)):
        totalSegmentAlpha = np.sum(segmentIsoformIndicatorMatrix[i, :] * Alpha)
        for j in range(len(isoformLengths)):
            if segmentIsoformIndicatorMatrix[i, j] == 1:
                observedValue = segmentCounts[i] * Alpha[j] / totalSegmentAlpha
                meanValue = readCount * Alpha[j] * effectiveSegmentLengths[i] / isoformLengths[j]
                #print("Observed Value: ", observedValue)
                #print("mean Value: ", meanValue)
                diffValue = (observedValue - meanValue)
                cost += (diffValue*diffValue)
    print("Cost: ", cost)
    return cost

def calculateLambda(readCount, effectiveSegmentLengths, isoformLengths):
    Lambda = np.zeros(size=(len(effectiveSegmentLengths), len(isoformLengths)))
    for i in range(len(effectiveSegmentLengths)):
        for j in range(len(isoformLengths)):
            Lambda = effectiveSegmentLengths[i] / isoformLengths[j] * readCount

def computeGradientAlpha(Alpha, segmentCounts, segmentIsoformIndicatorMatrix, Lambda):
    alphaGradients = np.zeros(size=len(Alpha))
    for i in range(len(segmentCounts)):
        totalSegmentAlpha = np.sum(segmentIsoformIndicatorMatrix[i, :] * Alpha)
        for j in range(len(Alpha)):
            if segmentIsoformIndicatorMatrix[i, j] == 1:
                alphaGradients[j] += 2*((Alpha[j] * (segmentCounts[i] / totalSegmentAlpha - Lambda[i, j])) * 
                    (segmentCounts[i] * (totalSegmentAlpha - Alpha)/(totalSegmentAlpha*totalSegmentAlpha) - Lambda[i, j]))
    return alphaGradients

def computeIsoformAbundances(pairedEndMode, countFileLines, readLength, fragmentLength, diffMax):
    geneCount = 0
    outputDict = dict()

    sumTPM = 0

    countFileLineIndexStart = 0
    countFileLineEnd = len(countFileLines)

    countFileLineIndex = countFileLineIndexStart
    if len(countFileLines) != 0:
        countFileLine = countFileLines[countFileLineIndex].strip()
    else:
        countFileLine = ""

    normSegCountsPlotList = []
    segmentLengthPlotList = []

    while countFileLine != "" and countFileLineIndex < countFileLineEnd:
        segmentIDs = []
        segmentCountsDict = dict()
        segmentCounts = []
        segmentLengths = []
        segmentIsoforms = dict()
        geneIsoforms = []

        splitLine = countFileLine.split("\t")

        if pairedEndMode:
            currGene = splitLine[4] #Get Current Gene
            while countFileLine != "" and countFileLineIndex < countFileLineEnd and splitLine[4] == currGene:
                segmentID1 = splitLine[0]
                segmentID2 = splitLine[1]
                if segmentID1 == segmentID2:
                    segmentID = "E_" + segmentID1 + "_" + segmentID2
                    segmentLengths.append(int(splitLine[5]))
                else:
                    segmentID = "J_" + segmentID1 + "_" + segmentID2
                    segmentLengths.append(fragmentLength - readLength)
                segmentIsoforms[segmentID] = splitLine[9].split(",")
                segmentIDs.append(segmentID)
                for isoform in segmentIsoforms[segmentID]:
                    if isoform not in geneIsoforms:
                        geneIsoforms.append(isoform)
                segmentCountsDict[segmentID] = int(splitLine[2])

                countFileLineIndex += 1
                if countFileLineIndex < countFileLineEnd:
                    countFileLine = countFileLines[countFileLineIndex].strip()
                    splitLine = countFileLine.split("\t")
        else:
            currGene = splitLine[3] #Get Current Gene
            while countFileLine != "" and countFileLineIndex < countFileLineEnd and splitLine[3] == currGene:
                segmentID = "E_" + splitLine[0]
                segmentIDs.append(segmentID)
                segmentCountsDict[segmentID] = int(splitLine[1])
                segmentLengths.append(int(splitLine[4]))

                segmentIsoforms[segmentID] = splitLine[6].split(",")
                for isoform in segmentIsoforms[segmentID]:
                    if isoform not in geneIsoforms:
                        geneIsoforms.append(isoform)

                countFileLineIndex += 1
                if countFileLineIndex < countFileLineEnd:
                    countFileLine = countFileLines[countFileLineIndex].strip()
                    splitLine = countFileLine.split("\t")

        for segmentID in segmentIDs: #Only get segment counts of certain segments
            segmentCounts.append(segmentCountsDict[segmentID])
        readCount = sum(segmentCounts)

        segmentIsoformIndicatorMatrix = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)), dtype="bool_")
        for i in range(len(segmentIDs)):
            segmentID = segmentIDs[i]
            segmentIsoformList = segmentIsoforms[segmentID]
            for isoform in segmentIsoformList:
                j = geneIsoforms.index(isoform)
                segmentIsoformIndicatorMatrix[i, j] = True

        effectiveSegmentLengths = np.zeros(shape=(len(segmentIDs)), dtype="int_")
        for i in range(len(segmentIDs)):
            if pairedEndMode:
                if (segmentIDs[i][0] == "E"):
                    effectiveSegmentLengths[i] = max(segmentLengths[i] - fragmentLength + 1, 1)
                else:
                    effectiveSegmentLengths[i] = segmentLengths[i]
            else:
                effectiveSegmentLengths[i] = max(segmentLengths[i] - readLength + 1, 1)

        isoformCounts = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
        isoformLengths = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
        for j in range(len(geneIsoforms)):
            for i in range(len(segmentIDs)):
                if(segmentIsoformIndicatorMatrix[i,j] and segmentIDs[i][0] == "E"):
                    isoformCounts[j] += segmentCounts[i]
                    isoformLengths[j] += segmentLengths[i] - readLength + 1
            isoformLengths[j] += (readLength - fragmentLength) #Make effective transcript length L - fragmentLength + 1
            #isoformLengths[j] = max(isoformLengths[j], 1)
        
        ############################################################################################################################################
        ## Find H for each segment
        ############################################################################################################################################
        
        segmentH = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
        for i in range(len(segmentIDs)):
            for j in range(len(geneIsoforms)):
                if(segmentIsoformIndicatorMatrix[i,j]):
                    segmentH[i, j] = 1.0/(isoformLengths[j])#1.0/isoformCounts[j]#float(segmentCounts[i])/(isoformCounts[j] * effectiveSegmentLengths[i])

        print("Time to Analyze Distribution", time.time() - startTime)

        for i in range(len(segmentIDs)):
            for j in range(len(geneIsoforms)):
                if(segmentIsoformIndicatorMatrix[i,j] and isoformCounts[j] != 0):
                    normSegCountsPlotList.append(segmentCounts[i]/float(isoformCounts[j]))
                    segmentLengthPlotList.append(log(1+effectiveSegmentLengths[i]))
    #     #####################################################################################################
    #     ## EM algorithm
    #     #####################################################################################################

        Lambda = calculateLambda(readCount, effectiveSegmentLengths, isoformLengths)

        Thetas = np.zeros(shape=len(geneIsoforms))
        Alpha = np.zeros(shape=len(geneIsoforms)) #This is Theta with ~ over it
        
        isoformNumBases = np.zeros(shape=len(geneIsoforms))
        isoformUniqueCountSum = np.zeros(shape=len(geneIsoforms))
        #If segment belongs to only one isoform, change value to isoform index
        segmentUniqueIsoform = np.full(shape=len(segmentIDs), fill_value=-1)
        for i in range(len(segmentIDs)):
            if np.sum(segmentIsoformIndicatorMatrix[i,:]) == 1: #If segment has only one isoform
                segmentUniqueIsoform[i] = np.where(segmentIsoformIndicatorMatrix[i,:] == 1)[0]
        for i in range(len(segmentIDs)):
            if segmentUniqueIsoform[i] != -1:
                isoformNumBases[segmentUniqueIsoform[i]] += effectiveSegmentLengths[i]
                isoformUniqueCountSum[segmentUniqueIsoform[i]] += segmentCounts[i]
        numIsoformsUniqueSegments = 0
        for j in range(len(geneIsoforms)):
            if isoformUniqueCountSum[j] != 0:
                Thetas[j] = float(isoformUniqueCountSum[j])/isoformNumBases[j]
                numIsoformsUniqueSegments += 1
        sumTheta = np.sum(Thetas)
        if sumTheta != 0:
            Thetas = np.divide(Thetas, sumTheta)
        Thetas = np.multiply(Thetas, float(numIsoformsUniqueSegments)/len(geneIsoforms))
        for j in range(len(geneIsoforms)):
            if isoformUniqueCountSum[j] == 0:
                Thetas[j] = 1.0/len(geneIsoforms)
        for j in range(len(geneIsoforms)):
            Alpha[j] = Thetas[j] * isoformLengths[j]
        sumAlpha = np.sum(Alpha)
        for j in range(len(geneIsoforms)):
            Alpha[j] /= sumAlpha
        oldAlpha = np.zeros(shape=len(geneIsoforms))

         #########################################################################################################
         ## iteration begins        
                
        if readCount != 0:
            numInitialThetas = 1
            ThetasLikelihoods = []
            for initialThetaIndex in range(numInitialThetas):

                diff = 1.0
                iterCount = 0
                likelihoodMatrix = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
                while diff > diffMax:
                    #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+str(diff)+"\t"+str(tmpTime))

                    #Calculate H
                    # for j in range(len(geneIsoforms)):
                    #     x = []
                    #     y = []
                    #     for i in range(len(segmentIDs)):
                    #         if(segmentIsoformIndicatorMatrix[i,j] == 1):
                    #             ThetaSegmentSum = np.sum(np.multiply(segmentIsoformIndicatorMatrix[i,:], Thetas))
                    #             x.append(effectiveSegmentLengths[i])
                    #             y.append(Thetas[j]/ThetaSegmentSum*segmentCounts[i])
                    #     x = np.reshape(np.array(x), newshape=(-1,1))
                    #     regr = linear_model.LinearRegression()
                    #     regr.fit(x, y)
                    #     yPred = regr.predict(x)
                    #     yPredTotal = sum(yPred)
                    #     yIndex = 0
                    #     for i in range(len(segmentIDs)):
                    #         if(segmentIsoformIndicatorMatrix[i,j] == 1):
                    #             segmentH[i,j] = yPred[yIndex]/yPredTotal
                    #             yIndex += 1

                    #Expectation Step
                    for i in range(len(segmentIDs)):
                        probNumerators = np.multiply(Alpha, segmentH[i])
                        probDenominator = np.sum(probNumerators)
                        if probDenominator != 0:
                            likelihoodMatrix[i] = np.multiply(probNumerators, float(segmentCounts[i]) / probDenominator)
                        else:
                            likelihoodMatrix[i] = np.multiply(probNumerators, 0)  

                    #Maximization Step
                    oldAlpha = Alpha
                    Alpha = np.multiply(np.sum(likelihoodMatrix, axis=0), 1/float(readCount))
                    
                    diffArray = np.absolute(np.subtract(Alpha, oldAlpha))
                    diff = diffArray.max()

                    sumAlpha = np.sum(Alpha)
        
                    Alpha = np.multiply(Alpha, 1.0 / sumAlpha)

                    Thetas = np.zeros(shape=(len(geneIsoforms),))
                    for j in range(len(Alpha)):
                        Thetas[j] = Alpha[j] / isoformLengths[j]
                    sumTheta = np.sum(Thetas)
                    Thetas = np.divide(Thetas, sumTheta)

                    iterCount += 1

                logliklihood = np.sum(np.log(likelihoodMatrix+1))
                ThetasLikelihoods.append([logliklihood, Thetas])
                previousAlpha = oldAlpha

                #New Thetas
                Thetas = np.random.dirichlet(np.ones(len(geneIsoforms)), size=1)[0]
                for j in range(len(geneIsoforms)):
                    Alpha[j] = Thetas[j] * isoformLengths[j]
                sumAlpha = np.sum(Alpha)
                for j in range(len(geneIsoforms)):
                    Alpha[j] /= sumAlpha
                oldAlpha = np.zeros(shape=len(geneIsoforms))

            #Ensure no duplicates
            uniqueLikelihoods = []
            uniqueThetasLikelihoods = []
            for ThetaLikelihood in ThetasLikelihoods:
                if ThetaLikelihood[0] not in uniqueLikelihoods:
                    uniqueLikelihoods.append(ThetaLikelihood[0])
                    uniqueThetasLikelihoods.append(ThetaLikelihood)

            maxLikelihoodThetas = max(uniqueThetasLikelihoods)[1]
            print(str(geneCount) + ": " + str(currGene)+"\t"+str(iterCount)+" iterations\tDone!")

            computeAlphaCost(previousAlpha, segmentCounts, readCount, segmentIsoformIndicatorMatrix, effectiveSegmentLengths, isoformLengths)
            print("Alpha True: ")
            AlphaReal = [0.34727334,  0.35028438,  0.04717297,  0.0267648,  0.04750753,  0.03178321, 0.06557377,  0.08364001]
            computeAlphaCost(AlphaReal, segmentCounts, readCount, segmentIsoformIndicatorMatrix, effectiveSegmentLengths, isoformLengths)
            isoformTPMs, isoformCounts = calculateIsoformTPMs(segmentIDs, geneIsoforms, segmentIsoformIndicatorMatrix, isoformLengths, segmentCounts, maxLikelihoodThetas)
        else:
            isoformTPMs = np.zeros(len(geneIsoforms))
            isoformCounts = np.zeros(len(geneIsoforms))
            maxLikelihoodThetas = np.full(len(geneIsoforms), fill_value=1.0/len(geneIsoforms))

        print("Isoforms ", geneIsoforms)
        print("Alphas: ", previousAlpha)

        effectiveTranscriptSegmentLengths = [[] for j in range(len(geneIsoforms))]
        for j in range(len(geneIsoforms)):
            for i in range(len(effectiveSegmentLengths)):
                if(segmentIsoformIndicatorMatrix[i,j]):
                    effectiveTranscriptSegmentLengths[j].append(segmentLengths[i])

        outputDict[currGene] = (geneIsoforms, isoformTPMs, isoformCounts, maxLikelihoodThetas, effectiveTranscriptSegmentLengths) 
        sumTPM += np.sum(isoformTPMs)
        
        #print("Output Dict: ", outputDict)
        #for j in range(len(Alpha)):
            #print(currGene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+isoformNames[i]+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(tmpTime))
            
        print("Time for EM", time.time() - startTime)
        geneCount += 1
        
    #Plot #1
    plt.scatter(segmentLengthPlotList, normSegCountsPlotList)
    plt.xlabel(r"$log(l_{eff})$")
    plt.ylabel(r"$\frac{Seg Count}{Total_t}$")
    plt.savefig("YanagiQuantFigures/figure_1_1")
    return sumTPM, outputDict

def arrayToString(array):
    arrayString = ""
    if len(array) != 0:
        for i in range(len(array)-1):
            arrayString += str(array[i]) + ","
        arrayString += str(array[len(array)-1])
    return arrayString

if __name__ == "__main__":
    countFile, readLength, fragmentLength, diffMax, numProcesses, outFile = getArgs()

    countFileHandle = open(countFile, 'rU')
    OUT = open(outFile, 'w')
    OUT.write("IsoformName\tGeneName\tCount\tTPMs\tRelativeAbundance\tSegmentLengths\n") ## Header of Results

    startTime = time.time()

    pairedEndMode = detectPairedEnd(countFileHandle, fragmentLength)

    if pairedEndMode:
        fragmentLength = calculatePairedFragmentLength(countFile, readLength)

    countFileLines = countFileHandle.readlines()
    countFileLength = len(countFileLines)
    
    countFileLineStarts, countFileLineEnds = partitionCountFile(countFileLines, countFileLength, pairedEndMode, numProcesses)

    pool = Pool(processes=numProcesses)
    #results = [pool.apply_async(computeIsoformAbundances, args=(countFileLines[countFileLineStarts[x]: countFileLineEnds[x]], readLength, diffMax,)) for x in range(numProcesses)]
    results = [computeIsoformAbundances(pairedEndMode, countFileLines[countFileLineStarts[0]:countFileLineEnds[0]], readLength, fragmentLength, diffMax,)]

    totalOutputDict = {}
    totalSumTPM = 0
    for result in results:
        sumTPM, outputDict = result#result.get()
        totalSumTPM += sumTPM
        totalOutputDict.update(outputDict) 

    sortedGeneIDs = sorted(totalOutputDict.keys())
    for geneID in sortedGeneIDs:
        geneOutput = totalOutputDict[geneID]
        for j in range(len(geneOutput[0])):
            transcriptID = geneOutput[0][j]
            isoformTPM = geneOutput[1][j]/totalSumTPM *1e6
            isoformCounts = geneOutput[2][j]
            Theta = geneOutput[3][j]
            effectiveSegmentLengthsString = arrayToString(geneOutput[4][j])
            OUT.write(transcriptID+"\t"+geneID+"\t"+str(isoformCounts)+"\t"+str(isoformTPM)+"\t"+str(Theta)+"\t"+effectiveSegmentLengthsString+"\n")

    countFileHandle.close()
    OUT.close()

    print(time.time() - startTime)