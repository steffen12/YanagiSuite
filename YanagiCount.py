import time
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
import itertools
import sys, pickle
import os

#Todo:
#1) Paired end resolve multimapping
#2) paired end frag
#3) Decide multimapped
#5) Check if paired end counts are valid
#from InputLoader import *

#python YanagiCount.py --single --stdin1 segments_ENSG00000100842.fa YanagiCountOutput 100
#python YanagiCount.py --paired --stdin1 segments_ENSG00000100842.fa YanagiCountOutput 100
#python YanagiCount.py --paired --stdin2 segments_ENSG00000100842.fa YanagiCountOutput 100

#python YanagiCount.py --single --bam segments_ENSG00000100842.fa SamOutput.sam YanagiCountOutput 100
#python YanagiCount.py --paired --bam segments_ENSG00000100842.fa SamOutput.sam SamOutput.sam YanagiCountOutput 100

#python ~/Penn_Internship/YanagiSuite/YanagiCount.py --single --kallisto ~/Penn_Internship/YanagiSuite/segments_ENSG00000100842.fa ~/Penn_Internship/YanagiSuite/kallisto_output_segments_ENSG00000100842 ~/Penn_Internship/YanagiSuite/YanagiCountOutput 100
#python ~/Penn_Internship/YanagiSuite/YanagiCount.py --paired --kallisto ~/Penn_Internship/YanagiSuite/segments_ENSG00000100842.fa ~/Penn_Internship/YanagiSuite/kallisto_output_segments_ENSG00000100842 ~/Penn_Internship/YanagiSuite/kallisto_output_segments_ENSG00000100842_2 ~/Penn_Internship/YanagiSuite/YanagiCountOutput 100

class Seg:
    ID = ""
    name = ""
    header = ""
    txs = set()

    length = 0
    segtype = "E"
    geneID = ""
    startLoc = 0
    strand = "+"

    def __init__(self, header):
        self.header = header

        ID, name, txs, segtype = header.split(" ")
        self.ID = ID
        self.name = name
        
        tokens = name.split(":")
        self.geneID = tokens[1]
        self.startLoc = long(tokens[2])
        self.strand = tokens[-1]
        
        self.segtype = segtype.split(":")[-1]
        self.txs = set(txs.split(":")[1].split(","))

    def __str__(self):
        return "%s:<%s>:%d" % (self.ID, self.name, self.length)
    def __repr__(self):
        return self.__str__()
        

class ReadAlignments:
    readID = ""
    readLen = 0
    aligns = [[], []]

    def __init__(self, readID, readLen):
        self.readID = readID
        self.readLen = readLen
        self.aligns = [[], []]

    def addAlignment(self, seg, flags, end_idx):
        self.aligns[end_idx].append((seg, flags))

    def __str__(self):
        return "%s:%d:[%s,%s]" % (self.readID, self.readLen, str(self.aligns[0]), str(self.aligns[1]))
    def __repr__(self):
        return self.__str__()
        
def load_SegmentsLib(seg_file):
    segs_dict = {}
    segIDs = []
    for record in SeqIO.parse(seg_file, "fasta"):
        seg = Seg(record.description)
        segs_dict[seg.ID] = seg
        segIDs.append(seg.ID)
        seg.length = len(record.seq)
    return(segs_dict, segIDs)

def usageExit():
    print("Usage: ")
    print("--single or --paired for reads")
    print("--stdin1 or --stdin2 or --bam or --kallisto for input mode")
    print("segment reference file")
    print("bam files (if --bam mode) or kallisto output directories (if --kallisto mode)")
    print("output directory")
    print("process count")
    sys.exit()

def read_alignment_SAM(f, read_aligns, readIDs, end_idx):  
    line = f.readline()
    #print(line)
    if not line:
        return ("", False) #Return when file is done
    if line[0] == "@": #If header line:
        return ("@", False)
    tokens = line.split("\t")
    readID = tokens[0]
    segID = tokens[2]

    flags = int(tokens[1])

    if segID == "*": #"*" represents a N/A
        return(readID, False)

    # Append to readIDs list
    if not readIDs[end_idx] or readID != readIDs[end_idx][-1]: #end_idx is file
        readIDs[end_idx].append(readID)
    
    pos = int(tokens[3])
    readLen = len(tokens[9])
    if readID in read_aligns:
        #if flags <= 16:
        read_aligns[readID].addAlignment(segID, flags, end_idx)
    else:
        align = ReadAlignments(readID, readLen)
        #if flags <= 16:
        align.addAlignment(segID, flags, end_idx)
        read_aligns[readID] = align
    return(readID, True)
    
def read_alignment_Kallisto(kallistoECFile, kallistoECReadsFile, read_aligns, totalReadIDs, segIDs, end_idx):  
    #end_idx is index of paired-end
    kallistoECFileLine = kallistoECFile.readline().strip()
    kallistoECReadsFileLine = kallistoECReadsFile.readline().strip()
    if not kallistoECFileLine and not kallistoECReadsFileLine:
        return (True) #Return when file is done
    if (not kallistoECFileLine and kallistoECReadsFileLine) or (kallistoECFileLine and not kallistoECReadsFileLine):
        print("Malformed Kallisto Files")
        usageExit()
    
    tokens = kallistoECReadsFileLine.split("\t")
    if len(tokens) == 1:
        ecID1 = tokens[0]
        readIDs = []
    else:
        ecID1 = tokens[0]
        readIDs = tokens[1].split(",")[:-1]

    tokens = kallistoECFileLine.split("\t")
    ecID2 = tokens[0]
    segmentIndexes = tokens[1].split(",")

    if ecID1 != ecID2:
        print("Malformed Kallisto Files")
        usageExit()

    # Append to totalReadIDs list
    
    if len(segmentIndexes) == 1:
        for readID in readIDs:
            if not totalReadIDs[end_idx] or readID != totalReadIDs[end_idx][-1]: #end_idx is file
                totalReadIDs[end_idx].append(readID)
        
            if readID not in read_aligns:
                align = ReadAlignments(readID, -1)
                read_aligns[readID] = align

            for segmentIndex in segmentIndexes:
                segID = segIDs[int(segmentIndex)]
                read_aligns[readID].addAlignment(segID, -1, end_idx)
    return(False)

def validOrientation(flags1, flags2):
    flag1Bin = bin(flags1)[2:]
    flag2Bin = bin(flags2)[2:]
    if(len(flag1Bin) < 5):
        flags1Orientation = 0
    else:
        flags1Orientation = 16*int(flag1Bin[-5])
    if(len(flag2Bin) < 5):
        flags2Orientation = 0
    else:
        flags2Orientation = 16*int(flag2Bin[-5])
    validOrientationTrue = ((flags1Orientation + flags2Orientation) == 16)
    return validOrientationTrue

def saveObj(obj, name):
    pickle.dump(obj, open(name+".p", "wb"))

def loadObj(name):
    return pickle.load(open(name+".p", "rb"))

def processSingleReads(readIDs, read_aligns, segs_dict, segID_counts):
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    c=0
    for readID in readIDs:
        aligns = read_aligns.pop(readID)
        mapped = False
        passed_txs = False
        c += len(aligns.aligns[0])
        for align in aligns.aligns[0]:
            segID = align[0]
            txs = segs_dict[segID].txs
            segID_counts[segID] += 1
            mapped = True
            passed_txs = True
            break
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    #print str(float(c)/len(readIDs)) + " alignments per read"
    return readsMapped

def mainSingleEnd(stdinTrue, stdinPair, pairedEndMode, kallistoMode):
    segmentReferenceFilename = sys.argv[3]  # path to the segments.fa reference 
    if stdinTrue:
        outputDir = sys.argv[4]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[5])    # How many reads alignments to store before processing 
    elif kallistoMode:
        kallistoOutputDir1 = sys.argv[4] 
        outputDir = sys.argv[5]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[6])    # How many reads alignments to store before processing 
        kallistoECFilename1 = kallistoOutputDir1 + "/pseudoalignments.ec"
        kallistoECReadsFilename1 = kallistoOutputDir1 + "/pseudoalignments.tsv"
    else:
        BAM1 = sys.argv[4]  
        outputDir = sys.argv[5]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[6])    # How many reads alignments to store before processing 
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputCountsFilename = outputDir + "/seg_counts.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    #DExons = load_disjointExons(seg_inDir)
    #txs2exs = load_Txs2Exs(seg_inDir, DExons)
    segs_dict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print "Done!"
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Reading Read alignments..."
    start_t = time.time()
    read_aligns = {}
    segID_counts = dict()
    for segID in segIDs:
        segID_counts[segID] = 0
    #TODO clean up based on sorted readIDs
    readsMapped = [0,0,0]
    if stdinTrue:
        f1 = sys.stdin
    elif kallistoMode:
        kallistoECFile1 = open(kallistoECFilename1)
        kallistoECReadsFile1 = open(kallistoECReadsFilename1)
    else:
        f1 = open(BAM1) 
    fout = open(outputCountsFilename, "w")
    done = False
    readIDs = [[], []]
    i = 0
    while not done:
        if not done: # not done with the first bam
            if kallistoMode:
                done = read_alignment_Kallisto(kallistoECFile1, kallistoECReadsFile1, read_aligns, readIDs, segIDs, 0)
            else:
                readID, aligned = read_alignment_SAM(f1, read_aligns, readIDs, 0)
                done = not readID

        # Process some records (if stored more than process_count reads in memory)
        if not (stdinTrue and pairedEndMode) and len(readIDs[0]) > process_count:
            process_readIDs = readIDs[0][:process_count]
            readIDs[0] = readIDs[0][process_count:]
            print i, len(read_aligns)
            start_t_2 = time.time()
            rMapped = processSingleReads(process_readIDs, read_aligns, segs_dict, segID_counts)
            readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
            elapsed = time.time() - start_t_2
            print "Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed
        i+=1
    if(stdinTrue and pairedEndMode):
        saveObj(read_aligns, "read_aligns")
    else:
        print i, len(read_aligns)
        # Process the rest of records
        process_readIDs = readIDs[0]
        rMapped = processSingleReads(process_readIDs, read_aligns, segs_dict, segID_counts)
        readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
        print "Processed:", len(process_readIDs), len(read_aligns), rMapped
        print "Done!"
        print "Mapped Reads:", readsMapped[0], "Unmapped Reads:", readsMapped[1]+readsMapped[2], \
              "Unmapped Due Txs:", readsMapped[1], "Unmapped Due Direction:", readsMapped[2]
        elapsed = time.time() - start_t
        print "Elapsed Time: ", elapsed

        print "Writing Segments Counts..."
        start_t = time.time()
        with open(outputCountsFilename, "w") as f:
            f.write("SEGID\tcount\tSEGTYPES\tGENE\tSEGLEN\tSEGSLoc\tTXS\n")
            
            for segID in sorted(segID_counts.iterkeys()):
                count = segID_counts[segID]
                seg = segs_dict[segID]
                types = seg.segtype
                line = "\t".join([seg.ID, str(count), types,
                                  seg.geneID,
                                  str(seg.length),
                                  str(seg.startLoc), ','.join(seg.txs)])
                f.write(line + "\n")
        elapsed = time.time() - start_t
        print "Elapsed Time: ", elapsed

#@profile
def processPairedReads(readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs):
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    c=0
    for readID in readIDs:
        aligns = read_aligns.pop(readID)
        mapped = False
        passed_txs = False
        c += len(aligns.aligns[0]+aligns.aligns[1])
        for align in itertools.product(aligns.aligns[0], aligns.aligns[1]):
            segID1, segID2 = align[0][0], align[1][0]
            pairkey = segID1 + "_" + segID2
            if pairkey in segPairs_txs:
                txs = segPairs_txs[pairkey]
            else:
                seg1 = segs_dict[segID1]
                seg2 = segs_dict[segID2]
                txs = (seg1.txs & seg2.txs)
                segPairs_txs[pairkey] = txs
            if kallistoMode:
                if len(txs) < 1:    # alignment between segs with no tx in common
                    newSegPairs_counts[pairkey] += 1    # Counted as a novel junction
                else:
                    if pairkey in segPairs_counts:
                        segPairs_counts[pairkey] += 1
                    else:
                        segPairs_counts[pairkey] = 1
                    mapped = True
            else:
                if len(txs) < 1 and validOrientation(align[0][1], align[1][1]): # alignment between segs with no tx in common
                    newsegPairs_counts[pairkey] += 1    # Counted as a novel junction
                elif not validOrientation(align[0][1], align[1][1]):  # not valid alignment
                    passed_txs = True
                else:
                    if pairkey in segPairs_counts:
                        segPairs_counts[pairkey] += 1
                    else:
                        segPairs_counts[pairkey] = 1
                    mapped = True
                    break
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print str(float(c)/len(readIDs)) + " alignments per read"
    return readsMapped

#@profile
def mainPairedEnd(stdinTrue, stdinPair, kallistoMode):
    segmentReferenceFilename = sys.argv[3]  # path to the segments.fa reference 
    if stdinTrue:
        outputDir = sys.argv[4]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[5])    # How many reads alignments to store before processing 
    elif kallistoMode:
        kallistoOutputDir1 = sys.argv[4]  
        kallistoOutputDir2 = sys.argv[5]
        outputDir = sys.argv[6]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[7])    # How many reads alignments to store before processing 
        kallistoECFilename1 = kallistoOutputDir1 + "/pseudoalignments.ec"
        kallistoECReadsFilename1 = kallistoOutputDir1 + "/pseudoalignments.tsv"
        kallistoECFilename2 = kallistoOutputDir2 + "/pseudoalignments.ec"
        kallistoECReadsFilename2 = kallistoOutputDir2 + "/pseudoalignments.tsv"
    else:
        BAM1 = sys.argv[4]  
        BAM2 = sys.argv[5]
        outputDir = sys.argv[6]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[7])    # How many reads alignments to store before processing 

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputCountsFilename = outputDir + "/seg_counts.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    #DExons = load_disjointExons(seg_inDir)
    #txs2exs = load_Txs2Exs(seg_inDir, DExons)
    segs_dict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print "Done!"
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Reading Read alignments..."
    start_t = time.time()
    read_aligns = {}
    segPairs_counts = dict()
    for segID in segIDs:
        segPairs_counts[segID+"_"+segID] = 0
    newsegPairs_counts = Counter()
    segPairs_txs = dict()
    for segID in segIDs:
        segPairs_txs[segID+"_"+segID] = segs_dict[segID].txs
    #TODO clean up based on sorted readIDs
    readsMapped = [0,0,0]
    if stdinTrue:
        read_aligns = loadObj("read_aligns")
        f2 = sys.stdin
        done = [True, False]
    elif kallistoMode:
        kallistoECFile1 = open(kallistoECFilename1)
        kallistoECReadsFile1 = open(kallistoECReadsFilename1)
        kallistoECFile2 = open(kallistoECFilename2)
        kallistoECReadsFile2 = open(kallistoECReadsFilename2)
        done = [False, False]
    else:
        f1 = open(BAM1)
        f2 = open(BAM2)
        done = [False, False]
    fout = open(outputCountsFilename, "w")
    readIDs = [[], []]
    i = 0
    while not all(done):
        if not done[0]: # not done with the first bam
            if kallistoMode:
                done[0] = read_alignment_Kallisto(kallistoECFile1, kallistoECReadsFile1, read_aligns, readIDs, segIDs, 0)
            else:
                readID, aligned = read_alignment_SAM(f1, read_aligns, readIDs, 0)
                done[0] = not readID
        if not done[1]: # not done with the second bam
            if kallistoMode:
                done[1] = read_alignment_Kallisto(kallistoECFile2, kallistoECReadsFile2, read_aligns, readIDs, segIDs, 1)
            else:
                readID, aligned = read_alignment_SAM(f2, read_aligns, readIDs, 1)
                done[1] = not readID

        # Process some records (if stored more than process_count reads in memory)
        if not stdinTrue:
            if len(readIDs[0]) > process_count and len(readIDs[1]) > process_count:
                process_readIDs = readIDs[0][:process_count]
                readIDs[0] = readIDs[0][process_count:]
                readIDs[1] = readIDs[1][process_count:]
                print i, len(read_aligns)
                start_t_2 = time.time()
                rMapped = processPairedReads(process_readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
                elapsed = time.time() - start_t_2
                print "Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed
        elif kallistoMode:
            if i == process_count:
                process_readIDsSet = set(readIDs[0]) & set(readIDs[1])
                process_readIDs = list(process_readIDsSet)
                readIDs[0] = list(set(readIDs[0]) - process_readIDsSet)
                readIDs[1] = list(set(readIDs[1]) - process_readIDsSet)
                print i, len(read_aligns)
                start_t_2 = time.time()
                newReadsMapped = processPairedReads(process_readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newSegPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+newReadsMapped[0], readsMapped[1]+newReadsMapped[1], readsMapped[2]+newReadsMapped[2]]
                elapsed = time.time() - start_t_2
                i = 0
                print "Processed:", len(process_readIDs), len(read_aligns), newReadsMapped, elapsed
        else:
            if len(readIDs[1]) > process_count:
                process_readIDs = readIDs[1][:process_count]
                readIDs[1] = readIDs[1][process_count:]
                print i, len(read_aligns)
                start_t_2 = time.time()
                rMapped = processPairedReads(process_readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
                elapsed = time.time() - start_t_2
                print "Stdin Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed    
        i+=1
    print i, len(read_aligns)
    # Process the rest of records
    if stdinTrue:
        process_readIDs = readIDs[1]
        rMapped = processPairedReads(process_readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
    else:
        process_readIDs = readIDs[0]
        rMapped = processPairedReads(process_readIDs, kallistoMode, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
    readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
    print "Processed:", len(process_readIDs), len(read_aligns), rMapped
    print "Done!"
    print "Mapped Reads:", readsMapped[0], "Unmapped Reads:", readsMapped[1]+readsMapped[2], \
          "Unmapped Due Txs:", readsMapped[1], "Unmapped Due Direction:", readsMapped[2]
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Writing Segments Counts..."
    start_t = time.time()
    with open(outputCountsFilename, "w") as f:
        f.write("SEG1ID\tSEG2ID\tcount\tSEGTYPES\tGENE\tSEG1LEN\tSEG2LEN\tSEG1SLoc\tSEG2SLoc\tTXS\n")
        
        for segPair in sorted(segPairs_counts.iterkeys()):
            count = segPairs_counts[segPair]
            segs = [segs_dict[segID] for segID in segPair.split("_")]
            types = segs[0].segtype
            types += segs[1].segtype if segs[0].ID != segs[1].ID else ""
            line = "\t".join([segs[0].ID, segs[1].ID, str(count), types,
                              segs[0].geneID,
                              str(segs[0].length), str(segs[1].length),
                              str(segs[0].startLoc), str(segs[1].startLoc), ','.join(segPairs_txs[segPair])])
            f.write(line + "\n")
    with open(outputCountsFilename+".newJuncs", "w") as f:
        f.write("SEG1ID\tSEG2ID\tcount\tSEGTYPES\tSEG1GENE\tSEG2GENE\tSEG1LEN\tSEG2LEN\tSEG1SLoc\tSEG2SLoc\n")
        
        for segPair in sorted(newsegPairs_counts.iterkeys()):
            count = newsegPairs_counts[segPair]
            segs = [segs_dict[segID] for segID in segPair.split("_")]
            types = segs[0].segtype
            types += segs[1].segtype if segs[0].ID != segs[1].ID else ""
            line = "\t".join([segs[0].ID, segs[1].ID, str(count), types,
                              segs[0].geneID, segs[1].geneID,
                              str(segs[0].length), str(segs[1].length),
                              str(segs[0].startLoc), str(segs[1].startLoc)])
            f.write(line + "\n")

    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed


if __name__ == "__main__":
    pairedEndMode = False
    if sys.argv[1] == "--single":
        pairedEndMode = False
    elif sys.argv[1] == "--paired":
        pairedEndMode = True
    else:
        print("Invalid single/paired end parameter")
        usageExit()

    stdinTrue = False
    kallistoMode = False
    stdinPair = 0
    if sys.argv[2] == "--stdin1":
        stdinTrue = True
        kallistoMode = False
        stdinPair = 1
    elif sys.argv[2] == "--stdin2":
        stdinTrue = True
        kallistoMode = False
        stdinPair = 2
    elif sys.argv[2] == "--bam":
        stdinTrue = False
        kallistoMode = False
        stdinPair = 0
    elif sys.argv[2] == "--kallisto":
        stdinTrue = False
        kallistoMode = True
        stdinPair = 0
    else:
        print("Invalid input mode parameter")
        usageExit()

    if pairedEndMode and (not stdinTrue or (stdinTrue and stdinPair == 2)):
        mainPairedEnd(stdinTrue, stdinPair, kallistoMode)
    else:
        mainSingleEnd(stdinTrue, stdinPair, pairedEndMode, kallistoMode)
    
