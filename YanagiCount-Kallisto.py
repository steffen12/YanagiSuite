import time
from collections import defaultdict
from Bio import SeqIO 
from collections import Counter
import itertools
import sys

#python ~/Penn_Internship/PennSeq/YanagiCount-Kallisto.py --single ~/Penn_Internship/PennSeq/segments_ENSG00000100842.fa ~/Penn_Internship/PennSeq/kallisto_output_segments_ENSG00000100842 ~/Penn_Internship/PennSeq/YanagiQuant 100
#python ~/Penn_Internship/PennSeq/YanagiCount-Kallisto.py --paired ~/Penn_Internship/PennSeq/segments_ENSG00000100842.fa ~/Penn_Internship/PennSeq/kallisto_output_segments_ENSG00000100842 ~/Penn_Internship/PennSeq/kallisto_output_segments_ENSG00000100842_2 ~/Penn_Internship/PennSeq/YanagiQuant 10000

#Todo:
#1) Stdin, paired-end
#2) paired end frag
#3) Fix load segmen read FASTA
#4) Ask about break
#5) Check valid alignment
#6) process read ids
#7) break at end of processReads
#8) process eff
#from InputLoader import *

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
        return "%s:%d:[%s,%d]" % (self.readID, self.readLen, self.aligns[0], self.aligns[1])
    def __repr__(self):
        return self.__str__()
        
def usageExit():
    print("Usage: ")
    print("--single or --paired for reads")
    print("segment reference file")
    print("kallisto output directories")
    print("output directory")
    print("process count")
    sys.exit()

def load_SegmentsLib(seg_file):
    segs_dict = {}
    segIDs = []
    for record in SeqIO.parse(seg_file, "fasta"):
        seg = Seg(record.description)
        segs_dict[seg.ID] = seg
        segIDs.append(seg.ID)
        seg.length = len(record.seq)
    return(segs_dict, segIDs)

def read_alignment(kallistoECFile, kallistoECReadsFile, read_aligns, totalReadIDs, segIDs, end_idx):  
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
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print str(float(c)/len(readIDs)) + " alignments per read"
    return readsMapped

def mainSingleEnd(segmentReferenceFilename, kallistoOutputDir1, outputDir, process_count):
    outputCountsFilename = outputDir + "/seg_counts.tsv"
    kallistoECFilename1 = kallistoOutputDir1 + "/pseudoalignments.ec"
    kallistoECReadsFilename1 = kallistoOutputDir1 + "/pseudoalignments.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    segs_dict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print "Done!"
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Estimate Fragment Lengths..."
    start_t = time.time()
    fl_params = [181, 62.01, 1]

    print "Reading Read alignments..."
    start_t = time.time()
    read_aligns = {}
    segID_counts = Counter()
    readsMapped = [0,0,0]
    
    kallistoECFile1 = open(kallistoECFilename1)
    kallistoECReadsFile1 = open(kallistoECReadsFilename1)
    outputCountsFile = open(outputCountsFilename, "w")
    
    done = False
    totalReadIDs = [[]]
    i = 0
    while not done:
        if not done: # not done with the first file
            done = read_alignment(kallistoECFile1, kallistoECReadsFile1, read_aligns, totalReadIDs, segIDs, 0)        
        
        # Process some records (if stored more than process_count reads in memory)
        if len(totalReadIDs[0]) > process_count:
            process_readIDs = totalReadIDs[0][:process_count]
            totalReadIDs[0] = totalReadIDs[0][process_count:]
            #print i, len(read_aligns)
            start_t_2 = time.time()
            newReadsMapped = processSingleReads(process_readIDs, read_aligns, segs_dict, segID_counts)
            readsMapped = [readsMapped[0]+newReadsMapped[0], readsMapped[1]+newReadsMapped[1], readsMapped[2]+newReadsMapped[2]]
            elapsed = time.time() - start_t_2
            print "Processed:", len(process_readIDs), len(read_aligns), newReadsMapped, elapsed
        i+=1
    #print i, len(read_aligns)
    # Process the rest of records
    process_readIDs = totalReadIDs[0]
    newReadsMapped = processSingleReads(process_readIDs, read_aligns, segs_dict, segID_counts)
    readsMapped = [readsMapped[0]+newReadsMapped[0], readsMapped[1]+newReadsMapped[1], readsMapped[2]+newReadsMapped[2]]
    print "Processed:", len(process_readIDs), len(read_aligns), newReadsMapped
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
def processPairedReads(readIDs, read_aligns, segs_dict, segPairs_counts, newSegPairs_counts, segPairs_txs):
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
            if len(txs) < 1:    # alignment between segs with no tx in common
                newSegPairs_counts[pairkey] += 1    # Counted as a novel junction
            else:
                segPairs_counts[pairkey] += 1
                mapped = True
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print str(float(c)/len(readIDs)) + " alignments per read"
    return readsMapped

#@profile
def mainPairedEnd(segmentReferenceFilename, kallistoOutputDir1, kallistoOutputDir2, outputDir, process_count):
    outputCountsFilename = outputDir + "/seg_counts.tsv"
    kallistoECFilename1 = kallistoOutputDir1 + "/pseudoalignments.ec"
    kallistoECReadsFilename1 = kallistoOutputDir1 + "/pseudoalignments.tsv"
    kallistoECFilename2 = kallistoOutputDir2 + "/pseudoalignments.ec"
    kallistoECReadsFilename2 = kallistoOutputDir2 + "/pseudoalignments.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    segs_dict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print "Done!"
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Estimate Fragment Lengths..."
    start_t = time.time()
    fl_params = [181, 62.01, 1]

    print "Reading Read alignments..."
    start_t = time.time()
    read_aligns = {}
    segPairs_counts = Counter()
    newSegPairs_counts = Counter()
    segPairs_txs = {}
    readsMapped = [0,0,0]
    
    kallistoECFile1 = open(kallistoECFilename1)
    kallistoECReadsFile1 = open(kallistoECReadsFilename1)
    kallistoECFile2 = open(kallistoECFilename2)
    kallistoECReadsFile2 = open(kallistoECReadsFilename2)
    outputCountsFile = open(outputCountsFilename, "w")
    
    done = [False, False]
    totalReadIDs = [[], []]
    i = 0
    while not all(done):
        if not done[0]: # not done with the first file
            done[0] = read_alignment(kallistoECFile1, kallistoECReadsFile1, read_aligns, totalReadIDs, segIDs, 0)
        if not done[1]: # not done with the second file
            done[1] = read_alignment(kallistoECFile2, kallistoECReadsFile2, read_aligns, totalReadIDs, segIDs, 1)
        
        # Process some records (if stored more than process_count reads in memory)
        
        # Process process_count records
        if i == process_count:
            process_readIDsSet = set(totalReadIDs[0]) & set(totalReadIDs[1])
            process_readIDs = list(process_readIDsSet)
            print(process_readIDs)
            totalReadIDs[0] = list(set(totalReadIDs[0]) - process_readIDsSet)
            totalReadIDs[1] = list(set(totalReadIDs[1]) - process_readIDsSet)
            print i, len(read_aligns)
            start_t_2 = time.time()
            newReadsMapped = processPairedReads(process_readIDs, read_aligns, segs_dict, segPairs_counts, newSegPairs_counts, segPairs_txs)
            readsMapped = [readsMapped[0]+newReadsMapped[0], readsMapped[1]+newReadsMapped[1], readsMapped[2]+newReadsMapped[2]]
            elapsed = time.time() - start_t_2
            i = 0
            print "Processed:", len(process_readIDs), len(read_aligns), newReadsMapped, elapsed
        i+=1
    #print i, len(read_aligns)
    # Process the rest of records
    process_readIDs = totalReadIDs[0]
    newReadsMapped = processPairedReads(process_readIDs, read_aligns, segs_dict, segPairs_counts, newSegPairs_counts, segPairs_txs)
    readsMapped = [readsMapped[0]+newReadsMapped[0], readsMapped[1]+newReadsMapped[1], readsMapped[2]+newReadsMapped[2]]
    print "Processed:", len(process_readIDs), len(read_aligns), newReadsMapped
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
        
        for segPair in sorted(newSegPairs_counts.iterkeys()):
            count = newSegPairs_counts[segPair]
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


#log_f = open("pairends_log.txt", "w")
#log_f2 = open("pairends_log2.txt", "w")
#log_f3 = open("pairends_log3.txt", "w")

if __name__ == "__main__":
    pairedEndMode = False
    if sys.argv[1] == "--single":
        pairedEndMode = False
    elif sys.argv[1] == "--paired":
        pairedEndMode = True
    else:
        print("Invalid single/paired end parameter")
        usageExit()

    segmentReferenceFilename = sys.argv[2]  # path to the segments.fa reference 

    if pairedEndMode:
        kallistoOutputDir1 = sys.argv[3]  
        kallistoOutputDir2 = sys.argv[4]
        outputDir = sys.argv[5]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[6])    # How many reads alignments to store before processing
        mainPairedEnd(segmentReferenceFilename, kallistoOutputDir1, kallistoOutputDir2, outputDir, process_count)
    else:
        kallistoOutputDir1 = sys.argv[3]  
        outputDir = sys.argv[4]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[5])    # How many reads alignments to store before processing
        mainSingleEnd(segmentReferenceFilename, kallistoOutputDir1, outputDir, process_count)
    
#log_f.close()
#log_f2.close()
#log_f3.close()
