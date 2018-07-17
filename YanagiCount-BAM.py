import time
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
import itertools
import sys, pickle
import os

#Todo:
#1) Stdin, paired-end
#2) paired end frag
#3) Decide multimapped
#4) Valid Orientation
#from InputLoader import *

#python YanagiCount-BAM.py --single --stdin1 segments_ENSG00000100842.fa YanagiCountOutput 100

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
    for record in SeqIO.parse(seg_file, "fasta"):
        seg = Seg(record.description)
        segs_dict[seg.ID] = seg
        seg.length = len(record.seq)
    return(segs_dict)

def usageExit():
    print("Usage: ")
    print("--single or --paired for reads")
    print("--stdin1 or --stdin2 or --bam for input mode")
    print("segment reference file")
    print("bam files (if --bam mode)")
    print("output directory")
    print("process count")
    sys.exit()

def read_alignment(f, read_aligns, readIDs, end_idx):  
    line = f.readline()
    if not line:
        return ("", False) #Return when file is done
    if line[0] == "@": #If header line:
        return ("@", False)
    tokens = line.split("\t")
    readID = tokens[0][:-2]
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
    

def validLength(fl_params, length):
    f_mean, f_sd, sd_factor = fl_params
    return True #((length < (f_mean + sd_factor * f_sd)) and (length > (f_mean - sd_factor * f_sd)))

def validOrientation(flags1, flags2):
    return True#flags1+flags2 == 16

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
    print str(float(c)/len(readIDs)) + " alignments per read"
    return readsMapped

def mainSingleEnd(stdinTrue, stdinPair, pairedEndMode):
    segmentReferenceFilename = sys.argv[3]  # path to the segments.fa reference 
    if not stdinTrue:
        BAM1 = sys.argv[4]  
        outputDir = sys.argv[5]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[6])    # How many reads alignments to store before processing 
    else:
        outputDir = sys.argv[4]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[5])    # How many reads alignments to store before processing 
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputCountsFilename = outputDir + "/seg_counts.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    #DExons = load_disjointExons(seg_inDir)
    #txs2exs = load_Txs2Exs(seg_inDir, DExons)
    segs_dict = load_SegmentsLib(segmentReferenceFilename)
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
    #TODO clean up based on sorted readIDs
    readsMapped = [0,0,0]
    if(stdinTrue):
        f1 = sys.stdin
    else:
        f1 = open(BAM1) 
    fout = open(outputCountsFilename, "w")
    done = False
    readIDs = [[], []]
    i = 0
    while not done:
        if not done: # not done with the first bam
            readID, aligned = read_alignment(f1, read_aligns, readIDs, 0)
            done = not readID

        # Process some records (if stored more than process_count reads in memory)
        if not stdinTrue and len(readIDs[0]) > process_count:
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
def processPairedReads(readIDs, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs):
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    c=0
    for readID in readIDs:
        aligns = read_aligns.pop(readID)
        mapped = False
        print(aligns)
        passed_txs = False
        c += len(aligns.aligns[0]+aligns.aligns[1])
        for align in itertools.product(aligns.aligns[0], aligns.aligns[1]):
            segID1, segID2 = align[0][0], align[1][0]
            pairkey = segID1 + "_" + segID2
            print(segID1, segID2)
            if pairkey in segPairs_txs:
                txs = segPairs_txs[pairkey]
            else:
                seg1 = segs_dict[segID1]
                seg2 = segs_dict[segID2]
                txs = (seg1.txs & seg2.txs)
                segPairs_txs[pairkey] = txs
            if len(txs) < 1:    # alignment between segs with no tx in common
                if validOrientation(align[0][1], align[1][1]):  # Valid alignment
                    newsegPairs_counts[pairkey] += 1    # Counted as a novel junction
                continue
            print(align)
            if not validOrientation(align[0][1], align[1][1]):  # not valid alignment
                passed_txs = True
                #print align[0][1], align[1][1]
                continue
            segPairs_counts[pairkey] += 1
            mapped = True
            break
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print (c*1.0)/len(readIDs)
    return readsMapped

#@profile
def mainPairedEnd(stdinTrue, stdinPair):
    segmentReferenceFilename = sys.argv[3]  # path to the segments.fa reference 
    if not stdinTrue:
        BAM1 = sys.argv[4]  
        BAM2 = sys.argv[5]
        outputDir = sys.argv[6]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[7])    # How many reads alignments to store before processing 
    else:
        outputDir = sys.argv[4]    # output Dir to dump the segcounts file
        process_count = int(sys.argv[5])    # How many reads alignments to store before processing 
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputCountsFilename = outputDir + "/seg_counts.tsv"

    print "Loading Segments Lib...",
    start_t = time.time()
    #DExons = load_disjointExons(seg_inDir)
    #txs2exs = load_Txs2Exs(seg_inDir, DExons)
    segs_dict = load_SegmentsLib(segmentReferenceFilename)
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
    newsegPairs_counts = Counter()
    segPairs_txs = {}
    #TODO clean up based on sorted readIDs
    readsMapped = [0,0,0]
    if(stdinTrue):
        read_aligns = loadObj("read_aligns")
        print(read_aligns)
        f2 = sys.stdin
        done = [True, False]
    else:
        f1 = open(BAM1)
        f2 = open(BAM2)
        done = [False, False]
    fout = open(outputCountsFilename, "w")
    readIDs = [[], []]
    i = 0
    while not all(done):
        if not done[0]: # not done with the first bam
            readID, aligned = read_alignment(f1, read_aligns, readIDs, 0)
            done[0] = not readID
        if not done[1]: # not done with the second bam
            readID, aligned = read_alignment(f2, read_aligns, readIDs, 1)
            done[1] = not readID

        # Process some records (if stored more than process_count reads in memory)
        if not stdinTrue:
            if len(readIDs[0]) > process_count and len(readIDs[1]) > process_count:
                process_readIDs = readIDs[0][:process_count]
                readIDs[0] = readIDs[0][process_count:]
                readIDs[1] = readIDs[1][process_count:]
                print i, len(read_aligns)
                start_t_2 = time.time()
                rMapped = processPairedReads(process_readIDs, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
                elapsed = time.time() - start_t_2
                print "Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed
        else:
            if len(readIDs[1]) > process_count:
                process_readIDs = readIDs[1][:process_count]
                readIDs[1] = readIDs[1][process_count:]
                print i, len(read_aligns)
                start_t_2 = time.time()
                rMapped = processPairedReads(process_readIDs, read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
                elapsed = time.time() - start_t_2
                print "Stdin Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed    
        i+=1
    print i, len(read_aligns)
    # Process the rest of records
    if stdinTrue:
        rMapped = processPairedReads(readIDs[1], read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
    else:
        rMapped = processPairedReads(readIDs[0], read_aligns, segs_dict, segPairs_counts, newsegPairs_counts, segPairs_txs)
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
            print(segPair)
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
    stdinPair = 0
    if sys.argv[2] == "--stdin1":
        stdinTrue = True
        stdinPair = 1
    elif sys.argv[2] == "--stdin2":
        stdinTrue = True
        stdinPair = 2
    elif sys.argv[2] == "--bam":
        stdinTrue = False
        stdinPair = 0
    else:
        print("Invalid input mode parameter")
        usageExit()

    if pairedEndMode and (not stdinTrue or (stdinTrue and stdinPair == 2)):
        mainPairedEnd(stdinTrue, stdinPair)
    else:
        mainSingleEnd(stdinTrue, stdinPair, pairedEndMode)
    
