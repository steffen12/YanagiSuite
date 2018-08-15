from collections import defaultdict
from collections import Counter

from Bio import SeqIO

import sys, pickle
import time

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

def load_SegmentsLib(seg_file):
    segs_dict = {}
    segIDs = []
    for record in SeqIO.parse(seg_file, "fasta"):
        seg = Seg(record.description)
        segs_dict[seg.ID] = seg
        segIDs.append(seg.ID)
        seg.length = len(record.seq)
    return(segs_dict, segIDs)

######################################
######################################

def saveObj(obj, name):
    pickle.dump(obj, open(name+".p", "wb"))

def loadObj(name):
    return pickle.load(open(name+".p", "rb"))

def readSAM(fin):
    readAligns = defaultdict(list)
    for line in fin:
        if line[0] == "@": #If header line:
            continue
        print(line)
        tokens = line.strip().split("\t")
        readID, flags, segID, pos = tokens[:4]
        #readLen = len(tokens[9])
        readAligns[readID].append((segID, flags))
    return(readAligns)

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
        
def processPairAligns(aligns1, aligns2, segsDict):
    segPairs_counts = Counter()
    newsegPairs_counts = Counter()
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    
    segPairsTxs = {}
    for readID in aligns1:
        a1 = aligns1[readID]
        a2 = aligns2[readID]
        mapped = False
        passed_txs = False
        for align in itertools.product(a1, a2):
            segID1, segID2 = align[0][0], align[1][0]
            pairkey = segID1 + "_" + segID2
            if pairkey in segPairsTxs:
                txs = segPairsTxs[pairkey]
            else:
                seg1 = segsDict[segID1]
                seg2 = segsDict[segID2]
                txs = (seg1.txs & seg2.txs)
                segPairsTxs[pairkey] = txs
                validOri = validOrientation(align[0][1], align[1][1])
                if len(txs) < 1 and validOri: # alignment between segs with no tx in common
                    newsegPairs_counts[pairkey] += 1    # Counted as a novel junction
                elif not validOri:  # not valid alignment
                    passed_txs = True
                else:
                    segPairs_counts[pairkey] += 1
                    mapped = True
                    break
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print("Mapped Reads:", readsMapped[0], "Unmapped Reads:", readsMapped[1]+readsMapped[2], \
          "Unmapped Due Txs:", readsMapped[1], "Unmapped Due Direction:", readsMapped[2])
    return(segPairs_counts, newsegPairs_counts)

endIdx = int(sys.argv[1])-1

fin = sys.stdin
aligns = readSAM(fin)
if endIdx == 0:
    saveObj(aligns, "aligns1")
else:

    segmentReferenceFilename = sys.argv[2]
    outputCountsFilename = sys.argv[3]
    
    aligns1 = loadObj("aligns1")

    print("Loading Segments Lib...")
    start_t = time.time()
    segsDict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print("Done!")
    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)

    print("Processing Alignments...")
    start_t = time.time()
    segPairs_counts, newsegPairs_counts = processAligns(aligns1, aligns, segsDict)
    print("Done!")
    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)

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
