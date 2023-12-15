import re
import regex
count = 0
seqCount = 0
fileName = input("Enter the fastq file name: ")
print("This script extracts sequences close to the specified length that lie between the Gibson Assembly sequences of pUC-GAPDH.")
print("Do NOT use it on any other plasmid backbone!")
insSize = int(input("Enter the size of the insert: "))
phred_cutoff = int(input("Enter the phred cutoff: "))
#insSize = int(input("Enter the size of the insert: "))
data = open(fileName, 'r')
outfile = open("all_inserts_q" +str(phred_cutoff) +".fasta", 'w')

def getInsert(seqn, gbSeq1, gbSeq2, seqCount, insSize, hasInsert):
    # adds insert read to outfile if it is an acceptable size (+/- 20bp)
    # generates the insert by splitting the sequence and checking if it is an acceptable size (+/- 20bp)
    first = seqn.split(gbSeq1)
    startIns = len(first[0]) + len(gbSeq1)
    second = first[1].split(gbSeq2)
    finishIns = startIns + len(second[0])
    insert = ""
    if (len(second[0]) > insSize - 21) and (len(second[0]) < insSize + 21):
        insert = second[0]
        hasInsert = True
    else: pass
    return hasInsert, startIns, finishIns, insert

hasInsert = False
startIns = 0
finishIns = 0

for line in data:
    # if the read length is at least the size of the backbone
    if count % 4 == 1 and len(line.strip()) > 2900:
        seqn = line.strip()
        # check for the Gibson Assembly sequences on the forward strand, allow for one mismatch
        if regex.search("CCACATCGCTCAGACAC{e<=1}", seqn) and regex.search("ACTGGCCGTCGTTTTAC{e<=1}", seqn):
            hasInsert, startIns, finishIns, insert = getInsert(seqn, regex.findall("CCACATCGCTCAGACAC{e<=1}", seqn)[0], regex.findall("ACTGGCCGTCGTTTTAC{e<=1}", seqn)[0], seqCount, insSize, hasInsert)
        # check for the Gibson Assembly sequences on the reverse strand
        elif regex.search("GTAAAACGACGGCCAGT{e<=1}", seqn) and regex.search("GTGTCTGAGCGATGTGG{e<=1}", seqn):
            hasInsert, startIns, finishIns, insert = getInsert(seqn, regex.findall("GTAAAACGACGGCCAGT{e<=1}", seqn)[0], regex.findall("GTGTCTGAGCGATGTGG{e<=1}", seqn)[0], seqCount, insSize, hasInsert)
        else: pass
    # get the Phred scores
    elif count % 4 == 3 and hasInsert:
        phred = line[startIns:finishIns]
        phredScore = 0
        for i in range(len(phred)):
            phredScore += (ord(phred[i]) - 33)
            phredAv = phredScore/len(phred)
        if phredAv >= phred_cutoff:
            seqCount +=1
            print(">" + str(seqCount) + "_q" + str(phredAv), file=outfile)
            print(insert, file=outfile)
            #print("+", file=outfile)
            #print(str(phred), file=outfile)
        else: pass
        hasInsert = False
    else: pass
    count +=1
data.close()
outfile.close()

# get the first 50 reads
data = open("all_inserts_q" +str(phred_cutoff) +".fasta", 'r')
outfile = open("first_50_inserts.fasta", 'w')
count = 0
for line in data:
    if count < 100:
        print(line.strip(), file=outfile)
    else: pass
    count +=1
data.close()
outfile.close()

print("Number of sequences extracted: ", seqCount)