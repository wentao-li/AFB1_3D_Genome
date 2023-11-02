#!/bin/python

import sys

inputFile=sys.argv[1]
out=sys.argv[2]
gap = int(sys.argv[3])

# read TAD/Loop/Gene...bed file and store the value in the dictionary
## format of bed belike: chr1       1115000 1245000 TAD1    1.2660695872919137      None

OP=open(inputFile,"r")
st=False
geneDict = {}
for lines in OP:
    lineList = lines.strip().split("\t")
    name = lineList[-3] ## ensure uniq name
    geneDict[name] = lineList

# Pairwise comparison of two genes, if the overlap is larger of the gap, then filter both of them.
comparedDict = geneDict
abandomList=[]
for key,value in geneDict.items():
    if key  in abandomList:
        continue
    chrom, start, end, name, score, strand = value
    start = int(start) - gap if int(start) - gap >=0 else 0
    end = int(end) + gap
    tag = False
    for key2,value2 in geneDict.items():
        if key2==key:
            continue
        chrom2, start2, end2, name2, score2, strand2 = value2
        start2,end2 = int(start2),int(end2)
        if chrom2==chrom and strand2==strand:
            if end2 < start or start2 > end:
                pass
            else:    
                tag = True
                abandomList.append(key2)
    if tag:
        abandomList.append(key)

# write results to a file.
OT=open(out,"w")
for k,v in geneDict.items():
    if k not in abandomList:
        v[1] = str(v[1])
        v[2] = str(v[2])
        OT.write('{v}\n'.format(v='\t'.join(v)))
OP.close()
OT.close()    
