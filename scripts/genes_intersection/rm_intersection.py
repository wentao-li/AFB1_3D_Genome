#!/bin/python
import sys
#inputFile="test.bed"
inputFile=sys.argv[1]
#inputFile="hg38_genes_and_expression_300exp_15kb.bed"
OP=open(inputFile,"r")

st=False
geneDict = {}
for lines in OP:
    lineList = lines.strip().split("\t")
    ## warning
    ##    if length(lineList)
    gene = lineList[-1]
    geneDict[gene] = lineList
    #print(geneDict)

gap = int(sys.argv[3])
comparedDict = geneDict
abandomList=[]
for key,value in geneDict.items():
    if key  in abandomList:
        continue
    print("1st")
    print(value)
    #chrom, start, end, gene, length, strand=value
    chrom, start, end, strand, length, gene=value
    start = int(start) - gap if int(start) - gap >=0 else 0
    end = int(end) + gap
    tag = False
    for key2,value2 in geneDict.items():
        if key2==key:
            continue
        #chrom2, start2, end2, gene2, length2, strand2=value2
        chrom2, start2, end2, strand2, length2, gene2=value2
        start2,end2 = int(start2),int(end2)
        if chrom2==chrom and strand2==strand:
            if end2<start or start2>end:
                pass
            else:    
                tag = True
                print("2st")
                print(value2)
                abandomList.append(key2)
                #abandomList.append(key)
                #break
    if tag:
        abandomList.append(key)
out=sys.argv[2]
OT=open(out,"w")
#OT=open("testout.bed2","w")
for k,v in geneDict.items():
    if k not in abandomList:
        start = v[1]
        #v[1] = int(start) - gap if int(start) - gap >=0 else 0
        v[1] = str(v[1])
        end = v[2]
        #v[2] = int(end) + gap
        v[2] = str(v[2])
        OT.write('{v}\n'.format(v='\t'.join(v)))
OP.close()
OT.close()    
