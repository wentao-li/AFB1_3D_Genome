#!/bin/python
from collections import OrderedDict
import argparse
parser = argparse.ArgumentParser(description='XRseq')
parser.add_argument('-sep', '--seperate_genelist',action='store_true', help='True or False')
parser.add_argument('-related', '--reads_genes_relationship', type=bool,default=False, help='True or False')
parser.add_argument('-s', '--sample_bed', help='sample bed file')
parser.add_argument('-g', '--gene_list', help='gene list file')
parser.add_argument('-gt', '--Gtype',action='store_true',  help='cal G number')
parser.add_argument('-ct', '--Ctype', action='store_true',help='cal C number')
parser.add_argument('-gco', '--GCout', help='GC output')
parser.add_argument('-fa', '--fasta', help='fasta')
parser.add_argument('-rpkgcm', '--RPKGCM', action='store_true',help='RPKGCM or RPKGM')
parser.add_argument('-gcf', '--GC_file', help='file contains GC number of a seq')
parser.add_argument('-inter', '--intersect_file', help='file contains read of genes')
parser.add_argument('-out', '--output', help='output')
parser.add_argument('-id', '--sampleid', help='sampleid')
parser.add_argument('-map', '--mapped_count', type=int, help='mapped read counts in total')
parser.add_argument('-bin_rpkgcm', '--bin_RPKGCM', action='store_true',help='')
parser.add_argument('-rpf', '--RP_file', help='file contains RPKGCM of a gene')
parser.add_argument('-tadf', '--TAD_file', help='tad file')
    

args = parser.parse_args()

def seperate_genelist(input_file):
    """
    if bin contain more than one gene, output to $mul_file
    else output to $one_file
    eg. input_file="hg38_genes_and_expression_uniqname_inTAD.bed"
        chr1    1114000 1114040 C1orf159_sub1   426 -   chr1_TAD1   1

    """ 
    gene_dict = OrderedDict()
    f = open(input_file, "r")
    mul_file = open(f"{input_file}_multigenes.bed", "w")
    one_file = open(f"{input_file}_onegene.bed", "w")
    for line in f.readlines() :
        line = line.strip()
        gr = line.split("\t")
        tag = "_".join(gr[6:8])
        if tag not in gene_dict.keys():
            gene_dict[tag] = [line]
        else:
            gene_dict[tag].append(line)

    f.close() 
    for k,v in gene_dict.items():
        if len(v) == 1:
            one_file.write(v[0]+"\n")
        else:
            for line in v:
                mul_file.write(line+"\n")
    mul_file.close()
    one_file.close()


def get_GC_number(fasta,G_type=True,C_type=True, output=None):
    print(str(G_type)+str(C_type)+"\n")
    f = open(fasta, "r")
    if output is None:
        output = f"{fasta}.gc_number.txt"
    out = open(output, "w")
    for line in f.readlines():
        if line.startswith(">"):
            line = line.strip().replace(">", "")
            out.write(line+"\t")
        else:
            g_count = 0
            c_count = 0
            if G_type:
                g_count = line.count("g") +  line.count("G") 
            if C_type:
                c_count = line.count("c") +  line.count("C") 
            gc_count = g_count + c_count
#                gc_count = line.count("g") + line.count("c") + line.count("G") + line.count("C")
            out.write(str(gc_count) + "\n")

    out.close()
    return output

def get_overlap_genes(sample_list, ) :
    #chr1    1114000 1114040 C1orf159_sub1   426     -       chr1_TAD1       1
    overlap_dict=OrderedDict()
    f = open(sample_list, "w" )
    for line in f.readlines():
        line = line.strip()
        gr = line.split("\t")
        tag = "_".join()

def reads_genes_relationship(sample_bed, gene_list):
    """
    eg. gene_list = "hg38_genes_and_expression_uniqname_inTAD.bed"
    """
    pass

def get_rpkgcm(gc_file, mapped_count, intersect_file, sampleid, output=None):
    nf = 10**9/float(mapped_count)
    # read count file
    if output is None:
        output = sampleid + "_rpkgcm.txt"
    gc_dict = {}
    f = open(gc_file, "r")
    for i in f.readlines():
        gr = i.strip().split("\t")
        gc_dict[gr[0]] = int(gr[-1])

    out = open(output, "w")
    with open(intersect_file) as f:
        for line in f:
            bed_line = line.strip().split("\t")
            region = bed_line[0] + ":" +  bed_line[1] + "-" +  bed_line[2] 
            if region in gc_dict.keys():
                gc_count = int(gc_dict[region])
                if gc_count == 0:
                    bed_line[-1] = str(0)
                else:
                    bed_line[-1] = str(float(bed_line[-1]) *nf/gc_count)
                out.write("\t".join(bed_line)+"\n")
    out.close()
    return output

def gene2bin_rpkgcm(intersect_file,rpkgcm_file, tad_file,output):
    outbin_dict = OrderedDict()
    bin_dict = OrderedDict()
    gene_num = OrderedDict()    
    ta = open(tad_file,"r")
    for line in ta.readlines():
        gr = line.strip().split("\t")
        tad = gr[-1]
        bi = gr[-2]
        if tad not in bin_dict.keys():
            bin_dict[tad] = {bi:0}
            gene_num[tad] = {bi:0}
            outbin_dict[tad] = {bi:gr}
        else:
            bin_dict[tad][bi] = 0
            gene_num[tad][bi] = 0
            outbin_dict[tad][bi] = gr
    ta.close()

    rp_dict=OrderedDict()
    rp = open(rpkgcm_file,"r")
    for line in rp.readlines():
        gr = line.strip().split()
        gene = gr[3]
        score = gr[-1]
        rp_dict[gene] = float(score)
    rp.close()

    f = open(intersect_file,"r")
    for line in f.readlines():
        gr = line.strip().split()
        gene = gr[3]
        tad = gr[-2]
        bi = gr[-3]
        if tad in bin_dict.keys() and bi in bin_dict[tad].keys():    
            bin_dict[tad][bi] = bin_dict[tad][bi] + rp_dict[gene]
            gene_num[tad][bi] = gene_num[tad][bi] + 1
    f.close()
    
    out = open(output, "w")
    for tad,v in bin_dict.items():
        for bi,score in v.items() :
            bin_rp = 0
            if gene_num[tad][bi] > 0:
                bin_rp = score/gene_num[tad][bi]
            out.write("\t".join(outbin_dict[tad][bi])+"\t"+str(bin_rp)+"\n")
    out.close()

if __name__=='__main__':
    if args.seperate_genelist and args.gene_list is not None:
        seperate_genelist(args.gene_list)

    if args.reads_genes_relationship and args.sample_bed is not None and args.gene_list is not None:
        reads_genes_relationship(sample_bed=args.sample_bed, gene_list = args.gene_list)

    if args.Gtype is True or args.Ctype is True:
        print(args.Gtype)
        get_GC_number(fasta=args.fasta,G_type=args.Gtype,C_type=args.Ctype, output=args.GCout)

    if args.RPKGCM is True :
        print(args.RPKGCM)
        get_rpkgcm(gc_file= args.GC_file, mapped_count=args.mapped_count, intersect_file=args.intersect_file, sampleid=args.sampleid, output=args.output)

    if args.bin_RPKGCM is True :
        print(args.bin_RPKGCM)
        gene2bin_rpkgcm(intersect_file = args.intersect_file ,rpkgcm_file=args.RP_file, tad_file = args.TAD_file , output = args.output)


