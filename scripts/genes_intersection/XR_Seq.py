import argparse
from os import listdir
from statistics import mean
import numpy as np
import pandas as pd


class BedLine:
    def __init__(self, line):
        columns = line.strip().split()
        self.chrom = columns[0]
        self.start = int(columns[1])
        self.end = int(columns[2])
        self.name = columns[3]
        try:
            self.score = int(columns[4])
        except ValueError:
            self.score = columns[4]
        self.strand = columns[5]
        try:
            self.count = float(columns[-1])
        except ValueError:
            pass

    def bed_line(self):
        """
        This part won't work properly with versions of Python earlier than 3.6 because by default and class attributes
        aren't ordered until Python 3.6.
        Be careful, Longleaf loads Python 3.5.1 by default, but has 3.6.8 if you don't load anything.
        """
        return "\t".join([str(value) for value in vars(self).values()]) + "\n"

    def sequence(self, assembly):

        query = {
            'genome': assembly,
            'chrom': self.chrom,  # might need to change chrom names
            'start': self.start,
            'end': self.end
        }

        # UCSC gives you the sequence on the + strand, so you need to keep it in mind for genes on - strand.
        r = requests.get('https://api.genome.ucsc.edu/getData/sequence?', params=query)
        seq = r.json()['dna'].upper()

        if self.strand == '+':
            return seq

        return seq.translate(string.maketrans("ACTG", "TGAC"))

    def __len__(self):
        return self.end - self.start


def rpkm(input_file,read_count_file,output_file):
#    read_count_file = name + "_total_mappedreads.txt"
    with open(read_count_file) as f:
        total_mapped_count = int(f.readline().strip())
        nf = 10 ** 9 / total_mapped_count

    intersect_list = [input_file]

    for intersect in intersect_list:
        print("Normalizing " + intersect)
        with open(intersect) as f:
            for line in f:
                line = BedLine(line)
                line.count = (line.count * nf) / len(line)
                with open(output_file, 'a') as d:
                    d.write(line.bed_line())

def get_rpkm_mean(bed_file, output, factor="rpkm"):
    rp = pd.read_csv(bed_file, sep="\t", header=None)
    rp.columns = ['chr', 'start', 'end', 'gene', 'bin', "strand", factor]
    rp.groupby(by='bin')[factor].agg([np.mean]).to_csv(output, sep="\t", header=None)


parser = argparse.ArgumentParser(description='Performs the Python half of the XR-Seq analysis.')
parser.add_argument('-i', '--inputfile', help='Name of the sample input file.')
parser.add_argument('-o', '--outputfile', help='Name of the sample output file.')
parser.add_argument('-t', '--totalread', help='Name of the total read count file.')
parser.add_argument('--rpkm', action='store_true', help='Perform RPKM normalization.')
parser.add_argument('--mean', action='store_true', help='get mean of RPKM file.')
args = parser.parse_args()

if args.rpkm:
    rpkm(args.inputfile,args.totalread,args.outputfile)
if args.mean:
    get_rpkm_mean(args.inputfile,args.outputfile)
