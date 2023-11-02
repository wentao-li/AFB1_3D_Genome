#!/bin/python
import sys

def get_rpkm_mean(output, bed_file, factor="rpkgm"):
    '''
    :param bed_file: "chr1	35313	35553	WASH7P	1	-	0.0"
    :return:
    '''

    import numpy as np
    import pandas as pd
    rp = pd.read_csv(bed_file, sep="\t", header=None)
    rp.columns = ['chr', 'start', 'end', 'name', 'bin', "strand", factor]
    rp.groupby(by='bin')[factor].agg([np.mean]).to_csv(output, sep="\t", header=None)
    return output

in_file=sys.argv[1]
out_file=sys.argv[2]
factor = sys.argv[3]
get_rpkm_mean(output = out_file, bed_file=in_file, factor=factor)
