import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ts = np.loadtxt("./results/rep1_beforeplot.noheader.TS.txt")

nts = np.loadtxt("./results/rep1_beforeplot.noheader.NTS.txt")
plt.plot(ts[:,:1],ts[:,1:2],color="red",linewidth=1,label='TS' )
plt.plot(nts[:,:1],nts[:,1:2],color="blue",linewidth=1,label="NTS" )
plt.legend(loc='upper right', frameon=False)
plt.show()

    # 新的github脚本算出来的，差别在于，这个脚本只拓展2kb，没有去重叠gene，分区150（up25+100+25down）
