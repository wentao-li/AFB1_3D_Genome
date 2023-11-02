from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy


ax1 = plt.subplot2grid((1,1), (0,0))

# ———————————————————————————————————————————————————————————————————————————————————————————————————————————

lw = 2
# bin26 
plt.axvline(26 ,color = "grey", linestyle='--', linewidth=lw) # 
# bin125
plt.axvline(125,color = "grey",linestyle='--', linewidth=lw) # 

# ———————————————————————————————————————————————————————————————————————————————————————————————————————————
linew = 3
msize = 2

## read
ts = np.loadtxt("HP4HA1_No50_sorted_TS_rpkm_mean.bed") ##  TS 
nts = np.loadtxt("HP4HA1_No50_sorted_NTS_rpkm_mean.bed") ## NTS

ts_color="red"
plt.plot(ts[:,:1],ts[:,1:2],marker='o',color=ts_color,linewidth=linew,markersize=msize, label='TS')

nts_color = "blue"
plt.plot(nts[:,:1],nts[:,1:2],marker='o',color=nts_color,linewidth=linew,markersize=msize, label="NTS" )

# ———————————————————————————————————————————————————————————————————————————————————————————————————————————

plt.legend(frameon=False,loc="upper right", fontsize= 18 ) 
plt.rcParams['figure.figsize']=(10,7.5) 
plt.rcParams['savefig.dpi'] = 600 
plt.rcParams['figure.dpi'] = 600 


plt.tick_params(labelsize=18)
ax1.set_ylabel ('AFB1-dG tXR-Seq RPKM', fontsize=22)


bwith = 1.5
bwith_up = 0.2
bwith_right = 0.2
TK = plt.gca()
TK.spines['bottom'].set_linewidth(bwith) 
TK.spines['left'].set_linewidth(bwith) 
TK.spines['top'].set_linewidth(bwith_up)
TK.spines['right'].set_linewidth(bwith_right)


# ——————————————————————————————————————————————————————————————————————————————————————————————————————————

ori_xticks =[0,26,125,150]
replaced_xticks = copy.deepcopy(ori_xticks)
replaced_xticks=['-6kb','TSS','TES','+6kb']
plt.xticks(ori_xticks,replaced_xticks)
labels = replaced_xticks
ax1.set_xticklabels(labels,) 
plt.tick_params(labelsize=18)

# ——————————————————————————————————————————————————————————————————————————————————————————————————————————

output="HP4HA1_mean_RPKM.png"
plt.savefig(output, dpi=600)
#plt.show()


