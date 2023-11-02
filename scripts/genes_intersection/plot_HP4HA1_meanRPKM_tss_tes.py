## click the Run button above to run following command



import matplotlib as mpl
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

## Arial font
from matplotlib.font_manager import FontProperties # 关于字体的设定使用的模块
arial_font = FontProperties(fname='/scratch/wl57926/wllab/wentao/font/arial.ttf') # 设定路径为自己安装的路径


ax1 = plt.subplot2grid((1,1), (0,0))

# —————————————————————————————————————————————————————设定灰色虚线——————————————————————————————————————————————————————

## bin26 和 bin125 上的灰色虚线的粗细
lw = 2
# 画bin26的灰色虚线（TAD的第一个bin）
plt.axvline(26 ,color = "grey", linestyle='--', linewidth=lw) # 
# 画bin125的灰色虚线 （TAD的最后一个bin）
plt.axvline(125,color = "grey",linestyle='--', linewidth=lw) # 

# —————————————————————————————————————————————————————设定TS和NTS——————————————————————————————————————————————————————
# TS 和 NTS 线条粗细
linew = 3

# TS 和 NTS 线条上的圆点的大小
msize = 2

## 读取 文件
ts = np.loadtxt("HP4HA1_No50_sorted_TS_rpkm_mean.bed") ##  TS 
nts = np.loadtxt("HP4HA1_No50_sorted_NTS_rpkm_mean.bed") ## NTS

## TS 线条颜色
ts_color="red"
## 画 TS
plt.plot(ts[:,:1],ts[:,1:2],marker='o',color=ts_color,linewidth=linew,markersize=msize, label='TS')

## NTS 线条颜色
nts_color = "blue"
## 画 NTS
plt.plot(nts[:,:1],nts[:,1:2],marker='o',color=nts_color,linewidth=linew,markersize=msize, label="NTS" )

# —————————————————————————————————————————————————————图片的整体设置——————————————————————————————————————————————————————

plt.legend(frameon=False,loc="upper right", fontsize= 18 ) #设置图例无边框，将图例放在左上角
plt.rcParams['figure.figsize']=(10,7.5) #图形大小
plt.rcParams['savefig.dpi'] = 600 #图片像素
plt.rcParams['figure.dpi'] = 600 #分辨率


plt.tick_params(labelsize=18)
ax1.set_ylabel ('AFB1-dG tXR-Seq RPKM', fontproperties=arial_font, fontsize=22)


bwith = 1.5 #边框宽度设置为2
bwith_up = 0.2
bwith_right = 0.2
TK = plt.gca()#获取边框
TK.spines['bottom'].set_linewidth(bwith) #下边框
TK.spines['left'].set_linewidth(bwith) # 左边框
TK.spines['top'].set_linewidth(bwith_up) # 上边框
TK.spines['right'].set_linewidth(bwith_right) # 右边框


# —————————————————————————————————————————————————————画x轴的标签—————————————————————————————————————————————————————

ori_xticks =[0,26,125,150]
replaced_xticks = copy.deepcopy(ori_xticks)
replaced_xticks=['-6kb','TSS','TES','+6kb']
plt.xticks(ori_xticks,replaced_xticks)
labels = replaced_xticks
ax1.set_xticklabels(labels, fontproperties=arial_font)
plt.tick_params(labelsize=18)

# —————————————————————————————————————————————————————设定图片title—————————————————————————————————————————————————————


# —————————————————————————————————————————————————————保存图片—————————————————————————————————————————————————————
output="HP4HA1_mean_RPKM.png"
plt.savefig(output, dpi=600)
#plt.show()


