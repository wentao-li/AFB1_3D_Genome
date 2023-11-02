def plot_mean_RPKM_quantile_1fig(q1_file,q2_file,q3_file,q4_file,nts_q1,nts_q2,nts_q3,nts_q4,out_file):

    import matplotlib as mpl
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import copy

    ax1 = plt.subplot2grid((1,1), (0,0))

    # —————————————————————————————————————————————————————设定灰色虚线——————————————————————————————————————————————————————

    ## bin26 and bin125 with grey dash line
    lw = 2
    # bin26（TAD 1st bin）
    plt.axvline(26 ,color = "grey", linestyle='--', linewidth=lw) # 
    # bin125 （TAD last bin）
    plt.axvline(125,color = "grey",linestyle='--', linewidth=lw) # 

    # —————————————————————————————————————————————————————TS and NTS——————————————————————————————————————————————————————
    # TS NTS
    linew = 3

    # TS  NTS
    msize = 1
# —————————————————————————————————————————————————————TS——————————————————————————————————————————————————————
    
    ## read files

    q1 = np.loadtxt(q1_file) ##  TS 
    q2 = np.loadtxt(q2_file) ## NTS
    q3 = np.loadtxt(q3_file) ## NTS
    q4 = np.loadtxt(q4_file) ## NTS
    


    q1_color="red"
    plt.plot(q1[:,:1],q1[:,1:2],marker='o',color=q1_color, alpha=0.1, linewidth=linew,markersize=msize, label='Q1')


    q2_color = "red"
    plt.plot(q2[:,:1],q2[:,1:2],marker='o',color=q2_color,alpha=0.4,linewidth=linew,markersize=msize, label="Q2" )
    
    q3_color="red"
    plt.plot(q3[:,:1],q3[:,1:2],marker='o',color=q3_color,alpha=0.6,linewidth=linew,markersize=msize, label='Q3')
    
    q4_color="red"
    plt.plot(q4[:,:1],q4[:,1:2],marker='o',color=q4_color,alpha=1,linewidth=linew,markersize=msize, label='Q4')
# —————————————————————————————————————————————————————NTS——————————————————————————————————————————————————————
    nq1 = np.loadtxt(nts_q1) ##  TS 
    nq2 = np.loadtxt(nts_q2) ## NTS
    nq3 = np.loadtxt(nts_q3) ## NTS
    nq4 = np.loadtxt(nts_q4) ## NTS
    
    marker = 'None'
    lstyle='dotted'
    
    nq1_color="blue"
    plt.plot(nq1[:,:1],nq1[:,1:2],linestyle=lstyle,marker=marker ,color=nq1_color,alpha=0.1,linewidth=linew,markersize=msize, label='Q1')


    nq2_color = "blue"
    plt.plot(nq2[:,:1],nq2[:,1:2],linestyle=lstyle,marker=marker ,color=nq2_color,alpha=0.4,linewidth=linew,markersize=msize, label="Q2" )
    
    nq3_color="blue"
    plt.plot(nq3[:,:1],nq3[:,1:2],linestyle=lstyle,marker=marker ,color=nq3_color,alpha=0.6,linewidth=linew,markersize=msize, label='Q3')
    
    nq4_color="blue"
    plt.plot(nq4[:,:1],nq4[:,1:2],linestyle=lstyle,marker=marker ,color=nq4_color,alpha=1,linewidth=linew,markersize=msize, label='Q4')

    
    
    # —————————————————————————————————————————————————————setting——————————————————————————————————————————————————————

    plt.legend(frameon=False,loc="upper right", fontsize= 18 ) 
    plt.rcParams['figure.figsize']=(10,7.5)
    plt.rcParams['savefig.dpi'] = 600 
    plt.rcParams['figure.dpi'] = 600 


    plt.tick_params(labelsize=18)
    ax1.set_ylabel ('AFB1-dG tXR-Seq RPKM',  fontsize=22)


    bwith = 1.5 
    bwith_up = 0.2
    bwith_right = 0.2
    TK = plt.gca()
    TK.spines['bottom'].set_linewidth(bwith) 
    TK.spines['left'].set_linewidth(bwith) 
    TK.spines['top'].set_linewidth(bwith_up) 
    TK.spines['right'].set_linewidth(bwith_right) 


    # —————————————————————————————————————————————————————x axis label—————————————————————————————————————————————————————

    ori_xticks =[0,26,125,150]
    replaced_xticks = copy.deepcopy(ori_xticks)
    replaced_xticks=['-6kb','TSS','TES','+6kb']
    plt.xticks(ori_xticks,replaced_xticks)
    labels = replaced_xticks
    ax1.set_xticklabels(labels)
    plt.tick_params(labelsize=18)

    # —————————————————————————————————————————————————————set_title—————————————————————————————————————————————————————


    # —————————————————————————————————————————————————————save fig—————————————————————————————————————————————————————
    output = out_file
    plt.savefig(output, dpi=600)
    plt.show()


TS_dir = "./quantile/TS/"
q1= TS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile1_bin150_TS_mean_RPKM.txt"
q2= TS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile2_bin150_TS_mean_RPKM.txt"
q3= TS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile3_bin150_TS_mean_RPKM.txt"
q4= TS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile4_bin150_TS_mean_RPKM.txt"

NTS_dir = "./quantile/NTS/"
nq1= NTS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile1_bin150_NTS_mean_RPKM.txt"
nq2= NTS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile2_bin150_NTS_mean_RPKM.txt"
nq3= NTS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile3_bin150_NTS_mean_RPKM.txt"
nq4= NTS_dir + "hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_quantile4_bin150_NTS_mean_RPKM.txt"

out = NTS_dir + "/../TSandNTS_mean_RPKM.png"


plot_mean_RPKM_quantile_1fig(q1,q2,q3,q4,nq1,nq2,nq3,nq4,out)

