import sys

def plot_distance_all(in_file, output):
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    ax1 = plt.subplot2grid((1,1), (0,0))

    ## read input file
    gr = np.loadtxt(in_file,delimiter="\t")
    x = gr[:,0:1]
    y = gr[:,1:2]
    linew = 2
    msize = 1.5
    plt.plot(x,y,color="#4B0082",marker='o', linewidth=linew, markersize=msize,)
    
    
    
    lw = 2
    plt.axvline(125,color = "#069AF3",linestyle='--', alpha=0.8,linewidth=lw)
    plt.axvline(26,color = "#069AF3",linestyle='--', alpha=0.8,linewidth=lw)
    #plt.legend(frameon=False,loc="upper right",fontsize='large') 
    plt.rcParams['figure.figsize']=(8,6) #figure_size
    plt.rcParams['savefig.dpi'] = 600 #figure_dpi
    plt.rcParams['figure.dpi'] = 600 #dpi
    #plt.ylabel('RPK(GC)M', fontsize=20)
    plt.ylim(0.8,1.8)
    ori_xticks =[9,76,141]
    replaced_xticks = copy.deepcopy(ori_xticks)
    replaced_xticks=['Inter-TAD','Intra-TAD','Inter-TAD',]


    plt.xticks(ori_xticks,replaced_xticks)
    plt.tick_params(labelsize=18)
    #
    plt.ylabel('AFB1-dG repair level (RPKGM)', fontsize=20)
    bwith = 0.2 #width
    bwithxy = 2 
    TK = plt.gca()
    TK.spines['bottom'].set_linewidth(bwithxy)
    TK.spines['left'].set_linewidth(bwithxy)
    TK.spines['top'].set_linewidth(bwith)
    TK.spines['right'].set_linewidth(bwith)

    # save figure
    plt.savefig(output, dpi=600)#
    

in_file=sys.argv[1]
out_file=sys.argv[2]
plot_distance_all(in_file=in_file, output=out_file)
