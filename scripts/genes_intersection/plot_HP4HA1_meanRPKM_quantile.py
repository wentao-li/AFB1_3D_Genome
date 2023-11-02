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

