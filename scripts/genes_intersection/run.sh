#!/bin/bash

module load BEDTools/2.30.0-GCC-11.3.0

# set sample bed file and output sample ID
sample_id=HP4HA1_No50_sorted
sample_bed=HP4HA1_No50_sorted.bed

# required files: gene bed file, reference genome,
prefix="hg38_genes_and_expression_uniqname_sorted"
gene_list ="hg38_genes_and_expression_uniqname_sorted.bed"
fasta=/work/wllab/public/00.ref/hg38/hg38_male/refined.hg38.fa


## criteria1: length >3kb
input_file = gene_list 
len=3000
out1=$prefix\_len$len\.bed
echo | awk -v gene_len="$len" '(($3-$2)>gene_len){print $0}' $gene_list > $out1

out2=$prefix\_len$len\_newformat.bed
awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$4}' $out1 > $out2


## criteria2: expand up_downstream to 6kb and remove  genes with intersection, then output genes with non-expanded region(because next step will expand it again, so need to merge thess two steps.)
gap=6000
out3=$prefix\_len$len\_newformat_rm$gap\.bed
python rm_intersection.py $out2 $out3 $gap

## criteria3: expand region to 6kb and split into 150 bins.
extend_len=6000
out4=$prefix\_len$len\_newformat_rm$gap\_extend$extend_len\_150bins\.bed
python divide_transcripts.py $out3 $out4 $extend_len


## get_count from sample_bed
total_read=$sample_id\_total_mappedreads.txt
grep -c "^"  $sample_bed > $total_read

# get TS/NTS
ts_bed=$sample_id\_TS.bed
nts_bed=$sample_id\_NTS.bed
bedtools intersect -c -a $out4 -b $sample_bed -wa -c -S -F 0.5 > $ts_bed
bedtools intersect -c -a $out4 -b $sample_bed -wa -c -s -F 0.5 > $nts_bed



input_file=$out4
minus_out=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_minus.bed
plus_out=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_plus.bed
awk '($6=="-"){print $0}' $input_file > $minus_out
awk '($6=="+"){print $0}' $input_file > $plus_out
##


#### getfa
output1=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_getfasta.fa
plus_out_fa=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_plus_getfasta.fa
minus_out_fa=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_minus_getfasta.fa
bedtools getfasta -fi $fasta -bed $plus_out -fo $plus_out_fa
bedtools getfasta -fi $fasta -bed $minus_out -fo $minus_out_fa

##
#### step2: cal gc for TS and NTS
TS_GCout=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_getfasta_TS_Gnumber.txt
NTS_GCout=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_getfasta_NTS_Gnumber.txt
all_GCout=hg38_genes_and_expression_uniqname_sorted_len3000_newformat_rm6000_extend6000_150bins_getfasta_all_Gnumber.txt
##
#####TS
python3 gene.py -fa $plus_out_fa  -ct -gco plusC_forTS.txt
python3 gene.py -fa $minus_out_fa -gt -gco minusG_forTS.txt
cat plusC_forTS.txt minusG_forTS.txt > $TS_GCout && rm plusC_forTS.txt minusG_forTS.txt
##
####### NTS
python3 gene.py -fa $plus_out_fa -gt  -gco plusG_forNTS.txt
python3 gene.py -fa $minus_out_fa  -ct -gco minusC_forNTS.txt
cat plusG_forNTS.txt minusC_forNTS.txt >$NTS_GCout && rm plusG_forNTS.txt minusC_forNTS.txt


## get_readcount
TS_count=$sample_id"_TS_intersect.txt"
NTS_count=$sample_id"_NTS_intersect.txt"
all_count=$sample_id"_all_intersect.txt"
#TS
bedtools intersect -a $input_file -b $sample_bed -wa -c -S -F 0.5 > $TS_count
#NTS
bedtools intersect -a $input_file -b $sample_bed -wa -c -s -F 0.5 > $NTS_count
#all
bedtools intersect -a $input_file -b $sample_bed -wa -c -F 0.5 > $all_count

## get_mapped_reads
mapped_reads=$(grep "^" $sample_bed | wc -l |  cut -d " " -f1)

## get RPKGM
#TS
TS_RPKGCM=$sample_id"_TS_RPKGCM.txt"
NTS_RPKGCM=$sample_id"_NTS_RPKGCM.txt"
all_RPKGCM=$sample_id"_all_RPKGCM.txt"
python3 gene.py --RPKGCM  -gcf $TS_GCout -inter $TS_count -out $TS_RPKGCM -map $mapped_reads
python3 gene.py --RPKGCM  -gcf $NTS_GCout -inter $NTS_count -out $NTS_RPKGCM -map $mapped_reads
python3 gene.py --RPKGCM  -gcf $all_GCout -inter $all_count -out $all_RPKGCM -map $mapped_reads


# get mean of RPKM
TS_meanRPKGCM=$sample_id\_TS_rpkm_mean.bed
NTS_meanRPKGCM=$sample_id\_NTS_rpkm_mean.bed
ml SciPy-bundle/2022.05-foss-2022a
python XR_Seq.py -i $TS_RPKGCM --mean -o $TS_meanRPKGCM
python XR_Seq.py -i $NTS_RPKGCM --mean -o $NTS_meanRPKGCM


