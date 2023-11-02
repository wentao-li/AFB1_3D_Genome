#!/bin/bash

## required files

ref_fa=refined.hg38.fa
sample_bed=HP4HA1_No50_sorted.bed
area_file=ENCFF018XKF.sorted.cut.TADnum.bed

###################################################### step1  #####################################################
## read area_file, eg. TAD.bed/Loop.bed/Gene.bed; Then pairwise comparison of two areas, if the overlap is larger of the gap, then filter both of them.

rm_output=ENCFF018XKF.sorted.cut.TADnum.filtered.bed
gap=0

python3 rm_intersection.py $area_file $rm_output $gap


##################################################### step2 #######################################################
## Calculate the areas that should be expanded upstream and downstream separately: if there is a TAD with another adjacent TAD, their distance/2 will be used as the extension length; If there are no adjacent TADs, the default extension is 1kb;
## Then seperate each area into 150bins.

distance=1000
bins_output=ENCFF018XKF.sorted.cut.TADnum.filtered.bins.bed

python3 mid.py $bins_output $rm_output $distance


###################################################### step3  #####################################################
## Extract the corresponding sequence using the bins file

fa_output=ENCFF018XKF.sorted.cut.TADnum.filtered.bins.fa

module load BEDTools/2.30.0-GCC-10.2.0
bedtools getfasta -fi $ref_fa -bed $bins_output -fo $fa_output


###################################################### step4  #####################################################
## Calulation G/g number of the corresponding sequence

Gnumber_output=ENCFF018XKF.sorted.cut.TADnum.filtered.bins.fa.Gnumber.txt

python3 gene.py -fa $fa_output  -ct -gt -gco $Gnumber_output


###################################################### step5  #####################################################
## Calculate reads that is intersected by the area file

reads_output=ENCFF018XKF_HP4HA1_intersect.txt

module load BEDTools/2.30.0-GCC-10.2.0
bedtools intersect -a $bins_output -b $sample_bed -wa -c -F 0.5 > $reads_output


###################################################### step6  #####################################################
## get total mapped reads from sample_bed file

mapped_reads=$(grep "^" $sample_bed | wc -l |  cut -d " " -f1)


###################################################### step7  #####################################################
## RPKGM

rpkgm_output=ENCFF018XKF_HP4HA1_intersect_rpkgcm.txt
python3 gene.py --RPKGCM  -gcf $all_GCout -inter $reads_output -out $rpkgm_output -map $mapped_reads


###################################################### step8  #####################################################
## get mean_RPKGM of 150bins
mean_output=ENCFF018XKF_HP4HA1_intersect_rpkgcm_mean.txt
python3 get_mean.py $$rpkgm_output $mean_output $factor

###################################################### step9  #####################################################


