#! /bin/bash

sample_list=()
basename_list=()

while IFS= read -r sample; do
    sample_list+=("$sample")
done < <(tail -n +2 H4K20me1_CnR_sample_sheet.txt | cut -f3)
echo ${sample_list[@]}

for file in ${sample_list}
do

tmp=${file##*/}
basename=${tmp%%"_dm6_trim_q5_dupsRemoved_allFrags.bed"}
basename_list+=(${basename})
echo ${basename_list[@]}

echo "Starting ${basename}..."
SRPMC=$(grep "${basename}" H4K20me1_CnR_sample_sheet.txt | cut -f9)

echo ${SRPMC}

module load bedtools

bedtools genomecov -i ${basename}_dm6_trim_q5_dupsRemoved_allFrags.bed -bga -g dm6.chrom.sizes -scale ${SRPMC} > ${basename}_dm6_trim_q5_dupsRemoved_allFrags_SRPMC.bg

module load ucsctools

wigToBigWig ${basename}_dm6_trim_q5_dupsRemoved_allFrags_SRPMC.bg dm6.chrom.sizes ${basename}_dm6_trim_q5_dupsRemoved_allFrags_SRPMC.bw

done