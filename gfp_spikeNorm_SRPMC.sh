#!/bin/bash

sample_list=()
basename_list=()

# Read sample names from the file and populate sample_list array
while IFS= read -r sample; do
    sample_list+=("$sample")
done < <(tail -n +2 GFP-l3mbt_CnR_sample_sheet.txt | cut -f3)

# Debugging: Print the entire sample_list
echo "Sample list: ${sample_list[@]}"

# Loop over each sample in the sample_list
for file in "${sample_list[@]}"
do
    # Extract the basename from the filename
    tmp=${file##*/}
    basename=${tmp%%"_dm6_trim_q30_dupsKept_allFrags.bed"}
    basename_list+=("$basename")
    
    # Debugging: Print the basename list
    echo "Basename list: ${basename_list[@]}"
    
    echo "Starting ${basename}..."
    SRPMC=$(grep "$basename" GFP-l3mbt_CnR_sample_sheet.txt | cut -f9)
    
    echo "SRPMC value: ${SRPMC}"
    
    module load bedtools

    bedtools genomecov -i "${basename}_dm6_trim_q30_dupsKept_allFrags.bed" -bga -g dm6.chrom.sizes -scale "$SRPMC" > "${basename}_dm6_trim_q30_dupsKept_allFrags_SRPMC.bg"

    module load ucsctools

    wigToBigWig "${basename}_dm6_trim_q30_dupsKept_allFrags_SRPMC.bg" dm6.chrom.sizes "${basename}_dm6_trim_q30_dupsKept_allFrags_SRPMC.bw"

done
