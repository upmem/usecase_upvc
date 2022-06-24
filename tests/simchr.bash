#!/bin/bash

# USAGE: simchr.bash chromNumber [rngSeed]
# If rngSeed is not specified, a random one will be chosen instead

# EDIT THESE TWO VARIABLES TO POINT TO THE PROPER BINARIES FOR VCF2DIPLOID AND ART_ILLUMINA
# You can obtain VCF2Diploid from github.com/moselhy/vcf2diploid, you need to run "make" after you clone it to create the jar file
# You can obtain ART from https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm, you just need to untar the binary package
VCF2DIPLOID="/home/ubuntu/work/genomee/vcf2diploid/vcf2diploid.jar"
ART_ILLUMINA="/home/ubuntu/work/genomee/art_bin_MountRainier/art_illumina"

# add in the $HOME/.bash_profile
# PATH="/Users/lavenier/Documents/Projets/UPMEM/upvc1.4/VCScripts:${PATH}"

wait_jobs() {
    for I in $(jobs -p)
    do
	echo "waiting $I"
	wait $I
    done
}

################## START OF SIMULATION ##################

# Get the user-defined seed from input
if [ $1 ]; then
	seed=$1
else
	seed=64189 # 0xfabd
fi

# If the common variants do not already exist on disk, download them
if [ ! -f "common_all_20170710.vcf.gz" ] && [ ! -f "common_all_20170710.vcf" ]; then
        echo "Get variants file form NCBI"
	wget "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/common_all_20170710.vcf.gz"
fi
# If they are not unzipped, unzip them
if [ ! -f "common_all_20170710.vcf" ]; then
        echo "unzip"
	gzip -d "common_all_20170710.vcf.gz"
fi

if [ ! -f "chrallvars.vcf" ]; then
    mv common_all_20170710.vcf chrallvars.vcf
fi

# Filter the chromosomes even more by choosing one out of 10 that are not in weak-spots of the reference genome
if [ ! -f "chrallvars_filtered.vcf" ]; then
        echo "Select variants"
	sel_var.py "chrall.fasta" "chrallvars.vcf" "chrallvars_filtered.vcf" "${seed}"
fi

# Split the reference genome into paternal/maternal chromosomes and insert the variants into them
echo "split Genome into paternal/maternal chromosomes"
java -jar $VCF2DIPLOID -id "chrall" -seed "${seed}" -nochains -chr "chrall.fasta" -vcf "chrallvars_filtered.vcf"

# Simulate reads for each reference genome
for paternal in *chrall_paternal.fa
do
    prefix=$(echo ${paternal} | sed 's/^\(.*\)_chrall_paternal.fa/\1/')
    echo "prefix=${prefix}"
    $ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 15 -rs "${seed}" -na -o "paternal_${prefix}_chrall_PE" -i "${paternal}" &
done

wait_jobs
rm -rf *chrall_paternal.fa

for maternal in *chrall_maternal.fa
do
    prefix=$(echo ${maternal} | sed 's/^\(.*\)_chrall_maternal.fa/\1/')
    echo "prefix=${prefix}"
    $ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 15 -rs "${seed}" -na -o "maternal_${prefix}_chrall_PE" -i "${maternal}" &
done

wait_jobs
rm -rf *chrall_maternal.fa

# Combine all the reads into one pair
echo "Combine reads into a single pair of files"
for ((ii=1; ii<=24; ii++))
do
    combineFastq.py "maternal_chr${ii}_chrall_PE" "paternal_chr${ii}_chrall_PE" "chr${ii}_chrall_PE" "${seed}" &
done

echo "first round!"
wait_jobs

rm -rf *.map
rm -rf maternal* paternal*


for ((ii=1; ii<=24; ii+=2))
do
    combineFastq.py "chr${ii}_chrall_PE" "chr$((${ii}+1))_chrall_PE" "chr${ii}_$((${ii}+1))_chrall_PE" "${seed}" &
done

echo "second round!"
wait_jobs

for ((ii=1; ii<=24; ii+=1))
do
    rm -rf chr${ii}_chrall_PE*
done

for ((ii=1; ii<=24; ii+=4))
do
    combineFastq.py "chr${ii}_$((${ii}+1))_chrall_PE" "chr$((${ii}+2))_$((${ii}+3))_chrall_PE" "chr${ii}_$((${ii}+3))_chrall_PE" "${seed}" &
done

echo "third round!"
wait_jobs

for ((ii=1; ii<=24; ii+=2))
do
    rm -rf chr${ii}_$((${ii}+1))_chrall_PE*
done

for ((ii=1; ii<=24; ii+=8))
do
    combineFastq.py "chr${ii}_$((${ii}+3))_chrall_PE" "chr$((${ii}+4))_$((${ii}+7))_chrall_PE" "chr${ii}_$((${ii}+7))_chrall_PE" "${seed}" &
done

echo "fourth round!"
wait_jobs

for ((ii=1; ii<=24; ii+=4))
do
    rm -rf chr${ii}_$((${ii}+3))_chrall_PE*
done

echo "fiveth round!"
combineFastq.py "chr1_8_chrall_PE" "chr9_16_chrall_PE" "chr1_16_chrall_PE" "${seed}"

rm -rf chr1_8_chrall_PE*
rm -rf chr9_16_chrall_PE*

echo "last round!"
combineFastq.py "chr1_16_chrall_PE" "chr17_24_chrall_PE" "chrall_PE" "${seed}"

rm -rf chr1_16_chrall_PE*
rm -rf chr17_24_chrall_PE*

mv "chrallvars_filtered.vcf" "chrallvars_ref.vcf"
# Print the output file names to the user
echo "Created chrall_PE1.fastq, chrall_PE2.fastq and chrallvars_ref.vcf"

################## END OF SIMULATION ##################
