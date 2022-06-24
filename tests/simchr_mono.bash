#!/bin/bash

set -x

# USAGE: simchr.bash chromNumber [rngSeed]
# If rngSeed is not specified, a random one will be chosen instead

# EDIT THESE TWO VARIABLES TO POINT TO THE PROPER BINARIES FOR VCF2DIPLOID AND ART_ILLUMINA
# You can obtain VCF2Diploid from github.com/moselhy/vcf2diploid, you need to run "make" after you clone it to create the jar file
# You can obtain ART from https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm, you just need to untar the binary package
GENOME_HOME="/home/upmemstaff/rjodin/work/usecases/genomee/tests"
VCF2DIPLOID="${GENOME_HOME}/vcf2diploid/vcf2diploid.jar"
ART_ILLUMINA="${GENOME_HOME}/art_bin_MountRainier/art_illumina"

# add in the $HOME/.bash_profile
# PATH="/Users/lavenier/Documents/Projets/UPMEM/upvc1.4/VCScripts:${PATH}"

wait_jobs() {
    for I in $(jobs -p)
    do
	wait $I
    done
}
################## START OF SIMULATION ##################

num=$1
# Get the user-defined seed from input
if [ $2 ]; then
	seed=$2
else
	seed=$RANDOM
fi

PATH=${GENOME_HOME}:$PATH python ${GENOME_HOME}/sel_var.py "chr${num}.fasta" "chr${num}vars.vcf" "chr${num}vars_filtered.vcf" "${seed}"

# Split the reference genome into paternal/maternal chromosomes and insert the variants into them
java -jar $VCF2DIPLOID -id "chr${num}" -seed "${seed}" -nochains -chr "chr${num}.fasta" -vcf "chr${num}vars_filtered.vcf"

# Simulate reads for each reference genome
for paternal in *chr${num}_paternal.fa
do
    prefix=$(echo ${paternal} | sed 's/^\(.*\)_[^_]*_paternal.fa/\1/')
    $ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 15 -rs "${seed}" -na -o "paternal_${prefix}_PE" -i "${paternal}" &
done

wait_jobs
rm -rf *chr${num}_paternal.fa

for maternal in *chr${num}_maternal.fa
do
    prefix=$(echo ${maternal} | sed 's/^\(.*\)_[^_]*_maternal.fa/\1/')
    $ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 15 -rs "${seed}" -na -o "maternal_${prefix}_PE" -i "${maternal}" &
done

wait_jobs
rm -rf *chr${num}_maternal.fa

PATH=${GENOME_HOME}:$PATH python ${GENOME_HOME}/combineFastq.py "maternal_chr${num}_PE" "paternal_chr${num}_PE" "chr${num}_PE" "${seed}"


mv "chr${num}vars_filtered.vcf" "chr${num}vars_ref.vcf"
# Print the output file names to the user
echo "Created chr${num}_PE1.fastq, chr${num}_PE2.fastq and chr${num}vars_ref.vcf"
rm -f maternal* paternal* *.map chr${num}vars.vcf

################## END OF SIMULATION ##################
