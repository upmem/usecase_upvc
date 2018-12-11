#!/bin/bash

# USAGE: simchr.bash chromNumber [rngSeed]
# If rngSeed is not specified, a random one will be chosen instead

# EDIT THESE TWO VARIABLES TO POINT TO THE PROPER BINARIES FOR VCF2DIPLOID AND ART_ILLUMINA
# You can obtain VCF2Diploid from github.com/moselhy/vcf2diploid, you need to run "make" after you clone it to create the jar file
# You can obtain ART from https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm, you just need to untar the binary package
VCF2DIPLOID="/home/rjodin/work/genomee/vcf2diploid/vcf2diploid.jar"
ART_ILLUMINA="/home/rjodin/work/genomee/art_bin_MountRainier/art_illumina"

# add in the $HOME/.bash_profile
# PATH="/Users/lavenier/Documents/Projets/UPMEM/upvc1.4/VCScripts:${PATH}"


if [ $(hostname -f | cut -d . -f2,3) == "genouest.org" ]; then
	source /local/env/envjava-1.8.0.sh
fi


################## START OF SIMULATION ##################

# Get the user-defined seed from input
if [ $2 ]; then
	seed=$2
else
	seed=$RANDOM
fi

# If the reference genome does not already exist on disk, download it
if [ ! -f "chr${1}.fasta" ]; then
        echo "get chromosome"
	wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chr${1}.fa.gz"
        echo "unzip"
	gzip -d "hs_ref_GRCh38.p12_chr${1}.fa.gz"
	mv "hs_ref_GRCh38.p12_chr${1}.fa" "chr${1}.fasta"
fi

# If the first line of the reference genome (header line) is not ">chr" then the chromosome number, make it so
if [ "$(head -n1 chr${1}.fasta)" != ">chr${1}" ]; then
#	sed -i "1s/.*/>chr${1}/g" "hs_ref_GRCh38.p12_chr${1}.fa"
# modif MAC OS X
        sed -i '' "1s/.*/>chr${1}/g" "chr${1}.fasta"
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

# Make the command that we will use with AWK to filter variants of the given chromosome
awkcmd='$1 ~ /^#/ || $1 == '
awkcmd+="$1 "

# If the variants are not already filtered, then filter them using AWK
if [ ! -f "chr${1}vars.vcf" ]; then
        echo "filter variant with AWK:" $awkcmd
	awk "$awkcmd" common_all_20170710.vcf > chr${1}vars.vcf
fi

# Filter the chromosomes even more by choosing one out of 10 that are not in weak-spots of the reference genome
if [ ! -f "chr${1}vars_filtered.vcf" ]; then
        echo "Select variants"
	sel_var.py "chr${1}.fasta" "chr${1}vars.vcf" "chr${1}vars_filtered.vcf" "${seed}"
fi

# Split the reference genome into paternal/maternal chromosomes and insert the variants into them
echo "split Genome into paternal/maternal chromosomes"
java -jar $VCF2DIPLOID -id "chr${1}" -seed "${seed}" -nochains -chr "chr${1}.fasta" -vcf "chr${1}vars_filtered.vcf"

# Simulate reads for each reference genome
$ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 25 -rs "${seed}" -na -o "paternal_chr${1}_PE" -i "chr${1}_chr${1}_paternal.fa"
$ART_ILLUMINA -m 400 -s 50 -l 120 -p -f 25 -rs "${seed}" -na -o "maternal_chr${1}_PE" -i "chr${1}_chr${1}_maternal.fa"

# Combine all the reads into one pair
echo "Combine reads into a single pair of files"
combineFastq.py "maternal_chr${1}_PE" "paternal_chr${1}_PE" "chr${1}_PE" "${seed}"

mv "chr${1}vars_filtered.vcf" "chr${1}vars_ref.vcf"
# Print the output file names to the user
echo "Created chr${1}.fasta, chr${1}_PE1.fastq, chr${1}_PE2.fastq and chr${1}vars_ref.vcf"

################## END OF SIMULATION ##################
