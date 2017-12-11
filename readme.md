Common work Dominique Lavenier (INRIA) and uPmem.

=========================================================
version 1.0
=========================================================
11-12-2017


Installation

- get the *.c and upvc.h files
- run the Makefile
- the ./upvc binary file mist be generated

Run

- ./upvc name_dataset

  the dataset includes the following files :
    - name_dataset.fasta      ==> reference genome
    - name_dataset_PE1.fastq  ==> file pair-read 1
    - name_dataset_PE2.fastq  ==> file pair-read 2

  the program generates the file:
    - name_dataset.vcf        ==> list of small variants

