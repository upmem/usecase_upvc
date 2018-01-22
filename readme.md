Common work Dominique Lavenier (INRIA) and uPmem.

=========================================================
version 1.0
=========================================================
11-12-2017
==========

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


==========
17-01-2018
==========

- correct a bug in get_genome (getgenome.c)
- lower NB_DPU to 128 (upvc.h)
- modify the function align (upvc_dpu.c)
   - replace ODPD by noDP
   - suppress decoding step

==========
22-01-2018
==========

- mixe the use of noDP et DDPD functions in the align function (upvc_dpu.c)

