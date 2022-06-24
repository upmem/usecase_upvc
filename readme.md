UPVC: Common work Dominique Lavenier (INRIA) and UPMEM.
=======================================================

UPVC performs small paired reads mapping (~150 nucleotides per read) onto a reference genome followed by variant calling on the Host CPU.
UPVC takes advantage of the high PIM parallelism to complete the mapping stage. It consists in finding all potential positions of short reads on the map of the reference genome loaded in memory. The CPU relieved from the data heavy computation is available to perform advanced dynamic computing and complex tasks such as alignment with variant calling. 

The input dataset is composed of an indexed human reference genome generated from a fasta text file, along with two fastq files containing the paired reads and some metadata. UPVC outputs a text file listing called variants in vcf format.

This program was developped by UPMEM team. Reach us at contact@upmem.com if you would like more details about this implementation (workflow structure, benchmarks, etc.).

Installation
------------

```
cmake -B build
cd build
make
```

Run on integration dataset
--------------------------

```
make run
make check
```

Run on other dataset
--------------------

Data-set must consist of 3 files:
  - ``<dataset_prefix>.fasta`` the reference genomee (or part of genomee)
  - ``<dataset_prefix>_PE1.fastq`` the PE1 of the input to compare to the reference
  - ``<dataset_prefix>_PE2.fastq`` the PE2 of the input to compare to the reference
  - ``<reference_vcf> the reference vcf output to check the quality of the computation

Run once to create the MRAM for the reference genomee to compare to:

```
./<path_to_build>/host/upvc -i <dataset_prefix> -n <number_of_virtual_dpus_during_execution> -g index
```

Then run:

```
./<path_to_build>/host/upvc -i <dataset_prefix> -g map  [-n <number_of_physical_dpus_available>]
```

If the number of physical dpus available is not specified, the program will try to alloc every dpus available at runtime.

Result are in ``<dataset_prefix>_upvc.vcf``

To check the quality of the results use:

```
python3 <path_to_build>/../tests/compareVCF.py <reference_vcf> <dataset_prefix>_upvc.vcf
```

Paper link
----------

https://hal.archives-ouvertes.fr/hal-03006764/document
