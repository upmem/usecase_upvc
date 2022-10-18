UPVC: Common work Dominique Lavenier (INRIA) and UPMEM.
=======================================================

Installation
------------

```
cmake -B build
cd build
make
```

Run on integration dataset
--------------------------

First generate the index with:

```
make mrams
```

You can then run on the integration dataset with:

```
make run
make check
```

Or, if not running on upmem DIMMs, you can instead run the integration dataset with:

```
make run_simu
make check
```

Run on other dataset
--------------------

Data-set must consist of 3 files:
  - ``<dataset_prefix>.fasta`` the reference genomee (or part of genomee)
  - ``<dataset_prefix>_PE1.fastq`` the PE1 of the input to compare to the reference
  - ``<dataset_prefix>_PE2.fastq`` the PE2 of the input to compare to the reference
  - ``<reference_vcf>`` the reference vcf output to check the quality of the computation

Run once to create the MRAM for the reference genomee to compare to:

```
./<path_to_build>/host/upvc -i <dataset_prefix> -n <number_of_virtual_dpus_during_execution> -g index
```

Then run:

```
./<path_to_build>/host/upvc -i <dataset_prefix> -g map  [-n <number_of_physical_dpus_available>]
```

If the number of physical dpus available is not specified, the program will try to alloc every dpus available at runtime.

Results are in ``<dataset_prefix>_upvc.vcf``

To check the quality of the results use:

```
python3 <path_to_build>/../tests/compareVCF.py <reference_vcf> <dataset_prefix>_upvc.vcf
```
