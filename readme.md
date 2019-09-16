Common work Dominique Lavenier (INRIA) and uPmem.
=================================================

Installation
------------

```
cd <genomee>
mkdir build
cd build
cmake ..
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

Data must consist of 3 files:
  - ``<dataset_prefix>.fasta`` the reference genomee (or part of genomee)
  - ``<dataset_prefix>_PE1.fastq`` the PE1 of the input to compare to the reference
  - ``<dataset_prefix>_PE2.fastq`` the PE2 of the input to compare to the reference

Run once to create the MRAM for the reference genomee to compare to:

```
./build/host/upvc -i <dataset_prefix> -d <number_of_virtual_dpus_during_execution> -g index
```

Then run:

```
./build/host/upvc -i <dataset_prefix> -g map -n <number_of_physical_dpus_available> -t <target=fsim|fpga> -b <dpu_binary>
```

Result are in ``<dataset_prefix>_upvc.vcf``

Paper link
----------

https://drive.google.com/open?id=1NJYijXnVzEpsVQdHgtMZ7RoL4d_BW-CJ