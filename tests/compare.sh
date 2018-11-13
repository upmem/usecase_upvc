#!/bin/bash
# $1: simu
# $2: fpga
# $3: chr22

gcc -O3 extract_res.c -o extract_res

rm -rf $1/* $2/* diff_$1_$2.txt
mkdir -p $1
mkdir -p $2

cp $3/$1.txt $1/

cp $3/$2.txt $2/

./extract_res $1/$1.txt
./extract_res $2/$2.txt

for I in $1/*\.txt\.*
do
    echo -ne "\rsorting ${I}     "
    cat $I | sort -nk 1 > $I.sorted
done

echo -ne "\n"

for I in $2/*\.txt\.*
do
    echo -ne "\rsorting ${I}     "
    cat $I | sort -nk 1 > $I.sorted
done

echo -ne "\n"

nb_file=$(ls -l $2/*\.txt\.*\.sorted | wc -l)

for ((i = 0; i < $nb_file; i++))
do
    echo -ne "\rcomparing ${i} / ${nb_file}"
    echo "###### ${i} ######" >> diff_$1_$2.txt
    diff -y --suppress-common-lines $1/*\.txt\.$i\.sorted $2/*\.txt\.$i\.sorted >> diff_$1_$2.txt
done
echo -ne "\nDONE!\n"

exit 0
