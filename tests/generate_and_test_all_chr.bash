#!/bin/bash

set -x

GENOMEE=$(realpath $1)
CHRALL=$(realpath $2)
RESULT=$(pwd)/result.txt

rm -rf ${RESULT}

for ((c = 1; c < 25; c++))
do
    echo $c >> $RESULT
    date >> $RESULT
    rm -rf chr${c}_640
    mkdir -p chr${c}_640
    pushd chr${c}_640
    mv ${CHRALL}/chr${c}.fasta .
    cat ${CHRALL}/chrallvars.vcf | grep -E "^${c}\>" > chr${c}vars.vcf
    date >> $RESULT
    time ${GENOMEE}/tests/simchr_mono.bash ${c}
    date >> $RESULT
    time ${GENOMEE}/build/host/upvc -i chr${c} -g index -n 640
    date >> $RESULT
    time ${GENOMEE}/build/host/upvc -i chr${c} -g map
    date >> $RESULT
    python ${GENOMEE}/tests/compareVCF.py *.vcf >> ${RESULT}
    popd
done
