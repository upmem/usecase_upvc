#!/bin/bash

acc=0
nb_line=0
cat $1 | while read line
do
    val=$(echo $line | sed 's/^\([^,]*\),.*/\1/')
    acc=$(($acc+$val))
    nb_line=$(($nb_line+1))
    moy=$(($acc/$nb_line))
    echo -ne "\r$acc $nb_line $moy"
done

echo -ne "\nfinish\n"
