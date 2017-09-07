#!/bin/bash 
##include <stdio.h>

# 1 default params

init="./ProjSFF -D "
LIST_FOLD=("Simulations/s96serie/" "Simulations/" "" "Simulations/semireal_1/" "Simulations/semireal_2/")

for fold in ${LIST_FOLD[*]}
  do
  LIST_SAMP=("../Samples/"$fold*.txt)
  #echo $LIST_SAMP
  #echo $fold
  for samp in ${LIST_SAMP[*]}
    do
    linecode=$init""$samp
    echo $linecode
    eval $linecode
  done
done

