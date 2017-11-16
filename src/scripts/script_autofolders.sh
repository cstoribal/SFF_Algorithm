#!/bin/bash 
##include <stdio.h>

# 1 default params

init="../ProjSFF -D "
LIST_FOLD=$1

for fold in ${LIST_FOLD[*]}
  do
  LIST_SAMP=("/media/sf_F_DRIVE/Stage/Samples/"$fold*.txt)
  #echo $LIST_SAMP
  #echo $fold
  for samp in ${LIST_SAMP[*]}
    do
    linecode=$init""$samp" "$2
    echo $linecode
    eval $linecode
  done
done
