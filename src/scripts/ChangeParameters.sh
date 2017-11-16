#!/bin/bash 
##include <stdio.h>

# 1 default params


init="./SFFDataMgmt -M "$1" -f "
LIST_FOLD=$2

for fold in ${LIST_FOLD[*]}
  do
  LIST_SAMP=("/media/sf_F_DRIVE/Stage/Samples/"$fold*.txt)
  #echo $LIST_SAMP
  #echo $fold
  for samp in ${LIST_SAMP[*]}
    do
    linecode=$init""$samp
    echo $linecode
    eval $linecode
  done
done

