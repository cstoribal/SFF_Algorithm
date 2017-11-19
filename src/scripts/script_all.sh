#!/bin/bash 
##include <stdio.h>

# 1 default params



LIST_FOLD=( "customseries/" "Semireal/semireal_mb4/" "Semireal/semireal_mb2/" "Semireal/semireal_mb1/" "Semireal/semireal_mb3/" "Simulations/simu4/" "Simulations/simu5/" "Simulations/s96serie/" "Simulations/")
# LIST_FOLD=( "customseries/") 
LIST_PARAMS=("./*.param.txt")
opt_param=" -optf 0 "

prog_param="./SFFDataMgmt -M "
prog_sff="../../ProjSFF "$opt_param" -D "

for fold in ${LIST_FOLD[*]}
  do
  subdir=${fold//\//_}
  ##echo $subdir
  ##mkdir $subdir
  ##cd $subdir
  LIST_SAMP=("/media/sf_F_DRIVE/Stage/Samples/"$fold*.txt)
  for samp in ${LIST_SAMP[*]}
    do
    for param in ${LIST_PARAMS[*]}
      do
      linecode=$prog_param""$param" -f "$samp
      echo $linecode
      #eval $linecode
      mkdir ${param%??????????}
      cd ${param%??????????}
      mkdir $subdir
      cd $subdir
      
      linecode=$prog_sff""$samp
      echo $linecode
      #eval $linecode
      cd ..
      cd ..
    done
  done
done
