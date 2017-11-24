#!/bin/bash 
##include <stdio.h>


convert_and_run() {
    prog_param="./SFFDataMgmt -M "
    prog_sff="../../ProjSFF "$2" -D "
    for param in $1
      do
      linecode=$prog_param""$param" -f "$3
      echo $linecode
      eval $linecode
      mkdir ${param%??????????}
      cd ${param%??????????}
      mkdir $4
      cd $4
      
      linecode=$prog_sff""$3
      echo $linecode
      eval $linecode
      cd ..
      cd ..
    done
}

export -f convert_and_run



LIST_FOLD=( "customseries/" "Semireal/semireal_mb4/" "Semireal/semireal_mb2/" "Semireal/semireal_mb1/" "Semireal/semireal_mb3/" "Simulations/simu4/" "Simulations/simu5/" "Simulations/s96serie/" "Simulations/")
LIST_PARAMS=("./*.param.txt")
opt_param=" -optf 0 "



for fold in ${LIST_FOLD[*]}
  do
  subdir=${fold//\//_}
  LIST_SAMP=("/data/data/Stage/Samples/"$fold*.txt)
  parallel convert_and_run ::: "${LIST_PARAMS[*]}" ::: $opt_param ::: ${LIST_SAMP[*]} ::: $subdir
done


