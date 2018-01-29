#!/bin/bash 
##include <stdio.h>


cp ./SFFDataMgmt /data/tmp/
cp ./ProjSFF     /data/tmp/

convert_and_run() {
    prog_param="/data/tmp/SFFDataMgmt -M "
    prog_sff="/data/tmp/ProjSFF -optf "$2" -D "
    #prog_param="./SFFDataMgmt -M "
    #prog_sff="../../ProjSFF -optf "$2" -D "
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
      cd ..
    done
}

export -f convert_and_run



#LIST_FOLD=( "customseries/" "Semireal/semireal_mb4/" "Semireal/semireal_mb2/" "Semireal/semireal_mb1/" "Semireal/semireal_mb3/" "Simulations/simu4/" "Simulations/simu5/" "Simulations/s96serie/" "Simulations/")
LIST_FOLD=( "Semireal/semireal_mb1/" "Semireal/semireal_mb4/")
#LIST_FOLD=( "customseries/" )
LIST_PARAMS=("./exp/*.param.txt")
opt_param=" 0 "

mkdir "./exp"
mkdir "./exp/log"


for fold in ${LIST_FOLD[*]}
  do
  subdir=${fold//\//_}
  subdir=${subdir%?}
  LIST_SAMP=("/data/data/Stage/Samples/"$fold*.txt)
  parallel -j 8 --resume-failed --joblog ./exp/log/$subdir.log convert_and_run ::: "${LIST_PARAMS[*]}" ::: $opt_param ::: ${LIST_SAMP[*]} ::: $subdir
done

cd ./exp
../concatenate.sh
cd ..
