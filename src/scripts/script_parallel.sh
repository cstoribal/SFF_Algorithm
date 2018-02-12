#!/bin/bash 
##include <stdio.h>

#Only if usefull, for instance if this location is RAM
cp ./SFFDataMgmt /data/tmp/
cp ./ProjSFF     /data/tmp/


### The actual program
run_command() {
    prog_param="/data/tmp/SFFDataMgmt -M "
    prog_sff="/data/tmp/ProjSFF -optf "$2" -D "

    linecode=$prog_param""$1" -f "$3
    echo $linecode
    eval $linecode
    mkdir ${1%??????????} 
    cd ${1%??????????}
    mkdir $4
    cd $4 #enter subfolders. Caution, it is actually 3 leveled-subfolder
    
    linecode=$prog_sff""$3
    echo $linecode
    eval $linecode
    cd ..
    cd ..
    cd ..
    
}
export -f run_command

# Logging the output of each sample (all parameters)
convert_and_run() {
    sampname=${3##*/}
    sampname=${sampname%?????????}
    parallel -j 1 --resume-failed --joblog ./exp/log/$4/$sampname.log run_command ::: $1 ::: $2 ::: $3 ::: $4
}
export -f convert_and_run



#LIST_FOLD=( "customseries/" "Semireal/semireal_mb4/" "Semireal/semireal_mb2/" "Semireal/semireal_mb1/" "Semireal/semireal_mb3/" "Simulations/simu4/" "Simulations/simu5/" "Simulations/s96serie/" "Simulations/")
LIST_FOLD=( "Semireal/semireal_mb4/" "Semireal/semireal_mb1/")
LIST_PARAMS=("./exp/*.param.txt")
opt_param=" 0 "

mkdir "./exp"
mkdir "./exp/log"

for fold in ${LIST_FOLD[*]}
  do
  subdir=${fold//\//_}
  subdir=${subdir%?}
  LIST_SAMP=("/data/data/Stage/Samples/"$fold*.txt)
  mkdir ./exp/log/$subdir
  parallel -j 6 --resume-failed --joblog ./exp/log/$subdir.log convert_and_run ::: "${LIST_PARAMS[*]}" ::: $opt_param ::: ${LIST_SAMP[*]} ::: $subdir
done

cd ./exp
../concatenate.sh
cd ..
