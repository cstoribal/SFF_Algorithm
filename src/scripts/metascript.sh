#!/bin/bash 
##include <stdio.h>

# 1 default params


init="./ChangeParameters.sh "
compute="../script_autofolders.sh"
LIST_FOLD=("Semireal/semireal_mb2/") ##Semireal/semireal_mb2/
LIST_PARAMS=("./*.param.txt")
opt_param=" -optf 0 "

for param in ${LIST_PARAMS[*]}
  do
    linecode=$init""$param' "'${LIST_FOLD[*]}'"'
    echo $linecode
    eval $linecode
    mkdir ${param%??????????}
    cd ${param%??????????}
    linecode=$compute' "'${LIST_FOLD[*]}'" "'$opt_param'"'
    echo $linecode
    eval $linecode
    cd ..
  done

