#!/bin/bash 
##include <stdio.h>

# 1 default params


init="../src/ChangeParameters.sh "
compute="../../src/script_autofolders.sh"
LIST_PARAMS=("./*.param.txt")

for param in ${LIST_PARAMS[*]}
  do
    linecode=$init""$param
    echo $linecode
    eval $linecode
    mkdir ${param%??????????}
    cd ${param%??????????}
    linecode=$compute
    echo $linecode
    eval $linecode
    cd ..
  done

