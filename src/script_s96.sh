#!/bin/bash 
##include <stdio.h>

# 1 default params

init="./ProjSFF -D ../../Samples/Simulations/s96serie/"
LIST_SAMP=("s96_cone04_data.txt" "s96_cone_disc360_data.txt" "s96_cone_strip180_data.txt" "s96_cos04_data.txt" "s96_cos_disc360_data.txt" "s96_cos_strip180_data.txt" "s96_plane2_disc360_data.txt" "s96_plane2_strip180_data.txt" "s96_plane2v04_data.txt" "s96_plane04_data.txt" "s96_plane_disc360_data.txt" "s96_plane_strip180_data.txt" "s96_sphere04_data.txt" "s96_sphere_disc360_data.txt" "s96_sphere_strip180_data.txt")

for samp in ${LIST_SAMP[*]}
  do
  linecode=$init""$samp
  eval $linecode
done

