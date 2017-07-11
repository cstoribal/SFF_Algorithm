#!/bin/bash 
##include <stdio.h>

# 1 default params

init="./ProjSFF -D ../../Samples/"
LIST_SAMP=("Simulations/sim_cone_01_data.txt" "Simulations/sim_cos_01_data.txt" "Simulations/sim_sphere_01_data.txt" "Simulations/sim_plane_01_data.txt" "Simulations/simu4/sim_cone_04_data.txt" "Simulations/simu4/sim_cos_04_data.txt" "Simulations/simu4/sim_plane2_04_data.txt" "Simulations/simu4/sim_plane_04_data.txt" "Simulations/simu4/sim_sphere_04_data.txt")

for samp in ${LIST_SAMP[*]}
  do
  linecode=$init""$samp
  eval $linecode
done

