#!/bin/bash

wdir="./inv"      # Working directory to be cleaned
deleteEnsemble=0  # 1 = delete ensemble of solutions (xmodels*.dat)

## ----------------------------------------------
## ----------------------------------------------
## Check working directory
if [ ! -d "$wdir" ]; then
    exit 1
fi

## ----------------------------------------------
## Delete temporary files
find "$wdir" -maxdepth 1 -type f -name "node*" -delete -print

## ----------------------------------------------
## Delete ensemble of solutions
if [ "$deleteEnsemble" -eq 1 ]; then
    find "$wdir" -maxdepth 1 -type f -name "xmodels*.dat" -delete -print
fi