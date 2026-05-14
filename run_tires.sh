#!/bin/bash

fileInp='input.para'
fileDat='data.para'

## ----------------------------------------------
## ----------------------------------------------
## Read working directory from fileInp
n=0
while read line; do
  n=$((n+1))
  # read number of zones
  if [ $n = 8 ]; then
    nzones=$line
  fi
  wline=$((33+(5*nzones)))
  # read working directory name
  if [ $n = $wline ]; then
    wdir=$line
  fi
done < $fileInp
[ "${wdir: -1}" != "/" ] && wdir=${wdir}'/'
echo WorkDir: $wdir

## ----------------------------------------------
## Prepare time stamp
tstamp=$(date +'%y%m%d_%H%M%S')
echo TimeStamp: $tstamp

## ----------------------------------------------
## Run the afterprocessing of Trans-dimensional Inversion
nohup ./tires $fileInp $fileDat > ${wdir}${tstamp}.log 2>&1 &

## ----------------------------------------------
## PID to file
echo PID: $!
echo $! > ${wdir}${tstamp}.pid