#!/bin/bash

fileInp='input.para'
fileDat='data.para'
nnodes=16

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
## Create working directory
if [ ! -d "$wdir" ]; then
  mkdir -p "$wdir"
fi
## Create log directory
if [ ! -d "${wdir}log" ]; then
  mkdir -p "${wdir}log"
fi

## ----------------------------------------------
## Copy input files to the working directory
cp $fileInp ${wdir}${fileInp##*/}
cp $fileDat ${wdir}${fileDat##*/}

## ----------------------------------------------
## Prepare time stamp
tstamp=$(date +'%y%m%d_%H%M%S')
echo TimeStamp: $tstamp

## ----------------------------------------------
## Run Trans-D Inversion
if [ $nnodes = 1 ]; then
  # Run TI on single CPU
  #./tiser input.para data.para
  nohup ./tiser $fileInp $fileDat > ${wdir}${tstamp}.log 2>&1 &
else
  # Run TI on multiple CPU with MPI
  #mpirun -np 12 ./timpi input.para data.para
  nohup mpirun -np $nnodes ./timpi $fileInp $fileDat > ${wdir}${tstamp}.log 2>&1 &
fi

## ----------------------------------------------
## PID to file
echo PID: $!
echo $! > ${wdir}${tstamp}.pid