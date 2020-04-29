#!/bin/bash


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

$DIR/convert_data $1 $2 $3 $4

ERR=$?

if [[ $ERR -eq 0 ]];
then

if [[ $# -eq 4 ]];
then
  OUTDIR=$4
  cat $OUTDIR/data_*.bin > $OUTDIR/data.bin
  rm $OUTDIR/data_*.bin
elif [[ $# -eq 3 ]];
then
  OUTDIR=$3
  cat $OUTDIR/data_*.bin > $OUTDIR/data.bin
  rm $OUTDIR/data_*.bin
fi

fi
