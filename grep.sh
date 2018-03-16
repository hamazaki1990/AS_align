#!/usr/bin/env bash

for FILE in 0 45 90 135
do
  cat DX1_1_B_${FILE}.list | while read line
  do
    grep "${line}" Btype_${FILE}.fa -A1
  done > "DX1_1_${FILE}_Btype.fa"
done
