#!/usr/bin/env bash

for chr in $(seq 1 22) X Y
do
  cat "chr${chr}_ASs_labels.csv" |awk '{if (substr($1,1,2)=="D1") print "B"}'|head
done
