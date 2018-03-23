#!/usr/bin/env bash

for chr in $(seq 1 22) X Y
do
  cat human_livingASs_labels.csv |grep "D${chr}Z"|awk 'BEGIN{print "AS"}{print $0}' > "chr${chr}_ASs_labels.csv"
done
