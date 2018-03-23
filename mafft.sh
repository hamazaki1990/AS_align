#!/bin/sh

for Start in 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150
do
  echo $Start
  mafft --auto AS_30bps${Start}.fa > AS_30bps${Start}_align.fasta
done
