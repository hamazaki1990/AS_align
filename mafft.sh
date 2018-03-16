#!/bin/sh

for Start in 0 45 90 135
do
  echo $Start
  mafft --auto DX1_1_${Start}_Btype.fa > DX1_1_B_${Start}_align.fasta
done
