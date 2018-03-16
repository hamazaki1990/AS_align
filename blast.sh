#!/bin/bash -f

for db in $(seq 1 4)
do
  for AS in $(seq 1 11)
  do
    echo $db
    echo $AS
    blastn -db "hg38_D1Z5_${db}" -query "D1Z5_AS${AS}.fa" -outfmt "6" -out "D1Z5_peri${db}_AS${AS}blast.tsv"
  done
done
