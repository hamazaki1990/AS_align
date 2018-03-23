#!/bin/bash -f

for db in 150
do
  for AS in 140
  do
    echo $db
    echo $AS
    blastn -subject "AS_30bps${db}.fa" -query "AS_30bps${AS}.fa"
  done >"AS_30bps_blast${AS}.txt"
done
