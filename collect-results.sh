#!/bin/bash
echo '' > all-results.txt
for f in $(ls */results/variants-comparison-freebayes.txt); do
  echo $f >> all-results.txt
  cat $f >> all-results.txt
  echo '' >> all-results.txt
done
