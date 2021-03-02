#!/bin/bash
set -e -x 
tail -n+11 ${1}.txt | awk '{print "0",$2,$1,$3,$4}' > ${1}.lgen &
tail -n+11 ${1}.txt | grep -P "\t${2}\t" | awk '{print $6,$1,"0",$7}' > ${1}.map &
tail -n+11 ${1}.txt | awk '{print "0",$2,"0","0","0","0"}' | uniq > ${1}.fam &

wait
sed -i 's/ -/ 0/g' ${1}.lgen &
wait
