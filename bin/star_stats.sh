#!/bin/bash

files=$(find ./logs -name "*Log.final.out")
declare -a stats=( "Number of input reads"  "Average input read length" "Uniquely mapped reads number" "Uniquely mapped reads %" "Average mapped length" );
out="star_stats.tab"

printf "%s\t" "ID" > $out
printf "%s\t" "${stats[@]}" >> $out
printf "\n" >> $out

for file in $files; do
	printf ${file##*/}"\t" >> $out
	for stat in  "${stats[@]}"; do  
		grep "$stat" $file | awk -F $'\t' '{ printf "%f", $2; printf "\t"} ' >> $out
	done
printf "\n" >> $out
done
