#!/bin/bash

mkdir scaffolds;

for D in ./S*; do
	if [ -d "$D" ]; then
		cd "$D"
		for S in ./S*; do
			if [ -d "$S" ]; then
				cd "$S"
				
				FASTA="scaffolds.fasta"
				FILE="$(echo $S | cut -d "-" -f2|cut -d "_" -f1,2)_scaffolds.fasta"
				#echo "$FASTA"
				#echo "$FILE"
				cp "$FASTA" ../../scaffolds/$FILE

				cd ..
				
			fi
		done
		cd ..
	fi
done
notify-send "Done!"