#!/bin/bash


for D in ./S*; do
	if [ -d "$D" ]; then
		cd "$D"
		for S in ./S*; do
			if [ -d "$S" ]; then
				cd "$S"
				prodigal \
				-a ${S##*/}N.faa \
				-d ${S##*/}N.cDNA.fna \
				-f gff \
				-i ./scaffolds.fasta \
				-m \
				-c \
				-t ${S##*/}N.training \
				-o ${S##*/N.gff}

				prodigal \
				-a ${S##*/}N.faa \
				-d ${S##*/}N.cDNA.fna \
				-f gff \
				-i ./scaffolds.fasta \
				-m \
				-c \
				-t ${S##*/}N.training \
				-o ${S##*/N.gff}
				cd ..
				
			fi
		done
		cd ..
	fi
done
notify-send "Done!"