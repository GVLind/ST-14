#!/bin/bash

#Spades algo - only one folder.
for S in ./S*; do
	if [ -d "$S" ]; then
		cd "$S"

		for f1 in *_trimmed_FP.fastq
		do
			notify-send "Start $f1"
			f2=${f1%%_trimmed_FP.fastq}"_trimmed_RP.fastq"

			/home/gvl/Programs/SPAdes-3.10.0-Linux/bin/spades.py \
			--pe1-1 $f1 \
			--pe1-2 $f2 \
			-k 21,33,55,77 \
			--careful \
			-t 32 \
			-o SpadesN_${f1/_trimmed*/}
		done

		cd ..
	fi
done

