#!/bin/bash


#in each subfolder:
#1) locate files with ending _1.fastq 
#2) change ending to _2.fastq to exactly match read paris.
#3) setup trimmomatic for each pair with input/output according to read paris.


for S in ./S*; do
#locate all subfolder
	if [ -d "$S" ]; then
		cd "$S"
		# if subfolder == folder: move into subfolder.
		for f1 in *_1.fastq
		do
			f2=${f1%%_1.fastq}"_2.fastq"
			# create variable f2  which is f1 with ending "_2.fastq"

			java -jar ~/bin/trimmomatic/trimmomatic.jar \
			PE -phred33 -trimlog \
			${f1%%_1.fastq}.log \
			$f1 \
			$f2 \
			${f1%%_1.fastq}_trimmed_FP.fastq \
			${f1%%_1.fastq}_trimmed_FU.fastq \
			${f1%%_1.fastq}_trimmed_RP.fastq \
			${f1%%_1.fastq}_trimmed_RU.fastq \
			LEADING:10 TRAILING:10 SLIDINGWINDOW:8:15 MINLEN:150 AVGQUAL:18 \
			ILLUMINACLIP:/home/gvl/bin/trimmomatic/adapters/TruSeq3-PE.fa:2:10:3:1
			echo ------------------------------------------------
		done

		cd ..
	fi
done
