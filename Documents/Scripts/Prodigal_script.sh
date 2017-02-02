#!/bin/bash

# add correct path.


for S in ./S*; do
	if [ -d "$S" ]; then
		cd "$S"
		mkdir prodigal_${S##./}
		cd ./prodigal*
		prodigal \
		-a ${S##./}.faa \
		-d ${S##./}cDNA.fna \
		-f gff \
		-i /home/gvl/Bioinformatics_Test/TestRun_Prodigal_script/prodgial_test/S5/SpadesN_S5_spades_output/scaffolds.fasta \
		-m -c \
		-t $S.training \
		-o $S.gff
		
		prodigal \
		-a $S_.faa \
		-d $S_cDNA.fna \
		-f gff \
		-i /home/gvl/Bioinformatics_Test/TestRun_Prodigal_script/prodgial_test/S5/SpadesN_S5_spades_output/scaffolds.fasta \
		-m -c \
		-t $S.training \
		-o $S.gff

		cd ../..
	fi
done

