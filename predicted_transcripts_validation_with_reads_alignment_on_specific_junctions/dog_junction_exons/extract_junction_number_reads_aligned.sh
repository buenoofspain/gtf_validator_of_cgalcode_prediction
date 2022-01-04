#! /bin/bash

#Developped by Nicolas Guillaudeux
#Dyliss team (INRIA Rennes Bretagne Atlantique/IRISA/Université de Rennes 1)

WORK_DIRECTORY=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/Tissue_analysed

for tissue_dir in ${WORK_DIRECTORY}/*
	do
		fullfilename=$(basename $tissue_dir)
		#echo ${fullfilename}
		cp ${tissue_dir}/output_*_sortByCoord.2To4Exons.out ../analysed/
		
	done

cd ${WORK_DIRECTORY}/analysed/

for reads_mapped in ${WORK_DIRECTORY}/analysed/*
	do
		fullfilename=$(basename $reads_mapped)
		filename=${fullfilename%.*}
		echo "${fullfilename::-26}" | cut -c8- > ${reads_mapped}_reads_mapped.txt
		if [ ! -f "junction_file.txt" ] && [ ! -s "junction_file.txt" ]; then
			touch junction_file.txt
			sed '1,2d' ${reads_mapped} | awk '{print $1}' > junction_file.txt
			#paste junction_file.txt ${reads_mapped}_reads_mapped.txt > nb_reads_mapped_tissue.csv
			#echo "fichier vide"
		fi
#		${reads_mapped} > ${reads_mapped}_reads_mapped.txt
		sed '1,3d' ${reads_mapped} | awk '{print $2}' >> ${reads_mapped}_reads_mapped.txt
		rm ${reads_mapped}
		#paste nb_reads_mapped_tissue.csv ${reads_mapped}_reads_mapped.txt >> nb_reads_mapped_tissue.csv
		#echo ${reads_mapped}
				#sed '1,2d' ${reads_mapped}
	done

paste * > nb_reads_mapped_tissue.csv
