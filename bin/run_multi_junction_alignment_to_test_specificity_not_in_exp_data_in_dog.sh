#! /bin/bash

#Developped by Nicolas Guillaudeux
#Dyliss team (INRIA Rennes Bretagne Atlantique/IRISA/Université de Rennes 1)
#09/2019

# This script run different part of the validation process: it extracts specific exons couples from predicted transcript by CG-alcode [Blanquart S., 2016]. Then, it align this exon couples with bam files and extract the correct alignment reads. After, it creates a file with the number of reads correctly mapped with the exon couples.
#Choose a $1 SPECIES and an junction extension in $2.

CODE=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/bin

#################################################
# To add a species name (e.g. "dog" or "clf", "mouse" or "mm", "human" or "hs")
SPECIES=$1
echo "The species choice is write: $SPECIES"

#################################################
# The localisation where BAM files where located.
#TISSUE=$3
TISSUE=/run/media/niguilla/Dyliss_NG/IGDR_data/bam #dog data

echo "All tissus are stored in: $TISSUE"

#################################################
# The repository where data are stored
#JUNCTION_TEST=$3
#JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/test_dog_without_all_new_annotation_junction/specific_junction_cga_not_in_ENS96_ENS98_UCSC_FEELNC.bed #dog data

#To the 253:
JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/test_dog_without_all_new_annotation_junction/analysis_253/specific_junction_cga_not_in_ens96_ens98_ucsc_igdr.bed

echo "Corresponding data are stored in : $JUNCTION_TEST"

#LIST_GENES=${JUNCTION_TEST}/list_312_tr_clf.txt #The list of transcripts used

#################################################
# The delta of the junction extension
#DELTA_EXTENSION=65
DELTA_EXTENSION=$2
echo "Delta chosen for the junction extension is: $DELTA_EXTENSION"

##################
# Main programm: #
##################
echo "----------------------"
echo ">>> START PROGRAMM <<<"
echo "----------------------"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-> $current_date_time"

DATE=`date +%Y-%m-%d`
OUTPUT_DIR="junction_validation_${SPECIES}_${DATE}"

mkdir ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

#python3 ${CODE}/extract_exon_junctions.py -p ${CGALCODE_PREDICTION_GTF} -r ${ENSEMBL_GTF} -t ${LIST_GENES} -d ${DELTA_EXTENSION} > output_run.txt

for tissue_file in ${TISSUE}/*
	do
		fullfilename=$(basename $tissue_file)
		#extension=${fullfilename##*.}
		filename=${fullfilename%.*}
		if [ $filename != "junction_validation_$1_${DATE}" ]
			then
				mkdir $filename
				echo ">>> WORK WITH ${filename}"
				echo "# STEP1: EXTRACTION OF SPECIFIC EXON COUPLES"
				echo "NOT REALIZED"

				cd $filename
				echo "# STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT"
				current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
				echo "-> $current_date_time"
				bedtools intersect -abam ${TISSUE}/${fullfilename} \
				-b ${JUNCTION_TEST} -wa -wb -split -bed \
				> output_BedToolsIntersect_specific_junction_$filename.bed
				
				echo "# STEP3: EXTRACTION OF READS CORRECTLY ALIGNED"
				current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
				echo "-> $current_date_time"
				python3 ${CODE}/exon_junction_validator.py \
				-a output_BedToolsIntersect_specific_junction_$filename.bed \
				-j ${JUNCTION_TEST} \
				> output_run_${filename}.txt

				awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' output_reads_mapped_on_specific_junctions.bed \
				> output_reads_mapped_on_specific_junctions_without_junction.bed
				cd ..
		fi
	done

#mkdir tissue
mkdir analysed

echo "# STEP4: COPY FILES WITH NUMBER OF READS MAPPED BY TISSUE IN A SAME DIRECTORY"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-> $current_date_time"
for tissue_dir in *
	do
		if [ -d ${tissue_dir} ]; then
			cd ${tissue_dir}
			if [ -f output_run*.txt ]; then
				cp output_run*.txt ../analysed/
			fi
			cd ..
		fi
	done

cd analysed/

echo "# STEP5: CREATE A FILE WITH THE READS NUMBER MAPPED WITH SPECIFIC EXON COUPLES"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-> $current_date_time"
for reads_mapped in *
	do
		echo "${reads_mapped::-4}" | cut -c12- > ${reads_mapped}_reads_mapped.txt
		if [ ! -f "junction_file.txt" ] && [ ! -s "junction_file.txt" ]; then
			touch junction_file.txt
			echo Junctions > junction_file.txt
			sed '1,3d' ${reads_mapped} | awk '{print $1}' >> junction_file.txt
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

echo "STEP6: CREATE A GRAPHIC PDF FILE"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-> $current_date_time"
python3 ${CODE}/nb_reads_mapped_graph.py \
-f nb_reads_mapped_tissue.csv \
-n 82

echo "--------------------"
echo ">>> End programm <<<"
echo "--------------------"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo ">>> $current_date_time"
