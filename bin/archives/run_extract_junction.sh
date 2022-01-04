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
# The repository where data are stored
#JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons
#JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_mouse_with_ensembl_data_in_ref
#JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_mouse_with_ensembl98_data_in_ref
JUNCTION_TEST=/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_human_with_ensembl98_data_in_ref

echo "Corresponding data are stored in : $JUNCTION_TEST"

# The CG-alcode gtf file prediction with all prediction information
#CGALCODE_PREDICTION_GTF=${JUNCTION_TEST}/dog_cgalcode_tot.gtf
#CGALCODE_PREDICTION_GTF=${JUNCTION_TEST}/mouse_cgalcode_tot.gtf
CGALCODE_PREDICTION_GTF=${JUNCTION_TEST}/human_cgalcode_tot.gtf

# The gtf file with all data used in CG-alcode programm with all known information from a reference source (e.g. Ensembl)
#ENSEMBL_GTF=${JUNCTION_TEST}/dog_ensembl_tot.gtf
#ENSEMBL_GTF=${JUNCTION_TEST}/mouse_ensembl_tot.gtf
ENSEMBL_GTF=${JUNCTION_TEST}/human_ensembl_tot.gtf

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

echo "# STEP: EXTRACTION OF SPECIFIC EXON COUPLES"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-> $current_date_time"
python3 ${CODE}/extract_specific_junctions.py \
-p ${CGALCODE_PREDICTION_GTF} \
-r ${ENSEMBL_GTF} \
-d ${DELTA_EXTENSION} \
> output_run.txt


echo "--------------------"
echo ">>> End programm <<<"
echo "--------------------"
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo ">>> $current_date_time"
