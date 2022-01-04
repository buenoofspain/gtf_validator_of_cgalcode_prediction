# Script descriptions

## gtf_CGalcode_generator_and_validator_predictions.py

This script convert genes and transcripts in CG-Alcode format to a GTF file. Example:
`python3 gtf_CGalcode_generator_and_validator_predictions.py`
`python3 ../../bin/gtf_CGalcode_generator_and_validator_predictions.py -rep ../../examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/ -c ../../examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_mm_withprediction.csv -set ../../examples/set_genes/set_253_mouse.txt -ref ../../examples/reference_gtf/mouse/ensembl/Ensembl96/Mus_musculus.GRCm38.96.gtf`

### arguments: 
- "`-gtf2bed`", "`--gtf2bedconvertion`", Give the Gtf to Bed script converter. By default, use the perl script, default = "dependencies/gtf2bed.pl".
- "`-rep`", "`--repository`", Enter the path repository where are located genes and transcripts files used by CG-Alcode programm, default = "examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/".
- "`-c`", "`--correspondance`", Give the file with the description relationship between genes and their transcripts, default = "examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_hs_withprediction.csv".
- "`-set`", "`--setgene`", Give a list of interested genes, default = "examples/set_genes/c135_human.txt".
- "`-ref`", "`--reference`", Give a reference gtf file to validate prediction realize with CG-Alcode programm, default = "examples/reference_gtf/human/ensembl/Ensembl98/Homo_sapiens.GRCh38.98.gtf".

### Information to the script execution with another data not Ensembl:
This script is executable with Ensembl data, for another data (UCSC, XBSeq, FEELnc), use the `human_cgalcode_tot.gtf` file in Ensembl results and convert in `.bed` format with `perl ../../dependencies/gtf2bed.pl ../output_human_with_ensembl96_data_in_ref_in_total_gene_140621/human_cgalcode_tot.gtf > human_cgalcode_tot.bed` (example). In the next step, use the `bedtools intersect` command : `bedtools intersect -a human_cgalcode_tot.bed -b ../../examples/reference_gtf/human/ucsc/UCSC_Human_all_mrna_without_utr.bed -wa -wb -f 1 -r -split > result_bedtools_intersect_test.out` (example).


## bed_without_utr.py

This script delete the UTR information of the file and return a file without this utr information. Example:
`python3 bed_without_utr.py -i file`

### arguments: 
- "`-i`", "`--input`", Give a file with utr information to obtain a file without utr.



## extract_specific_junctions.py (used in run_multi_junction_alignment.sh)

This script extract the specific exon junctions of prediction (in GTF format) from a reference GTF file. Example:
`python3 extract_specific_junctions.py -p ${CGALCODE_PREDICTION_GTF} -r ${ENSEMBL_GTF} -d ${DELTA_EXTENSION}`

### arguments:
- "`-p`", "`--pred`", Give an exon gtf file with exon from prediction to search specific exon from prediction.
- "`-r`", "`--ref`", Give an exon gtf file with exon from reference to search specific exon from prediction.
- "`-d`", "`--delta`", Give a delta score for the exon extension, default = 25.



## extract_specific_junctions_from_FEELnc.py

This script extract the specific exon junctions of prediction (in GTF format) from a reference GTF file and in consider a experimental file. Example, we search the specific exon junction of our prediction and we check if this junction are known in Ensembl release 90 like reference and we check that FEELnc has not this junctions. Example:
`python3 extract_specific_junctions_from_FEELnc.py -p dog_cgalcode_tot.gtf -r dog_ensembl_tot.gtf -e canfam3_cons_annot_TritouPuppetMasterChief_23-03-2016_lncClasse_geneBiotype_withEnsCds_withSilicoCds_withGoodGeneSource.gtf -d 65`

### arguments:
- "`-p`", "`--pred`", Give an exon gtf file with exon from prediction to search specific exon from prediction.
- "`-r`", "`--ref`", Give an exon gtf file with exon from reference to search specific exon from prediction.
- "`-e`", "`--exp`", Give an exon gtf file with exon from an experimental procedure to search specific exon from prediction.
- "`-d`", "`--delta`", Give a delta score for the exon extension, default = 65.



## extract_specific_junctions_from_gff-exp_data.py

This script extract the specific exon junctions of prediction (in GTF format) from a reference GTF and in consider a experimental GFF file. Example, we search the specific exon junction of our prediction and we check if this junction are known in Ensembl release 90 like reference and we check that XBSeq or UCSC has not this junctions. Example: 
`python3 extract_specific_junctions_from_gff_exp_data.py -p human_cgalcode_tot.gtf -r human_ensembl_tot.gtf -e XBSeq_test.gff -d 65`
`python3 extract_specific_junctions_from_gff_exp_data.py -p human_cgalcode_tot.gtf -r human_ensembl_tot.gtf -e UCSC_Human_all_mrna_without_utr.gff -d 65`

### arguments:
- "`-p`", "`--pred`", Give an exon gtf file with exon from prediction to search specific exon from prediction.
- "`-r`", "`--ref`", Give an exon gtf file with exon from reference to search specific exon from prediction.
- "`-e`", "`--exp`", Give an exon gtf file with exon from an experimental procedure to search specific exon from prediction.
- "`-d`", "`--delta`", Give a delta score for the exon extension, default = 65.



## exon_junction_validator.py (used in run_multi_junction_alignment.sh)

This script extract the set of reads correctly aligned. Example: 
`python3 exon_junction_validator.py -a output_BedToolsIntersect_specific_junction_$filename.bed -j ../specific_junction_cga.bed`

### arguments:
- "`-a`", "`--align`", Give the file with the output alignment obtain with BedTools intersect for example, in BED12 format.
- "`-j`", "`--junction`", Give the file with the exon-exon junctions in BED12 format.



## nb_reads_mapped_graph.py (used in run_multi_junction_alignment.sh)

This script creates barplot to illustrate the read number aligned for each specific junction. Example:
`python3 nb_reads_mapped_graph.py -f nb_reads_mapped_tissue.csv -n 82`

### arguments:
- "`-f`", "`--file`", Give the file with the number of reads aligned in exon-exon junctions.
- "`-n`", "`--nb`", Give a number of juncions for each figure.



## run_multi_junction_alignment.sh

This script is a pipeline to realise the extraction of specific junction and their validation with bam files. Example:
`bash run_multi_junction_alignment.sh dog 65`

### arguments:
- `$1`: the species name (e.g. dog, human, mouse).
- `$2`: the delta extension (e.g. 65).
- other: see the content.



## run_multi_junction_alignment_dog_after_feelnc.sh

This script is a derived from *run_multi_junction_alignment.sh*. It uses a file obtain from *extract_specific_junctions_from_FEELnc.py* and its utilisation is the same like *run_multi_junction_alignment.sh*



## run_multi_junction_alignment_to_test_specificity_not_in_exp_data_in_*species*.sh

This script is a derived from *run_multi_junction_alignment.sh*. It uses a file obtain from *extract_specific_junctions_from_FEELnc.py* and its utilisation is the same like *run_multi_junction_alignment.sh*



## final_confortation_junction_without_annotation_validations.py

This script creates files whith all results obtains with the validation with ne annotation and merge with the specific junction research results. The final file was use to creates a diagram with `nb_reads_mapped_graph.py` script. Results are stored in `results_validation`. No options are implemented (See inside the script). Example: `python3 ../bin/final_confortation_junction_without_annotation_validations.py`.



## nb_tr_predicted_already_found.py

This script gives the number of transcripts predicted but already found in another species. Example:
`python3 nb_tr_predicted_already_found.py -s1 ../examples/set_genes/set_253_human.txt -s2 ../examples/set_genes/set_253_mouse.txt -s3 ../examples/set_genes/set_253_dog.txt `



## extract_read_alignment_class.r

This script split the file "nb_reads_mapped_tissue.csv" into two files: "null_list.txt" and "nonull_list.txt" (and more if necessary). The first contains all junctions with no read alignment and the second contains all junctions with read alignment. To execute this script, open RStudio and load the set working directory and execute the script.



## visual_gene_and_transcript_model_with_utr.py

This script extract gtf information and show a gene and transcript model with UTR (if gtf contains this information). Example of execution: `python3 visual_gene_and_transcript_model_with_utr.py -f /home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/mouse/ensembl/Ensembl102/Mus_musculus.GRCm38.102.gtf -g ENSMUSG00000070780`

### arguments:
- `-g`: the Ensembl identifier of a gene tested
- `-f`: the gtf of the species tested


# archive and not_conserved contain others scripts not conserved.
