# Gtf_validator_of_cgalcode_prediction

A software that converts the output from cgalcode to gtf and validates the predictions with new annotations or to conforts validation by a research of specific new junctions and validates that by a read alignment.

## How to run ?

### Transform a CG-alcode output in GTF and predicted transcripts validation part with new annotation

Create a new repository in the `predicted_transcripts_validation_with_new_annotations` like `test` to create a work repository. Use the command line `$ python3 ../../bin/gtf_CGalcode_generator_and_validator_predictions.py` in the work repository to run the software with the exemple data (reference gtf use correspond to the gtf found in [Ensembl](ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz). The predicted transcripts were obtain with CG-Alcode programm (development by S. Blanquart and J.S. Varré [Blanquart S., 2006]) and complete by N. Guillaudeux. All data generated for our work are available in `predicted_transcripts_validation_with_new_annotations` directory with the `output_species_with_database&version_data_in_ref` in work repositories name. In this Readme, we presented a example based on a set of 135 genes analysis (all CDS conserved in bijective relation between human, mouse and dog species): `output_human_with_ensembl98_data_in_ref`. In results, there are the 253 genes set (all proteome conserved between human, mouse and dog species) and the total set of genes analysis (2,167 genes).

With this example, 19 predicted transcripts are validated on the 98 predicted transcripts by CG-Alcode software (see `result_bedtols_intersect.out` file and the `Sample output` section).

If you use this programm, it needs these parameters:
 - `-gtf2bed`, `--gtf2bedconvertion`: give the Gtf to Bed script converter. By default, use the perl script.
 - `-rep`, `--repository`: enter the path repository where are located genes and transcripts files used by CG-Alcode programm.
 - `-c`, `--correspondance`: give the file with the description relationship between genes and their transcripts.
 - `-set`, `--setgene`: give a list of interested genes.
 - `-ref`, `--reference`: give a reference gtf file to validate prediction realize with CG-Alcode programm.

The output results are stored in `result_bedtools_intersect_test.out` within the work repository.

Example with the Ensembl release 103: `python3 ../../bin/gtf_CGalcode_generator_and_validator_predictions.py -rep /home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ -c /home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_hs_withprediction.csv -set ../../examples/set_genes/set_253_human.txt -ref ../../examples/reference_gtf/human/ensembl/Ensembl103/Homo_sapiens.GRCh38.103.gtf`

Information to the script execution with another data not Ensembl: This script is executable with Ensembl data, for another data (UCSC, XBSeq, FEELnc), use the `human_cgalcode_tot.gtf` file in Ensembl results and convert in `.bed` format with `perl ../../dependencies/gtf2bed.pl ../output_human_with_ensembl96_data_in_ref_in_total_gene_140621/human_cgalcode_tot.gtf > human_cgalcode_tot.bed` (example). In the next step, use the `bedtools intersect` command : `bedtools intersect -a human_cgalcode_tot.bed -b ../../examples/reference_gtf/human/ucsc/UCSC_Human_all_mrna_without_utr.bed -wa -wb -f 1 -r -split > result_bedtools_intersect_test.out` (example).


### Validation of specific junction part with reads files in BAM format

Use the command line `$ bash ../../../bin/run_multi_junction_alignment.sh species delta bam_repository data_repository` to run the validation part with the bam files (the localisation of bam files is necessary in the script file). Chose the parameters `species` (e.g.: "human", "mouse", "dog"), `delta` of extension read length (e.g.: 65), `bam_repository` the repository with bam files and `data_repository` with data informations (species prediction information in cgalcode gtf file, species reference information same that ensembl in gtf format (cf. script for details)). Example: 
    `gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/dog_junction_exons$ bash ../../bin/run_multi_junction_alignment.sh dog 65 /.../IGDR_data/bam /.../gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/output_dog_with_ensembl98_data_in_ref`

### Final results

The script `final_confortation_junction_without_annotation_validations.py` finalise all results obtains in the two approach in combination. All this results are stored in `results_validation`. The three pictures are obtain with the `nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_dog_without_alignment.csv` script. Here are the command lines used:
- `python3 ../bin/final_confortation_junction_without_annotation_validations.py`
- `python3 ../bin/nb_reads_mapped_graph.py -f nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_dog_with_alignment.csv -n 11`
- `python3 ../bin/nb_reads_mapped_graph.py -f nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_mouse_with_alignment.csv -n 32`
- `python3 ../bin/nb_reads_mapped_graph.py -f nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_human_with_alignment.csv -n 30`


## Sample output

### Transform a CG-alcode output in GTF and validation part with new annotation

    `gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/test$ python3 ../../bin/    gtf_CGalcode_generator_and_validator_predictions.py`
    ########################################
    ########## PROGRAMM BEGINNING ##########
    
    >>> Directory containing genes and transcripts files: True
    >>> File containing genes and transcripts correspondance: True
    >>> File containing the list of genes of interest: True
    >>> GTF reference file to check predictions: True
    
    >>> Creation: a gtf file with all information (Ensembl, prediction)
    >>> Creation: gtf Ensembl files for each gene
    >>> Creation: a gtf file with all Ensembl information
    Number of empty gene without Ensembl ref: 0
    >>> Creation: gtf prediction files for each gene
    >>> Creation: a gtf file with all prediction information
    Number of empty gene without prediction: 70
    >>> Creation: gtf reference files for each gene
    Number of empty gene without reference: 0
    
    >>> Conservation of genes of interest
    Number of genes with predicted transcripts conserved: 135/135
    Number of ref genes conserved: 135/60622
    Number of ens genes conserved: 135/135
    Number of ref transcripts conserved: 1373
    Number of ens transcripts used by CGAlcode: 364
    Number of pred transcripts to check: 98
    
    >>> Convertion GTF to BED12
    >>> BedTools intersect analysis
    Number of validate transcript:
    19
    
    ############ PROGRAMM ENDING ############
    #########################################

### Validation part with read files in BAM format

    `gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/dog_junction_exons$ bash ../../bin/run_multi_junction_alignment.sh dog 65 /run/media/.../Dyliss_NG/IGDR_data/bam /home/.../gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/output_dog_with_ensembl98_data_in_ref`

    The species choice is write: dog
    All tissus are stored in: /run/media/niguilla/Dyliss_NG/IGDR_data/bam
    Corresponding data are stored in : /home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/output_dog_with_ensembl98_data_in_ref
    Delta chosen for the junction extension is: 65
    ----------------------
    >>> START PROGRAMM <<<
    ----------------------
    -> 2020-01-20 09:16:46
    >>> WORK WITH ADRENALGLAND_BEMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:16:46
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:16:46
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:21:07
    >>> WORK WITH BLOOD_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:21:08
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:21:08
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:24:08
    >>> WORK WITH BRAIN_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:24:08
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:24:08
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:27:11
    >>> WORK WITH CEREBELLUM_BESH
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:27:12
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:27:12
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:31:02
    >>> WORK WITH CEREBELLUM_GSMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:31:02
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:31:03
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:34:24
    >>> WORK WITH CORTEX_BESH
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:34:25
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:34:25
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:37:43
    >>> WORK WITH GUTCOLON_BEMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:37:43
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:37:44
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:41:49
    >>> WORK WITH H01_GRET
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:41:50
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:41:50
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:45:51
    >>> WORK WITH H01_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:45:52
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:45:52
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:49:35
    >>> WORK WITH H06_PODL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:49:35
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:49:35
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:53:54
    >>> WORK WITH HAIRFOLLICULE_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:53:55
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:53:55
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 09:58:14
    >>> WORK WITH HEART_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 09:58:15
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 09:58:15
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:00:16
    >>> WORK WITH JEJUNUM_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:00:17
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:00:17
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:04:20
    >>> WORK WITH KERATINOCYTE_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:04:20
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:04:20
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:08:19
    >>> WORK WITH KIDNEY_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:08:19
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:08:19
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:10:23
    >>> WORK WITH LIVER_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:10:23
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:10:23
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:12:26
    >>> WORK WITH LUNG_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:12:26
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:12:26
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:14:09
    >>> WORK WITH MAMMARYGLAND_GSMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:14:09
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:14:10
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:17:38
    >>> WORK WITH MUSCLE_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:17:39
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:17:39
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:19:30
    >>> WORK WITH NOSE01_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:19:31
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:19:31
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:24:15
    >>> WORK WITH NOSE02_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:24:16
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:24:16
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:29:46
    >>> WORK WITH NOSE03_LABR
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:29:46
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:29:47
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:36:27
    >>> WORK WITH OLFBULB_GSMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:36:27
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:36:28
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:40:18
    >>> WORK WITH OVARY_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:40:18
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:40:18
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:42:27
    >>> WORK WITH PANCREAS_BESH
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:42:28
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:42:28
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:47:34
    >>> WORK WITH RETINA_BOCL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:47:34
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:47:34
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:51:31
    >>> WORK WITH SKIN_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:51:31
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:51:31
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 10:55:51
    >>> WORK WITH SKIN_GSMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 10:55:51
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 10:55:52
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 11:00:11
    >>> WORK WITH SPINALCORD_GSMD
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 11:00:12
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 11:00:12
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 11:03:44
    >>> WORK WITH SPLEEN_BESH
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 11:03:44
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 11:03:44
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 11:08:00
    >>> WORK WITH TESTIS_BEGL
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 11:08:00
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 11:08:00
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 11:10:15
    >>> WORK WITH THYMUS_SLKI
    # STEP1: EXTRACTION OF SPECIFIC EXON COUPLES
    -> 2020-01-20 11:10:15
    # STEP2: ALIGNMENT OF BAM FILES WITH BEDTOOLS INTERSECT
    -> 2020-01-20 11:10:15
    # STEP3: EXTRACTION OF READS CORRECTLY ALIGNED
    -> 2020-01-20 11:14:16
    # STEP4: COPY FILES WITH NUMBER OF READS MAPPED BY TISSUE IN A SAME DIRECTORY
    -> 2020-01-20 11:14:17
    cp: 'output_run_ADRENALGLAND_BEMD.txt' et '../analysed/output_run_ADRENALGLAND_BEMD.txt' identifient le même fichier
    # STEP5: CREATE A FILE WITH THE READS NUMBER MAPPED WITH SPECIFIC EXON COUPLES
    -> 2020-01-20 11:14:17
    STEP6: CREATE A GRAPHIC PDF FILE
    -> 2020-01-20 11:14:17
    --------------------
    >>> End programm <<<
    --------------------
    >>> 2020-01-20 11:14:27

### Merge results

    niguilla@raptor:~/Documents/Software/gtf_validator_of_cgalcode_prediction/results_validation$ python3 ../bin/final_confortation_junction_without_annotation_validations.py 
    ########################################
    ########## PROGRAMM BEGINNING ##########
    ########## DOG ANALYSIS ##########
    >>>Total tr validated in dog: 185
    tr with reads alignments: 11
    tr without reads alignments: 16
    tr with reads alignments: 11
    tr without reads alignments: 16
    tr with and without reads alignments: 0
    ########## HUMAN ANALYSIS ##########
    >>>Total tr validated in human: 26
    tr with reads alignments: 27
    tr without reads alignments: 18
    tr with reads alignments: 26
    tr without reads alignments: 17
    tr with and without reads alignments: 1
    ########## MOUSE ANALYSIS ##########
    >>>Total tr validated in mouse: 32
    tr with reads alignments: 25
    tr without reads alignments: 37
    tr with reads alignments: 22
    tr without reads alignments: 34
    tr with and without reads alignments: 3
    ########### PROGRAMM ENDING ############

## Data tested:

We have used different references data to validate predicted transcript ('-ref'):

### New annotations

We had stored all GTF files used in `examples/reference_gtf/species/database/database&Version/file.gtf`.

- Ensembl data from the release 96 for [human](ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens), [mouse](ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus) and [dog](ftp://ftp.ensembl.org/pub/release-96/gtf/canis_familiaris/) stored in 06/2019.

- Ensembl data from the release 98 for [human](ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens), [mouse](ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus) and [dog](ftp://ftp.ensembl.org/pub/release-98/gtf/canis_familiaris/) stored in 10/2019.

- Ensembl data from the release 102 for [human](ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens), [mouse](ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus) and [dog](ftp://ftp.ensembl.org/pub/release-102/gtf/canis_lupus_familiaris/) stored in 10/2019.

- Ensembl data from the release 103 for [human](ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens) and [dog](ftp://ftp.ensembl.org/pub/release-103/gtf/canis_familiaris/) stored in 10/2019.

- XBSeq data availale [here](https://github.com/Liuy12/XBSeq_files/) for human and mouse and stored in 06/2019.

- UCSC data available for the three species and stored in 06/2019.

- IGDR data (from FEELnc study) available for the dog species and stored in 06/2019.

### Research of new specific informations

- Wuchet V., 2017 - FEELnc: a tool for long non-coding RNA annotation and its application to the dog transcriptome.
- Söllner J. F., 2017 -  An RNA-Seq atlas of gene expression in mouse and rat normal tissues.
- Wang D., 2019 - A deep proteome and transcriptome abundance atlas of 29 healthy human tissues.
