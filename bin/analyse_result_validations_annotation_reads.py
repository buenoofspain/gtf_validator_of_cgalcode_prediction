#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:24:41 2020

@author: niguilla
"""

def createListDetailsValidationByAnnotation(directoryPath, species, list_species_tr_validated):
    """
    FUNCTION DESCRIPTION:
        This function creates a list of predicted transcripts where a 
        correspondance was found in a new annotation.
    
    Input:
        - directoryPath: the path of the repository where is stored the file
                       "result_bedtools_intersect.out"
        - species: the species tested : human, mouse or dog
        - list_species_tr_validated: the list to store transcripts with a 
                                     correspondance in a new annotation
    
    Return:
        - list_species_tr_validated: the list to store transcripts with a 
                                     correspondance in a new annotation
    """
    
#    list_species_tr_validated = list()
    
    with open(directoryPath+"/result_bedtools_intersect_test.out", 'r') as path_d_f:
        lecture = path_d_f.readlines()
        for line in lecture:
            line = line.replace('\n', '').split('\t')
            try:
                tr = line[3]
                if not tr in list_species_tr_validated:
                    list_species_tr_validated.append(tr)
            except IndexError:
                pass

#    print("Number of predicted transcripts view in a new annotation:", len(list_species_tr_validated))

    return list_species_tr_validated
        

def createFileWithValidationByAnnotation(listTrValidatedByAnnotation, species):
    """
    FUNCTION DESCRIPTION:
        This functuion creates a file with the list of transcripts where a
        correspondance was found in a new annotation.
    
    Input:
        - listTrValidatedByAnnotation: the list to store transcripts with a 
                                       correspondance in a new annotation
        - species: the species tested : human, mouse or dog
    """
    with open("list_tr_validated_in_"+species+".txt", 'w') as list_tr_o:
        for tr_in_list in listTrValidatedByAnnotation:
            list_tr_o.write(tr_in_list+'\n')


def reduceJunctionReadsAlignmentResults(nbReadsMappedFileResult, species, listTrValidatedByNewAnnotation):
    """
    FUNCTION DESCRIPTION:
        This function reduced the file with the number of reads mapped on 
        specific junctiono if the junction does not belong to a predicted 
        transcript validated by a new annotation.
        
    Input:
        - nbReadsMappedFileResult: the file with all results of reads mapping 
                                   on junction tested
        - species: the species tested : human, mouse or dog
        - listTrValidatedByNewAnnotation: the list to store transcripts with a 
                                          correspondance in a new annotation
    """
    with open(nbReadsMappedFileResult, 'r') as mappedResult_i:
        with open("nb_reads_mapped_tissue_reduces_without_tr_validated_by_annotation_in_"+species+".csv", 'w') as mappedResult_o:
            lecture = mappedResult_i.readlines()
            for line in lecture:
                if line.startswith("Junctions"):
                    mappedResult_o.write(line)
                else:
                    line_split = line.replace('\n', '').split('\t')#
                    junc_tested = line_split[0]
                    tr_tested = junc_tested.split('JUNC')[0]
                    if not tr_tested in listTrValidatedByNewAnnotation:
                        mappedResult_o.write(line)
    
    
def splitFileResultWithMappedReadReduced(species):
    """
    """
    
    tr_with_reads_alignment = list()
    tr_without_reads_alignment = list()
    tr_with_and_without_alignment = list()
    
    with open("nb_reads_mapped_tissue_reduces_without_tr_validated_by_annotation_in_"+species+".csv", 'r') as f_i:
        with open("nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_"+species+"_with_alignment.csv", 'w') as f_o_yes:
            with open("nb_reads_mapped_tissue_reduces_without_tr_validates_by_annotation_in_"+species+"_without_alignment.csv", 'w') as f_o_no:
                lecture = f_i.readlines()
                for line in lecture:
                    if line.startswith("Junction"):
                        f_o_yes.write(line)
                        f_o_no.write(line)
                        
                    else:
                        line_split = line.replace('\n', '').split('\t')
                        total_reads_aligned = 0
                        for tissu in range(1, len(line_split)):
                            total_reads_aligned += int(line_split[tissu])
                        if total_reads_aligned == 0:
                            f_o_no.write(line)
                            if not line_split[0].split("JUNC")[0] in tr_without_reads_alignment:
                                tr_without_reads_alignment.append(line_split[0].split("JUNC")[0])
                        else:
                            f_o_yes.write(line)
                            if not line_split[0].split("JUNC")[0] in tr_with_reads_alignment:
                                tr_with_reads_alignment.append(line_split[0].split("JUNC")[0])
                        
    print("tr with reads alignments (possible without):", len(tr_with_reads_alignment))
    print("tr without reads alignments (possible with):", len(tr_without_reads_alignment))
    
    if "CGACAFTENST00000618462" in tr_with_reads_alignment:
        print(True)
    if "CGACAFTENST00000618462" in tr_without_reads_alignment:
        print(True)
    

    for tr in tr_with_reads_alignment:
        if tr in tr_without_reads_alignment:
#            print("TEST1", tr)
            tr_with_and_without_alignment.append(tr)
            tr_with_reads_alignment.remove(tr)
            tr_without_reads_alignment.remove(tr)
            
    for tr in tr_without_reads_alignment:
        if tr in tr_with_reads_alignment:
#            print("TEST2", tr)
            tr_with_and_without_alignment.append(tr)
            tr_with_reads_alignment.remove(tr)
            tr_without_reads_alignment.remove(tr)
            
    print("tr with reads alignments:", len(tr_with_reads_alignment))
    print("tr without reads alignments:", len(tr_without_reads_alignment))
    print("tr with and without reads alignments:", len(tr_with_and_without_alignment))

if __name__ == "__main__":
    
    print("########################################")
    print("########## PROGRAMM BEGINNING ##########")
    
    path = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/"
    
    print("########## DOG ANALYSIS ##########")
    sp = "dog"
    list_species_tr_validated = list()
    
    path_dir = path+"/dog_2167_180621/output_dog_with_ensembl96_data_in_ref_in_total_gene_150621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"/dog_2167_180621/output_dog_with_ensembl98_data_in_ref_in_total_gene_150621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"/dog_2167_180621/output_dog_with_ensembl102_data_in_ref_in_total_gene_150621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"/dog_2167_180621/output_dog_with_ensembl103_data_in_ref_in_total_gene_150621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"/dog_2167_180621/output_dog_with_feelnc_data_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"/dog_2167_180621/output_dog_with_ucsc_data_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    
    print(">>>Total tr validated in "+sp+":", len(list_species_tr_validated))
    
    createFileWithValidationByAnnotation(list_species_tr_validated, sp)
    
#    mappedResults_pathFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/dog_junction_exons/junction_validation_dog_2020-01-20/analysed/nb_reads_mapped_tissue.csv"
#    mappedResults_pathFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/dog_junction_exons/junction_validation_dog_2019-12-19_final/analysed/nb_reads_mapped_tissue.csv"
    mappedResults_pathFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/dog_junction_exons/junction_validation_dog_2021-06-25/analysed/nb_reads_mapped_tissue.csv"
    reduceJunctionReadsAlignmentResults(mappedResults_pathFile, sp, list_species_tr_validated)
    splitFileResultWithMappedReadReduced(sp)
    
    print("########## HUMAN ANALYSIS ##########")
    sp = "human"
    list_species_tr_validated = list()
    
#    path_dir = path+"output_human_with_ensembl_data_in_ref_new_gtf2bed"
#    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_ensembl96_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_ensembl98_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_ensembl102_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_ensembl103_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_ucsc_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"human_2167_genes_180621/output_human_with_xbseq_data_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    
    print(">>>Total tr validated in "+sp+":", len(list_species_tr_validated))
    
    createFileWithValidationByAnnotation(list_species_tr_validated, sp)
    
    mappedResults_pathFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/human_junction_exons/junction_validation_human_2021-06-25/analysed/nb_reads_mapped_tissue.csv"
    reduceJunctionReadsAlignmentResults(mappedResults_pathFile, sp, list_species_tr_validated)
    splitFileResultWithMappedReadReduced(sp)
    
    print("########## MOUSE ANALYSIS ##########")
    sp = "mouse"
    list_species_tr_validated = list()
    
    path_dir = path+"mouse_2167_180621/output_mouse_with_ensembl96_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"mouse_2167_180621/output_mouse_with_ensembl98_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"mouse_2167_180621/output_mouse_with_ensembl102_data_in_ref_in_total_gene_140621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"mouse_2167_180621/output_mouse_with_ucsc_data_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    path_dir = path+"mouse_2167_180621/output_mouse_with_xbseq_data_in_ref_in_total_gene_160621"
    list_species_tr_validated = createListDetailsValidationByAnnotation(path_dir, sp,list_species_tr_validated)
    
    
    print(">>>Total tr validated in "+sp+":", len(list_species_tr_validated))
    
    createFileWithValidationByAnnotation(list_species_tr_validated, sp)
    
    mappedResults_pathFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/mouse_junction_exons/junction_validation_mouse_2021-06-25/analysed/nb_reads_mapped_tissue.csv"
    reduceJunctionReadsAlignmentResults(mappedResults_pathFile, sp, list_species_tr_validated)
    splitFileResultWithMappedReadReduced(sp)
    
    
    print("########### PROGRAMM ENDING ############")
    print("########################################")
