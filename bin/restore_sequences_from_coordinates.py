#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:06:14 2020

@author: niguilla
"""

def gene_string(geneFile):
    """
    Function description:
        This function uses a gene file where all sequences are split with a 
        certain number of caracters and return all caracter in one string 
        sequence.
    
    Input:
        - geneFile: a gene file in fasta format.
    
    Return:
        joinSequences: a string sequence of the gene
    """
    
    geneDescriptionList = list()
    
    with open(geneFile, 'r') as g_o:
        lecture = g_o.readlines()
        
        joinSequences = ""
        for line in lecture:
            if line.startswith(">"):
                line = line.replace('\n', '').split(" ")
                geneId = line[1]
                species = line[2]
                annotation = line[4]
                startGenePosition = line[5]
                endGenePosition = line[6]
                strand = line[7]
                geneDescriptionList.append(geneId)
                geneDescriptionList.append(species)
                geneDescriptionList.append(startGenePosition)
                geneDescriptionList.append(endGenePosition)
                geneDescriptionList.append(strand)
                
            else:
                joinSequences += line.replace('\n', '')
            
        
        if strand == '-1':
            joinSequencesReversed = "".join(reversed(joinSequences))
            geneDescriptionList.append(joinSequencesReversed)
            geneDescriptionList.append("REVERSED")
        else:
            geneDescriptionList.append(joinSequences)
            geneDescriptionList.append("NOT_REVERSED")
            
#        geneDescriptionList.append(joinSequences)
#            
    return geneDescriptionList
                

def restore_sequences_from_coordinates(geneInformation, transcriptFile):
    """
    Function description:
        This function uses coordinates to reconstruct nucleotidic sequence from
        the sequences research.
        
    Input:
        - geneInformation: The set of all information of the gene and sequence.
        - transcriptFile: a file with all transcript information in CGalcode format.
    """

    geneId = geneInformation[0]
    species = geneInformation[1]
    startGenePos = int(geneInformation[2])
    endGenePos = int(geneInformation[3])
    strand = int(geneInformation[4])
    geneSequence = geneInformation[5]
    reversedState = geneInformation[6]
    
    with open(transcriptFile, 'r') as tr_o:
        lecture = tr_o.readlines()
#        print(lecture)
#        print(lecture[0].replace('\n', ''))
        seq = ""
        for transcriptDescription in (lecture[4].replace('\n', '').replace("'", '').replace("], [", ";").replace('[', '').replace(']', '').split(';')):
            eventDescription = transcriptDescription.split(',')
            startEventPos = int(eventDescription[0])
            endEventPos = int(eventDescription[1])
            typeEvent = eventDescription[2].replace(" ", "")
            
            if typeEvent == "EXON":
#                print(strand)
                if strand == -1:
#                    print(startEventPos, endEventPos, typeEvent)
                    seqExon = (geneSequence[startEventPos-startGenePos:endEventPos-startGenePos+1])
                    seqExonReversed = "".join(reversed(seqExon))
#                    print(seqExonReversed)
                    seq += seqExonReversed
                else:
#                    print(startEventPos, endEventPos, typeEvent)
                    seqExon = (geneSequence[startEventPos-startGenePos:endEventPos-startGenePos+1])
#                    print(seqExon)
                    seq += seqExon
    
        print(lecture[0].replace('\n', ''))
        print(seq)


if __name__ == "__main__":
    
#    geneFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/ENSCAFG00000004626.fasta"
#    geneInfo = gene_string(geneFile)
#    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/CGACAFTENST00000252485.fasta"
#    restore_sequences_from_coordinates(geneInfo, transcriptFile)


    
    geneFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENSG00000181722.fasta"
    geneInfo = gene_string(geneFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000357258.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000462705.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000481632.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000471418.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000464560.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/ENST00000393785.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    
    
    geneFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/ENSMUSG00000022708.fasta"
    geneInfo = gene_string(geneFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/ENSMUST00000114691.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/ENSMUST00000114690.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)
    
    geneFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/ENSCAFG00000010823.fasta"
    geneInfo = gene_string(geneFile)
    
    transcriptFile = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/CGACAFTENST00000357258.fasta"
    restore_sequences_from_coordinates(geneInfo, transcriptFile)