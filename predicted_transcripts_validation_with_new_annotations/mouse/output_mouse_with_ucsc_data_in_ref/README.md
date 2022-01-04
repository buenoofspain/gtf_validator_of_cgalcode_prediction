Use the bed format to prediction file and the bed format of the reference file w
ith bedtools intersect:
> bedtools intersect -a file_pred.bed -b file_ref.bed -wa -wb -f 1 -r -split > r
esult_bedtools_intersect.out
