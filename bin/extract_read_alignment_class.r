# read file
data = read.table("nb_reads_mapped_tissue.csv", header = T, sep='\t')
summary(data)

# sum by rows
data$sum = apply(data[,2:ncol(data)],1,sum)
summary(data)

# write files
write.table(data[data$sum==0,],"null_sum.csv", sep='\t')
write.table(data[data$sum>0,],"nonull_sum.csv", sep='\t')

data_null = read.table("null_sum.csv", h=T, sep='\t')
data_not_null = read.table("nonull_sum.csv", h=T, sep='\t')
summary(data_not_null)

#write.table(data_null$Junctions, "null_list.txt", row.names = F, col.names = F)
write(as.vector(data_null$Junctions), "null_list.txt")
write(as.vector(data_not_null$Junctions), "nonull_list.txt")
