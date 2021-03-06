#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("reshape2")
library("ggplot2")
library("cowplot")
library("reshape2")

####################################################################################
### HUMAN DATA

# read file
data_hs = read.table("nb_reads_mapped_tissue_hs.csv", header = T, sep='\t')
summary(data_hs)

# sum by rows
data_hs$sum = apply(data_hs[,2:ncol(data_hs)],1,sum)
summary(data_hs)

# write files
#write.table(data[data$sum==0,],"null_sum.csv", sep='\t')
write.table(data_hs[data_hs$sum>0,],"nonull_sum_hs.csv", sep='\t')

#data_null = read.table("null_sum.csv", h=T, sep='\t')
data_hs_not_null = read.table("nonull_sum_hs.csv", h=T, sep='\t')
summary(data_hs_not_null)

plot(density(data_hs_not_null$sum))
hist(data_hs_not_null$sum, freq = T)
?density
p <- ggplot(data_hs_not_null, aes(x=sum)) + 
  geom_density()

p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez l'humain") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")


boxplot(data_hs_not_null$sum)
hist(data_hs_not_null$sum, nclass=10000)
#hist(data_hs_not_null$sum, br = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300))


data_hs  = read.table("nonull_sum_hs.csv", sep='\t', header=TRUE)
summary(data_hs)
tableau_hs <- table(data_hs$sum)
#hist(tableau)
print(tableau_hs)
addmargins(tableau_hs)
prop.table(tableau_hs) 
write.table(tableau_hs,"tableau_hs_sum.csv", sep='\t')
plot(tableau_hs)

####################################################################################
### MOUSE DATA

# read file
data_mm = read.table("nb_reads_mapped_tissue_mm.csv", header = T, sep='\t')
summary(data_mm)

# sum by rows
data_mm$sum = apply(data_mm[,2:ncol(data_mm)],1,sum)
summary(data_mm)

# write files
write.table(data_mm[data_mm$sum>0,],"nonull_sum_mm.csv", sep='\t')

data_mm_not_null = read.table("nonull_sum_mm.csv", h=T, sep='\t')
summary(data_mm_not_null)

p <- ggplot(data_mm_not_null, aes(x=sum)) + 
  geom_density()
p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez la souris") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")

hist(data_mm_not_null$sum, nclass=922)
hist(data_mm_not_null$sum, br = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300))

data_mm  = read.table("nonull_sum_mm.csv", sep='\t', header=TRUE)
summary(data_mm)
tableau_mm <- table(data_mm$sum)
#hist(tableau)
print(tableau_mm)
addmargins(tableau_mm)
prop.table(tableau_mm) 

plot(tableau_mm)


####################################################################################
### DOG DATA

# read file
data_clf = read.table("nb_reads_mapped_tissue_clf.csv", header = T, sep='\t')
summary(data_clf)

# sum by rows
data_clf$sum = apply(data_clf[,2:ncol(data_clf)],1,sum)
summary(data_clf)

# write files
write.table(data_clf[data_clf$sum>0,],"nonull_sum_clf.csv", sep='\t')

data_clf_not_null = read.table("nonull_sum_clf.csv", h=T, sep='\t')
summary(data_clf_not_null)

p <- ggplot(data_clf_not_null, aes(x=sum)) + 
  geom_density()
p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez le chien") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")

hist(data_clf_not_null$sum)


data_clf  = read.table("nonull_sum_clf.csv", sep='\t', header=TRUE)
summary(data_clf)
tableau_clf <- table(data_clf$sum)
#hist(tableau)
print(tableau_clf)
addmargins(tableau_clf)
prop.table(tableau_clf) 

plot(tableau_clf)


plot_grid(tableau_hs, tableau_mm, tableau_clf, labels=c("hs", "mm", "clf"), ncol = 3, nrow = 1)

ggplot()




####################################################################################
### ALL DATA #File manually creates

data_not_null = read.table("nonull_sum_all.csv", h=T, sep='\t')
summary(data_not_null)

p <- ggplot(data_not_null, aes(x=species, y=sum, fill=species)) + 
  geom_boxplot()
p+scale_fill_manual(values=c("pink", "#56B4E9", "grey"))
p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées") +
  xlab("Espèces") + 
  ylab("Nombre de reads")


### BOX PLOT

#layout(matrix(1:3,1,3))

data_not_null = read.table("nonull_sum_all_v3_hs.csv", h=T, sep='\t')
summary(data_not_null)

p1 <- ggplot(data_not_null, aes(x=species, y=sum)) + 
  geom_boxplot() + 
  geom_boxplot(fill='#56B4E9', color="black") + 
  ggtitle("HumanChez l'humain") +
  xlab("") + 
  ylab("")

data_not_null = read.table("nonull_sum_all_v3_mm.csv", h=T, sep='\t')
summary(data_not_null)

p2 <- ggplot(data_not_null, aes(x=species, y=sum)) + 
  geom_boxplot() + 
  geom_boxplot(fill='grey', color="black") + 
  ggtitle("Chez la souris") +
  xlab("") + 
  ylab("")

data_not_null = read.table("nonull_sum_all_v3_clf.csv", h=T, sep='\t')
summary(data_not_null)

p3 <- ggplot(data_not_null, aes(x=species, y=sum)) + 
  geom_boxplot() + 
  geom_boxplot(fill='pink', color="black") + 
  ggtitle("Chez le chien") +
  xlab("") + 
  ylab("")

plot_grid(p1, p2, p3, labels=c("", "", ""), ncol = 3, nrow = 1)

#p+scale_fill_manual(values=c("pink", "#56B4E9", "grey"))
#geom_density()

data  = read.table("nonull_sum_all.csv", sep='\t', header=TRUE)
summary(data)
tableau <- table(data$sum, data$species)
#hist(tableau)
print(tableau)
addmargins(tableau)
prop.table(tableau) 
tableau2 <- addmargins(prop.table(tableau,2) )

plot(tableau)

write.table(tableau,"table_sum_all.csv", sep='\t')
data_lecture_table  = read.table("table_sum_all.csv", sep='\t', header=TRUE)
data_lecture_table
data_lecture_table$read = as.numeric(rownames(data_lecture_table))
data_lecture_table
?ggplot

dataMelted <- reshape2::melt(data_lecture_table, id.var='read') 
dataMelted
ggplot(dataMelted, aes(x = read, y=value, col=variable)) + geom_line()



####################################################################################
### REFERENCE HUMAN DATA (other "Session > Set Working Directory > Choose Directory... > reference" !!!)

# read file
data_ref_hs = read.table("nb_reads_mapped_tissue_hs_ref.csv", header = T, sep='\t')
summary(data_ref_hs)

# sum by rows
data_ref_hs$sum = apply(data_ref_hs[,2:ncol(data_ref_hs)],1,sum)
summary(data_ref_hs)

# write files
#write.table(data[data$sum==0,],"null_sum.csv", sep='\t')
write.table(data_ref_hs[data_ref_hs$sum>0,],"nonull_sum_hs.csv", sep='\t')

#data_null = read.table("null_sum.csv", h=T, sep='\t')
data_ref_hs_not_null = read.table("nonull_sum_hs.csv", h=T, sep='\t')
summary(data_ref_hs_not_null)

p <- ggplot(data_ref_hs_not_null, aes(x=sum)) + 
  geom_density()

p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez l'humain") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")





# read file
data = read.table("nb_reads_mapped_tissue_mm_ref.csv", header = T, sep='\t')
summary(data)

# sum by rows
data$sum = apply(data[,2:ncol(data)],1,sum)
summary(data)

# write files
#write.table(data[data$sum==0,],"null_sum.csv", sep='\t')
write.table(data[data$sum>0,],"nonull_sum_mm.csv", sep='\t')

#data_null = read.table("null_sum.csv", h=T, sep='\t')
data_not_null = read.table("nonull_sum_mm.csv", h=T, sep='\t')
summary(data_not_null)

p <- ggplot(data_not_null, aes(x=sum)) + 
  geom_density()

p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez l'humain") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")





# read file
data = read.table("nb_reads_mapped_tissue_clf_ref.csv", header = T, sep='\t')
summary(data)

# sum by rows
data$sum = apply(data[,2:ncol(data)],1,sum)
summary(data)

# write files
#write.table(data[data$sum==0,],"null_sum.csv", sep='\t')
write.table(data[data$sum>0,],"nonull_sum_clf.csv", sep='\t')

#data_null = read.table("null_sum.csv", h=T, sep='\t')
data_not_null = read.table("nonull_sum_clf.csv", h=T, sep='\t')
summary(data_not_null)

p <- ggplot(data_not_null, aes(x=sum)) + 
  geom_density()

p + ggtitle("Distribution du nombre de reads alignés sur les jonctions testées chez l'humain") +
  xlab("Nombre de reads alignés sur les jonctions testées") + 
  ylab("Fréquence")



data  = read.table("reference/nonull_sum_all_ref.csv", sep='\t', header=TRUE)
summary(data)
tableau <- table(data$sum, data$species)
#hist(tableau)
print(tableau)
addmargins(tableau)
prop.table(tableau) 
tableau2 <- addmargins(prop.table(tableau,2) )

write.table(tableau,"table_sum_all_ref.csv", sep='\t')




data_analysis_spec = read.table("table_sum_all.csv", sep = '\t')
plot(data_analysis_spec$human)
hist(data_mm_not_null$sum, nclass=922)
hist(data_mm_not_null$sum, br = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300))

