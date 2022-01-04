# Lecture du tableau de l'ensemble des blocs d'exons
data <- read.table("distri_exon_junction_len_all.csv", header = TRUE, sep = ",", dec=',')
head(data)
summary(data)

# Lecture du tableau de l'ensemble des blocs d'exons séparés
data_split <- read.table("distri_exon_junction_len_sep.csv", header = TRUE, sep = ",", dec=',')
head(data_split)
summary(data_split)

#####################################################################
# Répartition des blocs d'exons dans les couples d'adjacence totaux #
#####################################################################
main=paste("Répartition de la taille des", nrow(data), "blocs d'exons pour les", nrow(data)/2, "couples d'adjacence")
sub_main=paste("triés en fonction du brin + pour les", nrow(data)/2, "couples d'adjacence spécifiques")
xlab="Length of exon blocks"
ylab="Number of blocks"

graph_kn <- plot(table(data$All),
                 #col="grey",
                 #col = rainbow(4, start = 0.7, end = 0.9),
                 col = rgb((1:4)/4, 0, 0),
                 main=paste(main,'\n',sub_main),
                 #col = c("red", "blue"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,10),
                 xlim=c(min(data$All),max(data$All)),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn, y = table(data$tr_known), label = table(data$tr_known), pos = 3, cex = 0.9, col = "red")

##########################################################################
# Répartition des blocs d'exons dans les couples d'adjacence spécifiques #
##########################################################################

# Lecture du tableau de l'ensemble des 89 blocs d'exons séparés
data_split_89 <- read.table("distri_exon_89_junctions_len_sep.csv", header = TRUE, sep = ",", dec=',')
head(data_split_89)
summary(data_split_89)

#####################################
# Sur les blocs exoniques de gauche #
#####################################
main=paste("Répartition de la taille des exons situés à gauche de chacun des", nrow((data_split_89)))
sub_main=paste("couples d'adjacence triés en fonction du brin +")
xlab="Length of exon blocks"
ylab="Number of blocks"

graph_kn <- plot(table(data_split_89$Left),
                 #col="grey",
                 #col = rainbow(4, start = 0.7, end = 0.9),
                 col = rgb((1:4)/4, 0, 0),
                 main=paste(main,'\n',sub_main),
                 #col = c("red", "blue"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,4),
                 xlim=c(min(data_split_89$Left),max(data_split_89$Left)),
                 lwd=20,type="h",lend="butt",
                 las=1)


#####################################
# Sur les blocs exoniques de gauche #
#####################################
main=paste("Répartition de la taille des exons situés à droite de chacun des", nrow((data_split_89)))
sub_main=paste("couples d'adjacence triés en fonction du brin +")
xlab="Length of exon blocks"
ylab="Number of blocks"

graph_kn <- plot(table(data_split_89$Right),
                 #col="grey",
                 #col = rainbow(4, start = 0.7, end = 0.9),
                 col = rgb((1:4)/4, 0, 0),
                 main=paste(main,'\n',sub_main),
                 #col = c("red", "blue"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,4),
                 xlim=c(min(data_split_89$Right),max(data_split_89$Right)),
                 lwd=20,type="h",lend="butt",
                 las=1)


#####################################
# Sur les blocs exoniques confondus #
#####################################

# Lecture du tableau de l'ensemble des 89 blocs d'exons séparés
data_89 <- read.table("distri_exon_89_junctions_len_all.csv", header = TRUE, sep = ",", dec=',')
head(data_89)
summary(data_89)


main=paste("Répartition de la taille des", nrow(data_89), "blocs d'exons de chacun des", nrow(data_89)/2)
sub_main=paste("couples d'adjacence triés en fonction du brin +")
xlab="Length of exon blocks"
ylab="Number of blocks"

graph_kn <- plot(table(data_89$All_89),
                 #col="grey",
                 #col = rainbow(4, start = 0.7, end = 0.9),
                 col = rgb((1:4)/4, 0, 0),
                 main=paste(main,'\n',sub_main),
                 #col = c("red", "blue"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,8),
                 xlim=c(min(data_89$All_89),max(data_89$All_89)),
                 lwd=20,type="h",lend="butt",
                 las=1)

