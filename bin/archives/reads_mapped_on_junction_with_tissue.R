library(ggplot2)

# Lecture du tableau de l'ensemble des blocs d'exons
data <- read.table("nb_reads_mapped_tissue.csv", header = TRUE, sep = "\t", dec=',')
head(data)
summary(data)
data[1]
length(data)

cbind(data, total = rowSums(data[2:16]))

hist(data[2:16])

data_cp <- data[2:16]
summary(data_cp)

barplot(data$ADRENALGLAND, data$BLOOD, data$BRAIN, xlab(data$Junctions))

ggplot(data) +
  geom_col(aes(x = data$Junctions, y=c(data$ADRENALGLAND)))


ggplot(gather(data_cp, cols, value), aes(x = value)) + 
  geom_histogram(binwidth = 20) + facet_grid(.~cols)

for (col in 2:ncol(data)) {
  hist(data[,col], breaks=10)
}

ggplot(gather(mtcars), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

