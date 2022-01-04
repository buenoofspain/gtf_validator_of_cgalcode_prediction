library("vegan")
data <- read.csv("test.csv", sep='\t')
summary(data)

#?specaccum
#?rarecurve
specaccum(data)
md <- specaccum(data, method = "rarefaction") 
plot(data)
plot(md, ci.type="poly",col="darkblue", lwd=2, ci.lty=0, ci.col="lightblue",xlab="Nb_data (#read)", ylab="Nombre d'observations (#jonctions)") 

md
