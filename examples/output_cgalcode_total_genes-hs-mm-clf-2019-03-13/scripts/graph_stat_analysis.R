# Lecture du tableau de l'ensemble des gènes : 112
data <- read.table("relation_gene_known_predict_accumul.csv", header = TRUE, sep = ",", dec=',')
head(data)

# Ajouter les informations sur le pourcentage de tr connus / predits pour chaque gènes de chaque espèce
data$tr_known_percent <-(data$tr_known/(data$tr_known+data$tr_pred))*100
data$tr_pr_percent <- (data$tr_pred/(data$tr_known+data$tr_pred))*100
head(data)
#library(ggplot2)

# Tableau de répartition des données pour les plot
#table(data$tr_known)
#table(data$tr_pred)

# Lecture du tableau des gènes par rapport aux espèces : 112*3
data_sp <- read.table("relation_gene_known_predict_sp.csv", header = TRUE, sep = ",", dec=',')
head(data_sp)
summary(data_sp)

# Ajouter le pourcentage de tr connus / predits pour les triplets de gènes
data_sp$tr_known_percent <- ((data_sp$hs_tr_known+data_sp$mm_tr_known+data_sp$clf_tr_known)/(data_sp$hs_tr_known+data_sp$mm_tr_known+data_sp$clf_tr_known+data_sp$hs_tr_pred+data_sp$mm_tr_pr+data_sp$clf_tr_pr))*100
data_sp$tr_pr_percent <- ((data_sp$hs_tr_pred+data_sp$mm_tr_pr+data_sp$clf_tr_pr)/(data_sp$hs_tr_known+data_sp$mm_tr_known+data_sp$clf_tr_known+data_sp$hs_tr_pred+data_sp$mm_tr_pr+data_sp$clf_tr_pr))*100

# Ajouter le pourcentage de tr connus par espèces
data_sp$tr_know_percent_hs <-(data_sp$hs_tr_known/(data_sp$hs_tr_known+data_sp$hs_tr_pred))*100
data_sp$tr_know_percent_mm <-(data_sp$mm_tr_known/(data_sp$mm_tr_known+data_sp$mm_tr_pr))*100
data_sp$tr_know_percent_clf <-(data_sp$clf_tr_known/(data_sp$clf_tr_known+data_sp$clf_tr_pr))*100

data_sp$tr_pr_percent_hs <-(data_sp$hs_tr_pred/(data_sp$hs_tr_known+data_sp$hs_tr_pred))*100
data_sp$tr_pr_percent_mm <-(data_sp$mm_tr_pr/(data_sp$mm_tr_known+data_sp$mm_tr_pr))*100
data_sp$tr_pr_percent_clf <-(data_sp$clf_tr_pr/(data_sp$clf_tr_known+data_sp$clf_tr_pr))*100

head(data_sp)

#table(data_sp$hs_tr_known)
#table(data_sp$mm_tr_known)
#table(data_sp$clf_tr_known)

#table(data_sp$tr_known_percent)
#table(data_sp$tr_pr_percent_hs)

####
# Boîte à moutache
####

#summary(data$tr_known)
#boxplot(data_sp$hs_tr_known~data_sp$hs_tr_pred,varwidth = TRUE, notch = TRUE, outline = TRUE)
#boxplot(data_sp$hs_tr_known)

# Un seul graphique par figure
par(mfrow=c(1,1))

########################################################
# Distribution des transcrits connus selon les espèces #
########################################################

#Définition des couleurs
color1 <- c(4) #Définir la couleur
colort1 <- adjustcolor(color1, alpha.f = 0.3) #rendre la couleur transparente
color2 <- c(1)
colort2 <- adjustcolor(color2, alpha.f = 0.3)
color3 <- c(2)
colort3 <- adjustcolor(color3, alpha.f = 0.3)

#Affichez le premier jeu de données
plot(table(data_sp$hs_tr_known), col=colort1, lwd=20, xlim=c(0,10), ylim=c(0,100), ylab='')
#Appelez la fonction par() avec le paramètre new=TRUE
par(new=TRUE)

#Affichez le deuxième jeu de données sur le même graphe=
plot(table(data_sp$mm_tr_known), col=colort2, lwd=20, xlim=c(0,10), ylim=c(0,100), ylab='')
#Appelez la fonction par() avec le paramètre new=TRUE
par(new=TRUE)

#Affichez le troisième jeu de données sur le même graphe et définir les détails
plot(table(data_sp$clf_tr_known), col=colort3, lwd=20, xlim=c(0,10), ylim=c(0,100), #axes=T,
     main=paste("Distribution des transcrits connus selon les espèces"),
     ylab="Nombre de gènes",
     xlab="Nombre de transcrits")

# Ajouter une légende
legend(4, 100, legend=c("Transcrits humains", "Transcrits murins", "Transcrits chiens"),
       col=c(colort1, colort2, colort3), lty=1:1,lwd=rep(5,2), cex=1)

##################################################################
# OBJECTIF: VISUALISER LE NOMBRE DE TRANSCRITS (la distribution) #
##################################################################

## Définition d'une figure de 1*4 = 4 figures
par(mfrow=c(2,4))

################################################################
# Répartition des transcrits connus parmi l'ensemble des gènes #
################################################################
main=paste("Distribution des transcrits connus pour les", nrow(data), "gènes")
xlab="Nombre de transcrits connus"
ylab="Nombre de gènes"

graph_kn <- plot(table(data$tr_known),
                 # col="grey",
                 main=main,
                 col=9,
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,200),
                 xlim=c(0,10),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn, y = table(data$tr_known), label = table(data$tr_known), pos = 3, cex = 0.9, col = "red")

########################################################################
# Répartition des transcrits connus parmi l'ensemble des gènes humains #
########################################################################
graph_kn_hs <- plot(table(data_sp$hs_tr_known),
                 # col="grey",
                 main="Dans les gènes humains",
                 col=c("#66CCFF"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,100),
                 xlim=c(0,10),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn_hs, y = table(data_sp$hs_tr_known), label = table(data_sp$hs_tr_known), pos = 3, cex = 0.9, col = "black")

#######################################################################
# Répartition des transcrits connus parmi l'ensemble des gènes murins #
#######################################################################
graph_kn_mm <- plot(table(data_sp$mm_tr_known),
                 # col="grey",
                 main="Dans les gènes murins",
                 col=8,
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,100),
                 xlim=c(0,10),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn_mm, y = table(data_sp$mm_tr_known), label = table(data_sp$mm_tr_known), pos = 3, cex = 0.9, col = "black")

#######################################################################
# Répartition des transcrits connus parmi l'ensemble des gènes chiens #
#######################################################################
graph_kn_clf <- plot(table(data_sp$clf_tr_known),
                 # col="grey",
                 main="Dans les gènes canins",
                 col=c("#FF99FF"),
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,100),
                 xlim=c(0,10),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn_clf, y = table(data_sp$clf_tr_known), label = table(data_sp$clf_tr_known), pos = 1, cex = 0.9, col = "black")

## Définition d'une figure de 2*2 = 4 figures
#par(mfrow=c(1,4))

################################################################
# Répartition des transcrits prédits parmi l'ensemble des gènes #
################################################################
main=paste("Distribution des transcrits prédits pour les", nrow(data), "gènes")
xlab="Nombre de transcrits prédits"
ylab="Nombre de gènes"

graph_kn <- plot(table(data$tr_pred),
                 # col="grey",
                 main=main,
                 col=9,
                 xlab=xlab, ylab=ylab,
                 ylim=c(0,200),
                 xlim=c(0,10),
                 lwd=20,type="h",lend="butt",
                 las=1)
## Add text at top of bars
text(x = graph_kn, y = table(data$tr_pred), label = table(data$tr_pred), pos = 3, cex = 0.9, col = "red")

########################################################################
# Répartition des transcrits prédits parmi l'ensemble des gènes humains #
########################################################################
graph_kn_hs <- plot(table(data_sp$hs_tr_pred),
                    # col="grey",
                    main="Dans les gènes humains",
                    col=c("#66CCFF"),
                    xlab=xlab, ylab=ylab,
                    ylim=c(0,100),
                    xlim=c(0,10),
                    lwd=20,type="h",lend="butt",
                    las=1)
## Add text at top of bars
text(x = graph_kn_hs, y = table(data_sp$hs_tr_pred), label = table(data_sp$hs_tr_pred), pos = 3, cex = 0.9, col = "black")

#######################################################################
# Répartition des transcrits prédits parmi l'ensemble des gènes murins #
#######################################################################
graph_kn_mm <- plot(table(data_sp$mm_tr_pr),
                    # col="grey",
                    main="Dans les gènes murins",
                    col=8,
                    xlab=xlab, ylab=ylab,
                    ylim=c(0,100),
                    xlim=c(0,10),
                    lwd=20,type="h",lend="butt",
                    las=1)
## Add text at top of bars
text(x = graph_kn_mm, y = table(data_sp$mm_tr_pr), label = table(data_sp$mm_tr_pr), pos = 3, cex = 0.9, col = "black")

#######################################################################
# Répartition des transcrits pr"dits parmi l'ensemble des gènes chiens #
#######################################################################
graph_kn_clf <- plot(table(data_sp$clf_tr_pr),
                     # col="grey",
                     main="Dans les gènes canins",
                     col=c("#FF99FF"),
                     xlab=xlab, ylab=ylab,
                     ylim=c(0,100),
                     xlim=c(0,10),
                     lwd=20,type="h",lend="butt",
                     las=1)
## Add text at top of bars
text(x = graph_kn_clf, y = table(data_sp$clf_tr_pr), label = table(data_sp$clf_tr_pr), pos = 3, cex = 0.9, col = "black")







################################################################
# OBJECTIF: VISUALISER LES TAUX DE TRANSCRITS CONNUS / PREDITS #
################################################################
par(mfrow=c(2,4))

length((data_sp$tr_known_percent))

ylab="Nombre de gènes"
xlab="Taux de transcrits (en %)"
border="white"
##################################################################################
# Répartition des fréquences de transcrits connus par rapport au nombre de gènes #
##################################################################################
hist1 <- hist(data_sp$tr_known_percent, 
     breaks = c(0,10,20,30,40,50,60,70,80,90,100),
     axes=T,
     main=paste("Histogramme du taux de transcrits connus sur les", nrow(data_sp), "gène"),
     xlab=xlab, ylab=ylab, labels=T,
     xlim=c(0,100), ylim=c(0,60),
     col="black", border=border,
     las=1, lwd=1.5)
#axis(1, 0:100)  ## Axe des x allant de 0 à 10 en pas de 1
#axis(2, 0:50)  ## Axe des y allant de -1 à 1 en pas de 1
#abline(h=40, col="red")
#abline(h=30, col="blue")
#abline(h=20, col="green")
#abline(h=10, col="orange")

################################################################################################
# Répartition des fréquences de transcrits connus par rapport au nombre de gènes chez l"humain #
################################################################################################
hist1 <- hist(data_sp$tr_know_percent_hs, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes humains"),
              xlab=xlab, ylab=ylab, labels=T,
              xlim=c(0,100),ylim=c(0,120),
              col=c("#66CCFF"), border=border,
              las=1, lwd=1.5)

#################################################################################################
# Répartition des fréquences de transcrits connus par rapport au nombre de gènes chez la souris #
#################################################################################################
hist1 <- hist(data_sp$tr_know_percent_mm, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes murins"),
              xlab=xlab, ylab=ylab,
              xlim=c(0,100), ylim=c(0,120), labels=T,
              col=8, border=border,
              las=1, lwd=1.5)

################################################################################################
# Répartition des fréquences de transcrits connus par rapport au nombre de gènes chez le chien #
################################################################################################
hist1 <- hist(data_sp$tr_know_percent_clf, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes canins"),
              xlab=xlab, ylab=ylab, labels=T,
              xlim=c(0,100),
              ylim=c(0,120),
              col=c("#FF99FF"), border=border,
              las=1, lwd=1.5)

##################################################################################
# Répartition des fréquences de transcrits prédits par rapport au nombre de gènes #
##################################################################################
hist1 <- hist(data_sp$tr_pr_percent, 
     breaks = c(0,10,20,30,40,50,60,70,80,90,100),
     axes=T,
     main=paste("Histogramme du taux de transcrits prédits sur les", nrow(data_sp), "gène"),
     xlab=xlab, ylab=ylab, labels=T,
     xlim=c(0,100), ylim=c(0,60),
     col="black", border=border,
     las=1, lwd=1.5)

################################################################################################
# Répartition des fréquences de transcrits prédits par rapport au nombre de gènes chez l"humain #
################################################################################################
hist1 <- hist(data_sp$tr_pr_percent_hs, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes humains"),
              xlab=xlab, ylab=ylab, labels=T,
              xlim=c(0,100),ylim=c(0,120),
              col=c("#66CCFF"), border=border,
              las=1, lwd=1.5)

#################################################################################################
# Répartition des fréquences de transcrits prédits par rapport au nombre de gènes chez la souris #
#################################################################################################
hist1 <- hist(data_sp$tr_pr_percent_mm, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes murins"),
              xlab=xlab, ylab=ylab,
              xlim=c(0,100), ylim=c(0,120), labels=T,
              col=8, border=border,
              las=1, lwd=1.5)

################################################################################################
# Répartition des fréquences de transcrits prédits par rapport au nombre de gènes chez le chien #
################################################################################################
hist1 <- hist(data_sp$tr_pr_percent_clf, 
              breaks = c(0,10,20,30,40,50,60,70,80,90,100),
              axes=T,
              main=paste("Dans les gènes canins"),
              xlab=xlab, ylab=ylab, labels=T,
              xlim=c(0,100),
              ylim=c(0,120),
              col=c("#FF99FF"), border=border,
              las=1, lwd=1.5)





######################################################################################
# OBJECTIF: VISUALISER LES TAUX DE TRANSCRITS CONNUS / PREDITS DE MANIERE CUMULATIVE #
######################################################################################
par(mfrow=c(2,4))
length(data_sp$tr_known_percent_cat) #187
xlab="Taux de transcrits cumulé"
ylab="Taux de transcrits cumulé (en %)"

border="black"

#######################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes #
#######################################################################################
data_sp$tr_know_Pourc.classe <- cut(data_sp$tr_known_percent, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$tr_know_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$tr_know_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Taux de transcrits connus pour les",nrow(data),"gènes"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col="black", border="white",
                  las=1, lwd=1.5)

#effectif <- table(data_sp$tr_know_percent_hs_cat)
#effectif
#frequence <- effectif/187*100
#frequence
#effcum <- cumsum(effectif)
#effcum
#barplot(effcum, ylim = c(0,200))

#abline(h=60, col="red")
#abline(h=30, col="blue")
#abline(h=20, col="green")
#abline(h=10, col="orange")

####################################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes chez l'homme #
####################################################################################################
data_sp$hs_know_Pourc.classe <- cut(data_sp$tr_know_percent_hs, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$hs_know_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$hs_know_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes humains"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=c("#66CCFF"), border=border,
                  las=1, lwd=1.5)

######################################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes chez la souris #
######################################################################################################
data_sp$mm_know_Pourc.classe <- cut(data_sp$tr_know_percent_mm, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$mm_know_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$mm_know_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes murins"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=8, border=border,
                  las=1, lwd=1.5)
           
######################################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes chez le chien #
######################################################################################################
data_sp$clf_know_Pourc.classe <- cut(data_sp$tr_know_percent_clf, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$clf_know_Pourc.classe

#plot(data_sp$clf_know_Pourc.classe)

effectif <- table(data_sp$clf_know_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes canins"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=c("#FF99FF"), border=border,
                  las=1, lwd=1.5)

#######################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes #
#######################################################################################
data_sp$tr_pr_Pourc.classe <- cut(data_sp$tr_pr_percent, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$tr_pr_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$tr_pr_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Taux de transcrits connus pour les",nrow(data),"gènes"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col="black", border="white",
                  las=1, lwd=1.5)

####################################################################################################
# Répartition des taux cumulatifs de transcrits prédits par rapport au nombre de gènes chez l'homme #
####################################################################################################
data_sp$hs_pr_Pourc.classe <- cut(data_sp$tr_pr_percent_hs, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$hs_pr_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$hs_pr_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes humains"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=c("#66CCFF"), border=border,
                  las=1, lwd=1.5)

######################################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes chez la souris #
######################################################################################################
data_sp$mm_pr_Pourc.classe <- cut(data_sp$tr_pr_percent_mm, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$mm_pr_Pourc.classe

#plot(data_sp$hs_know_Pourc.classe)

effectif <- table(data_sp$mm_pr_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes murins"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=8, border=border,
                  las=1, lwd=1.5)

######################################################################################################
# Répartition des taux cumulatifs de transcrits connus par rapport au nombre de gènes chez le chien #
######################################################################################################
data_sp$clf_pr_Pourc.classe <- cut(data_sp$tr_pr_percent_clf, c(0,10,20,30,40,50,60,70,80,90,100), right=F,include.lowest = TRUE)
data_sp$clf_pr_Pourc.classe

#plot(data_sp$clf_know_Pourc.classe)

effectif <- table(data_sp$clf_pr_Pourc.classe)
effectif

effcum <- cumsum(effectif)
effcum

my_bar <- barplot(effcum,
                  #breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                  main=paste("Dans les gènes canins"),
                  #xlab=xlab, ylab=ylab, labels=T,
                  ylim=c(0,200),
                  #xlim = c(0,15), ylim = c(0,200),
                  col=c("#FF99FF"), border=border,
                  las=1, lwd=5)




# ######################################
# ########## POUBELLE DE TEST ##########
# ######################################
# 
# my_hist_test=plot((effcum), border=T, las=2)
# # Add the text 
# text(my_bar, data$average+0.4 , paste("n = ",data_sp$number,sep="") ,cex=1) 
# 
# x <- rnorm(100)
# plot(ecdf(x))
# 
# classes<-seq(180,230,10)
# classes
# effectifs<-c(2,26,35,60,7,0)
# effectifs
# taille<-rep(classes,effectifs)
# taille
# hist(taille,breaks=classes,right=FALSE,main="",xlab="Taille",ylab="Effectif",col="magenta")
# 
# x=rnorm(1000)
# plot(effcum,cumsum(effcum),type="b")
# 
# ?barplot
# 
# name= c("DD","with himself","with DC","with Silur" ,"DC","with himself","with DD","with Silur" ,"Silur","with himself","with DD","with DC" )
# average= sample(seq(1,10) , 12 , replace=T)
# number= sample(seq(4,39) , 12 , replace=T)
# data=data.frame(name,average,number)
# 
# # Basic Barplot
# my_bar=barplot(data$average , border=F , names.arg=data$name , las=2 , col=c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6) , rgb(0.3,0.9,0.4,0.6) ,  rgb(0.3,0.9,0.4,0.6)) , ylim=c(0,11) , main="" )
# abline(v=c(4.9 , 9.7) , col="grey")
# 
# 
# classes<-seq(180,230,10)
# classes
# effectifs<-c(2,26,35,60,7,0)
# effectifs
# taille<-rep(classes,effectifs)
# taille
# hist(taille,breaks=classes,right=FALSE,main="",xlab="Taille",ylab="Effectif",col="magenta")
# 
# ## données de l'histogramme
# (x <- hist(taille, breaks = classes, right = FALSE, plot = FALSE))
# ## nb de cm pour abcsisse: 2 cm pour 10 sinon largeur histo = 1 cm !!!
# xcm <- 2 * diff(range(x$breaks)) / 10
# ## nb de cm pour ordonnée
# ycm <- max(x$counts) / 10
# ## graphique
# par(pin = c(xcm / 2.54, ycm / 2.54))
# hist(taille, breaks = classes, right = FALSE, freq = TRUE,
#      main = "", xlab = "Taille", ylab = "Effectif", col = "magenta")
# 
# #Legende
# legend("topleft", legend = c("Alone","with Himself","With other genotype" ) , 
#        col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6) , rgb(0.3,0.9,0.4,0.6) ,  rgb(0.3,0.9,0.4,0.6)) , 
#        bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))
# 
#               
#               
# 
# #plot(table(data_sp$tr_known_percent))
# 
# # Nombre de triplets de gènes ayant des transcrits prédits entre 0 (compris) et 20 non compris : [0:20[
# length(0 <= data_sp$tr_pr_percent[data_sp$tr_pr_percent < 20])
# 
