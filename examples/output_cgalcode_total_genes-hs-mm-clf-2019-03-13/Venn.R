#install.packages('VennDiagram')
#install.packages("VennDiagram")
library(grid)
library(futile.logger)
library(VennDiagram)

#setwd("/home/Stage_M2/DATA/test_analysis_graph-hs-mm-clf-2018-05-24/")
#setwd("/home/Stage_M2/DATA/export_files_to_statistical_analysis-hs-mm-clf-2018-05-29/")
#setwd("/home/niguilla/Documents/test_donnees_stat_cliques_doublons_singleton_181011/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/test-hs-mm-clf-2018-10-11/")

#######################################################################
########################### SIGNAL ANALYSIS ###########################
#######################################################################

############################### TYPE ##################################
cliques <- readLines ("signal_repartition/clique_without_ambiguous.txt")
doublons  <- readLines ("signal_repartition/doublons.txt")
singletons  <- readLines ("signal_repartition/singletons.txt")

input<-list(clique=cliques,
            doublon=doublons,
            singleton=singletons)

n1 <- "Three species"
n2 <- "Two species"
n3 <- "One species"

color1 <- "brown2"
color2 <- "chocolate1"
color3 <- "darkgoldenrod1"

filename = "cliques_doublons_singletons_signal_en_without.png"

venn<-venn.diagram(input,
                   filename=filename,
                   height = 4000, width = 4000,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des composantes des graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.45,
                   cex = 3,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

############################### DOUBLONS ##################################
doublons_hsmm <- readLines ("signal_repartition/doublons_hsmm.txt")
doublons_hsclf <- readLines ("signal_repartition/doublons_hsclf.txt")
doublons_mmclf <- readLines ("signal_repartition/doublons_mmclf.txt")

input<-list(doublons_mmclf,
            doublons_hsclf,
            doublons_hsmm)

n1 <- "No human"
n2 <- "No mouse"
n3 <- "No dog"

color1 <- "skyblue1"
color2 <- "cornsilk3"
color3 <- "firebrick3"

filename = "doublons_signal.png"

venn<-venn.diagram(input,
                   filename=filename,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des doublons entre les graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.25,
                   cex = 2,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

############################### SINGLETONS ##################################
singletons_hs <- readLines ("signal_repartition/singletons_hs.txt")
singletons_mm <- readLines ("signal_repartition/singletons_mm.txt")
singletons_clf <- readLines ("signal_repartition/singletons_clf.txt")

input<-list(singletons_hs,
            singletons_mm,
            singletons_clf)

n1 <- "Human"
n2 <- "Mouse"
n3 <- "Dog"

color1 <- "skyblue1"
color2 <- "cornsilk3"
color3 <- "firebrick3"

filename = "singletons_signal.png"

####################### MAIN ########################
#jpeg("clique_doublon_singleton.png")
venn<-venn.diagram(input,
                   filename=filename,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des singletons entre les graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.25,
                   cex = 2,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

#######################################################################
######################### TRANSCRIPTS ANALYSIS ########################
#######################################################################

############################### TYPE ##################################
cliques <- readLines ("tr_repartition/clique.txt")
doublons  <- readLines ("tr_repartition/doublons.txt")
singletons  <- readLines ("tr_repartition/singletons.txt")

input<-list(clique=cliques,
            doublon=doublons,
            singleton=singletons)

n1 <- "Cliques"
n2 <- "Duplicates"
n3 <- "Singletons"

n1 <- "Three species"
n2 <- "Two species"
n3 <- "One species"

color1 <- "brown2"
color2 <- "chocolate1"
color3 <- "darkgoldenrod1"

filename = "cliques_doublons_singletons_tr_en.png"

# venn<-venn.diagram(input,
#                    filename=filename,
#                    #main = "Diagramme de Venn",
#                    #sub = "répartition des composantes des graphes (23/05/2018)",
#                    imagetype = "png",
#                    #main.cex = 1.5,
#                    cat.cex = 1.25,
#                    cex = 2,
#                    fill = c(color1, color2, color3),
#                    category.names = c(n1, n2, n3),
#                    col=c(color1, color2, color3),
#                    euler.d = FALSE,
#                    scaled=FALSE)

venn<-venn.diagram(input,
                   filename=filename,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des composantes des graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.25,
                   cex = 2,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

############################### DOUBLONS ##################################
doublons_hsmm <- readLines ("tr_repartition/doublons_hsmm.txt")
doublons_hsclf <- readLines ("tr_repartition/doublons_hsclf.txt")
doublons_mmclf <- readLines ("tr_repartition/doublons_mmclf.txt")

input<-list(doublons_mmclf,
            doublons_hsclf,
            doublons_hsmm)

n1 <- "No human"
n2 <- "No mouse"
n3 <- "No dog"

color1 <- "skyblue1"
color2 <- "cornsilk3"
color3 <- "firebrick3"

filename = "doublons_tr.png"

venn<-venn.diagram(input,
                   filename=filename,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des doublons entre les graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.25,
                   cex = 2,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

############################### SINGLETONS ##################################
singletons_hs <- readLines ("tr_repartition/singletons_hs.txt")
singletons_mm <- readLines ("tr_repartition/singletons_mm.txt")
singletons_clf <- readLines ("tr_repartition/singletons_clf.txt")

input<-list(singletons_hs,
            singletons_mm,
            singletons_clf)

n1 <- "Human"
n2 <- "Mouse"
n3 <- "Dog"

color1 <- "skyblue1"
color2 <- "cornsilk3"
color3 <- "firebrick3"

filename = "singletons_tr.png"

####################### MAIN ########################
#jpeg("clique_doublon_singleton.png")
venn<-venn.diagram(input,
                   filename=filename,
                   #main = "Diagramme de Venn",
                   #sub = "répartition des singletons entre les graphes (23/05/2018)",
                   imagetype = "png",
                   #main.cex = 1.5,
                   cat.cex = 1.10,
                   cex = 2,
                   fill = c(color1, color2, color3),
                   category.names = c(n1, n2, n3),
                   col=c(color1, color2, color3),
                   euler.d = FALSE,
                   scaled=FALSE)

                   #cat.dist = c(0.2,0.2,0.1,0.1),
#dev.off()
#?venn.diagram

