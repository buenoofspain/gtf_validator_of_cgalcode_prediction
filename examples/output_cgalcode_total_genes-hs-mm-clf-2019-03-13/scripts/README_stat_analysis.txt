# Réalisation des diagrammes de Venn

> Utiliser le script "Venn.R", il enregistrera directement les figures dans le répertoire courant.


# Réalisation des analyses statistiques

> Tout n'est pas automatisés, voici les étapes pour réaliser l'analyse

## Etape 1. Génération des fichiers utilisés :

> Utiliser le script "Proportion_tr_known_analysis.py" pour générer les fichiers nommés:
>> graph_tr_analysis_detail.csv
>> relation_gene_known_predict.csv

> Ce script à besoin des éléments suivant à lui transmettre en interne : 
>> benchmark : dossier où sont contenus les fichiers (REGLER SUR LE DOSSIER COURANT)
>> genes_list : la liste de gènes hs à analyser
>> graph_analysis : la liste d'information provenant de "graph_tr_analysis.csv" et ne contenant que les informations sur la liste d'intérêt (un tri des colonnnes et une suppression des informations non nécessaires est à réaliser)
>> tr_otr_relationship : correspond au fichier "corresponding_tr_otr.txt" (A NE PAS MODIFIER)
>> dog : lien vers le dossier ayant les informations sur les gènes/transcrits canins
>> mouse : lien vers le dossier ayant les informations sur les gènes/transcrits murins
>> human : lien vers le dossier ayant les informations sur les gènes/transcrits humains

> A partir du fichier "relation_gene_known_predict.csv", copier le et renommer la copie "relation_gene_known_predict_sp.csv".

> Modifier le premier fichier "relation_gene_known_predict.csv" de manière à regrouper tous les gènes en pemière colonne, le nombre de transcrits connus en deuxième colonne et le nombre de transcrits prédits en troisième colonne. Renommer les colonnes : gene, tr_known, tr_pred et renommer le fichier "relation_gene_known_predict_accumul.csv"

## Etape 2. Analyse

> Utiliser le script "graph_stat_analysis.R" qui utilise les fichier "relation_gene_known_predict_sp.csv" et "relation_gene_known_predict_accumul.csv". Ce script vous permet d'obtenir les figures d'analyses statistiques que vous pourrez exporter.

# L'ensemble des figures a été placé dans un dossier nommé "figures_analyses_stat". La taille d'exportation des diagrammes de Venne est de 1500 x 700.
