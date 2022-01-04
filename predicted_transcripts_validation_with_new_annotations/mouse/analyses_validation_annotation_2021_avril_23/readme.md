# Analyse facilité des résultats de validation sur de nouvelles annotations.

1. Extraire les fichiers résultats et renommer en result_bedtools_intersect_source_species.out
2. awk -F\\t '{print $4}' result_bedtools_intersect_source_species.out > result_via_source_species.txt
3. Utiliser le script analyses_validation.py pour avoir les infos + un fichier résultat au format csv




TODO: faire un csv pour la base de données avec les nouvelles données
