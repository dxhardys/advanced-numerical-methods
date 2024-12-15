# Nom du fichier de données
datafile = "results.txt"

# Nom du fichier de sortie PNG
outputfile = "plot.png"

# Titre du graphique
set title "NxN Matrix multiplication time"

# Étiquettes des axes
set ylabel "Time(s)"
set xlabel "Matrix Size N"

# Activer la grille
set grid

# Type de tracé et style
# Vous pouvez modifier le style, la couleur, etc., selon vos besoins
plot datafile using 1:2 with lines
 

# Exporter le graphique en tant qu'image PNG
set terminal pngcairo enhanced font 'Verdana,10'
set output outputfile
replot