Code crée par Ammi Haroun, Belyazid Ali & de Wouters Louise
Projet cpa - remise le 18 décembre 2019

Compilation : make

Exécution : ./projet proteine.fasta database.fasta

où :
- 'proteine.fasta' contient la protéine de base
- les trois fichiers de la database sont nommés 'database.fasta.pin' 'database.fasta.psq' et 'database.fasta.phr'
- le fichier texte pour écrire les informations et les résulats est nommé 'results.txt'
- le fichier de la matrice BLOSUM est nommé 'BLOSUM62'

Temps d'éxécution approximatif pour la protéine P00533 en supposant que le swipe prend 1 seconde : t = 70-75 secondes

NOTE :
Nous avons essayé d'optimiser notre code avec l'utilisation de thread. Malheureusement, nous n'avons pas réussi à insérer cette amélioration dans le code principal, structuré en différentes classes. En effet, la présence de classes complique l'utilisation des thread et nous n'avons pas eu le temps de travailler d'avantage dessus. Le code nommé "projet_opt.cpp" reprend cette optimisation. C'est un code peu structuré qui renvoie uniquement les 10 premiers meilleurs scores dans le terminal, avec leur place respective. 
Compilation : g++ -o projet_opt projet_opt.cpp -lpthread
Execution : ./projet2 proteine.fasta
où on suppose que 'proteine.fasta' contient la protéine de base, que les trois fichers de la database sont nommées 'database.fasta.pin' 'database.fasta.psq' et 'database.fasta.phr' et que le fichier de la matrice BLOSUM est nommé 'BLOSUM62'
