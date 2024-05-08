Ce projet a été réalisé par Hugo TIERES, Elise REBER, Maxime ORLHAC, Amira HAMAD et Naoufal Amallah.

Description

Ce projet comprend des implémentations en langage C de plusieurs algorithmes fondamentaux en théorie des nombres. Les algorithmes couverts sont les suivants :

Test de primalité de Solovay-Strassen
Test de primalité de Fermat
Test de primalité de Miller-Rabin
Algorithme du crible quadratique
Algorithme du Rho de Pollard
Algorithme du p-1 de Pollard
Algorithme d'exponentiation modulaire
Algorithme de génération d'entier aléatoire

Implémentations

Chaque algorithme est implémenté dans un fichier C distinct, comme suit :

Solovay-Strassen.c : Test de primalité de Solovay-Strassen
Fermat.c : Test de primalité de Fermat
Rabin-Miller.c : Test de primalité de Rabin-Miller
crible_quadratique.c : Algorithme du crible quadratique
Rho_de_Pollard.c : Algorithme du Rho de Pollard
p-1_de_Pollard.c : Algorithme du p-1 de Pollard
SQM.c : Algorithme d'exponentiation modulaire
Generateur.c : Algorithme de génération d'entier aléatoire

Utilisation

Chaque fichier contient une fonction main qui illustre l'utilisation de l'algorithme correspondant. Vous pouvez compiler chaque fichier séparément à l'aide de votre compilateur C préféré.

Voici comment compiler et exécuter chaque algorithme :

Test de primalité de Solovay-Strassen

gcc solovay_strassen.c -o solovay_strassen -lgmp
./solovay_strassen

Test de primalité de Fermat

gcc Fermat.c -o Fermat -lgmp
./Fermat

Test de primalité de Rabin-Miller

gcc Rabin-Miller.c -o Rabin-Miller -lgmp
./Rabin-Miller

Algorithme du crible quadratique

gcc crible_quadratique.c -o crible_quadratique -lgmp
./crible_quadratique

Algorithme du Rho de Pollard

gcc Rho_de_Pollard.c -o Rho_de_Pollard -lgmp
./Rho_de_Pollard

Algorithme du p-1 de Pollard

gcc p-1_de_Pollard.c -o p-1_de_Pollard -lgmp
./p-1_de_Pollard

Algorithme d'exponentiation modulaire

gcc SQM.c -o SQM -lgmp
./SQM

Algorithme de génération d'entier aléatoire

gcc Generateur.c -o Generateur -lgmp
./Generateur

Dépendances

Ce projet utilise la bibliothèque GMP (GNU Multiple Precision Arithmetic Library) pour effectuer des opérations sur de grands entiers.

Assurez-vous d'avoir GMP installé sur votre système avant de compiler et d'exécuter les programmes.