#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include "SQM.c"
#include "Generateur.c"

/*calcule le symbole de jacobi (a/n), n étant le nombre qu'on  veut déterminer premier.
source de l'algorithme : https://www-fourier.ujf-grenoble.fr/~lefourns/contenu/cours/Agreg2018/CourssymboleLegendre.pdf*/
int symbole_jacobi(mpz_t a, mpz_t n) {

    int parite;
    parite = mpz_even_p(n);
    if (mpz_cmp_ui(n, 0) <= 0 || parite) { // on vérifie que n est supérieur à zéro et impair
        return 0; } // si ce n'est pas le cas, le symbole de jacobi est nul

    
    mpz_t jacobi; // on crée une variable symbole de jacobi dans laquelle on stockera notre résultat et on l'initialise à 1
    mpz_init_set_str(jacobi, "1", 10);

    while (mpz_cmp_ui(a, 0) != 0) { 

        if (mpz_cmp(a,n)>=0) { // cas où a est supérieur à n
        
            /*on calcule le modulo de a par n 
            afin d'obtenir quelque chose de plus petit.
            si a est un multiple de n, le symbole de jacobi est égal à zéro.*/
            mpz_t reste;
            mpz_init(reste);
            mpz_mod(reste, a,n);
            mpz_set(a, reste); // on calcule le reste de a divisé par n et on le stocke dans la variable a
            if (mpz_cmp_ui(a,0) == 0) { // si a est un multiple de n, le symbole de jacobi est zéro
                return 0;
            }
        }
        while (mpz_even_p(a) != 0) { // si a = 2^k*a', on répète la boucle k fois jusqu'à ce que a = a' impair
            mpz_t quotient;
            mpz_init(quotient);
            mpz_div_ui(quotient, a, 2);
            mpz_set(a, quotient); // on divise a par deux pour obtenir a = 2^(k-1)*a'

            /* 
            on fait le changement de signe du symbole de jacobi : 
            
            */
            mpz_t reste;
            mpz_init(reste);
            mpz_mod_ui(reste, n, 8);
            if (mpz_cmp_ui(reste,3) == 0 | mpz_cmp_ui(reste,5) == 0) {
                mpz_neg(jacobi, jacobi);
            }        
            // on calcule jacobi = jacobi * -1^((n^2 -1)/8)
            mpz_t temp;
            mpz_t base;
            mpz_init(temp);
            mpz_init(temp);
            mpz_init_set_str(base, "-1", 10);
            // on initialise l'exposant
            mpz_t exp;
            mpz_init(exp);
            mpz_pow_ui(exp, n, 2);
            mpz_sub_ui(exp, exp, -1);
            mpz_div_ui(exp, exp, 8);
            // si l'exposant est pair : la base devient 1 et jacobi reste identique
            mpz_mod_ui(exp, exp, 2);
            if (mpz_odd_p(exp) == 1) { // si l'exposant est impair : la base reste -1 et jacobi change de signe
                mpz_neg(jacobi, jacobi);
            }
             
        }

        // quand a n'est plus un multiple de 2
        // on calcule jacobi = jacobi * -1^((a-1)(b-1)/8)
        mpz_t temp;
        mpz_t base;
        mpz_init(temp);
        mpz_init(base);
        mpz_init_set_str(base, "-1", 10);
        mpz_t a1;
        mpz_init(a1);
        mpz_init_set_str(a1, "-1", 10);
        mpz_add(a1, a, a1);
        mpz_t n1;
        mpz_init(n1);
        mpz_init_set_str(n1, "-1", 10);
        mpz_add(n1, n, n1);
        mpz_t mul;
        mpz_init(mul);
        mpz_mul(mul, a1, n1);
        mpz_t exp;
        mpz_init(exp);
        mpz_div_ui(exp, mul, 8);
        if (mpz_odd_p(exp) == 1) {
            mpz_neg(jacobi, jacobi);
        }

        // on intervertit les valeurs de a et n
        mpz_t echange;
        mpz_init(echange);
        mpz_set(echange, a);
        mpz_set(a, n);
        mpz_set(n, echange);
        // si le nouveau a est égal à zéro, c'est la fin de la boucle, sinon elle recommence
    }
    
    // maintenant que a = 0, il y a deux cas possibles
    if (mpz_cmp_ui(n,1) == 0) { // a et n sont premiers entre eux 
        int jacobi_int;
        if (mpz_cmp_ui(jacobi, 1) == 0) {
            jacobi_int = 1;
        }
        else if (mpz_cmp_ui(jacobi, -1) == 0) {
            jacobi_int = -1;
        }
        else {
            jacobi_int = 0;
        }
        //unsigned long int jacobi_int;
        return jacobi_int;// le symbole de jacobi est égal à 1 ou -1
    } 
    // a et n ne sont pas premiers entre eux 
    return 0; // le symbole de jacobi est égal à 0
}




/*
fonction effectuant le test de primalité de solovay-strassen
k fois sur un nombre n.
renvoie 1 si le nombre est premier et 0 sinon.
*/
int testSolovay (const mpz_t n, int k)
{	
	//Premier test si n est inférieur ou égal à 2
	if (mpz_cmp_ui(n, 2) < 0 || mpz_cmp_ui(n, 2) == 0) {
		return 0;
	}
	
	//on teste la parité de n
	mpz_t parite;
	mpz_init(parite);
	mpz_mod_ui(parite, n, 2);
	if ( mpz_cmp_ui(parite, 0) == 0 ) {
		// si n est pair, on renvoie zéro car n est composé
		return 0;
	}
	mpz_clear(parite);
	
	// Si n est impair
	// ici on mets en place une graine pour le générateur et on initialise les différentes variables
	mpz_t a, n_1, exposant, euler, reste;
	mpz_inits (a, n_1, exposant, euler, reste, NULL);
	gmp_randstate_t r_state;
	unsigned long int seed = time(NULL)*(rand()%100+1);
	gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);
	
	//On répète le test k fois 
	for (int i = 0; i < k ; i++)
	{
		// Générer a aléatoire dans [2, n-1]

		mpz_sub_ui(n_1, n, 1);
		int nbBits = mpz_sizeinbase(n, 2);
		//Boucle do while pour s'assurer que le nombre aléatoire est inférieur à n
		do { 
			mpz_urandomb(a, r_state, nbBits-1);
			mpz_add_ui(a, a, 1);
		} while ( mpz_cmp_ui(a,1) <= 0 || mpz_cmp(a,n_1) > 0);
		
		//on calcule le symbole de jacobi a/n
		int jacobi = mpz_jacobi(a, n);

		//int jacobi = symbole_jacobi(a, n);
		
		//on calcule le critère d'Euler = a^(n-1/2)
		mpz_cdiv_q_ui(exposant, n_1, 2);
		mpz_mod_ui(reste, n_1, 2);
		mpz_sub(exposant, exposant, reste);
		square_and_multiply(euler, a, exposant, n);
		
		if ( jacobi == 0 // le symbole de jacobi est égal à 0
			|| (jacobi == 1 && mpz_cmp_ui(euler, jacobi) != 0) // le symbole de jacobi est égal à 1 et le critère d'euler est différent du symbole de jacobi
			|| (jacobi == -1 && mpz_cmp(euler, n_1) != 0) ) // le symbole de jacobi est égal à -1 et le critère d'euler est différent de n-1 congrut à = -1 modulo n
		{
			// n n'est pas premier et on renvoie 0
			mpz_clears(a, n_1, exposant, euler, reste, NULL);
			gmp_randclear(r_state);
			return 0;
		}
	
	}
	// n est premier, on renvoie 1
	mpz_clears(a, n_1, exposant, euler, reste, NULL);
	gmp_randclear(r_state);
	return 1;
}

int main() {
	// test de la fonction testSolovay
	// avec un nombre premier
    char nombre_teste [] = "95 647 806 479 275 528 135 733 781 266 203 904 794 419 563 064 407";
    mpz_t constant;
    mpz_init_set_str(constant, nombre_teste, 10);
    int solovay = testSolovay(constant, 10);
	if (solovay == 1) {
		char message [] = "le nombre est premier. \n";
		printf("%s", message);
	}
	else {
		char message [] = "le nombre est composé. \n";
		printf("%s", message);
	}
    printf("nombre testé :%s,\n résultat du test : %d\n", nombre_teste, solovay);

	// avec un très grand nombre généré aléatoirement
	generateRandomNumber(constant);
	int solovay2 = testSolovay(constant, 10);
	if (solovay2 == 1) {
		char message [] = "le nombre est premier. \n";
		printf("%s", message);
	}
	else {
		char message [] = "le nombre est composé. \n";
		printf("%s", message);
	}
    gmp_printf("nombre testé :%Zd,\nrésultat du test : %d\n", constant, solovay2);

	// test de la fonction symbole_jacobi avec deux nombres premiers entre eux
	mpz_t a;
    mpz_t n;
    mpz_init_set_str(a, "15", 10);
    mpz_init_set_str(n, "2", 10);
	//on vérifie si symbole_jacobi renvoie la même chose que la fonction mpz_jacobi native à gmp
    printf("symbole de jacobi: %d \n", symbole_jacobi(a, n));
    if (symbole_jacobi(a, n) == mpz_jacobi(a, n)) {
        char texte [] = "correct";
        printf("%s", texte);
    }
    else {
        char texte [] = "incorrect";
        printf("%s", texte);
    }
	// ici, on voit en exécutant que la fonction symbole_jacobi ne fonctionne pas correctement

}