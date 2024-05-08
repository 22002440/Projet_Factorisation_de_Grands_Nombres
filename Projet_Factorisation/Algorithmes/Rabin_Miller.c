#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "SQM.c"
#include "Generateur.c"

int RabinMiller(const mpz_t n, int k)
{
	//Premier test si n <= 2
	if ( mpz_cmp_ui(n, 2) < 0 )
		return 0;
	if ( mpz_cmp_ui(n, 2) == 0 )
		return 2;
		
	mpz_t a, t, t_tmp, reste, y, nMoins1;
	mpz_inits (a, t, t_tmp, reste, y, nMoins1, NULL);
	int s = 0;
	gmp_randstate_t r_state;
	unsigned long int seed = time(NULL)*(rand()%100+1);
	gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);
	
	//Décomposition de n - 1 = 2^s*t
	//t = n - 1
	mpz_sub_ui(t, n, 1);
	//reste = t%2
	mpz_mod_ui(reste, t, 2);
	//tant que t est divisible par 2
	while( mpz_cmp_ui(reste, 0) == 0)
	{
		//t = t/2
		mpz_div_ui(t, t, 2);
		s++;
		//reste = t%2
		mpz_mod_ui(reste, t, 2);
	}

	//On répète le test k fois 
	for (int i = 0; i < k ; i++)
	{
		int iContinue = 0;
		mpz_sub_ui(nMoins1, n, 1);
		
		//Générer a aléatoire dans ]1, n[
		int nbBits = mpz_sizeinbase(n, 2);
		//Boucle do while pour sécurité, on sait jamais si le générateur aléatoire a un problème
		do {
			mpz_urandomb(a, r_state, nbBits-1);
			mpz_add_ui(a, a, 1);
			//gmp_printf("a = %Zd\n", a);
		} while ( mpz_cmp_ui(a,1) <= 0 || mpz_cmp(a,nMoins1) > 0);
		
		//y = a^t mod n
		mpz_set(t_tmp, t);
		mpz_powm(y, a, t_tmp, n);
		mpz_set(t_tmp, t);
		
		if (	(mpz_cmp_ui(y,1) != 0)	&&	(mpz_cmp(y,nMoins1) != 0)	)
		{
			for ( int j = 1 ; j <= (s-1) ; j++ )
			{
				//y = y^2 mod n
				mpz_powm_ui(y, y, 2, n);
				
				if (mpz_cmp_ui(y,1) == 0) {
					mpz_clears(a, t, t_tmp, reste, y, nMoins1, NULL);
					gmp_randclear(r_state);
					return 0;
				}
				if (mpz_cmp(y,nMoins1) == 0) {
					iContinue = 1;
					break;
				}
			}
			//Dans le cas où toute la boucle du j a été itérée 
			if(!iContinue) {
				mpz_clears(a, t, t_tmp, reste, y, nMoins1, NULL);
				gmp_randclear(r_state);
				return 0;
			}
		}
	}
	
	//Probablement premier
	mpz_clears(a, t, t_tmp, reste, y, nMoins1, NULL);
	gmp_randclear(r_state);
	return 1;
}

// Main function to drive primality testing
int main() {
    int test_round = 5;
    mpz_t nombre;
    mpz_init (nombre);
    printf("Enter integer to test for primality: ");
    gmp_scanf("%Zd", &nombre);
    if (RabinMiller(nombre, test_round))
        gmp_printf("\n %Zd is prime.\n", nombre);
    else
        gmp_printf("\n %Zd is not prime.\n", nombre);
    return 0;
}
