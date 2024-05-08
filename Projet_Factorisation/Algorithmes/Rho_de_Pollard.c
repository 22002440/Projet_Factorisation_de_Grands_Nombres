#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

// Fonction pour calculer le PGCD de deux nombres
void gcd(mpz_t result, mpz_t a, mpz_t b) {
    mpz_gcd(result, a, b);
}

// Fonction pour calculer le prochain élément dans la séquence
void next_element(mpz_t result, mpz_t x, mpz_t n) {
    mpz_mul(result, x, x);
    mpz_add_ui(result, result, 1);
    mpz_mod(result, result, n);
}

// Algorithme rho de Pollard pour la factorisation
void rho_pollard(mpz_t factor, mpz_t n) {
    mpz_t x, y, d;
    mpz_inits(x, y, d, NULL);

    // Choisir une valeur aléatoire pour x et y
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_urandomb(x, state, 256); // Génère un nombre aléatoire de 256 bits

    mpz_set(y, x);

    mpz_set_ui(factor, 1);

    // Boucle principale de l'algorithme rho de Pollard
    while (mpz_cmp_ui(factor, 1) == 0) {
        next_element(x, x, n);
        next_element(y, y, n);
        next_element(y, y, n);

        mpz_sub(d, x, y);
        mpz_abs(d, d);
        gcd(factor, d, n);
    }

    // Libérer la mémoire allouée
    mpz_clears(x, y, d, NULL);
}

int main() {
    mpz_t n, factor, quotient;
    mpz_inits(n, factor, quotient, NULL);

    // Nombre à factoriser (input de l'utilisateur)
    printf("Entrez un nombre à factoriser : ");
    char input[1000];
    scanf("%s", input);
    if (mpz_set_str(n, input, 10) != 0) {
        printf("Erreur lors de la lecture du nombre.\n");
        mpz_clears(n, factor, quotient, NULL);
        return 1;
    }

    // Copie du nombre à factoriser pour être utilisé comme quotient
    mpz_set(quotient, n);

    // Affichage des facteurs premiers
    printf("Facteurs de %s : ", mpz_get_str(NULL, 10, n));

    // Appel de l'algorithme rho de Pollard pour chaque facteur premier
    while (mpz_cmp_ui(quotient, 1) != 0) {
        rho_pollard(factor, quotient);
        // Affichage du facteur trouvé
        gmp_printf("%Zd ", factor);
        // Diviser n par le facteur trouvé pour répéter l'algorithme avec le quotient
        mpz_divexact(quotient, quotient, factor);
    }
    printf("\n");

    // Libérer la mémoire allouée
    mpz_clears(n, factor, quotient, NULL);

    return 0;
}
