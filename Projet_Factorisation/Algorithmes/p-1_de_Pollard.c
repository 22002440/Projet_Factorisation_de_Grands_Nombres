#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "SQM.c"
#include "Generateur.c"

// Fonction pour calculer le PGCD de deux nombres
void gcd(mpz_t result, mpz_t a, mpz_t b) {
    mpz_gcd(result, a, b);
}

// Algorithme p-1 de Pollard pour la factorisation
void p_minus_1_pollard(mpz_t factor, mpz_t n, unsigned long B) {
    mpz_t a, b, tmp, g;
    mpz_init(a);
    mpz_init(b);
    mpz_init(tmp);
    mpz_init(g);

    // Choisir une valeur aléatoire pour a
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_urandomm(a, state, n);

    // Initialiser b à a
    mpz_set(b, a);

    // Itérer jusqu'à B
    unsigned long i;
    for (i = 2; i < B; i++) {
        // Calculer a^b mod n
        square_and_multiply(tmp, a, b, n);

        // Calculer le PGCD(a^b - 1, n)
        mpz_sub_ui(tmp, tmp, 1);
        gcd(g, tmp, n);

        // Si le PGCD est non trivial, on a trouvé un facteur
        if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, n) != 0) {
            mpz_set(factor, g);
            break;
        }

        // Incrémenter b
        mpz_add_ui(b, b, 1);
    }

    // Libérer la mémoire allouée
    mpz_clears(a, b, tmp, g, NULL);
}

// Fonction pour tester si un nombre est probablement premier
int is_probably_prime(mpz_t n, int k) {
    return mpz_probab_prime_p(n, k);
}

int main() {
    mpz_t n, factor;
    mpz_inits(n, factor, NULL);

    // Demander à l'utilisateur d'entrer le nombre à factoriser
    printf("Entrez le nombre à factoriser : ");
    char input[100];
    fgets(input, sizeof(input), stdin);
    mpz_set_str(n, input, 10);

    // Bornes B pour l'algorithme p-1 de Pollard
    unsigned long B = 10000000; // Augmentation de la valeur de B

    // Test si n est probablement premier
    if (is_probably_prime(n, 25)) {
        printf("Le nombre est probablement premier : %s\n", mpz_get_str(NULL, 10, n));
        return 0;
    }

    // Exécuter l'algorithme p-1 de Pollard jusqu'à ce que le nombre à factoriser soit un nombre premier lui-même
    while (!is_probably_prime(n, 25)) {
        p_minus_1_pollard(factor, n, B);
        // Affichage du facteur trouvé
        gmp_printf("Facteur trouvé : %Zd\n", factor);
        // Diviser n par le facteur trouvé pour répéter l'algorithme avec le quotient
        mpz_divexact(n, n, factor);
    }

    // Affichage du dernier facteur trouvé (nombre premier)
    printf("Dernier facteur premier trouvé : %s\n", mpz_get_str(NULL, 10, n));

    // Libérer la mémoire allouée
    mpz_clears(n, factor, NULL);

    return 0;
}