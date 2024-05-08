#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>
#include "SQM.c"
#include "Generateur.c"


bool fermatTest(const mpz_t n, int k) {
    if (mpz_cmp_ui(n, 1) <= 0) {
        return false;
    }
    if (mpz_cmp_ui(n, 2) == 0) {
        return true;
    }

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    unsigned long seed = time(NULL);
    gmp_randseed_ui(state, seed);

    mpz_t a;
    mpz_init(a);
    mpz_t result;
    mpz_init(result);

    for (int i = 0; i < k; i++) {
        mpz_urandomm(a, state, n);
        if (mpz_cmp_ui(a, 0) == 0) {
            mpz_add_ui(a, a, 1);
        }

        square_and_multiply(result, a, n, n); // result = a^n mod n

        if (mpz_cmp(result, a) != 0) {
            mpz_clear(a);
            mpz_clear(result);
            gmp_randclear(state);
            return false; // n est compos�
        }
    }

    mpz_clear(a);
    mpz_clear(result);
    gmp_randclear(state);
    return true; // n est probablement premier
}

int main() {
    for (int i = 0; i < 10; i++) {
        mpz_t primeNumber;
        mpz_init(primeNumber);

        generateRandomNumber(primeNumber);

        printf("Nombre %d genere aleatoirement : ", i + 1);
        mpz_out_str(stdout, 10, primeNumber);
        printf("\n");

        if (fermatTest(primeNumber, 2)) {
            printf("Ce nombre est probablement premier.\n");
        } else {
            printf("Ce nombre est compose .\n");
        }

        mpz_clear(primeNumber);
    }

    return 0;
}
