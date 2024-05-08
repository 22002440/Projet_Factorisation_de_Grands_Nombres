#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>

void generateRandomNumber(mpz_t randomNumber) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    unsigned long seed = (unsigned long)(tv.tv_sec * 1000 + tv.tv_usec / 1000); // Convertir l'heure actuelle en millisecondes

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, seed); 
    mpz_urandomb(randomNumber, state, 2048);
    gmp_randclear(state);
}