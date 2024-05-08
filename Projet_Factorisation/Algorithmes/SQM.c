#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>

void square_and_multiply(mpz_t res, const mpz_t base, const mpz_t exp, const mpz_t m) {
    mpz_set_ui(res, 1); // Initialize result to 1
    mpz_t temp;
    mpz_init(temp);

    // Iterate through the bits of the exponent
    for (int i = mpz_sizeinbase(exp, 2) - 1; i >= 0; i--) {
        mpz_mul(res, res, res); // Square
        mpz_mod(res, res, m); // Modulo

        // If the current bit of the exponent is 1, multiply by base
        if (mpz_tstbit(exp, i)) {
            mpz_mul(temp, res, base);
            mpz_mod(res, temp, m); // Modulo
        }
    }

    mpz_clear(temp);
}