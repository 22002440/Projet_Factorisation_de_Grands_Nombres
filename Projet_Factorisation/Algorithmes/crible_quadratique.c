#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

// fonction calculant le polynôme Q(x) = x^2 - n
void Q(mpz_t result, mpz_t x, mpz_t n) {
    mpz_mul(x, x, x); // mets x au carré
    mpz_sub(x, x, n); // soustrait n à x^2
    mpz_set(result, x); // stocke le résultat dans result
}

// Fonction pour vérifier si un nombre est premier
bool estPremier(int n) {
    if (n <= 1) { // vérifie que le nombre est égal ou supérieur à 2
        return false;
    }
    for (int i = 2; i * i <= n; i++) { // crible d'érasthothène
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}

// Fonction pour remplir le tableau avec les nombres premiers inférieurs ou égaux à B
void remplirTableauPremiers(int tableau[], int B, int *taille) {
    for (int i = 2; i <= B; i++) {
        if (estPremier(i)) {
            tableau[(*taille)++] = i;
        }
    }
}

void remplirTableauQ(mpz_t tableau[], mpz_t n, int m, int *taille) {
    // stocke la racine carrée de n dans une variable x
    mpz_t x;
    mpz_init(x);
    mpz_sqrt(x, n);
    mpz_add_ui(x, x, 1);
    for (int i = 0; i < m; i++) {
        mpz_t temp;
        mpz_init(temp);
        mpz_set(temp, x);
        mpz_add_ui(temp, temp, i);
        mpz_t result;
        mpz_init(result);
        Q(result, temp, n);
        mpz_set(tableau[(*taille)++], result);        
        mpz_clear(result);
        mpz_clear(temp);
    }
    mpz_clear(x);
}

void factoriser(mpz_t tableau[], mpz_t n, int B, int m, int *taille) {
    // on commence par remplir un tableau avec les résultats des polynomes et calculer quels résultats sont B-friables

    int tableau_premiers[B];
    int taille_entiers = 0;

    remplirTableauPremiers(tableau_premiers, B, &taille_entiers);

    int taille_Q = 0;
    mpz_t tableau_polynomes[m];
    for (int i = 0; i < m; i++) {
        mpz_init(tableau_polynomes[i]);
    }
    remplirTableauQ(tableau_polynomes, n, m, &taille_Q);

    mpz_t tableau_initial[*taille];
    int taille_Q_2 = 0;
    for (int i = 0; i < *taille; i++) {
        mpz_init(tableau_initial[i]);
    }
    remplirTableauQ(tableau_initial, n, m, &taille_Q_2);


    for (int pass = 0; pass < 3; pass++) {

        for (int i = 0; i < taille_Q; i++) {

            for (int j = 0; j < taille_entiers; j++) {
                // Si x n'est pas divisible par le premier, passez au prochain premier
                if (!mpz_divisible_ui_p(tableau_polynomes[i], tableau_premiers[j])) {
                    continue;
                }
                mpz_t temp;
                mpz_init(temp);
                mpz_set(temp, tableau_polynomes[i]);
                mpz_div_ui(temp, temp, tableau_premiers[j]);
                mpz_set(tableau_polynomes[i], temp);
                mpz_clear(temp);
            }

        }
    }



    // on utilise ensuite les polynomes Q friables pour factoriser n
    gmp_printf("On factorise n = %Zd\n", n);

    for (int j = 0; j < taille_Q; j++) {
        if (mpz_cmp_ui(tableau_polynomes[j], 1) == 0) {
            mpz_t temp_x;
            mpz_init(temp_x);
            mpz_sqrt(temp_x, n);
            mpz_add_ui(temp_x, temp_x, j);
            mpz_t temp_q;
            mpz_init(temp_q);
            mpz_set_ui(temp_q, 0);
            Q(temp_q, temp_x, n);
            mpz_t pgc;
            mpz_init(pgc);
            mpz_set_ui(pgc, 0);
            mpz_gcd(pgc, n, temp_q);
            
            if (mpz_cmp(pgc, n) != 0 && mpz_cmp_ui(pgc, 1) != 0) {
                if (estPremier(mpz_get_ui(pgc))) {
                    gmp_printf("Un facteur premier de n trouvé : %Zd\n", pgc);
                    mpz_divexact(n, n, pgc); // Réduire n
                    gmp_printf("n réduit : %Zd\n", n);
                    if (estPremier(mpz_get_ui(n))) {
                        gmp_printf("dernier facteur de n : %Zd\n", n);
                        break;
                    }
                }
                else {
                    mpz_divexact(n, n, pgc); // Réduire n
                    mpz_t temp;
                    int count = 0;
                    mpz_init(temp);
                    mpz_set(temp, pgc);
                    int entier;
                    while (mpz_cmp_ui(temp, 1) > 0) {
                        for (int i = 0; i < taille_entiers; i++) {
                            mpz_t temp_2;
                            mpz_init(temp_2);
                            mpz_mod_ui(temp_2, temp, tableau_premiers[i]);
                            if (mpz_cmp_ui(temp_2, 0) == 0) {
                                entier = tableau_premiers[i];
                                mpz_div_ui(temp, temp, tableau_premiers[i]);
                                count++;
                            }
                        }
                    }
                    printf("Un facteur non premier de n trouvé : %d puissance %d \n", entier, count);
                    
                    gmp_printf("n réduit : %Zd\n", n);
                    if (estPremier(mpz_get_ui(n))) {
                        gmp_printf("dernier facteur de n : %Zd\n", n);
                        break;
                    }
                }
            }
    
            
            if (mpz_cmp_ui(n, 1) == 0) {
                printf("Tous les facteurs de n ont été trouvés.\n");
                break;
            }

        }
    }

}

int main() {
    mpz_t n;
    mpz_init(n);
    mpz_set_str(n, "108108", 10); // On choisit le nombre à factoriser
    int B = 50;
    int m = 1000;
    int taille = 0;
    mpz_t tableau[m];
    for (int i = 0; i < m; i++) {
        mpz_init(tableau[i]);
    }

    factoriser(tableau, n, B, m, &taille);
    
    return 0;
}