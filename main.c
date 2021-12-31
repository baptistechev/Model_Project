#ifndef MAIN_C
#define MAIN_C

#include "poly.c"
#include "karatsuba.c"
#include "toom.c"
#include "expe.c"

int main(int argc, char** argv){

    if(argc!=2) return 1;
    //Initialisation du corps modulo N avec l'entier passe en argument.
    N = (__uint32_t) atoi(argv[1]);

    //Initialisation de la matrice de vandermonde dans Z/NZ
    initVandermonde();

    /*** Test des algorithmes de multiplications ***/
    poly32_t a = allocate(10);
    poly32_t b = allocate(9);
    __uint32_t ac[] = {5%N,1%N,2%N,0,0,8%N,0,2%N,0,3%N};
    __uint32_t bc[] = {1%N,3%N,0,2%N,0,8%N,0,0,7%N};
    a->coeffs = ac;
    b->coeffs = bc;

    printf("Multiplication des polynomes suivants :\n");
    affichage(a);
    affichage(b);

    printf("\nMultiplication naive : ");
    affichage(prodPoly(a,b));

    printf("\nMultiplication de karatsuba : ");
    affichage(karatsuba(a,b));

    printf("\nMultiplication de toom-cook : ");
    affichage(toom3(a,b));

    /*** Test des performances ***/
    int size = 15000;

    poly32_t c = allocate(size);
    poly32_t d = allocate(size);

    srand(time(NULL));

    for(int i=0;i<size;i++){
        c->coeffs[i] = rand()%N; 
        d->coeffs[i] = rand()%N;
    }

    printf("\nMultiplication naive : ");
    timeProd(prodPoly,c,d);

    printf("\nMultiplication de karatsuba : ");
    timeProd(karatsuba,c,d);

    printf("\nMultiplication de toom-cook : ");
    timeProd(toom3,c,d);

    /*
    *   !! Se referer au rapport pour plus de details sur les performances !!
    */
    printf("\n Se referer au rapport pour plus de details sur les performances\n");
}

#endif