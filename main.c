#include "arit.c"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define T 100 //Threshold pour karatsuba, degre minimum pour utiliser Ka plutot que naif

poly32_t static inline timeProd(poly32_t (*f)(poly32_t,poly32_t), poly32_t a, poly32_t b){
    /* 
    *   Affiche le temps de calcul d'une fonction f calculant le produit de deux polynomes et renvoit la sortie de f
    */
    clock_t start, end;
    double cpu_time_used;
    poly32_t res;

    start = clock();
        res = (*f)(a,b);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time : %f\n", cpu_time_used);
    return res;
}

poly32_t* splitPoly(poly32_t p, __uint32_t k){
    /*
    *   Renvoit un pointeur sur un tableau contenant 2 polynomes de tailles respectives k et deg(p)-k+1
    */
    poly32_t p0 = allocate(k);
    poly32_t p1 = allocate(p->length-k);
    for(__uint32_t i=0;i<p0->length;i++) p0->coeffs[i] = p->coeffs[i];
    for(__uint32_t i=0;i<p1->length;i++) p1->coeffs[i] = p->coeffs[i+k];
    poly32_t* res = malloc(2*sizeof(poly32_t)); 
    res[0] = p0;
    res[1] = p1;
    return res;
}

poly32_t static inline karatsuba(poly32_t p, poly32_t q){
    /*
    *   Performe la multiplication de deux polynomes en utilisant l'algorithme de karatsuba 
    */

    __uint32_t d = min(p->length,q->length);

    // Si l'un des polynomes est de degre inferieur au parametre T
    // on utilise l'algorithme de multiplication naif.
    if(d-1 <= T) return prodPoly(p,q); 

    __uint32_t k = (__uint32_t)floor(d/2) + 1;

    poly32_t* res;
    poly32_t p0,p1,q0,q1;

    res = splitPoly(p,k);
    p0 = res[0];
    p1 = res[1];

    res = splitPoly(q,k);
    q0 = res[0];
    q1 = res[1];

    poly32_t a = karatsuba(p0,q0);
    poly32_t c = karatsuba(p1,q1);
    poly32_t b = subPoly( subPoly( karatsuba( addPoly(p0,p1), addPoly(q0,q1)) , a) , c);

    return addPoly(a, addPoly( increaseDegre(b, k), increaseDegre(c, 2*k) ));
}

int main(int argc, char** argv){
    if(argc!=2) return 1;
    N = (__uint32_t) atoi(argv[1]);

    // __uint32_t e = 4294967290;
    // __uint32_t f = 4294967290;

    // printf("%u\n",add(e,f));

    int size = 1000;

    poly32_t a = allocate(size);
    poly32_t b = allocate(size);

    srand(time(NULL));

    for(int i=0;i<size;i++){
        a->coeffs[i] = rand()%1000; 
        b->coeffs[i] = rand()%1000;
    }

    // __uint32_t ac[] = {2,4,3,2,1};
    // __uint32_t bc[] = {2,3,1,5,8,7,5,9};

    // a->coeffs = ac;
    // b->coeffs = bc;

    timeProd(prodPoly,a,b);
    timeProd(karatsuba,a,b);
}