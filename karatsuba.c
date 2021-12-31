/*
*   Fonctions intermediaires pour karatsuba
*   Algorithme de karatsuba pour la multiplication de deux polynomes a coefficients dans Z/NZ
*/
#ifndef KARA_C
#define KARA_C

#include "poly.c"

#define Tk 50 //Threshold pour karatsuba, degre minimum pour utiliser Ka plutot que naif

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

    // Si l'un des polynomes est de degre inferieur au parametre Tk
    // on utilise l'algorithme de multiplication naif.
    if(d-1 <= Tk) return prodPoly(p,q); 

    __uint32_t k = (__uint32_t)floor(d/2) + 1;

    // Decoupage des polynomes
    poly32_t* res;
    poly32_t p0,p1,q0,q1;

    res = splitPoly(p,k);
    p0 = res[0];
    p1 = res[1];

    res = splitPoly(q,k);
    q0 = res[0];
    q1 = res[1];

    free(res);

    //Calcul des termes intermediaires par un appel recursif
    poly32_t a = karatsuba(p0,q0);
    poly32_t c = karatsuba(p1,q1);
    poly32_t b = subPoly( subPoly( karatsuba( addPoly(p0,p1), addPoly(q0,q1)) , a) , c);

    deallocate(p0);
    deallocate(p1);
    deallocate(q0);
    deallocate(q1);

    return addPoly(a, addPoly( increaseDegre(b, k), increaseDegre(c, 2*k) ));
}

#endif