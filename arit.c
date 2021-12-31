/*
* Operation de bases pour l'arithmetique modulaire.
*/
#ifndef ARIT_C
#define ARIT_C

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/*-------------Macros des opérations de base------------------*/
#define add(a,b) (__uint32_t)(((__uint64_t)a+(__uint64_t)b)%N)
#define sub(a,b) (a>=b) ? ((a-b)%N) : (a-b+N)
#define prod(a,b) (__uint32_t)(((__uint64_t)a*(__uint64_t)b)%N)
#define max(a,b) a>b ? a : b
#define min(a,b) a<b ? a : b

/*--Variables globales--*/
__uint32_t N;

__uint32_t static inline modInverse(__uint32_t a){
    /*
    * Calcul l'inverse de a dans Z/NZ
    */
    __uint32_t m = N;
    __uint32_t y = 0, x = 1;
 
    if (m == 1)
        return 0;

    while (a > 1) {
        __uint32_t q = a / m;
        __uint32_t t = m;
        m = a % m, a = t;
        t = y;
        y = sub(x,prod(q,y));
        x = t;
    }
    if (x < 0)
        x = add(x,N);
    return x;
}

#endif