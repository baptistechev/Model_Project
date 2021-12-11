#include <stdio.h>
#include <stdlib.h>

#define aritc

#define add(a,b) (__uint32_t)(((__uint64_t)a+(__uint64_t)b)%N)
#define sub(a,b) (a>=b) ? ((a-b)%N) : (a-b+N)
#define prod(a,b) (__uint32_t)(((__uint64_t)a*(__uint64_t)b)%N)

#define max(a,b) a>b ? a : b
#define min(a,b) a<b ? a : b

__uint32_t N;

typedef struct{
    __uint32_t* coeffs;
    int length;
}poly32;

typedef poly32* poly32_t;

poly32_t allocate(long length){
    poly32_t p = malloc(sizeof(poly32));
    p->length = length;
    p->coeffs = malloc(sizeof(__uint32_t)*length);
    return p;
}

poly32_t static inline copyTo(poly32_t p, poly32_t q){
    for(int i=0;i<p->length;i++) 
        q->coeffs[i] = p->coeffs[i];
    return q;
}

int static inline modInverse(__uint32_t a)
{
    int m = N;
    int y = 0, x = 1;
 
    if (m == 1)
        return 0;
    while (a > 1) {
        int q = a / m;
        int t = m;
        m = a % m, a = t;
        t = y;
        y = sub(x,prod(q,y));
        x = t;
    }
    if (x < 0)
        x = add(x,N);
    return x;
}

poly32_t static inline addPoly(poly32_t p, poly32_t q){
    
    poly32_t r = allocate(max(p->length,q->length));
    r = p->length>q->length ? copyTo(p,r) : copyTo(q,r);
    int n = min(p->length,q->length);
    for (int i = 0; i < n; i++) 
        r->coeffs[i] = add(p->coeffs[i],q->coeffs[i]);
    return r;    

}

poly32_t static inline subPoly(poly32_t p, poly32_t q){
    
    poly32_t r = allocate(max(p->length,q->length));
    r = p->length>q->length ? copyTo(p,r) : copyTo(q,r);
    int n = min(p->length,q->length);
    for (int i = 0; i < n; i++) 
        r->coeffs[i] = sub(p->coeffs[i],q->coeffs[i]);
    return r;    

}

poly32_t static inline prodPoly(poly32_t p, poly32_t q){

    poly32_t r = allocate(p->length + q->length - 1);
    for (int i = 0; i < r->length; i++) r->coeffs[i] = 0;

    for (int i = 0; i < p->length; i++)
    {
        for (int j = 0; j < q->length; j++)
        {
            r->coeffs[i+j] = add(r->coeffs[i+j],prod(p->coeffs[i],q->coeffs[j]));
        }    
    }

    return r;

}

poly32_t static inline increaseDegre(poly32_t p, __uint32_t k){
    /*
    *   Augment de k le degre de chaque terme du polynome p.
    */
    poly32_t q = allocate(p->length+k);
    for(int i=0;i<p->length;i++)
        q->coeffs[i+k] = p->coeffs[i];
    return q;
}

poly32_t static inline constantMult(poly32_t p, int k){
    poly32_t q = allocate(p->length);
    for(int i=0;i<p->length;i++) q->coeffs[i] = p->coeffs[i]*k;
    return q;
}

void static inline affichage(poly32_t p){
    printf("%i",p->coeffs[0]);
    for (__uint32_t i = 1; i < p->length; i++)
    {
        printf(" + %uX^%u",p->coeffs[i],i);
    }
    printf("\n");
}




//evolution temps calcul degré
//impact modulo sur les performance (difference sans / avec)