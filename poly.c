/*
*  Definition des polynomes
*  Operation de bases pour les polynomes
*/
#ifndef POLY_C
#define POLY_C

#include "arit.c"

/*--Structure du polynome--*/
typedef struct{
    __uint32_t* coeffs;
    int length;
}poly32;
typedef poly32* poly32_t;


/*-------------Methodes de gestion memoire-------------*/
poly32_t allocate(long length){
    /*
    *   cree un polynome de taille length
    */
    poly32_t p = malloc(sizeof(poly32));
    p->length = length;
    p->coeffs = malloc(sizeof(__uint32_t)*length);
    return p;
}

void deallocate(poly32_t p){
    /*
    * Libere la memoire occupe par le polynome p
    */
    free(p->coeffs);
    free(p);
}

poly32_t static inline copyTo(poly32_t p, poly32_t q){
    /*
    *   Copie un polynome p dans un polynome q
    */
    for(int i=0;i<p->length;i++) 
        q->coeffs[i] = p->coeffs[i];
    return q;
}

/*---------------Opérations sur les polynomes---------------*/
poly32_t static inline addPoly(poly32_t p, poly32_t q){
    /*
    *   Additionne deux polynomes a coefficients dans Z/NZ
    */
    
    poly32_t r = allocate(max(p->length,q->length));
    r = p->length>q->length ? copyTo(p,r) : copyTo(q,r);
    int n = min(p->length,q->length);
    for (int i = 0; i < n; i++) 
        r->coeffs[i] = add(p->coeffs[i],q->coeffs[i]);
    return r;    

}

poly32_t static inline subPoly(poly32_t p, poly32_t q){
    /*
    *   Soustrait deux polynomes a coefficients dans Z/NZ
    */
    
    poly32_t r = allocate(max(p->length,q->length));
    r = p->length>q->length ? copyTo(p,r) : copyTo(q,r);
    int n = min(p->length,q->length);
    for (int i = 0; i < n; i++) 
        r->coeffs[i] = sub(p->coeffs[i],q->coeffs[i]);
    return r;    

}

poly32_t static inline prodPoly(poly32_t p, poly32_t q){
    /*
    *   Multiplie deux polynomes a coefficients dans Z/NZ
    */

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
    *   Augmente de k le degre de chaque terme du polynome p.
    */
    poly32_t q = allocate(p->length+k);
    for(int i=0;i<p->length;i++)
        q->coeffs[i+k] = p->coeffs[i];
    return q;
}

poly32_t static inline constantMult(poly32_t p, int k){
    /*
    *   Multiplie les coefficients du polynome p par k.
    */
    poly32_t q = allocate(p->length);
    for(int i=0;i<p->length;i++) q->coeffs[i] = prod(p->coeffs[i],k);
    return q;
}

void static inline affichage(poly32_t p){
    /*
    * Affiche le polynome p
    */
    printf("%i",p->coeffs[0]);
    for (__uint32_t i = 1; i < p->length; i++)
    {
        printf(" + %uX^%u",p->coeffs[i],i);
    }
    printf("\n");
}

int maxDeg (poly32_t* L,int n){
    /*
    *   Renvoit le degré maximum d'une liste de polynôme
    */
    int max=0;
    for(int i=0;i<n;i++){
        if(L[i]->length>max) max=L[i]->length;
    }
    return max;
}

#endif