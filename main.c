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

poly32_t* splitPoly3(poly32_t p, __uint32_t k){
    /*
    *   Renvoit un pointeur sur un tableau contenant 3 polynomes de tailles respectives k
    */
    poly32_t p0 = allocate(k);
    poly32_t p1 = allocate(k);
    poly32_t p2 = allocate(k);

    for(__uint32_t i=0;i<k;i++){
        p0->coeffs[i] = 0;
        p1->coeffs[i] = 0;
        p2->coeffs[i] = 0;
    } 

    for(__uint32_t i=0;i<k;i++) p0->coeffs[i] = p->coeffs[i];
    for(__uint32_t i=0;i<k;i++) p1->coeffs[i] = p->coeffs[i+k];
    for(__uint32_t i=0;i<p->length-2*k;i++) p2->coeffs[i] = p->coeffs[i+2*k];

    poly32_t* res = malloc(3*sizeof(poly32_t)); 
    res[0] = p0;
    res[1] = p1;
    res[2] = p2;
    return res;
}

int maxDeg (poly32_t* L,int n){
    int max=0;
    for(int i=0;i<n;i++){
        if(L[i]->length>max) max=L[i]->length;
    }
    return max;
}

poly32_t* interpol(poly32_t R0,poly32_t R1,poly32_t R2,poly32_t R3,poly32_t R4){

    poly32_t r0=R0;

    poly32_t*L1=malloc(5*sizeof(poly32_t));
    L1[0]=constantMult(R0,-1/2);L1[1]=constantMult(R1,1);L1[2]=constantMult(R2,-1/3);L1[3]=constantMult(R3,-1/6);L1[4]=constantMult(R4,2);
    poly32_t r1=allocate(maxDeg(L1,5));
    for(int i=0;i<r1->length;i++) r1->coeffs[i]=0;
    for(int i =0;i<5;i++) r1=addPoly(r1,L1[i]);

    poly32_t*L2=malloc(4*sizeof(poly32_t));
    L2[0]=constantMult(R0,-1);L2[1]=constantMult(R1,1/2);L2[2]=constantMult(R2,1/2);L2[3]=constantMult(R3,-1);
    poly32_t r2=allocate(maxDeg(L2,4));
    for(int i=0;i<r2->length;i++) r2->coeffs[i]=0;
    for(int i =0;i<4;i++) r2=addPoly(r2,L2[i]);

    poly32_t*L3=malloc(5*sizeof(poly32_t));
    L3[0]=constantMult(R0,1/2);L3[1]=constantMult(R1,-1/2);L3[2]=constantMult(R2,-1/6);L3[3]=constantMult(R3,1/6);L3[4]=constantMult(R4,-2);
    poly32_t r3=allocate(maxDeg(L3,5));
    for(int i=0;i<r3->length;i++) r3->coeffs[i]=0;
    for(int i =0;i<5;i++) r1=addPoly(r3,L3[i]);

    poly32_t r4=R4;

    poly32_t* lret=malloc(5*sizeof(poly32_t));
    lret[0]=r0;lret[1]=r1;lret[2]=r2;lret[3]=r3;lret[4]=r4;

    return lret;

}


poly32_t static inline toom3(poly32_t p, poly32_t q){
    /*
    *   Performe la multiplication de deux polynomes en utilisant l'algorithme de karatsuba 
    */

    __uint32_t k = max(floor(p->length/3),floor(q->length/3))+1;

    // if(d-1 <= T) return prodPoly(p,q); 

    // __uint32_t k = (__uint32_t)floor(d/2) + 1;

    poly32_t* res;
    poly32_t p0,p1,p2,q0,q1,q2;
    poly32_t p_v0,p_v1,p_v2,p_v3,p_v4,q_v0,q_v1,q_v2,q_v3,q_v4;
    poly32_t r_v0,r_v1,r_v2,r_v3,r_v4;

    res = splitPoly3(p,k);
    p0 = res[0];
    p1 = res[1];
    p2 = res[2];

    res = splitPoly3(q,k);
    q0 = res[0];
    q1 = res[1];
    q2 = res[2];

    p_v0 = p0;
    p_v1 = addPoly(p2,addPoly(p1,p0));
    p_v2 = addPoly(subPoly(p2,p1),p0);
    p_v3 = addPoly( constantMult(p2,4), addPoly( constantMult(p1,2), p0));
    p_v4 = p2;

    q_v0 = q0;
    q_v1 = addPoly(q2,addPoly(q1,q0));
    q_v2 = addPoly(subPoly(q2,q1),q0);
    q_v3 = addPoly( constantMult(q2,4), addPoly( constantMult(q1,2), q0));
    q_v4 = q2;

    r_v0 = prodPoly(p_v0,q_v0);
    r_v1 = prodPoly(p_v1,q_v1);
    r_v2 = prodPoly(p_v2,q_v2);
    r_v3 = prodPoly(p_v3,q_v3);
    r_v4 = prodPoly(p_v4,q_v4);

    affichage(r_v0);
    affichage(r_v1);
    affichage(r_v2);
    affichage(r_v3);
    affichage(r_v4);

    poly32_t* lInter=interpol(r_v0,r_v1,r_v2,r_v3,r_v4);

    printf("\ninterpolation\n");
    for(int i=0;i<5;i++){
        affichage(lInter[i]);
    }

    return p;
}


int main(int argc, char** argv){
    if(argc!=2) return 1;
    N = (__uint32_t) atoi(argv[1]);

    // __uint32_t e = 4294967290;
    // __uint32_t f = 4294967290;

    // printf("%u\n",add(e,f));

    int size = 1000;

    poly32_t a = allocate(10);
    poly32_t b = allocate(9);

    // srand(time(NULL));

    // for(int i=0;i<size;i++){
    //     a->coeffs[i] = rand()%1000; 
    //     b->coeffs[i] = rand()%1000;
    // }

    __uint32_t ac[] = {5,1,2,0,0,7,0,2,0,3};
    __uint32_t bc[] = {1,3,0,2,0,8,0,0,7};

    a->coeffs = ac;
    b->coeffs = bc;

    // timeProd(prodPoly,a,b);
    // timeProd(karatsuba,a,b);

    toom3(a,b);

}