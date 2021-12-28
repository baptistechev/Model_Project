#include "arit.c"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

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

    free(res);

    poly32_t a = karatsuba(p0,q0);
    poly32_t c = karatsuba(p1,q1);
    poly32_t b = subPoly( subPoly( karatsuba( addPoly(p0,p1), addPoly(q0,q1)) , a) , c);

    deallocate(p0);
    deallocate(p1);
    deallocate(q0);
    deallocate(q1);

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

poly32_t static inline toom3(poly32_t p, poly32_t q){
    /*
    *   Performe la multiplication de deux polynomes en utilisant l'algorithme de toom-cook 
    */

    __uint32_t k = max(floor(p->length/3),floor(q->length/3))+1;
    // if(d-1 <= T) return prodPoly(p,q); 

    // __uint32_t k = (__uint32_t)floor(d/2) + 1;
    
    if(k<3){
        return prodPoly(p,q);
    }

    poly32_t *res;
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

    free(res);

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


    r_v0 = toom3(p_v0,q_v0);
    r_v1 = toom3(p_v1,q_v1);
    r_v2 = toom3(p_v2,q_v2);
    r_v3 = toom3(p_v3,q_v3);
    r_v4 = toom3(p_v4,q_v4);

    deallocate(p0);
    deallocate(p1);
    deallocate(p2);
    deallocate(q0);
    deallocate(q1);
    deallocate(q2);

    //TEST temps d'execution interpol
   
    poly32_t* lInter=interpol(r_v0,r_v1,r_v2,r_v3,r_v4);

    poly32_t ret=allocate(p->length+q->length-1);
    for(int i=0;i<ret->length;i++) ret->coeffs[i]=0;

    for(int i=0;i<5;i++){
        ret=addPoly(ret,increaseDegre(lInter[i],i*k));
    }

    for(int i=0;i<5;i++) deallocate(lInter[i]);
    free(lInter);

    return ret;
}

int main(int argc, char** argv){

    if(argc!=2) return 1;
    N = (__uint32_t) atoi(argv[1]);

    //Initialisation of vandermonde matrix in current prime field
    initVandermonde();

    int size = 15000;

    poly32_t a = allocate(size);
    poly32_t b = allocate(size);

    srand(time(NULL));

    for(int i=0;i<size;i++){
        a->coeffs[i] = rand()%N; 
        b->coeffs[i] = rand()%N;
    }

    // __uint32_t ac[] = {5%N,1%N,2%N,0,0,8%N,0,2%N,0,3%N};
    // __uint32_t bc[] = {1%N,3%N,0,2%N,0,8%N,0,0,7%N};

    //a->coeffs = ac;
    //b->coeffs = bc;
    
    timeI=0;
    timeProd(prodPoly,a,b);
    timeProd(karatsuba,a,b);
    timeI=0;
    timeProd(toom3,a,b);
    printf("timeInterpol:%lf\n",timeI);

    deallocate(a);
    deallocate(b);
}