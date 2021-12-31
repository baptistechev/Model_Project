/*
*   Fonctions intermediaires pour toom-cook
*   Algorithme de Toom-Cook (toom-3) pour la multiplication de deux polynomes a coefficients dans Z/NZ
*/

#ifndef TOOM_C
#define TOOM_C

#include "karatsuba.c"

#define Ttc 250 //Threshold pour toom-3, degre minimum pour utiliser Toom-3 plutot que Kara

__uint32_t *vanMonMatrix;

void initVandermonde(){
    /*
    *   Initialise la matrice de vandermonde pour le corps Z/NZ
    */

    vanMonMatrix=malloc((5*5)*sizeof(__uint32_t));
    vanMonMatrix[0]=1;
    for (int i = 1; i < 5; i++)vanMonMatrix[i]=0;
    vanMonMatrix[5]=modInverse(sub(N,2));
    vanMonMatrix[6]=1;
    vanMonMatrix[7]=modInverse(sub(N,3));
    vanMonMatrix[8]=modInverse(sub(N,6%N));
    vanMonMatrix[9]=2;
    vanMonMatrix[10]=sub(N,1);
    vanMonMatrix[11]=modInverse(2);
    vanMonMatrix[12]=modInverse(2);
    vanMonMatrix[13]=0;
    vanMonMatrix[14]=sub(N,1);
    vanMonMatrix[15]=modInverse(2);
    vanMonMatrix[16]=modInverse(sub(N,2));
    vanMonMatrix[17]=modInverse(sub(N,6%N));
    vanMonMatrix[18]=modInverse(6%N);
    vanMonMatrix[19]=sub(N,2);
    for (int i = 20; i < 24; i++)vanMonMatrix[i]=0;
    vanMonMatrix[24]=1;
}

poly32_t* interpol(poly32_t R0,poly32_t R1,poly32_t R2,poly32_t R3,poly32_t R4){
    /*
    *   Retourne la liste de polynome r0...r4 correspondant au resultat de l'interpolation 
    *   des valuation R0...R4 avec la matrice de Vandermonde
    */

    poly32_t* lret=malloc(5*sizeof(poly32_t));
    poly32_t* L=malloc(5*sizeof(poly32_t));
    lret[0]=R0;
    lret[4]=R4;
    for (int i = 1; i < 4; i++)
    {

        L[0]=constantMult(R0,vanMonMatrix[i*5]);
        L[1]=constantMult(R1,vanMonMatrix[i*5+1]);
        L[2]=constantMult(R2,vanMonMatrix[i*5+2]);
        L[3]=constantMult(R3,vanMonMatrix[i*5+3]);
        L[4]=constantMult(R4,vanMonMatrix[i*5+4]);
        
        poly32_t r=allocate(maxDeg(L,5));
        for(int i=0;i<r->length;i++) r->coeffs[i]=0;   
        for(int i =0;i<5;i++) r=addPoly(r,L[i]);
        lret[i]=r;

    }

    return lret;
}

poly32_t* splitPoly3(poly32_t p, __uint32_t k){
    /*
    *   Renvoit un pointeur sur un tableau contenant 3 polynomes de tailles k
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

   int d=max(p->length,q->length);

    // Si les polynomes sont de degre inferieur au parametre Ttc
    // on utilise l'algorithme de multiplication de karatsuba.
    if(d<Ttc){
        return karatsuba(p,q);
    }

    __uint32_t k = max(floor(p->length/3),floor(q->length/3))+1;  
    
    //Decoupage des polynomes
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

    //Calcul des termes intermediaires
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

    //Interpolation
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

#endif