#include <stdio.h>
#include <stdlib.h>

#define aritc

/*-------------Macros des opérations de base------------------*/
#define add(a,b) (__uint32_t)(((__uint64_t)a+(__uint64_t)b)%N)
#define sub(a,b) (a>=b) ? ((a-b)%N) : (a-b+N)
#define prod(a,b) (__uint32_t)(((__uint64_t)a*(__uint64_t)b)%N)
#define max(a,b) a>b ? a : b
#define min(a,b) a<b ? a : b
/*-------------------------------------------------------------*/

/*--Variables globales--*/
__uint32_t N;
__uint32_t *vanMonMatrix;
/*----------------------*/

/*--Structure du polynome--*/
typedef struct{
    __uint32_t* coeffs;
    int length;
}poly32;
typedef poly32* poly32_t;
/*-------------------------*/

/*-------------Methodes de gestion memoire-------------*/
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
/*----------------------------------------------------*/

/*---------------Opérations sur les polynomes---------------*/
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
    /*
    *   Multiplie les coefficients du polynome p par k.
    */
    poly32_t q = allocate(p->length);
    for(int i=0;i<p->length;i++) q->coeffs[i] = prod(p->coeffs[i],k);
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
/*---------------------------------------------------*/

__uint32_t static inline modInverse(__uint32_t a)
{
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

/*__uint32_t static inline modInverse(__uint32_t a)
{
    int r1=a,r2=N,u1=0,v1=1,u2=1,v2=0;
    while(r2!=0)
    {
        int q=r1/r2;
        int r3=r1,u3=u1,v3=v1;
        r1 = r2, u1 = u2, v1 = v2;
        r2 = sub(r3,prod(q,r2)), u2 = sub(u3,prod(q,u2)), v2 = sub(v3,prod(q,v2));
    }

    return v1;
    
}*/

void initVandermonde(){
    /*
    *   Initialise la matrice de vandermonde pour le corp modulo N
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

// poly32_t* interpol(poly32_t R0,poly32_t R1,poly32_t R2,poly32_t R3,poly32_t R4){

//     /*
//     *   Retourne la liste de polynome r0...r4 correspondant au resultat de l'interpolation des valuation R0...R4 avec la matrice de Vandermonde
//     */
//     poly32_t r0=R0;

//     poly32_t*L1=malloc(5*sizeof(poly32_t));
//     L1[0]=constantMult(R0,modInverse(sub(N,2)));L1[1]=R1;L1[2]=constantMult(R2,modInverse(sub(N,3)));L1[3]=constantMult(R3,modInverse(sub(N,6%N)));L1[4]=constantMult(R4,2);
//     poly32_t r1=allocate(maxDeg(L1,5));
//     for(int i=0;i<r1->length;i++) r1->coeffs[i]=0;
//     for(int i =0;i<5;i++) r1=addPoly(r1,L1[i]);

//     poly32_t*L2=malloc(5*sizeof(poly32_t));
//     L2[0]=constantMult(R0,sub(N,1));L2[1]=constantMult(R1,modInverse(2));L2[2]=constantMult(R2,modInverse(2));L2[3]=constantMult(R3,0);L2[4]=constantMult(R4,sub(N,1));
//     poly32_t r2=allocate(maxDeg(L2,5));
//     for(int i=0;i<r2->length;i++) r2->coeffs[i]=0;
//     for(int i =0;i<5;i++) r2=addPoly(r2,L2[i]);

//     poly32_t*L3=malloc(5*sizeof(poly32_t));
//     L3[0]=constantMult(R0,modInverse(2));L3[1]=constantMult(R1,modInverse(sub(N,2)));L3[2]=constantMult(R2,modInverse(sub(N,6%N)));L3[3]=constantMult(R3,modInverse(6%N));L3[4]=constantMult(R4,sub(N,2));
//     poly32_t r3=allocate(maxDeg(L3,5));
//     for(int i=0;i<r3->length;i++) r3->coeffs[i]=0;
//     for(int i =0;i<5;i++) r3=addPoly(r3,L3[i]);

//     poly32_t r4=R4;

//     poly32_t* lret=malloc(5*sizeof(poly32_t));
//     lret[0]=r0;lret[1]=r1;lret[2]=r2;lret[3]=r3;lret[4]=r4;

//     return lret;

// }

//evolution temps calcul degré
//impact modulo sur les performance (difference sans / avec)