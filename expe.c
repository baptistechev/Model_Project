/*
*   Fonctions pour tester les performances des differents algorithmes.
*/

#ifndef EXPE_C
#define EXPE_C

#include <time.h>
#include "poly.c"

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

double static inline timeProd2(poly32_t (*f)(poly32_t,poly32_t), poly32_t a, poly32_t b){
    /* 
    *   Retourne le temps de calcul d'une fonction f calculant le produit de deux polynomes
    */
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    (*f)(a,b);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("Time : %f\n", cpu_time_used);
    return cpu_time_used;
}

// void executionTests(poly32_t (*f1)(poly32_t,poly32_t),poly32_t (*f2)(poly32_t,poly32_t),int maxDeg){
//     /*
//     *    Test des performances selon la valeur du threshold
//     */
//     FILE *tab = fopen("Tab.csv","w");
//     for (int i = 100; i <3000; i=i+100)
//     {   
//         thresh=i;
//         poly32_t a = allocate(20000);
//         poly32_t b = allocate(20000);

//         srand(time(NULL));

//         for(int j=0;j<20000;j++){
//             a->coeffs[i] = rand()%N; 
//             b->coeffs[i] = rand()%N;
//         }

//         //double t1=timeProd2(f1,a,b);
//         double t1=0;
//         double t2=timeProd2(f2,a,b);
//         fprintf(tab,"%d,%lf,%lf,%lf\n",i,t1,t2,t2-t1);
            
//     }

//     fclose(tab);   
// }

void executionTests2(poly32_t (*f1)(poly32_t,poly32_t),poly32_t (*f2)(poly32_t,poly32_t),int maxDeg){
    FILE *tab = fopen("Tab.csv","w");

    for (int i = 1; i <maxDeg; i=i+100)
    {   
        poly32_t a = allocate(i);
        poly32_t b = allocate(i);

        srand(time(NULL));

        for(int j=0;j<i;j++){
            a->coeffs[i] = rand()%N; 
            b->coeffs[i] = rand()%N;
        }

        //double t1=timeProd2(f1,a,b);
        double t1=timeProd2(f1,a,b);
        double t2=timeProd2(f2,a,b);
        fprintf(tab,"%d,%lf,%lf,%lf\n",i,t1,t2,t2-t1);
            
    }

    fclose(tab);
    
}

#endif