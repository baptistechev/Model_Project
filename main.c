#include "arit.c"

poly32_t static inline karatsuba(poly32_t p, poly32_t q){
    int T = max(p->length,q->length);
    if (T == 1) return p->coeffs[0]*q->coeffs[0];
    

}

int main(int argc, char** argv){
    if(argc!=2) return 1;
    N = (__uint32_t) atoi(argv[1]);

    __uint32_t e = 4294967290;
    __uint32_t f = 4294967290;

    printf("%u\n",add(e,f));

    poly32_t a = allocate(4);
    poly32_t b = allocate(2);

    __uint32_t ac[] = {5,4294967290,5,2};
    __uint32_t bc[] = {4,4294967290};

    a->coeffs = ac;
    b->coeffs = bc;

    affichage(a);
    affichage(b);

    affichage(addPoly(a,b));
    affichage(subPoly(a,b));
    affichage(prodPoly(a,b));

}