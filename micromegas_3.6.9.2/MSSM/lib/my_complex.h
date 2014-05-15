
#ifndef  __MY_COMPLEX__
#define  __MY_COMPLEX__

typedef struct{ double r;  double i; } my_complex;

const my_complex icomp;

extern void affichec(my_complex a);

extern my_complex prod(my_complex a,my_complex b);

extern my_complex conjug(my_complex a);

extern my_complex prodscal(double a, my_complex b);

extern my_complex diff(my_complex a, my_complex b);

extern my_complex somme(my_complex a, my_complex b);

#endif

