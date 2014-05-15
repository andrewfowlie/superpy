#include<stdio.h>
#include "my_complex.h"

const my_complex icomp ={0,1};

void affichec(my_complex a)
{
printf("%e +i %e \n",a.r,a.i); 
}

my_complex prod(my_complex a, my_complex b)
{
my_complex temp;
temp.r=a.r*b.r-a.i*b.i;
temp.i=a.r*b.i+b.r*a.i;
return temp;
}

my_complex conjug(my_complex a)
{
my_complex temp;
temp.r=a.r;
temp.i=-a.i;
return temp;
}

my_complex prodscal(double a,my_complex b)

{
my_complex temp;
temp.r=a*b.r;
temp.i=a*b.i;
return temp;
}

my_complex somme(my_complex a,my_complex b)
{
my_complex temp;
temp.r=a.r+b.r;
temp.i=a.i+b.i;
return temp;
}

my_complex diff(my_complex a,my_complex b)
{
my_complex temp;
temp.r=a.r-b.r;
temp.i=a.i-b.i;
return temp;
}
