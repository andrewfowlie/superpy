/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <math.h>
#include<stdlib.h>
#include"tools.h"
#include"simpson.h"

double divy_(double y)
{
    if (fabs(y) > 0.01) return  (1 - exp(-y)) / y;
    else   return 1 - y * (1 - y * (1 - y * (1 - y/5) / 4) / 3 ) / 2;
} /* divy_ */



double dinter_(double x, int n, double *xi, double *yi)
{
    int i;
    double x1, x2, x3, x4;

    /* Parameter adjustments */
    --yi;
    --xi;

    /* Function Body */
    if (x <= xi[2]) 
    {	x1 = xi[1];
	x2 = xi[2];
	x3 = xi[3];
	return  yi[1] * (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3)) + 
		yi[2] * (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3)) + 
		yi[3] * (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
    } else if (x > xi[n - 1]) 
    {	x1 = xi[n - 2];
	x2 = xi[n - 1];
	x3 = xi[n];
	return  yi[n - 2] * (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3)) +
		yi[n - 1] * (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3)) +
		yi[n    ] * (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
    } else 
    {	for(i = 2; x > xi[i]; i++) ;
	x1 = xi[i - 1];
	x2 = xi[i];
	x3 = xi[i + 1];
	x4 = xi[i + 2];
	return   yi[i-1]*(x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4)) + 
	         yi[i  ]*(x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4)) + 
	         yi[i+1]*(x-x1)*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4)) + 
	         yi[i+2]*(x-x1)*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3));
    }
} /* dinter_ */

static int gamma0_param;

static double gamma0_(double z)
{  return pow(z, gamma0_param-1.) * exp(-(z)); }


double gammai_(int n, double a)
{
    static double aa = -1.;

    double ret_val;

    int i;
    static double ea, eps, err, gmem[100];
    static int iscalk[100];

    if (a != aa) 
    {
	aa = a;
	ea = exp(-(a));
	for (i = 0; i < 100; ++i)    iscalk[i] = 0;
    }
    if (n <= 100 && iscalk[n - 1]) 
    {	ret_val = gmem[n - 1];
	return ret_val;
    }
    ea = exp(-(a));
    ret_val = 1 - ea;
    iscalk[0] = 1;
    gmem[0] = ret_val;
    err = 0.;
    for (i = 1; i < n; ++i) 
    {	ea = a * ea;
	ret_val = i * ret_val - ea;
	iscalk[i] = 1;
	gmem[i] = ret_val;
	err = i * err + ea * 1e-14f;
	if (err > ret_val * 1e-8) 
	{
	    gamma0_param = n;
	    eps = 1e-8;
	    ret_val = gauss345(gamma0_, 0., a, eps,NULL);
	    if (n <= 100){ iscalk[n - 1] = 1; gmem[n - 1] = ret_val;}
	    return ret_val;
	}
    }
    return ret_val;
} /* gammai_ */


static double (*conv_f1) (double);
static double (*conv_f2) (double);
static double conv_x,conv_b1, conv_b2,conv_c1, conv_c2;

static double conv_f(double x)
{
double z1,z2;
  if(x==0.) { z1=0.; z2=0.;}
  else {z1= 0.5*pow(x,1/conv_b1), z2=0.5*pow(x,1/conv_b2);}
  return  
  (*conv_f1)(conv_x*z1)*(*conv_f2)(conv_x*(1-z1))*pow(1-z1,conv_b2-1)*conv_c1 
 +(*conv_f2)(conv_x*z2)*(*conv_f1)(conv_x*(1-z2))*pow(1-z2,conv_b1-1)*conv_c2 ;
 }

double convol_(double (*f1)(double ),double (*f2)(double ), 
          double b1, double b2,  double x, double eps)
{
/*
convol_(x)=x^(1-b1-b2)*Int(from 0 to x, f1(z)*z^(b1-1)*f2(x-z)*(x-z)^(b2-1))
for the case x==0, conv=f1(0)*f2(0)*gamma_(b1)*gamma_(b2)/gamma_(b1+b2) 
*/
  conv_f1=f1;
  conv_f2=f2;
  conv_b1=b1;
  conv_b2=b2;
  conv_c1=pow(0.5,b1)/b1;
  conv_c2=pow(0.5,b2)/b2;
  conv_x=x;
 
  return   gauss345(&(conv_f),0.,1.,eps,NULL); 
  

 } /* convol_ */
