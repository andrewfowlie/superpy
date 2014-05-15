#include<stdio.h>
#include<math.h>

#include"logInt.h"
extern double simpson( double (*func)(double),double a,double b, double  eps);

#define eps 1.E-4


static double cn[3],cd[3],xa,xb;

static double jIntegrand(double y)
{ 
  double r,x;
  
  if(y==0|| y==1) return 0;
  x=xa+(xb-xa)*y*y*(3-2*y);
  r=(cn[0]+x*(cn[1]+x*cn[2]))/(cd[0]+x*(cd[1]+x*cd[2]));
  
  r=log(fabs(r))/x*(xb-xa)*y*6*(1-y);
  return r;
}

static double j5Integrand(double y)
{ 
  double r,x;
  
  if(y==0|| y==1) return 0;
  x=xa+(xb-xa)*y*y*(3-2*y);
  r=(cn[0]+x*(cn[1]+x*cn[2]))/(cd[0]+x*(cd[1]+x*cd[2]));
  
  r=log(fabs(r))/(x-(cd[0]-cn[0])/(cd[2]-cn[2]) )*(xb-xa)*y*6*(1-y);
  return r;
}




static double iIntegrand(double y)
{ 
  double r,x,n,d;
  
  if(y==0) return 0;
  x=y*y;
  n=(cn[0]+x*(cn[1]+x*cn[2]));  if(n==0) return 0;
  d=(cd[0]+x*(cd[1]+x*cd[2]));  if(d==0) return 0;
  
  r=2*log(fabs(n/d))/y;
  return r;
}


static void addzero(double *x, int*N, double x1)
{ int i;

  if(x1<0 || x1>1) return;

  for(i=*N; i>0; i--)
  {
     if(  x1<x[i-1]) x[i]=x[i-1]; else { x[i]=x1; (*N)++; return;} 
  } 
  (*N)++;
  x[0]=x1;
}


static void findzero(double *c,int*N, double*x)
{
    
  if(c[2]==0) 
  { 
    if(c[1]==0) return; 
    addzero(x,N,-c[0]/c[1]);
  } else
  { double x0=-c[1]/2/c[2], d= x0*x0 - c[0]/c[2];
    if(d<0) return;
    d=sqrt(d);
    addzero(x,N, x0+d);
    if(d) addzero(x,N, x0-d); 
  }
}

static double jInt(void)
{
  int i,n0=2;
  double s,x0[6]={0,1,0,0,0,0};

  findzero(cn,&n0,x0); 
  findzero(cd,&n0,x0);


  for(i=1,s=0;i<n0;i++)
  { double s1;
    xa=x0[i-1];xb=x0[i];
    s1=simpson(jIntegrand,0,1,eps);
    s+=s1;  
  }
  return s;
}


double  Jslog5(double An,double Bn,double Ad,double Bd,double C)
{
  cn[2]=An;  cn[1]=Bn; cn[0]=C;
  cd[2]=Ad;  cd[1]=Bd; cd[0]=C;
  return jInt();
}

double Jflog5(double r1,double r2,double r3,double r4,double r5)
{
  int i,n0=2;
  double s,x0[8]={0,1,0,0,0,0,0,0};

  cn[2]=-r1;  cn[1]=r1-r2+r3;  cn[0]=r2;
  cd[2]= r5;  cd[1]=(r3-r4)/2; cd[0]=(r3+r4)/2-r5;
    
  findzero(cn,&n0,x0); 
  findzero(cd,&n0,x0);

  addzero(x0,&n0,(cd[0]-cn[0])/(cd[2]-cn[2]));


  for(i=1,s=0;i<n0;i++)
  { double s1;
    xa=x0[i-1];xb=x0[i];
    if(xa<xb){ s1=simpson(j5Integrand,0,1,eps);s+=s1;}  
  }
  return s;
}

double J1(double a, double b)
{
  cn[2]=4*a; cn[1]=-4*a;cn[0]=b;
  cd[2]=0;   cd[1]=0;   cd[0]=b;
  return jInt();
}


double J2(double a, double b)
{
  cn[2]=-a; cn[1]= a+b-1; cn[0]=1;
  cd[2]= a; cd[1]=-a+b-1; cd[0]=1;
  return jInt();
}


double J3(double a, double b)
{
  cn[2]=-a; cn[1]= a+1-b;  cn[0]=b;
  cd[2]= a; cd[1]=-a+1-b;  cd[0]=b;
  return jInt();
}


/*
	
int main(void)
{
  double a=0.3, b=0.001, c=0.1;
 printf("J1(a,b)=%E\n", J1(a,b));
 printf("I1(a,b)=%E\n", I1(a,b));
 printf("J2(a,b)=%E\n", J2(a,b));
 printf("I2(a,b)=%E\n", I2(a,b));
 printf("J3(a,b)=%E\n", J3(a,b));
 printf("I3(a,b)=%E\n", I3(a,b));

for(;b>1.E-20; b/=2)
 printf("I_43(%e,%e,%e)= %E\n", a,b,c,I_43(a,b,1,c));

}

*/
