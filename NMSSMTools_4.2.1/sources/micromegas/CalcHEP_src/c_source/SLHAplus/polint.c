#include "SLHAplus.h"


static double  polintN(double x, int n,  double *xa, double *ya)
{  double z[20];
   int i,m;
   for(i=0;i<n;i++) z[i]=ya[i];
   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}

static int  leftXN(int n,int dim,  double * xa, double x)
{  int k1,k2,k3;

   k1=n/2;                         
   k2=dim-(n+1)/2-1;
   
   if(xa[0]< xa[dim-1])
   {              
     if(x<=xa[k1]) return 0;
     if(x>=xa[k2]) return dim-n;                   
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]>x)k2=k3; else k1=k3;
     }
   } else 
   {  
     if(x>=xa[k1]) return 0;
     if(x<=xa[k2]) return dim-n;
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]<x)k2=k3; else k1=k3;
     }
   }
   return k1+1-n/2;
}

static int  leftX(int dim, double * xa, double x)
{  int k1,k2,k3;
                                                                                
   if(x<=xa[0]) return 0;
   if(x>=xa[dim-3]) return dim-3;
                                                                                
   k1=0;
   k2=dim-3;
                                                                                
   while(k2-k1>1)
   { k3=(k1+k2)/2;
     if(xa[k3]>x)k2=k3; else k1=k3;
   }
   return k1;
}


double polint2(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(2,n, xa, x);
   return polintN(x,2,xa+shift, ya+shift);
}


double polint3(double x, int n,  double *xa, double *ya)
{ int shift=leftX(n, xa, x);
  double ar;
  ar=polintN(x,3,xa+shift, ya+shift);
  if(shift==0) return ar;
  shift--;
  return 0.5*( ar+ polintN(x,3,xa+shift, ya+shift));
}


double polint4(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(4,n, xa, x);
   return polintN(x,4,xa+shift, ya+shift);
}
