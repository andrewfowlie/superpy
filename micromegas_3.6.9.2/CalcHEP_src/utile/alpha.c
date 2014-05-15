#include<math.h>
#include<stdio.h>
#include"alpha.h"

double alpha(int nf, int odr, double lambda,  double dscale)
{

    double  d__2, d__4;

    double b0 = 11. -  (2./3.)*nf;
    double b1 = 51. - (19./3.)*nf;
    double b2 = 2857. - (5033./9.)*nf + (325./27.)*nf*nf;

    double rl = 2*log(dscale / lambda);
    double alpha0= 4*M_PI/(b0*rl);

    d__4 = log(rl) - .5;
    d__2= 2*b1/(b0*b0*rl);

    if(odr==1) return alpha0;
    else if(odr==2) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl));
    else if(odr==3) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl)
         + d__2*d__2 *(d__4*d__4 + b2*b0 /(8*b1*b1) - 1.25)  );
    else { fprintf(stderr,"Can not evaluate alpha in so large oder (%d).\n",odr);
           exit(1);
         }          
}

double findLambda(int nf,int odr, double alpha0, double M)
{ double l1=0.1, l2=0.3;
  double l,a,a1,a2;

  while((a1=alpha(nf,odr,l1,M)-alpha0) > 0)  l1*=0.7;

  while((a2=alpha(nf,odr,l2,M)-alpha0) < 0) l2*=1.3;

  do{ l=(l1*a2-l2*a1)/(a2-a1);
      a=alpha(nf,odr, l,M)-alpha0;
      if(fabs(a1)>fabs(a2)){ a1=a;l1=l;} else {a2=a;l2=l;}
    } while (fabs(a) > 0.00001*alpha0);
  return l;
}

int writeAlpha(FILE*f,int nf,int ordr,double lambda,int nfMx,
    double Mc,double Mb,double Mt,int N, double *q)
{
  int    nf3,nf4,nf5,nf6;
  double  L3, L4, L5, L6;
  double al;
  int i,j;
  
  nf5=nf; L5=lambda; 
  if(nf5==5) 
  {  nf4=4; 
     al=alpha(5, ordr,L5, Mb );
     if(ordr==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
     L4=findLambda(4,ordr, al ,Mb);
  }
  else {nf4=nf5;L4=lambda;}
  if(nf4==4) 
  {  nf3=3; 
     al=alpha(4, ordr,L4, Mc );
     if(ordr==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
      L3=findLambda(3,ordr,al ,Mc);
  }
  else {nf3=nf4;L3=lambda;}
  
  if(nfMx>=6) 
  {  nf6=6; 
     al=alpha(5, ordr,L5, Mt);
     if(ordr==3) al*=(1.-(11./72.)*pow(al/M_PI,2.));
     L6=findLambda(6,ordr, al ,Mt);
  }
  else {nf6=nf; L6=lambda;}            
 

fprintf(stderr,"L3=%f L4=%f L5=%f L6=%f\n",L3,L4,L5,L6);
  
  fprintf(f,"\n#Alpha\n");
  for(i=0,j=1;i<N;i++,j++)
  {  double al;
     double Q=q[i];
           if(Q<1.4)  al=alpha(nf3, ordr, L3, Q);
     else  if(Q<4.5)  al=alpha(nf4, ordr, L4, Q);
     else  if(Q<175.) al=alpha(nf5, ordr, L5, Q);
     else             al=alpha(nf6, ordr, L6, Q);
     fprintf(f," %.5E",al);
     if(j==10) {fprintf(f,"\n"); j=0;}
  }
  fprintf(f,"\n");
}
