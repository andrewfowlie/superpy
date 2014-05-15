#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include"pmodel.h"
                               
static double rc=34.53; 
static double mz=91.1884;

/*double MLZP(double c, double m)
{ double mr=m/2.405, ac=fabs(c+0.5);
  double exp1=exp(-rc*ac); 
   double r;
  if(c<-0.5) r= 2.*mr*exp1*sqrt(ac*(ac+1)); 
  if(c==-0.5) r=mr*sqrt(2./rc);
 if(c>-0.5&&c<=-0.3) r= mr*2*sqrt(c+0.5);
if(c>-0.3&&c<=-0.2) r= mr/0.1*(2*sqrt(0.2)*(-.2-c)+(c+.3)*M_PI*(0.8)/2.);
  if(c>-0.2) r=mr*M_PI*(1+c)/2.;
return r*1000.;
}  */ 



/*  double MLZP(double c, double m)
{ double mr=m/2.405, ac=fabs(c+0.5);
  double exp1=exp(-rc*ac), eps=1./sqrt(2.*rc);  
   double r;
  if(c<-0.52) r= 2.*mr*exp1*sqrt(ac*(ac+1));
if(c>=-0.52&&c<-0.5) r= mr/0.02*(-2.*exp1*sqrt(ac*(ac+1))*(c+.5)+sqrt(2./rc)*(c+0.52)); 
  if(c==-0.5) r=mr*sqrt(2./rc);
 if(c>-0.5&&c<=-0.47) r= mr/0.03*(2*sqrt(.03)*(c+0.5)-sqrt(2./rc)*(c+0.47) );
if(c>-0.47&&c<=-0.3) r= mr*2*sqrt(c+0.5);
if(c>-0.3&&c<=-0.2) r= mr/0.1*(2*sqrt(0.2)*(-.2-c)+(c+.3)*M_PI*(0.8)/2.);
  if(c>-0.2) r=mr*M_PI*(1+c)/2.;
return r*1000.;
}*/


 double MLZP(double c, double m)
{ double mr=m/2.405, ac=fabs(c+0.5);
  double exp1=exp(-rc*ac), eps=1./sqrt(2.*rc);  
  double r,x1=-0.54,x2=-0.46,x0=-0.5;
  double p0=((c-x1)*(c-x1)*(c-x2)*(c-x2))/((x0-x1)*(x0-x1)*(x0-x2)*(x0-x2));
  double p11= (c-x1)*(c-x0)*(c-x2)*(c-x2)/((x1-x0)*(x1-x2)*(x1-x2));
  double p21= (c-x2)*(c-x0)*(c-x1)*(c-x1)/((x2-x0)*(x2-x1)*(x2-x1));
  double p20a= (c-x0)*(c-x1)*(c-x1) - p21*( (x2-x1)*(x2-x1) +2*(x2-x0)*(x2-x1));
  double p20b= (x2-x0)*(x2-x1)*(x2-x1);
  double p20=p20a/p20b;
  double p10a=(c-x0)*(c-x2)*(c-x2) - p11*( (x1-x2)*(x1-x2) +2*(x1-x0)*(x1-x2)); 
  double p10b= (x1-x0)*(x1-x2)*(x1-x2); 
  double p10=p10a/p10b;  
  double FL=2*sqrt(x1*x1-x0*x0)*exp(rc*(x1-x0));
  double FL1=2*sqrt(x1*x1-0.25)*exp(-rc*(-x1-0.5))*(x1/(x1*x1-0.25) +rc);
  if(c<x1) r= 2.*mr*exp1*sqrt(ac*(ac+1));
  if(c>=x1&&c<x2)
{ r= mr*(sqrt(2./rc)*p0+FL*p10+FL1*p11+2*sqrt(x2+0.5)*p20+1/sqrt(x2+0.5)*p21);
  printf("cnu %8.3f m %8.3f \n", c,r);
}  if(c>x2&&c<=-0.3) r= mr*2*sqrt(c+0.5);
  if(c>-0.3&&c<=-0.2) r= mr/0.1*(2*sqrt(0.2)*(-.2-c)+(c+.3)*M_PI*(0.8)/2.);
  if(c>-0.2) r=mr*M_PI*(1+c)/2.;
  return r*1000.;
}

/* Invisible width of the Z in MeV*/
extern double zinv(double gzp, double m)
{double f;
 if(m<mz/2.)
{ f=gzp*gzp/24/M_PI/mz*(1-4.*m*m/mz/mz)*(mz*mz-2.*m*m)*1000.;
}
else f=0;
return f;
}   

