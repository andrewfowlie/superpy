#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include"pmodel.h"


static double ee;
static double sw2;
static double mt;
static double v;
static double g10;
static double mzp;

static double mm=2.405; /* ???????????*/


static double del=1.E-6;     /*to  avoid division by zero */
static double rc=34.53; 
/* problems with using SAVE as below, I do not know why */                             
static double cc,cd;
static double xx,xy;

double SAVE(double EE, double sw, double Mtop, double V, double G10, double MZP, double RC)
{ ee=EE; sw2=sw*sw; mt=Mtop/1000.; v=V/1000,g10=G10; mzp=MZP/1000; return 1;}

static double ftr(double z, double c)
{ double exp1=exp(rc), exp2=exp(rc*(1-2*c)), exp3=exp(rc*(1-2*(c+del))); 
  double  f=pow(z*exp1,2-c);
  if(c==1/2) return f*sqrt(exp1*(1-2*(c+del))/(exp3-1));
  else return f*sqrt(exp1*(1-2*c)/(exp2-1));
}   


static double bLKP(double c, double m)
{ double ac=fabs(c+0.5);
  double x=m*exp(-rc); 
  double jnu,ynu;
if(ac<0.01)
{jnu=j0(x);
 ynu=y0(x);
}
else{
  bessjy(x,ac,&jnu,&ynu);   
}
  return -jnu/ynu;
}   

static double NLKP(double c, double m)
{ double ac=fabs(c+0.5);
  double x=m*exp(-rc), acp; 
  double jnu,ynu,jnu1,ynu1,jnu2,ynu2;
  double f;
if(ac<.01)
{
jnu=j0(m);
 ynu=y0(m);
}
else
{bessjy(m,ac,&jnu,&ynu);
}
 if(c<-0.5){acp=ac+1;
 bessjy(x,acp,&jnu1,&ynu1);
}
else{acp=ac-1;}
if(acp<0.)
{
   if(ac<.01)
    {jnu1=-j1(x);
      ynu1=-y1(x);
     }
    else
    {
bessjy(x,-acp,&jnu2,&ynu2);
jnu1=cos(-acp*M_PI)*jnu2-sin(-acp*M_PI)*ynu2;
ynu1=sin(-acp*M_PI)*jnu2+cos(-acp*M_PI)*ynu2;
   }

}
 
 f= sqrt( exp(2*rc)/(2*rc)*( pow(jnu+bLKP(c,m)*ynu,2)-exp(-2*rc)*
pow(jnu1+bLKP(c,m)*ynu1,2)   ) );
  return f;
} 

static double fLKP(double z, double c, double m)
{ double ac=fabs(c+0.5);
  double jnu,ynu;
  double f;
if(ac<.01)
{
jnu=j0(m*z);
 ynu=y0(m*z);
}
else
{bessjy(m*z,ac,&jnu,&ynu);
 }
 f=pow(z*exp(rc),5/2.)*exp(rc/2)/sqrt(rc)/NLKP(c,m)*(jnu+bLKP(c,m)*ynu);

  return f;
} 

static double bgauge(double x)
{   double z=x*exp(-rc); 
  return -j1(z)/y1(z);
}

static double Ngauge(double x)
{ double z=x*exp(-rc); 
  return sqrt(0.5*(pow((j1(x)+bgauge(x)*y1(x)),2)-exp(-2*rc)*
  pow((j0(z)+bgauge(x)*y0(z)),2)));
} 

static double fgauge(double z, double x)
{ return exp(rc/2)*z/Ngauge(x)*(j1(z*x)+bgauge(x)*y1(z*x));
} 

static double rg1fac(double z)
{ double f=fgauge(z,mm)*ftr(z,cd)*fLKP(z,cc,xx)*pow(z*exp(rc),-4);
  return f;
} 

static double rgzpfac(double z)
{  double f=fgauge(z,mm)*fLKP(z,cc,xx)*fLKP(z,cc,xx)*pow(z*exp(rc),-4);
     return f;
}

static double rgzpfac2(double z)
{ double  xx1=MLZP(cc,mm)/1000.;
  double  xx2=MLZP(cd,mm)/1000.;
  double f=fgauge(z,mm)*fLKP(z,cc,xx1)*fLKP(z,cd,xx2)*pow(z*exp(rc),-4);
  return f;
}


static double rgzptfac(double z)
{ double f;
  f=fgauge(z,mm)*ftr(z,cc)*ftr(z,cc)*pow(z*exp(rc),-4);
  return f;
} 


double g1f(double c1,double c2)
{ double  low=exp(-rc),  eps=1.e-4;
  double f;
  cc=c1;
  cd=c2;
  xx=MLZP(cc,mm)/1000.;
  f=simps(rg1fac,low,1,eps)*exp(-rc/2);
   /*printf(" g1f %e rc %e\n",f,rc);
  */return f;
}

double g2f(double c1,double c2)
{ double  low=exp(-rc),  eps=1.e-4;
  double f;
  cc=c1;
  cd=c2;
  f=simps(rgzpfac2,low,1,eps)*exp(-rc/2);
   /*printf(" g2f %e cc %e cd %e \n",f,cc,cd);
  */return f;
}



double g6f(double c1)
{  double   low=exp(-rc), eps=1.e-4;  
  double f;
  cc=c1;
  f=  simps(rgzptfac,low,1,eps)*exp(-rc/2); 
   /*printf(" g6f %e \n",f);
  */return   f;
}


double g8f(double c1)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f;
  cc=c1;
  xx=MLZP(c1,mm)/1000.;
  f=simps(rgzpfac,low,high,eps)*exp(-rc/2);
 /*printf("g8f=%e  cc=%e xx=%e\n",f,cc,xx);
 */ return f;
}

static double fh(double z)
{
  return sqrt( 2*z*z*exp(rc)/(1-exp(-2*rc)));
}

static double theta(double x1, double x2, double x3)
{ return 1./( 1+pow(-x1*x1+x2*x2-x3*x3+sqrt(-4*x1*x1*x2*x2+
pow(-x1*x1-x2*x2-x3*x3,2)),2 )/(4.*x2*x2*x3*x3) );
}

static double ml(double c1)
{return M_PI*(c1+1.)/2.*mzp/2.405;
}

static double ryt(double z)
{  double f;
  f=ftr(z,cc)*ftr(z,cd)*fh(z)*pow(z*exp(rc),-4);;
 return f;
}

static double yt(double c1, double c2)
{double  low=exp(-rc),high=1, eps=1.e-4;
  double f;
  cc=c1;
  cd=c2;
  f=simps(ryt,low,high,eps)*sqrt(exp(-rc));
/*printf("yt= %e c1=%e c2=%e fh=%e\n",f,c1,c2,fh(0.999));
*/
  return f; 
}

static double rykk(double z)
{  double f;
  xx=MLZP(cc,mm)/1000.;
  xy=MLZP(cd,mm)/1000.;
  f=fLKP(z,cc,xx)*fLKP(z,cd,xy)*fh(z)*pow(z*exp(rc),-4);
return f;
}

static double ykk(double c1, double c2)
{double low=exp(-rc),high=1, eps=1.e-4;
  double f;
  cc=c1;
  cd=c2;
  f=simps(rykk,low,high,eps)*sqrt(exp(-rc));
/*printf("ykk= %e c1=%e c2=%e\n",f,c1,c2);
 */ return f; 
}

static double mlr(double c1,double c2, double c3, double c4)
{ 
  double f=mt*ykk(c1,c2)/yt(c3,c4);
  return f;
}


static double gzp2(double c1)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f;
  cc=c1;
  xx=MLZP(c1,mm)/1000.;
   
  f=simps(rgzpfac,low,high,eps)*exp(-rc/2);
  return f*sqrt(rc*5./2.)*g10/2.;
}

static double rgzph(double z)
{ double f=fgauge(z,mm)*fh(z)*fh(z)*pow(z*exp(rc),-1);
  return f;
}

double gztot(double c1, double c2, double c3, double c4)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f,f1,gznu;
  cc=c1;
  cd=c2;  
  gznu=ee/2./sqrt(sw2*(1-sw2)); 
  f1=gznu*simps(rgzph,low,high,eps)*sqrt(rc*exp(-rc))*g10/sqrt(10)*v*v/mzp/mzp
     *gzp2(c2);
  f=f1+ gznu*theta(MLZP(c2,mzp)/1000.,ml(c1),mlr(c1,c2,c3,c4)) ;
 
/* printf(" gztot %e mlzp%e  ml%e mlr%e theta=%e \n",f,MLZP(c2,mzp)/1000.,ml(c1),mlr(c1,c2,c3,c4),theta(MLZP(c2,mzp)/1000.,ml(c1),mlr(c1,c2,c3,c4)));
  */return f1;
}


double mixz(double m2)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f,gznu,mtev=m2/1000.;
  gznu=ee/2./sqrt(sw2*(1-sw2)); 
  f=gznu*simps(rgzph,low,high,eps)*sqrt(rc*exp(-rc))*g10/sqrt(10)*v*v/mtev/mtev;
  return f;
}


double ghlzp(double c1, double c2, double c3, double c4, double m1, double m2)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f,f1,gznu,gold;
  double Mtp=175;
  double sinL,cosR,delta,m3,X,test1,test2,test3,test4,l1,cosL,sinR,l2;  
  cc=c1;
  cd=c2;  
 /* Here c3 is cnur c4 cnul', c1=ctl c2=ctr*/
 
  if(c3>-0.5)
 { m3=Mtp*sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2));
  delta=pow((m1*m1+m2*m2+m3*m3),2)-4*m1*m1*m2*m2;
  X=-m1*m1+m2*m2-m3*m3+sqrt(delta);
  sinL=1./sqrt(1+X*X/(4.*m1*m1*m3*m3));
  cosR= -X/2/m2/m3/sqrt(1+X*X/(4.*m2*m2*m3*m3));
  f=sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2))*sinL*cosR;
  }
    else
 { m3=Mtp*sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2))/sqrt(2/(1-2*c3));
  delta=pow((m1*m1+m2*m2+m3*m3),2)-4*m1*m1*m2*m2;
  X=-m1*m1+m2*m2-m3*m3+sqrt(delta);
  sinL=1./sqrt(1+X*X/(4.*m1*m1*m3*m3));
  cosR= -X/2/m2/m3/sqrt(1+X*X/(4.*m2*m2*m3*m3));
  f=sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2))/sqrt(2/(1-2*c3))*sinL*cosR;
 }
  gold=2.*Mtp*(1-2*c3)/(1-2*c1)/(1-2*c2)*m1/(sqrt(2)*m2*m2);
/* printf(" ghlzp %e sinR^2 %e  mlr %e  old %e mass %e\n",f,(1-cosR*cosR),Mtp*sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2)),gold,m1);
printf(" m1= %e m2=%e  m3= %e \n",m1,m2,m3);
*/
/*   m3=Mtp*sqrt(2/(1-2*c1))*sqrt(2/(1-2*c2));
    delta=pow((m1*m1+m2*m2+m3*m3),2)-4*m1*m1*m2*m2; 
   X=-m1*m1+m2*m2-m3*m3+sqrt(delta);
   sinL=1./sqrt(1+X*X/(4.*m1*m1*m3*m3));
   cosR= X/2/m2/m3/sqrt(1+X*X/(4.*m2*m2*m3*m3));  
   cosL= -X/2/m1/m3/sqrt(1+X*X/(4.*m1*m1*m3*m3));
   sinR=1./sqrt(1+X*X/(4.*m2*m2*m3*m3));
   l1=sqrt(1/2.*((m1*m1+m2*m2+m3*m3)-sqrt(delta)));
    l2=sqrt(1/2.*((m1*m1+m2*m2+m3*m3)+sqrt(delta)));
   test1=m1*cosR-l1*cosL;
   test2=-m1*sinR+l2*sinL;
   test3=m3*cosR+m2*sinR-l1*sinL;
   test4=-m3*sinR+m2*cosR-l2*cosL;
    printf(" l1= %e l2=%e  delta=%e\n",l1,l2,sqrt(delta));
   printf(" test1= %e test2=%e  test3= %e test4 %e\n",test1,test2,test3,test4);
printf(" mat11= %e mat12=%e  mat13=%e mat14=%e\n",m1*cosR,-m1*sinR,m3*cosR+m2*sinR,-m3*sinR+m2*cosR);
printf(" mat11= %e mat12=%e  mat13=%e mat14=%e\n",l1*cosL,-l2*sinL,l1*sinL,l2*cosL);
*/
  return f;
}

double ghbr(double c1, double c2, double c3, double c4, double m1, double m2)
{ double  low=exp(-rc),high=1, eps=1.e-4;
  double f,f1,gznu,gold;
  double Mtp=175;
  double cosR,delta,m3,m4,X,X1,aL,aLp,ftl,fbr,ftr;  
  cc=c1;
  cd=c2;  
 /* Here c3 is cbr c4 cnul', c1=ctl c2=ctr*/
  ftl=sqrt(2/(1-2*c1));
  ftr=sqrt(2/(1-2*c2));
  fbr=sqrt(2/(1-2*c3));
  m3=Mtp*ftl*ftr;
  m4=Mtp*ftr;
  delta=pow((m1*m1+m2*m2+m3*m3+m4*m4),2)-4*(m1*m1+m4*m4)*m2*m2;
   X=m1*m1+m2*m2+m3*m3+m4*m4-sqrt(delta);
    X1=m1*m1-m2*m2+m3*m3+m4*m4-sqrt(delta);
   cosR= X1/2/m2/m3/sqrt(1+X1*X1/(4.*m2*m2*m3*m3));
   aL= 1./sqrt(1+m1*m1/m4/m4+m3*m3/m4/m4*X*X/X1/X1);
   aLp= aL*m3/m4*X/X1;
 if(c3>-0.5)
   {f=ftl*ftr*(aLp+aL/ftl)*cosR;
 }
    else
 { 
  f=ftl*ftr*(aLp+aL/ftl)*cosR/fbr;
 }
  
/*printf(" ghbr %e %e \n",f,cosR);*/
  return f;
}

