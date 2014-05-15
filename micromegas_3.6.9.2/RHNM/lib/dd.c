#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include"pmodel.h"


static double ee=0.31223;
static double sw=0.4730;
static double mz=91.1884;
static double mp=0.938; 
static double mn=0.939;
static double ftn=0.3776;
static double ftp=0.3755;


extern double sigmap(double gzp)
{ double f;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
 
f=mp*mp*ee*ee*(1-4.*sw2)*(1-4.*sw2)/64/M_PI/sw2/cw2/mz/mz/mz/mz*gzp*gzp*.389e9;
printf("sigmap= %.5e coeff %.5e \n",mz,f/gzp/gzp);
return f;
} 

extern double sigman(double gzp)
{ double f;
  double sw2=sw*sw;
  double cw2=(1-sw2);
  f=mn*mn*ee*ee/64/M_PI/sw2/cw2/mz/mz/mz/mz*gzp*gzp*.389e9;
printf("sigman= %.5e coeff %.5e \n",mz,f/gzp/gzp);
 return f;
} 

extern double sigmanh(double gzp, double gh, double mh, double mnu)
{ double f,fn,fp,s;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
   double cw=sqrt(cw2);
double mw=mz*cw;
fn=gzp-mn/mw*mz*mz/mh/mh*4.*cw*gh*ftn;
fp=gzp*(1.-4.*sw2)-mp/mw*mz*mz/mh/mh*4.*cw*gh*ftp;
f=mn*mn*mnu*mnu/(mn+mnu)/(mn+mnu)*ee*ee*(fn*fn)/64./M_PI/sw2/cw2/mz/mz/mz/mz*.389e9;
printf("fn= %.5e gzp=%.5e h=%.5e\n",fn,gzp,mn/mw*mz*mz/mh/mh*4.*cw*gh*ftn);
return f;
} 

extern double sigmaph(double gzp, double gh, double mh, double mnu)
{ double f,fn,fp;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
   double cw=sqrt(cw2);
double mw=mz*cw;
fn=gzp-mn/mw*mz*mz/mh/mh*4.*cw*gh*ftn;
fp=gzp*(1.-4.*sw2)-mp/mw*mz*mz/mh/mh*4.*cw*gh*ftp;
f=mn*mn*mnu*mnu/(mn+mnu)/(mn+mnu)*ee*ee*(fp*fp)/64./M_PI/sw2/cw2/mz/mz/mz/mz*.389e9;
return f;
} 


extern double sigmahe(double gzp)
{ double f;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
/* for Helium3*/
f=9.*mp*mp*ee*ee*(2.*(1-4.*sw2)-2)*(2*(1-4.*sw2)-2.)/64./M_PI/sw2/cw2/mz/mz/mz/mz*gzp*gzp*.389e9;
return f;
} 

extern double sigmage(double gzp, double gh, double mh, double mnu)
{ double f,fn,fp,s;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
   double cw=sqrt(cw2);
   double mw=mz*cw;
/* for Xenon131*/
fn=gzp-mn/mw*mz*mz/mh/mh*4.*cw*gh*ftn;
fp=gzp*(1.-4.*sw2)-mp/mw*mz*mz/mh/mh*4.*cw*gh*ftp;
s=mn*mn*mnu*mnu/(mn+mnu)/(mn+mnu)*ee*ee*(fn*fn)/64./M_PI/sw2/cw2/mz/mz/mz/mz*.389e9;
f=s*(41+fp/fn*32)*(41+fp/fn*32)/73/73;
return f;
} 

extern double sigmaxe(double gzp, double gh, double mh, double mnu)
{ double f,fn,fp,s;
  double sw2=sw*sw;
  double cw2=(1-sw2); 
   double cw=sqrt(cw2);
   double mw=mz*cw;
/* for Xenon131*/
fn=gzp-mn/mw*mz*mz/mh/mh*4.*cw*gh*ftn;
fp=gzp*(1.-4.*sw2)-mp/mw*mz*mz/mh/mh*4.*cw*gh*ftp;
/*f=131*mn*131*mn*mnu*mnu/(131*mn+mnu)/(131*mn+mnu)*ee*ee*(54*fp+77*fn)*(54*fp+77*fn)/64./M_PI/sw2/cw2/mz/mz/mz/mz*.389e9;
*/
s=mn*mn*mnu*mnu/(mn+mnu)/(mn+mnu)*ee*ee*(fn*fn)/64./M_PI/sw2/cw2/mz/mz/mz/mz*.389e9;
f=s*(77+fp/fn*54)*(77+fp/fn*54)/131/131;
return f;
} 


