

/*        SuSpect fortran code
c--------------------------------------------------
c   calculates leading one-loop SUSY delta_rho contributions of 3rd gen
c sfermions (plus leading two-loop QCD contributions) 
c  INPUT: MT, gmst(2), gmsb(2),gmstau(2),msn: top,stop,sbottom,
c  stau, stau neutrino masses and stop, sbottom, stau mixing angles
c  OUTPUT: drho = rho-1 
c--------------------------------------------------
*/

#include"../../sources/micromegas.h"
#include<math.h>
#define   su_fr(x,y)  ((x)+(y)-2*(x)*(y)/((x)-(y))*log((x)/(y)))

static double sq(double x) {return x*x;}

double  deltarho_(void)
{

      double mt,msn,thetat,thetab,thel;
      double GF=1.16639E-5;
      double ct,st,cb,ctau ,stau,cta2,sta2,ct2,st2,cb2,sb2,drho;
      double drhotau,drhotb,mt1,mt2,mb1,mb2,mta1,mta2,sb;

 
      mt=findValW("Mtp");
      thetat=atan2(findValW("Zt12"),findValW("Zt11"));
      thetab=atan2(findValW("Zb12"),findValW("Zb11"));
      thel=atan2(findValW("Zl12"),findValW("Zl11"));
/*printf("deltarho: %E %E %E %E\n",mt,thetat,thetab,thel);*/
      ct=cos(thetat);
      st=sin(thetat);
      cb=cos(thetab);
      sb=sin(thetab);
      ctau =cos(thel);
      stau =sin(thel);
      cta2=sq(ctau);
      sta2=sq(stau);
      ct2=sq(ct);
      st2=sq(st);
      cb2=sq(cb);
      sb2=sq(sb);

      mt1=sq(findValW("MSt1"));
      mt2=sq(findValW("MSt2"));
      mb1=sq(findValW("MSb1"));
      mb2=sq(findValW("MSb2"));
      mta1=sq(findValW("MSl1"));
      mta2=sq(findValW("MSl2"));
      msn=findValW("MSnl");

      drhotb= (ct2*(cb2*su_fr(mt1,mb1)+sb2*su_fr(mt1,mb2)) +
           st2*(cb2*su_fr(mt2,mb1)+sb2*su_fr(mt2,mb2)) -
           ct2*st2*su_fr(mt1,mt2)-cb2*sb2*su_fr(mb1,mb2));
      drhotau= -cta2*sta2*su_fr(mta1,mta2)+cta2*su_fr(mta1,sq(msn)) +
     sta2*su_fr(mta2,sq(msn));
      drho = 3*drhotb*(1. +2*0.12/3/M_PI*(1.+M_PI*M_PI/3))+drhotau;
      return  GF/(8* M_PI*M_PI* sqrt(2.))*drho;
}
