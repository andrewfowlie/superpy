/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include<math.h>
#include"4_vector.h"
#include"decay.h"

#include"nType.h"

/*
  bases: e0..e3 : v=V0*e0  ;  p=P0*e0 + P1*e1; a =A0*e0 + A1*e1 + A2*e2
 
                 e0=E0V*v ; e1=E1V*v + E1P*p; e2=E2V*v + E2P*p + E2A*a
 
   out1= E1OUT*e0 + POUT*( PARCOS*e1 + PARSIN*(e2*COSFI+ e3*SINFI)) 
       = E1OUT*E0V*v + POUT*( PARCOS*(E1V*v + E1P*p ) 
       + PARSIN*(COSFI*(E2V*v + E2P*p + E2A*a )  + SINFI*e3 )) 
       =(E1OUT*E0V + POUT*( PARCOS*E1V +PARSIN*COSFI*E2V))*v 
       +POUT*(PARCOS*E1P+PARSIN*COSFI*E2P)*p 
       +POUT*PARSIN*COSFI*E2A*a 
       +POUT*PARSIN*SINFI*e3 
*/



#define  a0    memDecay->a0
#define  a1    memDecay->a1
#define  a2    memDecay->a2
#define  n0    memDecay->n0
#define  n1    memDecay->n1
#define  n2    memDecay->n2
#define  n3    memDecay->n3
#define  p0    memDecay->p0
#define  p1    memDecay->p1
#define  v0    memDecay->v0
#define  aa    memDecay->aa
#define  pa    memDecay->pa
#define  va    memDecay->va
#define  oo    memDecay->oo
#define  po    memDecay->po
#define  pp    memDecay->pp
#define  vo    memDecay->vo
#define  vp    memDecay->vp
#define  vv    memDecay->vv
#define  e2a   memDecay->e2a
#define  am1   memDecay->am1
#define  am2   memDecay->am2
#define  e1p   memDecay->e1p
#define  e2p   memDecay->e2p
#define  e0v   memDecay->e0v
#define  e1v   memDecay->e1v
#define  e2v   memDecay->e2v
#define  cath  memDecay->cath
#define  cinp  memDecay->cinp
#define  ceps  memDecay->ceps
#define  cpol  memDecay->cpol
#define  pout  memDecay->pout
#define  e1out memDecay->e1out
#define  e2out memDecay->e2out
                                   
void decay_0(int nvin,  REAL  amm1, REAL  amm2, double *factor,double *Emax, iDecay*memDecay, REAL*V)
{
    REAL  summ, diffm;
    const  double c_pi__ = 16*M_PI*M_PI;
    REAL chY,shY;
    am1 = amm1;
    am2 = amm2;
    n0 = nvin;
    vv = vdot4(nvin, nvin,V);
    
if(vv<0)printf("decay0: vv=%e m1=%e m2=%e\n", vv,amm1, amm2);
    if(vv<0) { vv=0; fprintf(stderr,"decay_0: negative incoming mass\n");}
    v0 = sqrt(vv);
    e0v = 1 / v0;
/* Computing 2nd power */
    summ = am1 + am2;
    diffm = am1 - am2;
    if(v0<summ){pout=0;  fprintf(stderr,"decay_0: virtual decay\n");}
    else pout = sqrt((v0-summ)*(v0+ summ)*(v0-diffm)*(v0+diffm))/(v0*2);
    e1out = sqrt(am1*am1 + pout*pout);
    e2out = sqrt(am2*am2 + pout*pout);
    *factor = pout/(c_pi__ * (e1out + e2out));
    chY=V[4*nvin-4]/v0;
    shY=sqrt(chY*chY-1);
    Emax[0]= e1out*chY+pout*shY;
    Emax[1]= e2out*chY+pout*shY;    
}

void decay_1(int nvpole, REAL *hsum, REAL *hdif,iDecay*memDecay,REAL*V)
{
    n1 = nvpole;
    vp = vdot4(n0, nvpole,V);
    pp = vdot4(nvpole, nvpole,V);
    p0 = e0v * vp;
    p1 = sqrt(p0*p0 - pp);
    *hdif = p1 * 2 * pout;
    hsum[0] = pp + am1*am1 + p0 * 2 * e1out;
    hsum[1] = pp + am2*am2 + p0 * 2 * e2out;
}

void decay_2(int nvout1,  double * parcos,iDecay*memDecay,REAL*V)
{
    oo = vdot4(nvout1, nvout1,V);
    vo = vdot4(n0, nvout1,V);
    po = vdot4(n1, nvout1,V);
    *parcos = (po*vv - vp*vo)/sqrt((pp*vv - vp*vp) * (oo*vv - vo*vo));
}

void decay_3(int nvath, REAL parcos, REAL parfi, int nvout1, int nvout2,iDecay*memDecay,REAL*V)
{
    int i;
    REAL lparcos, lparsin,cosfi, sinfi;
/* *   v = NVIN; p = NVPOLE; a = NVATH; */
    va = vdot4(n0, nvath,V);
    pa = vdot4(n1, nvath,V);
    aa = vdot4(nvath, nvath,V);
    e1p = 1 / p1;
    e1v = -e1p * p0 * e0v;
    a0 = e0v * va;
    a1 = -(e1v * va + e1p * pa);
    a2 = sqrt(a0*a0 - a1*a1 - aa);
    e2a = 1 / a2;
    e2p = -e2a * a1 * e1p;
    e2v = -e2a * (a0 * e0v + a1 * e1v);
/* **** out state ***** */
    lparcos=parcos;
    if(lparcos>=1) {lparcos=1; lparsin=0;} 
    else if(lparcos<=-1){lparcos=-1;lparsin=0;}
    else   lparsin = sqrt((1 - lparcos) * (lparcos + 1));
    
    sinfi = sin((double)parfi);
    cosfi = cos((double)parfi);
    n3 = nvout1;
    n2 = nvath;
    eps4(n0, n1, n2, n3,V);
 
    cinp = pout*(lparcos*e1v + lparsin*cosfi*e2v);
    cpol = pout*(lparcos*e1p + lparsin*cosfi*e2p);
    cath = pout*lparsin*cosfi*e2a;
    ceps = pout*lparsin*sinfi/sqrt(-vdot4(n3, n3,V));

	
    for (i = 0; i <= 3; ++i) 
    {
      V[i+4*n3-4] = cinp*V[i+4*n0-4] + cpol*V[i+4*n1-4] 
                       + cath*V[i+4*n2-4] + ceps*V[i+4*n3-4];
    }


    for (i = 0; i <= 3; ++i) 
    {
       V[i+4*nvout2-4] = e2out*e0v*V[i+4*n0-4] - V[i+4*n3-4];
       V[i+4*nvout1-4] = e1out*e0v*V[i+4*n0-4] + V[i+4*n3-4];
    }   

    
/*
   if(am1!=0 && vdot4(nvout1, nvout1)<0) printf("nvout1^2=%E am1=%E\n",vdot4(nvout1, nvout1),am1);
   if(am2!=0 && vdot4(nvout2, nvout2)<0) printf("nvout2^2=%E am2=%E\n",vdot4(nvout2, nvout2),am2);
*/
}


