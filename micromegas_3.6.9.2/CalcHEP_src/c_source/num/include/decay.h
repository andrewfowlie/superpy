#ifndef __DECAY__
#define __DECAY__

#include"nType.h"

typedef   struct 
{
     REAL a0, a1, a2;
     int n0, n1,n2,n3;
     REAL p0, p1;
     REAL v0, aa, pa, va, oo, po, pp, vo, vp, vv, e2a, am1, am2, e1p, e2p, e0v, e1v, e2v;
     REAL  cath, cinp, ceps, cpol, pout, e1out, e2out;
}  iDecay;


extern  void decay_0(int nvin, REAL amm1, REAL amm2, double *factor, double *Emax,iDecay*memDecay,REAL*V);
extern  void decay_1(int nvpole, REAL *hsum, REAL *hdif,iDecay*memDecay,REAL*V);
extern  void decay_2(int nvout1, double *parcos,iDecay*memDecay,REAL*V);
extern  void decay_3(int nvath, REAL parcos, REAL parfi, int nvout1,int nvout2,iDecay*memDecay,REAL*V);
extern  void wrmom(int n);

#endif
