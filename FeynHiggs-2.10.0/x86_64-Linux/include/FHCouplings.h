#if 0
	FHCouplings.h
		human-readable indexing for the
		couplings, gammas, and gammasms arrays
		this file is part of FeynHiggs
		last modified 15 Nov 13 th

Note1: comments are real funny here because we want to include
       this file in both Fortran and C

Note2: for the same reason, the funny notation couplingS etc.
       is used because we have to remap array indices in C and
       Fortran does not care about caps.
#endif


#ifndef FHCOUPLINGS_H
#define FHCOUPLINGS_H

#define ncouplings 681
#define Roffset 472
#define Coupling(c)		couplingS(c)
#define LCoupling(c)		couplingS(c)
#define RCoupling(c)		couplingS(c+Roffset)

#define ncouplingsms 231
#define RSMoffset 108
#define CouplingSM(c)		couplingsmS(c)
#define LCouplingSM(c)		couplingsmS(c)
#define RCouplingSM(c)		couplingsmS(c+RSMoffset)

#define ngammas 978
#define BRoffset 491
#define GammaTot(h)		gammaS(h)
#define Gamma(c)		gammaS(c+4)
#define BR(c)			gammaS(c+BRoffset)

#define ngammasms 250
#define BRSMoffset 127
#define GammaSMTot(h)		gammasmS(h)
#define GammaSM(c)		gammasmS(c+4)
#define BRSM(c)			gammasmS(c+BRSMoffset)


#define H0VV(h,vv) h+3*vv-3
#if 0
  h  = 1..3	Higgs: h0, HH, A0
  vv = 1..5	vector-boson pair: gamma gamma, gamma Z, ZZ, WW, gg
#endif

#define H0FF(h,t,g1,g2) h+3*t+12*g1+36*g2-36
#if 0
  h  = 1..3	Higgs: h0, HH, A0
  t  = 1..4	fermion type: nu, e, u, d
  g1 = 1..3	fermion 1 generation
  g2 = 1..3	fermion 2 generation
#endif

#define HpFF(p,g1,g2) p+2*g1+6*g2+115
#if 0
  p  = 1..2	decay products: leptons, quarks
  g1 = 1..3	up-type fermion 1 generation
  g2 = 1..3	down-type fermion 2 generation
#endif

#define H0ChaCha(h,c1,c2) h+3*c1+6*c2+132
#if 0
  h  = 1..3	Higgs: h0, HH, A0
  c1 = 1..2	chargino 1
  c2 = 1..2	chargino 2
#endif

#define H0NeuNeu(h,n1,n2) h+3*n1+12*n2+138
#if 0
  h  = 1..3	Higgs: h0, HH, A0
  n1 = 1..4	neutralino 1
  n2 = 1..4	neutralino 2
#endif

#define HpNeuCha(n1,c2) n1+4*c2+197
#if 0
  n1 = 1..4	neutralino
  c2 = 1..2	chargino
#endif

#define H0HV(h,hv) h+3*hv+206
#if 0
  h  = 1..3	decaying Higgs: h0, HH, A0
  hv = 1..3	produced pair: h0-Z, HH-Z, A0-Z
#endif

#define HpHV(hv) hv+218
#if 0
  hv = 1..3	produced pair: h0-W, HH-W, A0-W
#endif

#define H0HH(h,h1,h2) h+3*h1+12*h2+206
#if 0
  h  = 1..3	decaying Higgs: h0, HH, A0
  h1 = 1..4	produced Higgs 1: h0, HH, A0, Hp
  h2 = 1..4	produced Higgs 2: h0, HH, A0, Hp
#endif

#define H0SfSf(h,s1,s2,t,g) h+3*s1+6*s2+12*t+48*g+200
#if 0
  h  = 1..3	Higgs: h0, HH, A0
  s1 = 1..2	sfermion 1
  s2 = 1..2	sfermion 2
  t  = 1..4	sfermion type: nu, e, u, d
  g = 1..3	common sfermion generation
#endif

#define HpSfSf(s1,s2,p,g1,g2) s1+2*s2+4*p+8*g1+24*g2+375
#if 0
  s1 = 1..2	sfermion 1
  s2 = 1..2	sfermion 2
  p  = 1..2	decay products: sleptons, squarks
  g1 = 1..3	up-type sfermion 1 generation
  g2 = 1..3	down-type sfermion 2 generation
#endif

#define tBF(bf) bf+485
#if 0
  bf = 1..2	W-b, H-b
#endif


#define nprodxs 52

#define bbh(h)		prodXS(h)
#define bbhSM(h)	prodXS(h+3)
#define btagbh(h)	prodXS(h+6)
#define btagbhSM(h)	prodXS(h+9)
#define ggh(h)		prodXS(h+12)
#define ggh2(h)		prodXS(h+15)
#define gghSM(h)	prodXS(h+18)
#define qqh(h)		prodXS(h+21)
#define qqhSM(h)	prodXS(h+24)
#define tth(h)		prodXS(h+27)
#define tthSM(h)	prodXS(h+30)
#define Wh(h)		prodXS(h+33)
#define WhSM(h)		prodXS(h+36)
#define Zh(h)		prodXS(h+39)
#define ZhSM(h)		prodXS(h+42)
#define StSth(h)	prodXS(h+45)
#define tHm		prodXS(49)
#define tHm2		prodXS(50)
#define tHm2lo		prodXS(51)
#define tHm2hi		prodXS(52)

#endif

