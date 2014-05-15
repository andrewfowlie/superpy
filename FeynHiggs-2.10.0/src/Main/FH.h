* FH.h
* global variable declarations
* this file is part of FeynHiggs
* last modified 5 Nov 13 th


#ifndef SignSq
#define SignSq(x) (x)*abs(x)
#define SignSqrt(x) sign(sqrt(abs(Re(x))),Re(x))
#define signbit(i) ibits(i,31,1)
#define Delta(i,j) signbit(ieor(i,j)-1)
#define Key(se) 2**(se-1)

#define LOOP(var,from,to,step) do var = from, to, step
#define ENDLOOP(var) enddo

#define isQ(t) ibits(t+1,2,1)
#define ifQ(t,i) iand(-isQ(t),i)
#define isB(t) ibits(t,2,1)
#define ifB(t,i) iand(-isB(t),i)
#define isMFV(t) ibits(not(fv),t,1)
#define ifMFV(t,i) iand(-isMFV(t),i)
#define Ncolor(t) ior(t-1,1)

#define tQ(t) Ncolor(t)
#define tU(t) (t + isQ(t))

* for encoding sfermion type in SfUpdate and Couplings.F:
#define X2(x1,x0) (x1)*16 + x0
#define X3(x2,x1,x0) (x2)*256 + X2(x1,x0)
#define X4(x3,x2,x1,x0) (x3)*4096 + X3(x2,x1,x0)
#define X5(x4,x3,x2,x1,x0) (x4)*65536 + X4(x3,x2,x1,x0)
#define nib4(x) ibits(x,16,4)
#define nib3(x) ibits(x,12,4)
#define nib2(x) ibits(x,8,4)
#define nib1(x) ibits(x,4,4)
#define nib0(x) ibits(x,0,4)

#define Mf(t,g) Sf(g,t)
#define Mf2(t,g) Sf(g+3,t)
#define MSf(s,t,g) Sf(s+2*(g)+4,t)
#define MASf(as,t) Sf(as+12,t)
#define MSf2(s,t,g) Sf(s+4*(g)+14,t)
#define MASf2(as,t) Sf(as+30,t)
#define USf2(s1,s2,t,g) Sf(s1+2*(s2)+4*(g)+36,t)
#define USf(s1,s2,t,g) CSf(s1+2*(s2)+4*(g)+21,t)
#define USf_flat(i,t) CSf(i+27,t)
#define UASf(as1,as2,t) CSf(as1+6*(as2)+33,t)
#define UASf_flat(i,t) CSf(i+39,t)
#define DSS2(s,t,g) Sf(s+2*(g)+148,t)
#define Kf(g1,g2,t) CSf(g1+3*(g2)+75,t)
#define Kfflat(i,t) CSf(i+78,t)
#define Deltaf(t,g) CSf(g+87,t)
#define CKM(g1,g2) CSf(g1+3*(g2)+84,1)
#define CKM_flat(i) CSf(i+87,1)
#define NSf 194

#define USfC(s1,s2,t,g) Conjugate(USf(s1,s2,t,g))
#define UASfC(as1,as2,t) Conjugate(UASf(as1,as2,t))
#define VChaC(c1,c2) Conjugate(VCha(c1,c2))
#define UChaC(c1,c2) Conjugate(UCha(c1,c2))
#define ZNeuC(n1,n2) Conjugate(ZNeu(n1,n2))
#define USdLC(s1,s2,g) Conjugate(USdL(s1,s2,g))
#define VChaLC(c1,c2) Conjugate(VChaL(c1,c2))
#define UChaLC(c1,c2) Conjugate(UChaL(c1,c2))
#define ZNeuLC(n1,n2) Conjugate(ZNeuL(n1,n2))
#define CKMC(g1,g2) Conjugate(CKM(g1,g2))
#define CKMinC(g1,g2) Conjugate(CKMin(g1,g2))
#define KfC(g1,g2,t) Conjugate(Kf(g1,g2,t))
#define AfC(t,g) Conjugate(Af(t,g))
#define XfC(t,g) Conjugate(Xf(t,g))
#define MUEC Conjugate(MUE)
#define M_3C Conjugate(M_3)

#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(i+3,j,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)

#define MSS2_LL(i,j) (MSS(nQ,i)*MSS(nQ,j))
#define MSS2_LR(i,j) (MSS(nQ,i)*MSS(nU,j))
#define MSS2_RL(i,j) (MSS(nU,i)*MSS(nQ,j))
#define MSS2_RR(i,j) (MSS(nU,i)*MSS(nU,j))

#define MSS02_LL(i,j) (MSS0(nQ,i)*MSS0(nQ,j))
#define MSS02_LR(i,j) (MSS0(nQ,i)*MSS0(nU,j))
#define MSS02_RL(i,j) (MSS0(nU,i)*MSS0(nQ,j))
#define MSS02_RR(i,j) (MSS0(nU,i)*MSS0(nU,j))
#endif


#include "const.h"
#include "FHCouplings.h"


* SM parameters

	ComplexType CKMin(3,3)
	RealType CKMlambda, CKMA, CKMrhobar, CKMetabar
	RealType Qf(4), MB_MT
	RealType MW, MW2, MZ, MZ2
	RealType SW, SW2, CW, CW2
	RealType invAlfaMZ, GF, vev
	RealType ELGF, AlfaGF, EL0, ELMZ, AlfaMZ, Alfat
	RealType GSMT, AlfasMT, AlfasMZ, AlfasDb

	common /smpara/
     &    CKMin, CKMlambda, CKMA, CKMrhobar, CKMetabar,
     &    Qf, MB_MT,
     &    MW, MW2, MZ, MZ2, CW, CW2, SW, SW2,
     &    invAlfaMZ, GF, vev,
     &    ELGF, AlfaGF, EL0, ELMZ, AlfaMZ, Alfat,
     &    GSMT, AlfasMT, AlfasMZ, AlfasDb

	RealType Alfa1L, Alfa2L, EL1L, EL2L
	equivalence (AlfaGF, Alfa1L, Alfa2L)
	equivalence (ELGF, EL1L, EL2L)

	RealType Alfas2L, GS2L
	equivalence (AlfasMT, Alfas2L)
	equivalence (GSMT, GS2L)

	RealType ME, ME2, MM, MM2, ML, ML2
	RealType MU, MU2, MC, MC2, MT, MT2
	RealType MD, MD2, MS, MS2, MB, MB2
	equivalence (Mf(2,1), ME), (Mf2(2,1), ME2)
	equivalence (Mf(2,2), MM), (Mf2(2,2), MM2)
	equivalence (Mf(2,3), ML), (Mf2(2,3), ML2)
	equivalence (Mf(3,1), MU), (Mf2(3,1), MU2)
	equivalence (Mf(3,2), MC), (Mf2(3,2), MC2)
	equivalence (Mf(3,3), MT), (Mf2(3,3), MT2)
	equivalence (Mf(4,1), MD), (Mf2(4,1), MD2)
	equivalence (Mf(4,2), MS), (Mf2(4,2), MS2)
	equivalence (Mf(4,3), MB), (Mf2(4,3), MB2)

	ComplexType CKMin_flat(3*3)
	equivalence (CKMin, CKMin_flat)


* MSSM parameters

	ComplexType UCha(2,2), VCha(2,2), ZNeu(4,4)
	ComplexType MSS2(3,3,5), deltaSf(6,6,2:4)
	ComplexType Xf(4,3), Af(4,3), Af0(2:4,3)
	ComplexType MUETB(2:4), MUE, M_1, M_2, M_3
	RealType MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	RealType MSS(5,3), MSS0(5,3)
	RealType DSf(2,5), QSf(2:4)
        RealType MHtree(4), MHtree2(4)
	RealType MGl, MGl2
	RealType CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	RealType CA, SA, CA2, SA2, C2A, S2A
	RealType CAB, SAB, CBA, SBA, CBA2, SBA2, SCB(2:4)
	RealType scalefactor
	integer inputmass

	common /mssmpara/
     &    UCha, VCha, ZNeu,
     &    MSS2, deltaSf,
     &    Xf, Af, Af0,
     &    MUETB, MUE, M_1, M_2, M_3,
     &    MCha, MCha2, MNeu, MNeu2,
     &    MSS, MSS0,
     &    DSf, QSf,
     &    MHtree, MHtree2,
     &    MGl, MGl2,
     &    CB, SB, TB, CB2, SB2, TB2, C2B, S2B,
     &    CA, SA, CA2, SA2, C2A, S2A,
     &    CAB, SAB, CBA, SBA, CBA2, SBA2, SCB,
     &    scalefactor,
     &    inputmass

	RealType Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	equivalence (MHtree(1), Mh0), (MHtree2(1), Mh02)
	equivalence (MHtree(2), MHH), (MHtree2(2), MHH2)
	equivalence (MHtree(3), MA0), (MHtree2(3), MA02)
	equivalence (MHtree(4), MHp), (MHtree2(4), MHp2)

	ComplexType MSS2_flat(3*3*5)
	equivalence (MSS2, MSS2_flat)

	RealType reimMUE(2), reMUE, imMUE
	equivalence (MUE, reimMUE)
	equivalence (reimMUE(1), reMUE)
	equivalence (reimMUE(2), imMUE)

	ComplexType UCha_flat(2*2), VCha_flat(2*2), ZNeu_flat(4*4)
	equivalence (UCha, UCha_flat)
	equivalence (VCha, VCha_flat)
	equivalence (ZNeu, ZNeu_flat)

	ComplexType deltaSf_flat(6*6*2)
	equivalence (deltaSf, deltaSf_flat)

* variants for large TB:

	ComplexType UChaL(2,2), VChaL(2,2)
	ComplexType ZNeuL(4,4), USdL(2,2,3)
	RealType MChaL(2), MNeuL(4)
	RealType MSdL(2,3), MSdL2(4,3)

	common /mssmparaLargeTB/
     &    UChaL, VChaL, ZNeuL, USdL,
     &    MChaL, MNeuL, MSdL, MSdL2


* Note: despite its name, sfermpara contains not only
* sfermion parameters, but all variables which need to be
* conserved during FHUncertainties.

* Sf(*,1) = Sneutrino				- set in Para.F
* Sf(*,2) = Slepton				- set in Para.F
* Sf(*,3) = Sup with MT(pole)			- set in Sfermions.F
* Sf(*,4) = Sdown with MB(MB)			- set in Sfermions.F
*
* Sf(*,5=bBR) = Sdown with MB(MB)/|1 + Db|	- set in Sfermions.F
*
* Sf(*,6=tT) = Sup with MT(MT)			- set in Sfermions.F
* Sf(*,7=bTR) = Sdown with MB(MT)/|1 + Db|	- set in Sfermions.F
* Sf(*,8=bTR0) = ditto but compatible with TLps	- set in TLShifts.F
*   (latter used for neutral Higgs masses only)
*
* Sf(*,9=tH) = Sup with MT(Mh) for Decays	- set in Couplings.F
* Sf(*,10=bH) = Sdown with MB(Mh) for Decays	- set in Couplings.F
* Sf(*,11=bHR) = Sdown with MB(Mh)/|1 + Db|	- set in Couplings.F

	integer SfSlots
	integer*8 SfIni
	parameter (SfSlots = 11, SfIni = O'12344344344')

	RealType Sf(NSf,SfSlots)

	common /sfermpara/ Sf

	RealType Sf_flat(NSf*SfSlots)
	ComplexType CSf(NSf/2,SfSlots)
	equivalence (Sf, Sf_flat, CSf)

	integer bBR, tT, bTR, bTR0, tH, bH, bHR
	parameter (bBR = 5, tT = 6, bTR = 7, bTR0 = 8)
	parameter (tH = 9, bH = 10, bHR = 11)


* Higgs results

	integer h0h0, HHHH, A0A0, HmHp
	integer h0HH, h0A0, HHA0
	integer G0G0, h0G0, HHG0, A0G0
	integer GmGp, HmGp
	integer semax
	integer cpeven, cpodd, goldstones
	parameter (h0h0 = 1, HHHH = 2, A0A0 = 3, HmHp = 4)
	parameter (h0HH = 5, h0A0 = 6, HHA0 = 7)
	parameter (G0G0 = 8, h0G0 = 9, HHG0 = 10, A0G0 = 11)
	parameter (GmGp = 12, HmGp = 13)
	parameter (semax = HmGp)
	parameter (cpeven = Key(h0h0) + Key(HHHH) + Key(h0HH))
	parameter (cpodd = Key(A0A0) + Key(h0A0) + Key(HHA0))
	parameter (goldstones = Key(G0G0) + Key(h0G0) +
     &    Key(HHG0) + Key(A0G0))

	integer se11, se12, se22, se1A, se2A, seAA
	parameter (se11 = 1, se12 = 2, se22 = 3)
	parameter (se1A = 4, se2A = 5, seAA = 6)

	integer NNeutral, NCharged, NHiggs
	parameter (NNeutral = 3)
	parameter (NCharged = 1)
	parameter (NHiggs = NNeutral + NCharged)
	ComplexType SAeff, XHiggs(0:NNeutral,0:NNeutral,0:2)
	RealType MHiggs(NHiggs), MHiggs2(0:NHiggs)

* renormalized self-energies & counter terms

	RealType dZ(semax), Msq(semax)
	ComplexType dMsq(semax)
	ComplexType seR(semax), dseR(semax), se2R(semax)
	ComplexType seAdd(semax)

	integer MaxVars, MaxSlots, Nvr
	parameter (MaxVars = 4)
	parameter (MaxSlots = 2**MaxVars)
	parameter (Nvr = 18)
	RealType monomial(MaxSlots)
	RealType vr(Nvr,MaxSlots)
	integer vdmb(MaxSlots)

	common /higgsdata/
     &    MHiggs, SAeff, XHiggs, MHiggs2,
     &    dZ, dMsq, Msq, seR, dseR, se2R, seAdd,
     &    monomial, vr, vdmb

	ComplexType vc(Nvr/2,MaxSlots)
	equivalence (vr, vc)

	ComplexType UHiggs(0:NNeutral,0:NNeutral)
	ComplexType ZHiggs(0:NNeutral,0:NNeutral)
	equivalence (XHiggs(0,0,1), UHiggs)
	equivalence (XHiggs(0,0,2), ZHiggs)

	integer NHiggsErr, NHiggsData
	parameter (NHiggsErr = NHiggs + 2 + 2*(NNeutral+1)**2*2)
	parameter (NHiggsData = NHiggsErr + NHiggs+1 +
     &    (3 + 3*2)*semax +
     &    MaxSlots + Nvr*MaxSlots)
	RealType HiggsData(NHiggsData)
	equivalence (MHiggs, HiggsData)


* couplings and widths

	ComplexType couplings(ncouplings)
	ComplexType couplingsms(ncouplingsms)
	RealType gammas(ngammas)
	RealType gammasms(ngammasms)
	RealType ratios(H0FF(3,4,3,3))
	RealType chSt1St1(H0SfSf(1,1,1,3,3):H0SfSf(3,1,1,3,3))

	common /coupdata/
     &    couplings, couplingsms, gammas, gammasms,
     &    ratios, chSt1St1

	RealType hggU(3,3), hggDRe(3,3), hggDIm(3,3)
	RealType hggSq(3), hgagaQ, hgagaSq

	common /kfactors/
     &    hggU, hggDRe, hggDIm, hggSq, hgagaQ, hgagaSq

	RealType hggU_flat(3*3)
	RealType hggDRe_flat(3*3), hggDIm_flat(3*3)
	equivalence (hggU, hggU_flat)
	equivalence (hggDRe, hggDRe_flat)
	equivalence (hggDIm, hggDIm_flat)

* flags

	integer mssmpart, fieldren, tanbren
	integer higgsmix, p2approx, looplevel
	integer runningMT, botResum, tlCplxApprox
	integer debuglevel, debugunit, fv
	integer uzint, uzext, mfeff
	integer tM1, tM2, bM, bM0, gM

* debuglevel = 0: no debug messages
*              1: dump setflags and setpara values
*              2: display Higgs mass matrix at p^2 = 0 and CTs
*              3: display search for zeros

	common /flags/
     &    mssmpart, fieldren, tanbren,
     &    higgsmix, p2approx, looplevel,
     &    runningMT, botResum, tlCplxApprox,
     &    debuglevel, debugunit, fv,
     &    uzint, uzext, mfeff,
     &    tM1, tM2, bM, bM0, gM


	integer flags_valid, sm_valid, para_valid, sf_valid
	integer tl_valid, higgs_valid, coup_valid, Ab_bad

	common /valids/
     &    flags_valid, sm_valid, para_valid, sf_valid,
     &    tl_valid, higgs_valid, coup_valid, Ab_bad


	character*1 cMSS(5), cAf(4)
	common /debug/ cMSS, cAf

