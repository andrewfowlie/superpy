* OutputSLHA.F
* fill SLHA data structure from FeynHiggs inputs and outputs
* this file is part of FeynHiggs
* last modified 1 Aug 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"
#include "version.h"


	subroutine FHOutputSLHA(error, slhadata, key)
	implicit none
	integer error, key

#include "FH.h"
#include "SLHA.h"
#include "PDG.h"
#include "looptools.h"
#define __SUBROUTINE__ "FHOutputSLHA"

	ComplexType slhadata(nslhadata)

	integer t, ty, nQ, nU, g, h, i, j, g1, g2
	integer vv, h1, h2, c1, c2, n1, n2, s1, s2, hv
	integer parent, ccb
	RealType x

	RealType DeltaMHiggs(NHiggs)
	ComplexType DeltaSAeff
	ComplexType DeltaUHiggs(NNeutral,NNeutral)
	ComplexType DeltaZHiggs(NNeutral,NNeutral)

	integer higgs(3), v_v(2,5), ferm(4,3), sferm(2,4,3)
	integer cha(2), neu(4)

#include "DeltaAlfaTop.h"

* Important: unlike the charged (s)fermions, PDG_Hp, PDG_W, and 
* PDG_chargino{1,2} refer to *positive* particles.  This differs
* from the FeynArts conventions adopted in FeynHiggs where
* particles are generally negative and antiparticles positive.

	data higgs / PDG_h0, PDG_HH, PDG_A0 /

	data v_v /
     &    PDG_photon, PDG_photon,
     &    PDG_photon, PDG_Z,
     &    PDG_Z, PDG_Z,
     &    PDG_W, -PDG_W,
     &    PDG_gluon, PDG_gluon /

	data ferm /
     &    PDG_nu_e, PDG_electron, PDG_up, PDG_down,
     &    PDG_nu_mu, PDG_muon, PDG_charm, PDG_strange,
     &    PDG_nu_tau, PDG_tau, PDG_top, PDG_bottom /

	data sferm /
     &    PDG_snu_e1, PDG_snu_e2, PDG_selectron1, PDG_selectron2,
     &    PDG_sup1, PDG_sup2, PDG_sdown1, PDG_sdown2,
     &    PDG_snu_mu1, PDG_snu_mu2, PDG_smuon1, PDG_smuon2,
     &    PDG_scharm1, PDG_scharm2, PDG_sstrange1, PDG_sstrange2,
     &    PDG_snu_tau1, PDG_snu_tau2, PDG_stau1, PDG_stau2,
     &    PDG_stop1, PDG_stop2, PDG_sbottom1, PDG_sbottom2 /

	data cha / -PDG_chargino1, -PDG_chargino2 /

	data neu / PDG_neutralino1, PDG_neutralino2,
     &    PDG_neutralino3, PDG_neutralino4 /


* R1, R2 (being nothing or 'c' = comment char) decide into how many
* blocks the SLHA data are redundantly written, this is important
* e.g. to keep SLHA1-only readers happy.

#define R1
#define R2


	error = 0
	CheckPara(error)
	CheckSf(error)
	if( error .ne. 0 ) return

	ModSel_Model = 0
	ModSel_GridPts = 1
	ModSel_Content = 0
	ModSel_RPV = 0
	if( higgsmix .lt. 3 ) then
	  ModSel_CPV = 0
	else
	  ModSel_CPV = 2
	endif
	ModSel_FV = ibits(-fv, 15, 1)

	SPInfo_NLines = 0
	DCInfo_NLines = 0
	call SLHAPutInfo(BlockSPInfo(1), 1, "FeynHiggs")
	call SLHAPutInfo(BlockSPInfo(1), 2, FHVERSION)
	call SLHAPutInfo(BlockSPInfo(1), 2, "built on "//FHBUILDDATE)

*** bit 0: RetrievePara

	if( .not. btest(key, 0) ) goto 1

	SMInputs_invAlfaMZ = invAlfaMZ - 1/Alfa0*DeltaAlfaTop(MT2, MZ2)
	SMInputs_AlfasMZ = AlfasMZ
	SMInputs_GF = GF

	SMInputs_Me = ME
	SMInputs_Mu = MU
	SMInputs_Md = MD
	SMInputs_Mmu = MM
	SMInputs_Mc = MC
	SMInputs_Ms = MS
	SMInputs_Mtau = ML
	SMInputs_Mt = MT
	SMInputs_Mb = MB

	SMInputs_MZ = MZ

	VCKMIn_lambda = CKMlambda
	VCKMIn_A = CKMA
	VCKMIn_rhobar = CKMrhobar
	VCKMIn_etabar = CKMetabar

	DMass_Q = scalefactor*MT

	MinPar_TB = TB
R2	ExtPar_TB = TB
R1	HMix_TB = TB

	Mass_MA0 = MA0
R2	ExtPar_MA0 = MA0
	Mass_MHp = MHp
R2	ExtPar_MHp = MHp

	HMix_MUE = MUE
R2	ExtPar_MUE = MUE

	MSoft_M1 = M_1
R2	ExtPar_M1 = M_1
	MSoft_M2 = M_2
R2	ExtPar_M2 = M_2
	MSoft_M3 = M_3
R1	Mass_MGl = MGl
R2	ExtPar_M3 = M_3

	do t = 2, 4
	  nQ = tQ(t)
R2	  ExtPar_MSS(1,nQ) = MSS0(nQ,1)
R1	  MSoft_MSS(1,nQ) = MSS0(nQ,1)
	  MSS2_MSS2(1,1,nQ) = SignSq(MSS0(nQ,1))
R2	  ExtPar_MSS(2,nQ) = MSS0(nQ,2)
R1	  MSoft_MSS(2,nQ) = MSS0(nQ,2)
	  MSS2_MSS2(2,2,nQ) = SignSq(MSS0(nQ,2))
R2	  ExtPar_MSS(3,nQ) = MSS0(nQ,3)
R1	  MSoft_MSS(3,nQ) = MSS0(nQ,3)
	  MSS2_MSS2(3,3,nQ) = SignSq(MSS0(nQ,3))
R2	  ExtPar_Q = QSf(t)
R1	  MSoft_Q = QSf(t)
	  MSS2_Q(nQ) = QSf(t)

	  nU = tU(t)
R2	  ExtPar_MSS(1,nU) = MSS0(nU,1)
R1	  MSoft_MSS(1,nU) = MSS0(nU,1)
	  MSS2_MSS2(1,1,nU) = SignSq(MSS0(nU,1))
R2	  ExtPar_MSS(2,nU) = MSS0(nU,2)
R1	  MSoft_MSS(2,nU) = MSS0(nU,2)
	  MSS2_MSS2(2,2,nU) = SignSq(MSS0(nU,2))
R2	  ExtPar_MSS(3,nU) = MSS0(nU,3)
R1	  MSoft_MSS(3,nU) = MSS0(nU,3)
	  MSS2_MSS2(3,3,nU) = SignSq(MSS0(nU,3))
	  MSS2_Q(nU) = QSf(t)

	  x = sqrt2/(vev*SCB(t))
	  ty = t + (bBR - 4)*isB(t)
R1	  Af_Af(1,1,t) = Af0(t,1)
	  Yf_Yf(1,1,t) = x*Mf(ty,1)
	  Tf_Tf(1,1,t) = Yf_Yf(1,1,t)*Af0(t,1)
R1	  Af_Af(2,2,t) = Af0(t,2)
	  Yf_Yf(2,2,t) = x*Mf(ty,2)
	  Tf_Tf(2,2,t) = Yf_Yf(2,2,t)*Af0(t,2)
R2	  ExtPar_Af(t) = Af0(t,3)
R1	  Af_Af(3,3,t) = Af0(t,3)
	  Yf_Yf(3,3,t) = x*Mf(ty,3)
	  Tf_Tf(3,3,t) = Yf_Yf(3,3,t)*Af0(t,3)
R1	  Af_Q(t) = QSf(t)
	  Yf_Q(t) = QSf(t)
	  Tf_Q(t) = QSf(t)

	  if( btest(fv, t) ) then
	    MSS2_MSS2(1,2,nQ) = deltaSf_LL(1,2,t)*MSS02_LL(1,2)
	    MSS2_MSS2(2,3,nQ) = deltaSf_LL(2,3,t)*MSS02_LL(2,3)
	    MSS2_MSS2(1,3,nQ) = deltaSf_LL(1,3,t)*MSS02_LL(1,3)

	    MSS2_MSS2(1,2,nU) = deltaSf_RR(1,2,t)*MSS02_RR(1,2)
	    MSS2_MSS2(2,3,nU) = deltaSf_RR(2,3,t)*MSS02_RR(2,3)
	    MSS2_MSS2(1,3,nU) = deltaSf_RR(1,3,t)*MSS02_RR(1,3)

	    Tf_Tf(1,2,t) = x*deltaSf_RL(1,2,t)*MSS02_RL(1,2)
	    Tf_Tf(2,3,t) = x*deltaSf_RL(2,3,t)*MSS02_RL(2,3)
	    Tf_Tf(1,3,t) = x*deltaSf_RL(1,3,t)*MSS02_RL(1,3)
	    Tf_Tf(2,1,t) = x*Conjugate(deltaSf_LR(1,2,t)*MSS02_LR(1,2))
	    Tf_Tf(3,2,t) = x*Conjugate(deltaSf_LR(2,3,t)*MSS02_LR(2,3))
	    Tf_Tf(3,1,t) = x*Conjugate(deltaSf_LR(1,3,t)*MSS02_LR(1,3))
	  endif
	enddo

*** bit 1: GetPara

1	if( .not. btest(key, 1) ) goto 2

	CheckSf(error)
	if( error .ne. 0 ) return

	Mass_MSf(1,1,1) = MSf(1,1,1)
	Mass_MSf(1,1,2) = MSf(1,1,2)
	Mass_MSf(1,1,3) = MSf(1,1,3)

	do t = 2, 4
	  ty = t + (bBR - 4)*isB(t)
	  if( btest(fv, t) ) then
	    Mass_MSf(1,t,1) = MASf(1,ty)
	    Mass_MSf(2,t,1) = MASf(2,ty)
	    Mass_MSf(1,t,2) = MASf(3,ty)
	    Mass_MSf(2,t,2) = MASf(4,ty)
	    Mass_MSf(1,t,3) = MASf(5,ty)
	    Mass_MSf(2,t,3) = MASf(6,ty)
	    do g2 = 1, 6
	      do g1 = 1, 6
	        ASfMix_UASf(g1,g2,t) = UASf(g1,g2,ty)
	      enddo
	    enddo
	  else
	    do g = 1, 3
	      Mass_MSf(1,t,g) = MSf(1,ty,g)
	      Mass_MSf(2,t,g) = MSf(2,ty,g)
	    enddo
* Note: the SfMix blocks exist only for g = 3
	    g = 3
	    SfMix_USf(1,1,t) = USf(1,1,ty,g)
	    SfMix_USf(2,1,t) = USf(2,1,ty,g)
	    SfMix_USf(1,2,t) = USf(1,2,ty,g)
	    SfMix_USf(2,2,t) = USf(2,2,ty,g)
	  endif
	enddo

	Mass_MCha(1) = MCha(1)
	Mass_MCha(2) = MCha(2)

	UMix_UCha(1,1) = UCha(1,1)
	UMix_UCha(2,1) = UCha(2,1)
	UMix_UCha(1,2) = UCha(1,2)
	UMix_UCha(2,2) = UCha(2,2)

	VMix_VCha(1,1) = VCha(1,1)
	VMix_VCha(2,1) = VCha(2,1)
	VMix_VCha(1,2) = VCha(1,2)
	VMix_VCha(2,2) = VCha(2,2)

	do i = 1, 4
	  Mass_MNeu(i) = MNeu(i)
	  if( Im(ZNeu(i,1)) .eq. 0 ) then
	    do j = 1, 4
	      NMix_ZNeu(i,j) = ZNeu(i,j)
	    enddo
	  else
	    Mass_MNeu(i) = -MNeu(i)
	    do j = 1, 4
	      NMix_ZNeu(i,j) = cI*ZNeu(i,j)
	    enddo
	  endif
	enddo

	Mass_MGl = MGl

*** bit 2: HiggsCorr

2	if( .not. btest(key, 2) ) goto 3

	call CheckHiggsCorr(error)
	if( error .ne. 0 ) return

	Mass_Mh0 = MHiggs(1)
	Mass_MHH = MHiggs(2)
	Mass_MA0 = MHiggs(3)
	Mass_MHp = MHiggs(4)

	if( SAeff .eq. -2 ) then
	  Alpha_Alpha = invalid
	else
	  Alpha_Alpha = asin(Re(SAeff))
	endif

	do i = 1, 3
	  do j = 1, 3
	    CVHMix_UH(i,j) = UHiggs(i,j)
	  enddo
	enddo

*** bit 3: Uncertainties

3	if( .not. btest(key, 3) ) goto 4

	call FHUncertainties(error,
     &    DeltaMHiggs, DeltaSAeff, DeltaUHiggs, DeltaZHiggs)
	if( error .ne. 0 ) return

	DMass_DeltaMh0 = DeltaMHiggs(1)
	DMass_DeltaMHH = DeltaMHiggs(2)
	DMass_DeltaMA0 = DeltaMHiggs(3)
	DMass_DeltaMHp = DeltaMHiggs(4)

	if( DeltaSAeff .eq. -2 ) then
	  DAlpha_DeltaAlpha = invalid
	else
	  DAlpha_DeltaAlpha = .5D0*abs(
     &      asin(min(Re(SAeff) + Re(DeltaSAeff), 1D0)) -
     &      asin(max(Re(SAeff) - Re(DeltaSAeff), -1D0)) )
	endif

*** bit 4+5: Couplings (bit 5 = fast)

4	if( .not. btest(key, 4) ) goto 6

	call CheckCouplings(error, ibits(key, 5, 1))
	if( error .ne. 0 ) return

	do h = 1, 3
	  parent = SLHANewDecay(slhadata, GammaTot(h), higgs(h))

	  do vv = 1, 5
	    call SLHAAddDecay(slhadata, BR(H0VV(h,vv)),
     &        parent, 2, v_v(1,vv), v_v(2,vv))
	  enddo

	  do t = 2, 4
	    do g2 = 1, 3
	      do g1 = 1, 3
	        call SLHAAddDecay(slhadata, BR(H0FF(h,t,g1,g2)),
     &            parent, 2, -ferm(t,g1), ferm(t,g2))
	      enddo
	    enddo
	  enddo

	  do c1 = 1, 2
	    do c2 = 1, 2
	      call SLHAAddDecay(slhadata, BR(H0ChaCha(h,c1,c2)),
     &          parent, 2, -cha(c1), cha(c2))
	    enddo
	  enddo

	  do n1 = 1, 4
	    do n2 = n1, 4
	      call SLHAAddDecay(slhadata, BR(H0NeuNeu(h,n1,n2)),
     &          parent, 2, neu(n1), neu(n2))
	    enddo
	  enddo

	  do hv = 1, 3
	    call SLHAAddDecay(slhadata, BR(H0HV(h,hv)),
     &        parent, 2, higgs(hv), PDG_Z)
	  enddo

	  do h1 = 1, 4
	    do h2 = 1, 4
	      call SLHAAddDecay(slhadata, BR(H0HH(h,h1,h2)),
     &          parent, 2, higgs(h1), higgs(h2))
	    enddo
	  enddo

	  do t = 2, 4
	    do g = 1, 3
	      do s1 = 1, 2
	        do s2 = 1, 2
	          call SLHAAddDecay(slhadata, BR(H0SfSf(h,s1,s2,t,g)),
     &              parent, 2, -sferm(s1,t,g), sferm(s2,t,g))
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

	parent = SLHANewDecay(slhadata, GammaTot(4), PDG_Hp)

	do t = 2, 4, 2
	  do g1 = 1, 3
	    do g2 = 1, 3
	      call SLHAAddDecay(slhadata, BR(HpFF(t/2,g1,g2)),
     &          parent, 2, ferm(t-1,g1), -ferm(t,g2))
	    enddo
	  enddo
	enddo

	do n1 = 1, 4
	  do c2 = 1, 2
	    call SLHAAddDecay(slhadata, BR(HpNeuCha(n1,c2)),
     &        parent, 2, neu(n1), -cha(c2))
	  enddo
	enddo

	do hv = 1, 3
	  call SLHAAddDecay(slhadata, BR(HpHV(hv)),
     &      parent, 2, higgs(hv), PDG_W)
	enddo

	do t = 2, 4, 2
	  do g1 = 1, 3
	    do g2 = 1, 3
	      do s1 = 1, 2
	        do s2 = 1, 2
	          call SLHAAddDecay(slhadata, BR(HpSfSf(s1,s2,t/2,g1,g2)),
     &              parent, 2, sferm(s1,t-1,g1), -sferm(s2,t,g2))
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

	parent = SLHANewDecay(slhadata,
     &    Gamma(tBF(1)) + Gamma(tBF(2)), PDG_top)

	call SLHAAddDecay(slhadata, BR(tBF(1)),
     &    parent, 2, PDG_W, PDG_bottom)
	call SLHAAddDecay(slhadata, BR(tBF(2)),
     &    parent, 2, PDG_Hp, PDG_bottom)

*** bit 6: Constraints

6	if( .not. btest(key, 6) ) goto 7

	call FHConstraints(error,
     &    PrecObs_gminus2mu,
     &    PrecObs_DeltaRho,
     &    PrecObs_MWMSSM, PrecObs_MWSM,
     &    PrecObs_SW2effMSSM, PrecObs_SW2effSM,
     &    PrecObs_EDMeTh, PrecObs_EDMn, PrecObs_EDMHg,
     &    ccb)
	if( error .ne. 0 ) return

	if( ccb .ne. 0 ) Warning("Colour-breaking minimum detected")

*** bit 7: Flavour

7	if( .not. btest(key, 7) ) return

	call FHFlavour(error,
     &    PrecObs_bsgammaMSSM, PrecObs_bsgammaSM,
     &    PrecObs_DeltaMsMSSM, PrecObs_DeltaMsSM,
     &    PrecObs_BsmumuMSSM, PrecObs_BsmumuSM)
	end

