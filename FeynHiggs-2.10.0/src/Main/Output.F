* Output.F
* run FH functions encoded in key and write output to unit
* this file is part of FeynHiggs
* last modified 16 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"

#define MASS_FMT	F16.8
#define TRIG_FMT	F14.8
#define GAMMA_FMT	G15.6

#define TAG		"| "
#define SKIP		"%"
#define CONT		" \\"/TAG,T19
#define EQU		T17,"= "
#define GAMMA_EQU	T25,"= "


	subroutine FHOutput(error, filename, key, sqrts)
	implicit none
	integer error, key
	character*(*) filename
	RealType sqrts

#include "FH.h"
#define __SUBROUTINE__ "FHOutput"

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer unit, i, j, h, hv, h1, h2, t, vv, p
	integer g, g1, g2, c1, c2, n1, n2, s1, s2

	RealType DeltaMHiggs(NHiggs)
	ComplexType DeltaSAeff
	ComplexType DeltaUHiggs(NNeutral,NNeutral)
	ComplexType DeltaZHiggs(NNeutral,NNeutral)

	RealType gm2
	RealType DeltaRho, MWMSSM, MWSM, SW2MSSM, SW2SM
	RealType edmeTh, edmn, edmHg
	RealType bsgMSSM, bsgSM
	RealType deltaMsMSSM, deltaMsSM
	RealType bsmumuMSSM, bsmumuSM
	integer ccb

	RealType prodxs(nprodxs)

	character*2 higgs(4)

	character*7 f_f(2:4,3)
	data f_f /
     &    "e-e", "u-u", "d-d",
     &    "mu-mu", "c-c", "s-s",
     &    "tau-tau", "t-t", "b-b" /

	character*10 f_fprime(2,3,3)
	data f_fprime /
     &    "nu_e-e",     "u-d",
     &    "nu_mu-e",    "c-d",
     &    "nu_tau-e",   "t-d",
     &    "nu_e-mu",    "u-s",
     &    "nu_mu-mu",   "c-s",
     &    "nu_tau-mu",  "t-s",
     &    "nu_e-tau",   "u-b",
     &    "nu_mu-tau",  "c-b",
     &    "nu_tau-tau", "t-b" /

	character*6 sfe(4,3)
	data sfe /
     &    "snue", "se", "su", "sd",
     &    "snumu", "smu", "sc", "ss",
     &    "snutau", "stau", "stop", "sbot" /

	character*12 v_v(5)
	data v_v / "gamma-gamma", "gamma-Z", "Z-Z", "W-W", "g-g" /


	if( filename .eq. "-" .or. filename .eq. "stdout" ) then
	  unit = 6
	else if( filename .eq. "stderr" ) then
	  unit = 0
	else
	  unit = 10
	  open(unit, file=filename, status="unknown", err=999)
	endif

	if( higgsmix .lt. 3 ) then
	  higgs(1) = "h0"
	  higgs(2) = "HH"
	  higgs(3) = "A0"
	else
	  higgs(1) = "h1"
	  higgs(2) = "h2"
	  higgs(3) = "h3"
	endif
	higgs(4) = "Hp"

	error = 0
	CheckPara(error)

*** bit 0: RetrievePara

	if( .not. btest(key, 0) ) goto 1

100	format(/"------------------- INPUT PARAMETERS -------------------")
	write(unit,100)

101	format(TAG, A, EQU, 4 MASS_FMT)

	write(unit,101) "invAlfaMZ", invAlfaMZ
	write(unit,101) "AlfasMZ", AlfasMZ
	write(unit,101) "GF", GF

	write(unit,101) "ME", ME
	write(unit,101) "MU", MU
	write(unit,101) "MD", MD
	write(unit,101) "MM", MM
	write(unit,101) "MC", MC
	write(unit,101) "MS", MS
	write(unit,101) "ML", ML
	write(unit,101) "MT", MT
	write(unit,101) "MB", MB

	write(unit,101) "MW", MW
	write(unit,101) "MZ", MZ

	write(unit,101) "CKMlambda", CKMlambda
	write(unit,101) "CKMA", CKMA
	write(unit,101) "CKMrhobar", CKMrhobar
	write(unit,101) "CKMetabar", CKMetabar

	write(unit,101) "scalefactor", scalefactor

	write(unit,101) "TB", TB
	write(unit,101) "MA0", MA0
	write(unit,101) "MHp", MHp

	write(unit,101) "MUE", MUE

	write(unit,101) "M_1", M_1
	write(unit,101) "M_2", M_2
	write(unit,101) "M_3", M_3

	write(unit,101) "M1SL", MSS(1,1)
	write(unit,101) "M1SE", MSS(2,1)
	write(unit,101) "M1SQ", MSS(3,1)
	write(unit,101) "M1SU", MSS(4,1)
	write(unit,101) "M1SD", MSS(5,1)
	write(unit,101) "Ae", Af(2,1)
	write(unit,101) "Au", Af(3,1)
	write(unit,101) "Ad", Af(4,1)

	write(unit,101) "M2SL", MSS(1,2)
	write(unit,101) "M2SE", MSS(2,2)
	write(unit,101) "M2SQ", MSS(3,2)
	write(unit,101) "M2SU", MSS(4,2)
	write(unit,101) "M2SD", MSS(5,2)
	write(unit,101) "Amu", Af(2,2)
	write(unit,101) "Ac", Af(3,2)
	write(unit,101) "As", Af(4,2)

	write(unit,101) "M3SL", MSS(1,3)
	write(unit,101) "M3SE", MSS(2,3)
	write(unit,101) "M3SQ", MSS(3,3)
	write(unit,101) "M3SU", MSS(4,3)
	write(unit,101) "M3SD", MSS(5,3)
	write(unit,101) "Atau", Af(2,3)
	write(unit,101) "At", Af(3,3)
	write(unit,101) "Ab", Af(4,3)

	write(unit,101) "deltaQLL12", deltaSf_LL(1,2,3)
	write(unit,101) "deltaQLL23", deltaSf_LL(2,3,3)
	write(unit,101) "deltaQLL13", deltaSf_LL(1,3,3)

	write(unit,101) "deltaULR12", deltaSf_LR(1,2,3)
	write(unit,101) "deltaULR23", deltaSf_LR(2,3,3)
	write(unit,101) "deltaULR13", deltaSf_LR(1,3,3)

	write(unit,101) "deltaURL12", deltaSf_RL(1,2,3)
	write(unit,101) "deltaURL23", deltaSf_RL(2,3,3)
	write(unit,101) "deltaURL13", deltaSf_RL(1,3,3)

	write(unit,101) "deltaURR12", deltaSf_RR(1,2,3)
	write(unit,101) "deltaURR23", deltaSf_RR(2,3,3)
	write(unit,101) "deltaURR13", deltaSf_RR(1,3,3)

	write(unit,101) "deltaDLR12", deltaSf_LR(1,2,4)
	write(unit,101) "deltaDLR23", deltaSf_LR(2,3,4)
	write(unit,101) "deltaDLR13", deltaSf_LR(1,3,4)

	write(unit,101) "deltaDRL12", deltaSf_RL(1,2,4)
	write(unit,101) "deltaDRL23", deltaSf_RL(2,3,4)
	write(unit,101) "deltaDRL13", deltaSf_RL(1,3,4)

	write(unit,101) "deltaDRR12", deltaSf_RR(1,2,4)
	write(unit,101) "deltaDRR23", deltaSf_RR(2,3,4)
	write(unit,101) "deltaDRR13", deltaSf_RR(1,3,4)

*** bit 1: GetPara

1	if( .not. btest(key, 1) ) goto 2

110	format(/"------------------ DERIVED PARAMETERS ------------------")
	write(unit,110)

	CheckSf(error)
	if( error .ne. 0 ) goto 500

111	format(TAG, A, EQU, 2 TRIG_FMT)

	write(unit,101) "MSnu_e", MSf(1,1,1)
	write(unit,101) "MSnu_mu", MSf(1,1,2)
	write(unit,101) "MSnu_tau", MSf(1,1,3)
	write(unit,101) "MSe", MSf(1,2,1), MSf(2,2,1)
	write(unit,101) "MSmu", MSf(1,2,2), MSf(2,2,2)
	write(unit,101) "MStau", MSf(1,2,3), MSf(2,2,3)
	write(unit,111) "cstau", USf(1,1,2,3)
	write(unit,111) "sstau", USf(1,2,2,3)
	write(unit,101) "MSup", MSf(1,3,1), MSf(2,3,1)
	write(unit,101) "MScharm", MSf(1,3,2), MSf(2,3,2)
	write(unit,101) "MStop", MSf(1,3,3), MSf(2,3,3)
	write(unit,111) "cstop", USf(1,1,3,3)
	write(unit,111) "sstop", USf(1,2,3,3)
	write(unit,101) "MSdown", MSf(1,bBR,1), MSf(2,bBR,1)
	write(unit,101) "MSstrange", MSf(1,bBR,2), MSf(2,bBR,2)
	write(unit,101) "MSbot", MSf(1,bBR,3), MSf(2,bBR,3)
	write(unit,111) "csbot", USf(1,1,bBR,3)
	write(unit,111) "ssbot", USf(1,2,bBR,3)

	write(unit,101) "MCha", MCha
	write(unit,101) "MNeu", MNeu
	write(unit,101) "Deltab", Deltaf(bTR,3)
	write(unit,101) "MGl", MGl

*** bit 2: HiggsCorr

2	if( .not. btest(key, 2) ) goto 3

120	format(/"--------------------- HIGGS MASSES ---------------------")
	write(unit,120)

	call CheckHiggsCorr(error)
	if( error .ne. 0 ) goto 500

121	format(TAG, A,A, EQU, MASS_FMT)
	write(unit,121) ("M", higgs(h), MHiggs(h), h = 1, 4)

	write(unit,111) "SAeff", SAeff

122	format(TAG, A, EQU, 3 TRIG_FMT,
     &                CONT, 3 TRIG_FMT,
     &                CONT, 3 TRIG_FMT)
	write(unit,122) "UHiggsRe",
     &    ((Re(UHiggs(i,j)), j = 1, 3), i = 1, 3)
	write(unit,122) "UHiggsIm",
     &    ((Im(UHiggs(i,j)), j = 1, 3), i = 1, 3)

	write(unit,122) "ZHiggsRe",
     &    ((Re(ZHiggs(i,j)), j = 1, 3), i = 1, 3)
	write(unit,122) "ZHiggsIm",
     &    ((Im(ZHiggs(i,j)), j = 1, 3), i = 1, 3)

*** bit 3: Uncertainties

3	if( .not. btest(key, 3) ) goto 4

130	format(/"---------------- ESTIMATED UNCERTAINTIES ---------------")
	write(unit,130)

	call FHUncertainties(error,
     &    DeltaMHiggs, DeltaSAeff, DeltaUHiggs, DeltaZHiggs)
	if( error .ne. 0 ) goto 500

	print 121, ("DeltaM", higgs(h), DeltaMHiggs(h), h = 1, 4)

	write(unit,111) "DeltaSAeff", DeltaSAeff

	write(unit,122) "DeltaUHiggsRe",
     &    ((Re(DeltaUHiggs(i,j)), j = 1, 3), i = 1, 3)
	write(unit,122) "DeltaUHiggsIm",
     &    ((Im(DeltaUHiggs(i,j)), j = 1, 3), i = 1, 3)

	write(unit,122) "DeltaZHiggsRe",
     &    ((Re(DeltaZHiggs(i,j)), j = 1, 3), i = 1, 3)
	write(unit,122) "DeltaZHiggsIm",
     &    ((Im(DeltaZHiggs(i,j)), j = 1, 3), i = 1, 3)

*** bit 4+5: Couplings (bit 5 = fast)

4	if( .not. btest(key, 4) ) goto 6

140	format(/"---------- DECAY WIDTHS AND BRANCHING RATIOS -----------")
	write(unit,140)

	call CheckCouplings(error, ibits(key, 5, 1))
	if( error .ne. 0 ) goto 500

141	format(TAG, "GammaTot-",A, GAMMA_EQU, 1P, GAMMA_FMT)
	write(unit,141) (higgs(h), GammaTot(h), h = 1, 4)
	write(unit,141) ("SM"//higgs(h), GammaSMTot(h), h = 1, 3)

#define GAMMA_BR_FMT GAMMA_EQU, 1P, 3 GAMMA_FMT

142	format(SKIP/SKIP,"  Channel", T29,A, T44,A,:,T59,A)
143	format(SKIP, TAG, A,"-",A, GAMMA_BR_FMT)
144	format(SKIP, TAG, A,"-",A,I1,"-",A,I1, GAMMA_BR_FMT)
145	format(SKIP, TAG, A,"-",A,"-",A, GAMMA_BR_FMT)

	do h = 1, 3
	  write(unit,142) "Gamma", "BR", "BRSM"

	  do vv = 1, 5
	    if( Gamma(H0VV(h,vv)) .ne. 0 .or.
     &          GammaSM(H0VV(h,vv)) .ne. 0 )
     &        write(unit,143) higgs(h), v_v(vv),
     &          Gamma(H0VV(h,vv)),
     &          BR(H0VV(h,vv)),
     &          BRSM(H0VV(h,vv))
	  enddo

	  do t = 2, 4
	    do g = 1, 3
	      if( Gamma(H0FF(h,t,g,g)) .ne. 0 .or.
     &            GammaSM(H0FF(h,t,g,g)) .ne. 0 )
     &          write(unit,143) higgs(h), f_f(t,g),
     &            Gamma(H0FF(h,t,g,g)),
     &            BR(H0FF(h,t,g,g)),
     &            BRSM(H0FF(h,t,g,g))
	    enddo
	  enddo

	  do c1 = 1, 2
	    do c2 = 1, 2
	      if( Gamma(H0ChaCha(h,c1,c2)) .ne. 0 )
     &          write(unit,144) higgs(h), "Cha",c1, "Cha",c2,
     &            Gamma(H0ChaCha(h,c1,c2)),
     &            BR(H0ChaCha(h,c1,c2))
	    enddo
	  enddo

	  do n1 = 1, 4
	    do n2 = n1, 4
	      if( Gamma(H0NeuNeu(h,n1,n2)) .ne. 0 )
     &          write(unit,144) higgs(h), "Neu",n1, "Neu",n2,
     &            Gamma(H0NeuNeu(h,n1,n2)),
     &            BR(H0NeuNeu(h,n1,n2))
	    enddo
	  enddo

	  do hv = 1, 3
	    if( Gamma(H0HV(h,hv)) .ne. 0 )
     &        write(unit,145) higgs(h), higgs(hv), "Z",
     &          Gamma(H0HV(h,hv)),
     &          BR(H0HV(h,hv))
	  enddo

	  do t = 1, 4
	    do g = 1, 3
	      do s1 = 1, 2
	        do s2 = s1, 2
	          if( Gamma(H0SfSf(h,s1,s2,t,g)) .ne. 0 )
     &              write(unit,144)
     &                higgs(h), Strip(sfe(t,g)),s1, Strip(sfe(t,g)),s2,
     &                Gamma(H0SfSf(h,s1,s2,t,g)),
     &                BR(H0SfSf(h,s1,s2,t,g))
	        enddo
	      enddo
	    enddo
	  enddo

	  do h1 = 1, 4
	    do h2 = h1, 4
	      if( Gamma(H0HH(h,h1,h2)) .ne. 0 )
     &          write(unit,145) higgs(h), higgs(h1), higgs(h2),
     &            Gamma(H0HH(h,h1,h2)), BR(H0HH(h,h1,h2))
	    enddo
	  enddo
	enddo

	write(unit,142) "Gamma", "BR"

	do g2 = 1, 3
	  do g1 = 1, 3
	    do p = 1, 2
	      if( Gamma(HpFF(p,g1,g2)) .ne. 0 )
     &          write(unit,143) higgs(4), f_fprime(p,g1,g2),
     &            Gamma(HpFF(p,g1,g2)),
     &            BR(HpFF(p,g1,g2))
	    enddo
	  enddo
	enddo

	do c2 = 1, 2
	  do n1 = 1, 4
	    if( Gamma(HpNeuCha(n1,c2)) .ne. 0 )
     &        write(unit,144) higgs(4), "Neu",n1, "Cha",c2,
     &          Gamma(HpNeuCha(n1,c2)),
     &          BR(HpNeuCha(n1,c2))
	  enddo
	enddo

	do hv = 1, 3
	  if( Gamma(HpHV(hv)) .ne. 0 )
     &      write(unit,145) higgs(4), higgs(hv), "W",
     &        Gamma(HpHV(hv)),
     &        BR(HpHV(hv))
	enddo

	do t = 2, 4, 2
	  do g2 = 1, 3
	    do g1 = 1, 3
	      do s1 = 1, 2
	        do s2 = 1, 2
	          if( Gamma(HpSfSf(s1,s2,t/2,g1,g2)) .ne. 0 )
     &              write(unit,144)
     &                higgs(4), Strip(sfe(t-1,g1)),s1, Strip(sfe(t,g2)),s2,
     &                Gamma(HpSfSf(s1,s2,t/2,g1,g2)),
     &                BR(HpSfSf(s1,s2,t/2,g1,g2))
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

	write(unit,142) "Gamma", "BR"

	write(unit,145) "t", "W", "b", Gamma(tBF(1)), BR(tBF(1))
	write(unit,145) "t", higgs(4), "b", Gamma(tBF(2)), BR(tBF(2))

150	format(SKIP/SKIP,
     &    "------------------ SELECTED COUPLINGS ------------------")
	write(unit,150)

#define COUPLING_FMT GAMMA_EQU, 1P, 3 GAMMA_FMT

151	format(SKIP, TAG, A,":",A,"-",A, COUPLING_FMT)
152	format(SKIP, TAG, A,":",A,"-",A,I1,"-",A,I1, COUPLING_FMT)
153	format(SKIP, TAG, A,":",A,"-",A,"-",A, COUPLING_FMT)

#define SMRatio(c,mssm,sm) Im(mssm(c))/Im(sm(c))

	do h = 1, 3
	  write(unit,142) "Re Coupling", "Im Coupling", "SM Ratio"

	  if( Coupling(H0VV(h,3)) .ne. 0 )
     &      write(unit,151) "C", higgs(h), v_v(3),
     &        Coupling(H0VV(h,3)),
     &        SMRatio(H0VV(h,3), Coupling, CouplingSM)

	  if( Coupling(H0VV(h,4)) .ne. 0 )
     &      write(unit,151) "C", higgs(h), v_v(4),
     &        Coupling(H0VV(h,4)),
     &        SMRatio(H0VV(h,4), Coupling, CouplingSM)

	  do t = 2, 4
	    do g = 3, 3
	      if( LCoupling(H0FF(h,t,g,g)) .ne. 0 .or.
     &            RCoupling(H0FF(h,t,g,g)) .ne. 0 ) then
	        write(unit,151) "CL", higgs(h), f_f(t,g),
     &            LCoupling(H0FF(h,t,g,g)),
     &            SMRatio(H0FF(h,t,g,g), LCoupling, LCouplingSM)
	        write(unit,151) "CR", higgs(h), f_f(t,g),
     &            RCoupling(H0FF(h,t,g,g)),
     &            SMRatio(H0FF(h,t,g,g), RCoupling, RCouplingSM)
	      endif
	    enddo
	  enddo

	  do c1 = 1, 2
	    do c2 = 1, 2
	      if( LCoupling(H0ChaCha(h,c1,c2)) .ne. 0 .or.
     &            RCoupling(H0ChaCha(h,c1,c2)) .ne. 0 ) then
	        write(unit,152) "CL", higgs(h), "Cha",c1, "Cha",c2,
     &            LCoupling(H0ChaCha(h,c1,c2))
	        write(unit,152) "CR", higgs(h), "Cha",c1, "Cha",c2,
     &            RCoupling(H0ChaCha(h,c1,c2))
	      endif
	    enddo
	  enddo

	  do n1 = 1, 4
	    do n2 = n1, 4
	      if( LCoupling(H0NeuNeu(h,n1,n2)) .ne. 0 .or.
     &            RCoupling(H0NeuNeu(h,n1,n2)) .ne. 0 ) then
	        write(unit,152) "CL", higgs(h), "Neu",n1, "Neu",n2,
     &            LCoupling(H0NeuNeu(h,n1,n2))
	        write(unit,152) "CR", higgs(h), "Neu",n1, "Neu",n2,
     &            RCoupling(H0NeuNeu(h,n1,n2))
	      endif
	    enddo
	  enddo

	  do hv = 1, 3
	    if( Coupling(H0HV(h,hv)) .ne. 0 )
     &        write(unit,153) "C", higgs(h), higgs(hv), "Z",
     &          Coupling(H0HV(h,hv))
	  enddo

	  do t = 1, 4
	    do g = 3, 3
	      do s1 = 1, 2
	        do s2 = s1, 2
	          if( Coupling(H0SfSf(h,s1,s2,t,g)) .ne. 0 )
     &              write(unit,152) "C", higgs(h),
     &                Strip(sfe(t,g)),s1, Strip(sfe(t,g)),s2,
     &                Coupling(H0SfSf(h,s1,s2,t,g))
	        enddo
	      enddo
	    enddo
	  enddo

	  do h1 = 1, 4
	    do h2 = h1, 4
	      if( Coupling(H0HH(h,h1,h2)) .ne. 0 )
     &          write(unit,153) "C", higgs(h), higgs(h1), higgs(h2),
     &            Coupling(H0HH(h,h1,h2))
	    enddo
	  enddo
	enddo

	write(unit,142) "Re Coupling", "Im Coupling"

	do g = 3, 3
	  do p = 1, 2
	    if( LCoupling(HpFF(p,g,g)) .ne. 0 .or.
     &          RCoupling(HpFF(p,g,g)) .ne. 0 ) then
	      write(unit,151) "CL", higgs(4), f_fprime(p,g,g),
     &          LCoupling(HpFF(p,g,g))
	      write(unit,151) "CR", higgs(4), f_fprime(p,g,g),
     &          RCoupling(HpFF(p,g,g))
	    endif
	  enddo
	enddo

	do c2 = 1, 2
	  do n1 = 1, 4
	    if( LCoupling(HpNeuCha(n1,c2)) .ne. 0 .or.
     &          RCoupling(HpNeuCha(n1,c2)) .ne. 0 ) then
	      write(unit,152) "CL", higgs(4), "Neu",n1, "Cha",c2,
     &          LCoupling(HpNeuCha(n1,c2))
	      write(unit,152) "CR", higgs(4), "Neu",n1, "Cha",c2,
     &          RCoupling(HpNeuCha(n1,c2))
	    endif
	  enddo
	enddo

	do hv = 1, 3
	  if( Coupling(HpHV(hv)) .ne. 0 )
     &      write(unit,153) "C", higgs(4), higgs(hv), "W",
     &        Coupling(HpHV(hv))
	enddo

	do t = 2, 4, 2
	  do g = 3, 3
	    do s1 = 1, 2
	      do s2 = 1, 2
	        if( Coupling(HpSfSf(s1,s2,t/2,g,g)) .ne. 0 )
     &            write(unit,152) "C", higgs(4),
     &              Strip(sfe(t-1,g)),s1, Strip(sfe(t,g)),s2,
     &              Coupling(HpSfSf(s1,s2,t/2,g,g))
	      enddo
	    enddo
	  enddo
	enddo

*** bit 6: Constraints

6	if( .not. btest(key, 6) ) goto 7

160	format(/"-------------------- EW CONSTRAINTS --------------------"
     &    /"Note: MW, SW2, and the EDMs are not fully checked yet!")
	write(unit,160)

	call FHConstraints(error, gm2,
     &    DeltaRho, MWMSSM, MWSM, SW2MSSM, SW2SM,
     &    edmeTh, edmn, edmHg, ccb)
	if( error .ne. 0 ) goto 500

161	format(TAG, A, EQU, 1P, GAMMA_FMT)
	write(unit,161) "DeltaRho", DeltaRho
	write(unit,161) "MWMSSM", MWMSSM
	write(unit,161) "MWSM", MWSM
	write(unit,161) "SW2effMSSM", SW2MSSM
	write(unit,161) "SW2effSM", SW2SM
	write(unit,161) "g-2", gm2
	write(unit,161) "EDMeTh", edmeTh
	write(unit,161) "EDMn", edmn
	write(unit,161) "EDMHg", edmHg

	if( ccb .ne. 0 ) Warning("Colour-breaking minimum detected")

*** bit 7: Flavour

7	if( .not. btest(key, 7) ) goto 8

170	format(/"----------------- FLAVOUR OBSERVABLES ------------------")
	write(unit,170)

	call FHFlavour(error,
     &    bsgMSSM, bsgSM,
     &    deltaMsMSSM, deltaMsSM,
     &    bsmumuMSSM, bsmumuSM)
	if( error .ne. 0 ) goto 500

	write(unit,161) "bsgamma", bsgMSSM
	write(unit,161) "bsgammaSM", bsgSM
	write(unit,161) "deltaMs", deltaMsMSSM
	write(unit,161) "deltaMsSM", deltaMsSM
	write(unit,161) "bsmumu", bsmumuMSSM
	write(unit,161) "bsmumuSM", bsmumuSM

*** Production XS

8	if( sqrts .le. 0 ) goto 500

180	format(/"-------------- PRODUCTION CROSS-SECTIONS ---------------")
	print 180

	call FHHiggsProd(error, sqrts, prodxs)
	if( error .ne. 0 ) goto 500

181	format(TAG, A, GAMMA_BR_FMT)
182	format(/"  Channel", T29,A, T44,A,:,T59,A)
183	format(TAG, A,"-",A, GAMMA_BR_FMT)

	write(unit,181) "prod:sqrts", sqrts

	write(unit,182) "xsection/fb", "SM xsection/fb"

	do h = 1, 3
	  write(unit,183) "prod:b-b",higgs(h), bbh(h), bbhSM(h)
	  write(unit,183) "prod:btag-b",higgs(h), btagbh(h), btagbhSM(h)
	  write(unit,183) "prod:g-g",higgs(h), ggh(h), gghSM(h)
	  write(unit,183) "prod:alt-g-g",higgs(h), ggh2(h), gghSM(h)
	  write(unit,183) "prod:q-q",higgs(h), qqh(h), qqhSM(h)
	  write(unit,183) "prod:t-t",higgs(h), tth(h), tthSM(h)
	  write(unit,183) "prod:W",  higgs(h), Wh(h), WhSM(h)
	  write(unit,183) "prod:Z",  higgs(h), Zh(h), ZhSM(h)
	  write(unit,183) "prod:St-St", higgs(h), StSth(h)
	enddo

	write(unit,182) "xsection/fb"
	write(unit,183) "prod:t", higgs(4), tHm
	write(unit,183) "prod:alt-t", higgs(4), tHm2
	write(unit,183) "prod:altlo-t", higgs(4), tHm2lo
	write(unit,183) "prod:althi-t", higgs(4), tHm2hi

500	if( unit .eq. 10 ) close(unit)
	return

999	Error(error, "cannot open file for writing")
	goto 500
	end

