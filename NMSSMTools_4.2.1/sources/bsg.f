*      Subroutine for BR(b -> s gamma)


* Literature:
*   - THEORETICAL FORMULAE:
* [1] K.Chetyrkin, M.Misiak, M.Munz,
*     'Weak Radiative B-Meson Decay Beyond Leading Logarithms'
*     Phys.Lett.B400:206-219,1997,
*     Erratum-ibid.B425:414,1998, e-Print: hep-ph/9612313
*
* [2] M.Ciuchini, G.Degrassi, P.Gambino, G.Giudice
*     'Next to Leading QCD Corrections to B -> Xs gamma: Standard Model
*     and Two-Higgs Doublet Model' Nucl.Phys.B534:3-20,1998,
*     e-Print: hep-ph/9806308
*
* [3] G.Degrassi, P.Gambino, G.Giudice
*     'B -> Xs gamma in Supersymmetry: Large Contributions Beyond the
*     Leading Order', JHEP 0012:009,2000, e-Print: hep-ph/0009337
*
* [4] P.Gambino, M.Misiak, 'Quark Mass Effects in B -> Xs gamma'
*     Nucl.Phys.B611:338-366,2001, e-Print: hep-ph/0104034
*
* [5] A.Buras, A.Czarnecki, M.Misiak, J.Urban,
*     'Completing the NLO QCD Calculation of B -> Xs gamma'
*     Nucl.Phys.B631:219-238,2002, e-Print: hep-ph/0203135
*
* [6] G.Belanger, F.Boudjema, A.Pukhov, A.Semenov,
*     'micrOMEGAs: Version 1.3', Appendix B,
*     .omput.Phys.Commun.174:577-604,2006, e-Print: hep-ph/0405253
*
* [7] T.Hurth, E.Lunghi, W.Porod,
*     'Untagged B -> Xs+d gamma CP asymmetry as a probe for New Physics'
*     Nucl.Phys.B704:56-74,2005, e-Print: hep-ph/0312260
*
* [8] M.Misiak et al.,
*     'Estimate of B(anti-B ---> X(s) gamma) at O(alpha(s)**2)'
*     Phys.Rev.Lett.98:022002,2007, e-Print: hep-ph/0609232
*
* [9] T.Becher and M.~Neubert,
*     Phys.\ Rev.\ Lett.\  {\bf 98} (2007) 022003 [arXiv:hep-ph/0610067].
*
* [10] A. J. Buras, P. H. Chankowski, J. Rosiek, L. Slawianowska,
*     ' Delta M(d, s), B0(d, s) ---> mu+ mu- and B ---> X(s) gamma
*      in supersymmetry at large tan beta.'
*     Nucl.Phys.B659:3,2003, e-Print: hep-ph/0210145
*
* [11] C.Bobeth, A.J.Buras, F.Kruger and J.Urban,
*     'QCD corrections to anti-B --> X/d,s nu anti-nu, anti-B/d,s --> l+ l-,
*      K--> pi nu anti-nu and K(L) --> mu+ mu- in the MSSM,'
*     Nucl.\ Phys.\  B {\bf 630} (2002) 87 [arXiv:hep-ph/0112305].
*
* [12] A.G. Akeroyd, S. Recksiegel,
*     ' The Effect of H+- on B+- ---> tau+- nu(tau) and
*       B+- ---> mu+- muon neutrino'
*     J.Phys.G29:2311-2317,2003, e-Print: hep-ph/0306037
*
*
*   - SOURCES FOR EXPERIMENTAL/LATTICE QCD DATA:
*
* [13] Y.~Amhis {\it et al.}  [Heavy Flavor Averaging Group
*      Collaboration], ``Averages of b-hadron, c-hadron, and tau-lepton
*      properties as of early 2012,'' arXiv:1207.1158 [hep-ex].
*
* [14] A. Abulencia et al. [CDF Collaboration],
*      Phys.\ Rev.\ Lett.\  {\bf 97} (2006) 242003
*      arXiv:hep-ex/0609040.
*
* [15] A. Gray et al. [HPQCD Collaboration],
*     ' The B meson decay constant from unquenched lattice QCD,'
*     Phys.\ Rev.\ Lett.\  {\bf 95} (2005) 212001 [arXiv:hep-lat/0507015].
*
* [16] E. Dalgic et al.,
*      Phys.\ Rev.\  D {\bf 76} (2007) 011501 [arXiv:hep-lat/0610104].
*
* [17] M. Okamoto,
*      'Full determination of the CKM matrix using recent results from lattice QCD,'
*      PoS {\bf LAT2005} (2006) 013 [arXiv:hep-lat/0510113].
*
* [18] P. Ball and R. Fleischer,
*      Eur.\ Phys.\ J.\  C {\bf 48} (2006) 413 [arXiv:hep-ph/0604249].
*
* [19]  The Heavy Flavor Averaging Group,
*       ``Averages of b-hadron, c-hadron, and tau-lepton Properties,''
*       arXiv:1010.1589 [hep-ex].
* [20] Upper and lower bounds on Bs -> mu+ mu- from 1211.2674 (LHCb)

      SUBROUTINE BSG(PAR,PROB)
      IMPLICIT NONE

      INTEGER I,J,K
      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
      DOUBLE PRECISION PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION TANB,COSB,SINB,au,ad,SST,SSB
      DOUBLE PRECISION ST(2),RST(2,2),SB(2),RSB(2,2),CCD(2),CCT(2,2)
      DOUBLE PRECISION FF1,FF2,FF3,FG1,FG2,FG3,ffh
      DOUBLE PRECISION fgh,esm,eh
      DOUBLE PRECISION asf,H2,C70SM,C70HIG
      DOUBLE PRECISION C80SM,C80HIG
      DOUBLE PRECISION C71HIG,C81HIG,dC7SM,dC8SM,dC7HIG
      DOUBLE PRECISION dC8HIG,C7CHARS,C8CHARS,C7CHAR,C8CHAR
      DOUBLE PRECISION XT,YT,XSQC(2),XSTC(2,2),PI,AAT,AAB,mu
      DOUBLE PRECISION xQg,xBg1,xBg2,xTg1,xTg2,xTneu(2,5),xBneu(2,5)
      DOUBLE PRECISION akk,etaS,epsb,epsbp,epst
      DOUBLE PRECISION MT0,MTH,etaH,fBssqBs,fBdsqBd,sqBBs
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION C70,C80,ALEM0,BRSL
      DOUBLE PRECISION eta,z,MB1S,etaB
      DOUBLE PRECISION delt,delt2,lndelt,lndeltp
      DOUBLE PRECISION C20b,C70b,C80b,C7EMb,C80S0,eta0
      DOUBLE PRECISION lambd2,gg1,gg2
      DOUBLE PRECISION ff11,ff12,ff17,ff18,ff22,ff27,ff28,ff77,ff78
      DOUBLE PRECISION ff88,Mch2
      DOUBLE PRECISION aa(8),bb(4),hh(8),h8(4),ee(8)
      DOUBLE PRECISION C7HIG,C8HIG,C41SM,C41HIG
      DOUBLE PRECISION dd(8),dt(8),dte(8),dim(8),da(8),db(8)
      DOUBLE PRECISION C10b,C70BSM,C80BSM
      DOUBLE PRECISION KC0,KT0,KBSM0,KC,KT,KCIM,KTIM,Ktot
      DOUBLE PRECISION KIM,KBSM,KBSMIM
      DOUBLE PRECISION sc0,VVtu,VVtuim,Vbsg,rmu
      DOUBLE PRECISION EPSew,HQET,BREMS,CCSL
      DOUBLE PRECISION af,bf,At1,A1,Ft1,F1,dAto,dA0,dFto,dF0
      DOUBLE PRECISION afim,bfim,sp2,aux
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION asmt,asmh,asmsusy,asc0,asmc,asmb
      DOUBLE PRECISION DC70BSM,DC80BSM,DKBSM
      DOUBLE PRECISION BRSG1,BRSG2
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION DSM,DH,Dchid,Dchis,DDP,ffp,ggp,gg0,S0
      DOUBLE PRECISION VVc,VVu,Vtb2,Vcs2,Vud2,VtdVtb2,VtsVtb2,runmb
      DOUBLE PRECISION MQU,MD,MS0,MC0,MB0,mmu,scR,scR2,runmass
      DOUBLE PRECISION epst0,epst1,epst2,epst3
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13
      DOUBLE PRECISION DMdexpmin,DMdexpMax,DMsexpmin,DMsexpMax,
     .       BRBMUMUexpMax,BRBMUMUexpMin,BRBTAUNUexpmin,
     .       BRBTAUNUexpMax,BRSGexpmin,BRSGexpMax
      DOUBLE PRECISION sigRLbs,sigRLbd,sigLRbs,sigLRbd,BB0,BB1
      DOUBLE PRECISION C2LR,C1SLL,C1SRR,Ca,Cs,Cp,C70S0,sgn
      DOUBLE PRECISION Vub2,fB,mBu,tauB,rh
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION DBRSGmax,DBRSGmin
      DOUBLE PRECISION dVub,dVtdVtb2,dVtsVtb2,dfB,dfBssqBs,dfBdsqBd
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .       DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .       BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllmax,CQ12,CQ22,Intpropa
      

      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
       COMMON/STSBSCALE/QSTSB
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .      DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .      BRBtaunumax,BRBtaunumin
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .      BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH

      pi=4d0*datan(1d0)

      AAT=PAR(12)
      AAB=PAR(13)
      mu=PAR(4)
      Mch2=PAR(21)
      gg1=dsqrt(g1)
      gg2=dsqrt(g2)
      
*       Alpha_s at various scales:
*   Charm Quark Mass:
       asmc=asf(mc)
*   Bottom Quark Mass:
       asmb=asf(mb)      
*   M_top:
       asmt=asf(MT)
*   Charged Higgs Mass:
       asmh=asf(CMASS)
*   Susy scale (squark masses):
       asmsusy=asf(dsqrt(QSTSB))

***********************************************************************
***********************************************************************
*      MASSES AND PARAMETERS
      
*       Matching Scale sc0 = m_top(m_top)(MSbar)
      MT0=MT/(1d0+4d0/(3d0*pi)*asmt+11d0/pi**2*asmt**2)
      sc0=MT0
      
*       Alphas at the matching scale:
      asc0=asf(sc0)

*       Trig. Functions of Beta
      TANB=PAR(3)
      sinb=tanb/dsqrt(1d0+tanb**2)
      cosb=sinb/tanb
      au=1d0/tanb
      ad=-tanb
      
*       M_top at the Charged Higgs Mass      Scale:
      mth=mt*(asmh/asmt)**(4d0/7d0)/(1d0+4d0/(3d0*pi)*asmt)
      
*       Stop Masses and Mixing Angles
      ST(1)=MST1
      ST(2)=MST2
      SST=DSQRT(1d0-CST**2)
      RST(1,1)=CST
      RST(2,2)=CST
      RST(1,2)=SST
      RST(2,1)=-SST

*       Sbottom Masses and Mixing Angles
      SB(1)=MSB1
      SB(2)=MSB2
      SSB=DSQRT(1d0-CSB**2)
      RSB(1,1)=CSB
      RSB(2,2)=CSB
      RSB(1,2)=SSB
      RSB(2,1)=-SSB

*       x Parameters (Squares of Mass Ratios)

      xt= (MT0/MW)**2              !(m_t/m_W)^2
      yt=(MTH/CMASS)**2            !(m_t/m_H+)^2
            
      do i=1,2
      xsqc(i)=(MUL/MCH(i))**2         !(m_Q/m_ch(i))^2
      do k=1,2
      xstc(k,i)=(ST(k)/MCH(i))**2       !(m_St(k)/m_ch(i))^2
      enddo
      enddo

      xQg=(MUL/MGL)**2              !(m_Q/m_gl)^2
      xBg1=(SB(1)/MGL)**2             !(m_Sb(1)/m_gl)^2
      xBg2=(SB(2)/MGL)**2             !(m_Sb(2)/m_gl)^2
      xTg1=(ST(1)/MGL)**2             !(m_St(1)/m_gl)^2
      xTg2=(ST(2)/MGL)**2             !(m_St(2)/m_gl)^2

      do i=1,5
      do k=1,2
      xTneu(k,i)=(ST(k)/MNEU(i))**2       !(m_St(k)/m_neu(i))^2
      xBneu(k,i)=(SB(k)/MNEU(i))**2       !(m_Sb(k)/m_neu(i))^2
      enddo
      enddo
      
*       .KM coefficients

      VVc=0.9879d0                 ! |V_cs.V_cb|/|V_ts.V_tb|
      VVu=0.4742d0                 ! |V_us.V_ub|/|V_ts.V_tb|
      Vtb2=(0.9991d0)**2               ! (V_tb)^2
      Vcs2=(0.97296d0)**2              ! (V_cs)^2
      Vud2=(0.97383d0)**2              ! (V_ud)^2
      
*       Since (see, e.g., [18]) V_ub(incl) differs considerably from
*       V_ub(excl), we allow for the large range
*    3.3 10^-3 < V_ub < 4.7 10^-3:
*      (V_ub)^2
      Vub2=(4.0d-3)**2
*      uncertainty on V_ub
      dVub=0.7d-3

*       From [18]: V_tb*V_td=(8.6 +/- 2.8) 10^-3  (2sigma)
*      (V_tb*V_td)^2
      VtdVtb2=(8.6d-3)**2
*      uncertainty on (V_tb*V_td)^2
      dVtdVtb2=2d0*dsqrt(VtdVtb2)*2.8d-3
*       From [18]: V_ts*V_tb=(41.3 +/- 1.4) 10^-3  (2sigma)
*      (Vtb.V_ts)^2
      VtsVtb2=(0.0413d0)**2
*      uncertainty on (Vtb.V_ts)^2
      dVtsVtb2=2d0*dsqrt(VtsVtb2)*1.4d-3

*      (V_ts.V_tb/V_cb)^2
      Vbsg=VtsVtb2/(0.042d0)**2

*      V_us.V_ub/V_ts.V_tb
      VVtu=-0.011d0
      VVtuim=0.0180d0

*      Stop/Sbottom Scale
      scR=dsqrt(QSTSB)
      scR2=QSTSB
      
       
*      Quark Masses at the Stop/Sbottom scale

      MS0=runmass(0.095d0,scR)
      MC0=runmass(1.25d0,scR)
      MB0=runmb(scR)
      MQU=runmass(0.002d0,scR)
      MD=runmass(0.005d0,scR)

*      Myon mass

      mmu=0.10566d0  !GeV
      
*      Hadronic parameters

*      From [16]:
*      f_Bs*sqrt(B_Bs)=(0.281 +/- 0.042) GeV  (2sigma)
      fBssqBs=0.281d0        ! GeV       f_Bs sqrt(B_Bs)
      dfBssqBs=0.042d0       ! GeV       uncertainty  (2sigma)

*      From [17]:
*      f_Bs*sqrt(B_Bs)/(f_Bd*sqrt(B_Bd))=1.216 +/- 0.0674 (2sigma)
*       so f_Bd*sqrt(B_Bd)=0.231*(1 +/- 0.1639) GeV  (2sigma)
      fBdsqBd=0.231d0          ! GeV       f_Bd*sqrt(B_Bd)
      dfBdsqBd=fBdsqBd*0.1639d0  ! GeV       uncertainty

*      From [15]:
*      f_B=(0.216 +/- 0.044) GeV (2sigma)
      fB=0.216d0             ! GeV       for B+ --> tau+ nu_tau
      dfB=0.044d0            ! GeV       uncertainty
      
*      Experimental data

*      Delta Md from [12], 2 sigma bounds:
*      0.499ps-1 < DMd=0.507ps-1 < 0.515ps-1
      DMdexpmin=0.499d0
      DMdexpMax=0.515d0

*      Delta Ms from [13], 2 sigma bounds:
*      17.633ps-1 < DMs=17.719ps-1 < 17.805ps-1
*     (Old, [14]: 17.53ps-1 < DMs=17.77ps-1 < 18.01ps-1)
      DMsexpmin=17.633d0
      DMsexpMax=17.805d0

*      From [20]: 2.0 10^-9 < BR(Bs->mu+mu-) < 4.7 10^-9 (95% C.L.)
      BRBMUMUexpmax=4.7d-9
      BRBMUMUexpmin=2.0d-9
      
*      BR(B+ -> tau+ nu) from [13], 2 sigma bounds:
*      1.07 10^-4 < BR(B+ -> tau+ nu)=1.67 10^-4 < 2.27 10^-4
      BRBTAUNUexpMax=2.27d-4
      BRBTAUNUexpmin=1.07d-4

*      BR(B -> Xs gamma) from [13], 2 sigma bounds:
*      3.04 10^-4 < BR(B -> Xs gamma)=3.55 10^-4 < 4.06 10^-4
      BRSGexpmin=3.04d-4
      BRSGexpMax=4.06d-4

***********************************************************************
*      EFFECTIVE NEUTRAL HIGGS COUPLINGS
*       -> EPSILON COEFFICIENTS, notation following [10]:
*
*       epsilontilde_J as in eq. (5.1) in [10] (up to a factor tanb)
*       with Delta m_d as in eq.(2.5), and Sigma as in App. (A.2)
*
*       epsilon_Y^(JI) as in eq. (5.1) in [10] (up to a factor yt^2
*       and a factor tanb), with (3.7) for lambda_0^(JI) and (3.53)
*       for the CKM matrix elements in terms of V^eff

*       a) epst1 = (epsilontilde as in[10])*tanb, epst2 (* tanb)
      epst2=asf(scR)/(3d0*pi)*(BB1(0d0,MGL**2,MDL**2,scR2)
     .      +BB1(0d0,MGL**2,MDL**2,scR2))              !gluinos

      aux=0d0
      do i=1,5
      aux=aux+2d0/H2Q*MNEU(i)                      !neutralinos
     .   *((gg1/3d0*N(i,1)-gg2*N(i,2))*N(i,3)
     .       *BB0(0d0,MNEU(i)**2,MDL**2,scR2)
     .    +DSQRT(2d0)/3d0*gg1*N(i,1)*N(i,3)
     .       *BB0(0d0,MNEU(i)**2,MDL**2,scR2))
     .       +((MS0/H2Q)**2*N(i,3)**2
     .   +1d0/2d0*(gg1/3d0*N(i,1)-gg2*N(i,2))**2)
     .       *BB1(0d0,MNEU(i)**2,MDL**2,scR2)
     .       +(2d0/9d0*gg1**2*N(i,1)**2
     .   +(MS0/H2Q)**2*N(i,3)**2)
     .       *BB1(0d0,MNEU(i)**2,MDL**2,scR2)
      enddo
      epst2=epst2+aux/(32d0*pi**2)

      aux=0d0
      do i=1,2
      aux=aux-2d0*MCH(i)*gg2/H2Q                    !charginos
     .     *U(i,2)*V(i,1)*BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .       +(MC/H1Q)**2*V(i,2)**2
     .     *BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .       +((MS0/H2Q)**2*U(i,2)**2
     .   +gg2**2*V(i,1)**2)*BB1(0d0,MCH(i)**2,MUL**2,scR2)
      enddo

      epst1=epst2+aux*Vud2/(32d0*pi**2)
      epst2=epst2+aux*Vcs2/(32d0*pi**2)
      
*       b) epst3 (* tanb)
      epst3=asf(scR)/(3d0*pi)
     . *(BB1(0d0,MGL**2,SB(1)**2,scR2)+BB1(0d0,MGL**2,SB(2)**2,scR2)
     . -2d0*MGL/MB0*SSB*CSB*(BB0(0d0,MGL**2,SB(1)**2,scR2)
     . -BB0(0d0,MGL**2,SB(2)**2,scR2)))                  !gluinos

      aux=0d0
      do i=1,5
      aux=aux                                   !neutralinos
     .   +2d0*((CSB*(gg1/3d0*N(i,1)-gg2*N(i,2))
     .     +HBQ*dsqrt(2d0)*SSB*N(i,3))
     .      *(gg1/3d0*SSB*N(i,1)+HBQ/dsqrt(2d0)*CSB*N(i,3))
     .      *BB0(0d0,MNEU(i)**2,SB(1)**2,scR2)
     .      -(SSB*(gg1/3d0*N(i,1)-gg2*N(i,2))
     .     -HBQ*dsqrt(2d0)*CSB*N(i,3))
     .      *(gg1/3d0*CSB*N(i,1)-HBQ/dsqrt(2d0)*SSB*N(i,3))
     .      *BB0(0d0,MNEU(i)**2,SB(2)**2,scR2))*MNEU(i)/MB0
     .   +((dsqrt(2d0)*gg1/3d0*SSB*N(i,1)+HBQ*CSB*N(i,3))**2
     .   +(CSB/dsqrt(2d0)*(gg1/3d0*N(i,1)-gg2*N(i,2))
     .              +HBQ*SSB*N(i,3))**2)
     .      *BB1(0d0,MNEU(i)**2,SB(1)**2,scR2)
     .   +((dsqrt(2d0)*gg1/3d0*CSB*N(i,1)-HBQ*SSB*N(i,3))**2
     .   +(SSB/dsqrt(2d0)*(gg1/3d0*N(i,1)-gg2*N(i,2))
     .              -HBQ*CSB*N(i,3))**2)
     .      *BB1(0d0,MNEU(i)**2,SB(2)**2,scR2)
      enddo
      epst3=epst3+1d0/(32d0*pi**2)*aux

      aux=0d0
      do i=1,2                                  !charginos
      aux=aux-2d0*MCH(i)/H2Q
     .    *(CST*U(i,2)*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     .       *BB0(0d0,MCH(i)**2,ST(1)**2,scR2)
     .    +SST*U(i,2)*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     .       *BB0(0d0,MCH(i)**2,ST(2)**2,scR2))
     .    +(HBQ**2*U(i,2)**2*CST**2+(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2)
     .       *BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .    +(HBQ**2*U(i,2)**2*SST**2+(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2)
     .       *BB1(0d0,MCH(i)**2,ST(2)**2,scR2)
      enddo
      epst3=epst3+Vtb2/(32d0*pi**2)*aux

*       .) epsY32 (*Yt^2 tanb), epsY31 (*Yt^2 tanb)
      epsY32=0d0
      do i=1,2
      epsY32=epsY32
     .      +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     .       *BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     .       *BB1(0d0,MCH(i)**2,ST(2)**2,scR2)
     .      -2d0*MCH(i)/H2Q*U(i,2)
     .      *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     .         *BB0(0d0,MCH(i)**2,ST(1)**2,scR2)
     .       +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     .         *BB0(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +(MS0/H2Q)**2*U(i,2)**2
     .       *(CST**2*BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +SST**2*BB1(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +VVc*((gg2*V(i,1))**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    -2d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     .            *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .    +(MC/H1Q)**2*V(i,2)**2
     .            *BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    +(MS0/H2Q)**2*U(i,2)**2
     .            *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY32=epsY32/(32d0*pi**2)
      
      epsY31=0d0
      do i=1,2
      epsY31=epsY31
     .      +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     .       *BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     .       *BB1(0d0,MCH(i)**2,ST(2)**2,scR2)
     .      -2d0*MCH(i)/H2Q*U(i,2)
     .       *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     .        *BB0(0d0,MCH(i)**2,ST(1)**2,scR2)
     .         +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     .        *BB0(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +(MD/H2Q)**2*U(i,2)**2
     .       *(CST**2*BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +SST**2*BB1(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +VVu*((gg2*V(i,1))**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    -2d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     .           *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .    +(MQU/H1Q)**2*V(i,2)**2
     .           *BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    +(MD/H2Q)**2*U(i,2)**2
     .          *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY31=epsY31/(32d0*pi**2)
      
*       d) epsY13 (*Yt^2 tanb), epsY23 (*Yt^2 tanb)
      epsY13=0d0
      do i=1,2
      epsY13=epsY13-2d0/H2Q*U(i,2)*MCH(i)
     .       *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     .    *BB0(0d0,MCH(i)**2,ST(1)**2,scR2)
     .       +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     .    *BB0(0d0,MCH(i)**2,ST(2)**2,scR2))
     .  +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     .    *BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .  +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     .    *BB1(0d0,MCH(i)**2,ST(2)**2,scR2)
     .  +HBQ/H2Q*MB0*U(i,2)**2
     .     *(CST**2*BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +SST**2*BB1(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +VVu*(gg2**2*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    +(MQU/H1Q)**2*V(i,2)**2
     .         *BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    -2d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     .         *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .    +HBQ/H2Q*MB0*U(i,2)**2
     .         *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY13=epsY13/(32d0*pi**2)

      epsY23=0d0
      do i=1,2
      epsY23=epsY23-2d0/H2Q*U(i,2)*MCH(i)
     .       *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     .    *BB0(0d0,MCH(i)**2,ST(1)**2,scR2)
     .       +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     .    *BB0(0d0,MCH(i)**2,ST(2)**2,scR2))
     .  +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     .    *BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .  +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     .    *BB1(0d0,MCH(i)**2,ST(2)**2,scR2)
     .  +HBQ/H2Q*MB0*U(i,2)**2
     .     *(CST**2*BB1(0d0,MCH(i)**2,ST(1)**2,scR2)
     .      +SST**2*BB1(0d0,MCH(i)**2,ST(2)**2,scR2))
     .      +VVc*(gg2**2*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    +(MC0/H1Q)**2*V(i,2)**2
     .         *BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    -2d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     .         *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .    +HBQ/H2Q*MB0*U(i,2)**2
     .         *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY23=epsY23/(32d0*pi**2)
      
*       e) epst0 (*tanb)
      epst0=epst3-Vtb2*epsY31


*      f) Couplings [X^s_RL]^JI as in eqs. (3.55) and (3.56), but
*       WITHOUT the S dependent mixing angles

      sigRLbs=MB0*epsY32/(H1Q*(1d0+epst0)*(1d0+epst3))
      
      sigRLbd=MB0*epsY31/(H1Q*(1d0+epst0)*(1d0+epst3))

      sigLRbs=MS0*epsY23/(H1Q*(1d0+epst0)*(1d0+epst3))
     .  *(1d0+epst3+(epst2-epst3)*epsY32/epsY23)/(1d0+epst2)

      sigLRbd=MD*epsY13/(H1Q*(1d0+epst0)*(1d0+epst3))
     .  *(1d0+epst3+(epst1-epst3)*epsY31/epsY13)/(1d0+epst1)


***********************************************************************
***********************************************************************

*     Towards DMs=m_Bs-m_Bbar_s and DMd=m_Bd-m_Bbar_d
*      (still following [10])

*     Box contributions to both DMs and DMd
*  - Standard Model
      DSM=S0(xt)/xt
      
*   - Charged Higgs
      DH=gg0(yt)/tanb**4+2d0*xt*(ffp(xt,(CMASS/MW)**2)
     .      +ggp(xt,xt,(CMASS/MW)**2)/4d0)/tanb**2

*   - Charginos
      do j=1,2
       do k=1,2
       CCT(j,k)=V(j,1)*RST(k,1)
     .       -V(j,2)*RST(k,2)*MT0/(dsqrt(2d0)*MW*sinb)
       enddo
      enddo
      
      aux=0d0                   !3rd family contribution
      do i=1,2
       do j=1,2
        aux=aux+1d0/MCH(j)**2
     .     *(CCT(j,1)**2*CCT(i,1)**2
     .       *ggp(xstc(1,j),xstc(1,j),(MCH(i)/MCH(j))**2)
     .      +CCT(j,2)**2*CCT(i,2)**2
     .  *ggp(xstc(2,j),xstc(2,j),(MCH(i)/MCH(j))**2))
        if(i.ne.j)aux=aux+1d0/MCH(j)**2
     .      *CCT(j,1)*CCT(i,1)*CCT(i,1)*CCT(j,2)
     .       *(ggp(xstc(1,j),xstc(2,j),(MCH(i)/MCH(j))**2)
     .  +ggp(xstc(2,j),xstc(1,j),(MCH(i)/MCH(j))**2))
       enddo
      enddo

      Dchid=MW**2/xt*aux
      Dchis=MW**2/xt*aux

      aux=0d0                   !3rd family/1st-2nd interference
      do i=1,2
       do j=1,2
        if(i.ne.j)aux=aux+1d0/MCH(j)**2
     .      *(V(j,1)*V(i,1)*CCT(i,1)*CCT(j,1)
     .      *ggp(xsqc(j),xstc(1,j),(MCH(i)/MCH(j))**2)
     .       +V(j,1)*V(i,1)*CCT(i,2)*CCT(j,2)
     .      *ggp(xsqc(j),xstc(2,j),(MCH(i)/MCH(j))**2)
     .       +V(i,1)*V(j,1)*CCT(j,1)*CCT(i,1)
     .      *ggp(xstc(1,j),xsqc(j),(MCH(i)/MCH(j))**2)
     .       +V(i,1)*V(j,1)*CCT(j,2)*CCT(i,2)
     .      *ggp(xstc(2,j),xsqc(j),(MCH(i)/MCH(j))**2))
       enddo
      enddo

      Dchid=Dchid+VVu*MW**2/xt*aux
      Dchis=Dchis+VVc*MW**2/xt*aux

      aux=0d0                   !1st-2nd family contribution
      do i=1,2
       do j=1,2
        aux=aux+1d0/MCH(j)**2*V(j,1)**2*V(j,2)**2
     .      *ggp(xsqc(j),xsqc(j),(MCH(i)/MCH(j))**2)
       enddo
      enddo

      Dchid=Dchid+VVu**2*MW**2/xt*aux
      Dchis=Dchis+VVc**2*MW**2/xt*aux

***********************************************************************
      
*      Double Penguin Contributions to DMd
*   - Wilson coefficients
      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-(5.2793d0)**2)/
     . dsqrt((SMASS(i)**2-(5.2793d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-(5.2793d0)**2)/
     . dsqrt((PMASS(i)**2-(5.2793d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C2LR=-(4d0*pi/(GF*MW))**2*sigRLbd*sigLRbd*aux

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-(5.2793d0)**2)/
     . dsqrt((SMASS(i)**2-(5.2793d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      epst=aux
      do i=1,2
      aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-(5.2793d0)**2)/
     . dsqrt((PMASS(i)**2-(5.2793d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C1SLL=-(4d0*pi/(GF*MW))**2*sigRLbd**2*aux/2d0

      C1SRR=-(4d0*pi/(GF*MW))**2*sigLRbd**2*aux/2d0

      DDP=(0.90d0*C2LR-0.37d0*(C1SLL+C1SRR))/xt

**********************************************************************      

*      Results for DMd
      etaB=0.551d0

      aux=(GF**2*MW**2/(6d0*pi**2)*etaB*5.2794d0*(fBdsqBd)**2
     .       *xt*VtdVtb2/(6.58211915d-13))

      DMd=aux*abs(DSM+DH+Dchid+DDP)


*      Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
*       * lattice QCD (sources [16,17]):
*       1.134 < fBd sqrt(BBd) / fBs sqrt(BBs)=1.216 < 1.298
*       0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
*      * CKM factor (source [18]):
*       6.2 10^-3 < VtbVtd=8.6 10^-3 < 11. 10^-3
*       * 30% on BSM (1st order QCD) contributions

*     First: 2 sigma "SM" (relative) error bars from CKM and
*       lattice uncertainties: (2sigma, added quadratiCALLy)
      DMdMax=(1d0+dsqrt((dVtdVtb2/VtdVtb2)**2
     .  +(2d0*dfBdsqBd/fBdsqBd)**2))
      DMdmin=(1d0-dsqrt((dVtdVtb2/VtdVtb2)**2
     .  +(2d0*dfBdsqBd/fBdsqBd)**2))

*      Total error bars, allowing for 30% theory error
*  on each BSM contribution:
      DMdmax=DMdMax*(DMd+aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))
      
      DMdmin=DMdmin*(DMd-aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))

*      Comparison with experimental data (source [12]):
*      (Recall: 2 sigma bounds: 0.499ps-1 < DMd < 0.515ps-1)

      prob(34)=0d0

      IF(DMdmin.GE.DMdexpMax)
     .     PROB(34)=DMdmin/DMdexpMax-1d0
      IF(DMdmax.LE.DMdexpmin)
     .     PROB(34)=DMdmax/DMdexpMin-1d0

**********************************************************************      
**********************************************************************      

*       Double Penguin Contributions to DMs
*   - Wilson coefficients
      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-(5.3696d0)**2)/
     . dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-(5.3696d0)**2)/
     . dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C2LR=-(4d0*pi/(GF*MW))**2*sigRLbs*sigLRbs*aux

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-(5.3696d0)**2)/
     . dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-(5.3696d0)**2)/
     . dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C1SLL=-(4d0*pi/(GF*MW))**2*sigRLbs**2*aux/2d0

      C1SRR=-(4d0*pi/(GF*MW))**2*sigLRbs**2*aux/2d0      
      
      DDP=(0.90d0*C2LR-0.37d0*(C1SLL+C1SRR))/xt

**********************************************************************      

*      Results for DMs      

      aux=GF**2*MW**2/(6d0*pi**2)*0.55d0*5.3696d0*(fBssqBs)**2
     .       *VtsVtb2*xt/(6.58211915d-13)
      
      DMs=aux*dabs(DSM+DH+Dchis+DDP)

*      Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
*      * lattice QCD (source [16]):
*       0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
*      * CKM factor (source [18]):
*       37.9 10^-3 < VtsVtb=41.3 10^-3 < 44.7 10^-3
*      * 30% on BSM (1st order QCD) contributions

*      First: 2 sigma "SM" (relative) error bars from CKM and
*       lattice uncertainties: (2sigma, added quadratiCALLy)
      DMsMax=(1d0+dsqrt((dVtsVtb2/VtsVtb2)**2
     .  +(2d0*dfBssqBs/fBssqBs)**2))
      DMsmin=(1d0-dsqrt((dVtsVtb2/VtsVtb2)**2
     .  +(2d0*dfBssqBs/fBssqBs)**2))
      
*      Total error bars, allowing for 30% theory error
*  on each BSM contribution:

      DMsmax=DMsMax*(DMs+aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))
      DMsmin=DMsmin*(DMs-aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))

*      Comparison with experimental data (source [14])
*       (Recall: 2 sigma bounds: 17.633ps-1 < DMs < 17.805ps-1)

      prob(33)=0d0

      IF(DMsmin.GE.DMsexpMax)
     .     PROB(33)=DMsmax/DMsexpMax-1d0
      IF(DMsmax.LE.DMsexpmin)
     .     PROB(33)=DMsmax/DMsexpMin-1d0

**********************************************************************      
**********************************************************************      

*      Contributions to BR(Bs -> mu+ mu-)    (see [10,11])
*       a) Wilson coefficients - SM contribution
      ca=1d0/4d0*(xt/(1d0-xt)+xt/(xt-1d0)**2*dlog(xt))
     .   -xt/8d0*((xt-6d0)/(xt-1d0)
     .   +(3d0*xt+2d0)/(xt-1d0)**2*dlog(xt))

*       b) Wilson coefficients - DP contribution
      aux=0d0
      do i=1,3
      aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)
     . *sgn(SMASS(i)**2-(5.3696d0)**2)/(cosb
     . *dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2))
      enddo
      cs=gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

      aux=0d0
      do i=1,2
      aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     .  *PCOMP(i,1)*tanb*sgn(PMASS(i)**2-(5.3696d0)**2)/
     . dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo
      cp=-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

*       .) Branching ratio      
      sqBBs=dsqrt(1.3d0)       !([10])       sqrt(B_Bs)

      aux=GF**2*ALEM0**2/(s2tw**2*64d0*pi**3)*(5.3696d0)**5
     .     *1.454d0/(6.58211915d-13)*(fBssqBs/sqBBs)**2
     .     *Vtsvtb2*dsqrt(1d0-4d0*mmu**2/5.3696d0**2)

      BRBMUMU=aux*
     .      ((1d0-4d0*mmu**2/(5.3696d0)**2)/(1d0+MS/Mb)**2*Cs**2
     .       +(Cp/(1d0+ms/mb)+2d0*mmu/(5.3696d0)**2*Ca)**2)

*      Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
*      * lattice QCD (source [16]):
*       0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
*      * CKM factor (source [18]):
*       37.9 10^-3 < VtsVtb=41.3 10^-3 < 44.7 10^-3
*      * 30% on BSM (1st order QCD) contributions

*       First: 2 sigma "SM" (relative) error bars from CKM and
*       lattice uncertainties: (2sigma, added quadratiCALLy)
      BRBMUMUMAX=(1d0+dsqrt((dVtsVtb2/VtsVtb2)**2
     .      +(2d0*dfBssqBs/fBssqBs)**2))
      BRBMUMUmin=(1d0-dsqrt((dVtsVtb2/VtsVtb2)**2
     .      +(2d0*dfBssqBs/fBssqBs)**2))
      
*      Total error bars, allowing for 30% theory error
*  on each BSM contribution:
      BRBMUMUmax=BRBMUMUMAX*(BRBMUMU+aux*0.3d0*
     .      ((1d0-4d0*mmu**2/(5.3696d0)**2)/(1d0+MS/Mb)**2*Cs**2
     .       +(Cp/(1d0+ms/mb)+2d0*mmu/(5.3696d0)**2*Ca)**2))
      BRBMUMUmin=BRBMUMUmin*(BRBMUMU-aux*0.3d0*
     .      ((1d0-4d0*mmu**2/(5.3696d0)**2)/(1d0+MS/Mb)**2*Cs**2
     .       +(Cp/(1d0+ms/mb)+2d0*mmu/(5.3696d0)**2*Ca)**2))

*      Comparison with experimental data (source [20])
*    (Recall 95% C.L.: 2.0 10^-9 < BR(Bs->mu+mu-) < 4.7 10^-9)

      prob(35)=0d0

      IF(BRBMUMUmin.GE.BRBMUMUexpMax)
     .     PROB(35)=BRBMUMUmin/BRBMUMUexpMax-1d0
      IF(BRBMUMUmax.LE.BRBMUMUexpMin)
     .     PROB(35)=BRBMUMUmax/BRBMUMUexpMin-1d0

**********************************************************************      
**********************************************************************      

*      Branching ratio BR(B+ -> tau+ nu_tau)
*                           following [12]
      tauB=1.638d-12/(6.58211915d-25)
      mBu=5.279d0

      rh=(1d0-(mBu/CMASS)**2*tanb**2/(1d0+epst0))**2
      aux=GF**2*mBu*MTAU**2/(8d0*pi)*(1d0-(MTAU/mBu)**2)**2
     .       *fB**2*Vub2*tauB
      BRBtaunu=aux*rh

*      Comparison with experimental data:
*   2 sigma bounds (source [19]):
*   0.89 10^-4 < BR(B+ -> tau+ nu) < 2.45 10^-4
*  hadronic parameter (source [14]; 2sigma):
*   0.172 GeV < fB=0.216 GeV <0.260 GeV
*  CKM factor (2sigma): 3.3 10^-3 < Vub=4.0 10^-3 < 4.7 10^-3

*  2 sigma (absolute) error bars from CKM and
*       lattice uncertainties:

      BRBtaunumax=(1d0+dfB/fB)**2*(1d0+dVub/dsqrt(Vub2))**2*BRBtaunu
      BRBtaunumin=(1d0-dfB/fB)**2*(1d0-dVub/dsqrt(Vub2))**2*BRBtaunu

      prob(36)=0d0

      IF(BRBtaunumin.GE.BRBTAUNUexpmax)
     .     PROB(36)=BRBtaunumin/BRBTAUNUexpmax-1d0
      IF(BRBtaunumax.LE.BRBTAUNUexpmin)
     .     PROB(36)=BRBtaunumax/BRBTAUNUexpmin-1d0

**********************************************************************      
**********************************************************************      
*  Contributions to BR(B->Xs gamma) at the Matching scale (sc0)
***********************************************************************

*   Towards SUSY Contributions including large tan(beta) effects from
*   charged Higgs couplings      (Used: [2], [3], [6])

*   First:  Effective Charged Higgs Yukawa Couplings
*    Epsilon_b, Epsilon'_b and Epsilon_t

*      a) epsilon_b
       epsb=-asmsusy*2d0/(3d0*pi)*
     .       ((mu-AAB/tanb)/MGL*H2(xBg1,xBg2)
     .    +1d0/(2d0*tanb)*(1d0-BB1(0d0,MGL**2,SB(1)**2,QSTSB)
     .    -BB1(0d0,MGL**2,SB(2)**2,QSTSB)))

       aux=0d0
       do j=1,2
       aux=aux+U(j,2)*H2(xstc(1,j),xstc(2,j))*V(j,2)/MCH(j)
       enddo
      
       epsb=epsb-(HTQ/(4d0*pi))**2*(AAT-mu/tanb)*aux
      
       aux=(CST/ST(1))**2*H2(Mch2**2/ST(1)**2,mu**2/ST(1)**2)
     .       +(SST/ST(2))**2*H2(Mch2**2/ST(2)**2,mu**2/ST(2)**2)
     .       +(CSB/SB(1))**2*H2(Mch2**2/SB(1)**2,mu**2/SB(1)**2)/2d0
     .       +(SSB/SB(2))**2*H2(Mch2**2/SB(2)**2,mu**2/SB(2)**2)/2d0
      
       epsb=epsb+ALEMMZ/(4d0*S2TW*pi)*mu*Mch2*aux      

*      b) epsilon'_b(t)
       epsbp=-asmsusy*2d0/(3d0*pi)*(mu-AAB/tanb)/MGL*
     .     (CST**2*(H2(xTg1,xBg2)*CSB**2+H2(xTg1,xBg1)*SSB**2)+
     .      SST**2*(H2(xTg2,xBg2)*CSB**2+H2(xTg2,xBg1)*SSB**2))

       aux=0d0
       do j=1,5
       aux=aux+1d0/MNEU(j)*N(j,4)*N(j,3)*
     .     (CST**2*(H2(xTneu(2,j),xBneu(1,j))*CSB**2+
     .      H2(xTneu(2,j),xBneu(2,j))*SSB**2)+
     .      SST**2*(H2(xTneu(1,j),xBneu(1,j))*CSB**2+
     .      H2(xTneu(1,j),xBneu(2,j))*SSB**2))
        enddo

      epsbp=epsbp+(HTQ/(4d0*pi))**2*(AAT-mu/tanb)*aux
      
      aux=(CST/ST(1))**2*H2(Mch2**2/ST(1)**2,mu**2/ST(1)**2)/2d0
     .       +(SST/ST(2))**2*H2(Mch2**2/ST(2)**2,mu**2/ST(2)**2)/2d0
     .       +(CSB/SB(1))**2*H2(Mch2**2/SB(1)**2,mu**2/SB(1)**2)
     .       +(SSB/SB(2))**2*H2(Mch2**2/SB(2)**2,mu**2/SB(2)**2)
      
      epsbp=epsbp+ALEMMZ/(4d0*S2TW*pi)*mu*Mch2*aux

*      C) epsilon_t(s)
       epst= -asmsusy*2d0/(3d0*pi)*(mu+AAT/tanb)/MGL*
     .     (CST**2*H2(xTg2,xQg)+SST**2*H2(xTg1,xQg))

            aux=0d0
       do,i=1,5
       aux=aux+N(i,4)*N(i,3)/MNEU(i)*(
     .    CST**2*CSB**2*H2(ST(1)**2/MNEU(i)**2,SB(2)**2/MNEU(i)**2)
     .   +CST**2*SSB**2*H2(ST(1)**2/MNEU(i)**2,SB(1)**2/MNEU(i)**2)
     .   +SST**2*CSB**2*H2(ST(2)**2/MNEU(i)**2,SB(2)**2/MNEU(i)**2)
     .   +SST**2*SSB**2*H2(ST(2)**2/MNEU(i)**2,SB(1)**2/MNEU(i)**2))

       enddo

        epst=epst+HBQ**2/(16d0*pi**2)*mu/tanb*aux
      
***********************************************************************
*      Chargino/Squark Contributions

*      1) Factors
       do j=1,2
        CCD(J)=U(J,2)*MW/(dsqrt(2d0)*COSB*MCH(J))
       do k=1,2
        CCT(j,k)=V(j,1)*RST(k,1)
     .    -V(j,2)*RST(k,2)*HTQ/dsqrt(g2)
       enddo
       enddo

*      2) Calculation of the Wilson Coefficients (SUSY Scale)
       C7CHARS=0d0
       C8CHARS=0d0
      
       akk=1d0/(1d0+epsb*tanb)
      
       do j=1,2
        C7CHARS=C7CHARS+2d0/3d0*V(J,1)**2*(MW/MUL)**2*FF1(XSQC(J))
     .       +akk*CCD(J)*V(J,1)*FF3(XSQC(J))

        C8CHARS=C8CHARS+2d0/3d0*V(J,1)**2*(MW/MUL)**2*FG1(XSQC(J))
     .       +akk*CCD(J)*V(J,1)*FG3(XSQC(J))

       do k=1,2  ! k is the stop index
         C7CHARS=C7CHARS-
     .       2d0/3d0*CCT(J,K)**2*(MW/ST(K))**2*FF1(XSTC(K,J))
     .       -akk*CCD(J)*CCT(J,K)*RST(K,1)*FF3(XSTC(K,J))

         C8CHARS=C8CHARS-
     .       2d0/3d0*CCT(J,K)**2*(MW/ST(K))**2*FG1(XSTC(K,J))
     .       -akk*CCD(J)*CCT(J,K)*RST(K,1)*FG3(XSTC(K,J))

       enddo
       enddo
      
*      3) Evolution from the SUSY scale to sc0

       etaS= asmsusy/asc0

       C7CHAR=etaS**(16d0/21d0)*C7CHARS+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARS
       C8CHAR=etaS**(14d0/21d0)*C8CHARS
      
***********************************************************************
*       Charged Higgs Contributions
*      1) Lowest Order:

      C70HIG=au**2/3d0*ff1(yt)-au*ad*ff2(yt)
      C80HIG=au**2/3d0*fg1(yt)-au*ad*fg2(yt)

*      2) Order alpha_s:

       C41HIG=au**2*eh(yt)
       C71HIG=ffh(yt,tanb)-4d0/9d0*C41HIG
       C81HIG=fgh(yt,tanb)-1d0/6d0*C41HIG
      
*      3) Large tan(beta) Corrections:

       dC7HIG=-tanb*(epsb+epst)*akk*ff2(yt)
       dC8HIG=-tanb*(epsb+epst)*akk*fg2(yt)
      
*      4) Evolution from M_Higgs to sc0:

       etaH=asmh/asc0      

       C7HIG=etaH**(16d0/21d0)*(C70HIG+dC7HIG)
     .       +8d0/3d0*(etaH**(2d0/3d0)-etaH**(16d0/21d0))
     .       *(C80HIG+dC8HIG)
       C8HIG=etaH**(2d0/3d0)*(C80HIG+dC8HIG)

***********************************************************************      
*       SM Contributions at the matching scale sc0 = m_top

*       1) Lowest Order:
       C70SM=ff1(xt)
       C80SM=fg1(xt)

*       2) Order alpha_s correction to C41:

       C41SM=esm(xt)-2d0/3d0+2d0/3d0*dlog(sc0**2/MW**2)

*      3) Large tan(beta) Corrections:

       dC7SM=tanb*(epsb-epsbp)*akk*ff2(xt)
       dC8SM=tanb*(epsb-epsbp)*akk*fg2(xt)

***********************************************************************      
*      Neutral Higgs Contribution to b -> s gamma
*      (following [10], but including the evolution from the scale
*      PMASS(1) to the matching scale sc0)

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)
     .   *(SCOMP(i,2)+SCOMP(i,1)*epst3/tanb)/SMASS(i)**2
      enddo
      do i=1,2
      aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)
     .       *(PCOMP(i,1)*sinb-PCOMP(i,1)*cosb*epst3/tanb)/PMASS(i)**2
      enddo
      
      C70S0=1d0/18d0*MW**2/g2*MB0/(H2Q*(1d0+epst3))*sigRLbs*aux      

*       - Running to mt
      IF(PMASS(1).le.sc0)THEN
       IF(PMASS(1).ge.mb)THEN
        eta0=asc0/asf(PMASS(1))
       ELSE
        eta0=asmb/asc0
       ENDIF
      ELSE
        eta0=1d0
      ENDIF
      C80S0=eta0**(14d0/23d0)*C70S0
      C70S0=eta0**(16d0/23d0)*C70S0+
     .     8d0/3d0*(eta0**(14d0/23d0)-eta0**(16d0/23d0))*C70S0

***********************************************************************      
*       Sums at the matching scale sc0 = m_top
      
       C70BSM=dC7SM+C7HIG+C7CHAR+C70S0
       C80BSM=dC8SM+C8HIG+C8CHAR+C80S0
                 
*  DC70BSM and DC80BSM are (conservative) error estimates:
*  Dominant errors: HIG: Order (alphas**2) ~ 10%,
*  CHAR: Order(alphas) ~ 30%
                 
       DC70BSM=.1d0*DABS(C7HIG)+.3d0*DABS(C7CHAR)+.3d0*dabs(C70S0)
       DC80BSM=.1d0*DABS(C8HIG)+.3d0*DABS(C8CHAR)+.3d0*dabs(C80S0)

       C70=C70SM+C70BSM
       C80=C80SM+C80BSM
      
***********************************************************************
***********************************************************************
*      Coefficients for BR(B->Xs gamma) at the Scale m_b
*       Used: [1], [4], [5], [7], [8]
***********************************************************************
*
*      Parameters:

      eta=asc0/asmb

*      A: LO COEFFICIENTS C10b, C20b, C70b, C80b
*  (Used in Bremsstrahlung and HQET Corrections)

*    1) "Magic Numbers"
*       a) Array: aa_i
      aa(1)=14d0/23d0
      aa(2)=16d0/23d0
      aa(3)=6d0/23d0
      aa(4)=-12d0/23d0
      aa(5)=0.4086d0
      aa(6)=-0.4230d0
      aa(7)=-0.8994d0
      aa(8)=0.1456d0

*       b) Array: bb_i
      do i=1,4
       j=4+i
       bb(i)=aa(j)
      enddo

*       .) Array: hh_i
      hh(1)=626126d0/272277d0
      hh(2)=-56281d0/51730d0
      hh(3)=-3d0/7d0
      hh(4)=-1d0/14d0
      hh(5)=-0.6494d0
      hh(6)=-0.0380d0
      hh(7)=-0.0185d0
      hh(8)=-0.0057d0

*       d) Array: h8_i
      h8(1)=-0.9135d0
      h8(2)=0.0873d0
      h8(3)=-0.0571d0
      h8(4)=0.0209d0
      
*    2) Wilson Coefficients

      C10b=eta**(aa(3))-eta**(aa(4))

      C20b=1.d0/3.d0*(2.d0*eta**(aa(3))+eta**(aa(4)))

      C70b=0.d0
      do i=1,8
       C70b=C70b+hh(i)*eta**(aa(i))
      enddo
      C70b=C70b+eta**(16d0/23d0)*C70+
     .   8d0/3d0*(eta**(14d0/23d0)-eta**(16d0/23d0))*C80

      C80b=0d0
      do i=1,4
       C80b=C80b+h8(i)*eta**(bb(i))
      enddo
      C80b=C80b+eta**(14d0/23d0)*(C80+313063d0/363036d0)

***********************************************************************
*      B: Lowest Order c-Quark (KC0) and t-Quark contributions (KT0)
*      (Only for Bremsstrahlung and HQET Corrections)

      KC0=0d0
      
      do i=1,8
       KC0=KC0+hh(i)*eta**(aa(i))
      enddo

      KC0=KC0-23d0/36d0*eta**(16d0/23d0)
     .     -8d0/9d0*(eta**(14d0/23d0)-eta**(16d0/23d0))

           KT0=eta**(4d0/23d0)*(C70SM+23d0/36d0)
     .   +8d0/3d0*(eta**(2d0/23d0)-eta**(4d0/23d0))
     .           *(C80SM+1d0/3d0)
           
      KBSM0=C70BSM*ETA**(4d0/23d0)+
     .    C80BSM*(ETA**(2d0/23d0)-ETA**(4d0/23d0))*8d0/3d0
            
***********************************************************************
*      C: Towards KC and KT including NLO in alpha_s

*    1) "Magic Numbers" needed for NLO
*       a) Array: e_i
      ee(1)=5.2620d0
      ee(2)=-8516d0/2217d0
      ee(3)=0d0
      ee(4)=0d0
      ee(5)=-1.9043d0
      ee(6)=-0.1008d0
      ee(7)=0.1216d0
      ee(8)=0.0183d0
      
*       b) Array: dd_i
      do i=1,8
       dd(i)=hh(i)
      enddo
      dd(1)=dd(1)-8d0/9d0
      dd(2)=dd(2)+8d0/9d0-23d0/36d0
      
*       c) Array: dt_i
      dt(1)=-17.6507d0
      dt(2)=11.3460d0
      dt(3)=2.4692d0
      dt(4)=-0.8056d0
      dt(5)=4.8898d0
      dt(6)=-0.2308d0
      dt(7)=-0.5290d0
      dt(8)=0.1994d0

*       d) Array: dte_i
      dte(1)=9.2746d0
      dte(2)=-6.9366d0
      dte(3)=-0.8740d0
      dte(4)=0.4218d0
      dte(5)=-2.7231d0
      dte(6)=0.4083d0
      dte(7)=0.1465d0
      dte(8)=0.0205d0
      
*       e) Array: dim_i
      dim(1)=0.4702d0
      dim(2)=0d0
      dim(3)=-0.4268d0
      dim(4)=-0.2222d0
      dim(5)=-0.9042d0
      dim(6)=0.1150d0
      dim(7)=-0.0975d0
      dim(8)=0.0115d0
      
*       f) Array: da_i
      da(1)=0d0
      da(2)=0d0
      da(3)=0.8571d0
      da(4)=0.6667d0
      da(5)=0.1298d0
      da(6)=0.1951d0
      da(7)=0.1236d0
      da(8)=0.0276d0
      
*       g) Array: db_i
      db(1)=0d0
      db(2)=0d0
      db(3)=0.8571d0
      db(4)=0.6667d0
      db(5)=0.2637d0
      db(6)=0.2906d0
      db(7)=-0.0611d0
      db(8)=-0.0171d0
      
*    2) Charm Quark Contribution KC, LO + NLO (SM only)

*    For z=MC/MB, we choose a value that reproduces the NNLO result
*    BRSG ~ 3.15 10^(-4) [8] in the SM:
      z=(0.307d0)**2

      KC=0d0
      do i=1,8
      KC=KC+eta**(aa(i))*(dd(i)*(1d0+46d0/3d0*aa(i)
     .    *asmb/(4d0*pi)*eta*dlog(sc0/MW))
     .    +asmb/(4d0*pi)*(dt(i)+dte(i)*eta
     .    +(1d0+VVtu)*(da(i)*af(z)+db(i)*bf(z))
     .    -VVtuim*(da(i)*afim(z)+db(i)*bfim(z))))
      enddo

*    3) Top Quark Contribution KT, LO + NLO (SM only)

      At1=A1(xt)
      Ft1=F1(xt)
      dAto=dA0(xt)
      dFto=dF0(xt)
      
      aux=0d0
      do i=1,8
       aux=aux+ee(i)*eta**(aa(i)+11d0/23d0)
      enddo
      
      KT=(1d0-2d0/9d0*asmb**2)*
     .      (eta**(4d0/23d0)*(C70SM+23d0/36d0)
     .       +8d0/3d0*(eta**(2d0/23d0)-eta**(4d0/23d0))
     .                 *(C80SM+1d0/3d0))
      KT=KT+asmb/(4d0*pi)*(aux*(C41SM+7d0/9d0)
     .      +eta**(4d0/23d0)*(eta*(-1d0/2d0*At1)
     .      -2d0*(12523d0/3174d0-7411d0/4761d0*eta
     .   -2d0/9d0*pi**2)
     .      *(C70SM+23d0/36d0)-8d0/3d0*eta*(-1d0/2d0*Ft1)
     .      -2d0*(-50092d0/4761d0+1110842d0/357075d0*eta
     .  +16d0/27d0*pi**2)*(C80SM+1d0/3d0))
     .      +eta**(2d0/23d0)*(8d0/3d0*eta*(-1d0/2d0*Ft1)
     .      -2d0*(2745458d0/357075d0-38890d0/14283d0*eta
     . -4d0/9d0*pi*(pi))*(C80SM+1d0/3d0)))
      KT=KT+asc0/pi*dlog(sc0/MT)*4*xt*
     .  (-1d0/2d0*eta**(4d0/23d0)*dAto
     .   +4d0/3d0*(eta**(4d0/23d0)-eta**(2d0/23d0))*dFto)
     .       +asmb/(4d0*pi)*((8d0/3d0*(C70SM+23d0/36d0)
     .     -64d0/9d0*(C80SM+1d0/3d0))*eta**(27d0/23d0)
     .     +32d0/9d0*eta**(25d0/23d0)*(C80SM+1d0/3d0))
     .                   *dlog(sc0/MT)

*       4) BSM Contribution to KT, LO + NLO

      KBSM=(1d0-2d0/9d0*asmb**2)*(eta**(4d0/23d0)*C70BSM
     .   +8d0/3d0*(eta**(2d0/23d0)-eta**(4d0/23d0))*C80BSM)

      KBSM=KBSM+asmb/(4d0*pi)*(aux*C41HIG
     .      +eta**(4d0/23d0)*(eta*C71HIG
     .      -2d0*(12523d0/3174d0-7411d0/4761d0*eta
     .     -2d0/9d0*pi**2)
     .      *C70BSM-8d0/3d0*eta*C81HIG
     .      -2d0*(-50092d0/4761d0+1110842d0/357075d0*eta
     .    +16d0/27d0*pi**2)*C80BSM)
     .      +eta**(2d0/23d0)*(8d0/3d0*eta*C81HIG
     .      -2d0*(2745458d0/357075d0-38890d0/14283d0*eta
     .    -4d0/9d0*pi*(pi))*C80BSM))

*   BSM error estimate from DC70BSM and DC80BSM:
      DKBSM=(1d0-2d0/9d0*asmb**2)*(eta**(4d0/23d0)*DC70BSM
     .   +8d0/3d0*(eta**(2d0/23d0)-eta**(4d0/23d0))*DC80BSM)

*       5) Imaginary Parts of KC, KT, KBSM

      KCIM=0d0
      do i=1,8
       KCIM=KCIM+eta**(aa(i))*(dim(i)*pi
     .   +(1d0+VVtu)*(da(i)*afim(z)+db(i)*bfim(z))
     .   +VVtuim*(da(i)*af(z)+db(i)*bf(z)))
      enddo
      KCIM=asmb/(4d0*pi)*KCIM

      KTIM=2d0*asmb/9d0*eta**(2d0/23d0)*(C80SM+1d0/3d0)

      KBSMIM=2d0*asmb/9d0*eta**(2d0/23d0)*C80BSM

**********************************************************************
***********************************************************************
*    Towards BR(b -> s gamma) including Electroweal corrections,
*      Bremsstrahlung Corrections and Heavy Quark Effective Theory
*      Corrections
*                           Used: [4], [7]
***********************************************************************

*   MB1S = "1S b-Quark Mass":
      MB1S=4.68d0
      
*       Ratio rmu = m_b(sc0)/MB1S:

      rmu=0.578d0*0.1185d0/alsmz*(MB1S/4.69d0)**0.23d0
     . *(MC*(1d0+asmc/pi*(-4d0/3d0))/1.25d0
     .  )**(-0.003d0)*(sc0/165d0)**(-0.08d0)*(MB/4.69d0)**(0.006d0)

***********************************************************************
*   A: Summing the LO and NLO Contributions

      Ktot=KC+rmu*(KT+KBSM)
      KIM=KCIM+rmu*(KTIM+KBSMIM)

***********************************************************************
*  B: ELECTROWEAK Corrections from
*   P.Gambino, U.Haisch,
*   'Complete Electroweak Matching for Radiative B Decays'
*   JHEP 0110:020,2001, e-Print: hep-ph/0109058

      C7EMb=(88d0/575d0*eta**(16d0/23d0)
     .      -40d0/69d0*eta**(-7d0/23d0)
     .      +32d0/75d0*eta**(-9d0/23d0))
     .   *(C70BSM+asc0/(4d0*pi)*C71HIG)
     .     +(640d0/1449d0*eta**(14d0/23d0)
     .     -704d0/1725d0*eta**(16d0/23d0)
     .     +32d0/1449d0*eta**(-7d0/23d0)
     .     -32d0/575d0*eta**(-9d0/23d0))
     .   *(C80BSM+asc0/(4d0*pi)*C81HIG)
*      Summation of the Electroweak Corrections, the SM Contribution
*  is put to 0.0071

*      EPSew=0.0071d0
*     .     +ALEMMZ/(4d0*pi)*(C7EMb-4d0*rmu*KBSM0*dlog(MZ/MB))

      EPSew=0.0071d0
     .     +ALEMMZ*(C7EMb/asmb-rmu*KBSM0*dlog(MZ/MB)/pi)


***********************************************************************
*    C: Bremsstrahlung Contribution with delta=0.352
* For the cut on the photonic energy, we take E_ph > 1.6 GeV, that is
* delt=0.352. As for the f_ij coefficients, we will use numerical
* values for ff_22 and ff_27.
* Formulae for delt=0.9 will be kept as comments.

      delt=0.352d0

*       a) Coefficients: f_ij
      lndelt=dlog(delt)
      lndeltp=dlog(1d0-delt)
      delt2=delt**2

*   If delt=0.9 we would have
*      ff22=0.107636d0-0.208484d0*dsqrt(z)-0.156146d0*z
*   Now: delt=0.325:
      ff22=-0.01445908619d0+.9889254405d0*dsqrt(z)-6.618735072d0*z
     .       +17.52172828d0*z*dsqrt(z)-6.886559723d0*z**2
     .       -60.56748284d0*z**2*dsqrt(z)+89.45278858d0*z**3
      ff22=ff22*((1d0+VVtu)**2+VVtuim**2)+(VVtu**2+VVtuim**2)
     .   *(1d0/3d0*(1d0-(1d0-delt)**3)+1d0/2d0*(1d0-delt)**2)
      ff11=1d0/36d0*ff22
      ff12=-1d0/3d0*ff22
      
*   If delt=0.9 we would have
*      ff27=-0.190805d0+0.948865d0*dsqrt(z)-0.787805d0*z
*   Now: delt=0.325:
      ff27=-0.04047347658d0-.4349867931d0*dsqrt(z)+2.997393353d0*z
     .       -2.340077892d0*z*dsqrt(z)-4.025847473d0*z**2
     .       +7.412298665d0*z**2*dsqrt(z)-13.15396549d0*z**3
      ff27=(1d0+VVtu)*ff27
      ff17=-1d0/6d0*ff27
      ff28=-1d0/3d0*ff27
      ff18=-1d0/6d0*ff28
      ff77=1d0/3d0*(10d0*delt+delt2-2d0/3d0*delt**3
     .     +delt*(delt-4d0)*lndelt)
     .   -1d0/3d0*(2d0*lndelt**2+7d0*lndelt+31d0/3d0)
      ff78=8d0/9d0*(sp2(1d0-delt)-pi**2/6d0-delt*lndelt
     .      +9d0/4d0*delt-delt2/4d0+delt**3/12d0)
      ff88=1d0/27d0*(4d0*sp2(1d0-delt)-2d0/3d0*pi**2
     .       +8d0*lndeltp-delt*(2d0+delt)*lndelt
     .       +7d0*delt+3d0*delt2-2d0/3d0*delt**3
     .       -2d0*(2d0*delt+delt2+4d0*lndeltp)*dlog(MB/0.95d0))

*       b) Summing the Bremsstrahlung Contributions
      BREMS=asmb/pi*(ff11*C10b**2+ff12*C10b*C20b+ff17*C10b*C70b
     .      +ff18*C10b*C80b+ff22*C20b**2+ff27*C20b*C70b
     .      +ff28*C20b*C80b+ff77*C70b**2+ff78*C70b*C80b
     .      +ff88*C80b**2)

***********************************************************************
*    D: Heavy Quark Effective Theory Corrections

      lambd2=0.12d0  ! = 1/2(M_B*^2-M_B^2)
      aux=MC*(1d0+asmc/pi*(-4d0/3d0))
      HQET=-lambd2/(9d0*aux**2)*(C20b-1d0/6d0*C10b)
     .                  *(KC0+rmu*(KT0+KBSM0))

***********************************************************************
*  E: Finally: The Branching Ratio B -> Xs gamma

*     0) Parameters
      BRSL=0.1061d0

      CCSL=0.580d0
      
*     1) Coefficient
      aux=6d0*ALEM0/(pi*CCSL)*Vbsg*BRSL
      
*     2) Result: BR(b -> Xs gamma)
      BRSG=aux*((Ktot+EPSew)**2+(KIM)**2+BREMS+HQET)
      
*      - Error estimate:
*     First: Variations of BRSG from BSM uncertainties DKBSM:
      BRSG1=aux*((Ktot+DKBSM+EPSew)**2+(KIM)**2+BREMS+HQET)
      BRSG2=aux*((Ktot-DKBSM+EPSew)**2+(KIM)**2+BREMS+HQET)
*     Second: Add the SM uncertainties + 0.23 / - 0.43 (from [8,9]):
      DBRSGmax=(max(BRSG1,BRSG2)-BRSG)+0.23D-4
      DBRSGmin=(min(BRSG1,BRSG2)-BRSG)-0.43D-4
      BRSGmax= BRSG+DBRSGmax
      BRSGmin= BRSG+DBRSGmin

*     Comparison with experimental data:
*      BR(B -> Xs gamma) from [13], 2 sigma bounds:
*      3.04 10^-4 < BR(B -> Xs gamma)=3.55 10^-4 < 4.06 10^-4

      prob(32)=0d0

      IF(BRSGmin.GE.BRSGexpMax)
     .     PROB(32)=BRSGmin/BRSGexpMax-1d0
      IF(BRSGmax.LE.BRSGexpmin)
     .     PROB(32)=BRSGmax/BRSGexpmin-1d0

***********************************************************************
*  F: Constraints from BR(B-->X_s mu+ mu-)

*     1) Wilson Coefficients
      C70=eta**(16/23)*(C70+58d0/135d0*(eta**(-10/23)-1d0)
     .    +29d0/189d0*(eta**(-28/23)-1))
      C70SM=eta**(16/23)*(C70SM
     .      +58d0/135d0*(eta**(-10/23)-1d0)
     .      +29d0/189d0*(eta**(-28/23)-1d0))

      aux=0d0
      do i=1,3
      do j=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .  *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     .       *Intpropa(SMASS(i),WIDTH(i),SMASS(j),WIDTH(j),
     .             1d0/mb1s**2,6d0/mb1s**2)
      enddo
      enddo
      CQ12=(gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
      CQ12=(5.3696d0/s2tw)**2*CQ12*aux
      aux=0d0
      do i=1,2
      do j=1,2
      aux=aux
     .      +(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)
     .      *(-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)
     .       *tanb**2
     .      *Intpropa(PMASS(i),WIDTH(3+i),PMASS(j)
     .    ,WIDTH(3+j),1d0/mb1s**2,6d0/mb1s**2)
      enddo
      enddo
      CQ22=(gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
      CQ22=(5.3696d0/s2tw)**2*CQ22*aux

*     2) Branching Ratio between 1GeV^2 and 6GeV^2
      BRBSll=1.59d-6+BRSL*4d0/CCSL*Vbsg*(ALEM0/(4d0*pi))**2
     .   *((8d0*dlog(6d0)-15d0/mb1s**2+215d0/(3d0*mb1s**6))
     .          *(C70**2-C70SM**2)
     .    +3d0/2d0*(cQ12+cQ22))

      BRBSllmax=(BRBSll-1.59d-6)*1.6d0+(1.59d0+.22d0)*1.d-6
      BRBSllmin=(BRBSll-1.59d-6)*.4d0+(1.59d0-.22d0)*1.d-6

*     3) Comparison with experiment
*     2 sigma bounds from Huber, Hurth, Lunghi, 0712.3009:
*   0.60 10^-6 < BR(B -> Xs l+l-)=1.60 10^-6 < 2.60 10^-6
      prob(40)=0d0
      IF(BRBSllmin.GE.2.6d-6)
     .     PROB(40)=BRBSllmin/2.6d-6-1d0
      IF(BRBSllmax.LE.0.6d-6)
     .     PROB(40)=1d0-BRBSllmax/0.6d-6

*     4) High M_{l+l-} region: 14.4 GeV^2 < M_{l+l-}^2 < mb^2
      aux=0d0
      do i=1,3
      do j=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .  *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     .   *Intpropa(SMASS(i),WIDTH(i),SMASS(j),WIDTH(j),14.4d0/mb1s**2
     .                ,1d0)
      enddo
      enddo
      CQ12=(gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
      CQ12=(5.3696d0/s2tw)**2*CQ12*aux
      
      aux=0d0
      do i=1,2
      do j=1,2
      aux=aux
     .      +(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)
     .      *(-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)
     .       *tanb**2
     .*Intpropa(PMASS(i),WIDTH(3+i),PMASS(j),WIDTH(3+j),14.4d0/mb1s**2
     .                   ,1d0)
      enddo
      enddo
      CQ22=(gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
      CQ22=(5.3696d0/s2tw)**2*CQ22*aux

      BRBSll=2.40d-7+BRSL*4d0/CCSL*Vbsg*(ALEM0/(4d0*pi))**2
     .   *((8d0*dlog(mb1s**2/14.4d0)-4d0*(1d0-14.4d0/mb1s**2)
     .    +4d0/3d0*(1d0-(14.4d0/mb1s**2)**3))*(C70**2-C70SM**2)
     .    +3d0/2d0*(cQ12+cQ22))

      BRBSllmax=(BRBSll-2.40d-7)*1.6d0+2.40d-7*(1d0+2d0*.29d0)
      BRBSllmin=(BRBSll-2.40d-7)*.4d0+2.40d-7*(1d0-2d0*.26d0)

*     2 sigma bounds from Huber, Hurth, Lunghi, 0712.3009:
*     2.0 10^-7 < BR(B -> Xs l+l-)=4.4 10^-7 < 6.8 10^-7

      IF(BRBSllmin.GE.6.8d-7)
     .     PROB(40)=PROB(40)+BRBSllmin/6.8d-7-1d0
      IF(BRBSllmax.LE.2.0d-7)
     .     PROB(40)=PROB(40)+1d0-BRBSllmax/2.0d-7

* For comparison: the Pole contribution only (from Hiller):
      BRBSLL=BRSL/CCSL*Vbsg*3d0*pi/(2d0*dsqrt(2d0)*GF)
     .      *mmu**2*PMASS(1)/(mb1s**8*WIDTH(4))
     .      *(mb1s**2-PMASS(1)**2)**2
     .   *sigRLbs**2*(-PCOMP(1,1)*cosb-PCOMP(1,1)*sinb*tanb)**2
     .   *(PCOMP(1,1)*tanb)**2*(5.3696d0/MB0)**2
      
      return
      END
      
      
**********************************************************************      
**********************************************************************
**********************************************************************

* MSbar 5/6-flavour evolution for the Strong Coupling Constant
      double precision function asf(x)
      
*       ALPHA_S (x) - NB it uses 5 flavors for x<mt and 6 if x>=mt
      implicit none
      double precision x,pi,asc,fn,b0,b1,vvv,b0t,b1t,ast
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      pi=4d0*datan(1d0)
      asc=ALSMZ
      
      fn=5d0
      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      vvv=1d0-b0*asc/(2d0*pi)*dlog(MZ/x)
      asf=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*dlog(vvv))

      if(x.gt.MT)then
       vvv=1d0-b0*asc/(2d0*pi)*dlog(MZ/MT)
       ast=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*dlog(vvv))
       b0t=b0-2d0/3d0
       b1t=b1-38d0/3d0
       vvv=1d0-b0t*ast/(2d0*pi)*dlog(MT/x)
       asf=ast/vvv*(1d0-b1t/b0t*ast/(4d0*pi*vvv)*dlog(vvv))
      endif
      return
      end
            
****************************************************************
* Dilogarithm
      double precision function sp2(x)
*       Li_2(x)
      implicit none
      double precision x,pi,f,bla
      external f
      pi=4d0*datan(1d0)
      bla=0d0
      if(x.ge.-1d0.and.x.le.0.5d0)then
       bla=f(x)
      else if(x.gt.0.5d0.and.x.lt.1d0)then
       bla=-f(1d0-x)+pi**2/6d0-dlog(x)*dlog(1d0-x)
      else if(x.lt.-1d0)then
       bla=-f(1d0/x)-pi**2/6d0-.5d0*(dlog(-x))**2
*      else if(x.gt.-500d0.and.x.lt.-10d0)then
*       bla=1.50616d0+0.153765d0*x-0.0000484249d0*x**2
*     .      -2.69934d-8*x**3
*     .      -1.97807d0*dlog(dabs(x))-0.0245271d0*x*dlog(x)
      else if(x.ge.1d0)then
      bla=0d0
*       write(6,*)'error in dilog',x
      endif
      sp2=bla
      return
      end

****************************************************************

      double precision function sgn(x)
      
      implicit none
      double precision x
      if(x.ge.0d0)then
      sgn=1d0
      else
      sgn=-1d0
      endif
      end
      
****************************************************************

      double precision function f(x)
      
      implicit none
      integer i
      double precision b(12),x,z,cCc,sum
      z=-dlog(1d0-x)
      b(1)=-.5d0
      b(2)=1d0/6d0
      b(3)=0d0
      b(4)=-1d0/30d0
      b(5)=0d0
      b(6)=1d0/42d0
      b(7)=0d0
      b(8)=-1d0/30d0
      b(9)=0d0
      b(10)=5d0/66d0
      b(11)=0d0
      b(12)=-691d0/2730d0
      cCc=z
      sum=z
      do i=1,12
      cCc=cCc*z/(i+1d0)
      sum=sum+b(i)*cCc
      enddo
      f=sum
      end
      
********************************************************************

      double precision function ff1(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      ff1=x*(7d0-5d0*x-8d0*x**2)*d/24d0
     .      +x**2*(3d0*x-2d0)*dg/4d0
      else
      ff1=-5d0/48d0
      endif
      return
      end

*********************************************************************

      double precision function fg1(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-2)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      fg1=x*(2d0+5d0*x-x**2)*d/8d0-3d0*x**2*dg/4d0
      else
      fg1=-1d0/16d0
      endif
      return
      end

*********************************************************************

      double precision function ff2(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-2)then
      d=1d0/(x-1d0)**2
      dg=dlog(x)/(x-1d0)**3
      ff2=x*(3d0-5d0*x)*d/12d0+x*(3d0*x-2d0)*dg/6d0
      else
      ff2=-7d0/36d0
      endif
      return
      end

*********************************************************************

      double precision function fg2(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-2)then
      d=1d0/(x-1d0)**2
      dg=dlog(x)/(x-1d0)**3
      fg2=x*(3d0-x)*d/4d0-x*dg/2d0
      else
      fg2=-1d0/6d0
      endif
      return
      end

*********************************************************************

      double precision function ff3(x)

      implicit none
      double precision x,ff1,ff2
      ff3=2d0/3d0*(1d0-1d0/x)*ff1(x)+ff2(x)+23d0/36d0
      return
      end
*********************************************************************

      double precision function fg3(x)

      implicit none
      double precision x,fg1,fg2
      fg3=2d0/3d0*(1d0-1d0/x)*fg1(x)+fg2(x)+1d0/3d0
      return
      end

*********************************************************************

      double precision function esm(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      esm=x*(-18d0+11d0*x+x**2)*d/12d0
     .      +x**2*(15d0-16d0*x+4d0*x**2)*dg/6d0
      esm=esm-2d0/3d0*dlog(x)
      else
      esm=43d0/72d0
      endif
      return
      end

*********************************************************************

      double precision function eh(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1.d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      eh=x*(16d0-29d0*x+7d0*x**2)*d/36d0+x*(3d0*x-2d0)
     .       *dg/6d0
      else
      eh=1d0/8d0
      endif
      return
      end
      
*********************************************************************

      double precision function af(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      af=16d0/9d0*((5d0/2d0-1d0/3d0*pi**2-3d0*1.2021d0
     .     +(5d0/2d0-3d0/4d0*pi**2)*lnx+lnx**2/4d0
     .       +lnx**3/12d0)*x
     .    +(7d0/4d0+2d0/3d0*pi**2-1d0/2d0*pi**2*lnx-lnx**2/4d0
     .     +lnx**3/12d0)*x2
     .    +(-7d0/6d0-pi**2/4d0+2d0*lnx-3d0/4d0*lnx**2)*x3)
      return
      end

*********************************************************************

      double precision function afim(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      afim=16d0/9d0*pi*((2d0-pi**2/6d0+lnx/2d0+lnx**2/2d0)*x
     .  +(1d0/2d0-pi**2/6d0-lnx+lnx**2/2d0)*x2+x3)
      return
      end

*********************************************************************

      double precision function bf(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      bf=-8d0/9d0*((-3d0+pi**2/6d0-lnx)*x
     .    +(1d0/2d0+pi**2-2d0*lnx-lnx**2/2d0)*x2
     .    +(-25d0/12d0-pi**2/9d0-19d0/18d0*lnx+2d0*lnx**2)*x3
     .    -2d0/3d0*pi**2*x**(3d0/2d0))
      return
      end

*********************************************************************

      double precision function bfim(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      bfim=-8d0/9d0*pi*(-x+(1d0-2d0*lnx)*x2
     .       +(-10d0/9d0+4d0/3d0*lnx)*x3)
      return
      end

*********************************************************************
* NLO t-quark contribution to K
      double precision function A1(x)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      spx=sp2(1d0-1d0/x)
      A1=(32d0*x4+244d0*x3-160d0*x2+16d0*x)/(9d0*dd4)*spx
     .   -(-774d0*x4-2826d0*x3+1994d0*x2-130d0*x+8d0)
     .            /(81d0*dd5)*lnx
     .   +(-94d0*x4-18665d0*x3+20682d0*x2-9113d0*x+2006d0)
     .            /(243d0*dd4)
      return
      end

*********************************************************************
* NLO t-quark contribution to K
      double precision function F1(x)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      spx=sp2(1d0-1d0/x)
      F1=(4d0*x4-40d0*x3-41d0*x2-x)/(3d0*dd4)*spx
     .   -(-144d0*x4+3177d0*x3+3661d0*x2+250d0*x-32d0)
     .   /(108d0*dd5)*lnx
     .   +(-247d0*x4+11890d0*x3+31779d0*x2-2966d0*x+1016d0)
     .   /(648d0*dd4)
      return
      end
      
*********************************************************************
* NLO t-quark contribution to K
      double precision function dA0(x)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      spx=sp2(1d0-1d0/x)
      dA0=x*(3d0*x2+5d0*x-4)/(2d0*dd5)*lnx
     .   +(-141d0*x2+48d0*x+21d0)/(36d0*dd4)
      return
      end
      
**********************************************************************
* NLO t-quark contribution to K
      double precision function dF0(x)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      spx=sp2(1d0-1d0/x)
      dF0=-3d0*x*(x+1d0)/dd5*lnx
     .   +(x2+10d0*x+1d0)/(2d0*dd4)
      return
      end
            
*********************************************************************

      double precision function ffh(x,tanb)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd3,dd4,dd5,lnx,ln2x,spx
      double precision sp2,sum1,sum2,tanb
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd3=(x-1d0)**3
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      ln2x=lnx**2
      spx=sp2(1d0-1d0/x)
      sum1=4d0*(-3d0+7d0*x-2d0*x2)/(3d0*dd3)*spx+
     .       (8d0-14d0*x-3d0*x2)/(3d0*dd4)*ln2x+
     .       2d0*(-3d0-x+12d0*x2-2d0*x3)/(3d0*dd4)*lnx
     .       +(7d0-13d0*x+2d0*x2)/dd3
      sum2=x*(18d0-37d0*x+8d0*x2)/dd4*spx
     .      +x*(-14d0+23d0*x+3d0*x2)/dd5*ln2x
     .      +(-50d0+251d0*x-174d0*x2-192d0*x3+21d0*x4)
     .      /(9d0*dd5)*lnx+
     .      (797d0-5436d0*x+7569d0*x2-1202d0*x3)/(108d0*dd4)
      ffh=-4d0/3d0*x*sum1+2d0/9d0*x/(tanb**2)*sum2
      return
      end

***********************************************************************

      double precision function fgh(x,tanb)
      
      implicit none
      double precision x,x2,x3,x4,dd2,dd3,dd4,dd5,lnx,ln2x,sp2
      double precision spx,sum1,sum2,tanb
      x2=x**2
      x3=x**3
      x4=x**4      
      dd2=(x-1d0)**2
      dd3=(x-1d0)**3
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      ln2x=lnx**2
      spx=sp2(1d0-1d0/x)
      sum1=(-36d0+25d0*x-17d0*x**2)/(2d0*dd3)*spx+(19d0+17d0*x)
     .     /(dd4)*ln2x+(-3d0-187d0*x+12d0*x2-14d0*x3)/(4d0*dd4)*lnx
     .      +3d0*(143d0-44d0*x+29d0*x2)/(8d0*dd3)
      sum2=x*(30d0-17d0*x+13d0*x2)/dd4*spx
     .      -x*(31d0+17d0*x)/dd5*ln2x
     .      +(-226d0+817d0*x+1353d0*x2+318d0*x3+42d0*x4)
     .       /(36d0*dd5)*lnx
     .      +(1130d0-18153d0*x+7650d0*x2-4451d0*x3)/(216d0*dd4)
      fgh=-x/3d0*sum1+x/(6d0*tanb**2)*sum2
      return
      end

**********************************************************************

* For effective Yukawa Couplings
      double precision function h2(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-y).gt.1.d-2)then
        if(dabs(x-1d0).lt.1.d-2)then
             h2=1d0/(-1d0+y)-y*dlog(y)/(-1d0+y)**2
         else
            if(dabs(y-1d0).lt.1.d-2)then
             h2=1d0/(-1d0+x) - x*dlog(x)/(-1d0+x)**2
            else
             h2=x*dlog(x)/((1d0-x)*(x-y))+
     .      y*dlog(y)/((1d0-y)*(-x+y))
            endif
         endif
      else
         if(dabs(x-1d0).lt.1.d-2)then
            h2=-1d0/2d0
         else
            h2=1d0/(1d0 - x) + dlog(x)/(-1d0 + x)**2
         endif
      endif
      return
      end

*******************************************************************

      double precision function BB0(x,y,z,t)

      implicit none
      double precision x,y,z,t,NMB0
      BB0=1d0-NMB0(x,y,z,t)
      return
      end

*******************************************************************

      double precision function BB1(x,y,z,t)

      implicit none
      double precision x,y,z,t,BB0
      if(dabs(y-z).lt.1.d-5)then
      BB1=1d0/2d0*BB0(x,y,z,t)
      else
      BB1=1d0/2d0*BB0(x,y,z,t)-(y+z)/(4d0*(y-z))
     .      +y*z/(2d0*(y-z)**2)*dlog(y/z)
      endif
      return
      end

*******************************************************************

           double precision function S0(x)

      implicit none
      double precision x
      S0=x*(1d0/4d0+9d0/(4d0*(1d0-x))-3d0/(2d0*(1d0-x)**2)
     .       -3d0*x**2*dlog(x)/(2d0*(1d0-x)**3))
      return
      end
      
*******************************************************************

      double precision function gg0(x)

      implicit none
      double precision x
      if(dabs(x-1d0).lt.1.d-5)then
      gg0=1d0/12d0
      else
      gg0=x*(x**2-1d0-2d0*x*dlog(x))/(4d0*(x-1d0)**3)
      endif
      return
      end
      
*******************************************************************

      double precision function ffp(x,y)

      implicit none
      double precision x,y
      if(dabs(y-x).lt.1.d-5)then
       if(dabs(y-1d0).lt.1.d-5)then
       ffp=1d0/6d0
       else
       ffp=(-1d0+x**2-2d0*x*dlog(x))
     .      /(-2d0*x+6d0*x**2+2d0*x**4-6d0*x**3)
       endif
      else
       if(dabs(y-1d0).lt.1.d-5)then
       ffp=(-2d0*x+x*dlog(x)+2d0+dlog(x))/(x-1d0)**2
       else
        if(dabs(x-1d0).lt.1.d-5)then
        ffp=(-1d0*y**2-2d0*y*dlog(y))/(2d0*(y-1d0)**3)
        else
        ffp=(x**2-y)*dlog(x)/((x-y)**2*(x-1d0)**2)
     .      -y*dlog(y)/((x-y)**2*(y-1d0))-1d0/((x-y)*(x-1d0))
        endif
       endif
      endif
      return
      end
      
*******************************************************************

      double precision function ggp(x,y,z)

      implicit none
      double precision x,y,z
      if(dabs(z-1d0).le.1.d-5)then
        if(dabs(x-y).le.1.d-5)then
         if(dabs(1d0-x).le.1.d-5)then
         ggp=1d0/3d0
         else
         ggp=(x**2-2d0*dlog(x)*x-1d0)/(3d0*x-1d0+x**3-3d0*x**2)
         endif
        else
         if(dabs(x-1d0).le.1.d-5)then
         ggp=(2d0*y**2*dlog(x)+4d0*y-1d0-3d0*y**2)
     .   /(-2d0+6d0*y+2d0*y**3-6d0*y**2)
         else if(dabs(1d0-y).le.1.d-5)then
         ggp=(2d0*x**2*dlog(x)+4d0*x-1d0-3d0*x**2)
     .   /(-2d0+6d0*x+2d0*x**3-6d0*x**2)
         else
         ggp=(x**2*dlog(x)-y**2*dlog(y)-dlog(y)*y**2*x**2-x**2
     .       -2d0*dlog(x)*y*x**2-x*y**2+2d0*dlog(y)*y**2*x+x**2*y
     .       +x**2*dlog(x)*y**2+x+y**2-y)
         endif
        endif
      else
        if(dabs(x-y).le.1.d-5)then
         if(dabs(x-1).le.1.d-5)then
         ggp=(-1d0-3d0*z**2+4d0*z+2d0*z**2*dlog(z))
     .  /(-2d0-6d0*z**2+2d0*z**3+6d0*z)
         else
         ggp=x/((x-z)*(x-1d0))*(1d0-(1d0/(x-1d0)
     .   +z/(x-z))*dlog(x))+z**2*dlog(z)/((x-z)**2*(z-1d0))
         endif
        else if(abs(x-z).le.1.d-5)then
         if(abs(y-1d0).le.1.d-5)then
         ggp=(x**2-2d0*dlog(x)*x-1d0)/(x**3-3d0*x**2+3d0*x-1d0)
         else
         ggp=x/((x-y)*(y-1d0))*(1d0-(1d0/(y-1d0)+y/(x-y))*dlog(x))
     .   +y**2*dlog(y)/((x-y)**2*(y-1d0))      
         endif
        else if(abs(y-z).le.1.d-5)then
         if(abs(x-1d0).le.1.d-5)then
         ggp=(y**2-2d0*dlog(y)*y-1d0)/(y**3-3d0*y**2+3d0*y-1d0)
         else
         ggp=y/((y-x)*(x-1d0))*(1d0-(1d0/(x-1d0)+x/(y-x))*dlog(y))
     .   +x**2*dlog(x)/((y-x)**2*(x-1d0))
         endif
        else
         if(dabs(y-1d0).le.1.d-5)then
         ggp=(-x*z**2+x**2*dlog(x)+x**2*dlog(x)*z**2
     .  +2d0*z**2*dlog(z)*x-z-z**2*dlog(z)-z**2*dlog(z)*x**2
     .  -x**2+x**2*z-2d0*x**2*dlog(x)*z+z**2+x)
     .       /(2d0*x*z**3+3d0*x**2*z-3d0*x*z**2-2d0*x**2-z**3
     .       +2d0*z**2+x-z+x**3*z**2-x**2*z**3-2d0*x**3*z+x**3)
         else if(dabs(x-1d0).le.1.d-5)then
         ggp=(-y*z**2+y**2*dlog(y)+y**2*dlog(y)*z**2
     .  +2d0*z**2*dlog(z)*y-z-z**2*dlog(z)-z**2*dlog(z)*y**2
     .  -y**2+y**2*z-2d0*y**2*dlog(y)*z+z**2+y)
     .       /(2d0*y*z**3+3d0*y**2*z-3d0*y*z**2-2d0*y**2-z**3
     .       +2d0*z**2+y-z+y**3*z**2-y**2*z**3-2d0*y**3*z+y**3)
         else
         ggp=1d0/(x-y)*(1d0/(x-z)*(x**2/(x-1d0)*dlog(x)-3d0/2d0*x
     .           -z**2/(z-1d0)*dlog(z)+3d0/2d0*z)
     .       -1d0/(y-z)*(y**2/(y-1d0)*dlog(y)-3d0/2d0*y
     .           -z**2/(z-1d0)*dlog(z)+3d0/2d0*z))
         endif
        endif
      endif
      return
      end
      
*******************************************************************

      double precision function RUNMass(y,x)

      implicit none
      double precision x,y
      DOUBLE PRECISION PI,nf,ge0,ge1,b0,b1,aux,aux1,asf
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      PI=4d0*DATAN(1d0)

*       * x: renormalization scale
*       * y: quark running mass at 2 GeV;
*      Running Masses at 2 GeV (PDG 2006):
*      - Mu= 1.5 to 3.0 MeV
*      - Md= 3. to 7. MeV
*      - Ms= 95 +/- 25 MeV
*      - Mc= 1.25 +/- 0.09 GeV
*      - Mb= 4.20 +/- 0.07 GeV
      
      
      aux=y

       nf=4d0
       ge0=8d0
       ge1=404d0/3d0-38d0/3d0*nf
       b0=11d0-2d0*nf/3d0
       b1=102d0-38d0/3d0*nf
       aux1=2d0
      
       if(x.gt.MBP)then
        aux=aux*(asf(MBP)/asf(MC))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(MBP)-asf(aux1)))

        nf=5d0
        ge0=8d0
        ge1=404d0/3d0-38d0/3d0*nf
        b0=11d0-2d0*nf/3d0
        b1=102d0-38d0/3d0*nf
        aux1=MBP
      
        if(x.gt.MT)then
         aux=aux*(asf(MT)/asf(MBP))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(MT)-asf(aux1)))

         nf=6d0
         ge0=8d0
         ge1=404d0/3d0-38d0/3d0*nf
         b0=11d0-2d0*nf/3d0
         b1=102d0-38d0/3d0*nf
         aux1=MT
        endif
      
       endif
      
      
      RUNMass=aux*(asf(x)/asf(aux1))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(x)-asf(aux1)))
      
      
      return
      end
      
*******************************************************************

* Additional function Intpropa for BR(B-->Xs l+l-):

      double precision function Intpropa(mH1,Gam1,mH2,Gam2,bd1,bd2)

      implicit none
      double precision mH1,Gam1,mH2,Gam2,bd1,bd2
      double precision mHmin,mHmax,Gam,mb1s,aux,fac,fca
      
      MB1S=4.68d0
      if(mH1.gt.mH2)then
      mHmax=mH1
      mHmin=mH2
      Gam=Gam2
      else
      mHmax=mH2
      mHmin=mH1
      Gam=Gam1
      endif
      
      fac=dsqrt(mHmin*mHmax)
      IF(fac.gt.150d0)then
       aux=(bd2-bd1)/2*(bd1+bd2-4d0/3d0*(bd2**2+bd1*bd2+bd1**2)
     .    +(bd2**3+bd2**2*bd1+bd2*bd1**2+bd1**3)/2d0)*(mb1s/fac)**4
      ELSE
       fac=(mb1s/mHmin)**2
       fca=(mb1s/mHmax)**2
       if(abs(mHmin-mHmax).gt.1.d-4)then
        if(mHmin.gt.5d0)then
        aux=(2d0*(bd2-bd1)
     .  +mHmin**2/(mHmin**2-mHmax**2)*(1-fac)**2
     .  *dlog((1d0-fac*bd2)**2/(1d0-fac*bd1)**2)/fac)/(2d0*fac)
     .      +(2d0*(bd2-bd1)
     .  +mHmax**2/(mHmax**2-mHmin**2)*(1-fca)**2
     .  *dlog((1d0-fca*bd2)**2/(1d0-fca*bd1)**2)/fca)/(2d0*fca)
     .      +(bd2-bd1)*(bd2+bd1-4d0)/2d0
        else
        aux=(2d0*(bd2-bd1)+mHmin**2*(mHmin**2-mHmax**2)
     .      /((mHmin**2-mHmax**2)**2+(mHmin*Gam)**2)
     .      *dlog(((1-fac*bd2)**2+(Gam/mHmin)**2)
     .      /((1-fac*bd1)**2+(Gam/mHmin)**2))/fac
     .      *((1-fac)**2+mHmax**2*Gam**2
     .      /(mHmin**2*(mHmin**2-mHmax**2))
     .      *(3d0+(1d0-fac*Gam**2/mb1s**2)*fac*fca-2d0*fac-2*fca)
     .                        ))/(2d0*fac)
        aux=aux+(2d0*(bd2-bd1)+mHmax**2*(mHmax**2-mHmin**2)
     .      /((mHmin**2-mHmax**2)**2+(mHmin*Gam)**2)
     .      *dlog((1-fca*bd2)**2/(1-fca*bd1)**2)/fca
     .      *(1-fca)**2)/(2d0*fca)
        aux=aux+Gam*mHmin*mHmax**2/fac
     .    /((mHmax**2-mHmin**2)**2+(mHmin*Gam)**2)
     . *(datan((mHmin**2-mb1s**2*bd2)/(mHmin*Gam))
     .      -datan((mHmin**2-mb1s**2*bd1)/(mHmin*Gam)))
     .   *(3d0-4d0*fac+fac**2+2d0*fca-2*(mHmin/mHmax)**2
     .     -(Gam/mb1s)**2*(2d0*fca+fac-2d0*fac*fca))
        aux=aux+(bd2-bd1)*(bd1+bd2-4d0)/2d0
        endif
       else
        if(mHmin.gt.5d0)then
        aux=(2d0+(1d0-fac)**2/((1d0-fac*bd2)*(1d0-fac*bd1)))
     .       *(bd2-bd1)+(3d0-4d0*fac+fac**2)/2d0
     .      *dlog((1d0-fac*bd2)**2/((1d0-fac*bd1)**2))/fac
        aux=aux/fac+(bd2-bd1)/2d0*(bd2+bd1-4d0)
       else
        aux=-(bd2-bd1)*(bd2+bd1-4d0)/2d0-(2d0*(bd2-bd1)
     . +(3d0-4d0*fac**2+(1d0+mHmin**2*Gam**2/mb1s**4)*fac**2)
     .  *dlog(((1-fac*bd2)**2+(Gam/mHmin)**2)
     .   /((1-fac*bd1)**2+(Gam/mHmin)**2))/(2d0*fac))/fac
     .      +(datan((mb1s**2*bd2-mHmin**2)/(mHmin*Gam))
     .       -datan((mb1s**2*bd1-mHmin**2)/(mHmin*Gam)))
     .       /(fac**3*Gam*mHmin/mb1s**2)
     .  *((1d0-fac)**2
     .     +(mHmin*Gam/mb1s**2)**2*fac**2*(-3d0+2d0*fac))
       endif
      endif
      ENDIF
      aux=aux/mb1s**4
      intpropa=aux
      return
      end
