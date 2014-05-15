      SUBROUTINE FILLGAMBR2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,NCMAX,NHC_H,SHC_H,CHC_H,STMASS_H
     . ,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     . ,MC_H,MN_H,NMNH,GAMBRN,NMCH,GAMBRC)
************************************************************************
*
*JSL 10/Jun/2009: Input arrays are extended to include SMPARA_H(NSMIN)
*                 and SSPARA_H(NSSIN)
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
* Input Array:
      REAL*8     SMPARA_H(NSMIN),SSPARA_H(NSSIN)
      INTEGER*8  IFLAG_H(NFLAG)
      REAL*8     HMASS_H(3)
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      REAL*8     MC_H(2),MN_H(4)
*-----------------------------------------------------------------------
* Output Array:
      REAL*8 GAMBRN(NMNH,3,3)   
      REAL*8 GAMBRC(NMCH,3)      
*-----------------------------------------------------------------------
* For integration
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
      EXTERNAL FVVS,FHVS
*-----------------------------------------------------------------------
* Local :
      INTEGER*8 ISMN,ISUSYN,ISMC,ISUSYC
      COMPLEX*16 GHV,GSS
      COMPLEX*16 XI
      COMPLEX*16 GF_UD,GS_UD,GP_UD
      COMPLEX*16 GF_CS,GS_CS,GP_CS
      COMPLEX*16 GF_TB,GS_TB,GP_TB
      COMPLEX*16 GF1,GS1,GP1
      COMPLEX*16 GF2,GS2,GP2
      COMPLEX*16 GF3,GS3,GP3
      COMPLEX*16 GF4,GS4,GP4
* Radiative corrections to Htt Hbb Yukawa couplings
      COMPLEX*16 HB_H,HT_H,CKB_H,CKT_H
      COMPLEX*16 RB_H,RT_H,CKBB_H,CKBT_H
*JSL 10/Jun/2009: Including threshold corrections
      COMPLEX*16 CKD_H,CKS_H
* Higgs-photon-photon
      COMPLEX*16 SPH,PPH
      COMPLEX*16 SPP(15),PPP(7)
*-----------------------------------------------------------------------
      ISMN   = IFLAG_H(20)
      ISUSYN = IFLAG_H(21)
      ISMC   = IFLAG_H(22)
      ISUSYC = IFLAG_H(23)
*      print*,ismn,isusyn,ismc,isusyc
      PI=2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
*
* << NEUTRAL HIGGS BOSON DECAYS INTO SM PARTICLES >>

      DO IH=1,3
      MH    = HMASS_H(IH)
*=======================================================================
      CALL QMASS_MH(MH,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH)
      CALL QMASS_MH(MH/2.D0,ASMH2,MTMH2,MBMH2,MCMH2
     .                           ,MSMH2,MUMH2,MDMH2)
*=======================================================================
*
*---> H_IH -> e+ e-     [IM= 1]
      IM   = 1
      IFF  = 1
      MJ   = ME_H
      MK   = ME_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> mu+ mu-   [IM= 2]
      IM   = 2
      IFF  = 4
      MJ   = MMU_H
      MK   = MMU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> tau+ tau- [IM= 3]
      IM   = 3
      IFF  = 7
      MJ   = MTAU_H
      MK   = MTAU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
                 
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> d d       [IM= 4]
      IM   = 4
      IFF  = 10
      MJ   = MDMH
      MK   = MDMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MDMH/MDMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> s s       [IM= 5]
      IM   = 5
      IFF  = 13
      MJ   = MSMH
      MK   = MSMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MSMH/MSMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> b b       [IM= 6]
      IM   = 6
      IFF  = 16
      MJ   = MBMH
      MK   = MBMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MBMH/MBMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*      print*,'SAME?',GAMBRN(IM,1,IH)
*recalculate taking into account running b-quark mass
*      IF((MH-2.D0*MBMH).GT.0.D0) THEN
*       BETABB=DSQRT(1.D0-4.D0*MBMH**2/MH**2)
*      ELSE
*       BETABB=0.D0
*      ENDIF
*      GBB=3.D0*(GW_H*MBMH/2.D0/MW_H)**2*MH*BETABB/8.D0/PI
*     .   *(BETABB**2*CDABS(NHC_H(IFF+1,IH))**2
*     .              +CDABS(NHC_H(IFF+2,IH))**2)
*      print*,'SAME?',GBB
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> u u       [IM= 7]
      IM   = 7
      IFF  = 19
      MJ   = MUMH
      MK   = MUMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MUMH/MUMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> c c       [IM= 8]
      IM   = 8
      IFF  = 22
      MJ   = MCMH
      MK   = MCMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MCMH/MCMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> t t       [IM= 9]
      IM   = 9
      IFF  = 25
      MJ   = MTMH
      MK   = MTMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MTMH/MTMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
* 
*      print*,'H',ih,' -> e   e   : ',gambrn(1,1,ih)
*      print*,'H',ih,' -> mu  mu  : ',gambrn(2,1,ih)
*      print*,'H',ih,' -> tau tau : ',gambrn(3,1,ih)
*      print*,'H',ih,' -> d   d   : ',gambrn(4,1,ih)
*      print*,'H',ih,' -> s   s   : ',gambrn(5,1,ih)
*      print*,'H',ih,' -> b   b   : ',gambrn(6,1,ih)
*      print*,'H',ih,' -> u   u   : ',gambrn(7,1,ih)
*      print*,'H',ih,' -> c   c   : ',gambrn(8,1,ih)
*      print*,'H',ih,' -> t   t   : ',gambrn(9,1,ih)
*
*---> H_IH -> W W       [IM=10]
      IM  = 10
      IVV = 70
      DV  = 2.D0
      CALL HVV(GF_H,DV,NHC_H(IVV,IH),MW_H,HMASS_H(IH),GAMBRN(IM,1,IH))
* into W + W^* 
      DVVS=2.D0
      IF( (HMASS_H(IH).GT.MW_H) .AND. 
     .    (HMASS_H(IH).LT.(2.D0*MW_H+GAMW_H)) ) THEN
       IF(DABS(HMASS_H(IH)-2.D0*MW_H).LT.GAMW_H) THEN
        DVVS=2.D0-(HMASS_H(IH)-2.D0*MW_H+GAMW_H)/2.D0/GAMW_H
       ENDIF
      EPSV=GAMW_H/MW_H
      OMEGAI=HMASS_H(IH)**2/MW_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*DV*DVVS
     . *CDABS(NHC_H(IVV,IH))**2
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES ! 3body
*
      CALL WW4BODY(HMASS_H(IH),WSWS)
      GAMBRN(IM,1,IH)=CDABS(NHC_H(IVV,IH))**2*WSWS ! 4body
*
      ENDIF
*
*---> H_IH -> Z Z       [IM=11]
      IM  = 11
      IVV = 70
      DV  = 1.D0
      CALL HVV(GF_H,DV,NHC_H(IVV,IH),MZ_H,HMASS_H(IH),GAMBRN(IM,1,IH))
* into Z + Z^* 
      DVVS=2.D0
      IF( (HMASS_H(IH).GT.MZ_H) .AND. 
     .    (HMASS_H(IH).LT.(2.D0*MZ_H+GAMZ_H)) ) THEN
       IF(DABS(HMASS_H(IH)-2.D0*MZ_H).LT.GAMZ_H) THEN
        DVVS=2.D0-(HMASS_H(IH)-2.D0*MZ_H+GAMZ_H)/2.D0/GAMZ_H
       ENDIF
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*DV*DVVS
     . *CDABS(NHC_H(IVV,IH))**2
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES ! 3body
*
      CALL ZZ4BODY(HMASS_H(IH),ZSZS)
      GAMBRN(IM,1,IH)=CDABS(NHC_H(IVV,IH))**2*ZSZS ! 4body
*
      ENDIF
*      print*,'H',ih,' -> W   W   : ',gambrn(10,1,ih)
*      print*,'H',ih,' -> Z   Z   : ',gambrn(11,1,ih)
*---> H_IH -> H1 Z      [IM=12]
      IM  = 12
      IF(IH.EQ.1) GHV=DCMPLX(0.D0,0.D0) ! Neglecting overall sign
      IF(IH.EQ.2) GHV=NHC_H(70,3)       
      IF(IH.EQ.3) GHV=NHC_H(70,2)
      CALL HHV(GF_H,GHV,HMASS_H(IH),HMASS_H(1),MZ_H,GAMBRN(IM,1,IH))
* H_IH -> H1 Z*
      IF( (HMASS_H(IH).GT.HMASS_H(1)) .AND. 
     .    (HMASS_H(IH).LT.(HMASS_H(1)+MZ_H+5.D0*GAMZ_H)) ) THEN
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=HMASS_H(1)**2/MZ_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*CDABS(GHV)**2
     . /8.D0/SQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*---> H_IH -> H2 Z      [IM=13]
      IM  = 13
      IF(IH.EQ.1) GHV=NHC_H(70,3)       ! Neglecting overall sign
      IF(IH.EQ.2) GHV=DCMPLX(0.D0,0.D0) 
      IF(IH.EQ.3) GHV=NHC_H(70,1)       
      CALL HHV(GF_H,GHV,HMASS_H(IH),HMASS_H(2),MZ_H,GAMBRN(IM,1,IH))
* H_IH -> H2 Z*
      IF( (HMASS_H(IH).GT.HMASS_H(2)) .AND. 
     .    (HMASS_H(IH).LT.(HMASS_H(2)+MZ_H+5.D0*GAMZ_H)) ) THEN
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=HMASS_H(2)**2/MZ_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*CDABS(GHV)**2
     . /8.D0/SQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*      print*,'H',ih,' -> H1  Z   : ',gambrn(12,1,ih)
*      print*,'H',ih,' -> H2  Z   : ',gambrn(13,1,ih)
*---> H_IH -> H1 H1     [IM=14]
      IM   = 14
      SYMF = 2.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(10),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(9),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(6),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(1),HMASS_H(1)
     .,GAMBRN(IM,1,IH))
*---> H_IH -> H1 H2     [IM=15]
      IM   = 15
      SYMF = 1.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(9),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(8),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(5),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(1),HMASS_H(2)
     .,GAMBRN(IM,1,IH))
*---> H_IH -> H2 H2     [IM=16]
      IM   = 16
      SYMF = 2.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(8),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(7),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(4),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(2),HMASS_H(2)
     .,GAMBRN(IM,1,IH))
*      print*,'H',ih,' -> H1  H1  : ',gambrn(14,1,ih)
*      print*,'H',ih,' -> H1  H2  : ',gambrn(15,1,ih)
*      print*,'H',ih,' -> H2  H2  : ',gambrn(16,1,ih)
*---> H_IH -> P P       [IM=17]
      CALL HPP(IH,ASMH2,MTMH2,MBMH2,MCMH2
     .        ,MC_H,MCH,HMASS_H,STMASS_H,SBMASS_H,STAUMASS_H
     .        ,NCMAX,NHC_H,SPH,PPH,SPP,PPP)
*      print*,'---> FILLGAMBR'
*      print*,'H-P-P with IH = ',ih
*      print*,'bottom,top     ',spp(1),spp(2)
*      print*,'charm,tau      ',spp(3),spp(4)
*      print*,'c.ino11,c.ino22',spp(5),spp(6)
*      print*,'stop11,stop22  ',spp(7),spp(8)
*      print*,'sbot11,sbot22  ',spp(9),spp(10)
*      print*,'ww,c.hc.h      ',spp(11),spp(12)
*      print*,'stau11,stau22  ',spp(13),spp(14)
*      print*,'S(H-p-p)       ',nhc_h(88,ih)
*      print*,'zero?          ',nhc_h(88,ih)-spp(15)
*      print*,'zero?          ',nhc_h(88,ih)-sph
*      print*,'bottom,top     ',ppp(1),ppp(2)
*      print*,'charm,tau      ',ppp(3),ppp(4)
*      print*,'c.ino11,c.ino22',ppp(5),ppp(6)
*      print*,'P(H-p-p)       ',nhc_h(89,ih)
*      print*,'zero?',nhc_h(89,ih)-ppp(7)
*      print*,'zero?',nhc_h(89,ih)-pph
*
      DJT = -ASMH2/PI
      DJSQ= 8.D0*ASMH2/3.D0/PI
*      print*,DJT,DJSQ
*      print*,SPH
      SPH=SPH+DJT*SPP(2)+DJSQ*(SPP(7)+SPP(8)+SPP(9)+SPP(10))
*      print*,SPH
*
*JSL[12/JUL/10] alpha(Q^2=0) used
      AEM_0=1.D0/137.D0
      GAMBRN(17,1,IH)=HMASS_H(IH)**3*AEM_0**2/256.D0/PI**3/V_H**2
     . *(CDABS(SPH)**2+CDABS(PPH)**2)
* ---> H_IH -> G G       [IM=18]
* running alpha_s(m_h) effect and K factor included
      IF(                        HMASS_H(IH).GT.MTMH) XNF=6.D0
      IF(HMASS_H(IH).LE.MTMH.AND.HMASS_H(IH).GT.MBMH) XNF=5.D0
      IF(HMASS_H(IH).LE.MBMH.AND.HMASS_H(IH).GT.MCMH) XNF=4.D0
      IF(HMASS_H(IH).LE.MCMH.AND.HMASS_H(IH).GT.MSMH) XNF=3.D0
      IF(HMASS_H(IH).LE.MSMH.AND.HMASS_H(IH).GT.MDMH) XNF=2.D0
      IF(HMASS_H(IH).LE.MDMH                        ) XNF=1.D0
      CKHG=1.D0+ASMH/PI*(95.D0/4.D0-7.D0/6.D0*XNF)
      CKAG=1.D0+ASMH/PI*(97.D0/4.D0-7.D0/6.D0*XNF)
      GAMBRN(18,1,IH)=HMASS_H(IH)**3*ASMH**2/32.D0/PI**3/V_H**2
     . *(CKHG*CDABS(NHC_H(84,IH))**2+CKAG*CDABS(NHC_H(85,IH))**2)
*      print*,IH,MH,ASMH,MBMH
*      print*,'H',ih,' -> P  P  : ',gambrn(17,1,ih)
*      print*,'H',ih,' -> G  G  : ',gambrn(18,1,ih)
*-----------------------------------------------------------------------
*
* << NEUTRAL HIGGS BOSON DECAYS INTO SUSY PARTICLES >>
 
*---> H_IH -> N1 N1       [IM=ISMN+1]
      IM   = ISMN+1
      IFF  = 28
      MJ   = MN_H(1)
      MK   = MN_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N2       [IM=ISMN+2]
      IM   = ISMN+2
      IFF  = 40
      MJ   = MN_H(1)
      MK   = MN_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N3       [IM=ISMN+3]
      IM   = ISMN+3
      IFF  = 43
      MJ   = MN_H(1)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N4       [IM=ISMN+4]
      IM   = ISMN+4
      IFF  = 46
      MJ   = MN_H(1)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N2       [IM=ISMN+5]
      IM   = ISMN+5
      IFF  = 31
      MJ   = MN_H(2)
      MK   = MN_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N3       [IM=ISMN+6]
      IM   = ISMN+6
      IFF  = 49
      MJ   = MN_H(2)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N4       [IM=ISMN+7]
      IM   = ISMN+7
      IFF  = 52
      MJ   = MN_H(2)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N3 N3       [IM=ISMN+8]
      IM   = ISMN+8
      IFF  = 34
      MJ   = MN_H(3)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N3 N4       [IM=ISMN+9]
      IM   = ISMN+9
      IFF  = 55
      MJ   = MN_H(3)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N4 N4       [IM=ISMN+10]
      IM   = ISMN+10
      IFF  = 37
      MJ   = MN_H(4)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> N1 N1 : ',gambrn(ismn+1,1,ih)
*      print*,'H',IH,' -> N1 N2 : ',gambrn(ismn+2,1,ih)
*      print*,'H',IH,' -> N1 N3 : ',gambrn(ismn+3,1,ih)
*      print*,'H',IH,' -> N1 N4 : ',gambrn(ismn+4,1,ih)
*      print*,'H',IH,' -> N2 N2 : ',gambrn(ismn+5,1,ih)
*      print*,'H',IH,' -> N2 N3 : ',gambrn(ismn+6,1,ih)
*      print*,'H',IH,' -> N2 N4 : ',gambrn(ismn+7,1,ih)
*      print*,'H',IH,' -> N3 N3 : ',gambrn(ismn+8,1,ih)
*      print*,'H',IH,' -> N3 N4 : ',gambrn(ismn+9,1,ih)
*      print*,'H',IH,' -> N4 N4 : ',gambrn(ismn+10,1,ih)
*
*---> H_IH -> C1+ C1-       [IM=ISMN+11]
      IM   = ISMN+11
      IFF  = 58
      MJ   = MC_H(1)
      MK   = MC_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C1+ C2-       [IM=ISMN+11]
      IM   = ISMN+12
      IFF  = 64
      MJ   = MC_H(1)
      MK   = MC_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C2+ C1-       [IM=ISMN+11]
      IM   = ISMN+13
      IFF  = 61
      MJ   = MC_H(2)
      MK   = MC_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C2+ C2-       [IM=ISMN+11]
      IM   = ISMN+14
      IFF  = 67
      MJ   = MC_H(2)
      MK   = MC_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> C1+ C1- : ',gambrn(ismn+11,1,ih)
*      print*,'H',IH,' -> C1+ C2- : ',gambrn(ismn+12,1,ih)
*      print*,'H',IH,' -> C2+ C1- : ',gambrn(ismn+13,1,ih)
*      print*,'H',IH,' -> C2+ C2- : ',gambrn(ismn+14,1,ih)
*
*---> H_IH -> ST1* ST1     [IM=ISMN+15]
      IM   = ISMN+15
      ISS  = 71
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(1),STMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST1* ST2     [IM=ISMN+16]
      IM   = ISMN+16
      ISS  = 73
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(1),STMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST1     [IM=ISMN+17]
      IM   = ISMN+17
      ISS  = 72
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(2),STMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST2     [IM=ISMN+18]
      IM   = ISMN+18
      ISS  = 74
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(2),STMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SB1* SB1     [IM=ISMN+19]
      IM   = ISMN+19
      ISS  = 75
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(1),SBMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> SB1* SB2     [IM=ISMN+20]
      IM   = ISMN+20
      ISS  = 77
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(1),SBMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SB2* SB1     [IM=ISMN+21]
      IM   = ISMN+21
      ISS  = 76
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(2),SBMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST2     [IM=ISMN+22]
      IM   = ISMN+22
      ISS  = 78
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(2),SBMASS_H(2),GAMBRN(IM,1,IH))
*JSL[01/SEP/05] Bug related STAUMASS_H fixed : Thanks to G. Belanger and S. Pukhov
*---> H_IH -> STAU1* STAU1 [IM=ISMN+23]
      IM   = ISMN+23
      ISS  = 79
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(1),STAUMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> STAU1* STAU2 [IM=ISMN+24]
      IM   = ISMN+24
      ISS  = 80
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(1),STAUMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> STAU2* STAU1 [IM=ISMN+25]
      IM   = ISMN+25
      ISS  = 81
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(2),STAUMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> STAU2* STAU2 [IM=ISMN+26]
      IM   = ISMN+26
      ISS  = 82
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(2),STAUMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SNU3* SNU3   [IM=ISMN+27]
      IM   = ISMN+27
      ISS  = 83
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SNU3MASS_H,SNU3MASS_H,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> ST1* ST1 : ',gambrn(ismn+15,1,ih)
*      print*,'H',IH,' -> ST1* ST2 : ',gambrn(ismn+16,1,ih)
*      print*,'H',IH,' -> ST2* ST1 : ',gambrn(ismn+17,1,ih)
*      print*,'H',IH,' -> ST2* ST2 : ',gambrn(ismn+18,1,ih)
*      print*,'H',IH,' -> SB1* SB1 : ',gambrn(ismn+19,1,ih)
*      print*,'H',IH,' -> SB1* SB2 : ',gambrn(ismn+20,1,ih)
*      print*,'H',IH,' -> SB2* SB1 : ',gambrn(ismn+21,1,ih)
*      print*,'H',IH,' -> SB2* SB2 : ',gambrn(ismn+22,1,ih)
*
      ENDDO ! IH
*-----------------------------------------------------------------------
* 
* << BRANCHING RATIOS OF NEUTRAL HIGGS BOSON >>
*
      DO IH=1,3
*
       GAMBRN(ISMN,1,IH)=0.D0
       DO IM=1,ISMN-1
       GAMBRN(ISMN,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       IXX=ISMN+ISUSYN
       GAMBRN(IXX,1,IH)=0.D0
       DO IM=ISMN+1,IXX-1
       GAMBRN(IXX,1,IH)=GAMBRN(IXX,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       GAMBRN(NMNH,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IXX,1,IH)
*
       DO IM=1,ISMN+ISUSYN+1
        GAMBRN(IM,2,IH)=GAMBRN(IM,1,IH)/GAMBRN(ISMN,1,IH)
        GAMBRN(IM,3,IH)=GAMBRN(IM,1,IH)/GAMBRN(NMNH,1,IH)
       ENDDO
*
      ENDDO ! IH
*-----------------------------------------------------------------------
*      IF(IFLAG_H(6).EQ.1) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,1)
*      IF(IFLAG_H(6).EQ.2) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,2)
*      IF(IFLAG_H(6).EQ.3) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,3)
*      IF(IFLAG_H(6).EQ.5) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,5)
*-----------------------------------------------------------------------
*
* << CHARGED HIGGS BOSON DECAYS INTO SM PARTICLES >>
*

*---> running alpha_s and b-quark mass at Charged Higgs Mass scale
*     : MH > MB^pole assumed
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      ASMZ    = ASMZ_H
      ASMT    = ASMT_H
      MH      = MCH
      IF(MH.LE.MTPOLE_H) THEN                               ! AS(MH)
       ASMH   = ASMZ/(1.D0+B5*ASMZ*DLOG(MH**2/MZ_H**2))
      ELSE
       ASMH   = ASMT/(1.D0+B6*ASMT*DLOG(MH**2/MTPOLE_H**2))
      ENDIF
      IF(MH.LE.MTPOLE_H) THEN                               ! MQ(MH)
       MTMH   = MTMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MBMH   = MBMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MCMH   = MCMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MSMH   = MSMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MUMH   = MUMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MDMH   = MDMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
      ELSE
       MTMH   = MTMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MBMH   = MBMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MCMH   = MCMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MSMH   = MSMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MUMH   = MUMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MDMH   = MDMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
      ENDIF
*      print*,'MCH,AS(MCH),=',MCH,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH
*
*JSL 10/Jun/2009: Including threshold corrections
*       R_123Q=SSPARA_H(22)
*       R_123D=SSPARA_H(24)
*       MQ12  =R_123Q*MQ3_H
*       MD12  =R_123D*MD3_H
*       Q12SQ =DMAX1(MQ12**2,MD12**2)
*       AS_M12=ASMT/(1.D0+B6*ASMT*DLOG(Q12SQ/MTPOLE_H**2))
*       CKS_H=2.D0*AS_M12/3.D0/PI*DCONJG(MU_H*M3_H)
*     .           *F_I(MD12**2,MQ12**2,CDABS(M3_H)**2)
       CKS_H=CAUX_H(11)  ! from FILLCOUPL2
       CKD_H=CKS_H
*       print*,'>>> FILLGAMBR2: cks_h',cks_h
*---> CH+ -> e+ nu     [IM= 1]
      IM   = 1
      IFF  = 1
      MJ   = ME_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> mu+ nu     [IM= 2]
      IM   = 2
      IFF  = 4
      MJ   = MMU_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> tau+ nu     [IM= 3]
      IM   = 3
      IFF  = 7
      MJ   = MTAU_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> u d         [IM= 4] 
      IM   = 4
      IFF  = 10
      GF_UD=DCMPLX(-GW_H*MUMH/DSQRT(2.D0)/MW_H,0.D0)
*      GS_UD=DCMPLX((1.D0/TB_H+MDMH/MUMH*TB_H)/2.D0,0.D0)
*      GP_UD=DCMPLX(0.D0,(1.D0/TB_H-MDMH/MUMH*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections
      GS_UD=(1.D0/TB_H+TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMH/MUMH)/2.D0
      GP_UD=(1.D0/TB_H-TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMH/MUMH)*XI/2.D0
*       print*,'FILLGAMBR2',gs_ud,chc_h(11)
*       print*,'FILLGAMBR2',gp_ud,chc_h(12)
      MJ   = MUMH
      MK   = MDMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,GF_UD,GS_UD,GP_UD,MCH,MJ,MK,GAMBRC(IM,1))
*K-factor
      GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
*---> CH+ -> c s         [IM= 5]
      IM   = 5
      IFF  = 13
      GF_CS=DCMPLX(-GW_H*MCMH/DSQRT(2.D0)/MW_H,0.D0)
*      GS_CS=DCMPLX((1.D0/TB_H+MSMH/MCMH*TB_H)/2.D0,0.D0)
*      GP_CS=DCMPLX(0.D0,(1.D0/TB_H-MSMH/MCMH*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections
      GS_CS=(1.D0/TB_H+TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMH/MCMH)/2.D0
      GP_CS=(1.D0/TB_H-TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMH/MCMH)*XI/2.D0
*       print*,'FILLGAMBR2',gs_cs,chc_h(14)
*       print*,'FILLGAMBR2',gp_cs,chc_h(15)
*       print*,'FILLGAMBR2',msmh,mcmh,msmh/mcmh
      MJ   = MCMH
      MK   = MSMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,GF_CS,GS_CS,GP_CS,MCH,MJ,MK,GAMBRC(IM,1))
*K-factor
      GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
*---> CH+ -> t b         [IM= 6]
      CALL RADNHTB(NFLAG,IFLAG_H,SBMASS_H,STMASS_H
     .            ,HB_H,HT_H,CKB_H,CKT_H)
*       print*,'FILLGAMBR2:CKB_H',ckb_h
      CALL RADCHTB(NFLAG,IFLAG_H,SBMIX_H,SBMASS_H,STMIX_H,STMASS_H
     .            ,RB_H,RT_H,CKBB_H,CKBT_H)
*      print*,'>> FILLGAMBR'
*      print*,DSQRT(2.D0)*MBMT_H/V_H/CB_H,DSQRT(2.D0)*MTMT_H/V_H/SB_H
*      print*,abs(hb_h),abs(ht_h)
*      print*,hb_h,ht_h
*      print*,ckb_h,ckt_h
*      print*,rb_h,rt_h
*      print*,ckbb_h,ckbt_h
      IM   = 6
      IFF  = 16
      XI   =DCMPLX(0.D0,1.D0)
      GF_TB=DCMPLX(-GW_H*MTMH/DSQRT(2.D0)/MW_H,0.D0)
      GS_TB=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMH/MTMH)/2.D0
      GP_TB=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMH/MTMH)*XI/2.D0
*      print*,mtmh,mtmt_h,mtmh/mtmt_h
*      print*,mbmh,mbmt_h,mbmh/mbmt_h
*      print*,chc_h(16),gf_tb,chc_h(16)*mtmh/mtmt_h
*      print*,chc_h(17),gs_tb
*      print*,chc_h(18),gp_tb
      MJ   = MTMH
      MK   = MBMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
*JSLEE 31/MAR/2005 : Three body decay has been included CH -> t* b -> W b b
*JSLEE 10/JUN/2009 : For LO QCD corrections, see, for example, hep-ph/9910339
*      print*,'delta_0',
*     . (1.D0-MW_H**2/MTPOLE_H**2)**2*(1.D0+2.D0*MW_H**2/MTPOLE_H**2) 
      GTWB=GW_H**2*MTPOLE_H**3/64.D0/PI/MW_H**2
     .    *( (1.D0-MW_H**2/MTPOLE_H**2)**2
     .      *(1.D0+2.D0*MW_H**2/MTPOLE_H**2) 
     .      -2.20D0*ASMT_H/PI )
      IF(MTPOLE_H.GT.MCH+MBMH) THEN
       GTHB=CDABS(GF_TB)**2*MTPOLE_H/16.D0/PI
     .    *(CDABS(GS_TB)**2+CDABS(GP_TB)**2)
     .    *(1.D0-MCH**2/MTPOLE_H**2)**2
      ELSE
       GTHB=0.D0
      ENDIF
*JSLEE 10/JUN/2009 : Store top-quark decay widths
       RAUX_H(50)=GTWB
       RAUX_H(51)=GTHB
       RAUX_H(52)=GTWB/(GTWB+GTHB)
       RAUX_H(53)=GTHB/(GTWB+GTHB)
*Top quark decay width
*      GTOP=GTWB+GTHB 
*H^\pm-b loop does not contrubute to the absorptive part of top self-energy
      GTOP=GTWB
*      print*,gtwb,gthb,gtop
      NSTEP=500
*
      IF(MCH.LE.MW_H+2.D0*MK) THEN
       GAMBRC(IM,1)=0.D0
      ELSEIF(MCH.LE.MJ+MK-2.D0) THEN
       CALL CHWBB(NSTEP,GW_H,GF_TB,GS_TB,GP_TB
     .           ,MCH,TB_H,MJ,MW_H,MK,GTOP,GCHWBB)
       GAMBRC(IM,1)=GCHWBB
      ELSEIF(MCH.LE.MJ+MK+2.D0) THEN
*
        MCH1=MJ+MK-3.D0
        ASMH1=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH1**2/MZ_H**2))
        MT1  = MTMT_H*(ASMH1/ASMT)**(1.D0/B5/PI)
        MB1  = MBMT_H*(ASMH1/ASMT)**(1.D0/B5/PI)
        GF1=DCMPLX(-GW_H*MT1/DSQRT(2.D0)/MW_H,0.D0)
        GS1=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB1/MT1)/2.D0
        GP1=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB1/MT1)*XI/2.D0
       CALL CHWBB(NSTEP,GW_H,GF1,GS1,GP1
     .           ,MCH1,TB_H,MT1,MW_H,MB1,GTOP,GCHWBB1)
*
        MCH2=MJ+MK-2.D0
        ASMH_2=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH2**2/MZ_H**2))
        MT2  = MTMT_H*(ASMH_2/ASMT)**(1.D0/B5/PI)
        MB2  = MBMT_H*(ASMH_2/ASMT)**(1.D0/B5/PI)
        GF2=DCMPLX(-GW_H*MT2/DSQRT(2.D0)/MW_H,0.D0)
        GS2=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB2/MT2)/2.D0
        GP2=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB2/MT2)*XI/2.D0
       CALL CHWBB(NSTEP,GW_H,GF2,GS2,GP2
     .           ,MCH2,TB_H,MT2,MW_H,MB2,GTOP,GCHWBB2)
*
        MCH3=MJ+MK+2.D0
        ASMH3=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH3**2/MZ_H**2))
        MT3  = MTMT_H*(ASMH3/ASMT)**(1.D0/B5/PI)
        MB3  = MBMT_H*(ASMH3/ASMT)**(1.D0/B5/PI)
        GF3=DCMPLX(-GW_H*MT3/DSQRT(2.D0)/MW_H,0.D0)
        GS3=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB3/MT3)/2.D0
        GP3=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB3/MT3)*XI/2.D0
       CALL HFF(CF,SYMF,DJK,GF3,GS3,GP3,MCH3,MT3,MB3,GCHWBB3)
        GCHWBB3=(1.D0+5.67D0*ASMH3/PI)*GCHWBB3
*
        MCH4=MJ+MK+3.D0
        ASMH4=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH4**2/MZ_H**2))
        MT4  = MTMT_H*(ASMH4/ASMT)**(1.D0/B5/PI)
        MB4  = MBMT_H*(ASMH4/ASMT)**(1.D0/B5/PI)
        GF4=DCMPLX(-GW_H*MT4/DSQRT(2.D0)/MW_H,0.D0)
        GS4=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB4/MT4)/2.D0
        GP4=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB4/MT4)*XI/2.D0
       CALL HFF(CF,SYMF,DJK,GF4,GS4,GP4,MCH4,MT4,MB4,GCHWBB4)
        GCHWBB4=(1.D0+5.67D0*ASMH4/PI)*GCHWBB4
*
*       print*,mch,mj,mk,gtop
*       print*,mj+mk-3.d0,gchwbb1
*       print*,mj+mk-2.d0,gchwbb2
*       print*,mj+mk+2.d0,gchwbb3
*       print*,mj+mk+3.d0,gchwbb4
*Linear extrapolation
       DMCH=MCH-(MJ+MK-2.D0)
       GAMBRC(IM,1)=GCHWBB2+(GCHWBB3-GCHWBB2)*DMCH/4.D0
      ELSE
       CALL HFF(CF,SYMF,DJK,GF_TB,GS_TB,GP_TB,MCH,MJ,MK,GAMBRC(IM,1))
       GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
      ENDIF
*      print*,mch,mj,mk,gtop
*      print*,mch,gambrc(im,1)
*      MJ   = MTMT_H
*      MK   = MBMT_H
*      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
*     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*      print*,gambrc(im,1)
*
*      print*,'CH+ -> e   nu  : ',gambrc(1,1)
*      print*,'CH+ -> mu  nu  : ',gambrc(2,1)
*      print*,'CH+ -> tau nu  : ',gambrc(3,1)
*      print*,'CH+ -> u   d   : ',gambrc(4,1)
*      print*,'CH+ -> c   s   : ',gambrc(5,1)
*      print*,'CH+ -> t   b   : ',gambrc(6,1)
*---> CH+ -> H1 W        [IM=7]
      IM  = 7
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,1),MCH,HMASS_H(1),MW_H,GAMBRC(IM,1))
* CH+ -> H1 W*
      IF( (MCH.GT.HMASS_H(1)) .AND. 
     .    (MCH.LT.(HMASS_H(1)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(1)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,1))**2 
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*---> CH+ -> H2 W        [IM=8]
      IM  = 8
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,2),MCH,HMASS_H(2),MW_H,GAMBRC(IM,1))
* CH+ -> H2 W*
      IF( (MCH.GT.HMASS_H(2)) .AND. 
     .    (MCH.LT.(HMASS_H(2)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(2)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,2))**2 
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*JSL: Added on Jun.06.2008
*---> CH+ -> H3 W        [IM=9]
      IM  = 9
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,3),MCH,HMASS_H(3),MW_H,GAMBRC(IM,1))
* CH+ -> H3 W*
      IF( (MCH.GT.HMASS_H(3)) .AND.
     .    (MCH.LT.(HMASS_H(3)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(3)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,3))**2
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*      print*,'CH+ -> H1   W   : ',gambrc(7,1)
*      print*,'CH+ -> H2   W   : ',gambrc(8,1)
*      print*,'CH+ -> H3   W   : ',gambrc(9,1)
*-----------------------------------------------------------------------
*
* << CHARGED HIGGS BOSON DECAYS INTO SUSY PARTICLES >>
 
*---> CH+ -> N1 C1+      [IM=ISMC+1]
      IM   = ISMC+1
      IFF  = 19
      MJ   = MN_H(1)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N2 C1+      [IM=ISMC+2]
      IM   = ISMC+2
      IFF  = 25
      MJ   = MN_H(2)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N3 C1+      [IM=ISMC+3]
      IM   = ISMC+3
      IFF  = 31
      MJ   = MN_H(3)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N4 C1+      [IM=ISMC+4]
      IM   = ISMC+4
      IFF  = 37
      MJ   = MN_H(4)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N1 C2+      [IM=ISMC+5]
      IM   = ISMC+5
      IFF  = 22
      MJ   = MN_H(1)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N2 C2+      [IM=ISMC+6]
      IM   = ISMC+6
      IFF  = 28
      MJ   = MN_H(2)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N3 C2+      [IM=ISMC+7]
      IM   = ISMC+7
      IFF  = 34
      MJ   = MN_H(3)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N4 C2+      [IM=ISMC+8]
      IM   = ISMC+8
      IFF  = 40
      MJ   = MN_H(4)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> stop1 sbottom1*   [IM=ISMC+9]
      IM   = ISMC+9
      ISS  = 43
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(1),SBMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> stop1 sbottom2*   [IM=ISMC+10]
      IM   = ISMC+10
      ISS  = 44
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(1),SBMASS_H(2),GAMBRC(IM,1))
*---> CH+ -> stop2 sbottom1*   [IM=ISMC+11]
      IM   = ISMC+11
      ISS  = 45
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(2),SBMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> stop2 sbottom2*   [IM=ISMC+12]
      IM   = ISMC+12
      ISS  = 46
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(2),SBMASS_H(2),GAMBRC(IM,1))
*---> CH+ -> snu3 stau1*       [IM=ISMC+13]
      IM   = ISMC+13
      ISS  = 47
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,SNU3MASS_H,STAUMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> snu3 stau2*       [IM=ISMC+14]
      IM   = ISMC+14
      ISS  = 48
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,SNU3MASS_H,STAUMASS_H(2),GAMBRC(IM,1))

*      print*,'CH+ -> N1   C1+   : ',gambrc(ismc+1,1)
*      print*,'CH+ -> N2   C1+   : ',gambrc(ismc+2,1)
*      print*,'CH+ -> N3   C1+   : ',gambrc(ismc+3,1)
*      print*,'CH+ -> N4   C1+   : ',gambrc(ismc+4,1)
*      print*,'CH+ -> N1   C2+   : ',gambrc(ismc+5,1)
*      print*,'CH+ -> N2   C2+   : ',gambrc(ismc+6,1)
*      print*,'CH+ -> N3   C2+   : ',gambrc(ismc+7,1)
*      print*,'CH+ -> N4   C2+   : ',gambrc(ismc+8,1)
*
*-----------------------------------------------------------------------
* 
* << BRANCHING RATIOS OF CHARGED HIGGS BOSON >>
*
       GAMBRC(ISMC,1)=0.D0
       DO IM=1,ISMC-1
       GAMBRC(ISMC,1)=GAMBRC(ISMC,1)+GAMBRC(IM,1)
       ENDDO
*
       IXX=ISMC+ISUSYC
       GAMBRC(IXX,1)=0.D0
       DO IM=ISMC+1,IXX-1
       GAMBRC(IXX,1)=GAMBRC(IXX,1)+GAMBRC(IM,1)
       ENDDO
*
       GAMBRC(NMCH,1)=GAMBRC(ISMC,1)+GAMBRC(IXX,1)
*
       DO IM=1,ISMC+ISUSYC+1
        GAMBRC(IM,2)=GAMBRC(IM,1)/GAMBRC(ISMC,1)
        GAMBRC(IM,3)=GAMBRC(IM,1)/GAMBRC(NMCH,1)
       ENDDO
*
*-----------------------------------------------------------------------
*      IF(IFLAG_H(6).EQ.4.OR.IFLAG_H(6).EQ.5) 
*     . CALL DUMP_CHDCY(ISMC,ISUSYC,NMCH,GAMBRC)
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE DUMP_CHDCY(ISMC,ISUSYC,NMCH,GAMBRC)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 ISMC,ISUSYC
      REAL*8 GAMBRC(NMCH,3)   
*
      print*,'---------------------------------------------------------'
*      print*,'Charged Higgs Boson Decays with '
*     .,'ISMC = ',ISMC,' : ISUSYC = ',ISUSYC
      print*,'               Charged Higgs Boson Decays'
      print*,'---------------------------------------------------------'
      print*,'DECAY MODE      [IM]   WIDTH[GeV]  BR[SM]      BR[TOTAL]'
      print*,'---------------------------------------------------------'
      WRITE(*,1)GAMBRC(1,1),GAMBRC(1,2),GAMBRC(1,3)
      WRITE(*,2)GAMBRC(2,1),GAMBRC(2,2),GAMBRC(2,3)
      WRITE(*,3)GAMBRC(3,1),GAMBRC(3,2),GAMBRC(3,3)
      WRITE(*,4)GAMBRC(4,1),GAMBRC(4,2),GAMBRC(4,3)
      WRITE(*,5)GAMBRC(5,1),GAMBRC(5,2),GAMBRC(5,3)
      WRITE(*,6)GAMBRC(6,1),GAMBRC(6,2),GAMBRC(6,3)
      WRITE(*,7)GAMBRC(7,1),GAMBRC(7,2),GAMBRC(7,3)
      WRITE(*,8)GAMBRC(8,1),GAMBRC(8,2),GAMBRC(8,3)
      WRITE(*,9)GAMBRC(9,1),GAMBRC(9,2),GAMBRC(9,3)
      WRITE(*,25)GAMBRC(25,1),GAMBRC(25,2),GAMBRC(25,3)
      WRITE(*,26)GAMBRC(26,1),GAMBRC(26,2),GAMBRC(26,3)
      WRITE(*,27)GAMBRC(27,1),GAMBRC(27,2),GAMBRC(27,3)
      WRITE(*,28)GAMBRC(28,1),GAMBRC(28,2),GAMBRC(28,3)
      WRITE(*,29)GAMBRC(29,1),GAMBRC(29,2),GAMBRC(29,3)
      WRITE(*,30)GAMBRC(30,1),GAMBRC(30,2),GAMBRC(30,3)
      WRITE(*,31)GAMBRC(31,1),GAMBRC(31,2),GAMBRC(31,3)
      WRITE(*,32)GAMBRC(32,1),GAMBRC(32,2),GAMBRC(32,3)
      WRITE(*,33)GAMBRC(33,1),GAMBRC(33,2),GAMBRC(33,3)
      WRITE(*,34)GAMBRC(34,1),GAMBRC(34,2),GAMBRC(34,3)
      WRITE(*,35)GAMBRC(35,1),GAMBRC(35,2),GAMBRC(35,3)
      WRITE(*,36)GAMBRC(36,1),GAMBRC(36,2),GAMBRC(36,3)
      WRITE(*,37)GAMBRC(37,1),GAMBRC(37,2),GAMBRC(37,3)
      WRITE(*,38)GAMBRC(38,1),GAMBRC(38,2),GAMBRC(38,3)
      WRITE(*,39)GAMBRC(39,1),GAMBRC(39,2),GAMBRC(39,3)
      WRITE(*,50)GAMBRC(50,1),GAMBRC(50,2),GAMBRC(50,3)
      WRITE(*,51)GAMBRC(51,1),GAMBRC(51,2),GAMBRC(51,3)
      print*,' '
      print*,'* Note : WIDTH=GAMBRC(IM,1), BR[SM]   =GAMBRC(IM,2) '
      print*,'                      and    BR[TOTAL]=GAMBRC(IM,3) '
      print*,'---------------------------------------------------------'
*
  1   FORMAT(1X,'CH+ -> e+   nu  [ 1]:',3(2X,E10.4))
  2   FORMAT(1X,'CH+ -> mu+  nu  [ 2]:',3(2X,E10.4))
  3   FORMAT(1X,'CH+ -> tau+ nu  [ 3]:',3(2X,E10.4))
  4   FORMAT(1X,'CH+ -> u    d   [ 4]:',3(2X,E10.4))
  5   FORMAT(1X,'CH+ -> c    s   [ 5]:',3(2X,E10.4))
  6   FORMAT(1X,'CH+ -> t    b   [ 6]:',3(2X,E10.4))
  7   FORMAT(1X,'CH+ -> H1   W   [ 7]:',3(2X,E10.4))
  8   FORMAT(1X,'CH+ -> H2   W   [ 8]:',3(2X,E10.4))
  9   FORMAT(1X,'CH+ -> H3   W   [ 9]:',3(2X,E10.4))
 25   FORMAT(1X,'CH+ TOTAL(SM)   [25]:',3(2X,E10.4))
 26   FORMAT(1X,'CH+ -> N1   C1+ [26]:',3(2X,E10.4))
 27   FORMAT(1X,'CH+ -> N2   C1+ [27]:',3(2X,E10.4))
 28   FORMAT(1X,'CH+ -> N3   C1+ [28]:',3(2X,E10.4))
 29   FORMAT(1X,'CH+ -> N4   C1+ [29]:',3(2X,E10.4))
 30   FORMAT(1X,'CH+ -> N1   C2+ [30]:',3(2X,E10.4))
 31   FORMAT(1X,'CH+ -> N2   C2+ [31]:',3(2X,E10.4))
 32   FORMAT(1X,'CH+ -> N3   C2+ [32]:',3(2X,E10.4))
 33   FORMAT(1X,'CH+ -> N4   C2+ [33]:',3(2X,E10.4))
 34   FORMAT(1X,'CH+ -> ST1 SB1* [34]:',3(2X,E10.4))
 35   FORMAT(1X,'CH+ -> ST1 SB2* [35]:',3(2X,E10.4))
 36   FORMAT(1X,'CH+ -> ST2 SB1* [36]:',3(2X,E10.4))
 37   FORMAT(1X,'CH+ -> ST2 SB2* [37]:',3(2X,E10.4))
 38   FORMAT(1X,'CH+ ->SNU3 STA1*[38]:',3(2X,E10.4))
 39   FORMAT(1X,'CH+ ->SNU3 STA2*[39]:',3(2X,E10.4))
 50   FORMAT(1X,'CH+ TOTAL(SUSY) [50]:',3(2X,E10.4))
 51   FORMAT(1X,'CH+ TOTAL       [51]:',3(2X,E10.4))
*
      RETURN
      END

      SUBROUTINE DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,IPRI)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 ISMN,ISUSYN
      REAL*8 GAMBRN(NMNH,3,3)   
*
      IF(IPRI.LE.3) IH=IPRI
      IC = 0
 999  CONTINUE
      IC = IC+1
      IF(IPRI.EQ.5) IH=IC
*
      print*,'---------------------------------------------------------'
*      print*,'Neutral Higgs Boson Decays with '
*     .,'ISMN = ',ISMN,' : ISUSYN = ',ISUSYN
      print*,'               Neutral Higgs Boson Decays'
      print*,'---------------------------------------------------------'
      print*,'DECAY MODE    [ IM]   WIDTH[GeV]  BR[SM]      BR[TOTAL]'
      print*,'---------------------------------------------------------'
       WRITE(*,1)IH,GAMBRN(1,1,ih),GAMBRN(1,2,ih),GAMBRN(1,3,ih)
       WRITE(*,2)IH,GAMBRN(2,1,ih),GAMBRN(2,2,ih),GAMBRN(2,3,ih)
       WRITE(*,3)IH,GAMBRN(3,1,ih),GAMBRN(3,2,ih),GAMBRN(3,3,ih)
       WRITE(*,4)IH,GAMBRN(4,1,ih),GAMBRN(4,2,ih),GAMBRN(4,3,ih)
       WRITE(*,5)IH,GAMBRN(5,1,ih),GAMBRN(5,2,ih),GAMBRN(5,3,ih)
       WRITE(*,6)IH,GAMBRN(6,1,ih),GAMBRN(6,2,ih),GAMBRN(6,3,ih)
       WRITE(*,7)IH,GAMBRN(7,1,ih),GAMBRN(7,2,ih),GAMBRN(7,3,ih)
       WRITE(*,8)IH,GAMBRN(8,1,ih),GAMBRN(8,2,ih),GAMBRN(8,3,ih)
       WRITE(*,9)IH,GAMBRN(9,1,ih),GAMBRN(9,2,ih),GAMBRN(9,3,ih)
       WRITE(*,10)IH,GAMBRN(10,1,ih),GAMBRN(10,2,ih),GAMBRN(10,3,ih)
       WRITE(*,11)IH,GAMBRN(11,1,ih),GAMBRN(11,2,ih),GAMBRN(11,3,ih)
       WRITE(*,12)IH,GAMBRN(12,1,ih),GAMBRN(12,2,ih),GAMBRN(12,3,ih)
       WRITE(*,13)IH,GAMBRN(13,1,ih),GAMBRN(13,2,ih),GAMBRN(13,3,ih)
       WRITE(*,14)IH,GAMBRN(14,1,ih),GAMBRN(14,2,ih),GAMBRN(14,3,ih)
       WRITE(*,15)IH,GAMBRN(15,1,ih),GAMBRN(15,2,ih),GAMBRN(15,3,ih)
       WRITE(*,16)IH,GAMBRN(16,1,ih),GAMBRN(16,2,ih),GAMBRN(16,3,ih)
       WRITE(*,17)IH,GAMBRN(17,1,ih),GAMBRN(17,2,ih),GAMBRN(17,3,ih)
       WRITE(*,18)IH,GAMBRN(18,1,ih),GAMBRN(18,2,ih),GAMBRN(18,3,ih)
       WRITE(*,19)IH,GAMBRN(19,1,ih),GAMBRN(19,2,ih),GAMBRN(19,3,ih)
       WRITE(*,50)IH,GAMBRN(50,1,ih),GAMBRN(50,2,ih),GAMBRN(50,3,ih)
       WRITE(*,51)IH,GAMBRN(51,1,ih),GAMBRN(51,2,ih),GAMBRN(51,3,ih)
       WRITE(*,52)IH,GAMBRN(52,1,ih),GAMBRN(52,2,ih),GAMBRN(52,3,ih)
       WRITE(*,53)IH,GAMBRN(53,1,ih),GAMBRN(53,2,ih),GAMBRN(53,3,ih)
       WRITE(*,54)IH,GAMBRN(54,1,ih),GAMBRN(54,2,ih),GAMBRN(54,3,ih)
       WRITE(*,55)IH,GAMBRN(55,1,ih),GAMBRN(55,2,ih),GAMBRN(55,3,ih)
       WRITE(*,56)IH,GAMBRN(56,1,ih),GAMBRN(56,2,ih),GAMBRN(56,3,ih)
       WRITE(*,57)IH,GAMBRN(57,1,ih),GAMBRN(57,2,ih),GAMBRN(57,3,ih)
       WRITE(*,58)IH,GAMBRN(58,1,ih),GAMBRN(58,2,ih),GAMBRN(58,3,ih)
       WRITE(*,59)IH,GAMBRN(59,1,ih),GAMBRN(59,2,ih),GAMBRN(59,3,ih)
       WRITE(*,60)IH,GAMBRN(60,1,ih),GAMBRN(60,2,ih),GAMBRN(60,3,ih)
       WRITE(*,61)IH,GAMBRN(61,1,ih),GAMBRN(61,2,ih),GAMBRN(61,3,ih)
       WRITE(*,62)IH,GAMBRN(62,1,ih),GAMBRN(62,2,ih),GAMBRN(62,3,ih)
       WRITE(*,63)IH,GAMBRN(63,1,ih),GAMBRN(63,2,ih),GAMBRN(63,3,ih)
       WRITE(*,64)IH,GAMBRN(64,1,ih),GAMBRN(64,2,ih),GAMBRN(64,3,ih)
       WRITE(*,65)IH,GAMBRN(65,1,ih),GAMBRN(65,2,ih),GAMBRN(65,3,ih)
       WRITE(*,66)IH,GAMBRN(66,1,ih),GAMBRN(66,2,ih),GAMBRN(66,3,ih)
       WRITE(*,67)IH,GAMBRN(67,1,ih),GAMBRN(67,2,ih),GAMBRN(67,3,ih)
       WRITE(*,68)IH,GAMBRN(68,1,ih),GAMBRN(68,2,ih),GAMBRN(68,3,ih)
       WRITE(*,69)IH,GAMBRN(69,1,ih),GAMBRN(69,2,ih),GAMBRN(69,3,ih)
       WRITE(*,70)IH,GAMBRN(70,1,ih),GAMBRN(70,2,ih),GAMBRN(70,3,ih)
       WRITE(*,71)IH,GAMBRN(71,1,ih),GAMBRN(71,2,ih),GAMBRN(71,3,ih)
       WRITE(*,72)IH,GAMBRN(72,1,ih),GAMBRN(72,2,ih),GAMBRN(72,3,ih)
       WRITE(*,73)IH,GAMBRN(73,1,ih),GAMBRN(73,2,ih),GAMBRN(73,3,ih)
       WRITE(*,74)IH,GAMBRN(74,1,ih),GAMBRN(74,2,ih),GAMBRN(74,3,ih)
       WRITE(*,75)IH,GAMBRN(75,1,ih),GAMBRN(75,2,ih),GAMBRN(75,3,ih)
       WRITE(*,76)IH,GAMBRN(76,1,ih),GAMBRN(76,2,ih),GAMBRN(76,3,ih)
       WRITE(*,77)IH,GAMBRN(77,1,ih),GAMBRN(77,2,ih),GAMBRN(77,3,ih)
       WRITE(*,100)IH,GAMBRN(100,1,ih),GAMBRN(100,2,ih),GAMBRN(100,3,ih)
       WRITE(*,101)IH,GAMBRN(101,1,ih),GAMBRN(101,2,ih),GAMBRN(101,3,ih)
      print*,' '
      IF(IH.EQ.1) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,1), BR[SM]   =GAMBRN(IM,2,1) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,1) '
      ELSEIF(IH.EQ.2) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,2), BR[SM]   =GAMBRN(IM,2,2) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,2) '
      ELSEIF(IH.EQ.3) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,3), BR[SM]   =GAMBRN(IM,2,3) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,3) '
      ENDIF
      IF(IPRI.EQ.5) THEN
       IF(IC.EQ.1.OR.IC.EQ.2) GOTO 999
      ENDIF
      print*,'---------------------------------------------------------'
*
  1   FORMAT(1X,'H',I1,' -> e    e  [  1]:',3(2X,E10.4))
  2   FORMAT(1X,'H',I1,' -> mu   mu [  2]:',3(2X,E10.4))
  3   FORMAT(1X,'H',I1,' -> tau  tau[  3]:',3(2X,E10.4))
  4   FORMAT(1X,'H',I1,' -> d    d  [  4]:',3(2X,E10.4))
  5   FORMAT(1X,'H',I1,' -> s    s  [  5]:',3(2X,E10.4))
  6   FORMAT(1X,'H',I1,' -> b    b  [  6]:',3(2X,E10.4))
  7   FORMAT(1X,'H',I1,' -> u    u  [  7]:',3(2X,E10.4))
  8   FORMAT(1X,'H',I1,' -> c    c  [  8]:',3(2X,E10.4))
  9   FORMAT(1X,'H',I1,' -> t    t  [  9]:',3(2X,E10.4))
 10   FORMAT(1X,'H',I1,' -> W    W  [ 10]:',3(2X,E10.4))
 11   FORMAT(1X,'H',I1,' -> Z    Z  [ 11]:',3(2X,E10.4))
 12   FORMAT(1X,'H',I1,' -> H1   Z  [ 12]:',3(2X,E10.4))
 13   FORMAT(1X,'H',I1,' -> H2   Z  [ 13]:',3(2X,E10.4))
 14   FORMAT(1X,'H',I1,' -> H1   H1 [ 14]:',3(2X,E10.4))
 15   FORMAT(1X,'H',I1,' -> H1   H2 [ 15]:',3(2X,E10.4))
 16   FORMAT(1X,'H',I1,' -> H2   H2 [ 16]:',3(2X,E10.4))
 17   FORMAT(1X,'H',I1,' -> ph   ph [ 17]:',3(2X,E10.4))
 18   FORMAT(1X,'H',I1,' -> gl   gl [ 18]:',3(2X,E10.4))
 19   FORMAT(1X,'H',I1,' -> Z    ph [ 19]:',3(2X,E10.4))
 50   FORMAT(1X,'H',I1,' TOTAL(SM)  [ 50]:',3(2X,E10.4))
 51   FORMAT(1X,'H',I1,' -> N1   N1 [ 51]:',3(2X,E10.4))
 52   FORMAT(1X,'H',I1,' -> N1   N2 [ 52]:',3(2X,E10.4))
 53   FORMAT(1X,'H',I1,' -> N1   N3 [ 53]:',3(2X,E10.4))
 54   FORMAT(1X,'H',I1,' -> N1   N4 [ 54]:',3(2X,E10.4))
 55   FORMAT(1X,'H',I1,' -> N2   N2 [ 55]:',3(2X,E10.4))
 56   FORMAT(1X,'H',I1,' -> N2   N3 [ 56]:',3(2X,E10.4))
 57   FORMAT(1X,'H',I1,' -> N2   N4 [ 57]:',3(2X,E10.4))
 58   FORMAT(1X,'H',I1,' -> N3   N3 [ 58]:',3(2X,E10.4))
 59   FORMAT(1X,'H',I1,' -> N3   N4 [ 59]:',3(2X,E10.4))
 60   FORMAT(1X,'H',I1,' -> N4   N4 [ 60]:',3(2X,E10.4))
 61   FORMAT(1X,'H',I1,' -> C1+  C1-[ 61]:',3(2X,E10.4))
 62   FORMAT(1X,'H',I1,' -> C1+  C2-[ 62]:',3(2X,E10.4))
 63   FORMAT(1X,'H',I1,' -> C2+  C1-[ 63]:',3(2X,E10.4))
 64   FORMAT(1X,'H',I1,' -> C2+  C2-[ 64]:',3(2X,E10.4))
 65   FORMAT(1X,'H',I1,' -> ST1* ST1[ 65]:',3(2X,E10.4))
 66   FORMAT(1X,'H',I1,' -> ST1* ST2[ 66]:',3(2X,E10.4))
 67   FORMAT(1X,'H',I1,' -> ST2* ST1[ 67]:',3(2X,E10.4))
 68   FORMAT(1X,'H',I1,' -> ST2* ST2[ 68]:',3(2X,E10.4))
 69   FORMAT(1X,'H',I1,' -> SB1* SB1[ 69]:',3(2X,E10.4))
 70   FORMAT(1X,'H',I1,' -> SB1* SB2[ 70]:',3(2X,E10.4))
 71   FORMAT(1X,'H',I1,' -> SB2* SB1[ 71]:',3(2X,E10.4))
 72   FORMAT(1X,'H',I1,' -> SB2* SB2[ 72]:',3(2X,E10.4))
 73   FORMAT(1X,'H',I1,' ->STA1*STA1[ 73]:',3(2X,E10.4))
 74   FORMAT(1X,'H',I1,' ->STA1*STA2[ 74]:',3(2X,E10.4))
 75   FORMAT(1X,'H',I1,' ->STA2*STA1[ 75]:',3(2X,E10.4))
 76   FORMAT(1X,'H',I1,' ->STA2*STA2[ 76]:',3(2X,E10.4))
 77   FORMAT(1X,'H',I1,' ->SNU3*SNU3[ 77]:',3(2X,E10.4))
100   FORMAT(1X,'H',I1,' TOTAL(SUSY)[100]:',3(2X,E10.4))
101   FORMAT(1X,'H',I1,' TOTAL      [101]:',3(2X,E10.4))
      RETURN
      END

      SUBROUTINE HSS(SF,V,G,MH,MHJ,MHK,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 G
*
      PI=2.D0*DASIN(1.D0)
*
      XKJ=MHJ**2/MH**2
      XKK=MHK**2/MH**2
      XLAM=(1.D0-XKJ-XKK)**2-4.D0*XKJ*XKK
      IF(MH.GT.(MHJ+MHK) .AND. XLAM.GT.0.D0) THEN
       GAM=SF*V**2*CDABS(G)**2/16.D0/PI/MH*DSQRT(XLAM)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HHV(GF,GV,MH,MHJ,MV,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GV
*
      PI=2.D0*DASIN(1.D0)
*
      XKV=MV**2/MH**2
      XKJ=MHJ**2/MH**2
      XLAM=(1.D0-XKJ-XKV)**2-4.D0*XKJ*XKV
      IF(MH.GT.(MHJ+MV) .AND. XLAM.GT.0.D0) THEN
       GAM=GF*MH**3/8.D0/DSQRT(2.D0)/PI*CDABS(GV)**2*XLAM*DSQRT(XLAM)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HVV(GF,DV,GV,MV,MH,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GV
*
      PI=2.D0*DASIN(1.D0)
*
      XKV=MV**2/MH**2
      BETASQ=1.D0-4.D0*XKV
      IF(BETASQ.GE.0.D0) THEN
       BETA=DSQRT(BETASQ)
       GAM=GF*CDABS(GV)**2*MH**3*DV/16.D0/DSQRT(2.D0)/PI*BETA
     .    *(1.D0-4.D0*XKV+12.D0*XKV**2)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HFF(CF,SF,DJK,GF,GS,GP,MH,MJ,MK,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GF,GS,GP
*
      PI=2.D0*DASIN(1.D0)
*
      XKJ=MJ**2/MH**2
      XKK=MK**2/MH**2
      XLAM=(1.D0-XKJ-XKK)**2-4.D0*XKJ*XKK
      IF(XLAM.LE.0.D0 .OR. MH.LE.(MJ+MK)) THEN
       GAM=0.D0
       RETURN
      ELSE
       GAM=SF**2/(1.D0+DJK)*CF*CDABS(GF)**2*MH*DSQRT(XLAM)*
     .    ((1.D0-XKJ-XKK)*(CDABS(GS)**2+CDABS(GP)**2)-
     .     2.D0*DSQRT(XKJ*XKK)*(CDABS(GS)**2-CDABS(GP)**2))/8.D0/PI
*
*       IF(MJ.EQ.MK) THEN
*         BETA=DSQRT(1.D0-4.D0*XKJ)
*         GAMP=SF**2/(1.D0+DJK)*CF*CDABS(GF)**2*MH*
*     .        BETA*(BETA**2*CDABS(GS)**2+CDABS(GP)**2)/8.D0/PI
*         print*,'same ?',gam,gamp
*       ENDIF
*
      ENDIF
*
      RETURN
      END

      SUBROUTINE BODE(FBODE,XIN,XOUT,NSTEP,YINT)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%                                                                      %
c%     BODE makes one-dimensional numerical integration of a            % 
c%     F(x), using Bode's rule:                                         %
c%                                                                      %
c%     Integral_x0^x2 F(x) dx  =   h/598752 * [ 16067 (f_0 + f_10)      %
c%     + 106300 (f_1 + f_9) - 48525 (f_2 + f_8) + 272400 (f_3 +f_7)     %
c%     - 260550 (f_4 + f_6) + 427368 f_5 ]                              %
c%                                                                      %
c%    FBODE defines the integrand  F(x)                                 %
c%    XIN  is the lower limit of the integral                           %
c%    XOUT is the upper limit of the integral                           %
c%    NSTEP determines the total number of steps                        %
c%    YINT contains the result of the integration                       %
c%                                                                      %
c%                                                                      %
c%    NOTES: 1. We always assume that XOUT > XIN                        %
c%           2. The integrand F(x) must be defined as an external       %
C%              function in the main program                            %
c%                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      IMPLICIT REAL*8(A,B,D-H,O-Z), COMPLEX*16(C)
C
CCC
C
      YINT = 0.D0
*      DELTA = (XOUT - XIN)/DREAL(NSTEP)
      DELTA = (XOUT - XIN)/DBLE(NSTEP)
      X10 = XIN
999   CONTINUE
      X0 = X10
       X1 = X10 + DELTA/10.D0
        X2 = X10 + DELTA/5.D0
         X3 = X10 + 3.D0*DELTA/10.D0
          X4 = X10 + 2.D0*DELTA/5.D0
           X5 = X10 + DELTA/2.D0
            X6 = X10 + 3.D0*DELTA/5.D0
             X7 = X10 + 7.D0*DELTA/10.D0
              X8 = X10 + 4.D0*DELTA/5.D0
               X9 = X10 + 9.D0*DELTA/10.D0
                X10 = X10 + DELTA
      IF(X10.GE.XOUT) GOTO 9999
      F0 = FBODE(X0)
       F1 = FBODE(X1)
        F2 = FBODE(X2)
         F3 = FBODE(X3)
          F4 = FBODE(X4)
           F5 = FBODE(X5)
            F6 = FBODE(X6)
             F7 = FBODE(X7)
              F8 = FBODE(X8)
               F9 = FBODE(X9)
                F10 = FBODE(X10)
      YINT = DELTA*( 16067.D0*(F0+F10) + 106300.D0*(F1+F9)
     #- 48525.D0*(F2+F8) + 272400.D0*(F3+F7) - 260550.D0*(F4+F6)
     #+ 427368.D0*F5  ) / 598752.D0  + YINT
      GOTO 999
9999  CONTINUE
      RETURN
      END
C
      REAL*8 FUNCTION FVVS(R)
************************************************************************
*
* Here, we used the integration method:
*
* I    =               \int_{x-}^{x+} dx f(x) 
*      = {G(x+)-G(x-)} \int_{0}^{1}   dR f(x)/g(x)
*
*  where      f(x) = F(x)/[(x-a)^2+b^2] 
*             g(x) = b^2/[(x-a)^2+b^2]
*             G(x) = pi/2*b+b*arctg[(x-a)/b]
*               x  = G^{-1}[G(x-)+R*{G(x+)-G(x-)}]
*        G^{-1}(y) = a+b*tg[(y-pi/2*b)/b]
*
*  Note : a=1 and b=EPSV 
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
*
      PI=2.D0*DASIN(1.D0)
      GXUP=PI/2.D0*EPSV+EPSV*DATAN((XUP-1.D0)/EPSV)
      GXDW=PI/2.D0*EPSV+EPSV*DATAN((XDW-1.D0)/EPSV)
      Y=GXDW+R*(GXUP-GXDW)
      XX=1.D0+EPSV*DTAN((Y-PI/2.D0*EPSV)/EPSV)
      XLAM=(1.D0-OMEGAI-XX)**2-4.D0*OMEGAI*XX
      IF(XLAM.LT.0.D0) XLAM=0.D0
      F_XX=DSQRT(XLAM)*(XLAM+12.D0*XX)
      FVVS=(GXUP-GXDW)*F_XX/EPSV**2
*
      RETURN
      END

      REAL*8 FUNCTION FHVS(R)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
*
      PI=2.D0*DASIN(1.D0)
      GXUP=PI/2.D0*EPSV+EPSV*DATAN((XUP-1.D0)/EPSV)
      GXDW=PI/2.D0*EPSV+EPSV*DATAN((XDW-1.D0)/EPSV)
      Y=GXDW+R*(GXUP-GXDW)
      XX=1.D0+EPSV*DTAN((Y-PI/2.D0*EPSV)/EPSV)
      XLAM=(XX-OMEGAI-OMEGAJ)**2-4.D0*OMEGAI*OMEGAJ
      IF(XLAM.LT.0.D0) XLAM=0.D0
      F_XX=DSQRT(XLAM)**3
      FHVS=(GXUP-GXDW)*F_XX/EPSV**2
*
      RETURN
      END

      SUBROUTINE CHWBB(NX1,GW,GTB,GS,GP,MCH,TANB,MT,MW,MB,GTOP,GCHWBB)
************************************************************************
*
* Calculate the three body decay rate H+ -> t* b^bar -> W^+ b b^bar
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
* For integration
      COMMON /CHWBB_RBODE/ MCH_CHWBB,TANB_CHWBB,MT_CHWBB,MW_CHWBB
     .                    ,MB_CHWBB,GTOP_CHWBB
      COMPLEX*16 GS_CHWBB,GP_CHWBB
      COMMON /CHWBB_CBODE/ GS_CHWBB,GP_CHWBB
      EXTERNAL FX1
*Local
      COMPLEX*16 GTB,GS,GP
      PI=2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
      MCH_CHWBB  = MCH
      TANB_CHWBB = TANB
      MT_CHWBB   = MT
      MW_CHWBB   = MW
      MB_CHWBB   = MB
      GTOP_CHWBB = GTOP
*
      GS_CHWBB   = GS
      GP_CHWBB   = GP
*
      XKW = MW**2/MCH**2
*
      X1D=0.D0
      X1U=1.D0-MW**2/MCH**2
*
      CALL BODE(FX1,X1D,X1U,NX1,RES)
*
      GCHWBB=3.D0*GW**2*CDABS(GTB)**2*MCH/512.D0/PI**3*RES
*
      RETURN
      END

      REAL*8 FUNCTION FX1(X1)
************************************************************************
*
* It returns a function FX1 as function of x1 after x2 integration.
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      COMMON /CHWBB_RBODE/ MCH_CHWBB,TANB_CHWBB,MT_CHWBB,MW_CHWBB
     .                    ,MB_CHWBB,GTOP_CHWBB
      COMPLEX*16 GS_CHWBB,GP_CHWBB
      COMMON /CHWBB_CBODE/ GS_CHWBB,GP_CHWBB
*Local
      COMPLEX*16 GS,GP,GL,GR
      COMPLEX*16 XI
*-----------------------------------------------------------------------
      XI   =DCMPLX(0.D0,1.D0)
*
      MCH  = MCH_CHWBB
      TANB = TANB_CHWBB
      MT   = MT_CHWBB
      MW   = MW_CHWBB
      MB   = MB_CHWBB
      GTOP = GTOP_CHWBB
      GS   = GS_CHWBB
      GP   = GP_CHWBB
* 
      XKT = MT**2/MCH**2
      XKW = MW**2/MCH**2
      XKB = MB**2/MCH**2
      XGT = GTOP**2/MCH**2
      GR0 = MB/MT*TANB
      GL0 = 1.D0/TANB
      GL  = GS-XI*GP
      GR  = GS+XI*GP
*      print*,'GR',gr0,dreal(gr),dimag(gr)
*      print*,'GL',gl0,dreal(gl),dimag(gl)
*
*The x2 intgration has been done by REDUCE:
*
      CALL GET_FN(XKT,XKW,XKB,XGT,X1,F0,F1,F2,F3)
*
      FXL2=-2.D0*XKB*XKT*F0+XKT*((1.D0-X1)*(F0-F1)/XKW+2.D0*X1*F0
     .     +2.D0*F1-3.D0*F0+2.D0*XKW*F0)
      FXR2=-2.D0*XKB**2*F0
     .     +XKB*((2.D0*XKW-2.D0*X1+3.D0)*F0
     .     +(-2.D0*F2-X1*F1+X1*F0+5.D0*F1-3.D0*F0)/XKW)
     .     +(F2+2.D0*X1*F1-2.D0*X1*F0-4.D0*F1+3.D0*F0-2.D0*XKW*F0)
     .     +(F3+X1*F2-3.D0*F2-2.D0*X1*F1+X1*F0+3.D0*F1-F0)/XKW
      FXLR=2.D0*XKB*F0-2.D0*XKW*F0-F1+F0+(F2-2.D0*F1+F0)/XKW
  
      FX1=CDABS(GL)**2*FXL2+CDABS(GR)**2*FXR2
     .   +2.D0*DSQRT(XKB*XKT)*DREAL(GL*DCONJG(GR))*FXLR
*
      RETURN
      END


      SUBROUTINE GET_FN(XKT,XKW,XKB,XGT,X1,F0,F1,F2,F3)
************************************************************************
*
*      /               x2^n
* Fn = | dx2  --------------------------
*      /      [(1-x2-xkt+xkb)^2+xkt*xgt]
*
*     x2^up =1-kw/(1-x1)   ; x2_down = 1-kw-x1
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      X2U=1.-XKW/(1.-X1)
      X2D=1.-XKW-X1
*-----------------------------------------------------------------------
*
* This is an out put of a REDUCE program 
*
*-----------------------------------------------------------------------
* f0
      f0=(sqrt(xkt)*sqrt(xgt)*(-atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))+atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))))/(
     . xgt*xkt)
* f1
      f1=(-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))*xkb+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)
     . /(sqrt(xkt)*sqrt(xgt)))*xkt-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-
     . xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))+2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-2.0*sqrt(xkt
     . )*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt+
     . 2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))+log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d
     . *xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+
     . xkt**2-2.0*xkt+1.0))*xgt*xkt)/(2.0*xgt*xkt)
* f2
      ans2=log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2
     . -2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt*xkt**2+log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-
     . 2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)
     . /(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*
     . xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt+log((x2u**2-
     . 2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0
     . *xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*
     . x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*
     . xgt*xkt-x2d*xgt*xkt+x2u*xgt*xkt
      ans1=sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))*xgt*xkt-sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb**2+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d
     . -xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*xkt-2.0*sqrt(xkt)*
     . sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-
     . sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt
     . )))*xkt**2+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkt-sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))-sqrt(xkt)*sqrt(xgt)*atan((x2u-
     . xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt+sqrt(xkt)*sqrt(xgt
     . )*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2-2.0*
     . sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt
     . )))*xkb*xkt+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb+sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**2-2.0*sqrt(xkt)*sqrt(xgt)
     . *atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt+sqrt(xkt)*
     . sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))+ans2
      f2=ans1/(xgt*xkt)
* f3
      ans4=6.0*log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d
     . *xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+
     . xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt+3.0*log((x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkt**3+3.0
     . *log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-
     . 2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+
     . 2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt*xkt-x2d**2*xgt*xkt-4.0*x2d*xgt*xkb*xkt+4.0*
     . x2d*xgt*xkt**2-4.0*x2d*xgt*xkt+x2u**2*xgt*xkt+4.0*x2u*xgt*xkb*
     . xkt-4.0*x2u*xgt*xkt**2+4.0*x2u*xgt*xkt
      ans3=log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2
     . -2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt**2*xkt**2+6.0*log((x2d**2-2.0*x2d*xkb+2.0*
     . x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*
     . xkt+1.0)/(x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**
     . 2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt**2+6.0*
     . log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0
     . *xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+2.0*
     . x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*
     . xkt+1.0))*xgt*xkt**2+3.0*log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-
     . 2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)
     . /(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*
     . xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb**2*xkt+ans4
      ans2=-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)
     . *sqrt(xgt)))*xgt*xkb*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt**2-6.0*sqrt(xkt)*sqrt
     . (xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt+
     . 2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))*xkb**3-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/
     . (sqrt(xkt)*sqrt(xgt)))*xkb**2*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan
     . ((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2+6.0*sqrt(xkt)
     . *sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*
     . xkt**2-12.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkb*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**3+6.0*sqrt(
     . xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*
     . xkt**2-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkt+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt
     . -1.0)/(sqrt(xkt)*sqrt(xgt)))+ans3
      ans1=6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))*xgt*xkb*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt**2+6.0*sqrt(xkt)*sqrt(
     . xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt-2.0
     . *sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(
     . xgt)))*xkb**3+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb**2*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan(
     . (x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2-6.0*sqrt(xkt)*
     . sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*
     . xkt**2+12.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkb*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb+2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**3-6.0*sqrt(
     . xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*
     . xkt**2+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkt-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt
     . -1.0)/(sqrt(xkt)*sqrt(xgt)))+ans2
      f3=ans1/(2.0*xgt*xkt)
*-----------------------------------------------------------------------
*END 
*-----------------------------------------------------------------------
      RETURN
      END



      SUBROUTINE FILLSMBRS(MHSM,GAMBRN,NMNH)
************************************************************************
*
*      RAUX_H(600)=MHSM
*
*--------------------------------------------------------
*      Width in GeV : BR         : Decay mode
*--------------------------------------------------------
*      RAUX_H(601)  : RAUX (651) : H -> e    e  [  1]
*      RAUX_H(602)  : RAUX (652) : H -> mu   mu [  2]
*      RAUX_H(603)  : RAUX (653) : H -> tau  tau[  3]
*      RAUX_H(604)  : RAUX (654) : H -> d    d  [  4]
*      RAUX_H(605)  : RAUX (655) : H -> s    s  [  5]
*      RAUX_H(606)  : RAUX (656) : H -> b    b  [  6]
*      RAUX_H(607)  : RAUX (657) : H -> u    u  [  7]
*      RAUX_H(608)  : RAUX (658) : H -> c    c  [  8]
*      RAUX_H(609)  : RAUX (659) : H -> t    t  [  9]
*      RAUX_H(610)  : RAUX (660) : H -> W    W  [ 10]
*      RAUX_H(611)  : RAUX (661) : H -> Z    Z  [ 11]
*      RAUX_H(617)  : RAUX (667) : H -> ph   ph [ 17]
*      RAUX_H(618)  : RAUX (668) : H -> gl   gl [ 18]
*      RAUX_H(619)  : RAUX (669) : H -> Z    ph [ 19]
*      RAUX_H(650)  : RAUX (700) : Total
*--------------------------------------------------------
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
* Input Array:
      REAL*8 GAMBRN(NMNH,3,3)
*-----------------------------------------------------------------------
* For integration
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
      EXTERNAL FVVS
*-----------------------------------------------------------------------
*Local
      REAL*8     SMBR(50,2)
      REAL*8     NC,I3F
      COMPLEX*16 HC_FSF,HC_F1
      COMPLEX*16 SGG,SPP
      COMPLEX*16 SPZ
      COMPLEX*16 HZP_C0,HZP_C2
*-----------------------------------------------------------------------
*
      DO I=1,50
       DO J=1,2
        SMBR(I,J)=0.D0
       ENDDO
      ENDDO
*
      PI = 2.D0*DASIN(1.D0)
*>>>>> Quarks masses at MH
      CALL QMASS_MH(MHSM,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH)
      CALL QMASS_MH(MHSM/2.D0,ASMH2,MTMH2,MBMH2,MCMH2
     .                             ,MSMH2,MUMH2,MDMH2)
*      print*,'>>> SMBRS :',mhsm,asmh,mtmh,mbmh,mcmh,msmh,mumh,mdmh
*      print*,'>>> SMBRS :',mhsm/2.d0,asmh2,mtmh2,mbmh2,mcmh2
*=======================================================================
*---> H_IH -> e+ e-     [IM= 1]
      IM   = 1
      GFF  = GW_H*ME_H/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = ME_H
      MK   = ME_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
      SMBR(IM,1)=HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> mu+ mu-   [IM= 2]
      IM   = 2
      GFF  = GW_H*MMU_H/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MMU_H
      MK   = MMU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
      SMBR(IM,1)=HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> tau+ tau-   [IM= 3]
      IM   = 3
      GFF  = GW_H*MTAU_H/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MTAU_H
      MK   = MTAU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
      SMBR(IM,1)=HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> d d       [IM= 4]
      IM   = 4
      GFF  = GW_H*MDMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MDMH
      MK   = MDMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> s s       [IM= 5]
      IM   = 5
      GFF  = GW_H*MSMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MSMH
      MK   = MSMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> b b       [IM= 6]
      IM   = 6
      GFF  = GW_H*MBMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MBMH
      MK   = MBMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*      AAA=3.D0*GF_H*MHSM/4.D0/DSQRT(2.D0)/PI*MBMH**2
*      print*,'hff_res,aaa
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*For consistency, we only keep the leading order
*inlcuding the higher-order QCD corrections
*      QCDCOR_1=5.67D0*ASMH/PI
*      QCDCOR_2=(35.94D0-1.36*5.D0)*(ASMH/PI)**2
*      QCDCOR_3=(164.14D0-25.77D0*5.D0+0.259D0*5.D0**2)*(ASMH/PI)**3
*      QCDCOR_T=(ASMH/PI)**2
*     .        *( 1.57D0-2.D0/3.D0*DLOG(MHSM**2/MTMH**2)
*     .          +(DLOG(MBMH**2/MHSM**2))**2/9.D0 )
**      print*,qcdcor_1,qcdcor_2,qcdcor_3,qcdcor_t
*      QCDCOR=QCDCOR_1+QCDCOR_2+QCDCOR_3+QCDCOR_T
*      SMBR(IM,1)=(1.D0+QCDCOR)*HFF_RES
*---> H_IH -> u u       [IM= 7]
      IM   = 7
      GFF  = GW_H*MUMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MUMH
      MK   = MUMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> c c       [IM= 8]
      IM   = 8
      GFF  = GW_H*MCMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MCMH
      MK   = MCMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> t t       [IM= 9]
      IM   = 9
      GFF  = GW_H*MTMH/2.D0/MW_H
      GSF  = 1.D0
      GPF  = 0.D0
      MJ   = MTMH
      MK   = MTMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF_SM(CF,SYMF,DJK,GFF,GSF,GPF,MHSM,MJ,MK,HFF_RES)
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      SMBR(IM,1)=(1.D0+5.67D0*ASMH/PI)*HFF_RES
*      print*,im,smbr(im,1)
*---> H_IH -> W W       [IM=10]
      IM  = 10
      DV  = 2.D0
      CALL HVV(GF_H,DV,DCMPLX(1.D0,0.D0),MW_H,MHSM,SMBR(IM,1))
* into W + W^*
      DVVS=2.D0
      IF( (MHSM.GT.MW_H) .AND.
     .    (MHSM.LT.(2.D0*MW_H+GAMW_H)) ) THEN
       IF(DABS(MHSM-2.D0*MW_H).LT.GAMW_H) THEN
        DVVS=2.D0-(MHSM-2.D0*MW_H+GAMW_H)/2.D0/GAMW_H
       ENDIF
      EPSV=GAMW_H/MW_H
      OMEGAI=MHSM**2/MW_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      SMBR(IM,1)=GF_H*MHSM**3*DV*DVVS
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES ! 3body
*
      CALL WW4BODY(MHSM,WSWS) 
      SMBR(IM,1)=WSWS ! 4body
*
      ENDIF
*
*      print*,im,smbr(im,1)
*---> H_IH -> Z Z       [IM=11]
      IM  = 11
      DV  = 1.D0
      CALL HVV(GF_H,DV,DCMPLX(1.D0,0.D0),MZ_H,MHSM,SMBR(IM,1))
* into Z + Z^*
      DVVS=2.D0
      IF( (MHSM.GT.MZ_H) .AND.
     .    (MHSM.LT.(2.D0*MZ_H+GAMZ_H)) ) THEN
       IF(DABS(MHSM-2.D0*MZ_H).LT.GAMZ_H) THEN
        DVVS=2.D0-(MHSM-2.D0*MZ_H+GAMZ_H)/2.D0/GAMZ_H
       ENDIF
      EPSV=GAMZ_H/MZ_H
      OMEGAI=MHSM**2/MZ_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      SMBR(IM,1)=GF_H*MHSM**3*DV*DVVS
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES  ! 3body
*
      CALL ZZ4BODY(MHSM,ZSZS)
      SMBR(IM,1)=ZSZS ! 4body
*
      ENDIF
*
*      print*,im,smbr(im,1)
*--->
      IM=12
      SMBR(IM,1)=0.D0
      IM=13
      SMBR(IM,1)=0.D0
      IM=14
      SMBR(IM,1)=0.D0
      IM=15
      SMBR(IM,1)=0.D0
      IM=16
      SMBR(IM,1)=0.D0
* ---> H_IH -> P P       [IM=17]
      IM=17
      QB=-1.D0/3.D0
      QT= 2.D0/3.D0
      QC= 2.D0/3.D0
      DJT = -ASMH2/PI
      AEM_0=1.D0/137.D0
      SPP=0.D0
     .   +2.D0*3.D0*QB**2*HC_FSF(MHSM**2/4.D0/MBMH2**2)*(1.D0+DJT)
     .   +2.D0*3.D0*QT**2*HC_FSF(MHSM**2/4.D0/MTMH2**2)*(1.D0+DJT)
     .   +2.D0*3.D0*QC**2*HC_FSF(MHSM**2/4.D0/MCMH2**2)*(1.D0+DJT)
     .   +2.D0*HC_FSF(MHSM**2/4.D0/MTAU_H**2)
     .   -HC_F1(MHSM**2/4.D0/MW_H**2)
      SMBR(IM,1)=MHSM**3*AEM_0**2/256.D0/PI**3/V_H**2*CDABS(SPP)**2
*      print*,im,smbr(im,1)
*      if(mhsm.eq.125.d0) then
*       print*,'MH = ',mhsm,' : AS(MH,MH/2)=',asmh,asmh2
*       print*,'POLE MASSES',mc_pole,mb_pole,mtpole_h
*       print*,'_at_POLE   ',raux_h(5),raux_h(2),mtmt_h
*       print*,'_at_MH     ',mcmh,mbmh,mtmh
*       print*,'_at_MH/2   ',mcmh2,mbmh2,mtmh2
*       print*,'_at_MT     ',mcmt_h,mbmt_h,mtmt_h
*      endif
* ---> H_IH -> G G       [IM=18]
      IM=18
      MB_POLE=RAUX_H(1)
      SGG=HC_FSF(MHSM**2/4.D0/MB_POLE**2)
     .   +HC_FSF(MHSM**2/4.D0/MTPOLE_H**2)
* running alpha_s(m_h) effect and K factor included
      XNF=5.D0
      CKHG=1.D0+ASMH/PI*(95.D0/4.D0-7.D0/6.D0*XNF)
      SMBR(IM,1)=MHSM**3*ASMH**2/32.D0/PI**3/V_H**2
     . *CKHG*CDABS(SGG)**2
*      print*,im,smbr(im,1)
*
* ---> H_IH -> Z p       [IM=19]
      IM=19
      SPZ=DCMPLX(0.D0,0.D0)
*W
      MX=MW_H
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MHSM**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
      SPZ=MZ_H**2/TW_H
     .   *(2.D0*( MHSM**2/MW_H**2*(1.D0-2.D0*CW_H**2)
     .             +2.D0*(1.D0-6.D0*CW_H**2))*HZP_C2
     .    +4.D0*(1.D0-4.D0*CW_H**2)*HZP_C0)
*      print*,spz
*top
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
      MX=MTMH2                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0-ASMH2/PI
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MHSM**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
      SPZ=SPZ
     .   +(2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .         *(HZP_C0+4.D0*HZP_C2)*QCD_COR)
*      print*,spz
*bottom
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
      MX=MBMH2                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0-ASMH2/PI
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MHSM**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
      SPZ=SPZ
     .   +(2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .         *(HZP_C0+4.D0*HZP_C2)*QCD_COR)
*tau
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
      MX=MTAU_H                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MHSM**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
      SPZ=SPZ
     .   +(2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .         *(HZP_C0+4.D0*HZP_C2)*QCD_COR)
*
      SMBR(IM,1)=AEM_H**2*MHSM**3/128.D0/PI**3/V_H**2
     .        *(1.D0-MZ_H**2/MHSM**2)**3
     .        *CDABS(SPZ)**2
      IF(SMBR(IM,1).LT.0.D0) SMBR(IM,1)=0.D0
*-----------------------------------------------------------------------
      GAMTOT=0.D0
      DO I=1,49
       GAMTOT=GAMTOT+SMBR(I,1)
      ENDDO
      SMBR(50,1)=GAMTOT
*
      DO I=1,50
       SMBR(I,2)=SMBR(I,1)/SMBR(50,1)
      ENDDO
*
      DO I=1,50
       IF(SMBR(I,1).GT.0.D0) THEN
*        WRITE(*,2) I,SMBR(I,1),SMBR(I,2),gambrn(i,1,1),gambrn(i,3,1)
       ENDIF
      ENDDO
*-----------------------------------------------------------------------
*
      RAUX_H(600)=MHSM
      DO I=1,50
       RAUX_H(600+I)=SMBR(I,1)
       RAUX_H(650+I)=SMBR(I,2)
*       print*,raux_h(600),i,raux_h(600+I),raux_h(650+I)
      ENDDO
*-----------------------------------------------------------------------
  2   FORMAT(1X,'[',I2,']',1X,2(2X,E10.4),1X,' : ',2(2X,E10.4))
*
      RETURN
      END


      SUBROUTINE HFF_SM(CF,SF,DJK,GF,GS,GP,MH,MJ,MK,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      PI=2.D0*DASIN(1.D0)
*
      XKJ=MJ**2/MH**2
      XKK=MK**2/MH**2
      XLAM=(1.D0-XKJ-XKK)**2-4.D0*XKJ*XKK
      IF(XLAM.LE.0.D0 .OR. MH.LE.(MJ+MK)) THEN
       GAM=0.D0
       RETURN
      ELSE
       GAM=SF**2/(1.D0+DJK)*CF*DABS(GF)**2*MH*DSQRT(XLAM)*
     .    ((1.D0-XKJ-XKK)*(DABS(GS)**2+DABS(GP)**2)-
     .     2.D0*DSQRT(XKJ*XKK)*(DABS(GS)**2-DABS(GP)**2))/8.D0/PI
*
      ENDIF
*
      RETURN
      END

      SUBROUTINE ZZ4BODY(MH,WIDTH)
************************************************************************
*
* The 4-body decays width in GeV for MH=50..200 GeV in the step of 0.5 GeV
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8     ZZ4BD(301)
*
       ZZ4BD(  1)=0.850582E-08
       ZZ4BD(  2)=0.915672E-08
       ZZ4BD(  3)=0.985130E-08
       ZZ4BD(  4)=0.105918E-07
       ZZ4BD(  5)=0.113811E-07
       ZZ4BD(  6)=0.122222E-07
       ZZ4BD(  7)=0.131180E-07
       ZZ4BD(  8)=0.140715E-07
       ZZ4BD(  9)=0.150861E-07
       ZZ4BD( 10)=0.161653E-07
       ZZ4BD( 11)=0.173123E-07
       ZZ4BD( 12)=0.185309E-07
       ZZ4BD( 13)=0.198247E-07
       ZZ4BD( 14)=0.211978E-07
       ZZ4BD( 15)=0.226549E-07
       ZZ4BD( 16)=0.242009E-07
       ZZ4BD( 17)=0.258400E-07
       ZZ4BD( 18)=0.275778E-07
       ZZ4BD( 19)=0.294196E-07
       ZZ4BD( 20)=0.313713E-07
       ZZ4BD( 21)=0.334374E-07
       ZZ4BD( 22)=0.356244E-07
       ZZ4BD( 23)=0.379386E-07
       ZZ4BD( 24)=0.403864E-07
       ZZ4BD( 25)=0.429732E-07
       ZZ4BD( 26)=0.457110E-07
       ZZ4BD( 27)=0.486055E-07
       ZZ4BD( 28)=0.516640E-07
       ZZ4BD( 29)=0.548935E-07
       ZZ4BD( 30)=0.583061E-07
       ZZ4BD( 31)=0.619064E-07
       ZZ4BD( 32)=0.657074E-07
       ZZ4BD( 33)=0.697151E-07
       ZZ4BD( 34)=0.739464E-07
       ZZ4BD( 35)=0.784127E-07
       ZZ4BD( 36)=0.831254E-07
       ZZ4BD( 37)=0.880906E-07
       ZZ4BD( 38)=0.933212E-07
       ZZ4BD( 39)=0.988377E-07
       ZZ4BD( 40)=0.104648E-06
       ZZ4BD( 41)=0.110773E-06
       ZZ4BD( 42)=0.117229E-06
       ZZ4BD( 43)=0.124021E-06
       ZZ4BD( 44)=0.131159E-06
       ZZ4BD( 45)=0.138677E-06
       ZZ4BD( 46)=0.146584E-06
       ZZ4BD( 47)=0.154913E-06
       ZZ4BD( 48)=0.163675E-06
       ZZ4BD( 49)=0.172894E-06
       ZZ4BD( 50)=0.182585E-06
       ZZ4BD( 51)=0.192764E-06
       ZZ4BD( 52)=0.203465E-06
       ZZ4BD( 53)=0.214750E-06
       ZZ4BD( 54)=0.226603E-06
       ZZ4BD( 55)=0.239033E-06
       ZZ4BD( 56)=0.252094E-06
       ZZ4BD( 57)=0.265838E-06
       ZZ4BD( 58)=0.280300E-06
       ZZ4BD( 59)=0.295476E-06
       ZZ4BD( 60)=0.311429E-06
       ZZ4BD( 61)=0.328197E-06
       ZZ4BD( 62)=0.345838E-06
       ZZ4BD( 63)=0.364373E-06
       ZZ4BD( 64)=0.383869E-06
       ZZ4BD( 65)=0.404346E-06
       ZZ4BD( 66)=0.425887E-06
       ZZ4BD( 67)=0.448507E-06
       ZZ4BD( 68)=0.472223E-06
       ZZ4BD( 69)=0.497125E-06
       ZZ4BD( 70)=0.523350E-06
       ZZ4BD( 71)=0.550794E-06
       ZZ4BD( 72)=0.579645E-06
       ZZ4BD( 73)=0.609919E-06
       ZZ4BD( 74)=0.641662E-06
       ZZ4BD( 75)=0.675053E-06
       ZZ4BD( 76)=0.710226E-06
       ZZ4BD( 77)=0.747268E-06
       ZZ4BD( 78)=0.786225E-06
       ZZ4BD( 79)=0.827367E-06
       ZZ4BD( 80)=0.870630E-06
       ZZ4BD( 81)=0.916272E-06
       ZZ4BD( 82)=0.964161E-06
       ZZ4BD( 83)=0.101444E-05
       ZZ4BD( 84)=0.106748E-05
       ZZ4BD( 85)=0.112351E-05
       ZZ4BD( 86)=0.118246E-05
       ZZ4BD( 87)=0.124519E-05
       ZZ4BD( 88)=0.131116E-05
       ZZ4BD( 89)=0.138180E-05
       ZZ4BD( 90)=0.145579E-05
       ZZ4BD( 91)=0.153417E-05
       ZZ4BD( 92)=0.161833E-05
       ZZ4BD( 93)=0.170842E-05
       ZZ4BD( 94)=0.180452E-05
       ZZ4BD( 95)=0.190790E-05
       ZZ4BD( 96)=0.201882E-05
       ZZ4BD( 97)=0.213846E-05
       ZZ4BD( 98)=0.227025E-05
       ZZ4BD( 99)=0.241252E-05
       ZZ4BD(100)=0.256420E-05
       ZZ4BD(101)=0.272911E-05
       ZZ4BD(102)=0.290874E-05
       ZZ4BD(103)=0.310882E-05
       ZZ4BD(104)=0.332339E-05
       ZZ4BD(105)=0.356117E-05
       ZZ4BD(106)=0.382301E-05
       ZZ4BD(107)=0.410785E-05
       ZZ4BD(108)=0.442227E-05
       ZZ4BD(109)=0.476255E-05
       ZZ4BD(110)=0.512637E-05
       ZZ4BD(111)=0.552088E-05
       ZZ4BD(112)=0.595573E-05
       ZZ4BD(113)=0.642919E-05
       ZZ4BD(114)=0.694596E-05
       ZZ4BD(115)=0.751900E-05
       ZZ4BD(116)=0.814332E-05
       ZZ4BD(117)=0.881451E-05
       ZZ4BD(118)=0.955077E-05
       ZZ4BD(119)=0.103378E-04
       ZZ4BD(120)=0.111734E-04
       ZZ4BD(121)=0.121067E-04
       ZZ4BD(122)=0.131155E-04
       ZZ4BD(123)=0.142124E-04
       ZZ4BD(124)=0.153946E-04
       ZZ4BD(125)=0.166506E-04
       ZZ4BD(126)=0.180254E-04
       ZZ4BD(127)=0.194691E-04
       ZZ4BD(128)=0.210544E-04
       ZZ4BD(129)=0.227392E-04
       ZZ4BD(130)=0.245560E-04
       ZZ4BD(131)=0.264500E-04
       ZZ4BD(132)=0.285173E-04
       ZZ4BD(133)=0.307460E-04
       ZZ4BD(134)=0.330939E-04
       ZZ4BD(135)=0.355500E-04
       ZZ4BD(136)=0.382294E-04
       ZZ4BD(137)=0.410736E-04
       ZZ4BD(138)=0.441123E-04
       ZZ4BD(139)=0.473098E-04
       ZZ4BD(140)=0.507233E-04
       ZZ4BD(141)=0.543765E-04
       ZZ4BD(142)=0.582581E-04
       ZZ4BD(143)=0.623807E-04
       ZZ4BD(144)=0.667002E-04
       ZZ4BD(145)=0.712900E-04
       ZZ4BD(146)=0.761755E-04
       ZZ4BD(147)=0.814138E-04
       ZZ4BD(148)=0.868620E-04
       ZZ4BD(149)=0.925556E-04
       ZZ4BD(150)=0.985075E-04
       ZZ4BD(151)=0.104757E-03
       ZZ4BD(152)=0.111521E-03
       ZZ4BD(153)=0.118648E-03
       ZZ4BD(154)=0.126220E-03
       ZZ4BD(155)=0.134090E-03
       ZZ4BD(156)=0.142410E-03
       ZZ4BD(157)=0.151149E-03
       ZZ4BD(158)=0.160439E-03
       ZZ4BD(159)=0.170129E-03
       ZZ4BD(160)=0.180195E-03
       ZZ4BD(161)=0.190729E-03
       ZZ4BD(162)=0.201654E-03
       ZZ4BD(163)=0.213263E-03
       ZZ4BD(164)=0.225459E-03
       ZZ4BD(165)=0.238224E-03
       ZZ4BD(166)=0.251500E-03
       ZZ4BD(167)=0.265387E-03
       ZZ4BD(168)=0.279882E-03
       ZZ4BD(169)=0.295029E-03
       ZZ4BD(170)=0.310757E-03
       ZZ4BD(171)=0.327739E-03
       ZZ4BD(172)=0.345837E-03
       ZZ4BD(173)=0.364552E-03
       ZZ4BD(174)=0.384113E-03
       ZZ4BD(175)=0.404260E-03
       ZZ4BD(176)=0.425609E-03
       ZZ4BD(177)=0.447618E-03
       ZZ4BD(178)=0.470769E-03
       ZZ4BD(179)=0.494803E-03
       ZZ4BD(180)=0.519790E-03
       ZZ4BD(181)=0.545948E-03
       ZZ4BD(182)=0.573361E-03
       ZZ4BD(183)=0.602287E-03
       ZZ4BD(184)=0.632715E-03
       ZZ4BD(185)=0.663497E-03
       ZZ4BD(186)=0.695449E-03
       ZZ4BD(187)=0.728838E-03
       ZZ4BD(188)=0.763956E-03
       ZZ4BD(189)=0.800822E-03
       ZZ4BD(190)=0.838719E-03
       ZZ4BD(191)=0.878352E-03
       ZZ4BD(192)=0.920337E-03
       ZZ4BD(193)=0.964928E-03
       ZZ4BD(194)=0.101164E-02
       ZZ4BD(195)=0.105894E-02
       ZZ4BD(196)=0.110763E-02
       ZZ4BD(197)=0.115822E-02
       ZZ4BD(198)=0.121074E-02
       ZZ4BD(199)=0.126526E-02
       ZZ4BD(200)=0.132207E-02
       ZZ4BD(201)=0.138124E-02
       ZZ4BD(202)=0.144412E-02
       ZZ4BD(203)=0.151021E-02
       ZZ4BD(204)=0.157828E-02
       ZZ4BD(205)=0.164925E-02
       ZZ4BD(206)=0.172343E-02
       ZZ4BD(207)=0.180110E-02
       ZZ4BD(208)=0.188243E-02
       ZZ4BD(209)=0.196891E-02
       ZZ4BD(210)=0.205941E-02
       ZZ4BD(211)=0.215312E-02
       ZZ4BD(212)=0.225055E-02
       ZZ4BD(213)=0.235295E-02
       ZZ4BD(214)=0.245811E-02
       ZZ4BD(215)=0.256838E-02
       ZZ4BD(216)=0.268448E-02
       ZZ4BD(217)=0.280553E-02
       ZZ4BD(218)=0.293367E-02
       ZZ4BD(219)=0.306743E-02
       ZZ4BD(220)=0.320776E-02
       ZZ4BD(221)=0.335427E-02
       ZZ4BD(222)=0.350712E-02
       ZZ4BD(223)=0.366633E-02
       ZZ4BD(224)=0.383411E-02
       ZZ4BD(225)=0.401039E-02
       ZZ4BD(226)=0.419512E-02
       ZZ4BD(227)=0.439104E-02
       ZZ4BD(228)=0.459913E-02
       ZZ4BD(229)=0.481751E-02
       ZZ4BD(230)=0.504531E-02
       ZZ4BD(231)=0.528769E-02
       ZZ4BD(232)=0.554218E-02
       ZZ4BD(233)=0.581075E-02
       ZZ4BD(234)=0.609595E-02
       ZZ4BD(235)=0.639346E-02
       ZZ4BD(236)=0.671161E-02
       ZZ4BD(237)=0.705609E-02
       ZZ4BD(238)=0.742485E-02
       ZZ4BD(239)=0.781808E-02
       ZZ4BD(240)=0.823843E-02
       ZZ4BD(241)=0.868774E-02
       ZZ4BD(242)=0.916536E-02
       ZZ4BD(243)=0.967067E-02
       ZZ4BD(244)=0.102136E-01
       ZZ4BD(245)=0.108035E-01
       ZZ4BD(246)=0.114384E-01
       ZZ4BD(247)=0.121222E-01
       ZZ4BD(248)=0.128671E-01
       ZZ4BD(249)=0.136764E-01
       ZZ4BD(250)=0.145695E-01
       ZZ4BD(251)=0.155576E-01
       ZZ4BD(252)=0.166527E-01
       ZZ4BD(253)=0.178631E-01
       ZZ4BD(254)=0.192278E-01
       ZZ4BD(255)=0.207598E-01
       ZZ4BD(256)=0.225137E-01
       ZZ4BD(257)=0.245002E-01
       ZZ4BD(258)=0.267336E-01
       ZZ4BD(259)=0.293474E-01
       ZZ4BD(260)=0.324839E-01
       ZZ4BD(261)=0.362050E-01
       ZZ4BD(262)=0.406818E-01
       ZZ4BD(263)=0.460497E-01
       ZZ4BD(264)=0.524157E-01
       ZZ4BD(265)=0.599736E-01
       ZZ4BD(266)=0.684090E-01
       ZZ4BD(267)=0.776280E-01
       ZZ4BD(268)=0.875789E-01
       ZZ4BD(269)=0.978273E-01
       ZZ4BD(270)=0.108060E+00
       ZZ4BD(271)=0.118138E+00
       ZZ4BD(272)=0.128045E+00
       ZZ4BD(273)=0.137692E+00
       ZZ4BD(274)=0.147087E+00
       ZZ4BD(275)=0.156279E+00
       ZZ4BD(276)=0.165170E+00
       ZZ4BD(277)=0.173801E+00
       ZZ4BD(278)=0.182226E+00
       ZZ4BD(279)=0.190449E+00
       ZZ4BD(280)=0.198573E+00
       ZZ4BD(281)=0.206546E+00
       ZZ4BD(282)=0.214355E+00
       ZZ4BD(283)=0.222053E+00
       ZZ4BD(284)=0.229573E+00
       ZZ4BD(285)=0.236998E+00
       ZZ4BD(286)=0.244475E+00
       ZZ4BD(287)=0.251628E+00
       ZZ4BD(288)=0.258893E+00
       ZZ4BD(289)=0.266171E+00
       ZZ4BD(290)=0.273106E+00
       ZZ4BD(291)=0.280102E+00
       ZZ4BD(292)=0.287080E+00
       ZZ4BD(293)=0.294060E+00
       ZZ4BD(294)=0.301037E+00
       ZZ4BD(295)=0.307934E+00
       ZZ4BD(296)=0.314711E+00
       ZZ4BD(297)=0.321461E+00
       ZZ4BD(298)=0.328308E+00
       ZZ4BD(299)=0.335163E+00
       ZZ4BD(300)=0.341782E+00
       ZZ4BD(301)=0.348398E+00
*
       IMH=INT(2.D0*MH-100.D0)+1
       MHGRID=0.5D0*DBLE(IMH-1)+50.D0
       WIDTH=ZZ4BD(IMH)+(ZZ4BD(IMH+1)-ZZ4BD(IMH))*(MH-MHGRID)/0.5D0
*       print*,'>>> ZZ4BODY',imh,mh,mhgrid,width,zz4bd(imh),zz4bd(imh+1)
*
      RETURN
      END

      SUBROUTINE WW4BODY(MH,WIDTH)
************************************************************************
*
* The 4-body decays width in GeV for MH=50..200 GeV in the step of 0.5 GeV
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8     WW4BD(301)
*
       WW4BD(  1)=0.280840E-07
       WW4BD(  2)=0.302768E-07
       WW4BD(  3)=0.326207E-07
       WW4BD(  4)=0.351252E-07
       WW4BD(  5)=0.378012E-07
       WW4BD(  6)=0.406583E-07
       WW4BD(  7)=0.437070E-07
       WW4BD(  8)=0.469579E-07
       WW4BD(  9)=0.504232E-07
       WW4BD( 10)=0.541150E-07
       WW4BD( 11)=0.580512E-07
       WW4BD( 12)=0.622446E-07
       WW4BD( 13)=0.667077E-07
       WW4BD( 14)=0.714580E-07
       WW4BD( 15)=0.765105E-07
       WW4BD( 16)=0.818829E-07
       WW4BD( 17)=0.875916E-07
       WW4BD( 18)=0.936611E-07
       WW4BD( 19)=0.100115E-06
       WW4BD( 20)=0.106975E-06
       WW4BD( 21)=0.114252E-06
       WW4BD( 22)=0.121977E-06
       WW4BD( 23)=0.130179E-06
       WW4BD( 24)=0.138882E-06
       WW4BD( 25)=0.148122E-06
       WW4BD( 26)=0.157913E-06
       WW4BD( 27)=0.168280E-06
       WW4BD( 28)=0.179269E-06
       WW4BD( 29)=0.190912E-06
       WW4BD( 30)=0.203256E-06
       WW4BD( 31)=0.216338E-06
       WW4BD( 32)=0.230182E-06
       WW4BD( 33)=0.244827E-06
       WW4BD( 34)=0.260330E-06
       WW4BD( 35)=0.276789E-06
       WW4BD( 36)=0.294174E-06
       WW4BD( 37)=0.312547E-06
       WW4BD( 38)=0.331998E-06
       WW4BD( 39)=0.352609E-06
       WW4BD( 40)=0.374401E-06
       WW4BD( 41)=0.397457E-06
       WW4BD( 42)=0.421855E-06
       WW4BD( 43)=0.447696E-06
       WW4BD( 44)=0.475048E-06
       WW4BD( 45)=0.503971E-06
       WW4BD( 46)=0.534628E-06
       WW4BD( 47)=0.567032E-06
       WW4BD( 48)=0.601224E-06
       WW4BD( 49)=0.637409E-06
       WW4BD( 50)=0.675716E-06
       WW4BD( 51)=0.716078E-06
       WW4BD( 52)=0.758810E-06
       WW4BD( 53)=0.803904E-06
       WW4BD( 54)=0.851600E-06
       WW4BD( 55)=0.902215E-06
       WW4BD( 56)=0.955866E-06
       WW4BD( 57)=0.101278E-05
       WW4BD( 58)=0.107318E-05
       WW4BD( 59)=0.113737E-05
       WW4BD( 60)=0.120519E-05
       WW4BD( 61)=0.127691E-05
       WW4BD( 62)=0.135307E-05
       WW4BD( 63)=0.143411E-05
       WW4BD( 64)=0.152016E-05
       WW4BD( 65)=0.161208E-05
       WW4BD( 66)=0.171038E-05
       WW4BD( 67)=0.181519E-05
       WW4BD( 68)=0.192635E-05
       WW4BD( 69)=0.204680E-05
       WW4BD( 70)=0.217713E-05
       WW4BD( 71)=0.231771E-05
       WW4BD( 72)=0.246956E-05
       WW4BD( 73)=0.263545E-05
       WW4BD( 74)=0.281899E-05
       WW4BD( 75)=0.301994E-05
       WW4BD( 76)=0.323738E-05
       WW4BD( 77)=0.347613E-05
       WW4BD( 78)=0.374224E-05
       WW4BD( 79)=0.404000E-05
       WW4BD( 80)=0.436379E-05
       WW4BD( 81)=0.473075E-05
       WW4BD( 82)=0.513678E-05
       WW4BD( 83)=0.558691E-05
       WW4BD( 84)=0.608131E-05
       WW4BD( 85)=0.661854E-05
       WW4BD( 86)=0.720828E-05
       WW4BD( 87)=0.786556E-05
       WW4BD( 88)=0.859165E-05
       WW4BD( 89)=0.940109E-05
       WW4BD( 90)=0.103051E-04
       WW4BD( 91)=0.112863E-04
       WW4BD( 92)=0.123725E-04
       WW4BD( 93)=0.135459E-04
       WW4BD( 94)=0.148239E-04
       WW4BD( 95)=0.162487E-04
       WW4BD( 96)=0.178280E-04
       WW4BD( 97)=0.195412E-04
       WW4BD( 98)=0.213702E-04
       WW4BD( 99)=0.233980E-04
       WW4BD(100)=0.255751E-04
       WW4BD(101)=0.279666E-04
       WW4BD(102)=0.305374E-04
       WW4BD(103)=0.332682E-04
       WW4BD(104)=0.362552E-04
       WW4BD(105)=0.395346E-04
       WW4BD(106)=0.430143E-04
       WW4BD(107)=0.466768E-04
       WW4BD(108)=0.507103E-04
       WW4BD(109)=0.550262E-04
       WW4BD(110)=0.597024E-04
       WW4BD(111)=0.646305E-04
       WW4BD(112)=0.699552E-04
       WW4BD(113)=0.757009E-04
       WW4BD(114)=0.818506E-04
       WW4BD(115)=0.883687E-04
       WW4BD(116)=0.953488E-04
       WW4BD(117)=0.102815E-03
       WW4BD(118)=0.110801E-03
       WW4BD(119)=0.119154E-03
       WW4BD(120)=0.127888E-03
       WW4BD(121)=0.137295E-03
       WW4BD(122)=0.147435E-03
       WW4BD(123)=0.158232E-03
       WW4BD(124)=0.169725E-03
       WW4BD(125)=0.181874E-03
       WW4BD(126)=0.194585E-03
       WW4BD(127)=0.208257E-03
       WW4BD(128)=0.222613E-03
       WW4BD(129)=0.237768E-03
       WW4BD(130)=0.253377E-03
       WW4BD(131)=0.270150E-03
       WW4BD(132)=0.287937E-03
       WW4BD(133)=0.306540E-03
       WW4BD(134)=0.326072E-03
       WW4BD(135)=0.346613E-03
       WW4BD(136)=0.368367E-03
       WW4BD(137)=0.391158E-03
       WW4BD(138)=0.415243E-03
       WW4BD(139)=0.441622E-03
       WW4BD(140)=0.468966E-03
       WW4BD(141)=0.497777E-03
       WW4BD(142)=0.527699E-03
       WW4BD(143)=0.559566E-03
       WW4BD(144)=0.592597E-03
       WW4BD(145)=0.627538E-03
       WW4BD(146)=0.663901E-03
       WW4BD(147)=0.702113E-03
       WW4BD(148)=0.742420E-03
       WW4BD(149)=0.785160E-03
       WW4BD(150)=0.830171E-03
       WW4BD(151)=0.875797E-03
       WW4BD(152)=0.924147E-03
       WW4BD(153)=0.975064E-03
       WW4BD(154)=0.102919E-02
       WW4BD(155)=0.108462E-02
       WW4BD(156)=0.114312E-02
       WW4BD(157)=0.120548E-02
       WW4BD(158)=0.127261E-02
       WW4BD(159)=0.134197E-02
       WW4BD(160)=0.141306E-02
       WW4BD(161)=0.148671E-02
       WW4BD(162)=0.156412E-02
       WW4BD(163)=0.164433E-02
       WW4BD(164)=0.172835E-02
       WW4BD(165)=0.181732E-02
       WW4BD(166)=0.191279E-02
       WW4BD(167)=0.201210E-02
       WW4BD(168)=0.211508E-02
       WW4BD(169)=0.222369E-02
       WW4BD(170)=0.233768E-02
       WW4BD(171)=0.245787E-02
       WW4BD(172)=0.258697E-02
       WW4BD(173)=0.272287E-02
       WW4BD(174)=0.286403E-02
       WW4BD(175)=0.301290E-02
       WW4BD(176)=0.316737E-02
       WW4BD(177)=0.332929E-02
       WW4BD(178)=0.350087E-02
       WW4BD(179)=0.367986E-02
       WW4BD(180)=0.387283E-02
       WW4BD(181)=0.407470E-02
       WW4BD(182)=0.428658E-02
       WW4BD(183)=0.451010E-02
       WW4BD(184)=0.474337E-02
       WW4BD(185)=0.499037E-02
       WW4BD(186)=0.525137E-02
       WW4BD(187)=0.552841E-02
       WW4BD(188)=0.582664E-02
       WW4BD(189)=0.614204E-02
       WW4BD(190)=0.647351E-02
       WW4BD(191)=0.682703E-02
       WW4BD(192)=0.719897E-02
       WW4BD(193)=0.759696E-02
       WW4BD(194)=0.801843E-02
       WW4BD(195)=0.846803E-02
       WW4BD(196)=0.895513E-02
       WW4BD(197)=0.948285E-02
       WW4BD(198)=0.100536E-01
       WW4BD(199)=0.106685E-01
       WW4BD(200)=0.113303E-01
       WW4BD(201)=0.120401E-01
       WW4BD(202)=0.128017E-01
       WW4BD(203)=0.136275E-01
       WW4BD(204)=0.145282E-01
       WW4BD(205)=0.155056E-01
       WW4BD(206)=0.165856E-01
       WW4BD(207)=0.177732E-01
       WW4BD(208)=0.190978E-01
       WW4BD(209)=0.205777E-01
       WW4BD(210)=0.222455E-01
       WW4BD(211)=0.241230E-01
       WW4BD(212)=0.262740E-01
       WW4BD(213)=0.287543E-01
       WW4BD(214)=0.316432E-01
       WW4BD(215)=0.349403E-01
       WW4BD(216)=0.388647E-01
       WW4BD(217)=0.437079E-01
       WW4BD(218)=0.496242E-01
       WW4BD(219)=0.570217E-01
       WW4BD(220)=0.661209E-01
       WW4BD(221)=0.772500E-01
       WW4BD(222)=0.901905E-01
       WW4BD(223)=0.104641E+00
       WW4BD(224)=0.120473E+00
       WW4BD(225)=0.136587E+00
       WW4BD(226)=0.152627E+00
       WW4BD(227)=0.168258E+00
       WW4BD(228)=0.183523E+00
       WW4BD(229)=0.198282E+00
       WW4BD(230)=0.212659E+00
       WW4BD(231)=0.226472E+00
       WW4BD(232)=0.239843E+00
       WW4BD(233)=0.252848E+00
       WW4BD(234)=0.265524E+00
       WW4BD(235)=0.278024E+00
       WW4BD(236)=0.290136E+00
       WW4BD(237)=0.302141E+00
       WW4BD(238)=0.313803E+00
       WW4BD(239)=0.325452E+00
       WW4BD(240)=0.336978E+00
       WW4BD(241)=0.348132E+00
       WW4BD(242)=0.359155E+00
       WW4BD(243)=0.370110E+00
       WW4BD(244)=0.380949E+00
       WW4BD(245)=0.391782E+00
       WW4BD(246)=0.402589E+00
       WW4BD(247)=0.413393E+00
       WW4BD(248)=0.423989E+00
       WW4BD(249)=0.434519E+00
       WW4BD(250)=0.445062E+00
       WW4BD(251)=0.455737E+00
       WW4BD(252)=0.466043E+00
       WW4BD(253)=0.476357E+00
       WW4BD(254)=0.486718E+00
       WW4BD(255)=0.496939E+00
       WW4BD(256)=0.507283E+00
       WW4BD(257)=0.517688E+00
       WW4BD(258)=0.527896E+00
       WW4BD(259)=0.538414E+00
       WW4BD(260)=0.548685E+00
       WW4BD(261)=0.558930E+00
       WW4BD(262)=0.569419E+00
       WW4BD(263)=0.579868E+00
       WW4BD(264)=0.590358E+00
       WW4BD(265)=0.600943E+00
       WW4BD(266)=0.611299E+00
       WW4BD(267)=0.621573E+00
       WW4BD(268)=0.632484E+00
       WW4BD(269)=0.642815E+00
       WW4BD(270)=0.653535E+00
       WW4BD(271)=0.664527E+00
       WW4BD(272)=0.674944E+00
       WW4BD(273)=0.685469E+00
       WW4BD(274)=0.695985E+00
       WW4BD(275)=0.706770E+00
       WW4BD(276)=0.717824E+00
       WW4BD(277)=0.728629E+00
       WW4BD(278)=0.739567E+00
       WW4BD(279)=0.750568E+00
       WW4BD(280)=0.761360E+00
       WW4BD(281)=0.772817E+00
       WW4BD(282)=0.783958E+00
       WW4BD(283)=0.795133E+00
       WW4BD(284)=0.806302E+00
       WW4BD(285)=0.817605E+00
       WW4BD(286)=0.829310E+00
       WW4BD(287)=0.840600E+00
       WW4BD(288)=0.851834E+00
       WW4BD(289)=0.863588E+00
       WW4BD(290)=0.875493E+00
       WW4BD(291)=0.886807E+00
       WW4BD(292)=0.898821E+00
       WW4BD(293)=0.911072E+00
       WW4BD(294)=0.922992E+00
       WW4BD(295)=0.935154E+00
       WW4BD(296)=0.947205E+00
       WW4BD(297)=0.959565E+00
       WW4BD(298)=0.971867E+00
       WW4BD(299)=0.984006E+00
       WW4BD(300)=0.996043E+00
       WW4BD(301)=0.100847E+01
*
       IMH=INT(2.D0*MH-100.D0)+1
       MHGRID=0.5D0*DBLE(IMH-1)+50.D0
       WIDTH=WW4BD(IMH)+(WW4BD(IMH+1)-WW4BD(IMH))*(MH-MHGRID)/0.5D0
*       print*,'>>> WW4BODY',imh,mh,mhgrid,width,ww4bd(imh),ww4bd(imh+1)
*
      RETURN
      END


      SUBROUTINE H2ZGAMMA_V0(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC,HZP_WIDTH,HZP_WIDTH_SM)
************************************************************************
*
*  SPZ[IM,IH]=S[H_IH -> Z\gamma]  and 
*  PPZ[IM,IH]=P[H_IH -> Z\gamma]
*
*  with
*
*    IM=1   : top
*    IM=2   : bottom
*    IM=3   : tau
*    IM=6   : chargino1-chargino1
*    IM=7   : chargino1-chargino2
*    IM=8   : chargino2-chargino1
*    IM=9   : chargino2-chargino2
*    IM=11  : W
*    IM=12  : H^\pm
*    IM=13  : stop1-stop1
*    IM=14  : stop1-stop2
*    IM=15  : stop2-stop1
*    IM=16  : stop2-stop2
*    IM=17  : sbottom1-sbottom1
*    IM=18  : sbottom1-sbottom2
*    IM=19  : sbottom2-sbottom1
*    IM=20  : sbottom2-sbottom2
*    IM=21  : stau1-stau1
*    IM=22  : stau1-stau2
*    IM=23  : stau2-stau1
*    IM=24  : stau2-stau2
*    IM=30  : TOTAL
*
* OUTPUT: HZP_WIDTH[IH]=Gamma(H_IH-> Z \gamma) in GeV
*         HZP_WIDTH_SM[IH]=corresponding SM width
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
*ARRAYS:
      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
*
      INTEGER*8 IFLAG_H(NFLAG)
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)  ! 100 = NCMAX
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)  
      REAL*8 GAMBRC(NMCH,3)    
*-----------------------------------------------------------------------
*Local Vairables:
*
*Higgs-Photon-Z Couplings
      COMPLEX*16 SPZ_SM(30,3)
      COMPLEX*16 SPZ(30,3),PPZ(30,3)
*
      REAL*8     NC,I3F
*
      COMPLEX*16 HZP_FW
      COMPLEX*16 HZP_C0,HZP_C2
      COMPLEX*16 HZP_FTAU,HZP_GTAU
*
      COMPLEX*16 C_I1,C_I2,C_A1,C_A12
      COMPLEX*16 CX,CRES
      COMPLEX*16 F_122,G_122
      COMPLEX*16 C0_12,C0_21,C1_12,C1_21,C2_12,C2_21
      COMPLEX*16 GZ_12,GS_12,GP_12
      COMPLEX*16 XI,CTMP1,CTMP2
*
*=======================================================================
*OUTPUT
      REAL*8     HZP_WIDTH(3),HZP_WIDTH_SM(3)
*=======================================================================
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
*      print*,'    >>>>> H2ZGAMMA: START <<<<<'
*
*=======================================================================
*Initialize
      DO I=1,30
       DO J=1,3
        SPZ(I,J)=DCMPLX(0.D0,0.D0)
        SPZ_SM(I,J)=DCMPLX(0.D0,0.D0)
        PPZ(I,J)=DCMPLX(0.D0,0.D0)
*        print*,i,j,spz(i,j),ppz(i,j)
       ENDDO 
      ENDDO
*
*=======================================================================
*
      DO IH=1,3
*
      MH    = HMASS_H(IH)
*
      CALL QMASS_MH(MH,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH)
*      print*,MH,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH
      CALL QMASS_MH(MH/2.D0,ASMH2,MTMH2,MBMH2,MCMH2
     .                           ,MSMH2,MUMH2,MDMH2)
*-----------------------------------------------------------------------
*I=11 : W
      IM=11
      MX=MW_H
*
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MH**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
*      print*,ih,hmass_h(ih),hzp_c0,hzp_c2
*
      HZP_FW=2.D0*( MH**2/MW_H**2*(1.D0-2.D0*CW_H**2)
     .             +2.D0*(1.D0-6.D0*CW_H**2))*HZP_C2
     .      +4.D0*(1.D0-4.D0*CW_H**2)*HZP_C0
      SPZ(IM,IH)=MZ_H**2/TW_H*NHC_H(70,IH)*HZP_FW
      SPZ_SM(IM,IH)=MZ_H**2/TW_H*HZP_FW
*      print*,ih,im,hmass_h(ih),spz(im,ih),spz_sm(im,ih),'<= SM'
*
*Djouadi's Review I1 and I2
      C_I1=4.D0*MX**2*HZP_C2
      C_I2=-MX**2*HZP_C0
      C_A1=CW_H*( 4.D0*(3.D0-SW_H**2/CW_H**2)*C_I2
     .           +( (1.D0+MH**2/2.D0/MX**2)*SW_H**2/CW_H**2
     .             -(5.D0+MH**2/2.D0/MX**2) )*C_I1 )
*      print*,'Djouaid Review p.94:',c_a1,-4.6d0+0.3d0*mh**2/mx**2
*-----------------------------------------------------------------------
*I=1 : top
      IM=1
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
      MX=MTMH2                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0-ASMH2/PI
*
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MH**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
*      print*,ih,hmass_h(ih),hzp_c0,hzp_c2
*
      SPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(26,IH)*(HZP_C0+4.D0*HZP_C2)*QCD_COR
      SPZ_SM(IM,IH)
     .          =2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *(HZP_C0+4.D0*HZP_C2)*QCD_COR
      PPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(27,IH)*HZP_C0*QCD_COR
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*      print*,ih,im,hmass_h(ih),ppz(im,ih)
*
*Djouadi's Review I1 and I2
      C_I1=4.D0*MX**2*HZP_C2
      C_I2=-MX**2*HZP_C0
      C_A12=C_I1-C_I2
*      print*,'Djouaid Review p.94:',c_a12
*     .      ,nc*qf*(2.d0*i3f-4.d0*qf*sw_h**2)/3.d0/cw_h,0.3d0
*-----------------------------------------------------------------------
*I=2 : bottom
      IM=2
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
      MX=MBMH2                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0-ASMH2/PI
*
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MH**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
*      print*,ih,hmass_h(ih),hzp_c0,hzp_c2
*
      SPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(17,IH)*(HZP_C0+4.D0*HZP_C2)*QCD_COR
      SPZ_SM(IM,IH)
     .          =2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *(HZP_C0+4.D0*HZP_C2)*QCD_COR
      PPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(18,IH)*HZP_C0*QCD_COR
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*      print*,ih,im,hmass_h(ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=3 : tau
      IM=3
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
      MX=MTAU_H                ! same as in Higgs-gamma-gamma
      QCD_COR=1.D0
*
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MH**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
*      print*,ih,hmass_h(ih),hzp_c0,hzp_c2
*
      SPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(8,IH)*(HZP_C0+4.D0*HZP_C2)*QCD_COR
      SPZ_SM(IM,IH)
     .          =2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *(HZP_C0+4.D0*HZP_C2)*QCD_COR
      PPZ(IM,IH)=2.D0*QF*NC*MX**2*(I3F-2.D0*SW_H**2*QF**2)/SW_H/CW_H
     .          *NHC_H(9,IH)*HZP_C0*QCD_COR
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*      print*,ih,im,hmass_h(ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=12 : charged Higgs
      IM=12
      MX=MCH 
      QCD_COR=1.D0
*
      TAUZ=4.D0*MX**2/MZ_H**2
      TAUH=4.D0*MX**2/MH**2
      CALL HZP_C0C2(MX,TAUZ,TAUH,HZP_C0,HZP_C2)
      SPZ(IM,IH)=-V_H**2/2.D0/CW_H/SW_H*NHC_H(86,IH)
     .          *4.D0*HZP_C2
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*
*      print*,'DDHK, below eq.(19) (it should be -1/24 MCH^2):'
*     .      ,hzp_c2,-1.d0/24.d0/mx**2
*
*      if(ih.eq.1) then
*      sa=-omix_h(1,1)  ! O_{phi_1,1)=-sin(alpha)
*      ca= omix_h(2,1)  ! O_{phi_2,1)= cos(alpha)
*      print*,'DDHK, eq.(24):',-2.d0*nhc_h(86,ih)/gw_h**2
*     .      ,(sb_h+ca-cb_h*sa)
*     .      +(cb_h**2-sb_h**2)*(sb_h*ca+cb_h*sa)/2.d0/cw_h**2
*      endif
*
*-----------------------------------------------------------------------
*I=6 : charggino1-chargino1
      IM=6
      I1=1
      I2=1
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=MC_H(I1)
      M2=MC_H(I2)
      IF(I1.EQ.I2) THEN
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0-CW_H**2
      ELSE
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0
      ENDIF
      IF(I1.EQ.1.AND.I2.EQ.1) THEN
       GS_12=NHC_H(59,IH)
       GP_12=NHC_H(60,IH)
      ELSEIF(I1.EQ.1.AND.I2.EQ.2) THEN
       GS_12=NHC_H(62,IH)
       GP_12=NHC_H(63,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.1) THEN
       GS_12=NHC_H(65,IH)
       GP_12=NHC_H(66,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.2) THEN
       GS_12=NHC_H(68,IH)
       GP_12=NHC_H(69,IH)
      ELSE
       print*,'>>H2ZGAMMA: Error in H-charino-chargino coupling'
      ENDIF
*      print*,ih,im,gs_12,gp_12
*
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*      print*,ih,im,m1,m2,f_122,g_122,c2_12
*
      SPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*F_122*GZ_12*GS_12
      PPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*G_122*GZ_12*GP_12*XI
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=7 : charggino1-chargino2
      IM=7
      I1=1
      I2=2
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=MC_H(I1)
      M2=MC_H(I2)
      IF(I1.EQ.I2) THEN
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0-CW_H**2
      ELSE
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0
      ENDIF
      IF(I1.EQ.1.AND.I2.EQ.1) THEN
       GS_12=NHC_H(59,IH)
       GP_12=NHC_H(60,IH)
      ELSEIF(I1.EQ.1.AND.I2.EQ.2) THEN
       GS_12=NHC_H(62,IH)
       GP_12=NHC_H(63,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.1) THEN
       GS_12=NHC_H(65,IH)
       GP_12=NHC_H(66,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.2) THEN
       GS_12=NHC_H(68,IH)
       GP_12=NHC_H(69,IH)
      ELSE
       print*,'>>H2ZGAMMA: Error in H-charino-chargino coupling'
      ENDIF
*      print*,ih,im,gs_12,gp_12
*
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      SPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*F_122*GZ_12*GS_12
      PPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*G_122*GZ_12*GP_12*XI
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=8 : charggino2-chargino1
      IM=8
      I1=2
      I2=1
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=MC_H(I1)
      M2=MC_H(I2)
      IF(I1.EQ.I2) THEN
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0-CW_H**2
      ELSE
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0
      ENDIF
      IF(I1.EQ.1.AND.I2.EQ.1) THEN
       GS_12=NHC_H(59,IH)
       GP_12=NHC_H(60,IH)
      ELSEIF(I1.EQ.1.AND.I2.EQ.2) THEN
       GS_12=NHC_H(62,IH)
       GP_12=NHC_H(63,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.1) THEN
       GS_12=NHC_H(65,IH)
       GP_12=NHC_H(66,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.2) THEN
       GS_12=NHC_H(68,IH)
       GP_12=NHC_H(69,IH)
      ELSE
       print*,'>>H2ZGAMMA: Error in H-charino-chargino coupling'
      ENDIF
*      print*,ih,im,gs_12,gp_12
*
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      SPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*F_122*GZ_12*GS_12
      PPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*G_122*GZ_12*GP_12*XI
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=9 : charggino2-chargino2
      IM=9
      I1=2
      I2=2
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=MC_H(I1)
      M2=MC_H(I2)
      IF(I1.EQ.I2) THEN
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0-CW_H**2
      ELSE
       GZ_12=( UL_H(I1,2)*DCONJG(UL_H(I2,2))
     .        +UR_H(I1,2)*DCONJG(UR_H(I2,2)) )/4.D0
      ENDIF
      IF(I1.EQ.1.AND.I2.EQ.1) THEN
       GS_12=NHC_H(59,IH)
       GP_12=NHC_H(60,IH)
      ELSEIF(I1.EQ.1.AND.I2.EQ.2) THEN
       GS_12=NHC_H(62,IH)
       GP_12=NHC_H(63,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.1) THEN
       GS_12=NHC_H(65,IH)
       GP_12=NHC_H(66,IH)
      ELSEIF(I1.EQ.2.AND.I2.EQ.2) THEN
       GS_12=NHC_H(68,IH)
       GP_12=NHC_H(69,IH)
      ELSE
       print*,'>>H2ZGAMMA: Error in H-charino-chargino coupling'
      ENDIF
*      print*,ih,im,gs_12,gp_12
*
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      SPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*F_122*GZ_12*GS_12
      PPZ(IM,IH)=-MZ_H**2*CW_H/SW_H
     .          *2.D0*DSQRT(2.D0)*M1/MW_H*G_122*GZ_12*GP_12*XI
*
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih)
*-----------------------------------------------------------------------
*I=13: stop1-stop1
      IM=13
      I1=1
      I2=1
      INHC=71  ! H_IH-stop1^*-stop1
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STMASS_H(I1)
      M2=STMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=14: stop1-stop2
      IM=14
      I1=1
      I2=2
      INHC=72  ! H_IH-stop1^*-stop2
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STMASS_H(I1)
      M2=STMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=15: stop2-stop1
      IM=15
      I1=2
      I2=1
      INHC=73  ! H_IH-stop2^*-stop1
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STMASS_H(I1)
      M2=STMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=16: stop2-stop2
      IM=16
      I1=2
      I2=2
      INHC=74  ! H_IH-stop2^*-stop1
      QF=2.D0/3.D0
      NC=3.D0
      I3F=1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STMASS_H(I1)
      M2=STMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STMIX_H(1,I2))*STMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=17: sbottom1-sbottom1
      IM=17
      I1=1
      I2=1
      INHC=75  ! H_IH-sbottom1^*-sbottom1
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=SBMASS_H(I1)
      M2=SBMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=18: sbottom1-sbottom2
      IM=18
      I1=1
      I2=2
      INHC=76  ! H_IH-sbottom1^*-sbottom2
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=SBMASS_H(I1)
      M2=SBMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=19: sbottom2-sbottom1
      IM=19
      I1=2
      I2=1
      INHC=77  ! H_IH-sbottom2^*-sbottom1
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
*     
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=SBMASS_H(I1)
      M2=SBMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=20: sbottom2-sbottom2
      IM=20
      I1=2
      I2=2
      INHC=78  ! H_IH-sbottom2^*-sbottom1
      QF=-1.D0/3.D0
      NC=3.D0
      I3F=-1.D0/2.D0
*     
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=SBMASS_H(I1)
      M2=SBMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(SBMIX_H(1,I2))*SBMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=21: stau1-stau1
      IM=21
      I1=1
      I2=1
      INHC=79  ! H_IH-stau1^*-stau1
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STAUMASS_H(I1)
      M2=STAUMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=22: stau1-stau2
      IM=22
      I1=1
      I2=2
      INHC=80  ! H_IH-stau1^*-stau2
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STAUMASS_H(I1)
      M2=STAUMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=23: stau2-stau1
      IM=23
      I1=2
      I2=1
      INHC=81  ! H_IH-stau2^*-stau1
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
*     
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STAUMASS_H(I1)
      M2=STAUMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*I=24: stau2-stau2
      IM=24
      I1=2
      I2=2
      INHC=82  ! H_IH-stau2^*-stau2
      QF=-1.D0
      NC=1.D0
      I3F=-1.D0/2.D0
*
      MUREM=100.D0 ! renormalization scale, results should be indept. of this
      M1=STAUMASS_H(I1)
      M2=STAUMASS_H(I2)
      CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
     .        ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*
      IF(I1.EQ.I2) THEN
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)-QF*SW_H**2
      ELSE
       GZ_12=I3F*DCONJG(STAUMIX_H(1,I2))*STAUMIX_H(1,I1)
      ENDIF
*
      SPZ(IM,IH)=MZ_H**2*NC*QF
     .          *(-4.D0*V_H**2)/MZ_H**2/CW_H/SW_H
     .          *NHC_H(INHC,IH)*GZ_12*C2_12
*
*      print*,ih,im,hmass_h(ih),spz(im,ih)
*-----------------------------------------------------------------------
*Total
      IM=30
      DO I=1,29
       SPZ(IM,IH)=SPZ(IM,IH)+SPZ(I,IH)
       PPZ(IM,IH)=PPZ(IM,IH)+PPZ(I,IH)
       SPZ_SM(IM,IH)=SPZ_SM(IM,IH)+SPZ_SM(I,IH)
      ENDDO
*      print*,ih,im,hmass_h(ih),spz(im,ih),ppz(im,ih),spz_sm(im,ih)
*      print*,spz(im,ih),ppz(im,ih),';',spz_sm(im,ih)
*
*-----------------------------------------------------------------------
*Widths
*      print*,AEM_H**2*MH**3/128.D0/PI**3/V_H**2 
*      print*,AEM_H*GF_H**2*MW_H**2*SW_H**2*MH**3/64.D0/PI**4
      WIDTH_SM=AEM_H**2*MH**3/128.D0/PI**3/V_H**2
     .        *(1.D0-MZ_H**2/MH**2)**3
     .        *CDABS(SPZ_SM(30,IH))**2
      IF(WIDTH_SM.LT.0.D0) WIDTH_SM=0.D0
*
*Note: H2~H_SM with "run.h2sm"
*
      WIDTH=AEM_H**2*MH**3/128.D0/PI**3/V_H**2
     .     *(1.D0-MZ_H**2/MH**2)**3
     .     *(CDABS(SPZ(30,IH))**2+CDABS(PPZ(30,IH))**2)
      IF(WIDTH.LT.0.D0) WIDTH=0.D0
*
*      print*,ih,mh,width,width_sm
*
      HZP_WIDTH(IH)=WIDTH
      HZP_WIDTH_SM(IH)=WIDTH_SM
*-----------------------------------------------------------------------
*
      ENDDO ! IH
*
*      print*,smpara_h(20),hmass_h(2)
*      raux_h(650)=raux_h(650)+raux_h(619)
*      do i=1,50
*       raux_h(650+i)=raux_h(600+i)/raux_h(650)
*      enddo
**      do i=1,50
**       if(raux_h(600+i).gt.0.d0) print*,i,raux_h(600+i)
**     .                                   ,raux_h(650+i)*1.d2
**      enddo
*       print*,'ee :',raux_h(601),raux_h(651)*1.d2
*       print*,'mm :',raux_h(602),raux_h(652)*1.d2
*       print*,'tt :',raux_h(603),raux_h(653)*1.d2
*       print*,'dd :',raux_h(604),raux_h(654)*1.d2
*       print*,'ss :',raux_h(605),raux_h(655)*1.d2
*       print*,'bb :',raux_h(606),raux_h(656)*1.d2
*       print*,'uu :',raux_h(607),raux_h(657)*1.d2
*       print*,'cc :',raux_h(608),raux_h(658)*1.d2
*       print*,'WW :',raux_h(610),raux_h(660)*1.d2
*       print*,'ZZ :',raux_h(611),raux_h(661)*1.d2
*       print*,'pp :',raux_h(617),raux_h(667)*1.d2
*       print*,'gg :',raux_h(618),raux_h(668)*1.d2
*       print*,'Zp :',raux_h(619),raux_h(669)*1.d2
*       print*,'Tot:',raux_h(650),raux_h(700)*1.d2
*=======================================================================
*
*      print*,'    >>>>> H2ZGAMMA: E N D <<<<<'
*
*--> test of f(tau) and g(tau) at large tau
*
*      tau=100.d0
*      print*,'f(tau):',hzp_ftau(tau),1.d0/tau+1.d0/3.d0/tau**2
*      print*,'g(tau):',hzp_gtau(tau),1.d0-1.d0/3.d0/tau
*
*--> test of DILOG function
*
*      CX=DCMPLX(1.D0,0.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,pi**2/6.d0
*      CX=DCMPLX(-1.D0,0.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,-pi**2/12.d0
*      CX=DCMPLX(0.D0,0.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,0.d0
*      CX=DCMPLX(0.5D0,0.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,pi**2/12.d0-(dlog(2.d0))**2/2.d0
*      CX=DCMPLX(2.D0,0.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,pi**2/4.d0-dcmplx(0.d0,1.d0)*pi*dlog(2.d0)
***
*      CATALAN=0.9159656D0 ! Catalan constant
*      CX=DCMPLX(0.D0,1.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres,dcmplx(-pi**2/48.d0,CATALAN) 
*      CX=DCMPLX(1.D0,-1.D0)
*      CALL HZP_DILOG(CX,CRES)
*      print*,'DiLog:',cx,cres
*     .,dcmplx(3.d0*pi**2/16.d0-pi**2/8.d0,-CATALAN-pi*dlog(2.d0)/4.d0) 
*
*
*--> test of C functions
*
*       MUREM=100.D0
*       MH   =125.D0
*       M2   = 20.D0
*       M1   = 15.D0
*       print*,'mu_rem=',murem,mz_h,mh,m1,m2
*       CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
*     .         ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*       print*,'C0(m1,m2,m2)=',c0_12
*       print*,'C0(m2,m1,m1)=',c0_21
*       print*,'C1(m1,m2,m2)=',c1_12
*       print*,'C1(m2,m1,m1)=',c1_21
*       print*,'C2(m1,m2,m2)=',c2_12
*       print*,'C2(m2,m1,m1)=',c2_21
**       
*       do im2=1,3
*        if(im2.eq.1) m2=140.d0
*        if(im2.eq.2) m2=115.d0
*        if(im2.eq.3) m2= 85.d0
*        do im1=1,200
*         m1=dble(im1)
*         if(m1.gt.m2) goto 97
*         CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
*     .           ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*         X1=m1
*         X2=m2
*         Y1=dreal(c0_12)
*         Y2=dimag(c0_12)
*         Y3=dreal(c0_21)
*         Y4=dimag(c0_21)
*         Y5=0.d0
*         Y6=0.D0
*         Y7=0.D0
*         Y8=0.D0
*         WRITE(10,10) X1,X2,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
*        enddo
* 97     continue
*       enddo
*       print*,f_122,g_122,c2_12
*       MUREM=200.D0
*       M1   =200.D0
*       M2   =100.D0
*       print*,'mu_rem=',murem,mz_h,mh,m1,m2
*       CALL HZP_FGC(MUREM,MZ_H,MH,M1,M2,F_122,G_122
*     .         ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
*       print*,f_122,g_122,c2_12
*=======================================================================
 10   FORMAT(1X,10(2X,E12.6,2X))
*
      RETURN
      END



      SUBROUTINE HZP_C0C2(M,TZ,TH,C0,C2)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 C0_DDHK,C2_DDHK
      COMPLEX*16 C0,C2
      COMPLEX*16 HZP_FTAU,HZP_GTAU
*
      C0_DDHK=-2.D0*TZ*TH/(TZ-TH)*(HZP_FTAU(TZ)-HZP_FTAU(TH))/4.D0/M**2
      C2_DDHK=
     .   ( TZ*TH/2.D0/(TZ-TH)
     .    +TZ*TH**2/2.D0/(TZ-TH)**2
     .      *(   TZ*(HZP_FTAU(TZ)-HZP_FTAU(TH))
     .        +2.D0*(HZP_GTAU(TZ)-HZP_GTAU(TH)) ) )/4.D0/M**2
*
      C0= C0_DDHK
      C2= C2_DDHK
*
      RETURN
      END


      COMPLEX*16 FUNCTION HZP_FTAU(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 XI
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      IF(X.LE.1.D0.AND.X.GT.0.D0) THEN
       HZP_FTAU=-0.25D0*(DLOG( (1.D0+DSQRT(1.D0-X))
     .                        /(1.D0-DSQRT(1.D0-X)) ) -XI*PI)**2
      ELSEIF(X.GT.1.D0) THEN
       HZP_FTAU=DASIN(DSQRT(1.D0/X))**2
      ELSE
        RETURN
      ENDIF
*
      RETURN
      END


      COMPLEX*16 FUNCTION HZP_GTAU(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 XI
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      IF(X.LE.1.D0.AND.X.GT.0.D0) THEN
       HZP_GTAU=0.5D0*DSQRT(1.D0-X)
     .               *(DLOG( (1.D0+DSQRT(1.D0-X))
     .                      /(1.D0-DSQRT(1.D0-X)) ) -XI*PI)
      ELSEIF(X.GT.1.D0) THEN
       HZP_GTAU=DSQRT(X-1.D0)*DASIN(DSQRT(1.D0/X))
      ELSE
        RETURN
      ENDIF
*
      RETURN
      END

      SUBROUTINE HZP_DILOG(Z_IN,CRES_OUT)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON /HZP_DILOG_BODE/ Z_MAG,Z_RE,Z_IM
      EXTERNAL HZP_DILOG_RE,HZP_DILOG_IM
*
      COMPLEX*16 Z_IN,CRES_OUT
      COMPLEX*16 XI,CLOG_Z_IN
      COMPLEX*16 Z_C,CRES,CF_DIFF1,CF_DIFF2
*
      COMPLEX*16 HZP_CLOG1,HZP_CLOG2
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      IF(CDABS(Z_IN).GT.1.D0) THEN
       Z_C=1.D0/Z_IN
      ELSE
       Z_C=Z_IN
      ENDIF
*
      Z_MAG=CDABS(Z_C)
      Z_RE =DREAL(Z_C)
      Z_IM =DIMAG(Z_C) 
*      print*,'>> HZP_DILOG:',z_in,' ==> ',z_mag,z_re,z_im
*
      IF(Z_MAG.EQ.0.D0) THEN
       CRES_OUT=DCMPLX(0.D0,0.D0)
       RETURN
      ENDIF
*
      EPS=1.0D-6
      NSTEP=500
      CALL BODE(HZP_DILOG_RE,EPS,1.D0-EPS,NSTEP,RES_RE)
      CALL BODE(HZP_DILOG_IM,EPS,1.D0-EPS,NSTEP,RES_IM)
*
      CRES=DCMPLX(RES_RE,RES_IM)
*
      IF(CDABS(Z_IN).GT.1.D0) THEN
        CLOG_Z_IN=HZP_CLOG2(Z_IN)
        CF_DIFF1=PI**2/3.D0-(CLOG_Z_IN)**2/2.D0-XI*PI*CLOG_Z_IN
        CLOG_Z_IN=HZP_CLOG2(-Z_IN)
        CF_DIFF2=-PI**2/6.D0-(CLOG_Z_IN)**2/2.D0
*        print*,z_in,CF_DIFF1,CF_DIFF2
*
*       IF(CDABS(CF_DIFF1).LT.CDABS(CF_DIFF2)) THEN
*        CRES_OUT=-CRES+CF_DIFF1
*       ELSE
*        CRES_OUT=-CRES+CF_DIFF2
*       ENDIF
*
*CF_DIFF2 seems to be a correct choice in our convention
       CRES_OUT=-CRES+CF_DIFF2
      ELSE
       CRES_OUT=CRES
      ENDIF
*
      RETURN
      END

      REAL*8 FUNCTION HZP_DILOG_RE(T)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HZP_DILOG_BODE/ Z_MAG,Z_RE,Z_IM
*
      COMPLEX*16 CZ,HZP_CLOG1,HZP_CLOG2
*
      CZ=1.D0-DCMPLX(Z_RE,Z_IM)*T
      HZP_DILOG_RE=-DREAL(HZP_CLOG2(CZ))/T
*
      RETURN
      END

      REAL*8 FUNCTION HZP_DILOG_IM(T)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HZP_DILOG_BODE/ Z_MAG,Z_RE,Z_IM
*
      COMPLEX*16 CZ,HZP_CLOG1,HZP_CLOG2
*
      CZ=1.D0-DCMPLX(Z_RE,Z_IM)*T
      HZP_DILOG_IM=-DIMAG(HZP_CLOG2(CZ))/T
*
      RETURN
      END

      
      SUBROUTINE HZP_FGC(MU,MZ,MH,M1,M2,F_122,G_122
     .         ,C0_12,C0_21,C1_12,C1_21,C2_12,C2_21)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 B0H_12,B0Z_12,B0H_21,B0Z_21
      COMPLEX*16 C0_12,C0_21
      COMPLEX*16 C1_12,C1_21
      COMPLEX*16 C2_12,C2_21
      COMPLEX*16 F_122,G_122
*
      COMPLEX*16 XI
      COMPLEX*16 CLAM_H,CLAM_Z
      COMPLEX*16 CARG_HP,CARG_HM,CARG_ZP,CARG_ZM
      COMPLEX*16 CRES_HP,CRES_HM,CRES_ZP,CRES_ZM
*
      COMPLEX*16 HZP_CLOG2
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      XLAM_H=(MH**4+M1**4+M2**4-2.D0*MH**2*M1**2
     .       -2.D0*MH**2*M2**2-2.D0*M1**2*M2**2)
      IF(XLAM_H.GT.0.D0) THEN
       CLAM_H=DCMPLX(DSQRT(XLAM_H),0.D0)
      ELSE
       CLAM_H=DCMPLX(0.D0,DSQRT(-XLAM_H))
      ENDIF
*
      XLAM_Z=(MZ**4+M1**4+M2**4-2.D0*MZ**2*M1**2
     .       -2.D0*MZ**2*M2**2-2.D0*M1**2*M2**2)
      IF(XLAM_Z.GT.0.D0) THEN
       CLAM_Z=DCMPLX(DSQRT(XLAM_Z),0.D0)
      ELSE
       CLAM_Z=DCMPLX(0.D0,DSQRT(-XLAM_Z))
      ENDIF
*
*One should not ignore the i\epsilon terms from the loop masses
      EPS_POLE=-1.D-10
      CLAM_H=CDSQRT(MH**4
     .            +(M1**2+XI*EPS_POLE)**2
     .            +(M2**2+XI*EPS_POLE)**2
     .            -2.D0*MH**2*(M1**2+XI*EPS_POLE)
     .            -2.D0*MH**2*(M2**2+XI*EPS_POLE)
     .            -2.D0*(M1**2+XI*EPS_POLE)*(M2**2+XI*EPS_POLE))
      CLAM_Z=CDSQRT(MZ**4
     .            +(M1**2+XI*EPS_POLE)**2
     .            +(M2**2+XI*EPS_POLE)**2
     .            -2.D0*MZ**2*(M1**2+XI*EPS_POLE)
     .            -2.D0*MZ**2*(M2**2+XI*EPS_POLE)
     .            -2.D0*(M1**2+XI*EPS_POLE)*(M2**2+XI*EPS_POLE))
*      print*,clam_h,clam_z
*
      CARG_HP=2.D0*MH**2/(M2**2-M1**2+MH**2+CLAM_H)
      CARG_HM=2.D0*MH**2/(M2**2-M1**2+MH**2-CLAM_H)
      CARG_ZP=2.D0*MZ**2/(M2**2-M1**2+MZ**2+CLAM_Z)
      CARG_ZM=2.D0*MZ**2/(M2**2-M1**2+MZ**2-CLAM_Z)
      CALL HZP_DILOG(CARG_HP,CRES_HP)
      CALL HZP_DILOG(CARG_HM,CRES_HM)
      CALL HZP_DILOG(CARG_ZP,CRES_ZP)
      CALL HZP_DILOG(CARG_ZM,CRES_ZM)
      C0_12=(-CRES_HP-CRES_HM+CRES_ZP+CRES_ZM)/(MH**2-MZ**2)
*      print*,'C0(m1,m2,m2)=',c0_12
*      print*,carg_hp,cres_hp
*      print*,carg_hm,cres_hm
*      print*,carg_zp,cres_zp
*      print*,carg_zm,cres_zm
*
      CARG_HP=2.D0*MH**2/(M1**2-M2**2+MH**2+CLAM_H)
      CARG_HM=2.D0*MH**2/(M1**2-M2**2+MH**2-CLAM_H)
      CARG_ZP=2.D0*MZ**2/(M1**2-M2**2+MZ**2+CLAM_Z)
      CARG_ZM=2.D0*MZ**2/(M1**2-M2**2+MZ**2-CLAM_Z)
      CALL HZP_DILOG(CARG_HP,CRES_HP)
      CALL HZP_DILOG(CARG_HM,CRES_HM)
      CALL HZP_DILOG(CARG_ZP,CRES_ZP)
      CALL HZP_DILOG(CARG_ZM,CRES_ZM)
      C0_21=(-CRES_HP-CRES_HM+CRES_ZP+CRES_ZM)/(MH**2-MZ**2)
*      print*,'C0(m2,m1,m1)=',c0_21
*      print*,carg_hp,cres_hp
*      print*,carg_hm,cres_hm
*      print*,carg_zp,cres_zp
*      print*,carg_zm,cres_zm
*
*--> test of ComplexLOG function
*
*      print*,hzp_clog(dcmplx(1.d0,0.d0))
*      print*,hzp_clog(dcmplx(0.d0,1.d0))
*      print*,hzp_clog(dcmplx(1.d0,1.d0)),dlog(2.d0)/2.d0,pi/4.d0
*
      B0H_12= 2.D0-DLOG(M1*M2/MU**2)
     .      +(M1**2-M2**2)/MH**2*DLOG(M2/M1)
     .      +CLAM_H/MH**2
     .       *HZP_CLOG2((M1**2+M2**2-MH**2+CLAM_H)/2.D0/M1/M2)
      B0Z_12= 2.D0-DLOG(M1*M2/MU**2)
     .      +(M1**2-M2**2)/MZ**2*DLOG(M2/M1)
     .      +CLAM_Z/MZ**2
     .       *HZP_CLOG2((M1**2+M2**2-MZ**2+CLAM_Z)/2.D0/M1/M2)
      C1_12=(B0H_12-B0Z_12)/(MZ**2-MH**2)-C0_12
*      print*,'C1(m1,m2,m2)=',c1_12
*
      B0H_21= 2.D0-DLOG(M2*M1/MU**2)
     .      +(M2**2-M1**2)/MH**2*DLOG(M1/M2)
     .      +CLAM_H/MH**2
     .       *HZP_CLOG2((M2**2+M1**2-MH**2+CLAM_H)/2.D0/M2/M1)
      B0Z_21= 2.D0-DLOG(M2*M1/MU**2)
     .      +(M2**2-M1**2)/MZ**2*DLOG(M1/M2)
     .      +CLAM_Z/MZ**2
     .       *HZP_CLOG2((M2**2+M1**2-MZ**2+CLAM_Z)/2.D0/M2/M1)
      C1_21=(B0H_21-B0Z_21)/(MZ**2-MH**2)-C0_21
*      print*,'C2(m2,m1,m1)=',c1_21
*
      A0_1=M1**2*(1.D0-DLOG(M1**2/MU**2))
      A0_2=M2**2*(1.D0-DLOG(M2**2/MU**2))
*
      C2_12=(M1**2-M2**2-MZ**2)/2.D0/(MZ**2-MH**2)**2
     .     *(B0H_12-B0Z_12)
     .     +1.D0/2.D0/(MZ**2-MH**2)/MH**2
     .     *(MH**2+2.D0*M2**2*MH**2*C0_12
     .       +(M2**2-M1**2)*B0H_12+A0_1-A0_2)
*      print*,'C2(m1,m2,m2)=',c2_12
*
      C2_21=(M2**2-M1**2-MZ**2)/2.D0/(MZ**2-MH**2)**2
     .     *(B0H_21-B0Z_21)
     .     +1.D0/2.D0/(MZ**2-MH**2)/MH**2
     .     *(MH**2+2.D0*M1**2*MH**2*C0_21
     .       +(M1**2-M2**2)*B0H_21+A0_2-A0_1)
*      print*,'C2(m2,m1,m1)=',c2_21
*
      F_122=-2.D0*(C0_12+C1_12+2.D0*C2_12+2.D0*C2_21-C1_21)
      G_122=-2.D0*(C0_12+C1_12+C1_21)
*
*      print*,'>> HZP_FGC:',f_122,g_122,c2_12
*
      RETURN
      END


      COMPLEX*16 FUNCTION HZP_CLOG1(CX)
***********************************************************************
*
* Arg(CX)=PHI=[0..2pi]
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 CX
      COMPLEX*16 XI
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      IF(DREAL(CX).EQ.0.D0) THEN
       IF(DIMAG(CX).GT.0.D0) PHI=PI/2.D0
       IF(DIMAG(CX).LT.0.D0) PHI=3.D0*PI/2.D0
       GOTO 99
      ENDIF
*
      IF(DIMAG(CX).EQ.0.D0) THEN
       IF(DREAL(CX).GT.0.D0) PHI=0.D0
       IF(DREAL(CX).LT.0.D0) PHI=PI
       GOTO 99
      ENDIF
*
      PHI=DATAN(DIMAG(CX)/DREAL(CX))
      IF(DREAL(CX).GT.0.D0.AND.DIMAG(CX).GT.0.D0) PHI=PHI
      IF(DREAL(CX).LT.0.D0.AND.DIMAG(CX).GT.0.D0) PHI=PHI+PI
      IF(DREAL(CX).LT.0.D0.AND.DIMAG(CX).LT.0.D0) PHI=PHI+PI
      IF(DREAL(CX).GT.0.D0.AND.DIMAG(CX).LT.0.D0) PHI=PHI+2.D0*PI
      IF(PHI.LT.0.D0) print*,'Error in HZP_CLOG1: Arg(Z) <0 ',phi
*
 99   CONTINUE
      HZP_CLOG1=DLOG(CDABS(CX))+XI*PHI
*
      RETURN
      END

      COMPLEX*16 FUNCTION HZP_CLOG2(CX)
***********************************************************************
*
* Arg(CX)=PHI=[-pi..pi]
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 CX
      COMPLEX*16 XI
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      IF(DREAL(CX).EQ.0.D0) THEN
       IF(DIMAG(CX).GT.0.D0) PHI=PI/2.D0
       IF(DIMAG(CX).LT.0.D0) PHI=-PI/2.D0
       GOTO 99
      ENDIF
*
      IF(DIMAG(CX).EQ.0.D0) THEN
       IF(DREAL(CX).GT.0.D0) PHI=0.D0
       IF(DREAL(CX).LT.0.D0) PHI=PI
       GOTO 99
      ENDIF
*
      PHI=DATAN(DIMAG(CX)/DREAL(CX))
      IF(DREAL(CX).GT.0.D0.AND.DIMAG(CX).GT.0.D0) PHI=PHI
      IF(DREAL(CX).LT.0.D0.AND.DIMAG(CX).GT.0.D0) PHI=PHI+PI
      IF(DREAL(CX).LT.0.D0.AND.DIMAG(CX).LT.0.D0) PHI=PHI-PI
      IF(DREAL(CX).GT.0.D0.AND.DIMAG(CX).LT.0.D0) PHI=PHI
      IF(DABS(PHI).GT.PI) print*
     .                  ,'Error in HZP_CLOG2: |Arg(Z)| > pi ',phi
*
 99   CONTINUE
      HZP_CLOG2=DLOG(CDABS(CX))+XI*PHI
*
      RETURN
      END



      SUBROUTINE UPDATEGAMBR_0(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC,HZP_WIDTH,HZP_WIDTH_SM)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
*ARRAYS:
      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
*
      INTEGER*8 IFLAG_H(NFLAG)
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)  ! 100 = NCMAX
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)  
      REAL*8 GAMBRC(NMCH,3)    
*-----------------------------------------------------------------------
*Local Vairables:
*
      INTEGER*8 ISMN,ISUSYN,ISMC,ISUSYC
      REAL*8    HZP_WIDTH(3),HZP_WIDTH_SM(3)
*=======================================================================
*      do iwb=1,3
*      print*,'INITL',iwb,gambrn(19,iwb,1)
*      print*,'INITL',iwb,gambrn(19,iwb,2)
*      print*,'INITL',iwb,gambrn(19,iwb,3)
*      enddo
*
      GAMBRN(19,1,1)=HZP_WIDTH(1)
      GAMBRN(19,1,2)=HZP_WIDTH(2)
      GAMBRN(19,1,3)=HZP_WIDTH(3)
*=======================================================================
      ISMN   = IFLAG_H(20)
      ISUSYN = IFLAG_H(21)
      ISMC   = IFLAG_H(22)
      ISUSYC = IFLAG_H(23)
*      print*,ismn,isusyn
*
* << BRANCHING RATIOS OF NEUTRAL HIGGS BOSON >>
*
      DO IH=1,3
*
       GAMBRN(ISMN,1,IH)=0.D0
       DO IM=1,ISMN-1
       GAMBRN(ISMN,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       IXX=ISMN+ISUSYN
       GAMBRN(IXX,1,IH)=0.D0
       DO IM=ISMN+1,IXX-1
       GAMBRN(IXX,1,IH)=GAMBRN(IXX,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       GAMBRN(NMNH,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IXX,1,IH)
*
       DO IM=1,ISMN+ISUSYN+1
        GAMBRN(IM,2,IH)=GAMBRN(IM,1,IH)/GAMBRN(ISMN,1,IH)
        GAMBRN(IM,3,IH)=GAMBRN(IM,1,IH)/GAMBRN(NMNH,1,IH)
       ENDDO
*
*Refill GAMBRN(IM,2,IH) with the corresponding SM BRs
*
      CALL FILLSMBRS(HMASS_H(IH),GAMBRN,NMNH)
      DO IM=1,ISMN
       GAMBRN(IM,2,IH)=RAUX_H(650+IM)
      ENDDO
      DO IM=ISMN+1,ISMN+ISUSYN
       GAMBRN(IM,2,IH)=0.D0
      ENDDO
       GAMBRN(ISMN+ISUSYN+1,2,IH)=1.D0
*
      ENDDO ! IH
*
*      do iwb=1,3
*      print*,'FINAL',iwb,gambrn(19,iwb,1)
*      print*,'FINAL',iwb,gambrn(19,iwb,2)
*      print*,'FINAL',iwb,gambrn(19,iwb,3)
*      enddo
*
*
* << BRANCHING RATIOS OF CHARGED HIGGS BOSON >>
*
*
*Refill GAMBRC(IM,2) with the corresponding SM BRs
*
       DO IM=1,ISMC+ISUSYC+1
        GAMBRC(IM,2)=0.D0
       ENDDO
*
*      print*,'>>> UPDATEGAMBR ... done',ismn,isusyn
*
      IF(IFLAG_H(6).EQ.1) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,1)
      IF(IFLAG_H(6).EQ.2) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,2)
      IF(IFLAG_H(6).EQ.3) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,3)
      IF(IFLAG_H(6).EQ.5) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,5)
*
      IF(IFLAG_H(6).EQ.4.OR.IFLAG_H(6).EQ.5) 
     . CALL DUMP_CHDCY(ISMC,ISUSYC,NMCH,GAMBRC)
*=======================================================================
*
      RETURN
      END
