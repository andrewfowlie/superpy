      SUBROUTINE FILLMUON(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     . ,MCH,HMASS,OMIX
     . ,STMASS,STMIX,SBMASS,SBMIX,STAUMASS,STAUMIX,SNU3MASS
     . ,M_C,C_L,C_R,M_N,N_N,NCMAX,NHC,SHC,CHC)
************************************************************************
*
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
*Input arrays
      REAL*8 SMPARA(NSMIN),SSPARA(NSSIN)
      INTEGER*8 IFLAG(NFLAG)
*
      REAL*8     MCH,HMASS(3),OMIX(3,3)
*
      REAL*8     STMASS(2),SBMASS(2),STAUMASS(2),SNU3MASS
      COMPLEX*16 STMIX(2,2),SBMIX(2,2),STAUMIX(2,2)
*
      REAL*8     M_C(2),M_N(4)
      COMPLEX*16 C_L(2,2),C_R(2,2),N_N(4,4)
*
      COMPLEX*16 NHC(NCMAX,3)
      REAL*8     SHC(NCMAX)
      COMPLEX*16 CHC(NCMAX)
*
      COMPLEX*16 A_MU
*-----------------------------------------------------------------------
*Local arrays
      COMPLEX*16 XK_MU
      COMPLEX*16 H_MU
      COMPLEX*16 XRL,XLR
      REAL*8     SMUMASS(2)
      COMPLEX*16 SMUMIX(2,2)
*Couplings at Mt^pole
      COMPLEX*16 CGL,CGR,H_F,SFMIX(2,2)
*     chargino(I)-muon-senutrino_mu
      COMPLEX*16 GL_CI_MU_SNMU(2),  GR_CI_MU_SNMU(2)
*     neutralino(I)-muon-smuon(J)
      COMPLEX*16 GL_NI_MU_SMUJ(4,2),GR_NI_MU_SMUJ(4,2)
*Electric EDMs of muon: [d^E_mu/e]
      REAL*8     DEOE_MU(20)
      DATA       NEDM_SUB/20/
      REAL*8     MDM_MU(20)
      DATA       NMDM_SUB/20/
*-----------------------------------------------------------------------
      EXTERNAL EDM_A,EDM_B
      EXTERNAL MDM_A,MDM_B
* For integration
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
      EXTERNAL F0_HiggsEDM,F_HiggsEDM,G_HiggsEDM
*
      NX_EDM  = 1000   ! Number of calling SUBROUTINE BODE
      EPS_EDM = 1.D-6  ! integration region of 2 loop functions:
*                        [EPS_EDM ; (1-EPS_EDM)]
      NX  = NX_EDM
      EPS = EPS_EDM
*
      X1D=EPS
      X1U=1.D0-EPS
*-----------------------------------------------------------------------
      PI=2.D0*DASIN(1.D0)
      AEM=1.D0/137.D0
      GEVTOCM=1.97326968D-14
*
      A_MU=SSPARA(29)*DCMPLX(DCOS(SSPARA(30)/180.D0*PI)
     .                      ,DSIN(SSPARA(30)/180.D0*PI))
*      print*,'A_MU =',a_mu
*
      I_DEOE_MU    =360 ! d^E_mu/e [cm]
      I_MDM_MU     =380 ! a_mu
*-----------------------------------------------------------------------
*      print*,a_mu,atau_h
      R_123Q=SSPARA(22)
      R_123U=SSPARA(23)
      R_123D=SSPARA(24)
      R_123L=SSPARA(25)
      R_123E=SSPARA(26)
*
      MQ1=R_123Q*MQ3_H
      MQ2=R_123Q*MQ3_H
      MQ3=       MQ3_H
      MU1=R_123U*MU3_H
      MU2=R_123U*MU3_H
      MU3=       MU3_H
      MD1=R_123D*MD3_H
      MD2=R_123D*MD3_H
      MD3=       MD3_H
      ML1=R_123L*ML3_H
      ML2=R_123L*ML3_H
      ML3=       ML3_H
      ME1=R_123E*ME3_H
      ME2=R_123E*ME3_H
      ME3=       ME3_H
*      print*,mq1,mq2,mq3
*      print*,mu1,mu2,mu3
*      print*,md1,md2,md3
*      print*,ml1,ml2,ml3
*      print*,me1,me2,me3
*-----------------------------------------------------------------------
*Smuon mixing
*
      H_MU0=DCMPLX(DSQRT(2.D0)*MMU_H/V_H/CB_H,0.D0)
      XLL  = ML2**2+DABS(MMU_H)**2
     .      +(CB_H**2-SB_H**2)*MZ_H**2*(SW_H**2-1.D0/2.D0)
      XRR  = ME2**2+DABS(MMU_H)**2
     .      -(CB_H**2-SB_H**2)*MZ_H**2*SW_H**2
      XRL  = H_MU0*V_H*CB_H*(A_MU-DCONJG(MU_H)*SB_H/CB_H)/SQRT(2.D0)
      XLR  = DCONJG(XRL)

      PHI=DATAN(DIMAG(XRL)/DREAL(XRL))
      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SMUMASS,SMUMIX)
      IF(SMUMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative smuon mass squared'
        RETURN
      ENDIF

*EW threshold corrections: see Ellis, Lee, Pilaftsis, hep-ph/0404167
      SNU2MASS=DSQRT(ML2**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0)
      IF(IFLAG(10).EQ.0) THEN
       XK_MU=DCONJG(MU_H)*AEM_H/4.D0/PI*
     . (
     .  -DCONJG(M2_H)/SW_H**2*
     .   (F_I(SNU2MASS**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .   +F_I(SMUMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(1,1))**2/2.D0
     .   +F_I(SMUMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(1,2))**2/2.D0)  ! M2_H
     .  +DCONJG(M1_H)/CW_H**2*
     .   (F_I(SMUMASS(1)**2,SMUMASS(2)**2,CDABS(M1_H)**2)
     .   +F_I(SMUMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(1,1))**2/2.D0
     .   +F_I(SMUMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(1,2))**2/2.D0
     .   -F_I(SMUMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(2,1))**2
     .   -F_I(SMUMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SMUMIX(2,2))**2)       ! M1_H
     . )
       H_MU=DCMPLX(DSQRT(2.D0)*MMU_H/V_H/CB_H,0.D0)/(1.D0+XK_MU*TB_H)
      ELSE
       H_MU=DCMPLX(DSQRT(2.D0)*MMU_H/V_H/CB_H,0.D0)
      ENDIF
*      print*,h_mu0,h_mu,xk_mu

*      print*,'Smuon Sector:'
*      WRITE(*,1) SMUMASS(1),SMUMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SMUMIX(1,1)),DIMAG(SMUMIX(1,1))
*     .          ,DREAL(SMUMIX(1,2)),DIMAG(SMUMIX(1,2))
*      WRITE(*,4) DREAL(SMUMIX(2,1)),DIMAG(SMUMIX(2,1))
*     .          ,DREAL(SMUMIX(2,2)),DIMAG(SMUMIX(2,2))
*      print*,' '
*-----------------------------------------------------------------------
*Couplings:
*     chargino(I)-muon-snutrino_mu
*      COMPLEX*16 GL_CI_MU_SNMU(2),  GR_CI_MU_SNMU(2)
*      print*,'chargino(I)-muon-sneutrino_mu:'
      DO I=1,2
       GL_CI_MU_SNMU(I)=-GW_H*C_R(I,1)
       GR_CI_MU_SNMU(I)=DCONJG(H_MU)*C_L(I,2)
       CGL=GL_CI_MU_SNMU(I)
       CGR=GR_CI_MU_SNMU(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
*     neutralino(I)-muon-smuon(J)
*      COMPLEX*16 GL_NI_MU_SMUJ(4,2),GR_NI_MU_SMUJ(4,2)
      I_A =3
      T_3F=-1.0/2.D0
      Q_F =-1.D0
      H_F =H_MU
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SMUMIX(I,J)
       ENDDO
      ENDDO

      DO I=1,4
       DO J=1,2
        GL_NI_MU_SMUJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_N(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_N(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_N(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_MU_SMUJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_N(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_N(I,I_A)*DCONJG(SFMIX(1,J))
        CGL=GL_NI_MU_SMUJ(I,J)
        CGR=GR_NI_MU_SMUJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*Dumping couplings
      IF(IFLAG(19).EQ.2) THEN
      print*,'---------------------------------------------------------'
      print*,' Ino-couplings to fermion and sfermion needed for EDM :'
      print*,'                            G_R, G_L, and Im[G_R^* G_L]'
      print*,'---------------------------------------------------------'
      print*,'chargino(I)-muon-sneutrino_mu:'
      DO I=1,2
       CGL=GL_CI_MU_SNMU(I)
       CGR=GR_CI_MU_SNMU(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'neutralino(I)-muon-smuon(J):'
      DO I=1,4
       DO J=1,2
        CGL=GL_NI_MU_SMUJ(I,J)
        CGR=GR_NI_MU_SMUJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'---------------------------------------------------------'
      ENDIF ! IF(IFLAG(19).EQ.2) THEN
*-----------------------------------------------------------------------
*Electric EDMs of muon: [d^E_f/e]
*      REAL*8     DEOE_MU(20)
************************************************************************
*Convention of the array DEOE 
*      1 = total = 2+3+4+5
*                  2 = chargino contribution
*                  3 = neutralino contribution
*                  4 = gluino contribution
*                  5 = Higgs contribution = 11+12+13+14+15+16+17+18
*                                           11 = top contrubution 
*                                           12 = bottom contrubution 
*                                           13 = stop contrubution 
*                                           14 = sbottom contrubution 
*                                           15 = tau contrubution 
*                                           16 = stau contrubution 
*                                           17 = chargino contrubution 
*                                           18 = charged Higgs 
************************************************************************
*
*      print*, NEDM_SUB
*Initialize
       DO I=1,NEDM_SUB
        DEOE_MU(I)     =0.D0
       ENDDO
*
*Chargino contributions:
*=======================
      IEDM_SUB=2
*
      MSNMU=DSQRT(ML2**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0)

      DEOE_MU(IEDM_SUB)=-1.D0/16.D0/PI**2
     .*( M_C(1)/MSNMU**2
     .                *DIMAG(DCONJG(GR_CI_MU_SNMU(1))*GL_CI_MU_SNMU(1))
     .                *EDM_A(M_C(1)**2/MSNMU**2)
     .  +M_C(2)/MSNMU**2
     .                *DIMAG(DCONJG(GR_CI_MU_SNMU(2))*GL_CI_MU_SNMU(2))
     .                *EDM_A(M_C(2)**2/MSNMU**2) )

*
*Neutralino contributions:
*=========================
      IEDM_SUB=3

      Q_SF=-1.D0
      DEOE_MU(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_N(1)/SMUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(1,1))*GL_NI_MU_SMUJ(1,1))
     .             *Q_SF*EDM_B(M_N(1)**2/SMUMASS(1)**2)
     .  +M_N(1)/SMUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(1,2))*GL_NI_MU_SMUJ(1,2))
     .             *Q_SF*EDM_B(M_N(1)**2/SMUMASS(2)**2)
     .  +M_N(2)/SMUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(2,1))*GL_NI_MU_SMUJ(2,1))
     .             *Q_SF*EDM_B(M_N(2)**2/SMUMASS(1)**2)
     .  +M_N(2)/SMUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(2,2))*GL_NI_MU_SMUJ(2,2))
     .             *Q_SF*EDM_B(M_N(2)**2/SMUMASS(2)**2)
     .  +M_N(3)/SMUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(3,1))*GL_NI_MU_SMUJ(3,1))
     .             *Q_SF*EDM_B(M_N(3)**2/SMUMASS(1)**2)
     .  +M_N(3)/SMUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(3,2))*GL_NI_MU_SMUJ(3,2))
     .             *Q_SF*EDM_B(M_N(3)**2/SMUMASS(2)**2)
     .  +M_N(4)/SMUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(4,1))*GL_NI_MU_SMUJ(4,1))
     .             *Q_SF*EDM_B(M_N(4)**2/SMUMASS(1)**2)
     .  +M_N(4)/SMUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_MU_SMUJ(4,2))*GL_NI_MU_SMUJ(4,2))
     .             *Q_SF*EDM_B(M_N(4)**2/SMUMASS(2)**2) )

*
*Gluino contributions:
*=====================
      IEDM_SUB=4

      Q_SF=-1.D0
      DEOE_MU(IEDM_SUB)= 0.D0

*
*Two-loop Higgs contributions the muon EDM : d^E_mu/e
*====================================================
      IEDM_SUB=5

      IEDM_SUB_HIGGS=11 ! top
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MTMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTOP)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTOP)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*MMU_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(6,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(5,IH))*DREAL(NHC(27,IH))*GTOP)
      ENDDO

      IEDM_SUB_HIGGS=12 ! bottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MBMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FBOT)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GBOT)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*MMU_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(6,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(5,IH))*DREAL(NHC(18,IH))*GBOT)
      ENDDO

      IEDM_SUB_HIGGS=13 ! stop
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=STMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST1)
       Z_HiggsEDM=STMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST2)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MMU_H/32.D0/PI**3*DREAL(NHC(6,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
      ENDDO

      IEDM_SUB_HIGGS=14 ! sbottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=SBMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB1)
       Z_HiggsEDM=SBMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB2)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MMU_H/32.D0/PI**3*DREAL(NHC(6,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
      ENDDO

      IEDM_SUB_HIGGS=15 ! tau
      DO IH=1,3
       Z_HiggsEDM=MTAU_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTAU)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTAU)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . -AEM**2*MMU_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(6,IH))*DREAL(NHC(8,IH))*FTAU
     .   +DREAL(NHC(5,IH))*DREAL(NHC(9,IH))*GTAU)
      ENDDO

      IEDM_SUB_HIGGS=16 ! stau
      DO IH=1,3
       Z_HiggsEDM=STAUMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU1)
       Z_HiggsEDM=STAUMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU2)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . +AEM*MMU_H/32.D0/PI**3*DREAL(NHC(6,IH))/HMASS(IH)**2
     .  *DREAL(NHC(79,IH)*FSTAU1+NHC(82,IH)*FSTAU2)
      ENDDO

      IEDM_SUB_HIGGS=17 ! chargino
      DO IH=1,3
       Z_HiggsEDM=M_C(1)**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FC1)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GC1)
       Z_HiggsEDM=M_C(2)**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FC2)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GC2)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . -AEM**2*MMU_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(DREAL(NHC(6,IH))*DREAL(NHC(59,IH))*FC1
     .   +DREAL(NHC(5,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*MMU_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(DREAL(NHC(6,IH))*DREAL(NHC(68,IH))*FC2
     .   +DREAL(NHC(5,IH))*DREAL(NHC(69,IH))*GC2)
      ENDDO

      IEDM_SUB_HIGGS=18 ! charged Higgs
      DO IH=1,3
       Z_HiggsEDM=MCH**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FCH)
       DEOE_MU(IEDM_SUB_HIGGS)=DEOE_MU(IEDM_SUB_HIGGS)
     . +AEM*MMU_H/32.D0/PI**3*DREAL(NHC(6,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
      ENDDO

*Flip the signs of the fermionic contributions to Barr-Zee EDM
      DEOE_MU(11)=-DEOE_MU(11)
      DEOE_MU(12)=-DEOE_MU(12)
      DEOE_MU(15)=-DEOE_MU(15)
      DEOE_MU(17)=-DEOE_MU(17)

*      do i=11,18
*       print*,i,deoe_mu(i)
*      enddo
      DEOE_MU(IEDM_SUB)=DEOE_MU(11)+DEOE_MU(12)+DEOE_MU(13)+DEOE_MU(14)
     .                 +DEOE_MU(15)+DEOE_MU(16) 
     .                 +DEOE_MU(17)+DEOE_MU(18)
*
*TOTAL:
*======
      IEDM_SUB=1
      DEOE_MU(IEDM_SUB)=DEOE_MU(2)+DEOE_MU(3)+DEOE_MU(4)+DEOE_MU(5)
*
*STORE:Muon Elecron EDM:
      RAUX_H(I_DEOE_MU+0)=DEOE_MU(1)*GEVTOCM
      RAUX_H(I_DEOE_MU+1)=DEOE_MU(2)*GEVTOCM
      RAUX_H(I_DEOE_MU+2)=DEOE_MU(3)*GEVTOCM
      RAUX_H(I_DEOE_MU+3)=DEOE_MU(4)*GEVTOCM
      RAUX_H(I_DEOE_MU+4)=DEOE_MU(5)*GEVTOCM
*
      IF(IFLAG(19).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'The Electric EDM of muon in units of [cm]: (d^mu/e)/[cm]'
      print*,'---------------------------------------------------------'
      write(*,8) '(d^E_mu/e)/cm[Total]:',deoe_mu(1)*GEVTOCM
      write(*,8) '(d^E_mu/e)/cm[C,N,Gl,H]:'
     .      ,deoe_mu(2)*GEVTOCM,deoe_mu(3)*GEVTOCM
     .      ,deoe_mu(4)*GEVTOCM,deoe_mu(5)*GEVTOCM
      ENDIF ! IF(IFLAG(19).EQ.1) THEN
*-----------------------------------------------------------------------
*MDM of muon: [a_mu]
*      REAL*8     MDM_MU(20)
************************************************************************
*Convention of the array MDM
*      1 = total = 2+3+4+5
*                  2 = chargino contribution
*                  3 = neutralino contribution
*                  4 = gluino contribution
*                  5 = Higgs contribution = 11+12+13+14+15+16+17+18
*                                           11 = top contrubution 
*                                           12 = bottom contrubution 
*                                           13 = stop contrubution 
*                                           14 = sbottom contrubution 
*                                           15 = tau contrubution 
*                                           16 = stau contrubution 
*                                           17 = chargino contrubution 
*                                           18 = charged Higgs 
************************************************************************
*      print*, NMDM_SUB
*Initialize
       DO I=1,NMDM_SUB
        MDM_MU(I)     =0.D0
       ENDDO
*
*Chargino contributions:
*=======================
      IMDM_SUB=2
*
      MSNMU=DSQRT(ML2**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0)
*      print*,msnmu

      MDM_MU(IMDM_SUB)=MMU_H**2/8.D0/PI**2/MSNMU**2
     .*( (CDABS(GR_CI_MU_SNMU(1))**2+CDABS(GL_CI_MU_SNMU(1))**2)
     .   *MDM_A(M_C(1)**2/MSNMU**2)
     .  +(CDABS(GR_CI_MU_SNMU(2))**2+CDABS(GL_CI_MU_SNMU(2))**2)
     .   *MDM_A(M_C(2)**2/MSNMU**2)
     .  -M_C(1)/MMU_H*DREAL(DCONJG(GR_CI_MU_SNMU(1))*GL_CI_MU_SNMU(1))
     .               *EDM_A(M_C(1)**2/MSNMU**2)
     .  -M_C(2)/MMU_H*DREAL(DCONJG(GR_CI_MU_SNMU(2))*GL_CI_MU_SNMU(2))
     .               *EDM_A(M_C(2)**2/MSNMU**2)  
     . )
*      print*,mdm_mu(imdm_sub)
*
*Neutralino contributions:
*=======================
      IMDM_SUB=3

      MDM_MU(IMDM_SUB)=MMU_H**2/8.D0/PI**2/SMUMASS(1)**2
     .*(-(CDABS(GR_NI_MU_SMUJ(1,1))**2+CDABS(GL_NI_MU_SMUJ(1,1))**2)
     .   *MDM_B(M_N(1)**2/SMUMASS(1)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(2,1))**2+CDABS(GL_NI_MU_SMUJ(2,1))**2)
     .   *MDM_B(M_N(2)**2/SMUMASS(1)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(3,1))**2+CDABS(GL_NI_MU_SMUJ(3,1))**2)
     .   *MDM_B(M_N(3)**2/SMUMASS(1)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(4,1))**2+CDABS(GL_NI_MU_SMUJ(4,1))**2)
     .   *MDM_B(M_N(4)**2/SMUMASS(1)**2)
     .  -M_N(1)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(1,1))
     .                            *GL_NI_MU_SMUJ(1,1))
     .               *EDM_B(M_N(1)**2/SMUMASS(1)**2)
     .  -M_N(2)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(2,1))
     .                            *GL_NI_MU_SMUJ(2,1))
     .               *EDM_B(M_N(2)**2/SMUMASS(1)**2)  
     .  -M_N(3)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(3,1))
     .                            *GL_NI_MU_SMUJ(3,1))
     .               *EDM_B(M_N(3)**2/SMUMASS(1)**2)  
     .  -M_N(4)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(4,1))
     .                            *GL_NI_MU_SMUJ(4,1))
     .               *EDM_B(M_N(4)**2/SMUMASS(1)**2)  
     . )
     .                +MMU_H**2/8.D0/PI**2/SMUMASS(2)**2
     .*(-(CDABS(GR_NI_MU_SMUJ(1,2))**2+CDABS(GL_NI_MU_SMUJ(1,2))**2)
     .   *MDM_B(M_N(1)**2/SMUMASS(2)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(2,2))**2+CDABS(GL_NI_MU_SMUJ(2,2))**2)
     .   *MDM_B(M_N(2)**2/SMUMASS(2)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(3,2))**2+CDABS(GL_NI_MU_SMUJ(3,2))**2)
     .   *MDM_B(M_N(3)**2/SMUMASS(2)**2)
     .  -(CDABS(GR_NI_MU_SMUJ(4,2))**2+CDABS(GL_NI_MU_SMUJ(4,2))**2)
     .   *MDM_B(M_N(4)**2/SMUMASS(2)**2)
     .  -M_N(1)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(1,2))
     .                            *GL_NI_MU_SMUJ(1,2))
     .               *EDM_B(M_N(1)**2/SMUMASS(2)**2)
     .  -M_N(2)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(2,2))
     .                            *GL_NI_MU_SMUJ(2,2))
     .               *EDM_B(M_N(2)**2/SMUMASS(2)**2)
     .  -M_N(3)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(3,2))
     .                            *GL_NI_MU_SMUJ(3,2))
     .               *EDM_B(M_N(3)**2/SMUMASS(2)**2)
     .  -M_N(4)/MMU_H*DREAL(DCONJG(GR_NI_MU_SMUJ(4,2))
     .                            *GL_NI_MU_SMUJ(4,2))
     .               *EDM_B(M_N(4)**2/SMUMASS(2)**2)
     . )
*      print*,mdm_mu(imdm_sub)
*
*Gluino contributions:
*=====================
      IMDM_SUB=4
      MDM_MU(IMDM_SUB)=0.D0
*      print*,mdm_mu(imdm_sub)

*Two-loop Higgs contributions to the muon MDM : a_mu
*====================================================
      IMDM_SUB=5

      IMDM_SUB_HIGGS=11 ! top
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MTMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTOP)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTOP)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*MMU_H**2/4.D0/PI**2/SW_H**2/MW_H**2
     .   *(-DREAL(NHC(5,IH))*DREAL(NHC(26,IH))*FTOP
     .     +DREAL(NHC(6,IH))*DREAL(NHC(27,IH))*GTOP)
      ENDDO

      IMDM_SUB_HIGGS=12 ! bottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MBMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FBOT)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GBOT)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*MMU_H**2/4.D0/PI**2/SW_H**2/MW_H**2
     .  *(-DREAL(NHC(5,IH))*DREAL(NHC(17,IH))*FBOT
     .    +DREAL(NHC(6,IH))*DREAL(NHC(18,IH))*GBOT)
      ENDDO

      IMDM_SUB_HIGGS=13 ! stop
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=STMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST1)
       Z_HiggsEDM=STMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST2)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     .-3.D0*AEM*QQ**2*MMU_H**2/16.D0/PI**3*DREAL(NHC(5,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
      ENDDO

      IMDM_SUB_HIGGS=14 ! sbottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=SBMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB1)
       Z_HiggsEDM=SBMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB2)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     .-3.D0*AEM*QQ**2*MMU_H**2/16.D0/PI**3*DREAL(NHC(5,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
      ENDDO

      IMDM_SUB_HIGGS=15 ! tau
      DO IH=1,3
       Z_HiggsEDM=MTAU_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTAU)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTAU)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     . -AEM**2*MMU_H**2/4.D0/PI**2/SW_H**2/MW_H**2
     .  *(-DREAL(NHC(5,IH))*DREAL(NHC(8,IH))*FTAU
     .    +DREAL(NHC(6,IH))*DREAL(NHC(9,IH))*GTAU)
*          print*,ih,DREAL(NHC(5,IH))*DREAL(NHC(8,IH)),FTAU
*     .             ,DREAL(NHC(6,IH))*DREAL(NHC(9,IH)),GTAU
*     .          ,MDM_MU(IMDM_SUB_HIGGS)
      ENDDO

      IMDM_SUB_HIGGS=16 ! stau
      DO IH=1,3
       Z_HiggsEDM=STAUMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU1)
       Z_HiggsEDM=STAUMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU2)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     . -AEM*MMU_H**2/16.D0/PI**3*DREAL(NHC(5,IH))/HMASS(IH)**2
     .  *DREAL(NHC(79,IH)*FSTAU1+NHC(82,IH)*FSTAU2)
      ENDDO

      IMDM_SUB_HIGGS=17 ! chargino
      DO IH=1,3
       Z_HiggsEDM=M_C(1)**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FC1)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GC1)
       Z_HiggsEDM=M_C(2)**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FC2)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GC2)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     . -AEM**2*MMU_H**2/2.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(-DREAL(NHC(5,IH))*DREAL(NHC(59,IH))*FC1
     .    +DREAL(NHC(6,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*MMU_H**2/2.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(-DREAL(NHC(5,IH))*DREAL(NHC(68,IH))*FC2
     .    +DREAL(NHC(6,IH))*DREAL(NHC(69,IH))*GC2)
      ENDDO

      IMDM_SUB_HIGGS=18 ! charged Higgs
      DO IH=1,3
       Z_HiggsEDM=MCH**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FCH)
       MDM_MU(IMDM_SUB_HIGGS)=MDM_MU(IMDM_SUB_HIGGS)
     . -AEM*MMU_H**2/16.D0/PI**3*DREAL(NHC(5,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
      ENDDO

*Flip the signs of the fermionic contributions to Barr-Zee MDM
      MDM_MU(11)=-MDM_MU(11)
      MDM_MU(12)=-MDM_MU(12)
      MDM_MU(15)=-MDM_MU(15)
      MDM_MU(17)=-MDM_MU(17)

      MDM_MU(IMDM_SUB)=MDM_MU(11)+MDM_MU(12)+MDM_MU(15)  ! fermion
     .                +MDM_MU(13)+MDM_MU(14)+MDM_MU(16)  ! sfermion
     .                +MDM_MU(17)                        ! chargino
     .                +MDM_MU(18)                        ! charged Higgs


*      print*,mdm_mu(imdm_sub)
*
*TOTAL:
*======
      IMDM_SUB=1
      MDM_MU(IMDM_SUB)=MDM_MU(2)+MDM_MU(3)+MDM_MU(4)+MDM_MU(5)
*
*QED correction; RG evoulution from the SUSY scale to m_mu
*We assume SUSY scale \sim 100 GeV
      QED_LOG=1.D0-4.D0*AEM_H/PI*DLOG(1.D2/MMU_H) ! about 0.93
*      print*,'QED-log.',qed_log
*STORE:Muon MDM:
      RAUX_H(I_MDM_MU+0)=QED_LOG*MDM_MU(1)
      RAUX_H(I_MDM_MU+1)=QED_LOG*MDM_MU(2)
      RAUX_H(I_MDM_MU+2)=QED_LOG*MDM_MU(3)
      RAUX_H(I_MDM_MU+3)=QED_LOG*MDM_MU(4)
      RAUX_H(I_MDM_MU+4)=QED_LOG*MDM_MU(5)
*
*details of (a_mu)_SUSY Barr-Zee : Non-official
      do i=11,18
       raux_h(i_mdm_mu+i)=qed_log*mdm_mu(i)
      enddo
*
      IF(IFLAG(19).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'The SUSY MDM of muon:'
      print*,'---------------------------------------------------------'
      write(*,8) 'a_mu[Total]:',qed_log*mdm_mu(1)
      write(*,8) 'a_mu[C,N,Gl,H]:'
     .      ,qed_log*mdm_mu(2),qed_log*mdm_mu(3)
     .      ,qed_log*mdm_mu(4),qed_log*mdm_mu(5)
      print*,'---------------------------------------------------------'
      ENDIF ! IF(IFLAG(19).EQ.1) THEN
*-----------------------------------------------------------------------
 1    FORMAT(2X,'Masses',1X,E10.4,1X,E10.4,' GeV')
 3    FORMAT(2X,'[L] /','(',E10.4,1X,E10.4,') '
     .                 ,'(',E10.4,1X,E10.4,')',' \\')
 4    FORMAT(2X,'[R] \\','(',E10.4,1X,E10.4,') '
     .                  ,'(',E10.4,1X,E10.4,')',' /')
 5    FORMAT(4X,I1,' (',E10.4,1X,E10.4,') ','(',E10.4,1X,E10.4,') '
     .      ,E10.4)
 6    FORMAT(2X,I1,1X,I1,' (',E10.4,1X,E10.4,') '
     .      ,'(',E10.4,1X,E10.4,') ',E10.4)
 7    FORMAT(2X,A11,1X,2(1X,E10.4,1X))
 8    FORMAT(2X,A24,1X,4(1X,E10.4,1X))
 9    FORMAT(2X,A20,1X,E10.4)
 10   FORMAT(2X,A27,1X,E10.4)
 11   FORMAT(2X,A18,1X,2(1X,E10.4,1X))
*
      RETURN
      END

      REAL*8 FUNCTION MDM_A(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF(DABS(X-1.D0).LT.1.D-8) THEN
       MDM_A=1.D0/24.D0
       RETURN
      ENDIF

      MDM_A=1.D0/12.D0/(1.D0-X)**3*(2.D0+5.D0*X-X**2
     .      +6.D0*X*DLOG(X)/(1.D0-X))
*
      RETURN
      END

      REAL*8 FUNCTION MDM_B(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF(DABS(X-1.D0).LT.1.D-8) THEN
       MDM_B=1.D0/24.D0
       RETURN
      ENDIF

      MDM_B=1.D0/12.D0/(1.D0-X)**3*(1.D0-5.D0*X-2.D0*X**2
     .      -6.D0*X**2*DLOG(X)/(1.D0-X))
*
      RETURN
      END
