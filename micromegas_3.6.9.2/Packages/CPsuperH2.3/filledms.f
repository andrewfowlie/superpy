      SUBROUTINE FILLEDMS(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     . ,MCH,HMASS,OMIX
     . ,STMASS,STMIX,SBMASS,SBMIX,STAUMASS,STAUMIX,SNU3MASS
     . ,M_C,C_L,C_R,M_N,N_N,NCMAX,NHC,SHC,CHC)
************************************************************************
*
*When zq>0.1, user SHOULD provide the subroutine for the calculation of 
*the loop function H for the gluino contribution to the Weinberg operator 
*in the form of
*
*      SUBROUTINE EDM_HH_USER(STMASS,SBMASS,M3,MBMT,MTMT
*     .                      ,EDM_HH_STOP,EDM_HH_SBOT)
*
*
*The function H:
*
*               1   /1   /1   /1               N1 N2
* H(z1,z2,zq)= ---  | dx | du | dy x (1-x) u --------
*               2   /0   /0   /0                D^4
*
* N1 = u (1-x) + zq x (1-x) (1-u) - 2 u x [z1 y + z2 (1-y)]
* N2 = (1-x)^2 (1-u)^2 + u^2 - x^2 (1-u)^2/9
* D  = u (1-x) + zq x (1-x) (1-u) + u x [z1 y + z2 (1-y)]
*[Refs.: Dai, Dykstra, Leigh, Paban, Dicus, PLB237(1990)216]
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
      COMPLEX*16 A_E,A_U,A_C,A_D,A_S
*-----------------------------------------------------------------------
*Local arrays
      COMPLEX*16 XK_E,EGD(3,3),EGU(3,3)
      COMPLEX*16 H_E,H_U,H_C,H_D,H_S
      COMPLEX*16 XRL,XLR,H_Q,A_Q
      REAL*8     SEMASS(2),SUMASS(2),SCMASS(2),SDMASS(2),SSMASS(2)
      COMPLEX*16 SEMIX(2,2),SUMIX(2,2),SCMIX(2,2),SDMIX(2,2),SSMIX(2,2)
*Couplings at Mt^pole
      COMPLEX*16 CGL,CGR,H_F,SFMIX(2,2)
*     chargino(I)-electron-snutrino_e
      COMPLEX*16 GL_CI_E_SNE(2),  GR_CI_E_SNE(2)
*     chargino(I)-up quark-sdown(J)
      COMPLEX*16 GL_CI_U_SDJ(2,2),GR_CI_U_SDJ(2,2)
*     chargino(I)-down quark-sup(J)
      COMPLEX*16 GL_CI_D_SUJ(2,2),GR_CI_D_SUJ(2,2)
*     chargino(I)-strange quark-scharm(J)
      COMPLEX*16 GL_CI_S_SCJ(2,2),GR_CI_S_SCJ(2,2)
*     neutralino(I)-electron-selectron(J)
      COMPLEX*16 GL_NI_E_SEJ(4,2),GR_NI_E_SEJ(4,2)
*     neutralino(I)-up quark-sup(J)
      COMPLEX*16 GL_NI_U_SUJ(4,2),GR_NI_U_SUJ(4,2)
*     neutralino(I)-down quark-sdown(J)
      COMPLEX*16 GL_NI_D_SDJ(4,2),GR_NI_D_SDJ(4,2)
*     neutralino(I)-strange quark-sstrange(J)
      COMPLEX*16 GL_NI_S_SSJ(4,2),GR_NI_S_SSJ(4,2)
*     gluino-up quark-sup(I)
      COMPLEX*16 GL_GL_U_SUI(2),  GR_GL_U_SUI(2)
*     gluino-down quark-sdown(I)
      COMPLEX*16 GL_GL_D_SDI(2),  GR_GL_D_SDI(2)
*     gluino-strange quark-sstrange(I)
      COMPLEX*16 GL_GL_S_SSI(2),  GR_GL_S_SSI(2)
*     gluino-top quark-stop(I)
      COMPLEX*16 GL_GL_T_STI(2),  GR_GL_T_STI(2)
*     gluino-bottom quark-sbottom(I)
      COMPLEX*16 GL_GL_B_SBI(2),  GR_GL_B_SBI(2)
*Electric EDMs of electron and up, down, and strange  quarks: [d^E_f/e]
      REAL*8     DEOE_E(20),DEOE_U(20),DEOE_D(20),DEOE_S(20)
*Chromo-electric EDMs of up quark and down quark: [d^C_q]
      REAL*8     DC_U(20),DC_D(20)
*Purely-gluonic dimension six Weinberg operator: [d^G]
      REAL*8     DG_WEINBERG(20)
      DATA       NEDM_SUB/20/
*Four-fermion couplings at M_F(P) and/or 1 GeV when M_F(P)<1 GeV
      COMPLEX*16 C4_F_FP
*-----------------------------------------------------------------------
      EXTERNAL EDM_A,EDM_B,EDM_C,EDM_H
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
      A_E=SSPARA(27)*DCMPLX(DCOS(SSPARA(28)/180.D0*PI)
     .                     ,DSIN(SSPARA(28)/180.D0*PI))
      A_U=SSPARA(31)*DCMPLX(DCOS(SSPARA(32)/180.D0*PI)
     .                     ,DSIN(SSPARA(32)/180.D0*PI))
      A_C=SSPARA(33)*DCMPLX(DCOS(SSPARA(34)/180.D0*PI)
     .                     ,DSIN(SSPARA(34)/180.D0*PI))
      A_D=SSPARA(35)*DCMPLX(DCOS(SSPARA(36)/180.D0*PI)
     .                     ,DSIN(SSPARA(36)/180.D0*PI))
      A_S=SSPARA(37)*DCMPLX(DCOS(SSPARA(38)/180.D0*PI)
     .                     ,DSIN(SSPARA(38)/180.D0*PI))
*      print*,'A_E = ',a_e
*      print*,'A_U = ',a_u
*      print*,'A_C = ',a_c
*      print*,'A_D = ',a_d
*      print*,'A_S = ',a_s
*
*Basically 10 slots for each EDM except C^4f_ff'/m_f(')
      I_DEOE_E     =200 ! d^E_e/e [cm]
      I_DEOE_U     =210 ! d^E_u/e [cm]
      I_DEOE_D     =220 ! d^E_d/e [cm]
      I_DEOE_S     =230 ! d^E_s/e [cm]
      I_DC_U       =240 ! d^C_u   [cm]
      I_DC_D       =250 ! d^C_d   [cm]
      I_DC_S       =400 ! d^C_s   [cm]
      I_DG_WEINBERG=260 ! d^G     [cm/GeV]
      I_CSPP       =270 ! C_S, C_P C'_P [cm/GeV]
      I_C4FOM      =280 ! C4_ff'/m_f(') [cm/GeV^2]: needs 20 slots
      I_TL         =300 ! d^Tl [e cm]
      I_N1         =310 ! d^n  [e cm]: Chiral Quark Model
      I_N2         =320 ! d^n  [e cm]: Parton Quark Model
      I_N3         =330 ! d^n  [e cm]: QCD Sum Rule Technique
      I_HG         =340 ! d^Hg [e cm]
      I_DEUT       =350 ! d^D  [e cm]
*-----------------------------------------------------------------------
*      print*,a_e,a_u,a_c,a_d,a_s
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
      ASMT=ASMT_H
      GSMT=2.D0*DSQRT(PI*ASMT)
*Stop an Sbottom scales
      BT=(11.D0-2.D0*6.D0/3.D0)/(4.D0*PI)
       QQB=MQ3**2+MBMT_H**2
       QBB=MD3**2+MBMT_H**2
      QB2=DMAX1(QQB,QBB)
       QQT=MQ3**2+MTPOLE_H**2
       QTT=MU3**2+MTPOLE_H**2
      QT2=DMAX1(QQT,QTT)
*      print*,'QB2,QT2=',qb2,qt2
*AS(Stop,Sbottom)
      AS_MSB=ASMT/(1.D0+BT*ASMT*DLOG(QB2/MTPOLE_H**2))
      AS_MST=ASMT/(1.D0+BT*ASMT*DLOG(QT2/MTPOLE_H**2))
*      print*,'As(Mt^pole),As(M_sbottom),AS(M_stop)=',asmt,as_msb,as_mst
*
*The gluino contribution to the threshold corrections \Delta_d*sqrt(2)/v_2
*and \Delta_u*sqrt(2)/v_1
      DO I=1,3
       DO J=1,3
        EGD(I,J)=DCMPLX(0.D0,0.D0)
        EGU(I,J)=DCMPLX(0.D0,0.D0)
       ENDDO
      ENDDO
      EGD(1,1)=2.D0*AS_MSB/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MD1**2,MQ1**2,CDABS(M3_H)**2)
      EGD(2,2)=2.D0*AS_MSB/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MD2**2,MQ2**2,CDABS(M3_H)**2)
      EGD(3,3)=2.D0*AS_MSB/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MD3**2,MQ3**2,CDABS(M3_H)**2)
      EGU(1,1)=2.D0*AS_MST/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MU1**2,MQ1**2,CDABS(M3_H)**2)
      EGU(2,2)=2.D0*AS_MST/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MU2**2,MQ2**2,CDABS(M3_H)**2)
      EGU(3,3)=2.D0*AS_MST/3.D0/PI*DCONJG(MU_H*M3_H)
     .         *F_I(MU3**2,MQ3**2,CDABS(M3_H)**2)
*      print*,'EGD(1,1): ',egd(1,1)
*      print*,'EGD(2,2): ',egd(2,2)
*      print*,'EGD(3,3): ',egd(3,3)
*      print*,'EGU(1,1): ',egu(1,1),egd(1,1)/as_msb*as_mst
*      print*,'EGU(2,2): ',egu(2,2),egd(2,2)/as_msb*as_mst
*      print*,'EGU(3,3): ',egu(3,3),egd(3,3)/as_msb*as_mst
*
*      print*,h_e
      IF(IFLAG(10).EQ.0) THEN
       H_U=DCMPLX(DSQRT(2.D0)*MUMT_H/V_H/SB_H,0.D0)/(1.D0+EGU(1,1)/TB_H)
       H_C=DCMPLX(DSQRT(2.D0)*MCMT_H/V_H/SB_H,0.D0)/(1.D0+EGU(2,2)/TB_H)
*       print*,h_u,h_c
       H_D=DCMPLX(DSQRT(2.D0)*MDMT_H/V_H/CB_H,0.D0)/(1.D0+TB_H*EGD(1,1))
       H_S=DCMPLX(DSQRT(2.D0)*MSMT_H/V_H/CB_H,0.D0)/(1.D0+TB_H*EGD(2,2))
*      print*,h_d,h_s
      ELSE
       H_U=DCMPLX(DSQRT(2.D0)*MUMT_H/V_H/SB_H,0.D0)
       H_C=DCMPLX(DSQRT(2.D0)*MCMT_H/V_H/SB_H,0.D0)
*       print*,h_u,h_c
       H_D=DCMPLX(DSQRT(2.D0)*MDMT_H/V_H/CB_H,0.D0)
       H_S=DCMPLX(DSQRT(2.D0)*MSMT_H/V_H/CB_H,0.D0)
*       print*,h_d,h_s
      ENDIF

*-----------------------------------------------------------------------
*Selectron mixing
*
      H_E0=DCMPLX(DSQRT(2.D0)*ME_H/V_H/CB_H,0.D0)
      XLL = ML1**2+DABS(ME_H)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(SW_H**2-1.D0/2.D0)
      XRR = ME1**2+DABS(ME_H)**2
     .     -(CB_H**2-SB_H**2)*MZ_H**2*SW_H**2
      XRL = H_E0*V_H*CB_H*(A_E-DCONJG(MU_H)*SB_H/CB_H)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SEMASS,SEMIX)
      IF(SEMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative selectron mass squared'
        RETURN
      ENDIF

*EW threshold corrections: see Ellis, Lee, Pilaftsis, hep-ph/0404167
      SNU1MASS=DSQRT(ML1**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0)
      IF(IFLAG(10).EQ.0) THEN
       XK_E=DCONJG(MU_H)*AEM_H/4.D0/PI*
     . (
     .  -DCONJG(M2_H)/SW_H**2*
     .   (F_I(SNU1MASS**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .   +F_I(SEMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(1,1))**2/2.D0
     .   +F_I(SEMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(1,2))**2/2.D0)  ! M2_H
     .  +DCONJG(M1_H)/CW_H**2*
     .   (F_I(SEMASS(1)**2,SEMASS(2)**2,CDABS(M1_H)**2)
     .   +F_I(SEMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(1,1))**2/2.D0
     .   +F_I(SEMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(1,2))**2/2.D0
     .   -F_I(SEMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(2,1))**2
     .   -F_I(SEMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *CDABS(SEMIX(2,2))**2)       ! M1_H
     . )
       H_E=DCMPLX(DSQRT(2.D0)*ME_H/V_H/CB_H,0.D0)/(1.D0+XK_E*TB_H)
      ELSE
       H_E=DCMPLX(DSQRT(2.D0)*ME_H/V_H/CB_H,0.D0)
      ENDIF
*      print*,h_e0,h_e,xk_e

*Some iteration... may not needed
*      IF(IFLAG(10).EQ.0) THEN
*
*       DO III=0,9
*        XLL = ML1**2+DABS(ME_H)**2
*     .       +(CB_H**2-SB_H**2)*MZ_H**2*(SW_H**2-1.D0/2.D0)
*        XRR = ME1**2+DABS(ME_H)**2
*     .       -(CB_H**2-SB_H**2)*MZ_H**2*SW_H**2
*        XRL = H_E*V_H*CB_H*(A_E-DCONJG(MU_H)*SB_H/CB_H)/SQRT(2.D0)
*        XLR = DCONJG(XRL)
*
*        CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SEMASS,SEMIX)
*        IF(SEMASS(1).LE.0.D0) THEN
*          print*,'ERROR!: Negative selectron mass squared'
*          RETURN
*        ENDIF
*        XK_E=DCONJG(MU_H)*AEM_H/4.D0/PI*
*     .  (
*     .   -DCONJG(M2_H)/SW_H**2*
*     .    (F_I(SNU1MASS**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
*     .    +F_I(SEMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(1,1))**2/2.D0
*     .    +F_I(SEMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(1,2))**2/2.D0)  ! M2_H
*     .   +DCONJG(M1_H)/CW_H**2*
*     .    (F_I(SEMASS(1)**2,SEMASS(2)**2,CDABS(M1_H)**2)
*     .    +F_I(SEMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(1,1))**2/2.D0
*     .    +F_I(SEMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(1,2))**2/2.D0
*     .    -F_I(SEMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(2,1))**2
*     .    -F_I(SEMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
*     .     *CDABS(SEMIX(2,2))**2)       ! M1_H
*     .  )
*        H_E=DCMPLX(DSQRT(2.D0)*ME_H/V_H/CB_H,0.D0)/(1.D0+XK_E*TB_H)
*        print*,iii,h_e,xk_e
*       ENDDO ! III
*
*      ENDIF

*      print*,'Selectron Sector:'
*      WRITE(*,1) SEMASS(1),SEMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SEMIX(1,1)),DIMAG(SEMIX(1,1))
*     .          ,DREAL(SEMIX(1,2)),DIMAG(SEMIX(1,2))
*      WRITE(*,4) DREAL(SEMIX(2,1)),DIMAG(SEMIX(2,1))
*     .          ,DREAL(SEMIX(2,2)),DIMAG(SEMIX(2,2))
*      print*,' '
*

*
      T_3=1.D0/2.D0
      Q_Q=2.D0/3.D0
      V_Q=V_H*SB_H
      R_Q=CB_H/SB_H
      M_Q=MUMT_H
      H_Q=H_U
      A_Q=A_U

      XLL = MQ1**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MU1**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SUMASS,SUMIX)
      IF(SUMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative sup mass squared'
        RETURN
      ENDIF

*      print*,'Sup Sector:'
*      WRITE(*,1) SUMASS(1),SUMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SUMIX(1,1)),DIMAG(SUMIX(1,1))
*     .          ,DREAL(SUMIX(1,2)),DIMAG(SUMIX(1,2))
*      WRITE(*,4) DREAL(SUMIX(2,1)),DIMAG(SUMIX(2,1))
*     .          ,DREAL(SUMIX(2,2)),DIMAG(SUMIX(2,2))
*      print*,' '
*
*Scharm mixing
*
      T_3=1.D0/2.D0
      Q_Q=2.D0/3.D0
      V_Q=V_H*SB_H
      R_Q=CB_H/SB_H
      M_Q=MCMT_H
      H_Q=H_C
      A_Q=A_C

      XLL = MQ2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MU2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SCMASS,SCMIX)
      IF(SCMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative scharm mass squared'
        RETURN
      ENDIF

*      print*,'Scharm Sector:'
*      WRITE(*,1) SCMASS(1),SCMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SCMIX(1,1)),DIMAG(SCMIX(1,1))
*     .          ,DREAL(SCMIX(1,2)),DIMAG(SCMIX(1,2))
*      WRITE(*,4) DREAL(SCMIX(2,1)),DIMAG(SCMIX(2,1))
*     .          ,DREAL(SCMIX(2,2)),DIMAG(SCMIX(2,2))
*      print*,' '
*
*Sdown mixing
*
      T_3=-1.D0/2.D0
      Q_Q=-1.D0/3.D0
      V_Q=V_H*CB_H
      R_Q=SB_H/CB_H
      M_Q=MDMT_H
      H_Q=H_D
      A_Q=A_D

      XLL = MQ1**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MD1**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SDMASS,SDMIX)
      IF(SDMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative sdown mass squared'
        RETURN
      ENDIF

*      print*,'Sdown Sector:'
*      WRITE(*,1) SDMASS(1),SDMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SDMIX(1,1)),DIMAG(SDMIX(1,1))
*     .          ,DREAL(SDMIX(1,2)),DIMAG(SDMIX(1,2))
*      WRITE(*,4) DREAL(SDMIX(2,1)),DIMAG(SDMIX(2,1))
*     .          ,DREAL(SDMIX(2,2)),DIMAG(SDMIX(2,2))
*      print*,' '
*
*Sstrange mixing
*
      T_3=-1.D0/2.D0
      Q_Q=-1.D0/3.D0
      V_Q=V_H*CB_H
      R_Q=SB_H/CB_H
      M_Q=MSMT_H
      H_Q=H_S
      A_Q=A_S

      XLL = MQ2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MD2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SSMASS,SSMIX)
      IF(SSMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative sstrange mass squared'
        RETURN
      ENDIF

*      print*,'Sstrange Sector:'
*      WRITE(*,1) SSMASS(1),SSMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SSMIX(1,1)),DIMAG(SSMIX(1,1))
*     .          ,DREAL(SSMIX(1,2)),DIMAG(SSMIX(1,2))
*      WRITE(*,4) DREAL(SSMIX(2,1)),DIMAG(SSMIX(2,1))
*     .          ,DREAL(SSMIX(2,2)),DIMAG(SSMIX(2,2))
*      print*,' '
*-----------------------------------------------------------------------
*Couplings:
*     chargino(I)-electron-snutrino_e
*      COMPLEX*16 GL_CI_E_SNE(2),  GR_CI_E_SNE(2)
*      print*,'chargino(I)-electron-snutrino_e:'
      DO I=1,2
       GL_CI_E_SNE(I)=-GW_H*C_R(I,1)
       GR_CI_E_SNE(I)=DCONJG(H_E)*C_L(I,2)
       CGL=GL_CI_E_SNE(I)
       CGR=GR_CI_E_SNE(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
*     chargino(I)-up quark-sdown(J)
*      COMPLEX*16 GL_CI_U_SDJ(2,2),GR_CI_U_SDJ(2,2)
*      print*,'chargino(I)-up quark-sdown(J):'
      DO I=1,2
       DO J=1,2
        GL_CI_U_SDJ(I,J)=-GW_H*DCONJG(C_L(I,1))*DCONJG(SDMIX(1,J))
     .                    +H_D*DCONJG(C_L(I,2))*DCONJG(SDMIX(2,J))
        GR_CI_U_SDJ(I,J)=DCONJG(H_U*C_R(I,2)*SDMIX(1,J))
        CGL=GL_CI_U_SDJ(I,J)
        CGR=GR_CI_U_SDJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     chargino(I)-down quark-sup(J)
*      COMPLEX*16 GL_CI_D_SUJ(2,2),GR_CI_D_SUJ(2,2)
*      print*,'chargino(I)-down quark-sup(J):'
      DO I=1,2
       DO J=1,2
        GL_CI_D_SUJ(I,J)=-GW_H*C_R(I,1)*DCONJG(SUMIX(1,J))
     .                    +H_U*C_R(I,2)*DCONJG(SUMIX(2,J))
        GR_CI_D_SUJ(I,J)=DCONJG(H_D)*C_L(I,2)*DCONJG(SUMIX(1,J))
        CGL=GL_CI_D_SUJ(I,J)
        CGR=GR_CI_D_SUJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     chargino(I)-strange quark-scharm(J)
*      COMPLEX*16 GL_CI_S_SCJ(2,2),GR_CI_S_SCJ(2,2)
*      print*,'chargino(I)-strange quark-scharm(J):'
      DO I=1,2
       DO J=1,2
        GL_CI_S_SCJ(I,J)=-GW_H*C_R(I,1)*DCONJG(SCMIX(1,J))
     .                    +H_C*C_R(I,2)*DCONJG(SCMIX(2,J))
        GR_CI_S_SCJ(I,J)=DCONJG(H_S)*C_L(I,2)*DCONJG(SCMIX(1,J))
        CGL=GL_CI_S_SCJ(I,J)
        CGR=GR_CI_S_SCJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     neutralino(I)-electron-selectron(J)
*      COMPLEX*16 GL_NI_E_SEJ(4,2),GR_NI_E_SEJ(4,2)
      I_A =3
      T_3F=-1.0/2.D0
      Q_F =-1.D0
      H_F =H_E
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SEMIX(I,J)
       ENDDO
      ENDDO

*      print*,'neutralino(I)-electron-selectron(J):'
      DO I=1,4
       DO J=1,2
        GL_NI_E_SEJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_N(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_N(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_N(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_E_SEJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_N(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_N(I,I_A)*DCONJG(SFMIX(1,J))
        CGL=GL_NI_E_SEJ(I,J)
        CGR=GR_NI_E_SEJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     neutralino(I)-up quark-sup(J)
*      COMPLEX*16 GL_NI_U_SUJ(4,2),GR_NI_U_SUJ(4,2)
      I_A =4
      T_3F=1.0/2.D0
      Q_F =2.D0/3.D0
      H_F =H_U
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SUMIX(I,J)
       ENDDO
      ENDDO

*      print*,'neutralino(I)-up quark-sup(J):'
      DO I=1,4
       DO J=1,2
        GL_NI_U_SUJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_N(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_N(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_N(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_U_SUJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_N(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_N(I,I_A)*DCONJG(SFMIX(1,J))
        CGL=GL_NI_U_SUJ(I,J)
        CGR=GR_NI_U_SUJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     neutralino(I)-down quark-sdown(J)
*      COMPLEX*16 GL_NI_D_SDJ(4,2),GR_NI_D_SDJ(4,2)
      I_A =3
      T_3F=-1.0/2.D0
      Q_F =-1.D0/3.D0
      H_F =H_D
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SDMIX(I,J)
       ENDDO
      ENDDO

*      print*,'neutralino(I)-down quark-sdown(J):'
      DO I=1,4
       DO J=1,2
        GL_NI_D_SDJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_N(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_N(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_N(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_D_SDJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_N(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_N(I,I_A)*DCONJG(SFMIX(1,J))
        CGL=GL_NI_D_SDJ(I,J)
        CGR=GR_NI_D_SDJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     neutralino(I)-strange quark-sstrange(J)
*      COMPLEX*16 GL_NI_S_SSJ(4,2),GR_NI_S_SSJ(4,2)
      I_A =3
      T_3F=-1.0/2.D0
      Q_F =-1.D0/3.D0
      H_F =H_S
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SSMIX(I,J)
       ENDDO
      ENDDO

*      print*,'neutralino(I)-strange quark-sstrange(J):'
      DO I=1,4
       DO J=1,2
        GL_NI_S_SSJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_N(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_N(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_N(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_S_SSJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_N(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_N(I,I_A)*DCONJG(SFMIX(1,J))
        CGL=GL_NI_S_SSJ(I,J)
        CGR=GR_NI_S_SSJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      PHI_3=SSPARA(10)/180.D0*PI ! in Radian
*      print*,gsmt,phi_3,DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))
*     gluino-up quark-sup(I)
*      COMPLEX*16 GL_GL_U_SUI(2),  GR_GL_U_SUI(2)
*      print*,'gluino-up quark-sup(I):'
      DO I=1,2
       GL_GL_U_SUI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(SUMIX(1,I))
       GR_GL_U_SUI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(SUMIX(2,I))
       CGL=GL_GL_U_SUI(I)
       CGR=GR_GL_U_SUI(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
*     gluino-down quark-sdown(I)
*      COMPLEX*16 GL_GL_D_SDI(2),  GR_GL_D_SDI(2)
*      print*,'gluino-down quark-sdown(I):'
      DO I=1,2
       GL_GL_D_SDI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(SDMIX(1,I))
       GR_GL_D_SDI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(SDMIX(2,I))
       CGL=GL_GL_D_SDI(I)
       CGR=GR_GL_D_SDI(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
*     gluino-strange quark-sstrange(I)
*      COMPLEX*16 GL_GL_S_SSI(2),  GR_GL_S_SSI(2)
*      print*,'gluino-strange quark-sstrange(I):'
      DO I=1,2
       GL_GL_S_SSI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(SSMIX(1,I))
       GR_GL_S_SSI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(SSMIX(2,I))
       CGL=GL_GL_S_SSI(I)
       CGR=GR_GL_S_SSI(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
*     gluino-top quark-stop(I)
*      COMPLEX*16 GL_GL_T_STI(2),  GR_GL_T_STI(2)
      DO I=1,2
       GL_GL_T_STI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(STMIX(1,I))
       GR_GL_T_STI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(STMIX(2,I))
      ENDDO
*     gluino-bottom quark-sbottom(I)
*      COMPLEX*16 GL_GL_B_SBI(2),  GR_GL_B_SBI(2)
      DO I=1,2
       GL_GL_B_SBI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(SBMIX(1,I))
       GR_GL_B_SBI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(SBMIX(2,I))
      ENDDO
*Dumping couplings
      IF(IFLAG(18).EQ.3) THEN
      print*,'---------------------------------------------------------'
      print*,' Ino-couplings to fermion and sfermion needed for EDMs:'
      print*,'                            G_R, G_L, and Im[G_R^* G_L]'
      print*,'---------------------------------------------------------'
      print*,'chargino(I)-electron-snutrino_e:'
      DO I=1,2
       CGL=GL_CI_E_SNE(I)
       CGR=GR_CI_E_SNE(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'chargino(I)-up quark-sdown(J):'
      DO I=1,2
       DO J=1,2
        CGL=GL_CI_U_SDJ(I,J)
        CGR=GR_CI_U_SDJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'chargino(I)-down quark-sup(J):'
      DO I=1,2
       DO J=1,2
        CGL=GL_CI_D_SUJ(I,J)
        CGR=GR_CI_D_SUJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'chargino(I)-strange quark-scharm(J):'
      DO I=1,2
       DO J=1,2
        CGL=GL_CI_S_SCJ(I,J)
        CGR=GR_CI_S_SCJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'neutralino(I)-electron-selectron(J):'
      DO I=1,4
       DO J=1,2
        CGL=GL_NI_E_SEJ(I,J)
        CGR=GR_NI_E_SEJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'neutralino(I)-up quark-sup(J):'
      DO I=1,4
       DO J=1,2
        CGL=GL_NI_U_SUJ(I,J)
        CGR=GR_NI_U_SUJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'neutralino(I)-down quark-sdown(J):'
      DO I=1,4
       DO J=1,2
        CGL=GL_NI_D_SDJ(I,J)
        CGR=GR_NI_D_SDJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'neutralino(I)-strange quark-sstrange(J):'
      DO I=1,4
       DO J=1,2
        CGL=GL_NI_S_SSJ(I,J)
        CGR=GR_NI_S_SSJ(I,J)
        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
      print*,'gluino-up quark-sup(I):'
      DO I=1,2
       CGL=GL_GL_U_SUI(I)
       CGR=GR_GL_U_SUI(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'gluino-down quark-sdown(I):'
      DO I=1,2
       CGL=GL_GL_D_SDI(I)
       CGR=GR_GL_D_SDI(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'gluino-strange quark-sstrange(I):'
      DO I=1,2
       CGL=GL_GL_S_SSI(I)
       CGR=GR_GL_S_SSI(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'gluino-top quark-stop(I):'
      DO I=1,2
       CGL=GL_GL_T_STI(I)
       CGR=GR_GL_T_STI(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'gluino-bottom quark-sbottom(I):'
      DO I=1,2
       CGL=GL_GL_B_SBI(I)
       CGR=GR_GL_B_SBI(I)
       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
*
      print*,'---------------------------------------------------------'
      ENDIF ! IF(IFLAG(18).EQ.3) THEN
*-----------------------------------------------------------------------
*Electric EDMs of electron and up, down, and strange  quarks: [d^E_f/e]
*      REAL*8     DEOE_E(20),DEOE_U(20),DEOE_D(20),DEOE_S(20)
*      REAL*8     DC_U(20),DC_D(20)
************************************************************************
*Convention of the arrays DEOE and DC
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
        DEOE_E(I)     =0.D0
        DEOE_U(I)     =0.D0
        DEOE_D(I)     =0.D0
        DEOE_S(I)     =0.D0
        DC_U(I)       =0.D0
        DC_D(I)       =0.D0
       ENDDO
*
*      print*,EDM_A(0.1D0),EDM_A(0.5D0),EDM_A(0.9D0),EDM_A(1.0D0)
*     .                   ,EDM_A(1.1D0),EDM_A(2.0D0),EDM_A(10.D0)
*      print*,EDM_B(0.1D0),EDM_B(0.5D0),EDM_B(0.9D0),EDM_B(1.0D0)
*     .                   ,EDM_B(1.1D0),EDM_B(2.0D0),EDM_B(10.D0)
*      print*,EDM_C(0.1D0),EDM_C(0.5D0),EDM_C(0.9D0),EDM_C(1.0D0)
*     .                   ,EDM_C(1.1D0),EDM_C(2.0D0),EDM_C(10.D0)
*
*Chargino contributions:
*=======================
      IEDM_SUB=2
*
      MSNE=DSQRT(ML1**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0)
*       print*,msne,snu3mass

      DEOE_E(IEDM_SUB)=-1.D0/16.D0/PI**2
     .*( M_C(1)/MSNE**2*DIMAG(DCONJG(GR_CI_E_SNE(1))*GL_CI_E_SNE(1))
     .                *EDM_A(M_C(1)**2/MSNE**2)
     .  +M_C(2)/MSNE**2*DIMAG(DCONJG(GR_CI_E_SNE(2))*GL_CI_E_SNE(2))
     .                *EDM_A(M_C(2)**2/MSNE**2) )

      DEOE_U(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_C(1)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(1,1))*GL_CI_U_SDJ(1,1))
     .             *( EDM_A(M_C(1)**2/SDMASS(1)**2)
     .               -EDM_B(M_C(1)**2/SDMASS(1)**2)/3.D0)
     .  +M_C(1)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(1,2))*GL_CI_U_SDJ(1,2))
     .             *( EDM_A(M_C(1)**2/SDMASS(2)**2)
     .               -EDM_B(M_C(1)**2/SDMASS(2)**2)/3.D0)
     .  +M_C(2)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(2,1))*GL_CI_U_SDJ(2,1))
     .             *( EDM_A(M_C(2)**2/SDMASS(1)**2)
     .               -EDM_B(M_C(2)**2/SDMASS(1)**2)/3.D0)
     .  +M_C(2)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(2,2))*GL_CI_U_SDJ(2,2))
     .             *( EDM_A(M_C(2)**2/SDMASS(2)**2)
     .               -EDM_B(M_C(2)**2/SDMASS(2)**2)/3.D0) )

      DEOE_D(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_C(1)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(1,1))*GL_CI_D_SUJ(1,1))
     .             *(     -EDM_A(M_C(1)**2/SUMASS(1)**2)
     .               +2.D0*EDM_B(M_C(1)**2/SUMASS(1)**2)/3.D0)
     .  +M_C(1)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(1,2))*GL_CI_D_SUJ(1,2))
     .             *(     -EDM_A(M_C(1)**2/SUMASS(2)**2)
     .               +2.D0*EDM_B(M_C(1)**2/SUMASS(2)**2)/3.D0)
     .  +M_C(2)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(2,1))*GL_CI_D_SUJ(2,1))
     .             *(     -EDM_A(M_C(2)**2/SUMASS(1)**2)
     .               +2.D0*EDM_B(M_C(2)**2/SUMASS(1)**2)/3.D0)
     .  +M_C(2)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(2,2))*GL_CI_D_SUJ(2,2))
     .             *(     -EDM_A(M_C(2)**2/SUMASS(2)**2)
     .               +2.D0*EDM_B(M_C(2)**2/SUMASS(2)**2)/3.D0) )

      DEOE_S(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_C(1)/SCMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(1,1))*GL_CI_S_SCJ(1,1))
     .             *(     -EDM_A(M_C(1)**2/SCMASS(1)**2)
     .               +2.D0*EDM_B(M_C(1)**2/SCMASS(1)**2)/3.D0)
     .  +M_C(1)/SCMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(1,2))*GL_CI_S_SCJ(1,2))
     .             *(     -EDM_A(M_C(1)**2/SCMASS(2)**2)
     .               +2.D0*EDM_B(M_C(1)**2/SCMASS(2)**2)/3.D0)
     .  +M_C(2)/SCMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(2,1))*GL_CI_S_SCJ(2,1))
     .             *(     -EDM_A(M_C(2)**2/SCMASS(1)**2)
     .               +2.D0*EDM_B(M_C(2)**2/SCMASS(1)**2)/3.D0)
     .  +M_C(2)/SCMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(2,2))*GL_CI_S_SCJ(2,2))
     .             *(     -EDM_A(M_C(2)**2/SCMASS(2)**2)
     .               +2.D0*EDM_B(M_C(2)**2/SCMASS(2)**2)/3.D0) )

      DC_U(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( M_C(1)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(1,1))*GL_CI_U_SDJ(1,1))
     .             *EDM_B(M_C(1)**2/SDMASS(1)**2)
     .  +M_C(1)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(1,2))*GL_CI_U_SDJ(1,2))
     .             *EDM_B(M_C(1)**2/SDMASS(2)**2)
     .  +M_C(2)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(2,1))*GL_CI_U_SDJ(2,1))
     .             *EDM_B(M_C(2)**2/SDMASS(1)**2)
     .  +M_C(2)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_U_SDJ(2,2))*GL_CI_U_SDJ(2,2))
     .             *EDM_B(M_C(2)**2/SDMASS(2)**2) )

      DC_D(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( M_C(1)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(1,1))*GL_CI_D_SUJ(1,1))
     .             *EDM_B(M_C(1)**2/SUMASS(1)**2)
     .  +M_C(1)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(1,2))*GL_CI_D_SUJ(1,2))
     .             *EDM_B(M_C(1)**2/SUMASS(2)**2)
     .  +M_C(2)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(2,1))*GL_CI_D_SUJ(2,1))
     .             *EDM_B(M_C(2)**2/SUMASS(1)**2)
     .  +M_C(2)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_D_SUJ(2,2))*GL_CI_D_SUJ(2,2))
     .             *EDM_B(M_C(2)**2/SUMASS(2)**2) )
*
*Neutralino contributions:
*=========================
      IEDM_SUB=3

      Q_SF=-1.D0
      DEOE_E(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_N(1)/SEMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(1,1))*GL_NI_E_SEJ(1,1))
     .             *Q_SF*EDM_B(M_N(1)**2/SEMASS(1)**2)
     .  +M_N(1)/SEMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(1,2))*GL_NI_E_SEJ(1,2))
     .             *Q_SF*EDM_B(M_N(1)**2/SEMASS(2)**2)
     .  +M_N(2)/SEMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(2,1))*GL_NI_E_SEJ(2,1))
     .             *Q_SF*EDM_B(M_N(2)**2/SEMASS(1)**2)
     .  +M_N(2)/SEMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(2,2))*GL_NI_E_SEJ(2,2))
     .             *Q_SF*EDM_B(M_N(2)**2/SEMASS(2)**2)
     .  +M_N(3)/SEMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(3,1))*GL_NI_E_SEJ(3,1))
     .             *Q_SF*EDM_B(M_N(3)**2/SEMASS(1)**2)
     .  +M_N(3)/SEMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(3,2))*GL_NI_E_SEJ(3,2))
     .             *Q_SF*EDM_B(M_N(3)**2/SEMASS(2)**2)
     .  +M_N(4)/SEMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(4,1))*GL_NI_E_SEJ(4,1))
     .             *Q_SF*EDM_B(M_N(4)**2/SEMASS(1)**2)
     .  +M_N(4)/SEMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_E_SEJ(4,2))*GL_NI_E_SEJ(4,2))
     .             *Q_SF*EDM_B(M_N(4)**2/SEMASS(2)**2) )

      Q_SF=2.D0/3.D0
      DEOE_U(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_N(1)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(1,1))*GL_NI_U_SUJ(1,1))
     .             *Q_SF*EDM_B(M_N(1)**2/SUMASS(1)**2)
     .  +M_N(1)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(1,2))*GL_NI_U_SUJ(1,2))
     .             *Q_SF*EDM_B(M_N(1)**2/SUMASS(2)**2)
     .  +M_N(2)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(2,1))*GL_NI_U_SUJ(2,1))
     .             *Q_SF*EDM_B(M_N(2)**2/SUMASS(1)**2)
     .  +M_N(2)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(2,2))*GL_NI_U_SUJ(2,2))
     .             *Q_SF*EDM_B(M_N(2)**2/SUMASS(2)**2)
     .  +M_N(3)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(3,1))*GL_NI_U_SUJ(3,1))
     .             *Q_SF*EDM_B(M_N(3)**2/SUMASS(1)**2)
     .  +M_N(3)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(3,2))*GL_NI_U_SUJ(3,2))
     .             *Q_SF*EDM_B(M_N(3)**2/SUMASS(2)**2)
     .  +M_N(4)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(4,1))*GL_NI_U_SUJ(4,1))
     .             *Q_SF*EDM_B(M_N(4)**2/SUMASS(1)**2)
     .  +M_N(4)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(4,2))*GL_NI_U_SUJ(4,2))
     .             *Q_SF*EDM_B(M_N(4)**2/SUMASS(2)**2) )

      Q_SF=-1.D0/3.D0
      DEOE_D(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_N(1)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(1,1))*GL_NI_D_SDJ(1,1))
     .             *Q_SF*EDM_B(M_N(1)**2/SDMASS(1)**2)
     .  +M_N(1)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(1,2))*GL_NI_D_SDJ(1,2))
     .             *Q_SF*EDM_B(M_N(1)**2/SDMASS(2)**2)
     .  +M_N(2)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(2,1))*GL_NI_D_SDJ(2,1))
     .             *Q_SF*EDM_B(M_N(2)**2/SDMASS(1)**2)
     .  +M_N(2)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(2,2))*GL_NI_D_SDJ(2,2))
     .             *Q_SF*EDM_B(M_N(2)**2/SDMASS(2)**2)
     .  +M_N(3)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(3,1))*GL_NI_D_SDJ(3,1))
     .             *Q_SF*EDM_B(M_N(3)**2/SDMASS(1)**2)
     .  +M_N(3)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(3,2))*GL_NI_D_SDJ(3,2))
     .             *Q_SF*EDM_B(M_N(3)**2/SDMASS(2)**2)
     .  +M_N(4)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(4,1))*GL_NI_D_SDJ(4,1))
     .             *Q_SF*EDM_B(M_N(4)**2/SDMASS(1)**2)
     .  +M_N(4)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(4,2))*GL_NI_D_SDJ(4,2))
     .             *Q_SF*EDM_B(M_N(4)**2/SDMASS(2)**2) )

      Q_SF=-1.D0/3.D0
      DEOE_S(IEDM_SUB)= 1.D0/16.D0/PI**2
     .*( M_N(1)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(1,1))*GL_NI_S_SSJ(1,1))
     .             *Q_SF*EDM_B(M_N(1)**2/SSMASS(1)**2)
     .  +M_N(1)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(1,2))*GL_NI_S_SSJ(1,2))
     .             *Q_SF*EDM_B(M_N(1)**2/SSMASS(2)**2)
     .  +M_N(2)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(2,1))*GL_NI_S_SSJ(2,1))
     .             *Q_SF*EDM_B(M_N(2)**2/SSMASS(1)**2)
     .  +M_N(2)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(2,2))*GL_NI_S_SSJ(2,2))
     .             *Q_SF*EDM_B(M_N(2)**2/SSMASS(2)**2)
     .  +M_N(3)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(3,1))*GL_NI_S_SSJ(3,1))
     .             *Q_SF*EDM_B(M_N(3)**2/SSMASS(1)**2)
     .  +M_N(3)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(3,2))*GL_NI_S_SSJ(3,2))
     .             *Q_SF*EDM_B(M_N(3)**2/SSMASS(2)**2)
     .  +M_N(4)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(4,1))*GL_NI_S_SSJ(4,1))
     .             *Q_SF*EDM_B(M_N(4)**2/SSMASS(1)**2)
     .  +M_N(4)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(4,2))*GL_NI_S_SSJ(4,2))
     .             *Q_SF*EDM_B(M_N(4)**2/SSMASS(2)**2) )

      DC_U(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( M_N(1)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(1,1))*GL_NI_U_SUJ(1,1))
     .             *EDM_B(M_N(1)**2/SUMASS(1)**2)
     .  +M_N(1)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(1,2))*GL_NI_U_SUJ(1,2))
     .             *EDM_B(M_N(1)**2/SUMASS(2)**2)
     .  +M_N(2)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(2,1))*GL_NI_U_SUJ(2,1))
     .             *EDM_B(M_N(2)**2/SUMASS(1)**2)
     .  +M_N(2)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(2,2))*GL_NI_U_SUJ(2,2))
     .             *EDM_B(M_N(2)**2/SUMASS(2)**2)
     .  +M_N(3)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(3,1))*GL_NI_U_SUJ(3,1))
     .             *EDM_B(M_N(3)**2/SUMASS(1)**2)
     .  +M_N(3)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(3,2))*GL_NI_U_SUJ(3,2))
     .             *EDM_B(M_N(3)**2/SUMASS(2)**2)
     .  +M_N(4)/SUMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(4,1))*GL_NI_U_SUJ(4,1))
     .             *EDM_B(M_N(4)**2/SUMASS(1)**2)
     .  +M_N(4)/SUMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_U_SUJ(4,2))*GL_NI_U_SUJ(4,2))
     .             *EDM_B(M_N(4)**2/SUMASS(2)**2) )

      DC_D(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( M_N(1)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(1,1))*GL_NI_D_SDJ(1,1))
     .             *EDM_B(M_N(1)**2/SDMASS(1)**2)
     .  +M_N(1)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(1,2))*GL_NI_D_SDJ(1,2))
     .             *EDM_B(M_N(1)**2/SDMASS(2)**2)
     .  +M_N(2)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(2,1))*GL_NI_D_SDJ(2,1))
     .             *EDM_B(M_N(2)**2/SDMASS(1)**2)
     .  +M_N(2)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(2,2))*GL_NI_D_SDJ(2,2))
     .             *EDM_B(M_N(2)**2/SDMASS(2)**2)
     .  +M_N(3)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(3,1))*GL_NI_D_SDJ(3,1))
     .             *EDM_B(M_N(3)**2/SDMASS(1)**2)
     .  +M_N(3)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(3,2))*GL_NI_D_SDJ(3,2))
     .             *EDM_B(M_N(3)**2/SDMASS(2)**2)
     .  +M_N(4)/SDMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(4,1))*GL_NI_D_SDJ(4,1))
     .             *EDM_B(M_N(4)**2/SDMASS(1)**2)
     .  +M_N(4)/SDMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_D_SDJ(4,2))*GL_NI_D_SDJ(4,2))
     .             *EDM_B(M_N(4)**2/SDMASS(2)**2) )
*
*Gluino contributions:
*=====================
      IEDM_SUB=4

      Q_SF=-1.D0
      DEOE_E(IEDM_SUB)= 0.D0

      Q_SF=2.D0/3.D0
      DEOE_U(IEDM_SUB)= 1.D0/3.D0/PI**2
     .*( CDABS(M3_H)/SUMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_U_SUI(1))*GL_GL_U_SUI(1))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SUMASS(1)**2)
     .  +CDABS(M3_H)/SUMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_U_SUI(2))*GL_GL_U_SUI(2))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SUMASS(2)**2) )

      Q_SF=-1.D0/3.D0
      DEOE_D(IEDM_SUB)= 1.D0/3.D0/PI**2
     .*( CDABS(M3_H)/SDMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_D_SDI(1))*GL_GL_D_SDI(1))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SDMASS(1)**2)
     .  +CDABS(M3_H)/SDMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_D_SDI(2))*GL_GL_D_SDI(2))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SDMASS(2)**2) )

      Q_SF=-1.D0/3.D0
      DEOE_S(IEDM_SUB)= 1.D0/3.D0/PI**2
     .*( CDABS(M3_H)/SSMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_S_SSI(1))*GL_GL_S_SSI(1))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SSMASS(1)**2)
     .  +CDABS(M3_H)/SSMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_S_SSI(2))*GL_GL_S_SSI(2))
     .                  *Q_SF*EDM_B(CDABS(M3_H)**2/SSMASS(2)**2) )

      DC_U(IEDM_SUB)= -GSMT/8.D0/PI**2
     .*( CDABS(M3_H)/SUMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_U_SUI(1))*GL_GL_U_SUI(1))
     .                  *EDM_C(CDABS(M3_H)**2/SUMASS(1)**2)
     .  +CDABS(M3_H)/SUMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_U_SUI(2))*GL_GL_U_SUI(2))
     .                  *EDM_C(CDABS(M3_H)**2/SUMASS(2)**2) )

      DC_D(IEDM_SUB)= -GSMT/8.D0/PI**2
     .*( CDABS(M3_H)/SDMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_D_SDI(1))*GL_GL_D_SDI(1))
     .                  *EDM_C(CDABS(M3_H)**2/SDMASS(1)**2)
     .  +CDABS(M3_H)/SDMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_D_SDI(2))*GL_GL_D_SDI(2))
     .                  *EDM_C(CDABS(M3_H)**2/SDMASS(2)**2) )
*
*Two-loop Higgs contributions: the contribution to
*the electron EDMs        : d^E_e/e, d^E_u/e, d^E_d/e, d^E_s/e
*the chromo-electric CEDMs: d^C_u and d^C_d
      IEDM_SUB=5

*(d^E_e/e)^H, (d^E_u/e)^H, (d^E_d/e)^H, (d^E_s/e)^H:
*---------------------------------------------------

      IEDM_SUB_HIGGS=11 ! top
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MTMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTOP)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTOP)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*ME_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(3,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(2,IH))*DREAL(NHC(27,IH))*GTOP)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*MUMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(21,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(20,IH))*DREAL(NHC(27,IH))*GTOP)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*MDMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(12,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(11,IH))*DREAL(NHC(27,IH))*GTOP)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     .  -3.D0*AEM**2*QQ**2*MSMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(15,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(14,IH))*DREAL(NHC(27,IH))*GTOP)
      ENDDO

      IEDM_SUB_HIGGS=12 ! bottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=MBMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FBOT)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GBOT)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*ME_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(3,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(2,IH))*DREAL(NHC(18,IH))*GBOT)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*MUMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(21,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(20,IH))*DREAL(NHC(18,IH))*GBOT)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*MDMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(12,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(11,IH))*DREAL(NHC(18,IH))*GBOT)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . -3.D0*AEM**2*QQ**2*MSMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(15,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(14,IH))*DREAL(NHC(18,IH))*GBOT)
      ENDDO

      IEDM_SUB_HIGGS=13 ! stop
      QQ=2.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=STMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST1)
       Z_HiggsEDM=STMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST2)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*ME_H/32.D0/PI**3*DREAL(NHC(3,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MUMT_H/32.D0/PI**3*DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MDMT_H/32.D0/PI**3*DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MSMT_H/32.D0/PI**3*DREAL(NHC(15,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
      ENDDO

      IEDM_SUB_HIGGS=14 ! sbottom
      QQ=-1.D0/3.D0
      DO IH=1,3
       Z_HiggsEDM=SBMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB1)
       Z_HiggsEDM=SBMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB2)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*ME_H/32.D0/PI**3*DREAL(NHC(3,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MUMT_H/32.D0/PI**3*DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MDMT_H/32.D0/PI**3*DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . +3.D0*AEM*QQ**2*MSMT_H/32.D0/PI**3*DREAL(NHC(15,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
      ENDDO

      IEDM_SUB_HIGGS=15 ! tau
      DO IH=1,3
       Z_HiggsEDM=MTAU_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTAU)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTAU)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . -AEM**2*ME_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(3,IH))*DREAL(NHC(8,IH))*FTAU
     .   +DREAL(NHC(2,IH))*DREAL(NHC(9,IH))*GTAU)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . -AEM**2*MUMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(21,IH))*DREAL(NHC(8,IH))*FTAU
     .   +DREAL(NHC(20,IH))*DREAL(NHC(9,IH))*GTAU)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . -AEM**2*MDMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(12,IH))*DREAL(NHC(8,IH))*FTAU
     .   +DREAL(NHC(11,IH))*DREAL(NHC(9,IH))*GTAU)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . -AEM**2*MSMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(15,IH))*DREAL(NHC(8,IH))*FTAU
     .   +DREAL(NHC(14,IH))*DREAL(NHC(9,IH))*GTAU)
      ENDDO

      IEDM_SUB_HIGGS=16 ! stau
      DO IH=1,3
       Z_HiggsEDM=STAUMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU1)
       Z_HiggsEDM=STAUMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSTAU2)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . +AEM*ME_H/32.D0/PI**3*DREAL(NHC(3,IH))/HMASS(IH)**2
     .  *DREAL(NHC(79,IH)*FSTAU1+NHC(82,IH)*FSTAU2)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . +AEM*MUMT_H/32.D0/PI**3*DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(79,IH)*FSTAU1+NHC(82,IH)*FSTAU2)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . +AEM*MDMT_H/32.D0/PI**3*DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(79,IH)*FSTAU1+NHC(82,IH)*FSTAU2)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . +AEM*MSMT_H/32.D0/PI**3*DREAL(NHC(15,IH))/HMASS(IH)**2
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
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . -AEM**2*ME_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(DREAL(NHC(3,IH))*DREAL(NHC(59,IH))*FC1
     .   +DREAL(NHC(2,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*ME_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(DREAL(NHC(3,IH))*DREAL(NHC(68,IH))*FC2
     .   +DREAL(NHC(2,IH))*DREAL(NHC(69,IH))*GC2)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . -AEM**2*MUMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(DREAL(NHC(21,IH))*DREAL(NHC(59,IH))*FC1
     .   +DREAL(NHC(20,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*MUMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(DREAL(NHC(21,IH))*DREAL(NHC(68,IH))*FC2
     .   +DREAL(NHC(20,IH))*DREAL(NHC(69,IH))*GC2)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . -AEM**2*MDMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(DREAL(NHC(12,IH))*DREAL(NHC(59,IH))*FC1
     .   +DREAL(NHC(11,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*MDMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(DREAL(NHC(12,IH))*DREAL(NHC(68,IH))*FC2
     .   +DREAL(NHC(11,IH))*DREAL(NHC(69,IH))*GC2)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . -AEM**2*MSMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(1)
     .  *(DREAL(NHC(15,IH))*DREAL(NHC(59,IH))*FC1
     .   +DREAL(NHC(14,IH))*DREAL(NHC(60,IH))*GC1)
     . -AEM**2*MSMT_H/4.D0/DSQRT(2.D0)/PI**2/SW_H**2/MW_H/M_C(2)
     .  *(DREAL(NHC(15,IH))*DREAL(NHC(68,IH))*FC2
     .   +DREAL(NHC(14,IH))*DREAL(NHC(69,IH))*GC2)
      ENDDO

*JSL[2009.Feb.27]: Charged-Higgs loops included
      IEDM_SUB_HIGGS=18 ! charged Higgs
      DO IH=1,3
       Z_HiggsEDM=MCH**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FCH)
       DEOE_E(IEDM_SUB_HIGGS)=DEOE_E(IEDM_SUB_HIGGS)
     . +AEM*ME_H/32.D0/PI**3*DREAL(NHC(3,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
       DEOE_U(IEDM_SUB_HIGGS)=DEOE_U(IEDM_SUB_HIGGS)
     . +AEM*MUMT_H/32.D0/PI**3*DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
       DEOE_D(IEDM_SUB_HIGGS)=DEOE_D(IEDM_SUB_HIGGS)
     . +AEM*MDMT_H/32.D0/PI**3*DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
       DEOE_S(IEDM_SUB_HIGGS)=DEOE_S(IEDM_SUB_HIGGS)
     . +AEM*MSMT_H/32.D0/PI**3*DREAL(NHC(15,IH))/HMASS(IH)**2
     .  *DREAL(NHC(86,IH)*FCH)
      ENDDO

*JSL[2009.Apr.22]: Flip the signs of the fermionic contributions to Barr-Zee EDM
      DEOE_E(11) =-DEOE_E(11)
      DEOE_E(12) =-DEOE_E(12)
      DEOE_E(15) =-DEOE_E(15)
      DEOE_E(17) =-DEOE_E(17)
       DEOE_U(11)=-DEOE_U(11)
       DEOE_U(12)=-DEOE_U(12)
       DEOE_U(15)=-DEOE_U(15)
       DEOE_U(17)=-DEOE_U(17)
      DEOE_D(11) =-DEOE_D(11)
      DEOE_D(12) =-DEOE_D(12)
      DEOE_D(15) =-DEOE_D(15)
      DEOE_D(17) =-DEOE_D(17)
       DEOE_S(11)=-DEOE_S(11)
       DEOE_S(12)=-DEOE_S(12)
       DEOE_S(15)=-DEOE_S(15)
       DEOE_S(17)=-DEOE_S(17)


*      do i=11,17
*       print*,i,deoe_e(i),deoe_u(i),deoe_d(i),deoe_s(i)
*      enddo
      DEOE_E(IEDM_SUB)=DEOE_E(11)+DEOE_E(12)+DEOE_E(13)+DEOE_E(14)
     .                +DEOE_E(15)+DEOE_E(16) ! remove this for old result in HiggsEDM
     .                +DEOE_E(17)+DEOE_E(18)
      DEOE_U(IEDM_SUB)=DEOE_U(11)+DEOE_U(12)+DEOE_U(13)+DEOE_U(14)
     .                +DEOE_U(15)+DEOE_U(16)
     .                +DEOE_U(17)+DEOE_U(18)
      DEOE_D(IEDM_SUB)=DEOE_D(11)+DEOE_D(12)+DEOE_D(13)+DEOE_D(14)
     .                +DEOE_D(15)+DEOE_D(16)
     .                +DEOE_D(17)+DEOE_D(18)
      DEOE_S(IEDM_SUB)=DEOE_S(11)+DEOE_S(12)+DEOE_S(13)+DEOE_S(14)
     .                +DEOE_S(15)+DEOE_S(16)
     .                +DEOE_S(17)+DEOE_S(18)
*JSL[2008.Sep.24]: taking account the relative electric charge to electron -Q_f
      DEOE_U(IEDM_SUB)=-2.D0/3.D0*DEOE_U(IEDM_SUB) 
      DEOE_D(IEDM_SUB)=DEOE_D(IEDM_SUB)/3.D0
      DEOE_S(IEDM_SUB)=DEOE_S(IEDM_SUB)/3.D0

*(d^C_u)^H, (d^C_d)^H:
*---------------------

*      print*,asmt,gsmt,mumt_h

      IEDM_SUB_HIGGS=11 ! top
      DO IH=1,3 
       Z_HiggsEDM=MTMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FTOP)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GTOP)
       DC_U(IEDM_SUB_HIGGS)=DC_U(IEDM_SUB_HIGGS)
     .  -(ASMT*GSMT/2.D0)*AEM*MUMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(21,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(20,IH))*DREAL(NHC(27,IH))*GTOP)
       DC_D(IEDM_SUB_HIGGS)=DC_D(IEDM_SUB_HIGGS)
     .  -(ASMT*GSMT/2.D0)*AEM*MDMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC(12,IH))*DREAL(NHC(26,IH))*FTOP
     .    +DREAL(NHC(11,IH))*DREAL(NHC(27,IH))*GTOP)
      ENDDO

      IEDM_SUB_HIGGS=12 ! bottom
      DO IH=1,3 
       Z_HiggsEDM=MBMT_H**2/HMASS(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX,FBOT)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX,GBOT)
       DC_U(IEDM_SUB_HIGGS)=DC_U(IEDM_SUB_HIGGS)
     . -(ASMT*GSMT/2.D0)*AEM*MUMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(21,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(20,IH))*DREAL(NHC(18,IH))*GBOT)
       DC_D(IEDM_SUB_HIGGS)=DC_D(IEDM_SUB_HIGGS)
     . -(ASMT*GSMT/2.D0)*AEM*MDMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC(12,IH))*DREAL(NHC(17,IH))*FBOT
     .   +DREAL(NHC(11,IH))*DREAL(NHC(18,IH))*GBOT)
      ENDDO

      IEDM_SUB_HIGGS=13 ! stop
      DO IH=1,3 
       Z_HiggsEDM=STMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST1)
       Z_HiggsEDM=STMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FST2)
       DC_U(IEDM_SUB_HIGGS)=DC_U(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MUMT_H/32.D0/PI**3
     .  *DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
       DC_D(IEDM_SUB_HIGGS)=DC_D(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MDMT_H/32.D0/PI**3
     .  *DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(71,IH)*FST1+NHC(74,IH)*FST2)
      ENDDO

      IEDM_SUB_HIGGS=14 ! sbottom
      DO IH=1,3
       Z_HiggsEDM=SBMASS(1)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB1)
       Z_HiggsEDM=SBMASS(2)**2/HMASS(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX,FSB2)
       DC_U(IEDM_SUB_HIGGS)=DC_U(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MUMT_H/32.D0/PI**3
     .  *DREAL(NHC(21,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
       DC_D(IEDM_SUB_HIGGS)=DC_D(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MDMT_H/32.D0/PI**3
     .  *DREAL(NHC(12,IH))/HMASS(IH)**2
     .  *DREAL(NHC(75,IH)*FSB1+NHC(78,IH)*FSB2)
      ENDDO

*JSL[2009.Apr.22]: Flip the signs of the fermionic contributions to Barr-Zee EDM
      DC_U(11) =-DC_U(11)
      DC_U(12) =-DC_U(12)
       DC_D(11)=-DC_D(11)
       DC_D(12)=-DC_D(12)

      DC_U(IEDM_SUB)=DC_U(11)+DC_U(12)+DC_U(13)+DC_U(14)
      DC_D(IEDM_SUB)=DC_D(11)+DC_D(12)+DC_D(13)+DC_D(14)
*JSL[2008.Sep.24]: taking account the relative sign to
*                  the Higgs-mediated electron EDM
      DC_U(IEDM_SUB)=-DC_U(IEDM_SUB)
      DC_D(IEDM_SUB)=-DC_D(IEDM_SUB)
*
*JSL[2010.Mar.25]: More EDMs
*
      CALL MORE_EDMS(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     . ,MCH,HMASS,OMIX
     . ,STMASS,STMIX,SBMASS,SBMIX,STAUMASS,STAUMIX,SNU3MASS
     . ,M_C,C_L,C_R,M_N,N_N,NCMAX,NHC,SHC,CHC)
*
*TOTAL:
*======
      IEDM_SUB=1
      DEOE_E(IEDM_SUB)=DEOE_E(2)+DEOE_E(3)+DEOE_E(4)+DEOE_E(5)
      DEOE_U(IEDM_SUB)=DEOE_U(2)+DEOE_U(3)+DEOE_U(4)+DEOE_U(5)
      DEOE_D(IEDM_SUB)=DEOE_D(2)+DEOE_D(3)+DEOE_D(4)+DEOE_D(5)
      DEOE_S(IEDM_SUB)=DEOE_S(2)+DEOE_S(3)+DEOE_S(4)+DEOE_S(5)
        DC_U(IEDM_SUB)=  DC_U(2)+  DC_U(3)+  DC_U(4)+  DC_U(5)
        DC_D(IEDM_SUB)=  DC_D(2)+  DC_D(3)+  DC_D(4)+  DC_D(5)
*JSL[2010.Mar.25]: Adding the more EDMs
      DEOE_E(IEDM_SUB)=DEOE_E(IEDM_SUB)
     .                +(RAUX_H(I_DEOE_E+5)+RAUX_H(I_DEOE_E+6))/GEVTOCM
      DEOE_D(IEDM_SUB)=DEOE_D(IEDM_SUB)
     .                +(RAUX_H(I_DEOE_D+5)+RAUX_H(I_DEOE_D+6))/GEVTOCM
      DEOE_S(IEDM_SUB)=DEOE_S(IEDM_SUB)
     .                +(RAUX_H(I_DEOE_S+5)+RAUX_H(I_DEOE_S+6))/GEVTOCM
*
*STORE:Electric Elecron EDM:
      RAUX_H(I_DEOE_E+0)=DEOE_E(1)*GEVTOCM
      RAUX_H(I_DEOE_E+1)=DEOE_E(2)*GEVTOCM
      RAUX_H(I_DEOE_E+2)=DEOE_E(3)*GEVTOCM
      RAUX_H(I_DEOE_E+3)=DEOE_E(4)*GEVTOCM
      RAUX_H(I_DEOE_E+4)=DEOE_E(5)*GEVTOCM
*STORE:Electric Up-quark EDM:
      RAUX_H(I_DEOE_U+0)=DEOE_U(1)*GEVTOCM
      RAUX_H(I_DEOE_U+1)=DEOE_U(2)*GEVTOCM
      RAUX_H(I_DEOE_U+2)=DEOE_U(3)*GEVTOCM
      RAUX_H(I_DEOE_U+3)=DEOE_U(4)*GEVTOCM
      RAUX_H(I_DEOE_U+4)=DEOE_U(5)*GEVTOCM
*STORE:Electric Down-quark EDM:
      RAUX_H(I_DEOE_D+0)=DEOE_D(1)*GEVTOCM
      RAUX_H(I_DEOE_D+1)=DEOE_D(2)*GEVTOCM
      RAUX_H(I_DEOE_D+2)=DEOE_D(3)*GEVTOCM
      RAUX_H(I_DEOE_D+3)=DEOE_D(4)*GEVTOCM
      RAUX_H(I_DEOE_D+4)=DEOE_D(5)*GEVTOCM
*STORE:Electric Strange-quark EDM:
      RAUX_H(I_DEOE_S+0)=DEOE_S(1)*GEVTOCM
      RAUX_H(I_DEOE_S+1)=DEOE_S(2)*GEVTOCM
      RAUX_H(I_DEOE_S+2)=DEOE_S(3)*GEVTOCM
      RAUX_H(I_DEOE_S+3)=DEOE_S(4)*GEVTOCM
      RAUX_H(I_DEOE_S+4)=DEOE_S(5)*GEVTOCM
*STORE:Chromo-Electric Up-quark EDM:
      RAUX_H(I_DC_U+0)=DC_U(1)*GEVTOCM
      RAUX_H(I_DC_U+1)=DC_U(2)*GEVTOCM
      RAUX_H(I_DC_U+2)=DC_U(3)*GEVTOCM
      RAUX_H(I_DC_U+3)=DC_U(4)*GEVTOCM
      RAUX_H(I_DC_U+4)=DC_U(5)*GEVTOCM
*STORE:Chromo-Electric Down-quark EDM:
      RAUX_H(I_DC_D+0)=DC_D(1)*GEVTOCM
      RAUX_H(I_DC_D+1)=DC_D(2)*GEVTOCM
      RAUX_H(I_DC_D+2)=DC_D(3)*GEVTOCM
      RAUX_H(I_DC_D+3)=DC_D(4)*GEVTOCM
      RAUX_H(I_DC_D+4)=DC_D(5)*GEVTOCM
*
      IF(IFLAG(18).EQ.2) THEN
      print*,'---------------------------------------------------------'
      print*,'The Electric EDMs of particles in cm: e, u, d, s: '
      print*,'---------------------------------------------------------'
      write(*,8) 'd^E_e/e[Total]:',deoe_e(1)*GEVTOCM
      write(*,8) 'd^E_u/e[Total]:',deoe_u(1)*GEVTOCM
      write(*,8) 'd^E_d/e[Total]:',deoe_d(1)*GEVTOCM
      write(*,8) 'd^E_s/e[Total]:',deoe_s(1)*GEVTOCM
      write(*,8) 'd^E_e/e[C,N,Gl,H]:'
     .      ,deoe_e(2)*GEVTOCM,deoe_e(3)*GEVTOCM
     .      ,deoe_e(4)*GEVTOCM,deoe_e(5)*GEVTOCM
      write(*,8) 'd^E_u/e[C,N,Gl,H]:'
     .      ,deoe_u(2)*GEVTOCM,deoe_u(3)*GEVTOCM
     .      ,deoe_u(4)*GEVTOCM,deoe_u(5)*GEVTOCM
      write(*,8) 'd^E_d/e[C,N,Gl,H]:'
     .      ,deoe_d(2)*GEVTOCM,deoe_d(3)*GEVTOCM
     .      ,deoe_d(4)*GEVTOCM,deoe_d(5)*GEVTOCM
      write(*,8) 'd^E_s/e[C,N,Gl,H]:'
     .      ,deoe_s(2)*GEVTOCM,deoe_s(3)*GEVTOCM
     .      ,deoe_s(4)*GEVTOCM,deoe_s(5)*GEVTOCM
      print*,'---------------------------------------------------------'
      print*,'The MORE Electric EDMs of particles in cm: e, d, s: '
      print*,'---------------------------------------------------------'
      write(*,8) 'd^E_e/e[WH,WW]:'
     .          ,RAUX_H(I_DEOE_E+5),RAUX_H(I_DEOE_E+6)
      write(*,8) 'd^E_d/e[WH,WW]:'
     .          ,RAUX_H(I_DEOE_D+5),RAUX_H(I_DEOE_D+6)
      write(*,8) 'd^E_s/e[WH,WW]:'
     .          ,RAUX_H(I_DEOE_S+5),RAUX_H(I_DEOE_S+6)
      print*,'---------------------------------------------------------'
      print*,'The Chromo-Electric EDMs of particles in cm: u, d, s: '
      print*,'---------------------------------------------------------'
      write(*,8) 'd^C_u  [Total]:',dc_u(1)*GEVTOCM
      write(*,8) 'd^C_d  [Total]:',dc_d(1)*GEVTOCM
      write(*,8) 'd^C_s  [Total]:',RAUX_H(I_DC_S+0)
      write(*,8) 'd^C_u  [C,N,Gl,H]:' 
     .      ,dc_u(2)*GEVTOCM,dc_u(3)*GEVTOCM
     .      ,dc_u(4)*GEVTOCM,dc_u(5)*GEVTOCM
      write(*,8) 'd^C_d  [C,N,Gl,H]:'
     .      ,dc_d(2)*GEVTOCM,dc_d(3)*GEVTOCM
     .      ,dc_d(4)*GEVTOCM,dc_d(5)*GEVTOCM
      write(*,8) 'd^C_s  [C,N,Gl,H]:'
     .          ,RAUX_H(I_DC_S+1),RAUX_H(I_DC_S+2)
     .          ,RAUX_H(I_DC_S+3),RAUX_H(I_DC_S+4)
      ENDIF ! IF(IFLAG(18).EQ.2) THEN
*-----------------------------------------------------------------------
*Purely-gluonic D-6 Weinberg operator @ Electro-weak (or mt^pole) scale
*
       DO I=1,NEDM_SUB
        DG_WEINBERG(I)=0.D0
       ENDDO
*Higggs:
*=======
*       Z_HiggsEDM=1.D-6
*       CALL BODE(EDM_H,X1D,X1U,NX,HXXX)
*       print*,Z_HiggsEDM,HXXX,' simeq 0.0625 = 1/16 ?'
      IEDM_WEINBERG=2
      DO IH=1,3
       Z_HiggsEDM=HMASS(IH)**2/MTMT_H**2
       CALL BODE(EDM_H,X1D,X1U,NX,HTOP)
*        print*,Z_HiggsEDM,HTOP
       Z_HiggsEDM=HMASS(IH)**2/MBMT_H**2
       CALL BODE(EDM_H,X1D,X1U,NX,HBOT)
*        print*,Z_HiggsEDM,HBOT
       DG_WEINBERG(IEDM_WEINBERG)=DG_WEINBERG(IEDM_WEINBERG)
     .   +4.D0*DSQRT(2.D0)*GF_H*GSMT**3/(4.D0*PI)**4
     .   *(DREAL(NHC(26,IH))*DREAL(NHC(27,IH))*HTOP
     .    +DREAL(NHC(17,IH))*DREAL(NHC(18,IH))*HBOT)
      ENDDO
*
*Gluino:
*=======
      IEDM_WEINBERG=3
*If you are providing the subroutine for the calculation of the 
*loop function H(msq1^2/|M3|^2,msq2^2/|M3|^2,mq^2/|M3|^2) for 
*the gluino contribution to the Weinberg operator, use:
*      CALL EDM_HH_USER(STMASS,SBMASS,CDABS(M3_H),MBMT_H,MTMT_H
*     .                ,EDM_HH_STOP,EDM_HH_SBOT)
*If not, use:
      CALL EDM_HH_USER0(STMASS,SBMASS,CDABS(M3_H),MBMT_H,MTMT_H
     .                 ,EDM_HH_STOP,EDM_HH_SBOT)
*
*      print*,'H(stop), H(sbottom)',edm_hh_stop,edm_hh_sbot

      DG_WEINBERG(IEDM_WEINBERG)=
     .-3.D0/2.D0/PI*(GSMT/4.D0/PI/CDABS(M3_H))**3
     .*( 
     .   MTMT_H*(STMASS(1)**2/CDABS(M3_H)**2 ! stop
     .          *DIMAG(DCONJG(GR_GL_T_STI(1))*GL_GL_T_STI(1))
     .          +STMASS(2)**2/CDABS(M3_H)**2
     .          *DIMAG(DCONJG(GR_GL_T_STI(2))*GL_GL_T_STI(2))
     .          )*EDM_HH_STOP
     .  +MBMT_H*(SBMASS(1)**2/CDABS(M3_H)**2 ! sbottom
     .          *DIMAG(DCONJG(GR_GL_B_SBI(1))*GL_GL_B_SBI(1))
     .          +SBMASS(2)**2/CDABS(M3_H)**2
     .          *DIMAG(DCONJG(GR_GL_B_SBI(2))*GL_GL_B_SBI(2))
     .          )*EDM_HH_SBOT
     . )
*
*TOTAL:
*======
      DG_WEINBERG(1)=DG_WEINBERG(2)+DG_WEINBERG(3)
*
*STORE:Purely-gluonic D-6 Weinberg operator:
      RAUX_H(I_DG_WEINBERG+0)=DG_WEINBERG(1)*GEVTOCM
      RAUX_H(I_DG_WEINBERG+1)=DG_WEINBERG(2)*GEVTOCM
      RAUX_H(I_DG_WEINBERG+2)=DG_WEINBERG(3)*GEVTOCM
*
      IF(IFLAG(18).EQ.2) THEN
      print*,'---------------------------------------------------------'
      print*,'Purely-gluonic D-6 Weinberg operator in cm/GeV: '
      print*,'---------------------------------------------------------'
      write(*,11) 'd^G       [Total]:',dg_weinberg(1)*GEVTOCM
      write(*,11) 'd^G[Higgs,Gluino]:',dg_weinberg(2)*GEVTOCM
     .                                ,dg_weinberg(3)*GEVTOCM
      ENDIF ! IF(IFLAG(18).EQ.2) THEN
*-----------------------------------------------------------------------
*Four-fermion couplings at M_F(P) and/or 1 GeV when M_F(P)<1 GeV
*      COMPLEX*16 C4_F_FP
*
*Quark masses at arbitrary scale RSCALE:
*
      MTPOLE=MTPOLE_H
      MBPOLE=RAUX_H(1)
      MCPOLE=RAUX_H(4)
*      print*,mtpole,mbpole,mcpole
      AS_MT =ASMT_H
      AS_MZ =ASMZ_H
      AS_MB =RAUX_H(3)
      AS_MC =RAUX_H(6)
*      print*,as_mt,as_mz,as_mb,as_mc

      RSCALE=MTPOLE
      CALL MQ_RUN(RSCALE,MTPOLE,MZ_H,MBPOLE,MCPOLE
     .           ,MTMT_H,MBMT_H,MCMT_H,MSMT_H,MUMT_H,MDMT_H
     .           ,AS_MT,AS_MZ,AS_MB,AS_MC
     .           ,MT_MTPOLE,MB_MTPOLE,MC_MTPOLE
     .                               ,MS_MTPOLE,MU_MTPOLE,MD_MTPOLE)
*      print*,MT_MTPOLE,MB_MTPOLE,MC_MTPOLE,MS_MTPOLE,MU_MTPOLE,MD_MTPOLE

      RSCALE=MBPOLE
      CALL MQ_RUN(RSCALE,MTPOLE,MZ_H,MBPOLE,MCPOLE
     .           ,MTMT_H,MBMT_H,MCMT_H,MSMT_H,MUMT_H,MDMT_H
     .           ,AS_MT,AS_MZ,AS_MB,AS_MC
     .           ,MT_MBPOLE,MB_MBPOLE,MC_MBPOLE
     .                               ,MS_MBPOLE,MU_MBPOLE,MD_MBPOLE)
*      print*,MT_MBPOLE,MB_MBPOLE,MC_MBPOLE,MS_MBPOLE,MU_MBPOLE,MD_MBPOLE

      RSCALE=MCPOLE
      CALL MQ_RUN(RSCALE,MTPOLE,MZ_H,MBPOLE,MCPOLE
     .           ,MTMT_H,MBMT_H,MCMT_H,MSMT_H,MUMT_H,MDMT_H
     .           ,AS_MT,AS_MZ,AS_MB,AS_MC
     .           ,MT_MCPOLE,MB_MCPOLE,MC_MCPOLE
     .                               ,MS_MCPOLE,MU_MCPOLE,MD_MCPOLE)
*      print*,MT_MCPOLE,MB_MCPOLE,MC_MCPOLE,MS_MCPOLE,MU_MCPOLE,MD_MCPOLE

      RSCALE=1.D0
      CALL MQ_RUN(RSCALE,MTPOLE,MZ_H,MBPOLE,MCPOLE
     .           ,MTMT_H,MBMT_H,MCMT_H,MSMT_H,MUMT_H,MDMT_H
     .           ,AS_MT,AS_MZ,AS_MB,AS_MC
     .           ,MT_1,MB_1,MC_1,MS_1,MU_1,MD_1)
*      print*,MT_1,MB_1,MC_1,MS_1,MU_1,MD_1
*
*C4_de
      N_F =10
      N_FP=1
      M_F =MD_1
      M_FP=ME_H
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_DE=DREAL(C4_F_FP)
*      print*,'C4_de:',c4_f_fp,c4_de
*C4_se
      N_F =13
      N_FP=1
      M_F =MS_1
      M_FP=ME_H
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_SE=DREAL(C4_F_FP)
*      print*,'C4_se:',c4_f_fp,c4_se
*C4_ed
      N_F =1
      N_FP=10
      M_F =ME_H
      M_FP=MD_1
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_ED=DREAL(C4_F_FP)
*      print*,'C4_ed:',c4_f_fp,c4_ed
*C4_es
      N_F =1
      N_FP=13
      M_F =ME_H
      M_FP=MS_1
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_ES=DREAL(C4_F_FP)
*      print*,'C4_es:',c4_f_fp,c4_es
*C4_eb
      N_F =1
      N_FP=16
      M_F =ME_H
      M_FP=MB_MBPOLE
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_EB=DREAL(C4_F_FP)
*      print*,'C4_eb:',c4_f_fp,c4_eb
*C4_ec
      N_F =1
      N_FP=22
      M_F =ME_H
      M_FP=MC_MCPOLE
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_EC=DREAL(C4_F_FP)
*      print*,'C4_ec:',c4_f_fp,c4_ec
*C4_et
      N_F =1
      N_FP=25
      M_F =ME_H
      M_FP=MT_MTPOLE
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_ET=DREAL(C4_F_FP)
*      print*,'C4_et:',c4_f_fp,c4_et
*C4_dd
      N_F =10
      N_FP=10
      M_F =MD_1
      M_FP=MD_1
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_DD=DREAL(C4_F_FP)
*      print*,'C4_dd:',c4_f_fp,c4_dd
*C4_sd
      N_F =13
      N_FP=10
      M_F =MS_1
      M_FP=MD_1
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_SD=DREAL(C4_F_FP)
*      print*,'C4_sd:',c4_f_fp,c4_sd
*C4_bd
      N_F =16
      N_FP=10
      M_F =MB_MBPOLE
      M_FP=MD_1
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_BD=DREAL(C4_F_FP)
*      print*,'C4_bd:',c4_f_fp,c4_bd
*C4_db
      N_F =10
      N_FP=16
      M_F =MD_1
      M_FP=MB_MBPOLE
      CALL FOUR_FERMION(N_F,N_FP,V_H,M_F,M_FP,HMASS,NCMAX,NHC,C4_F_FP)
      C4_DB=DREAL(C4_F_FP)
*      print*,'C4_db:',c4_f_fp,c4_db
*
*STORE:C_4f
      RAUX_H(I_C4FOM+ 0)=C4_de/MD_1*GEVTOCM
      RAUX_H(I_C4FOM+ 1)=C4_se/MS_1*GEVTOCM
      RAUX_H(I_C4FOM+ 2)=C4_ed/MD_1*GEVTOCM
      RAUX_H(I_C4FOM+ 3)=C4_es/MS_1*GEVTOCM
      RAUX_H(I_C4FOM+ 4)=C4_eb/MB_MBPOLE*GEVTOCM
      RAUX_H(I_C4FOM+ 5)=C4_ec/MC_MCPOLE*GEVTOCM
      RAUX_H(I_C4FOM+ 6)=C4_et/MT_MTPOLE*GEVTOCM
      RAUX_H(I_C4FOM+ 7)=C4_dd/MD_1*GEVTOCM
      RAUX_H(I_C4FOM+ 8)=C4_sd/MS_1*GEVTOCM
      RAUX_H(I_C4FOM+ 9)=C4_bd/MB_MBPOLE*GEVTOCM
      RAUX_H(I_C4FOM+10)=C4_db/MB_MBPOLE*GEVTOCM
*
      IF(IFLAG(18).EQ.2) THEN
      print*,'---------------------------------------------------------'
      print*,' Four-fermion couplings needed for EDMs in cm/GeV^2:'
      print*,'---------------------------------------------------------'
      write(*,7) 'C4_de/m_d: ',C4_de/MD_1*GEVTOCM
      write(*,7) 'C4_se/m_s: ',C4_se/MS_1*GEVTOCM
      write(*,7) 'C4_ed/m_d: ',C4_ed/MD_1*GEVTOCM
      write(*,7) 'C4_es/m_s: ',C4_es/MS_1*GEVTOCM
      write(*,7) 'C4_eb/m_b: ',C4_eb/MB_MBPOLE*GEVTOCM
      write(*,7) 'C4_ec/m_c: ',C4_ec/MC_MCPOLE*GEVTOCM
      write(*,7) 'C4_et/m_t: ',C4_et/MT_MTPOLE*GEVTOCM
      write(*,7) 'C4_dd/m_d: ',C4_dd/MD_1*GEVTOCM
      write(*,7) 'C4_sd/m_s: ',C4_sd/MS_1*GEVTOCM
      write(*,7) 'C4_bd/m_b: ',C4_bd/MB_MBPOLE*GEVTOCM
      write(*,7) 'C4_db/m_b: ',C4_db/MB_MBPOLE*GEVTOCM
      ENDIF ! IF(IFLAG(18).EQ.2) THEN
*-----------------------------------------------------------------------
*CS_TOT, CP_TOT, CPP_TOT
      XKAPPA=0.5D0 ! \pm 0.25
      X_TOP=1.D0
      X_BOT=1.D0-XKAPPA*0.25D0
*      x_bot=1.d0 ! for old result in HiggsEDM
*
      CS_G =0.D0
      DO IH=1,3
       CS_G =CS_G
     .      +0.1D0*ME_H/V_H**2*DREAL(NHC(3,IH))/HMASS(IH)**2
     .      *( 2.D0/3.D0*X_TOP*DREAL(NHC(26,IH))              ! top
     .        +2.D0/3.D0*X_BOT*DREAL(NHC(17,IH))              ! bottom
     .        -V_H**2/12.D0*( DREAL(NHC(71,IH))/STMASS(1)**2  ! stop
     .                       +DREAL(NHC(74,IH))/STMASS(2)**2
     .                       +DREAL(NHC(75,IH))/SBMASS(1)**2  ! sbottom
     .                       +DREAL(NHC(78,IH))/SBMASS(2)**2 ) )
      ENDDO
      CS_4F=29.D-3*C4_DE/MD_1+XKAPPA*220.D-3*C4_SE/MS_1
      CS_TOT=CS_4F+CS_G
*      cs_tot=cs_g*128.d0/137.d0 ! for old result in HiggsEDM
*
      CP_4F=-375.D-3*( C4_EC/MC_MCPOLE+C4_ES/MS_1
     .                 +C4_ET/MT_MTPOLE+C4_EB/MB_MBPOLE )
      CP_TOT=CP_4F
*
      CPP_4F=-806.D-3*C4_ED/MD_1
     .       -181.D-3*( C4_EC/MC_MCPOLE+C4_ES/MS_1
     .                  +C4_ET/MT_MTPOLE+C4_EB/MB_MBPOLE )
      CPP_TOT=CPP_4F
*
*STORE:C's
      RAUX_H(I_CSPP+0)=CS_TOT*GEVTOCM
      RAUX_H(I_CSPP+1)=CP_TOT*GEVTOCM
      RAUX_H(I_CSPP+2)=CPP_TOT*GEVTOCM
*
      IF(IFLAG(18).EQ.2) THEN
      print*,'---------------------------------------------------------'
      print*,'   C_S, C_P, and C_P^prime in cm/GeV and in 1/GeV^2: '
      print*,'---------------------------------------------------------'
      write(*,7) 'C_S      : ',CS_tot*GEVTOCM,CS_tot
      write(*,7) 'C_P      : ',CP_tot*GEVTOCM,CP_tot
      write(*,7) 'C_P^prime: ',CPP_tot*GEVTOCM,CPP_tot
      print*,'---------------------------------------------------------'
      ENDIF ! IF(IFLAG(18).EQ.2) THEN
*-----------------------------------------------------------------------
*Thallium EDM d^Tl/(e cm)
      DTL=-585.D0*DEOE_E(1)*GEVTOCM-8.5D-13*CS_TOT
*STORE
      RAUX_H(I_TL+0)=DTL
      RAUX_H(I_TL+1)=-585.D0*DEOE_E(1)*GEVTOCM
      RAUX_H(I_TL+2)=-8.5D-13*CS_TOT
      IF(IFLAG(18).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'       Thallium EDM in units of [e cm]: d^Tl/[e cm] '
      print*,'---------------------------------------------------------'
      write(*,9) 'd^Tl/(e cm) [Total]=',dtl
      print*,' Each contribution to d^Tl from'
      write(*,9) '            [d^E_e]=',-585.D0*DEOE_E(1)*GEVTOCM
      write(*,9) '            [C_S  ]=',-8.5D-13*CS_TOT
      ENDIF ! IF(IFLAG(18).EQ.1) THEN
*      print*,' '
*      print*,' Each contribution to d^Tl through d^E_e from'
*      write(*,10)'[Charginos   ]=',-585.D0*DEOE_E(2)*GEVTOCM
*      write(*,10)'[Neutralinos ]=',-585.D0*DEOE_E(3)*GEVTOCM
*      write(*,10)'[Gluinos     ]=',-585.D0*DEOE_E(4)*GEVTOCM
*      write(*,10)'[Higgs bosons]=',-585.D0*DEOE_E(5)*GEVTOCM
*-----------------------------------------------------------------------
*Neutron EDM d^n/(e cm)
*
*(1) Chiral Quark Model
      ETAE=1.53D0
      ETAC=3.4D0
      ETAG=3.4D0
      XLAM=1.19D0  ! in GeV
      DOE_U=ETAE*DEOE_U(1)+ETAC/4.D0/PI*DC_U(1)
     .     +ETAG*XLAM/4.D0/PI*DG_WEINBERG(1)
      DOE_D=ETAE*DEOE_D(1)+ETAC/4.D0/PI*DC_D(1)
     .     +ETAG*XLAM/4.D0/PI*DG_WEINBERG(1)
      DN_CQM=(4.D0*DOE_D-DOE_U)/3.D0*GEVTOCM
*(2) Parton Quark Model
      SPIN_U=-0.508D0
      SPIN_D= 0.746D0
      SPIN_S=-0.226D0
      DN_PQM=ETAE*(SPIN_U*DEOE_U(1)
     .            +SPIN_D*DEOE_D(1)
     .            +SPIN_S*DEOE_S(1))*GEVTOCM
*(3) QCD sum rule
      DN_QCD1=1.4D0*(DEOE_D(1)-0.25D0*DEOE_U(1))*GEVTOCM  ! d^E 
      DN_QCD2=1.1D0*(DC_D(1)+0.5D0*DC_U(1))/GSMT*GEVTOCM  ! d^C
       DG_EWTO1GEV=ETAG/0.4D0
      DN_QCD3=DG_EWTO1GEV*(20.D-3*DG_WEINBERG(1))*GEVTOCM ! d^G(1GeV)
      DN_QCD4=(2.6D-3*(C4_bd+0.75D0*C4_db)/MB_MBPOLE)*GEVTOCM ! C_bd
      DN_QCD=DN_QCD1+DN_QCD2+DN_QCD3+DN_QCD4
*
*STORE
      RAUX_H(I_N1+0)=DN_CQM
      RAUX_H(I_N1+1)=ETAE*(4.D0*DEOE_D(1)-DEOE_U(1))/3.D0*GEVTOCM
      RAUX_H(I_N1+2)=ETAC/4.D0/PI*(4.D0*DC_D(1)-DC_U(1))/3.D0*GEVTOCM
      RAUX_H(I_N1+3)=ETAG*XLAM/4.D0/PI*DG_WEINBERG(1)*GEVTOCM
      RAUX_H(I_N2+0)=DN_PQM
      RAUX_H(I_N2+1)=ETAE*SPIN_U*DEOE_U(1)*GEVTOCM
      RAUX_H(I_N2+2)=ETAE*SPIN_D*DEOE_D(1)*GEVTOCM
      RAUX_H(I_N2+3)=ETAE*SPIN_S*DEOE_S(1)*GEVTOCM
      RAUX_H(I_N3+0)=DN_QCD
      RAUX_H(I_N3+1)=DN_QCD1
      RAUX_H(I_N3+2)=DN_QCD2
      RAUX_H(I_N3+3)=DN_QCD3
      RAUX_H(I_N3+4)=DN_QCD4
*
      IF(IFLAG(18).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'       Neutron  EDM in units of [e cm]: d^n/[e cm] '
      print*,'---------------------------------------------------------'
      print*,'(1) Chiral Quark Model'
      write(*,9) 'd^n/(e cm)  [Total]=',dn_cqm
      print*,' Each contribution to d^n from'
      write(*,10)'[d^E_u & d^E_d]=',ETAE*
     .           (4.D0*DEOE_D(1)-DEOE_U(1))/3.D0*GEVTOCM
      write(*,10)'[d^C_u & d^C_d]=',ETAC/4.D0/PI*
     .           (4.D0*DC_D(1)-DC_U(1))/3.D0*GEVTOCM
      write(*,10)'[ Weinberg-6D ]=',ETAG*XLAM/4.D0/PI*
     .           DG_WEINBERG(1)*GEVTOCM
      print*,' '
      print*,'(2) Parton Quark Model'
      write(*,9) 'd^n/(e cm)  [Total]=',dn_pqm
      print*,' Each contribution to d^n from'
      write(*,10)'[d^E_u        ]=',ETAE*SPIN_U*DEOE_U(1)*GEVTOCM
      write(*,10)'[d^E_d        ]=',ETAE*SPIN_D*DEOE_D(1)*GEVTOCM
      write(*,10)'[d^E_s        ]=',ETAE*SPIN_S*DEOE_S(1)*GEVTOCM
      print*,' '
      print*,'(3) QCD sum rule technique'
      write(*,9) 'd^n/(e cm)  [Total]=',dn_qcd
      print*,' Each contribution to d^n from'
      write(*,10)'[d^E_u & d^E_d]=',DN_QCD1
      write(*,10)'[d^C_u & d^C_d]=',DN_QCD2
      write(*,10)'[ Weinberg-6D ]=',DN_QCD3
      write(*,10)'[ C_bd & C_db ]=',DN_QCD4
      ENDIF ! IF(IFLAG(18).EQ.1) THEN
*-----------------------------------------------------------------------
*Mercury EDM d^Hg/(e cm)
*
      D_HG1=1.D-2*(DEOE_E(1))*GEVTOCM
      D_HG2=7.D-3*(DC_U(1)-DC_D(1))/GSMT*GEVTOCM
      D_HG3=-1.4D-5*( 0.5D0*C4_dd/MD_1
     .               +3.3D0*XKAPPA*C4_sd/MS_1
     .               +(1.D0-0.25D0*XKAPPA)*C4_bd/MB_MBPOLE )*GEVTOCM
      D_HG4=3.5D-3*CS_TOT*GEVTOCM
      D_HG5=4.0D-4*(CP_TOT-0.2D0*CPP_TOT)*GEVTOCM
      D_HG=D_HG1+D_HG2+D_HG3+D_HG4+D_HG5
*STORE
      RAUX_H(I_HG+0)=D_HG
      RAUX_H(I_HG+1)=D_HG1
      RAUX_H(I_HG+2)=D_HG2
      RAUX_H(I_HG+3)=D_HG3
      RAUX_H(I_HG+4)=D_HG4
      RAUX_H(I_HG+5)=D_HG5
*
      IF(IFLAG(18).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'       Mercury EDM in units of [e cm]: d^Hg/[e cm] '
      print*,'---------------------------------------------------------'
      write(*,9) 'd^Hg/(e cm) [Total]=',d_hg
      print*,' Each contribution to d^Hg from'
      write(*,10)'[d^E_e        ]=',D_HG1
      write(*,10)'[d^C_u & d^C_d]=',D_HG2
      write(*,10)'[C_4f         ]=',D_HG3
      write(*,10)'[C_S          ]=',D_HG4
      write(*,10)'[C_P & C_P^pr ]=',D_HG5
      ENDIF ! IF(IFLAG(18).EQ.1) THEN
*-----------------------------------------------------------------------
*Deuteron EDM d^D/(e cm)
*
      D_D1=0.5D0*(DEOE_U(1)+DEOE_D(1))*GEVTOCM
      D_D2=-(5.D0+0.6D0)*(DC_U(1)-DC_D(1))/GSMT*GEVTOCM
     .     -0.2D0*(DC_U(1)+DC_D(1))/GSMT*GEVTOCM
      D_D3=1.0D-2*( 0.5D0*C4_dd/MD_1
     .             +3.3D0*XKAPPA*C4_sd/MS_1
     .             +(1.D0-0.25D0*XKAPPA)*C4_bd/MB_MBPOLE )*GEVTOCM
       DG_EWTO1GEV=ETAG/0.4D0
      D_D4=DG_EWTO1GEV*(20.D-3*DG_WEINBERG(1))*GEVTOCM 
      D_D=D_D1+D_D2+D_D3+D_D4
*STORE
      RAUX_H(I_DEUT+0)=D_D
      RAUX_H(I_DEUT+1)=D_D1
      RAUX_H(I_DEUT+2)=D_D2
      RAUX_H(I_DEUT+3)=D_D3
      RAUX_H(I_DEUT+4)=D_D4
*
      IF(IFLAG(18).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'       Deuteron EDM in units of [e cm]: d^D/[e cm] '
      print*,'---------------------------------------------------------'
      write(*,9) 'd^D/(e cm)  [Total]=',d_d
      print*,' Each contribution to d^D from'
      write(*,10)'[d^E_u & d^E_d]=',D_D1
      write(*,10)'[d^C_u & d^C_d]=',D_D2
      write(*,10)'[C_4f         ]=',D_D3
      write(*,10)'[ Weinberg-6D ]=',D_D4
      print*,'---------------------------------------------------------'
      ENDIF ! IF(IFLAG(18).EQ.1) THEN
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
 8    FORMAT(2X,A18,1X,4(1X,E10.4,1X))
 9    FORMAT(2X,A20,1X,E10.4)
 10   FORMAT(2X,A27,1X,E10.4)
 11   FORMAT(2X,A18,1X,2(1X,E10.4,1X))
*
      RETURN
      END

      SUBROUTINE SFERMION_MIXING(XLL,XRR,XRL,XLR,SFMASS,SFMIX)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      COMPLEX*16 XRL,XLR
      REAL*8     SFMASS(2)
      COMPLEX*16 SFMIX(2,2)
*-----------------------------------------------------------------------
*Local
      COMPLEX*16 D11,D22,D12,D21
      COMPLEX*16 M1SQ,M2SQ,D121,D122
*-----------------------------------------------------------------------
      PI=2.D0*DASIN(1.D0)
*      print*,'sfermion mass squared in TeV^2:'
*      print*,'/',XLL/1.D6,XLR/1.D6,'\\'
*      print*,'\\',XRL/1.D6,XRR/1.D6,'/'
*
      DELTA = DSQRT((XLL-XRR)**2+4.D0*CDABS(XRL)**2)
      XMAVG = (XLL+XRR)/2.D0

      SFMASS(1)=DSQRT(DABS(XMAVG-DELTA/2.D0))
      IF((XMAVG-DELTA/2.D0).LT.0.D0) SFMASS(1)=-SFMASS(1)
      SFMASS(2)=DSQRT(XMAVG+DELTA/2.D0)

      IF(CDABS(XRL).EQ.0.D0) XRL=DCMPLX(1.D-10,0.D0)
      PHI=DATAN(DIMAG(XRL)/DREAL(XRL))
*      print*,'ACK',DIMAG(XRL)/DREAL(XRL),DTAN(PHI)
      THT_ABS = DATAN(-(SFMASS(1)**2-XLL)/CDABS(XRL))
      IF (DREAL(XRL).LT.0.D0) THT= THT_ABS
      IF (DREAL(XRL).GT.0.D0) THT=-THT_ABS
      IF(SFMASS(1).GT.0.D0) THEN
*      print*,'ACK',DCOS(PHI),DCOS(THT),THT_ABS
      IF(DCOS(PHI).LT.0.D0 .OR. DCOS(THT).LT.0.D0 .OR. THT_ABS.LT.0.D0)
     .   THEN
         print*,'ERROR in <SFERMION_MIXING> !!!'
     .         ,DCOS(PHI),DCOS(THT),THT_ABS
      ENDIF
      ENDIF

*     SFMIX(alpha,i)
*      print*,'ACK',PHI,DCOS(PHI),DSIN(PHI)
*      print*,'ACK',THT,DCOS(THT),DSIN(THT)
      SFMIX(1,1)=DCMPLX(DCOS(THT),0.D0)
      SFMIX(1,2)=DCMPLX(-DSIN(THT)*DCOS(PHI),DSIN(THT)*DSIN(PHI))
      SFMIX(2,1)=DCMPLX(DSIN(THT)*DCOS(PHI),DSIN(THT)*DSIN(PHI))
      SFMIX(2,2)=DCMPLX(DCOS(THT),0.D0)
*Check                 A           A
      D11=DCONJG(SFMIX(1,1))*SFMIX(1,1)
     .   +DCONJG(SFMIX(2,1))*SFMIX(2,1)
      D22=DCONJG(SFMIX(1,2))*SFMIX(1,2)
     .   +DCONJG(SFMIX(2,2))*SFMIX(2,2)
      D12=DCONJG(SFMIX(1,1))*SFMIX(1,2)
     .   +DCONJG(SFMIX(2,1))*SFMIX(2,2)
      D21=DCONJG(SFMIX(1,2))*SFMIX(1,1)
     .   +DCONJG(SFMIX(2,2))*SFMIX(2,1)
*      print*,'1,1,0,0?',d11,d22,d12,d21
*Check                  A      AB       B
      M1SQ=DCONJG(SFMIX(1,1))*XLL*SFMIX(1,1)
     .    +DCONJG(SFMIX(1,1))*XLR*SFMIX(2,1)
     .    +DCONJG(SFMIX(2,1))*XRL*SFMIX(1,1)
     .    +DCONJG(SFMIX(2,1))*XRR*SFMIX(2,1)
      M2SQ=DCONJG(SFMIX(1,2))*XLL*SFMIX(1,2)
     .    +DCONJG(SFMIX(1,2))*XLR*SFMIX(2,2)
     .    +DCONJG(SFMIX(2,2))*XRL*SFMIX(1,2)
     .    +DCONJG(SFMIX(2,2))*XRR*SFMIX(2,2)
      D121=DCONJG(SFMIX(1,1))*XLL*SFMIX(1,2)
     .    +DCONJG(SFMIX(1,1))*XLR*SFMIX(2,2)
     .    +DCONJG(SFMIX(2,1))*XRL*SFMIX(1,2)
     .    +DCONJG(SFMIX(2,1))*XRR*SFMIX(2,2)
      D122=DCONJG(SFMIX(1,2))*XLL*SFMIX(1,1)
     .    +DCONJG(SFMIX(1,2))*XLR*SFMIX(2,1)
     .    +DCONJG(SFMIX(2,2))*XRL*SFMIX(1,1)
     .    +DCONJG(SFMIX(2,2))*XRR*SFMIX(2,1)
*      print*,'Sfermion_Mixing : All zer0s [in units of GeV^2] ?'
*      write(*,4) DREAL(M1SQ-SFMASS(1)**2),DIMAG(M1SQ-SFMASS(1)**2)
*     .          ,DREAL(M2SQ-SFMASS(2)**2),DIMAG(M2SQ-SFMASS(2)**2)
*      write(*,4) DREAL(D121),DIMAG(D121),DREAL(D122),DIMAG(D122)
*-----------------------------------------------------------------------
 4    FORMAT(2X,4(1X,E10.4,1X))
*
      RETURN
      END

      SUBROUTINE MQ_RUN(SQRTS,MT_POLE,MZ,MB_POLE,MC_POLE
     .                 ,MT_MT,MB_MT,MC_MT,MS_MT,MU_MT,MD_MT
     .                 ,AS_MT,AS_MZ,AS_MB,AS_MC
     .                 ,MT_S,MB_S,MC_S,MS_S,MU_S,MD_S)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      PI      = 2.D0*DASIN(1.D0)
      B3      = (11.D0-2.D0/3.D0*3.D0)/4.D0/PI
      B4      = (11.D0-2.D0/3.D0*4.D0)/4.D0/PI
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
*
*Quark masses at mb^pole and mc^pole
      MT_MB = MT_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MB_MB = MB_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MC_MB = MC_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MS_MB = MS_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MU_MB = MU_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MD_MB = MD_MT*(AS_MB/AS_MT)**(1.D0/B5/PI)
*      print*,'at mb^pole:',mt_mb,mb_mb,mc_mb,ms_mb,mu_mb,md_mb
      MT_MC = MT_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MB_MC = MB_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MC_MC = MC_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MS_MC = MS_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MU_MC = MU_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MD_MC = MD_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
*      print*,'at mc^pole:',mt_mc,mb_mc,mc_mc,ms_mc,mu_mc,md_mc
*
*AS(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MT_POLE) THEN
       AS_S = AS_MT/(1.D0+B6*AS_MT*DLOG(SQRTS**2/MT_POLE**2))
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MT_POLE .AND. SQRTS.GT.MB_POLE ) THEN
       AS_S = AS_MZ/(1.D0+B5*AS_MZ*DLOG(SQRTS**2/MZ**2))
*  mc^pole < ss <=mb^pole
      ELSEIF(SQRTS.LE.MB_POLE .AND. SQRTS.GT.MC_POLE  ) THEN
       AS_S = AS_MB/(1.D0+B4*AS_MB*DLOG(SQRTS**2/MB_POLE**2))
*            ss <=mc^pole
      ELSEIF(SQRTS.LE.MC_POLE) THEN
       AS_S = AS_MC/(1.D0+B3*AS_MC*DLOG(SQRTS**2/MC_POLE**2))
      ELSE
       print*,'SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*      print*,'ACKKKKKK',sqrts,as_s
*MQ(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MT_POLE) THEN
       MT_S = MT_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
       MB_S = MB_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
       MC_S = MC_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
       MS_S = MS_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
       MU_S = MU_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
       MD_S = MD_MT*(AS_S/AS_MT)**(1.D0/B6/PI)
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MT_POLE .AND. SQRTS.GT.MB_POLE ) THEN
       MT_S = MT_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
       MB_S = MB_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
       MC_S = MC_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
       MS_S = MS_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
       MU_S = MU_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
       MD_S = MD_MT*(AS_S/AS_MT)**(1.D0/B5/PI)
*  mc^pole < ss <=mb^pole
      ELSEIF(SQRTS.LE.MB_POLE .AND. SQRTS.GT.MC_POLE  ) THEN
       MT_S = MT_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MB_S = MB_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MC_S = MC_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MS_S = MS_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MU_S = MU_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MD_S = MD_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
*            ss <=mc^pole
      ELSEIF(SQRTS.LE.MC_POLE) THEN
       MT_S = MT_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MB_S = MB_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MC_S = MC_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MS_S = MS_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MU_S = MU_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MD_S = MD_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
      ELSE
       print*,'SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*
*      print*,'SQRTS,AS(SQRTS) =',SQRTS,AS_S
*      print*,' > MQ(SQRTS)    :',MT_S,MB_S,MC_S,MS_S,MU_S,MD_S
*-----------------------------------------------------------------------
*
      RETURN
      END

      SUBROUTINE FOUR_FERMION(N_F,N_FP,V,M_F,M_FP,MH,NCMAX,NHC,C4_F_FP)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      COMPLEX*16 C4_F_FP
      COMPLEX*16 NHC(NCMAX,3)
      REAL*8     MH(3)
*
*      print*,n_f,n_fp
*      print*,v,m_f,m_fp,mh(1),mh(2),mh(3)
*
      N_S=N_F +1
      N_P=N_FP+2
      C4_F_FP=M_F*M_FP/V**2*
     .       ( NHC(N_S,1)*NHC(N_P,1)/MH(1)**2
     .        +NHC(N_S,2)*NHC(N_P,2)/MH(2)**2
     .        +NHC(N_S,3)*NHC(N_P,3)/MH(3)**2 )
*-----------------------------------------------------------------------
*
      RETURN
      END

      REAL*8 FUNCTION EDM_A(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF(DABS(X-1.D0).LT.1.D-8) THEN
       EDM_A=-1.D0/3.D0
       RETURN
      ENDIF

      EDM_A=1.D0/2.D0/(1.D0-X)**2*(3.D0-X+2.D0*DLOG(X)/(1.D0-X))
*
      RETURN
      END

      REAL*8 FUNCTION EDM_B(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF(DABS(X-1.D0).LT.1.D-8) THEN
       EDM_B=1.D0/6.D0
       RETURN
      ENDIF

      EDM_B=1.D0/2.D0/(1.D0-X)**2*(1.D0+X+2.D0*X*DLOG(X)/(1.D0-X))
*
      RETURN
      END

      REAL*8 FUNCTION EDM_C(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF(DABS(X-1.D0).LT.1.D-8) THEN
       EDM_C=19.D0/18.D0
       RETURN
      ENDIF

      EDM_C=1.D0/6.D0/(1.D0-X)**2
     . *(10.D0*X-26.D0+2.D0*X*DLOG(X)/(1.D0-X)-18.D0*DLOG(X)/(1.D0-X))
*
      RETURN
      END

      REAL*8 FUNCTION EDM_H(X1)
************************************************************************
*
* EDM function of the purely gluonic dimension six Weinberg operator
* (Neutral-Higgs contribution)
* Ref.: D. A. Dicus, PRD41(1990)999
*
*           1  /1   /1              u^3 x^3 (1-x)
*  h(k) =  --- | dx | du -----------------------------------
*           4  /0   /0    [ x (1 - u x) + k (1-u) (1-x) ]^2
*
*
* The "u" integration has been done by MAPLE:
*
*                1   /1              u^3 x^3 (1-x)
*  EDM_H(x) =   ---  | du -----------------------------------
*                4   /0    [ x (1 - u x) + k (1-u) (1-x) ]^2
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
      complex*16 s1,s2,t0
      complex*16 xlog
      external   xlog
*
      X  = X1
      XK = Z_HiggsEDM
*MAPLE output
*
* [1] fxu:=1/4*u^3*x^3*(1-x)/(x*(1-u*x)+xk*(1-u)*(1-x))^2;
* [2] fx:=int(fxu,u=0..1);
* [3] fortran(fx,mode=double,precision=double);
*
* The argument of LOG can be negative depending on XK. As a prescirption,
* we take "DLOG" -> "xlog": see the COMPLEX FUNCTION XLOG.
*
*
      s1 = -x**2*(-6.D0*xk**3*x**2-6.D0*x**5*xk+27.D0*xk**2*x**2+8.D0*xk
     #*x**4-11.D0*xk**2*x+6.D0*xk**3*x-6.D0*x**3*xlog(x*(-1.D0+x))-4.D0*
     #x**4-10.D0*xk*x**2+24.D0*xlog(x*(-1.D0+x))*xk*x**3-21.D0*x**3*xk**
     #2-6.D0*x*xlog(x*(-1.D0+x))*xk**2-18.D0*xlog(x*(-1.D0+x))*x**3*xk**
     #2-12.D0*x**4*xlog(x*(-1.D0+x))*xk-12.D0*x**2*xlog(x*(-1.D0+x))*xk+
     #6.D0*xlog(x*(-1.D0+x))*x**4*xk**2+18.D0*xlog(x*(-1.D0+x))*x**2*xk*
     #*2+x**6+6.D0*x**4*xlog(x*(-1.D0+x))-2.D0*xk**3+2.D0*xk**3*x**3+3.D
     #0*x**5+8.D0*x**3*xk+5.D0*x**4*xk**2-2.D0*x**3)/(2.D0*xk*x**2-2.D0*
     #xk**2*x-2.D0*x**3*xk+xk**2*x**2+xk**2+x**4)/(x**2+xk-xk*x)**2.D0/8
     #.D0
      s2 = x**3*(x**3*xk**2+3.D0*x**3*xlog(-x-xk+xk*x)*xk**2-6.D0*xlog(-
     #x-xk+xk*x)*xk*x**3-2.D0*x**3*xk+x**3+3.D0*x**3*xlog(-x-xk+xk*x)-9.
     #D0*x**2*xlog(-x-xk+xk*x)*xk**2-3.D0*xk**2*x**2+12.D0*xlog(-x-xk+xk
     #*x)*xk*x**2+4.D0*xk*x**2-3.D0*x**2*xlog(-x-xk+xk*x)-x**2+9.D0*xlog
     #(-x-xk+xk*x)*x*xk**2+3.D0*xk**2*x-6.D0*x*xlog(-x-xk+xk*x)*xk-2.D0*
     #xk*x-3.D0*xlog(-x-xk+xk*x)*xk**2-xk**2)/(x**2+xk-xk*x)**4.D0/4.D0
      t0 = s1+s2
*
*      if(dimag(t0).gt.1.d-6) print*,'ERROR!!! EDM_H',t0
*      if(dimag(t0)/dreal(t0).gt.1.d-4) print*,'WARNING!!! EDM_H',t0
      if(dimag(t0)/dreal(t0).gt.1.d-3) print*,'WARNING!!! EDM_H',t0
      EDM_H=DREAL(t0)
*
      RETURN
      END

      COMPLEX*16 FUNCTION XLOG(X1)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 XI
*
      PI=2.D0*DASIN(1.D0)
      XI=DCMPLX(0.D0,1.D0)
*
      IF(X1.GT.0.D0) XLOG=DLOG(X1)
      IF(X1.LT.0.D0) XLOG=DLOG(-X1)+XI*PI
*
      RETURN
      END

      SUBROUTINE EDM_HH_USER0(STMASS,SBMASS,M3,MBMT,MTMT
     .                       ,EDM_HH_STOP,EDM_HH_SBOT)
************************************************************************
*
*This is the approximated expression of H(z1,z2,zq) valid ONLY when zq
*is small, zq<0.1 or less
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      REAL*8     STMASS(2),SBMASS(2)
*
      REAL*8   EDM_HH_ZS
      EXTERNAL EDM_HH_ZS
*
      COMMON /EDM_HH_BASES/ Z1_BASES,Z2_BASES,ZQ_BASES
*-----------------------------------------------------------------------
      Z1=STMASS(1)**2/M3**2
      Z2=STMASS(2)**2/M3**2
      ZQ=MTMT**2/M3**2
      IF(ZQ.GT.1.0D-1) 
     . print*,'WARNING: H(z1,z2,zt) estimation is incorrect!!!',zq
       Z1_BASES=Z1
       Z2_BASES=Z2
       ZQ_BASES=ZQ
      EDM_HH_STOP=EDM_HH_ZS((Z1+Z2)/2.D0)/ZQ
*
      Z1=SBMASS(1)**2/M3**2
      Z2=SBMASS(2)**2/M3**2
      ZQ=MBMT**2/M3**2
      IF(ZQ.GT.1.0D-1) 
     . print*,'WARNING: H(z1,z2,zb) estimation is incorrect!!!',zq
       Z1_BASES=Z1
       Z2_BASES=Z2
       ZQ_BASES=ZQ
      EDM_HH_SBOT=EDM_HH_ZS((Z1+Z2)/2.D0)/ZQ
*
      RETURN
      END

      REAL*8 FUNCTION EDM_HH_ZS(X1)
************************************************************************
*
*     EDM_HH_ZS= limit [ZQ*H(z1,z2,zq)]
*                ZQ->0
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /EDM_HH_BASES/ Z1_BASES,Z2_BASES,ZQ_BASES
*
      ZS=X1  ! X1=(Z2_BASES+Z1_BASES)/2
      DZ=Z2_BASES-Z1_BASES
*
      IF(DABS(1.D0-ZS).LT.1.D-6) THEN
       H0=5.D0/108.D0
       H1=11.D0/1080.D0
      ELSE
       H0=1.D0/18.D0/(1.D0-ZS)**4
     .   *( 2.D0*(1.D0-ZS)*(1.D0+11.D0*ZS)
     .    -(1.D0-16.D0*ZS-9.D0*ZS**2)*DLOG(ZS) )
       H1=1.D0/108.D0/(1.D0-ZS)**6
     .   *( (1.D0-ZS)*(1.D0+7.D0*ZS+295.D0*ZS**2+177.D0*ZS**3)
     .    +6.D0*ZS**2*(21.D0+50.D0*ZS+9.D0*ZS**2)*DLOG(ZS) )
      ENDIF
*
      EDM_HH_ZS=H0+DZ**2/4.D0/ZS**2*H1
*
      RETURN
      END


      REAL*8 FUNCTION F0_HiggsEDM(X1)
************************************************************************
*
*  F(z)=x*(1-x)/[z-x*(1-x)]*ln[x*(1-x)/z]
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
*
      Z=Z_HiggsEDM
      X=X1
*
      IF(X.EQ.0.D0 .OR. X.EQ.1.D0) THEN
       F0_HiggsEDM=0.D0
      ELSE
       F0_HiggsEDM=X*(1.D0-X)/(Z-X*(1.D0-X))*DLOG(X*(1.D0-X)/Z)
      ENDIF
*
      RETURN
      END

      REAL*8 FUNCTION F_HiggsEDM(X1)
************************************************************************
*
*  f(z)=z/2*[1-2*x*(1-x)]/[x*(1-x)-z]*ln[x*(1-x)/z]
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
*
      Z=Z_HiggsEDM
      X=X1
*
      F_HiggsEDM=Z/2.D0*(1.D0-2.D0*X*(1.D0-X))/(X*(1.D0-X)-Z)
     .       *DLOG(X*(1.D0-X)/Z)
*
      RETURN
      END

      REAL*8 FUNCTION G_HiggsEDM(X1)
************************************************************************
*
*  g(z)=z/2*1/[x*(1-x)-z]*ln[x*(1-x)/z]
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
*
      Z=Z_HiggsEDM
      X=X1
*
      G_HiggsEDM=Z/2.D0/(X*(1.D0-X)-Z)*DLOG(X*(1.D0-X)/Z)
*
      RETURN
      END


      SUBROUTINE MORE_EDMS(
     . NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
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
*-----------------------------------------------------------------------
*Local Vairables:
*
      EXTERNAL EDM_B,EDM_C
* For integration
      COMMON /HiggsEDM_BODE/ Z_HiggsEDM
      EXTERNAL F0_HiggsEDM,F_HiggsEDM,G_HiggsEDM
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
      EXTERNAL EDM_JWH1,EDM_JWH2,EDM_JWH3,EDM_JWH4,EDM_JWW
*
      COMPLEX*16 XI
      REAL*8     FWH1(4,2),FWH2(4,2),FWH3(4,2),FWH4(4,2),FWW(4,2) 
      COMPLEX*16 G_FFP,GS_FFP,GP_FFP
      COMPLEX*16 GS_IJ(4,2),GP_IJ(4,2),GL_IJ(4,2),GR_IJ(4,2)
*
      REAL*8 DEOE_E(2),DEOE_D(2),DEOE_S(2)    
*
      COMPLEX*16 EGD22,EGU22
      COMPLEX*16 H_C,H_S,A_C,A_S
      COMPLEX*16 XRL,XLR,H_Q,A_Q
      REAL*8     SCMASS(2),SSMASS(2)
      COMPLEX*16 SCMIX(2,2),SSMIX(2,2)
*
*Couplings at Mt^pole
      COMPLEX*16 CGL,CGR,H_F,SFMIX(2,2)
*     chargino(I)-strange quark-scharm(J)
      COMPLEX*16 GL_CI_S_SCJ(2,2),GR_CI_S_SCJ(2,2)
*     neutralino(I)-strange quark-sstrange(J)
      COMPLEX*16 GL_NI_S_SSJ(4,2),GR_NI_S_SSJ(4,2)
*     gluino-strange quark-sstrange(I)
      COMPLEX*16 GL_GL_S_SSI(2),  GR_GL_S_SSI(2)
*
      REAL*8     DC_S(20)
*=======================================================================
*      print*,'>>>>> AURUN: START <<<<<'
*=======================================================================
      PI=2.D0*DASIN(1.D0)
      XI=DCMPLX(0.D0,1.D0)
      AEM=1.D0/137.D0
      GEVTOCM=1.97326968D-14
      I_DEOE_E = 200 ! d^E_e/e [cm]
      I_DEOE_D = 220 ! d^E_d/e [cm]
      I_DEOE_S = 230 ! d^E_s/e [cm]
      I_DC_D   = 250 ! d^C_d   [cm]
      I_DC_S   = 400 ! d^C_s   [cm]
*
*Loop functions
      NX_EDM  = 1000   ! Number of calling SUBROUTINE BODE
      EPS_EDM = 1.D-6  ! integration region of 2 loop functions:
*                        [EPS_EDM ; (1-EPS_EDM)]
      X1D=EPS_EDM
      X1U=1.D0-EPS_EDM
*
      DO I=1,4
       DO J=1,2
        RI=MN_H(I)**2/MCH**2
        RJ=MC_H(J)**2/MCH**2
        RW=MW_H**2/MCH**2
         CALL BODE(EDM_JWH1,X1D,X1U,NX_EDM,RES_JWH1)
         CALL BODE(EDM_JWH2,X1D,X1U,NX_EDM,RES_JWH2)
         CALL BODE(EDM_JWH3,X1D,X1U,NX_EDM,RES_JWH3)
         CALL BODE(EDM_JWH4,X1D,X1U,NX_EDM,RES_JWH4)
*         print*,ri,rj,rw,res_jwh1,res_jwh2,res_jwh3,res_jwh4
         FWH1(I,J)=RES_JWH1
         FWH2(I,J)=RES_JWH2
         FWH3(I,J)=RES_JWH3
         FWH4(I,J)=RES_JWH4
        RI=MN_H(I)**2/MW_H**2
        RJ=MC_H(J)**2/MW_H**2
         CALL BODE(EDM_JWW,X1D,X1U,NX_EDM,RES_JWW)
*         print*,ri,rj,res_jww
         FWW(I,J)=RES_JWW
*         print*,i,j,fwh1(i,j),fwh2(i,j),fwh3(i,j),fwh4(i,j),fww(i,j)
       ENDDO
      ENDDO
*
*Charged Higgs-neutralino_I-chargino_J couplings:
      GS_IJ(1,1) =CHC_H(20)
       GP_IJ(1,1)=CHC_H(21)
      GS_IJ(1,2) =CHC_H(23)
       GP_IJ(1,2)=CHC_H(24)
      GS_IJ(2,1) =CHC_H(26)
       GP_IJ(2,1)=CHC_H(27)
      GS_IJ(2,2) =CHC_H(29)
       GP_IJ(2,2)=CHC_H(30)
      GS_IJ(3,1) =CHC_H(32)
       GP_IJ(3,1)=CHC_H(33)
      GS_IJ(3,2) =CHC_H(35)
       GP_IJ(3,2)=CHC_H(36)
      GS_IJ(4,1) =CHC_H(38)
       GP_IJ(4,1)=CHC_H(39)
      GS_IJ(4,2) =CHC_H(41)
       GP_IJ(4,2)=CHC_H(42)
*W boson-neutralino_I-chargino_J couplings:
      DO I=1,4
       DO J=1,2
        GL_IJ(I,J)=N_H(I,3)*DCONJG(UL_H(J,2))
     .            +DSQRT(2.D0)*N_H(I,2)*DCONJG(UL_H(J,1))
        GR_IJ(I,J)=-DCONJG(N_H(I,4))*DCONJG(UR_H(J,2))
     .            +DSQRT(2.D0)*DCONJG(N_H(I,2))*DCONJG(UR_H(J,1))
*        print*,i,j,gs_ij(i,j),gp_ij(i,j)
*        print*,i,j,gl_ij(i,j),gr_ij(i,j)
       ENDDO
      ENDDO
*
*-----------------------------------------------------------------------
*For electron
      G_FFP =CHC_H(1)
      GS_FFP=CHC_H(2)
      GP_FFP=CHC_H(3)
*      print*,g_ffp,gs_ffp,gp_ffp
*
      DEOE_E(1)=0.D0 ! gamma-W-H      
      DO I=1,4
       DO J=1,2
        DEOE_E(1)=DEOE_E(1)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J))) 
     .        )*FWH1(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J))) 
     .        )*FWH2(I,J)*MN_H(I)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J))) 
     .        )*FWH3(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J))) 
     .        )*FWH4(I,J)*MN_H(I)
*        print*,i,j,deoe_e(1)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/64.D0/PI**2/SW_H**4/MCH**2
     .       *(-DSQRT(2.D0)*G_FFP/GW_H)*GEVTOCM
*      print*,fac_ffp
      DEOE_E(1)=FAC_FFP*DEOE_E(1)
*      print*,'(d^E_e/e)^WH [cm] = ',deoe_e(1)
*-----
      DEOE_E(2)=0.D0 ! gamma-W-W
      DO I=1,4
       DO J=1,2
        DEOE_E(2)=DEOE_E(2)
     .  +DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J)) 
     .        )*FWW(I,J)*MN_H(I)*MC_H(J)
*        print*,i,j
*     .  ,DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J)) ),FWW(I,J),MN_H(I)*MC_H(J)
*        print*,i,j,deoe_e(2)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/32.D0/PI**2/SW_H**4/MW_H**4*ME_H*GEVTOCM
*      print*,fac_ffp
      DEOE_E(2)=FAC_FFP*DEOE_E(2)
*      print*,'(d^E_e/e)^WW [cm] = ',deoe_e(2)
*
*STORE:
      RAUX_H(I_DEOE_E+5)=DEOE_E(1)
      RAUX_H(I_DEOE_E+6)=DEOE_E(2)
*      print*,(raux_h(I_DEOE_E+i),i=0,6)
*-----------------------------------------------------------------------
*For down quark
      G_FFP =CHC_H(10)
      GS_FFP=CHC_H(11)
      GP_FFP=CHC_H(12)
*      print*,g_ffp,gs_ffp,gp_ffp
*      print*,gs_ffp+xi*gp_ffp,tb_h*mdmt_h/mumt_h
*
      DEOE_D(1)=0.D0 ! gamma-W-H
      DO I=1,4
       DO J=1,2
        DEOE_D(1)=DEOE_D(1)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J)))
     .        )*FWH1(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J)))
     .        )*FWH2(I,J)*MN_H(I)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J)))
     .        )*FWH3(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J)))
     .        )*FWH4(I,J)*MN_H(I)
*        print*,i,j,deoe_d(1)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/64.D0/PI**2/SW_H**4/MCH**2
     .       *(-DSQRT(2.D0)*G_FFP/GW_H)*GEVTOCM
*      print*,fac_ffp,(-DSQRT(2.D0)*G_FFP/GW_H),mumt_h/mw_h
      DEOE_D(1)=FAC_FFP*DEOE_D(1)
*      print*,'(d^E_d/e)^WH [cm] = ',deoe_d(1)
*-----
      DEOE_D(2)=0.D0 ! gamma-W-W
      DO I=1,4
       DO J=1,2
        DEOE_D(2)=DEOE_D(2)
     .  +DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J))
     .        )*FWW(I,J)*MN_H(I)*MC_H(J)
*        print*,i,j
*     .  ,DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J))
*     ),FWW(I,J),MN_H(I)*MC_H(J)
*        print*,i,j,deoe_d(2)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/32.D0/PI**2/SW_H**4/MW_H**4*MDMT_H*GEVTOCM
*      print*,fac_ffp
      DEOE_D(2)=FAC_FFP*DEOE_D(2)
*      print*,'(d^E_d/e)^WW [cm] = ',deoe_d(2)
*
*STORE:
      RAUX_H(I_DEOE_D+5)=DEOE_D(1)
      RAUX_H(I_DEOE_D+6)=DEOE_D(2)
*      print*,(raux_h(I_DEOE_D+i),i=0,6)
*-----------------------------------------------------------------------
*For strange quark
      G_FFP =CHC_H(13)
      GS_FFP=CHC_H(14)
      GP_FFP=CHC_H(15)
*      print*,g_ffp,gs_ffp,gp_ffp
*      print*,gs_ffp+xi*gp_ffp,tb_h*msmt_h/mcmt_h
*
      DEOE_S(1)=0.D0 ! gamma-W-H
      DO I=1,4
       DO J=1,2
        DEOE_S(1)=DEOE_S(1)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J)))
     .        )*FWH1(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J)))
     .        )*FWH2(I,J)*MN_H(I)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GR_IJ(I,J)-GL_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GR_IJ(I,J)+GL_IJ(I,J)))
     .        )*FWH3(I,J)*MC_H(J)
     .  +DIMAG( (GS_FFP+XI*GP_FFP)
     .         *(   DCONJG(GS_IJ(I,J))*(GL_IJ(I,J)-GR_IJ(I,J))
     .          +XI*DCONJG(GP_IJ(I,J))*(GL_IJ(I,J)+GR_IJ(I,J)))
     .        )*FWH4(I,J)*MN_H(I)
*        print*,i,j,deoe_s(1)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/64.D0/PI**2/SW_H**4/MCH**2
     .       *(-DSQRT(2.D0)*G_FFP/GW_H)*GEVTOCM
*      print*,fac_ffp,(-DSQRT(2.D0)*G_FFP/GW_H),mumt_h/mw_h
      DEOE_S(1)=FAC_FFP*DEOE_S(1)
*      print*,'(d^E_s/e)^WH [cm] = ',deoe_s(1)
*-----
      DEOE_S(2)=0.D0 ! gamma-W-W
      DO I=1,4
       DO J=1,2
        DEOE_S(2)=DEOE_S(2)
     .  +DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J))
     .        )*FWW(I,J)*MN_H(I)*MC_H(J)
*        print*,i,j
*     .  ,DIMAG( GL_IJ(I,J)*DCONJG(GR_IJ(I,J))
*     ),FWW(I,J),MN_H(I)*MC_H(J)
*        print*,i,j,deoe_s(2)
       ENDDO
      ENDDO
      FAC_FFP=AEM_H**2/32.D0/PI**2/SW_H**4/MW_H**4*MSMT_H*GEVTOCM
*      print*,fac_ffp
      DEOE_S(2)=FAC_FFP*DEOE_S(2)
*      print*,'(d^E_s/e)^WW [cm] = ',deoe_s(2)
*
*STORE:
      RAUX_H(I_DEOE_S+5)=DEOE_S(1)
      RAUX_H(I_DEOE_S+6)=DEOE_S(2)
*      print*,(raux_h(I_DEOE_S+i),i=0,6)
*-----------------------------------------------------------------------
*=======================================================================
*
* Calculation of the strage-quark CEDM: DC_S
*
      A_C=SSPARA_H(33)*DCMPLX(DCOS(SSPARA_H(34)/180.D0*PI)
     .                       ,DSIN(SSPARA_H(34)/180.D0*PI))
      A_S=SSPARA_H(37)*DCMPLX(DCOS(SSPARA_H(38)/180.D0*PI)
     .                       ,DSIN(SSPARA_H(38)/180.D0*PI))
*      print*,a_c,a_s
      R_123Q=SSPARA_H(22)
      R_123U=SSPARA_H(23)
      R_123D=SSPARA_H(24)
*
      MQ2=R_123Q*MQ3_H
      MU2=R_123U*MU3_H
      MD2=R_123D*MD3_H
*      print*,mq2,mu2,md2
*
      ASMT=ASMT_H
      GSMT=2.D0*DSQRT(PI*ASMT)
*Stop an Sbottom scales
      BT=(11.D0-2.D0*6.D0/3.D0)/(4.D0*PI)
       QQB=MQ3**2+MBMT_H**2
       QBB=MD3**2+MBMT_H**2
      QB2=DMAX1(QQB,QBB)
       QQT=MQ3**2+MTPOLE_H**2
       QTT=MU3**2+MTPOLE_H**2
      QT2=DMAX1(QQT,QTT)
*      print*,'QB2,QT2=',qb2,qt2
*AS(Stop,Sbottom)
      AS_MSB=ASMT/(1.D0+BT*ASMT*DLOG(QB2/MTPOLE_H**2))
      AS_MST=ASMT/(1.D0+BT*ASMT*DLOG(QT2/MTPOLE_H**2))
*      print*,'As(Mt^pole),As(M_sbottom),AS(M_stop)=',asmt,as_msb,as_mst
*
*The gluino contribution to the threshold corrections
*\Delta_d*sqrt(2)/v_2
*and \Delta_u*sqrt(2)/v_1
      EGD22=2.D0*AS_MSB/3.D0/PI*DCONJG(MU_H*M3_H)
     .          *F_I(MD2**2,MQ2**2,CDABS(M3_H)**2)
      EGU22=2.D0*AS_MST/3.D0/PI*DCONJG(MU_H*M3_H)
     .          *F_I(MU2**2,MQ2**2,CDABS(M3_H)**2)
*      print*,'EGD22: ',egd22
*      print*,'EGU22: ',egu22,egd22/as_msb*as_mst
*
      IF(IFLAG_H(10).EQ.0) THEN
       H_C=DCMPLX(DSQRT(2.D0)*MCMT_H/V_H/SB_H,0.D0)/(1.D0+EGU22/TB_H)
       H_S=DCMPLX(DSQRT(2.D0)*MSMT_H/V_H/CB_H,0.D0)/(1.D0+TB_H*EGD22)
*       print*,h_c,h_s
*       print*,DSQRT(2.D0)*MCMT_H/V_H/SB_H
*       print*,DSQRT(2.D0)*MSMT_H/V_H/CB_H
      ELSE
       H_C=DCMPLX(DSQRT(2.D0)*MCMT_H/V_H/SB_H,0.D0)
       H_S=DCMPLX(DSQRT(2.D0)*MSMT_H/V_H/CB_H,0.D0)
*       print*,h_c,h_s
      ENDIF
*
*Scharm mixing
*
      T_3=1.D0/2.D0
      Q_Q=2.D0/3.D0
      V_Q=V_H*SB_H
      R_Q=CB_H/SB_H
      M_Q=MCMT_H
      H_Q=H_C
      A_Q=A_C

      XLL = MQ2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MU2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SCMASS,SCMIX)
      IF(SCMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative scharm mass squared'
        RETURN
      ENDIF
*
*      print*,'Scharm Sector:'
*      WRITE(*,1) SCMASS(1),SCMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SCMIX(1,1)),DIMAG(SCMIX(1,1))
*     .          ,DREAL(SCMIX(1,2)),DIMAG(SCMIX(1,2))
*      WRITE(*,4) DREAL(SCMIX(2,1)),DIMAG(SCMIX(2,1))
*     .          ,DREAL(SCMIX(2,2)),DIMAG(SCMIX(2,2))
*      print*,' '
*
*Sstrange mixing
*
      T_3=-1.D0/2.D0
      Q_Q=-1.D0/3.D0
      V_Q=V_H*CB_H
      R_Q=SB_H/CB_H
      M_Q=MSMT_H
      H_Q=H_S
      A_Q=A_S

      XLL = MQ2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*(T_3-Q_Q*SW_H**2)
      XRR = MD2**2+DABS(M_Q)**2
     .     +(CB_H**2-SB_H**2)*MZ_H**2*Q_Q*SW_H**2
      XRL = H_Q*V_Q*(A_Q-DCONJG(MU_H)*R_Q)/SQRT(2.D0)
      XLR = DCONJG(XRL)

      CALL SFERMION_MIXING(XLL,XRR,XRL,XLR,SSMASS,SSMIX)
      IF(SSMASS(1).LE.0.D0) THEN
        print*,'ERROR!: Negative sstrange mass squared'
        RETURN
      ENDIF
*
*      print*,'Sstrange Sector:'
*      WRITE(*,1) SSMASS(1),SSMASS(2)
*      print*,'                [1]                     [2]'
*      WRITE(*,3) DREAL(SSMIX(1,1)),DIMAG(SSMIX(1,1))
*     .          ,DREAL(SSMIX(1,2)),DIMAG(SSMIX(1,2))
*      WRITE(*,4) DREAL(SSMIX(2,1)),DIMAG(SSMIX(2,1))
*     .          ,DREAL(SSMIX(2,2)),DIMAG(SSMIX(2,2))
*      print*,' '
*
*-- couplings:
*
*     chargino(I)-strange quark-scharm(J)
*      print*,'chargino(I)-strange quark-scharm(J):'
      DO I=1,2
       DO J=1,2
        GL_CI_S_SCJ(I,J)=-GW_H*UR_H(I,1)*DCONJG(SCMIX(1,J))
     .                    +H_C*UR_H(I,2)*DCONJG(SCMIX(2,J))
        GR_CI_S_SCJ(I,J)=DCONJG(H_S)*UL_H(I,2)*DCONJG(SCMIX(1,J))
*        CGL=GL_CI_S_SCJ(I,J)
*        CGR=GR_CI_S_SCJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     neutralino(I)-strange quark-sstrange(J)
      I_A =3
      T_3F=-1.0/2.D0
      Q_F =-1.D0/3.D0
      H_F =H_S
      DO I=1,2
       DO J=1,2
        SFMIX(I,J)=SSMIX(I,J)
       ENDDO
      ENDDO

*      print*,'neutralino(I)-strange quark-sstrange(J):'
      DO I=1,4
       DO J=1,2
        GL_NI_S_SSJ(I,J)=
     .-DSQRT(2.D0)*GW_H*T_3F*DCONJG(N_H(I,2))*DCONJG(SFMIX(1,J))
     .-DSQRT(2.D0)*GW_H*TW_H*(Q_F-T_3F)*DCONJG(N_H(I,1))
     .                                 *DCONJG(SFMIX(1,J))
     .-H_F*DCONJG(N_H(I,I_A))*DCONJG(SFMIX(2,J))
        GR_NI_S_SSJ(I,J)=
     . DSQRT(2.D0)*GW_H*TW_H*Q_F*N_H(I,1)*DCONJG(SFMIX(2,J))
     .-DCONJG(H_F)*N_H(I,I_A)*DCONJG(SFMIX(1,J))
*        CGL=GL_NI_S_SSJ(I,J)
*        CGR=GR_NI_S_SSJ(I,J)
*        write(*,6) I,J,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
       ENDDO
      ENDDO
*
*     gluino-strange quark-sstrange(I)
      PHI_3=SSPARA_H(10)/180.D0*PI ! in Radian
*      print*,'gluino-strange quark-sstrange(I):'
      DO I=1,2
       GL_GL_S_SSI(I)=-GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0),-DSIN(PHI_3/2.D0))*DCONJG(SSMIX(1,I))
       GR_GL_S_SSI(I)= GSMT/DSQRT(2.D0)
     . *DCMPLX(DCOS(PHI_3/2.D0), DSIN(PHI_3/2.D0))*DCONJG(SSMIX(2,I))
*       CGL=GL_GL_S_SSI(I)
*       CGR=GR_GL_S_SSI(I)
*       write(*,5) I,CGR,CGL,DIMAG(DCONJG(CGR)*CGL)
      ENDDO
************************************************************************
*Convention of the arrays DEOE and DC
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
       NEDM_SUB=20
*Initialize
       DO I=1,NEDM_SUB
        DC_S(I)       =0.D0
       ENDDO
*
*Chargino contributions:
*=======================
      IEDM_SUB=2
      DC_S(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( MC_H(1)/SCMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(1,1))*GL_CI_S_SCJ(1,1))
     .             *EDM_B(MC_H(1)**2/SCMASS(1)**2)
     .  +MC_H(1)/SCMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(1,2))*GL_CI_S_SCJ(1,2))
     .             *EDM_B(MC_H(1)**2/SCMASS(2)**2)
     .  +MC_H(2)/SCMASS(1)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(2,1))*GL_CI_S_SCJ(2,1))
     .             *EDM_B(MC_H(2)**2/SCMASS(1)**2)
     .  +MC_H(2)/SCMASS(2)**2
     .             *DIMAG(DCONJG(GR_CI_S_SCJ(2,2))*GL_CI_S_SCJ(2,2))
     .             *EDM_B(MC_H(2)**2/SCMASS(2)**2) )
*      print*,iedm_sub,DC_S(IEDM_SUB)*gevtocm
*
*Neutralino contributions:
*=========================
      IEDM_SUB=3
      DC_S(IEDM_SUB)= GSMT/16.D0/PI**2
     .*( MN_H(1)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(1,1))*GL_NI_S_SSJ(1,1))
     .             *EDM_B(MN_H(1)**2/SSMASS(1)**2)
     .  +MN_H(1)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(1,2))*GL_NI_S_SSJ(1,2))
     .             *EDM_B(MN_H(1)**2/SSMASS(2)**2)
     .  +MN_H(2)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(2,1))*GL_NI_S_SSJ(2,1))
     .             *EDM_B(MN_H(2)**2/SSMASS(1)**2)
     .  +MN_H(2)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(2,2))*GL_NI_S_SSJ(2,2))
     .             *EDM_B(MN_H(2)**2/SSMASS(2)**2)
     .  +MN_H(3)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(3,1))*GL_NI_S_SSJ(3,1))
     .             *EDM_B(MN_H(3)**2/SSMASS(1)**2)
     .  +MN_H(3)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(3,2))*GL_NI_S_SSJ(3,2))
     .             *EDM_B(MN_H(3)**2/SSMASS(2)**2)
     .  +MN_H(4)/SSMASS(1)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(4,1))*GL_NI_S_SSJ(4,1))
     .             *EDM_B(MN_H(4)**2/SSMASS(1)**2)
     .  +MN_H(4)/SSMASS(2)**2
     .             *DIMAG(DCONJG(GR_NI_S_SSJ(4,2))*GL_NI_S_SSJ(4,2))
     .             *EDM_B(MN_H(4)**2/SSMASS(2)**2) )
*      print*,iedm_sub,DC_S(IEDM_SUB)*gevtocm
*
*Gluino contributions:
*=====================
      IEDM_SUB=4
      DC_S(IEDM_SUB)= -GSMT/8.D0/PI**2
     .*( CDABS(M3_H)/SSMASS(1)**2
     .                  *DIMAG(DCONJG(GR_GL_S_SSI(1))*GL_GL_S_SSI(1))
     .                  *EDM_C(CDABS(M3_H)**2/SSMASS(1)**2)
     .  +CDABS(M3_H)/SSMASS(2)**2
     .                  *DIMAG(DCONJG(GR_GL_S_SSI(2))*GL_GL_S_SSI(2))
     .                  *EDM_C(CDABS(M3_H)**2/SSMASS(2)**2) )
*      print*,DIMAG(DCONJG(GR_GL_S_SSI(1))*GL_GL_S_SSI(1))
*      print*,DIMAG(DCONJG(GR_GL_S_SSI(2))*GL_GL_S_SSI(2))
*      print*,iedm_sub,DC_S(IEDM_SUB)*gevtocm
*
*Two-loop Higgs contributions: the contribution to d^C_s
*=============================
      IEDM_SUB=5
*
      IEDM_SUB_HIGGS=11 ! top
      DO IH=1,3
       Z_HiggsEDM=MTMT_H**2/HMASS_H(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX_EDM,FTOP)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX_EDM,GTOP)
       DC_S(IEDM_SUB_HIGGS)=DC_S(IEDM_SUB_HIGGS)
     .  -(ASMT*GSMT/2.D0)*AEM*MSMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .   *(DREAL(NHC_H(15,IH))*DREAL(NHC_H(26,IH))*FTOP
     .    +DREAL(NHC_H(14,IH))*DREAL(NHC_H(27,IH))*GTOP)
      ENDDO
*      print*,iedm_sub_higgs,DC_S(IEDM_SUB_HIGGS)*gevtocm

      IEDM_SUB_HIGGS=12 ! bottom
      DO IH=1,3
       Z_HiggsEDM=MBMT_H**2/HMASS_H(IH)**2
       CALL BODE(F_HiggsEDM,X1D,X1U,NX_EDM,FBOT)
       CALL BODE(G_HiggsEDM,X1D,X1U,NX_EDM,GBOT)
       DC_S(IEDM_SUB_HIGGS)=DC_S(IEDM_SUB_HIGGS)
     . -(ASMT*GSMT/2.D0)*AEM*MSMT_H/8.D0/PI**2/SW_H**2/MW_H**2
     .  *(DREAL(NHC_H(15,IH))*DREAL(NHC_H(17,IH))*FBOT
     .   +DREAL(NHC_H(14,IH))*DREAL(NHC_H(18,IH))*GBOT)
      ENDDO
*      print*,iedm_sub_higgs,DC_S(IEDM_SUB_HIGGS)*gevtocm

      IEDM_SUB_HIGGS=13 ! stop
      DO IH=1,3
       Z_HiggsEDM=STMASS_H(1)**2/HMASS_H(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX_EDM,FST1)
       Z_HiggsEDM=STMASS_H(2)**2/HMASS_H(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX_EDM,FST2)
       DC_S(IEDM_SUB_HIGGS)=DC_S(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MSMT_H/32.D0/PI**3
     .  *DREAL(NHC_H(15,IH))/HMASS_H(IH)**2
     .  *DREAL(NHC_H(71,IH)*FST1+NHC_H(74,IH)*FST2)
      ENDDO
*      print*,iedm_sub_higgs,DC_S(IEDM_SUB_HIGGS)*gevtocm

      IEDM_SUB_HIGGS=14 ! sbottom
      DO IH=1,3
       Z_HiggsEDM=SBMASS_H(1)**2/HMASS_H(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX_EDM,FSB1)
       Z_HiggsEDM=SBMASS_H(2)**2/HMASS_H(IH)**2
       CALL BODE(F0_HiggsEDM,X1D,X1U,NX_EDM,FSB2)
       DC_S(IEDM_SUB_HIGGS)=DC_S(IEDM_SUB_HIGGS)
     . +(ASMT*GSMT/2.D0)*MSMT_H/32.D0/PI**3
     .  *DREAL(NHC_H(15,IH))/HMASS_H(IH)**2
     .  *DREAL(NHC_H(75,IH)*FSB1+NHC_H(78,IH)*FSB2)
      ENDDO
*      print*,iedm_sub_higgs,DC_S(IEDM_SUB_HIGGS)*gevtocm
*for the sign flips, see the subroutine FILLEDMS
       DC_S(11)=-DC_S(11)
       DC_S(12)=-DC_S(12)
      DC_S(IEDM_SUB)=DC_S(11)+DC_S(12)+DC_S(13)+DC_S(14)
       DC_S(IEDM_SUB)=-DC_S(IEDM_SUB)
*      print*,iedm_sub,DC_S(IEDM_SUB)*gevtocm
*
*TOTAL:
*======
      IEDM_SUB=1
        DC_S(IEDM_SUB)=DC_S(2)+DC_S(3)+DC_S(4)+DC_S(5)
*      print*,iedm_sub,DC_S(IEDM_SUB)*gevtocm
*STORE:Chromo-Electric Down-quark EDM:
      RAUX_H(I_DC_S+0)=DC_S(1)*GEVTOCM
      RAUX_H(I_DC_S+1)=DC_S(2)*GEVTOCM
      RAUX_H(I_DC_S+2)=DC_S(3)*GEVTOCM
      RAUX_H(I_DC_S+3)=DC_S(4)*GEVTOCM
      RAUX_H(I_DC_S+4)=DC_S(5)*GEVTOCM
*      print*,(RAUX_H(I_DC_S+i),i=0,4)
*      print*,(RAUX_H(I_DC_D+i)*msmt_h/mdmt_h,i=0,4)
*=======================================================================
*-------- printout:
*
*      IF(IFLAG_H(18).EQ.2) THEN
*      print*,'---------------------------------------------------------'
*      print*,'The MORE Electric EDMs of particles in cm: e, d, s: '
*      print*,'---------------------------------------------------------'
*      write(*,8) 'd^E_e/e[WH,WW]:'
*     .          ,RAUX_H(I_DEOE_E+5),RAUX_H(I_DEOE_E+6)
*      write(*,8) 'd^E_d/e[WH,WW]:'
*     .          ,RAUX_H(I_DEOE_D+5),RAUX_H(I_DEOE_D+6)
*      write(*,8) 'd^E_s/e[WH,WW]:'
*     .          ,RAUX_H(I_DEOE_S+5),RAUX_H(I_DEOE_S+6)
*      print*,'---------------------------------------------------------'
*      print*,'The Chromo-Electric EDMs of the particles in cm: s: '
*      print*,'---------------------------------------------------------'
*      write(*,8) 'd^C_s  [Total]:',RAUX_H(I_DC_S+0)
*      write(*,8) 'd^C_s  [C,N,Gl,H]:'
*     .          ,RAUX_H(I_DC_S+1),RAUX_H(I_DC_S+2)
*     .          ,RAUX_H(I_DC_S+3),RAUX_H(I_DC_S+4)
*      print*,'---------------------------------------------------------'
*      ENDIF ! IF(IFLAG_H(18).EQ.2) THEN
*=======================================================================
*      print*,'>>>>> AURUN: E N D <<<<<'
*=======================================================================
 99   CONTINUE
*
 1    FORMAT(2X,'Masses',1X,E10.4,1X,E10.4,' GeV')
 3    FORMAT(2X,'[L] /','(',E10.4,1X,E10.4,') '
     .                 ,'(',E10.4,1X,E10.4,')',' \\')
 4    FORMAT(2X,'[R] \\','(',E10.4,1X,E10.4,') '
     .                  ,'(',E10.4,1X,E10.4,')',' /')
 5    FORMAT(4X,I1,' (',E10.4,1X,E10.4,') ','(',E10.4,1X,E10.4,') '
     .      ,E10.4)
 6    FORMAT(2X,I1,1X,I1,' (',E10.4,1X,E10.4,') '
     .      ,'(',E10.4,1X,E10.4,') ',E10.4)
 8    FORMAT(2X,A18,1X,4(1X,E10.4,1X))
*
      RETURN
      END


      REAL*8 FUNCTION EDM_JWH1(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
*
      A=RW
      B=RJ/(1.D0-X)+RI/X
*
      IF(DABS(A).LT.1.D-6) THEN
       XJAB=DLOG(B)/(B-1.D0)
      ELSEIF(DABS(B).LT.1.D-6) THEN
       XJAB=DLOG(A)/(A-1.D0)
      ELSEIF(DABS(A-1.D0).LT.1.D-6 .AND. DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=1.D0/2.D0
      ELSEIF(DABS(A-1.D0).LT.1.D-6) THEN
       XJAB=(B*DLOG(B)-B+1.D0)/(B-1.D0)**2
      ELSEIF(DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=(A*DLOG(A)-A+1.D0)/(A-1.D0)**2
      ELSEIF(DABS(A-B).LT.1.D-6) THEN
       XJAB=(-DLOG(A)+A-1.D0)/(A-1.D0)**2
      ELSE 
       XJAB=(A*DLOG(A)/(A-1.D0)-B*DLOG(B)/(B-1.D0))/(A-B)
      ENDIF
*
      EDM_JWH1=XJAB*X**2/(1.D0-X)
*
      RETURN
      END

      REAL*8 FUNCTION EDM_JWH2(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
*
      A=RW
      B=RJ/(1.D0-X)+RI/X
*
      IF(DABS(A).LT.1.D-6) THEN
       XJAB=DLOG(B)/(B-1.D0)
      ELSEIF(DABS(B).LT.1.D-6) THEN
       XJAB=DLOG(A)/(A-1.D0)
      ELSEIF(DABS(A-1.D0).LT.1.D-6 .AND. DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=1.D0/2.D0
      ELSEIF(DABS(A-1.D0).LT.1.D-6) THEN
       XJAB=(B*DLOG(B)-B+1.D0)/(B-1.D0)**2
      ELSEIF(DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=(A*DLOG(A)-A+1.D0)/(A-1.D0)**2
      ELSEIF(DABS(A-B).LT.1.D-6) THEN
       XJAB=(-DLOG(A)+A-1.D0)/(A-1.D0)**2
      ELSE
       XJAB=(A*DLOG(A)/(A-1.D0)-B*DLOG(B)/(B-1.D0))/(A-B)
      ENDIF
*
      EDM_JWH2=XJAB*(1.D0-X)
*
      RETURN
      END

      REAL*8 FUNCTION EDM_JWH3(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
*
      A=RW
      B=RJ/(1.D0-X)+RI/X
*
      IF(DABS(A).LT.1.D-6) THEN
       XJAB=DLOG(B)/(B-1.D0)
      ELSEIF(DABS(B).LT.1.D-6) THEN
       XJAB=DLOG(A)/(A-1.D0)
      ELSEIF(DABS(A-1.D0).LT.1.D-6 .AND. DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=1.D0/2.D0
      ELSEIF(DABS(A-1.D0).LT.1.D-6) THEN
       XJAB=(B*DLOG(B)-B+1.D0)/(B-1.D0)**2
      ELSEIF(DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=(A*DLOG(A)-A+1.D0)/(A-1.D0)**2
      ELSEIF(DABS(A-B).LT.1.D-6) THEN
       XJAB=(-DLOG(A)+A-1.D0)/(A-1.D0)**2
      ELSE
       XJAB=(A*DLOG(A)/(A-1.D0)-B*DLOG(B)/(B-1.D0))/(A-B)
      ENDIF
*
      EDM_JWH3=XJAB*X/(1.D0-X)
*
      RETURN
      END

      REAL*8 FUNCTION EDM_JWH4(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
*
      A=RW
      B=RJ/(1.D0-X)+RI/X
*
      IF(DABS(A).LT.1.D-6) THEN
       XJAB=DLOG(B)/(B-1.D0)
      ELSEIF(DABS(B).LT.1.D-6) THEN
       XJAB=DLOG(A)/(A-1.D0)
      ELSEIF(DABS(A-1.D0).LT.1.D-6 .AND. DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=1.D0/2.D0
      ELSEIF(DABS(A-1.D0).LT.1.D-6) THEN
       XJAB=(B*DLOG(B)-B+1.D0)/(B-1.D0)**2
      ELSEIF(DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=(A*DLOG(A)-A+1.D0)/(A-1.D0)**2
      ELSEIF(DABS(A-B).LT.1.D-6) THEN
       XJAB=(-DLOG(A)+A-1.D0)/(A-1.D0)**2
      ELSE
       XJAB=(A*DLOG(A)/(A-1.D0)-B*DLOG(B)/(B-1.D0))/(A-B)
      ENDIF
*
      EDM_JWH4=XJAB
*
      RETURN
      END


      REAL*8 FUNCTION EDM_JWW(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMMON /More_HiggsEDM_BODE/ RI,RJ,RW
*
      A=0.D0
      B=RJ/(1.D0-X)+RI/X
*
      IF(DABS(A).LT.1.D-6) THEN
       XJAB=DLOG(B)/(B-1.D0)
      ELSEIF(DABS(B).LT.1.D-6) THEN
       XJAB=DLOG(A)/(A-1.D0)
      ELSEIF(DABS(A-1.D0).LT.1.D-6 .AND. DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=1.D0/2.D0
      ELSEIF(DABS(A-1.D0).LT.1.D-6) THEN
       XJAB=(B*DLOG(B)-B+1.D0)/(B-1.D0)**2
      ELSEIF(DABS(B-1.D0).LT.1.D-6) THEN
       XJAB=(A*DLOG(A)-A+1.D0)/(A-1.D0)**2
      ELSEIF(DABS(A-B).LT.1.D-6) THEN
       XJAB=(-DLOG(A)+A-1.D0)/(A-1.D0)**2
      ELSE
       XJAB=(A*DLOG(A)/(A-1.D0)-B*DLOG(B)/(B-1.D0))/(A-B)
      ENDIF
*
      EDM_JWW=XJAB/(1.D0-X)
*
      RETURN
      END


      SUBROUTINE HGRA_EDMS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
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
*-----------------------------------------------------------------------
*Local Vairables:
*
*=======================================================================
*EDMs of Mercury d^Hg/(e cm)
*
      CALL MORE_HG_EDMS(G_PINN_0,G_PINN_1
     .                 ,D_HG,D_HG_1,D_HG_2,D_HG_3,D_HG_4)
*      print*,'>>> HGRA_EDMS: ',D_HG,D_HG_1,D_HG_2,D_HG_3,D_HG_4
      RAUX_H(410)=D_HG
      RAUX_H(411)=D_HG_1
      RAUX_H(412)=D_HG_2
      RAUX_H(413)=D_HG_3
      RAUX_H(414)=D_HG_4
*----
*[d^Ra(225)/(e cm)]^IV
      GEVTOCM=1.97326968D-14   ! GEVTOCM  = GeV*cm
      FAC_0=-8.7D-2
      FAC_1=3.5D-1
*for d_Ra225/(e cm)
      FAC_0=FAC_0*GEVTOCM
      FAC_1=FAC_1*GEVTOCM
      D_RA_225=(FAC_0*G_PINN_0+FAC_1*G_PINN_1)
*      print*,'>>> HGRA_EDMS: ',d_ra_225
      RAUX_H(420)=D_RA_225
*----
      IF(IFLAG_H(18).EQ.1) THEN
*      print*,'---------------------------------------------------------'
      print*,'     More Mercury and Radium EDMs in units of [e cm]'
      print*,'---------------------------------------------------------'
      write(*,9) 'd^Hg    /(e cm) [Total]=',RAUX_H(410)
      write(*,9) 'd^Hg^I  /(e cm) [Total]=',RAUX_H(411)
      write(*,9) 'd^Hg^II /(e cm) [Total]=',RAUX_H(412)
      write(*,9) 'd^Hg^III/(e cm) [Total]=',RAUX_H(413)
      write(*,9) 'd^Hg^IV /(e cm) [Total]=',RAUX_H(414)
      write(*,9) 'd^Ra    /(e cm) [Total]=',RAUX_H(420)
      print*,'---------------------------------------------------------'
      print*,' '
      ENDIF ! IF(IFLAG(18).EQ.1) THEN
*=======================================================================
 9    FORMAT(2X,A25,1X,E10.4)
*
      RETURN
      END


      SUBROUTINE MORE_HG_EDMS(G_PINN_0,G_PINN_1
     .                       ,D_HG,D_HG_1,D_HG_2,D_HG_3,D_HG_4)
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
*Local Vairables:
*
*=======================================================================
*
      PI=2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
      ASMT=ASMT_H
      GSMT=2.D0*DSQRT(PI*ASMT)
      GEVTOCM=1.97326968D-14   ! GEVTOCM  = GeV*cm
      XKAPPA=0.5D0 ! \pm 0.25
*
*Mercury EDM d^Hg/(e cm)
*
      DEOE_E =RAUX_H(200)  ! d^E_e/e in cm
      DC_U   =RAUX_H(240)  ! d^C_u   in cm
      DC_D   =RAUX_H(250)  ! d^C_d   in cm
      C4DDD  =RAUX_H(287)  ! C_dd/md in cm/GeV^2
      C4SDS  =RAUX_H(288)  ! C_sd/ms in cm/GeV^2
      C4BDB  =RAUX_H(289)  ! C_bd/mb in cm/GeV^2
      CS     =RAUX_H(270)  ! C_S in cm/GeV
      CP     =RAUX_H(271)  ! C_P in cm/GeV
      CPP    =RAUX_H(272)  ! C_P^\prime in cm/GeV
*----
*[d^Hg/(e cm)] in CPsuperH
      D_HG1=1.D-2*(DEOE_E)
      D_HG2=7.D-3*(DC_U-DC_D)/GSMT
      D_HG3=-1.4D-5*( 0.5D0*C4DDD
     .               +3.3D0*XKAPPA*C4SDS
     .               +(1.D0-0.25D0*XKAPPA)*C4BDB )
      D_HG4=3.5D-3*CS
      D_HG5=4.0D-4*(CP-0.2D0*CPP)
      D_HG=D_HG1+D_HG2+D_HG3+D_HG4+D_HG5
*      print*,d_hg1,d_hg2,d_hg3,d_hg4,d_hg5
*      print*,'[d^Hg/(e cm)]^CPsH :',d_hg
*----
*The unit-less pion-N-N couplings:
      G_PINN_1_DC=2.D+14*(DC_U-DC_D)/GSMT                ! unit-less
      G_PINN_1_4F=-8.D-3*( 0.5D0*C4DDD
     .                    +3.3D0*XKAPPA*C4SDS
     .                    +(1.D0-0.25D0*XKAPPA)*C4BDB )  ! in GeV*cm
      G_PINN_1=G_PINN_1_DC+G_PINN_1_4F/GEVTOCM           ! unit-less
*      G_PINN_0=0.2D0*(DC_U+DC_D)/(DC_U-DC_D)*G_PINN_1_DC ! unit-less
      G_PINN_0=0.2D0*(DC_U+DC_D)*2.D14/GSMT ! unit-less
*      print*,g_pinn_0,g_pinn_1_dc
*----
*[d^Hg/(e cm)]^I
      FAC_0=0.D0
      FAC_1=1.8D-3
*for d_Hg/(e cm)
      FAC_0=FAC_0*GEVTOCM  
      FAC_1=FAC_1*GEVTOCM  
*      cor1=7.d-3/(1.8D-3*GEVTOCM*2.d14)  ! corrections to math CPsuperH
*      cor2=1.4d-5/(1.8D-3*8.d-3)
*      print*,cor1,cor2
*      G_PINN_1=cor1*G_PINN_1_DC+cor2*G_PINN_1_4F/GEVTOCM           
      D_HG23=(FAC_0*G_PINN_0+FAC_1*G_PINN_1)
      D_HG_1=D_HG1+D_HG23+D_HG4+D_HG5
*      print*,'[d^Hg/(e cm)]^I    :',d_hg_1
*----
*[d^Hg/(e cm)]^II
      FAC_0=7.6D-6
      FAC_1=1.0D-3
*for d_Hg/(e cm)
      FAC_0=FAC_0*GEVTOCM  
      FAC_1=FAC_1*GEVTOCM  
      D_HG23=(FAC_0*G_PINN_0+FAC_1*G_PINN_1)
      D_HG_2=D_HG1+D_HG23+D_HG4+D_HG5
*      print*,'[d^Hg/(e cm)]^II   :',d_hg_2
*----
*[d^Hg/(e cm)]^III
      FAC_0=1.3D-4
      FAC_1=1.4D-3
*for d_Hg/(e cm)
      FAC_0=FAC_0*GEVTOCM
      FAC_1=FAC_1*GEVTOCM
      D_HG23=(FAC_0*G_PINN_0+FAC_1*G_PINN_1)
      D_HG_3=D_HG1+D_HG23+D_HG4+D_HG5
*      print*,'[d^Hg/(e cm)]^III  :',d_hg_3
*----
*[d^Hg/(e cm)]^IV
      FAC_0=3.1D-4
      FAC_1=9.5D-5
*for d_Hg/(e cm)
      FAC_0=FAC_0*GEVTOCM
      FAC_1=FAC_1*GEVTOCM
      D_HG23=(FAC_0*G_PINN_0+FAC_1*G_PINN_1)
      D_HG_4=D_HG1+D_HG23+D_HG4+D_HG5
*      print*,'[d^Hg/(e cm)]^IV   :',d_hg_4
*----
*=======================================================================
*
      RETURN
      END
