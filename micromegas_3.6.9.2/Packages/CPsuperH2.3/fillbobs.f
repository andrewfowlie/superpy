      SUBROUTINE FILLBOBS(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     .                   ,MH,OH,MCH,M_C,C_L,C_R,STMASS,STMIX)
************************************************************************
*
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
**Input Arrays
      REAL*8     SMPARA(NSMIN),SSPARA(NSSIN)
      INTEGER*8  IFLAG(NFLAG)
      REAL*8     MH(3),OH(3,3)
      REAL*8     M_C(2),STMASS(2)
      COMPLEX*16 C_L(2,2),C_R(2,2),STMIX(2,2)
*-----------------------------------------------------------------------
*Local 
      REAL*8     DELTA(3,3)
      COMPLEX*16 CDELTA(3,3)
      COMPLEX*16 XI,S13C,CTMP,C33(3,3)
      COMPLEX*16 CKM(3,3), CKMD(3,3), CKM_CKMD(3,3) ! CKM matrix
      COMPLEX*16 EGD(3,3), EHD(3,3) ! Diagonal E_G and E_H for Down-type quarks
      COMPLEX*16 EGU(3,3), EHU(3,3) ! Diagonal E_G and E_H for Up-type quarks
      COMPLEX*16 RD(3,3), RDI(3,3), RD_RDI(3,3)    ! Diagonal R and R^-1
      COMPLEX*16 CKMD_RDI(3,3),CKMD_RDI_CKM(3,3)   ! V^dagger R^-1 V
      COMPLEX*16 GL_NH(3,3,3),GR_NH(3,3,3)   ! g^{L,R}_{H[IH] d^bar[I] d[j]}
      COMPLEX*16 GS_NH(3,3,3),GP_NH(3,3,3)   ! g^{S,P}_{H[IH] d^bar[I] d[j]}
      COMPLEX*16 GL_CH(3,3),GR_CH(3,3)       ! g^{L,R}_{CH    d^bar[I] u[j]}
      REAL*8     GS_HLL(3),GP_HLL(3)         ! H-lepton-lepton
      COMPLEX*16 CS,CP,FS,FP,FA              ! B -> ll
      COMPLEX*16 C1SLL,C1SRR,C2LR_DP,C2LR_2HDM,BBAMP ! DM_B
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
      DO I=1,3
       DO J=1,3
        IF(I.EQ.J) THEN
         DELTA(I,J)=1.D0
         CDELTA(I,J)=DCMPLX(1.D0,0.D0)
        ELSE
         DELTA(I,J)=0.D0
         CDELTA(I,J)=DCMPLX(0.D0,0.D0)
        ENDIF
       ENDDO
*       print*,(delta(i,k),k=1,3)
      ENDDO
*-----------------------------------------------------------------------
*===> B and K Parameters:
*-CKM:PDG
      XL     = SMPARA(16)
      A      = SMPARA(17)
      RB     = SMPARA(18)
      EB     = SMPARA(19)
*-Bs
      TAU_BS=1.461D-12/6.58211915D-25 ! in 1/GeV
      M_BS  =5.3696D0                 ! in GeV
*      F_BS  =0.260D0                  ! in GeV, see hep-ph/0606034
      F_BS  =0.230D0                  ! in GeV some smaller value
      RBF_BS=0.294D0                  ! sqrt(B)*f_B in GeV
      ETA_BS=0.55D0
*-Bd
      TAU_BD=1.536D-12/6.58211915D-25 ! in 1/GeV
      M_BD  =5.2794D0                 ! in GeV
*      F_BD  =0.216D0                  ! in GeV, see hep-ph/0606034
      F_BD  =0.200D0                  ! in GeV some smaller value
      RBF_BD=0.244D0                  ! sqrt(B)*f_B in GeV
      ETA_BD=0.55D0
*-Bu
      M_BU=5.279D0  ! in GeV
*-----------------------------------------------------------------------
*===> CKM matrix:PDG
      S12 =XL
      S23 =A*XL**2
      S13C=A*XL**3*(RB+XI*EB)*DSQRT(1.D0-A**2*XL**4)/DSQRT(1.D0-XL**2)
     .    /(1.D0-A**2*XL**4*(RB+XI*EB))
      C12=DSQRT(1.D0-S12**2)
      C23=DSQRT(1.D0-S23**2)
      C13=DSQRT(1.D0-CDABS(S13C)**2)

      CKM(1,1)= C12*C13
      CKM(1,2)= S12*C13
      CKM(1,3)= DCONJG(S13C)
      CKM(2,1)=-S12*C23-C12*S23*S13C
      CKM(2,2)= C12*C23-S12*S23*S13C
      CKM(2,3)= S23*C13
      CKM(3,1)= S12*S23-C12*C23*S13C
      CKM(3,2)=-C12*S23-S12*C23*S13C
      CKM(3,3)= C23*C13
*
*NO CKM MIXING
*      DO I=1,3
*       DO J=1,3
*        CKM(I,J)=CDELTA(I,J)
*       ENDDO
*      ENDDO
*
      DO I=1,3
       DO J=1,3
        CKMD(I,J)=DCONJG(CKM(J,I))
       ENDDO
      ENDDO
*
      DO I=1,3
       DO J=1,3
        CKM_CKMD(I,J)=CKM(I,1)*CKMD(1,J)
     .               +CKM(I,2)*CKMD(2,J)
     .               +CKM(I,3)*CKMD(3,J)
       ENDDO
      ENDDO
*
*      print*,'V_CKM :'
*      write(*,11) (dreal(ckm(1,j)),dimag(ckm(1,j)),j=1,3)
*      write(*,12) (dreal(ckm(2,j)),dimag(ckm(2,j)),j=1,3)
*      write(*,13) (dreal(ckm(3,j)),dimag(ckm(3,j)),j=1,3)
*      print*,'V_CKM^dagger :'
*      write(*,11) (dreal(ckmd(1,j)),dimag(ckmd(1,j)),j=1,3)
*      write(*,12) (dreal(ckmd(2,j)),dimag(ckmd(2,j)),j=1,3)
*      write(*,13) (dreal(ckmd(3,j)),dimag(ckmd(3,j)),j=1,3)
*      print*,'V_CKM * V_CKM^dagger :'
*      write(*,11) (dreal(ckm_ckmd(1,j)),dimag(ckm_ckmd(1,j)),j=1,3)
*      write(*,12) (dreal(ckm_ckmd(2,j)),dimag(ckm_ckmd(2,j)),j=1,3)
*      write(*,13) (dreal(ckm_ckmd(3,j)),dimag(ckm_ckmd(3,j)),j=1,3)
*-----------------------------------------------------------------------
      R_123Q=SSPARA(22)
      R_123U=SSPARA(23)
      R_123D=SSPARA(24)
*===> EGD and EGU
      MQ1=R_123Q*MQ3_H
      MQ2=R_123Q*MQ3_H
      MQ3=       MQ3_H
      MU1=R_123U*MU3_H
      MU2=R_123U*MU3_H
      MU3=       MU3_H
      MD1=R_123D*MD3_H
      MD2=R_123D*MD3_H
      MD3=       MD3_H
*      print*,'M123:',mq1,mq2,mq3
*      print*,'M123:',mu1,mu2,mu3
*      print*,'M123:',md1,md2,md3
*
      ASMT=ASMT_H
*Stop an Sbottom scales
      BT=(11.D0-2.D0*6.D0/3.D0)/(4.D0*PI)
       QQB=MQ3**2+MBMT_H**2
       QBB=MD3**2+MBMT_H**2
      QB2=DMAX1(QQB,QBB)
       QQT=MQ3**2+MTPOLE_H**2
       QTT=MU3**2+MTPOLE_H**2
      QT2=DMAX1(QQT,QTT)
*      print*,qb2,qt2
*AS(Stop,Sbottom)
      AS_MSB=ASMT/(1.D0+BT*ASMT*DLOG(QB2/MTPOLE_H**2))
      AS_MST=ASMT/(1.D0+BT*ASMT*DLOG(QT2/MTPOLE_H**2))
*      print*,'As(Mt^pole),As(M_sbottom),AS(M_stop)=',asmt,as_msb,as_mst
*YT(M_Stop),YB(M_Sbot)
      YT_MT=DSQRT(2.D0)*MTMT_H/V_H/SB_H
      YB_MT=DSQRT(2.D0)*MBMT_H/V_H/CB_H
      BHT=1.D0/(16.D0*PI**2)
     .   *(9.D0*YT_MT**2/2.D0+YB_MT**2/2.D0-32.D0*PI*ASMT)
      BHB=1.D0/(16.D0*PI**2)
     .   *(9.D0*YB_MT**2/2.D0+YT_MT**2/2.D0-32.D0*PI*ASMT)
*      print*,bht,bhb
      YT_MST=YT_MT*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
      YB_MSB=YB_MT*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
*      print*,'YT(MT^pole):YT(M_Stop) =',yt_mt,yt_mst
*      print*,'YB(MT^pole):YB(M_Sbot) =',yb_mt,yb_msb
*EGD and EGU
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
*--
*      print*,'EGD(1,1): ',egd(1,1)
*      print*,'EGD(2,2): ',egd(2,2)
*      print*,'EGD(3,3): ',egd(3,3)
*      print*,'EGU(1,1): ',egu(1,1),egd(1,1)/as_msb*as_mst
*      print*,'EGU(2,2): ',egu(2,2),egd(2,2)/as_msb*as_mst
*      print*,'EGU(3,3): ',egu(3,3),egd(3,3)/as_msb*as_mst
*
*EHD assuming A_u=A_c=A_t
*EHU assuming A_d=A_s=A_b
      YT=YT_MST ! Yukawa couplings at stop mass scale
      YC=DSQRT(2.D0)*MCMT_H/V_H/SB_H
      YU=DSQRT(2.D0)*MUMT_H/V_H/SB_H 
*      print*,yu**2,yc**2,yt**2
      YB=YB_MSB ! Yukawa couplings at sbottom mass scale
      YS=DSQRT(2.D0)*MSMT_H/V_H/CB_H
      YD=DSQRT(2.D0)*MDMT_H/V_H/CB_H 
*      print*,yb**2,ys**2,yd**2
      DO I=1,3
       DO J=1,3
        EHD(I,J)=DCMPLX(0.D0,0.D0)
        EHU(I,J)=DCMPLX(0.D0,0.D0)
       ENDDO
      ENDDO
      EHD(1,1)=1.D0/16.D0/PI**2*DCONJG(MU_H*AT_H)*YU**2
     .         *F_I(MU1**2,MQ1**2,CDABS(MU_H)**2)
      EHD(2,2)=1.D0/16.D0/PI**2*DCONJG(MU_H*AT_H)*YC**2
     .         *F_I(MU2**2,MQ2**2,CDABS(MU_H)**2)
      EHD(3,3)=1.D0/16.D0/PI**2*DCONJG(MU_H*AT_H)*YT**2
     .         *F_I(MU3**2,MQ3**2,CDABS(MU_H)**2)
      EHU(1,1)=1.D0/16.D0/PI**2*DCONJG(MU_H*AB_H)*YD**2
     .         *F_I(MD1**2,MQ1**2,CDABS(MU_H)**2)
      EHU(2,2)=1.D0/16.D0/PI**2*DCONJG(MU_H*AB_H)*YS**2
     .         *F_I(MD2**2,MQ2**2,CDABS(MU_H)**2)
      EHU(3,3)=1.D0/16.D0/PI**2*DCONJG(MU_H*AB_H)*YB**2
     .         *F_I(MD3**2,MQ3**2,CDABS(MU_H)**2)
*      print*,'EHD(1,1)[before 2HDM]: ',ehd(1,1)
*      print*,'EHD(2,2)[before 2HDM]: ',ehd(2,2)
*      print*,'EHD(3,3)[before 2HDM]: ',ehd(3,3)
*      print*,'EHU(1,1)[before 2HDM]: ',ehu(1,1)
*      print*,'EHU(2,2)[before 2HDM]: ',ehu(2,2)
*      print*,'EHU(3,3)[before 2HDM]: ',ehu(3,3)
*--
*2HDM Contribution: See 'B observables' paper by Ellis, Lee, Pilaftsis
*
*      print*,'MCH          ',raux_h(10)
*      print*,'Q_t^2        ',raux_h(11),qt2
*      print*,'Q_b^2        ',raux_h(12),qb2
*      print*,'Q_tb^2       ',raux_h(13)
*      print*,'v1(mt^pol)   ',raux_h(14)
*      print*,'v2(mt^pol)   ',raux_h(18)
*      print*,'|ht0(mt^pol)|',raux_h(22),dsqrt(2.d0)*mtmt_h/v_h/sb_h
*      print*,'|hb0(mt^pol)|',raux_h(23),dsqrt(2.d0)*mbmt_h/v_h/cb_h
*      print*,'|ht(mt^pol)| ',raux_h(24)
*      print*,'|hb(mt^pol)| ',raux_h(27)
*      print*,'tanb(mt^pol) ',raux_h(18)/raux_h(14)
*      print*,'MCH          ',raux_h(10)
*      print*,'MA^2         ',raux_h(30)
*      print*,'RePi_{H^+H^-}',raux_h(31)
*      print*,'L4 v^2/2     ',raux_h(32)
*      print*,'L4           ',raux_h(33)
*      print*,'L1           ',raux_h(34)
*      print*,'L2           ',raux_h(35)
*      print*,'L34=L3+L4    ',raux_h(36)
      REM12SQ=CB_H*SB_H*(RAUX_H(10)**2-RAUX_H(32)+RAUX_H(31))
*      print*,'ReM12^2      ',rem12sq
      X1=3.D0/8.D0/PI**2*(  RAUX_H(27)**2*MBMT_H**2
     .                    *(DLOG(MBMT_H**2/MTPOLE_H**2)-1.D0) )
      MU1_BARSQ=-(TB_H*REM12SQ+RAUX_H(34)*RAUX_H(14)**2
     .                        +RAUX_H(36)*RAUX_H(18)**2/2.D0+X1)
*      print*,'Mu1_bar^2    ',MU1_BARSQ,X1
      X2=3.D0/8.D0/PI**2*(  RAUX_H(24)**2*MTMT_H**2
     .                    *(DLOG(MTMT_H**2/MTPOLE_H**2)-1.D0) )
      MU2_BARSQ=-(REM12SQ/TB_H+RAUX_H(35)*RAUX_H(18)**2
     .                        +RAUX_H(36)*RAUX_H(14)**2/2.D0+X2)
*      print*,'Mu2_bar^2    ',MU2_BARSQ,X2
      MU1_1SQ=-3.D0/16.D0/PI**2
     .*(RAUX_H(24)**2*CDABS(MU_H)**2*DLOG(RAUX_H(11)/MTPOLE_H**2)
     . +RAUX_H(27)**2*CDABS(AB_H)**2*DLOG(RAUX_H(12)/MTPOLE_H**2))
      MHDSQ=-MU1_BARSQ-CDABS(MU_H)**2+MU1_1SQ
*      print*,'M_Hd^2       ',MHDSQ,cdabs(mu_h)**2,mu1_1sq
      MU2_1SQ=-3.D0/16.D0/PI**2
     .*(RAUX_H(24)**2*CDABS(AT_H)**2*DLOG(RAUX_H(11)/MTPOLE_H**2)
     . +RAUX_H(27)**2*CDABS(MU_H)**2*DLOG(RAUX_H(12)/MTPOLE_H**2))
      MHUSQ=-MU2_BARSQ-CDABS(MU_H)**2+MU2_1SQ
*      print*,'M_Hu^2       ',MHUSQ,cdabs(mu_h)**2,mu2_1sq
      FAC_2HDM=1.D0/16.D0/PI**2*REM12SQ/(MHDSQ-MHUSQ)
     . *DLOG(DABS( (MHDSQ+CDABS(MU_H)**2)/(MHUSQ+CDABS(MU_H)**2) ))
*      FAC_2HDM=0.D0
*      print*,'>>>>>>>>>>  FAC_2HDM     ',fac_2hdm
*
      EHD(1,1)=EHD(1,1)+FAC_2HDM*YU**2
      EHD(2,2)=EHD(2,2)+FAC_2HDM*YC**2
      EHD(3,3)=EHD(3,3)+FAC_2HDM*YT**2
      EHU(1,1)=EHU(1,1)+FAC_2HDM*YD**2
      EHU(2,2)=EHU(2,2)+FAC_2HDM*YS**2
      EHU(3,3)=EHU(3,3)+FAC_2HDM*YB**2
*      print*,'EHD(1,1)[after  2HDM]: ',ehd(1,1)
*      print*,'EHD(2,2)[after  2HDM]: ',ehd(2,2)
*      print*,'EHD(3,3)[after  2HDM]: ',ehd(3,3)
*      print*,'EHU(1,1)[after  2HDM]: ',ehu(1,1)
*      print*,'EHU(2,2)[after  2HDM]: ',ehu(2,2)
*      print*,'EHU(3,3)[after  2HDM]: ',ehu(3,3)
*      print*,'EHD(3,3)/YT**2: ',ehd(3,3)/yt**2
*-----------------------------------------------------------------------
*      print*,'EGD(1,1)+EHD(1,1)[after 2HDM]: ',egd(1,1)+ehd(1,1)
*      print*,'EGD(2,2)+EHD(2,2)[after 2HDM]: ',egd(2,2)+ehd(2,2)
*      print*,'EGD(3,3)+EHD(3,3)[after 2HDM]: ',egd(3,3)+ehd(3,3)
*===>R and R^-1 
      DO I=1,3
       DO J=1,3
        RD(I,J)=DELTA(I,J)+TB_H*(EGD(I,J)+EHD(I,J))
*        RD(I,J)=DELTA(I,J)+TB_H*(EGD(I,J))
*        RD(I,J)=DELTA(I,J)+TB_H*(EHD(I,J))
       ENDDO
      ENDDO
*
      CTMP=RD(1,1)*RD(2,2)*RD(3,3)-RD(1,1)*RD(2,3)*RD(3,2)
     .    -RD(2,1)*RD(1,2)*RD(3,3)+RD(2,1)*RD(1,3)*RD(3,2)
     .    +RD(3,1)*RD(1,2)*RD(2,3)-RD(3,1)*RD(1,3)*RD(2,2)
      IF(CDABS(CTMP).EQ.0.D0) THEN
       PRINT*,'ERROR to get R^-1'
       STOP
      ENDIF
*
      RDI(1,1)= (RD(2,2)*RD(3,3)-RD(2,3)*RD(3,2))/CTMP
      RDI(1,2)=-(RD(1,2)*RD(3,3)-RD(1,3)*RD(3,2))/CTMP
      RDI(1,3)= (RD(1,2)*RD(2,3)-RD(1,3)*RD(2,2))/CTMP
      RDI(2,1)=-(RD(2,1)*RD(3,3)-RD(2,3)*RD(3,1))/CTMP
      RDI(2,2)= (RD(1,1)*RD(3,3)-RD(1,3)*RD(3,1))/CTMP
      RDI(2,3)=-(RD(1,1)*RD(2,3)-RD(1,3)*RD(2,1))/CTMP
      RDI(3,1)= (RD(2,1)*RD(3,2)-RD(2,2)*RD(3,1))/CTMP
      RDI(3,2)=-(RD(1,1)*RD(3,2)-RD(1,2)*RD(3,1))/CTMP
      RDI(3,3)= (RD(1,1)*RD(2,2)-RD(1,2)*RD(2,1))/CTMP
*
      DO I=1,3
       DO J=1,3
        RD_RDI(I,J)=RD(I,1)*RDI(1,J)+RD(I,2)*RDI(2,J)+RD(I,3)*RDI(3,J)
       ENDDO
      ENDDO
*      print*,'R :'
*      write(*,11) (dreal(rd(1,j)),dimag(rd(1,j)),j=1,3)
*      write(*,12) (dreal(rd(2,j)),dimag(rd(2,j)),j=1,3)
*      write(*,13) (dreal(rd(3,j)),dimag(rd(3,j)),j=1,3)
*      print*,'R^-1 :'
*      write(*,11) (dreal(rdi(1,j)),dimag(rdi(1,j)),j=1,3)
*      write(*,12) (dreal(rdi(2,j)),dimag(rdi(2,j)),j=1,3)
*      write(*,13) (dreal(rdi(3,j)),dimag(rdi(3,j)),j=1,3)
*      print*,'R*R^-1 :'
*      write(*,11) (dreal(rd_rdi(1,j)),dimag(rd_rdi(1,j)),j=1,3)
*      write(*,12) (dreal(rd_rdi(2,j)),dimag(rd_rdi(2,j)),j=1,3)
*      write(*,13) (dreal(rd_rdi(3,j)),dimag(rd_rdi(3,j)),j=1,3)
*Checking R^-1*(EG+EU)=(1-R^-1)/tanb
      DO I=1,3
       DO J=1,3
        C33(I,J)=RDI(I,1)*(EGD(1,J)+EHD(1,J))
     .          +RDI(I,2)*(EGD(2,J)+EHD(2,J))
     .          +RDI(I,3)*(EGD(3,J)+EHD(3,J))
     .          -(DELTA(I,J)-RDI(I,J))/TB_H
       ENDDO
      ENDDO
*      print*,'R^-1*(EG+EU)-(1-R^-1)/tanb = 0:'
*      write(*,11) (dreal(c33(1,j)),dimag(c33(1,j)),j=1,3)
*      write(*,12) (dreal(c33(2,j)),dimag(c33(2,j)),j=1,3)
*      write(*,13) (dreal(c33(3,j)),dimag(c33(3,j)),j=1,3)
*-----------------------------------------------------------------------
*===>V^dagger R^-1 V:
      DO I=1,3
       DO J=1,3
        CKMD_RDI(I,J)=CKMD(I,1)*RDI(1,J)
     .               +CKMD(I,2)*RDI(2,J)
     .               +CKMD(I,3)*RDI(3,J)
       ENDDO
      ENDDO
      DO I=1,3
       DO J=1,3
        CKMD_RDI_CKM(I,J)=CKMD_RDI(I,1)*CKM(1,J)
     .                   +CKMD_RDI(I,2)*CKM(2,J)
     .                   +CKMD_RDI(I,3)*CKM(3,J)
       ENDDO
      ENDDO
*
*      print*,'V^dagger*R^-1*V :'
*      write(*,11) (dreal(ckmd_rdi_ckm(1,j))
*     .            ,dimag(ckmd_rdi_ckm(1,j)),j=1,3)
*      write(*,12) (dreal(ckmd_rdi_ckm(2,j))
*     .            ,dimag(ckmd_rdi_ckm(2,j)),j=1,3)
*      write(*,13) (dreal(ckmd_rdi_ckm(3,j))
*     .            ,dimag(ckmd_rdi_ckm(3,j)),j=1,3)
**
*      do i=1,3
*       do j=1,3
*        c33(i,j)=ckmd(i,1)*rdi(1,1)*ckm(1,j)
*     .          +ckmd(i,2)*rdi(2,2)*ckm(2,j)
*     .          +ckmd(i,3)*rdi(3,3)*ckm(3,j)
*       enddo
*      enddo
*      write(*,11) (dreal(c33(1,j)),dimag(c33(1,j)),j=1,3)
*      write(*,12) (dreal(c33(2,j)),dimag(c33(2,j)),j=1,3)
*      write(*,13) (dreal(c33(3,j)),dimag(c33(3,j)),j=1,3)
*-----------------------------------------------------------------------
*===>GL_NH(H[IH],d^bar[I],d[J]) and GR_NH(H_IH,d^bar[I],d[J])
*      print*,mh(1),mh(2),mh(3)
*      print*,(oh(1,k),k=1,3)
*      print*,(oh(2,k),k=1,3)
*      print*,(oh(3,k),k=1,3)

      DO IH=1,3
       DO I=1,3
        DO J=1,3
         GL_NH(IH,I,J)=            CKMD_RDI_CKM(I,J)*OH(1,IH)/CB_H
     .               +(DELTA(I,J)-CKMD_RDI_CKM(I,J))*OH(2,IH)/SB_H
     .    -XI*(DELTA(I,J)-CKMD_RDI_CKM(I,J)/CB_H**2)*OH(3,IH)/TB_H
        ENDDO
       ENDDO
      ENDDO
*
      DO IH=1,3
       DO I=1,3
        DO J=1,3
         GR_NH(IH,I,J)=DCONJG(GL_NH(IH,J,I))
        ENDDO
       ENDDO
      ENDDO
*
      DO IH=1,3
       DO I=1,3
        DO J=1,3
         GS_NH(IH,I,J)=   (GL_NH(IH,I,J)+GR_NH(IH,I,J))/2.D0
         GP_NH(IH,I,J)=XI*(GL_NH(IH,I,J)-GR_NH(IH,I,J))/2.D0
        ENDDO
       ENDDO
      ENDDO
*
*      do ih=1,3
*      print*,'---------------------------------------------------------'
*      print*,'GS for IH =',ih,' :  O[phi_1, IH]/CB = ',oh(1,ih)/cb_h
*      write(*,11) (dreal(gs_nh(ih,1,j)),dimag(gs_nh(ih,1,j)),j=1,3)
*      write(*,12) (dreal(gs_nh(ih,2,j)),dimag(gs_nh(ih,2,j)),j=1,3)
*      write(*,13) (dreal(gs_nh(ih,3,j)),dimag(gs_nh(ih,3,j)),j=1,3)
*      print*,'GP for IH =',ih,' : -O[a, IH]*TB = ',-oh(3,ih)*tb_h
*      write(*,11) (dreal(gp_nh(ih,1,j)),dimag(gp_nh(ih,1,j)),j=1,3)
*      write(*,12) (dreal(gp_nh(ih,2,j)),dimag(gp_nh(ih,2,j)),j=1,3)
*      write(*,13) (dreal(gp_nh(ih,3,j)),dimag(gp_nh(ih,3,j)),j=1,3)
*      enddo
*      print*,'---------------------------------------------------------'
*-----------------------------------------------------------------------
*GL_CH(3,3),GR_CH(3,3)       ! g^{L,R}_{CH    d^bar[I] u[J]}
      DO I=1,3
       DO J=1,3
        GL_CH(I,J)=1.D0/TB_H*CKMD(I,J)*(1.D0-RDI(J,J)/CB_H**2)
        GR_CH(I,J)=-1.D0/TB_H*CKMD(I,J)
     .             +CKMD(I,J)*DCONJG(EGU(J,J)+EHU(J,J))
       ENDDO
      ENDDO
*
*      print*,'GL for CH [Tree]:'
*      write(*,11) (dreal(-tb_h*ckmd(1,j))
*     .            ,dimag(-tb_h*ckmd(1,j)),j=1,3)
*      write(*,12) (dreal(-tb_h*ckmd(2,j))
*     .            ,dimag(-tb_h*ckmd(2,j)),j=1,3)
*      write(*,13) (dreal(-tb_h*ckmd(3,j))
*     .            ,dimag(-tb_h*ckmd(3,j)),j=1,3)
*      print*,'GL for CH :'
*      write(*,11) (dreal(gl_ch(1,j)),dimag(gl_ch(1,j)),j=1,3)
*      write(*,12) (dreal(gl_ch(2,j)),dimag(gl_ch(2,j)),j=1,3)
*      write(*,13) (dreal(gl_ch(3,j)),dimag(gl_ch(3,j)),j=1,3)
*      print*,'GR for CH :'
*      write(*,11) (dreal(gr_ch(1,j)),dimag(gr_ch(1,j)),j=1,3)
*      write(*,12) (dreal(gr_ch(2,j)),dimag(gr_ch(2,j)),j=1,3)
*      write(*,13) (dreal(gr_ch(3,j)),dimag(gr_ch(3,j)),j=1,3)
*-----------------------------------------------------------------------
*===>Br(Bs->mu mu):
      BS_RES=1.D0
*
*CS, CP, and C10 
      DO IH=1,3
       GS_HLL(IH)= OH(1,IH)/CB_H
       GP_HLL(IH)=-OH(3,IH)*TB_H
      ENDDO 
*
      AEM_0=AEM_H
      CS=2.D0*PI*MMU_H/AEM_0/CKM(3,3)/DCONJG(CKM(3,2))
     .  *(GR_NH(1,2,3)*GS_HLL(1)/(MH(1)**2-BS_RES*M_BS**2)
     .   +GR_NH(2,2,3)*GS_HLL(2)/(MH(2)**2-BS_RES*M_BS**2)
     .   +GR_NH(3,2,3)*GS_HLL(3)/(MH(3)**2-BS_RES*M_BS**2))
*      print*,'C_S[Bs->mu mu]',CS,GR_NH(1,2,3),GR_NH(2,2,3),GR_NH(3,2,3)
*JSL 05/FEB/2009
      CAUX_H(169)=CS
*
      CP=XI*2.D0*PI*MMU_H/AEM_0/CKM(3,3)/DCONJG(CKM(3,2))
     .  *(GR_NH(1,2,3)*GP_HLL(1)/(MH(1)**2-BS_RES*M_BS**2)
     .   +GR_NH(2,2,3)*GP_HLL(2)/(MH(2)**2-BS_RES*M_BS**2)
     .   +GR_NH(3,2,3)*GP_HLL(3)/(MH(3)**2-BS_RES*M_BS**2))
*      print*,'C_P[Bs->mu mu]',CP
*JSL 05/FEB/2009
      CAUX_H(170)=CP
*
      C10=-4.221D0
*JSL 05/FEB/2009
      CAUX_H(171)=C10
*
*FS, FP, and FA
      FS=-XI/2.D0*M_BS**2*F_BS*CS*MBMT_H/(MBMT_H+MSMT_H)
      FP=-XI/2.D0*M_BS**2*F_BS*CP*MBMT_H/(MBMT_H+MSMT_H)
      FA=-XI/2.D0*F_BS*C10

*For the SM prediction
*      FS=DCMPLX(0.D0,0.D0)
*      FP=DCMPLX(0.D0,0.D0)
*
      BSMM=GF_H**2*AEM_0**2/16.D0/PI**3*M_BS*TAU_BS
     .    *CDABS(CKM(3,3)*DCONJG(CKM(3,2)))**2
     .    *DSQRT(1.D0-4.D0*MMU_H**2/M_BS**2)
     .    *( (1.D0-4.D0*MMU_H**2/M_BS**2)*CDABS(FS)**2
     .      +CDABS(FP+2.D0*MMU_H*FA)**2 )
      BSMM_SM=GF_H**2*AEM_0**2/16.D0/PI**3*M_BS*TAU_BS
     .    *CDABS(CKM(3,3)*DCONJG(CKM(3,2)))**2
     .    *DSQRT(1.D0-4.D0*MMU_H**2/M_BS**2)
     .    *CDABS(2.D0*MMU_H*FA)**2 
*      print*,'> BSMM*10^7, BSMM[SM]*10^7 = ',bsmm*1.D7,bsmm_sm*1.D7
      RAUX_H(130)=BSMM*1.0D7
*-----------------------------------------------------------------------
*===>Br(Bd->tau tau):
      BD_RES=1.D0
*
*CS, CP, and C10
      DO IH=1,3
       GS_HLL(IH)= OH(1,IH)/CB_H
       GP_HLL(IH)=-OH(3,IH)*TB_H
      ENDDO
*
      AEM_0=AEM_H
      CS=2.D0*PI*MTAU_H/AEM_0/CKM(3,3)/DCONJG(CKM(3,1))
     .  *(GR_NH(1,1,3)*GS_HLL(1)/(MH(1)**2-BD_RES*M_BD**2)
     .   +GR_NH(2,1,3)*GS_HLL(2)/(MH(2)**2-BD_RES*M_BD**2)
     .   +GR_NH(3,1,3)*GS_HLL(3)/(MH(3)**2-BD_RES*M_BD**2))
*      print*,'C_S[Bd->tau tau]',CS
*
      CP=XI*2.D0*PI*MTAU_H/AEM_0/CKM(3,3)/DCONJG(CKM(3,1))
     .  *(GR_NH(1,1,3)*GP_HLL(1)/(MH(1)**2-BD_RES*M_BD**2)
     .   +GR_NH(2,1,3)*GP_HLL(2)/(MH(2)**2-BD_RES*M_BD**2)
     .   +GR_NH(3,1,3)*GP_HLL(3)/(MH(3)**2-BD_RES*M_BD**2))
*      print*,'C_P[Bd->tau tau]',CP
*
      C10=-4.221D0
*
*FS, FP, and FA
      FS=-XI/2.D0*M_BD**2*F_BD*CS*MBMT_H/(MBMT_H+MDMT_H)
      FP=-XI/2.D0*M_BD**2*F_BD*CP*MBMT_H/(MBMT_H+MDMT_H)
      FA=-XI/2.D0*F_BD*C10
*      print*,cdabs(fs),cdabs(fp)

*For the SM prediction
*      FS=DCMPLX(0.D0,0.D0)
*      FP=DCMPLX(0.D0,0.D0)
*
      BDTT=GF_H**2*AEM_0**2/16.D0/PI**3*M_BD*TAU_BD
     .    *CDABS(CKM(3,3)*DCONJG(CKM(3,1)))**2
     .    *DSQRT(1.D0-4.D0*MTAU_H**2/M_BD**2)
     .    *( (1.D0-4.D0*MTAU_H**2/M_BD**2)*CDABS(FS)**2
     .      +CDABS(FP+2.D0*MTAU_H*FA)**2 )
      BDTT_SM=GF_H**2*AEM_0**2/16.D0/PI**3*M_BD*TAU_BD
     .    *CDABS(CKM(3,3)*DCONJG(CKM(3,1)))**2
     .    *DSQRT(1.D0-4.D0*MTAU_H**2/M_BD**2)
     .    *CDABS(2.D0*MTAU_H*FA)**2
*      print*,'> BDTT*10^7, BDTT[SM]*10^7 = ',bdtt*1.D7,bdtt_sm*1.D7
      RAUX_H(131)=BDTT*1.0D7
*-----------------------------------------------------------------------
*\Delta M_Bd
      IQ     =1
      MQMT   =MDMT_H
      RBF    =RBF_BD/0.230D0
      AMP_FAC=1711.D0 ! in 1/ps
      ETA_B  =ETA_BD/0.55D0
*
      C1SLL=-16.D0*PI**2*MBMT_H**2/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GL_NH(1,3,IQ)*GL_NH(1,3,IQ)/MH(1)**2
     .   +GL_NH(2,3,IQ)*GL_NH(2,3,IQ)/MH(2)**2
     .   +GL_NH(3,3,IQ)*GL_NH(3,3,IQ)/MH(3)**2)
*      print*,16.D0*PI**2,MBMT_H**2,DSQRT(2.D0),GF_H,MW_H**2
*      print*,GL_NH(1,3,IQ)*GL_NH(1,3,IQ),MH(1)**2
*      print*,GL_NH(2,3,IQ)*GL_NH(2,3,IQ)
*     .      +GL_NH(3,3,IQ)*GL_NH(3,3,IQ),MH(3)**2
      C1SRR=-16.D0*PI**2*MQMT**2/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GR_NH(1,3,IQ)*GR_NH(1,3,IQ)/MH(1)**2
     .   +GR_NH(2,3,IQ)*GR_NH(2,3,IQ)/MH(2)**2
     .   +GR_NH(3,3,IQ)*GR_NH(3,3,IQ)/MH(3)**2)
      C2LR_DP=-32.D0*PI**2*MBMT_H*MQMT/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GL_NH(1,3,IQ)*GR_NH(1,3,IQ)/MH(1)**2
     .   +GL_NH(2,3,IQ)*GR_NH(2,3,IQ)/MH(2)**2
     .   +GL_NH(3,3,IQ)*GR_NH(3,3,IQ)/MH(3)**2)
      C2LR_2HDM=-2.D0*MBMT_H*MQMT/MW_H**2
     .  *(DCONJG(CKM(3,3))*CKM(3,IQ))**2*TB_H**2
*      print*,'M_BD',cdabs(c1sll),cdabs(c1srr)
*      print*,'M_BD',cdabs(c2lr_dp),cdabs(c2lr_2hdm)
      BBAMP=AMP_FAC*RBF**2*ETA_B
     .      *(0.88D0*(C2LR_DP+C2LR_2HDM)-0.52D0*(C1SLL+C1SRR))
*
      DMBD_SUSY=2.D0*CDABS(BBAMP)
*      print*,'M_BD',dmbd_susy
      CAUX_H(150)=BBAMP
      RAUX_H(132)=DMBD_SUSY
*-----------------------------------------------------------------------
*\Delta M_Bs
      IQ     =2
      MQMT   =MSMT_H
      RBF    =RBF_BS/0.265D0
      AMP_FAC=2310.D0 ! in 1/ps
      ETA_B  =ETA_BS/0.55D0
*
      C1SLL=-16.D0*PI**2*MBMT_H**2/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GL_NH(1,3,IQ)*GL_NH(1,3,IQ)/MH(1)**2
     .   +GL_NH(2,3,IQ)*GL_NH(2,3,IQ)/MH(2)**2
     .   +GL_NH(3,3,IQ)*GL_NH(3,3,IQ)/MH(3)**2)
      C1SRR=-16.D0*PI**2*MQMT**2/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GR_NH(1,3,IQ)*GR_NH(1,3,IQ)/MH(1)**2
     .   +GR_NH(2,3,IQ)*GR_NH(2,3,IQ)/MH(2)**2
     .   +GR_NH(3,3,IQ)*GR_NH(3,3,IQ)/MH(3)**2)
      C2LR_DP=-32.D0*PI**2*MBMT_H*MQMT/DSQRT(2.D0)/GF_H/MW_H**2
     .  *(GL_NH(1,3,IQ)*GR_NH(1,3,IQ)/MH(1)**2
     .   +GL_NH(2,3,IQ)*GR_NH(2,3,IQ)/MH(2)**2
     .   +GL_NH(3,3,IQ)*GR_NH(3,3,IQ)/MH(3)**2)
      C2LR_2HDM=-2.D0*MBMT_H*MQMT/MW_H**2
     .  *(DCONJG(CKM(3,3))*CKM(3,IQ))**2*TB_H**2
*      print*,'M_BS',cdabs(c1sll),cdabs(c1srr)
*      print*,'M_BS',cdabs(c2lr_dp),cdabs(c2lr_2hdm)
      BBAMP=AMP_FAC*RBF**2*ETA_B
     .      *(0.88D0*(C2LR_DP+C2LR_2HDM)-0.52D0*(C1SLL+C1SRR))
*
      DMBS_SUSY=2.D0*CDABS(BBAMP)
*      print*,'M_BS',dmbs_susy
      CAUX_H(151)=BBAMP
      RAUX_H(133)=DMBS_SUSY
*-----------------------------------------------------------------------
* B -> tau nu
      MCH=RAUX_H(10)
*       print*,'B -> tau nu',m_bu,mch
*       print*,'B -> tau nu',tb_msusy,ckm(1,3),rdi(1,1)
*       print*,'B -> tau nu',gl_ch(3,1),ckm(1,3)
*       print*,'B -> tau nu',dconjg(gl_ch(3,1))/ckm(1,3)
      RBTAUNU=CDABS( 1.D0
     .       +(M_BU/MCH)**2*TB_H*DCONJG(GL_CH(3,1))/CKM(1,3) )**2
*      print*,'RBtau',RBTAUNU
*     .      ,M3_H/CDABS(M3_H),TB_H,DCONJG(GL_CH(3,1))/CKM(1,3)
      RAUX_H(134)=RBTAUNU
*-----------------------------------------------------------------------
*b -> s gamma
      CALL BTOSGAM(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     .            ,MCH,M_C,C_L,C_R,STMASS,STMIX
     .            ,RDI,GL_NH,GR_NH,GS_NH,GP_NH,GL_CH,GR_CH)
*-----------------------------------------------------------------------
      IF(IFLAG(16).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'                     B Observables '
      print*,'---------------------------------------------------------'
*      write(*,21) bsmm*1.d7,bsmm_sm*1.d7
*      write(*,22) bdtt*1.d7,bdtt_sm*1.d7
      write(*,21) bsmm*1.d7
      write(*,26) raux_h(135)
      write(*,25) rbtaunu
      write(*,22) bdtt*1.d7
      write(*,27) raux_h(136)
      write(*,23) dmbd_susy
      write(*,24) dmbs_susy
      print*,'---------------------------------------------------------'
      print*,' '
      ENDIF
*-----------------------------------------------------------------------
  11  format(1x,'/',3(1x,'(',e10.4,1x,e10.4,')',1x),'\\')
  12  format(1x,'|',3(1x,'(',e10.4,1x,e10.4,')',1x),'|')
  13  format(1x,'\\',3(1x,'(',e10.4,1x,e10.4,')',1x),'/')
  14  format(2x,3(1x,e10.4,1x))
*  21  format(1x,'B(B_s->mu  mu ) x 10^7 = ',e10.4,1x,'[SM = ',e10.4,']')
*  22  format(1x,'B(B_d->tau tau) x 10^7 = ',e10.4,1x,'[SM = ',e10.4,']')
  21  format(2x,'B(B_s -> mu  mu )   x 10^7 = ',e10.4)
  22  format(2x,'B(B_d -> tau tau)   x 10^7 = ',e10.4)
  23  format(2x,'Delta M [B_d] (SUSY)       = ',e10.4,1x,'[1/ps]')
  24  format(2x,'Delta M [B_s] (SUSY)       = ',e10.4,1x,'[1/ps]')
  25  format(2x,'B(B_u -> tau nu)/B(SM)     = ',e10.4,1x)
  26  format(2x,'B(B   -> X_s gamma) x 10^4 = ',e10.4)
  27  format(2x,'ACP(B -> X_s gamma) x 10^2 = ',e10.4,1x,'[%]')
*
      RETURN
      END


      SUBROUTINE BTOSGAM(NSMIN,NSSIN,SMPARA,SSPARA,NFLAG,IFLAG
     .                  ,MCH,M_C,C_L,C_R,STMASS,STMIX
     .                  ,RDI,GL_NH,GR_NH,GS_NH,GP_NH,GL_CH,GR_CH)
************************************************************************
*
* Following conventions of microOMEGAs:V1.3 manual (CPC 174(2006)577)
* and references there in
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
**Input Arrays
      REAL*8     SMPARA(NSMIN),SSPARA(NSSIN)
      INTEGER*8  IFLAG(NFLAG)
      REAL*8     M_C(2),STMASS(2)
      COMPLEX*16 C_L(2,2),C_R(2,2),STMIX(2,2)
      COMPLEX*16 RDI(3,3)                    ! Diagonal R_d^-1
      COMPLEX*16 GL_NH(3,3,3),GR_NH(3,3,3)   ! g^{L,R}_{H[IH] d^bar[I] d[j]}
      COMPLEX*16 GS_NH(3,3,3),GP_NH(3,3,3)   ! g^{S,P}_{H[IH] d^bar[I] d[j]}
      COMPLEX*16 GL_CH(3,3),GR_CH(3,3)       ! g^{L,R}_{CH    d^bar[I] u[j]}
*-----------------------------------------------------------------------
      COMMON /FIJ_BODE/ DGAM_FIJ,Z_FIJ
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*-----------------------------------------------------------------------
      COMPLEX*16 B2SG_G7_H,B2SG_G8_H
      COMPLEX*16 B2SG_D7_H,B2SG_D8_H
* 
      EXTERNAL XLI2
      EXTERNAL B2SG_F22,B2SG_F27
      EXTERNAL B2SG_F7_1,B2SG_F8_1
      EXTERNAL B2SG_F7_2,B2SG_F8_2
      EXTERNAL B2SG_F7_3,B2SG_F8_3
      EXTERNAL B2SG_G7_1,B2SG_G8_1
      EXTERNAL B2SG_G7_H,B2SG_G8_H
      EXTERNAL B2SG_D7_H,B2SG_D8_H
      EXTERNAL B2SG_E
      EXTERNAL B2SG_EH
      EXTERNAL B2SG_GG
*-----------------------------------------------------------------------
*Local 
      COMPLEX*16 XI
      COMPLEX*16 S13C,CKM(3,3),CKMDAG(3,3)
      COMPLEX*16 C7_0_MW,C8_0_MW,C7_1_MW,C8_1_MW
      COMPLEX*16 C7_0_MW_CH,C8_0_MW_CH,C7_1_MW_CH,C8_1_MW_CH
      COMPLEX*16 C7_0_MW_CINO,C8_0_MW_CINO
      COMPLEX*16 C7_0_MS_CINO,C8_0_MS_CINO
      COMPLEX*16 C2_0,C7_0,C8_0,C7_1,C7_EM
      COMPLEX*16 AU_C,AD_C
      REAL*8     K22,K27,K28,K77,K78,K88
      REAL*8     K_NLO,K_NLO_0,K_NLO_1,K_NLO_EM
      REAL*8     H_MAGIC(8),A_MAGIC(8),H8_MAGIC(4),B_MAGIC(4)
      REAL*8     E_MAGIC(8),F_MAGIC(8),G_MAGIC(8)
      REAL*8     XJI(2,2)
      COMPLEX*16 CF_1,CF_2
      COMPLEX*16 CKMDAG_RDI(3,3)
*-----------------------------------------------------------------------
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*-CKM
      XL     = SMPARA(16)
      A      = SMPARA(17)
      RB     = SMPARA(18)
      EB     = SMPARA(19)
*
      S12 =XL
      S23 =A*XL**2
      S13C=A*XL**3*(RB+XI*EB)*DSQRT(1.D0-A**2*XL**4)/DSQRT(1.D0-XL**2)
     .    /(1.D0-A**2*XL**4*(RB+XI*EB))
      C12=DSQRT(1.D0-S12**2)
      C23=DSQRT(1.D0-S23**2)
      C13=DSQRT(1.D0-CDABS(S13C)**2)
      CKM(1,1)= C12*C13
      CKM(1,2)= S12*C13
      CKM(1,3)= DCONJG(S13C)
      CKM(2,1)=-S12*C23-C12*S23*S13C
      CKM(2,2)= C12*C23-S12*S23*S13C
      CKM(2,3)= S23*C13
      CKM(3,1)= S12*S23-C12*C23*S13C
      CKM(3,2)=-C12*S23-S12*C23*S13C
      CKM(3,3)= C23*C13
      DO I=1,3
       DO J=1,3
        CKMDAG(I,J)=DCONJG(CKM(J,I))
       ENDDO
      ENDDO
*      print*,'V_CKM :'
*      write(*,11) (dreal(ckm(1,j)),dimag(ckm(1,j)),j=1,3)
*      write(*,12) (dreal(ckm(2,j)),dimag(ckm(2,j)),j=1,3)
*      write(*,13) (dreal(ckm(3,j)),dimag(ckm(3,j)),j=1,3)
*      print*,'|V_CKM| :'
*      write(*,21) (cdabs(ckm(1,j)),j=1,3)
*      write(*,22) (cdabs(ckm(2,j)),j=1,3)
*      write(*,23) (cdabs(ckm(3,j)),j=1,3)
*-some constants
      MB_POLE=RAUX_H(1)
      MB_MB  =RAUX_H(2)
      AS_MB  =RAUX_H(3)
      MC_POLE=RAUX_H(4)
      MC_MC  =RAUX_H(5)
      AS_MC  =RAUX_H(6)
      AS_MT = ASMT_H      
      B5    = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      MC_MB = MCMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
*
      Z0     =(MC_POLE/MB_POLE)**2
      FZ0    =1.D0-8.D0*Z0+8.D0*Z0**3-Z0**4-12.D0*Z0**2*DLOG(Z0)
      DEL_GAM=1.0D0/3.D0 ! Photon-energy cut-off = 1.6 GeV for mb=4.8
*
*For the common bloak /FIJ_BODE/
      DGAM_FIJ=DEL_GAM
*To capture part of NNLO correction?:
*The smaller MU_C => The larger running charm-quark mass MC_R => The smaller BR
      MU_C=MC_POLE
      MC_R=MC_MC
*      MU_C=MB_POLE
*      MC_R=MC_MB
      Z_FIJ   =(MC_R/MB_POLE)**2 ! rather than Z0, see Ref.[41]
*-----------------------------------------------------------------------
      RMU     = MB_POLE
      B4      = (11.D0-2.D0/3.D0*4.D0)/4.D0/PI
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      IF(RMU.GT.MB_POLE) THEN
       AS_RMU = AS_MB/(1.D0+B5*AS_MB*DLOG(RMU**2/MB_POLE**2))
      ELSEIF(RMU.LT.MB_POLE) THEN
       AS_RMU = AS_MB/(1.D0+B4*AS_MB*DLOG(RMU**2/MB_POLE**2))
      ELSEIF(RMU.EQ.MB_POLE) THEN
       AS_RMU = AS_MB
      ELSE
       print*,'ERROR: AS_RMU'
      ENDIF
*-----------------------------------------------------------------------
*Test of Dilogarithmic function: XLI2
*      print*,'Test of Dilogarithmic function: XLI2'
*      print*,' 1',xli2(2.0d0),' same? ',pi**2/4.d0
*      print*,' 2',xli2(1.d0),' same? ',pi**2/6.d0
*      print*,' 3',xli2(0.5d0),' same? ',pi**2/12.d0-(dlog(2.d0))**2/2.d0
*      print*,' 4',xli2(0.0d0),' same? ',0.d0
*      print*,' 5',xli2(-1.d0),' same? ',-pi**2/12.d0
*      print*,' 7',xli2(3.0d0)+xli2(-3.d0),' same? ',xli2(9.d0)/2.d0
*      print*,' 8',xli2(1.5d0)+xli2(-1.5d0),' same? ',xli2(2.25d0)/2.d0
*      print*,' 9',xli2(0.7d0)+xli2(-0.7d0),' same? ',xli2(0.49d0)/2.d0
*      print*,'10',xli2(0.3d0)+xli2(-0.3d0),' same? ',xli2(0.09d0)/2.d0
*      print*,'11',xli2(0.5d0)+xli2(-0.5d0),' same? ',xli2(0.25d0)/2.d0
*      print*,'12',xli2(2.0d0)+xli2(-2.0d0),' same? ',xli2(4.00d0)/2.d0
*      print*,'13',xli2(9.0d0)+xli2(-9.0d0),' same? ',xli2(81.00d0)/2.d0
*-----------------------------------------------------------------------
*f_ij:
*
      NSTEP=500
      CALL BODE(B2SG_F22,1.D-9,1.D0,NSTEP,F22)
      CALL BODE(B2SG_F27,1.D-9,1.D0,NSTEP,F27)
      F28=-F27/3.D0
      F77=(10.D0*DEL_GAM+DEL_GAM**2-2.D0*DEL_GAM**3/3.D0
     .    +DEL_GAM*(DEL_GAM-4.D0)*DLOG(DEL_GAM))/3.D0
      F78=8.D0/9.D0*(XLI2(1.D0-DEL_GAM)-PI**2/6.D0-DEL_GAM*DLOG(DEL_GAM)
     .              +9.D0*DEL_GAM/4.D0-DEL_GAM**2/4.D0+DEL_GAM**3/12.D0)
      F88=(4.D0*XLI2(1.D0-DEL_GAM)-2.D0*PI**2/3.D0
     .    +8.D0*DLOG(1.D0-DEL_GAM)
     .    -DEL_GAM*(2.D0+DEL_GAM)*DLOG(DEL_GAM)+7.D0*DEL_GAM
     .    +3.D0*DEL_GAM**2-2.D0*DEL_GAM**3/3.D0
     .    -2.D0*(2.D0*DEL_GAM+DEL_GAM**2
     .          +4.D0*DLOG(1.D0-DEL_GAM))*DLOG(50.D0) ! mb/ms = 50
     .     )/27.D0
*      print*,f22,0.107636D0-0.208484D0*dsqrt(z_fij)-0.156146D0*z_fij
*      print*,f27,-0.190805D0+0.948865D0*dsqrt(z_fij)-0.787805D0*z_fij
*      print*,f77,f78,f88
*      print*,'fij:',f22,f27,f28,f77,f78,f88
*-----------------------------------------------------------------------
*k_ij
      SUDAKOV=DEXP(-2.D0*AS_RMU/3.D0/PI
     .             *((DLOG(DEL_GAM))**2+7.D0/2.D0*DLOG(DEL_GAM)) )
*      print*,sudakov
      K22=AS_RMU/PI*F22
      K28=AS_RMU/PI*F28
      K88=AS_RMU/PI*F88
*      print*,k22,k28,k88
      G77 =32.D0/3.D0
      G27 =416.D0/81.D0
      G87 =-32.D0/9.D0
      XL2 =0.12D0
      R7  =-10.D0/3.D0-8.D0*PI**2/9.D0
      RER8=44.D0/9.D0-8.D0*PI**2/27.D0
      RER2=-4.987D0+12.78D0*(DSQRT(Z_FIJ)-0.22D0)
      BKZ =3.672D0-4.14D0*(DSQRT(Z_FIJ)-0.22D0)
*
      K77=SUDAKOV*(1.D0
     .   +AS_RMU/2.D0/PI*(R7+G77*DLOG(MB_POLE/RMU)-16.D0/3.D0)
     .   +((1.D0-Z0)**4/FZ0-1.D0)*6.D0*XL2/MB_POLE**2)
     .   +AS_RMU/PI*F77+SUDAKOV*AS_RMU/2.D0/PI*BKZ ! Taking bar(mu)_b = mu_b
*      print*,AS_RMU/2.D0/PI*(R7+G77*DLOG(MB_POLE/RMU)-16.D0/3.D0)
*      print*,((1.D0-Z0)**4/FZ0-1.D0)*6.D0*XL2/MB_POLE**2
*      print*,AS_RMU/PI*F77,SUDAKOV*AS_RMU/2.D0/PI*BKZ 
*      print*,k77
      K27=SUDAKOV*(AS_RMU/2.D0/PI*(RER2+G27*DLOG(MB_POLE/RMU))
     .            -XL2/9.D0/MB_POLE**2/Z0)+AS_RMU/PI*F27
*      print*,AS_RMU/2.D0/PI*(RER2+G27*DLOG(MB_POLE/RMU))
*      print*,-XL2/9.D0/MB_POLE**2/Z0,AS_RMU/PI*F27
*      print*,'K_27[b -> s gam]:',dsqrt(z_fij),k27
      K78=SUDAKOV*AS_RMU/2.D0/PI*(RER8+G87*DLOG(MB_POLE/RMU))
     .   +AS_RMU/PI*F78
*      print*,SUDAKOV*AS_RMU/2.D0/PI*(RER8+G87*DLOG(MB_POLE/RMU))
*      print*,AS_RMU/PI*F78
*      print*,k78
*      print*,'kij:',k22,k27,k28,k77,k78,k88
*-----------------------------------------------------------------------
*MAGIC numbers
      H_MAGIC(1) = 626126.D0/272277.D0
      H_MAGIC(2) =-56281.D0/51730.D0
      H_MAGIC(3) =-3.D0/7.D0
      H_MAGIC(4) =-1.D0/14.D0
      H_MAGIC(5) =-0.6494D0
      H_MAGIC(6) =-0.0380D0
      H_MAGIC(7) =-0.0186D0
      H_MAGIC(8) =-0.0057D0
*
      A_MAGIC(1) = 14.0/23.D0
      A_MAGIC(2) = 16.0/23.D0
      A_MAGIC(3) = 6.0/23.D0
      A_MAGIC(4) =-12.D0/23.D0
      A_MAGIC(5) = 0.4086D0
      A_MAGIC(6) =-0.4230D0
      A_MAGIC(7) =-0.8994D0
      A_MAGIC(8) = 0.1456D0
*
      H8_MAGIC(1)=-0.9135D0
      H8_MAGIC(2)= 0.0873D0
      H8_MAGIC(3)=-0.0571D0
      H8_MAGIC(4)= 0.0209D0
*
      B_MAGIC(1) = 0.4086D0
      B_MAGIC(2) =-0.4230D0
      B_MAGIC(3) =-0.8994D0
      B_MAGIC(4) = 0.1456D0
*
      E_MAGIC(1) = 4661194.D0/816831.D0
      E_MAGIC(2) =-8516.D0/2217.D0
      E_MAGIC(3) = 0.D0
      E_MAGIC(4) = 0.D0
      E_MAGIC(5) =-1.9043D0
      E_MAGIC(6) =-0.1008D0
      E_MAGIC(7) = 0.1216D0
      E_MAGIC(8) = 0.0183D0
*
      F_MAGIC(1) =-17.3023D0
      F_MAGIC(2) = 8.5027D0
      F_MAGIC(3) = 4.5508D0
      F_MAGIC(4) = 0.7519D0
      F_MAGIC(5) = 2.0040D0
      F_MAGIC(6) = 0.7476D0
      F_MAGIC(7) =-0.5385D0
      F_MAGIC(8) = 0.0914D0
*
      G_MAGIC(1) = 14.8088D0
      G_MAGIC(2) =-10.8090D0
      G_MAGIC(3) =-0.8740D0
      G_MAGIC(4) = 0.4218D0
      G_MAGIC(5) =-2.9347D0
      G_MAGIC(6) = 0.3971D0
      G_MAGIC(7) = 0.1600D0
      G_MAGIC(8) = 0.0225D0
*-----------------------------------------------------------------------
*The SM contribution: C7_[0,1]_MW_SM and C8_[0,1]_MW_SM
      AS_MT = ASMT_H      
      AS_MZ = ASMZ_H
      AS_MW = AS_MB/(1.D0+B5*AS_MB*DLOG(MW_H**2/MB_POLE**2))
      AS_160= AS_MB/(1.D0+B5*AS_MB*DLOG(160.D0**2/MB_POLE**2))
*      print*,as_mt,as_mz,as_mw
*      VS_MW =1.D0-23.D0/3.D0*AS_MZ/2.D0/PI*DLOG(MZ_H/MW_H)
*      AS_MW = AS_MZ/VS_MW
*     .       *(1.D0-116.D0/23.D0*AS_MZ/4.D0/PI*DLOG(VS_MW)/VS_MW)
*      print*,as_mw
*
      MT_MT=MTMT_H
      MT_MW=MT_MT*(AS_MW/AS_MT)**(12.D0/23.D0)*( 1.D0
     . +AS_MT/4.D0/PI*(8.D0/2.D0/(23.D0/3.D0))
     .               *(1012.D0/9.D0/8.D0-116.D0/3.D0/(23.D0/3.D0))
     .               *(AS_MW/AS_MT-1.D0) )
*      print*,mt_mt,mt_mw
      XTW=(MT_MW/MW_H)**2
*      print*,xtw
      C7_0_MW_SM=B2SG_F7_1(XTW)
      C8_0_MW_SM=B2SG_F8_1(XTW)
*      print*,XTW,C7_0_MW_SM,C8_0_MW_SM
      C7_1_MW_SM=B2SG_G7_1(XTW)
      C8_1_MW_SM=B2SG_G8_1(XTW)
*      print*,XTW,C7_1_MW_SM,C8_1_MW_SM
*      print*,'C7(MW)[SM]^0',C7_0_MW_SM
*      print*,'C8(MW)[SM]^0',C8_0_MW_SM
*      print*,'  C7(MW)[SM]^1',C7_1_MW_SM
*      print*,'  C8(MW)[SM]^1',C8_1_MW_SM
*-----------------------------------------------------------------------
*Charged Higgs 
*
      AU_C=DCONJG(GR_CH(3,3))/CKM(3,3)*GR_CH(2,3)/CKMDAG(2,3)
      AD_C=DCONJG(GL_CH(3,3))/CKM(3,3)*GR_CH(2,3)/CKMDAG(2,3)
*For COMMON /C78_1_CH/ AU_CH,AD_CH
      AU_CH=AU_C
      AD_CH=AD_C
*      print*,'>>> MAIN',tb_h,au_ch,ad_ch
*
      XTH=(MT_MW/MCH)**2
*      print*,mch,xth
      C7_0_MW_CH=AU_C*B2SG_F7_1(XTH)/3.D0
     .          +AD_C*B2SG_F7_2(XTH)
      C8_0_MW_CH=AU_C*B2SG_F8_1(XTH)/3.D0
     .          +AD_C*B2SG_F8_2(XTH)
*      print*,B2SG_F7_2(XTH),B2SG_F8_2(XTH)
      C7_1_MW_CH=B2SG_G7_H(XTH)+B2SG_D7_H(XTH)*DLOG(MW_H**2/MCH**2)
     .          -4.D0/9.D0*B2SG_EH(XTH)
*      print*,DLOG(MW_H**2/MCH**2)
*      print*,'==========> G7_H ',B2SG_G7_H(XTH)
*      print*,'==========> D7_H ',B2SG_D7_H(XTH)
*      print*,'==========> E_H  ',B2SG_EH(XTH)
      C8_1_MW_CH=B2SG_G8_H(XTH)+B2SG_D8_H(XTH)*DLOG(MW_H**2/MCH**2)
     .          -1.D0/6.D0*B2SG_EH(XTH)
*      print*,'==========> G8_H ',B2SG_G8_H(XTH)
*      print*,'==========> D8_H ',B2SG_D8_H(XTH)
*
*      print*,'C7(MW)[CH]^0',C7_0_MW_CH
*      print*,'C8(MW)[CH]^0',C8_0_MW_CH
*      print*,'  C7(MW)[CH]^1',C7_1_MW_CH
*      print*,'  C8(MW)[CH]^1',C8_1_MW_CH
*-----------------------------------------------------------------------
*Charginos
*      print*,m_c(1),m_c(2)
*      print*,c_l(1,1),c_l(1,2)
*      print*,c_l(2,1),c_l(2,2)
*      print*,c_r(1,1),c_r(1,2)
*      print*,c_r(2,1),c_r(2,2)
*      print*,stmass(1),stmass(2)
*      print*,stmix(1,1),stmix(1,2)
*      print*,stmix(2,1),stmix(2,2)
      MST=DSQRT(STMASS(1)*STMASS(2))   ! stop-mass scale
      R_123Q=SSPARA(22)
      MSQ=R_123Q*MQ3_H                 ! common sup_L- and scharm_L-mass scale
*      print*,stmass(1),stmass(2),mst,msq
*
      MT_POLE=MTPOLE_H
      B6     = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      AS_MST = AS_MT/(1.D0+B6*AS_MT*DLOG(MST**2/MT_POLE**2))
      MT_MST = MT_MT*(AS_MST/AS_MT)**(1.D0/B6/PI)
*      print*,as_mz,as_mt,as_mst,mt_mt,mt_mst
*      print*,1.d0/b6/pi,12.d0/21.d0
*
      XQ1=(MSQ/M_C(1))**2
      XQ2=(MSQ/M_C(2))**2
      DO J=1,2
       DO I=1,2
        XJI(J,I)=(STMASS(J)/M_C(I))**2
       ENDDO
      ENDDO
*      print*,xq1,xq2,xji(1,1),xji(1,2),xji(2,1),xji(2,2)
      DO I=1,3
       DO J=1,3
        CKMDAG_RDI(I,J)=DCONJG(CKM(1,I))*RDI(1,J)
     .                 +DCONJG(CKM(2,I))*RDI(2,J)
     .                 +DCONJG(CKM(3,I))*RDI(3,J)
       ENDDO
      ENDDO
      CF_1 =(DCONJG(CKMDAG_RDI(3,1))*DCONJG(CKM(1,2))
     .      +DCONJG(CKMDAG_RDI(3,2))*DCONJG(CKM(2,2)))
     .             /CB_H/(CKM(3,3)*DCONJG(CKM(3,2)))
      CF_2 =DCONJG(CKMDAG_RDI(3,3))/CB_H/CKM(3,3)
*C7_0_MS_CINO
      C7_0_MS_CINO=2.D0/3.D0*MW_H**2/MSQ**2*CDABS(C_R(1,1))**2
     .             *B2SG_F7_1(XQ1)
     .            +2.D0/3.D0*MW_H**2/MSQ**2*CDABS(C_R(2,1))**2
     .             *B2SG_F7_1(XQ2)
     .            -CF_1
     .             *C_L(1,2)*DCONJG(C_R(1,1))*MW_H/DSQRT(2.D0)/M_C(1)
     .             *B2SG_F7_3(XQ1)
     .            -CF_1
     .             *C_L(2,2)*DCONJG(C_R(2,1))*MW_H/DSQRT(2.D0)/M_C(2)
     .             *B2SG_F7_3(XQ2)
*       print*,0,0,c7_0_ms_cino
      DO I=1,2
       DO J=1,2
        C7_0_MS_CINO=C7_0_MS_CINO
     .  -2.D0/3.D0*CDABS( C_R(I,1)*DCONJG(STMIX(1,J))
     .   -MT_MST/DSQRT(2.D0)/SB_H/MW_H*C_R(I,2)*DCONJG(STMIX(2,J)) )**2
     .   *MW_H**2/STMASS(J)**2*B2SG_F7_1(XJI(J,I))
     .  +CF_2*( -C_L(I,2)*DCONJG(C_R(I,1))*MW_H/DSQRT(2.D0)/M_C(I)
     .           *CDABS(STMIX(1,J))**2
     .          +DCONJG(STMIX(1,J))*STMIX(2,J)
     .           *C_L(I,2)*DCONJG(C_R(I,2))*MT_MST/2.D0/SB_H/M_C(I)
     .        )*B2SG_F7_3(XJI(J,I))
*       print*,i,j,c7_0_ms_cino
       ENDDO
      ENDDO
*
*C8_0_MS_CINO
      C8_0_MS_CINO=2.D0/3.D0*MW_H**2/MSQ**2*CDABS(C_R(1,1))**2
     .             *B2SG_F8_1(XQ1)
     .            +2.D0/3.D0*MW_H**2/MSQ**2*CDABS(C_R(2,1))**2
     .             *B2SG_F8_1(XQ2)
     .            -CF_1
     .             *C_L(1,2)*DCONJG(C_R(1,1))*MW_H/DSQRT(2.D0)/M_C(1)
     .             *B2SG_F8_3(XQ1)
     .            -CF_1
     .             *C_L(2,2)*DCONJG(C_R(2,1))*MW_H/DSQRT(2.D0)/M_C(2)
     .             *B2SG_F8_3(XQ2)
*       print*,0,0,c8_0_ms_cino
      DO I=1,2
       DO J=1,2
        C8_0_MS_CINO=C8_0_MS_CINO
     .  -2.D0/3.D0*CDABS( C_R(I,1)*DCONJG(STMIX(1,J))
     .   -MT_MST/DSQRT(2.D0)/SB_H/MW_H*C_R(I,2)*DCONJG(STMIX(2,J)) )**2
     .   *MW_H**2/STMASS(J)**2*B2SG_F8_1(XJI(J,I))
     .  +CF_2*( -C_L(I,2)*DCONJG(C_R(I,1))*MW_H/DSQRT(2.D0)/M_C(I)
     .           *CDABS(STMIX(1,J))**2
     .          +DCONJG(STMIX(1,J))*STMIX(2,J)
     .           *C_L(I,2)*DCONJG(C_R(I,2))*MT_MST/2.D0/SB_H/M_C(I)
     .        )*B2SG_F8_3(XJI(J,I))
*       print*,i,j,XJI(J,I),B2SG_F8_3(XJI(J,I))
*       print*,i,j,c8_0_ms_cino
       ENDDO
      ENDDO
*
*      print*,'C7(MS)[Ch]^0',C7_0_MS_CINO
*      print*,'C8(MS)[Ch]^0',C8_0_MS_CINO
*
*C7_0_MW_CINO,C8_0_MW_CINO
      ETA_S=AS_MST/AS_MW
*      print*,as_mw,as_mz,as_mt,as_mst,eta_s
      C7_0_MW_CINO=ETA_S**(16.D0/3.D0/7.D0)*C7_0_MS_CINO
     .            +8.D0/3.D0*(ETA_S**(14.D0/3.D0/7.D0)
     .                       -ETA_S**(16.D0/3.D0/7.D0))*C8_0_MS_CINO
      C8_0_MW_CINO=ETA_S**(14.D0/3.D0/7.D0)*C8_0_MS_CINO
*      print*,'C7(MW)[Ch]^0',C7_0_MW_CINO
*      print*,'C8(MW)[Ch]^0',C8_0_MW_CINO
*
      DO ISM=0,2  ! 0 = SM, 1=SM+C.Higgs, 2=SM+C.Higgs+C.ino
      IF(ISM.EQ.0) THEN
       ICH=0
       ICN=0
      ELSEIF(ISM.EQ.1) THEN
       ICH=1
       ICN=0
      ELSEIF(ISM.EQ.2) THEN
       ICH=1
       ICN=1
      ENDIF
*-----------------------------------------------------------------------
*C7_[0,1]_MW and C8_[0,1]_MW at MW
      C7_0_MW=C7_0_MW_SM+DBLE(ICH)*C7_0_MW_CH+DBLE(ICN)*C7_0_MW_CINO
      C8_0_MW=C8_0_MW_SM+DBLE(ICH)*C8_0_MW_CH+DBLE(ICN)*C8_0_MW_CINO
      C7_1_MW=C7_1_MW_SM+DBLE(ICH)*C7_1_MW_CH
      C8_1_MW=C8_1_MW_SM+DBLE(ICH)*C8_1_MW_CH
*      print*,'C7(MW)^0',C7_0_MW
*      print*,'C8(MW)^0',C8_0_MW
*      print*,'C7(MW)^1',C7_1_MW
*      print*,'C8(MW)^1',C8_1_MW
      CALL MAGIC_MISIAK(ISM,AS_MW,AS_160
     .         ,C7_0_MW_SM,C8_0_MW_SM,C7_1_MW_SM,C8_1_MW_SM
     .         ,C7_0_MW,C8_0_MW,C7_1_MW,C8_1_MW)
*-----------------------------------------------------------------------
*C2_0, C7_0 and C8_0 at RMU
      ETA=AS_MW/AS_RMU
*      print*,eta
      MAGIC_HA=H_MAGIC(1)*ETA**(A_MAGIC(1))
      DO I=2,8
       MAGIC_HA=MAGIC_HA+H_MAGIC(I)*ETA**(A_MAGIC(I))
      ENDDO
      MAGIC_H8B=H8_MAGIC(1)*ETA**(B_MAGIC(1))
      DO I=2,4
       MAGIC_H8B=MAGIC_H8B+H8_MAGIC(I)*ETA**(B_MAGIC(I))
      ENDDO
*      print*,magic_ha,magic_h8b
*
      C2_0=(ETA**(-12.D0/23.D0)+ETA**(6.D0/23.D0))/2.D0
      C7_0=ETA**(16.D0/23.D0)*C7_0_MW
     .    +8.D0/3.D0*(ETA**(14.D0/23.D0)-ETA**(16.D0/23.D0))*C8_0_MW
     .    +MAGIC_HA
      C8_0=ETA**(14.D0/23.D0)*(C8_0_MW+313063.D0/363036.D0)+MAGIC_H8B
*      print*,'C2_0 (mu_b) ',ich,icn,c2_0
*      print*,'C7_0 (mu_b) ',ich,icn,c7_0
*      print*,'C8_0 (mu_b) ',ich,icn,c8_0
*JSL 05/FEB/2009
      IF(ISM.EQ.0) THEN
       CAUX_H(160)=C2_0
       CAUX_H(161)=C7_0
       CAUX_H(162)=C8_0
      ELSEIF(ISM.EQ.1) THEN
       CAUX_H(163)=C2_0
       CAUX_H(164)=C7_0
       CAUX_H(165)=C8_0
      ELSEIF(ISM.EQ.2) THEN
       CAUX_H(166)=C2_0
       CAUX_H(167)=C7_0
       CAUX_H(168)=C8_0
      ENDIF
*-----------------------------------------------------------------------
*C7_1 at RMU
*      print*,eta
*      print*,xtw,b2sg_e(xtw)
       MAGIC_EFGA=( E_MAGIC(1)*ETA*B2SG_E(XTW)
     .             +F_MAGIC(1)+G_MAGIC(1)*ETA )*ETA**(A_MAGIC(1))
      DO I=2,8
       MAGIC_EFGA=MAGIC_EFGA
     .           +( E_MAGIC(I)*ETA*B2SG_E(XTW)
     .             +F_MAGIC(I)+G_MAGIC(I)*ETA )*ETA**(A_MAGIC(I))
      ENDDO
*      print*,magic_efga
*
      C7_1=ETA**(39.D0/23.D0)*C7_1_MW
     .    +8.D0/3.D0*(ETA**(37.D0/23.D0)-ETA**(39.D0/23.D0))*C8_1_MW
     .    +( 297664.D0/14283.D0*ETA**(16.D0/23.D0)
     .      -7164416.D0/357075.D0*ETA**(14.D0/23.D0)
     .      +256868.D0/14283.D0*ETA**(37.D0/23.D0)
     .      -6698884.D0/357075.D0*ETA**(39.D0/23.D0) )*C8_0_MW
     .      +37208.D0/4761.D0
     .        *(ETA**(39.D0/23.D0)-ETA**(16.D0/23.D0))*C7_0_MW
     .      +MAGIC_EFGA
*      print*,'C7_1 (mu_b) ',c7_1
*-----------------------------------------------------------------------
*C7_EM at RMU
*      print*,eta
      C7_EM=( 32.D0/75.D0*ETA**(-9.D0/23.D0)
     .       -40.D0/69.D0*ETA**(-7.D0/23.D0)
     .       +88.D0/575.D0*ETA**(16.D0/23.D0) )*C7_0_MW
     .     +(-32.D0/575.D0*ETA**(-9.D0/23.D0)
     .       +32.D0/1449.D0*ETA**(-7.D0/23.D0)
     .       +640.D0/1449.D0*ETA**(14.D0/23.D0)
     .       -704.D0/1725.D0*ETA**(16.D0/23.D0) )*C8_0_MW
     .     -190.D0/8073.D0*ETA**(-35.D0/23.D0)
     .     -359.D0/3105.D0*ETA**(-17.D0/23.D0)
     .     +4276.D0/121095.D0*ETA**(-12.D0/23.D0)
     .     +350531.D0/1009125.D0*ETA**(-9.D0/23.D0)
     .     +2.D0/4347.D0*ETA**(-7.D0/23.D0)
     .     -5956.D0/15525.D0*ETA**(6.D0/23.D0)
     .     +38380.D0/169533.D0*ETA**(14.D0/23.D0)
     .     -748.D0/8625.D0*ETA**(16.D0/23.D0)
*      print*,'C7_em(mu_b) ',c7_em
*-----------------------------------------------------------------------
*==>K_NLO_0
      K_NLO_0=K22*DREAL(C2_0*DCONJG(C2_0))
     .       +K27*DREAL(C2_0*DCONJG(C7_0))
     .       +K28*DREAL(C2_0*DCONJG(C8_0))
     .       +K77*DREAL(C7_0*DCONJG(C7_0))
     .       +K78*DREAL(C7_0*DCONJG(C8_0))
     .       +K88*DREAL(C8_0*DCONJG(C8_0))
*      print*,'K_NLO[0] ',k_nlo_0
*
*==>K_NLO_1
*      print*,sudakov,as_rmu
      K_NLO_1=SUDAKOV*AS_RMU/2.D0/PI*DREAL(C7_1*DCONJG(C7_0))
*      print*,'K_NLO[1] ',k_nlo_1
*
*==>K_NLO_EM
      AEM_LOW=1.D0/137.036D0
      XKSL_EM_1=2.D0*AS_RMU/PI*DLOG(MW_H/RMU)
      XKSL_EM_2=12.D0/23.D0*(1.D0/ETA-1.D0)
*      print*,xksl_em_1,xksl_em_2
      XKSL_EM=XKSL_EM_1
      K_NLO_EM=SUDAKOV*AEM_LOW/AS_RMU*(
     .         2.D0*DREAL(C7_EM*DCONJG(C7_0))
     .        -XKSL_EM*CDABS(C7_0)**2 )
*      print*,'K_NLO[em]',k_nlo_em
*
*==>K_NLO
      K_NLO=K_NLO_0+K_NLO_1+K_NLO_EM
*      print*,'K_NLO    ',k_nlo
*-----------------------------------------------------------------------
      B2CEN =0.1045D0
      BR_FAC=6.D0*AEM_LOW/PI/FZ0*B2CEN
     .      *CDABS(DCONJG(CKM(3,2))*CKM(3,3)/CKM(2,3))**2
      BR_NLO_0=BR_FAC*K_NLO_0
      BR_NLO  =BR_FAC*K_NLO
*      print*,'CKM factor       ',CDABS( DCONJG(CKM(3,2))
*     .                                 *CKM(3,3)/CKM(2,3) )
      IF (ISM.EQ.0) BR_SM      =BR_NLO
      IF (ISM.EQ.1) BR_SM_CH   =BR_NLO
      IF (ISM.EQ.2) BR_SM_CH_CN=BR_NLO
*-----------------------------------------------------------------------
*ACP:Kagan & Neubert PRD58(1998)094012
*      print*,DEL_GAM,Z_FIJ,dsqrt(z_fij)
      VZ=(5.D0+DLOG(Z_FIJ)+(DLOG(Z_FIJ))**2-Pi**2/3.D0)
     .  +((DLOG(Z_FIJ))**2-PI**2/3.D0)*Z_FIJ
     .  +(28.D0/9.D0-4.D0/3.D0*DLOG(Z_FIJ))*Z_FIJ**2
*       print*,vz
      BZ=B2SG_GG(1.D0)-B2SG_GG(1.D0-DEL_GAM)
*       print*,B2SG_GG(1.D0),B2SG_GG(1.D0-DEL_GAM),bz
*      print*,as_mt,as_mz,as_mb
      A27=AS_MB*(40.D0/81.D0-8.D0*Z_FIJ/9.D0*(VZ+BZ))
      A87=-4.D0/9.D0*AS_MB
      A28=8.D0/27.D0*AS_MB*Z_FIJ*BZ
*      print*,a27,a87,a28
      ACP=( A27*DIMAG(C2_0*DCONJG(C7_0))
     .     +A87*DIMAG(C8_0*DCONJG(C7_0))
     .     +A28*DIMAG(C2_0*DCONJG(C8_0))
     .    )/CDABS(C7_0)**2
*      print*,'ACP1',A27*DIMAG(C2_0*DCONJG(C7_0))/CDABS(C7_0)**2
*      print*,'ACP2',A87*DIMAG(C8_0*DCONJG(C7_0))/CDABS(C7_0)**2
      IF (ISM.EQ.0) ACP_SM      =ACP
      IF (ISM.EQ.1) ACP_SM_CH   =ACP
      IF (ISM.EQ.2) ACP_SM_CH_CN=ACP
*-----------------------------------------------------------------------
      ENDDO ! ISM
*
      RAUX_H(135)=BR_SM_CH_CN*1.0D4
      RAUX_H(136)=ACP_SM_CH_CN*1.0D2
*-----------------------------------------------------------------------
      IF(IFLAG(17).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,'                     B -> X_s gamma'
       write(*,31) DEL_GAM,(1.D0-DEL_GAM)*MB_POLE/2.D0
      print*,'---------------------------------------------------------'
       print*,' b-q masses [GeV]  (pole,       @mb^pole,   @mt^pole):'
       write(*,32) MB_POLE,MB_MB,MBMT_H
       print*,' c-q masses [GeV]  (pole,       @mc^pole,   @mb^pole):'
       write(*,33) MC_POLE,MC_MC,MC_MB
       write(*,34) RMU,MU_C
      print*,'---------------------------------------------------------'
       write(*,35) BR_SM_CH_CN*1.d4
       write(*,36) BR_SM_CH*1.d4
       write(*,37) BR_SM*1.d4
       write(*,38) ACP_SM_CH_CN*1.d2
      print*,'---------------------------------------------------------'
      print*,' '
      ENDIF
*-----------------------------------------------------------------------
  11  format(1x,'/',3(1x,'(',e10.4,1x,e10.4,')',1x),'\\')
  12  format(1x,'|',3(1x,'(',e10.4,1x,e10.4,')',1x),'|')
  13  format(1x,'\\',3(1x,'(',e10.4,1x,e10.4,')',1x),'/')
  14  format(2x,3(1x,'(',e10.4,1x,e10.4,')',1x))
  21  format(1x,'/',3(1x,'(',e10.4,')',1x),'\\')
  22  format(1x,'|',3(1x,'(',e10.4,')',1x),'|')
  23  format(1x,'\\',3(1x,'(',e10.4,')',1x),'/')
  31  format(2x,'delta and E_gamma^cut [GeV]:',2(1x,e10.4,1x))
  32  format(20x,3(1x,e10.4,1x))
  33  format(20x,3(1x,e10.4,1x))
  34  format(2x,'mu_b and mu_c  [GeV]       :',2(1x,e10.4,1x))
  35  format(2x,'BR  x 10^4:',1x,e10.4,1x,'(SM+Charged Higgs+Chargino)')
  36  format(13x,'[',e10.4,1x,'(SM+Charged Higgs)]')
  37  format(13x,'[',e10.4,1x,'(SM)]')
  38  format(2x,'ACP x 10^2:',1x,e10.4,1x,'%')
*
      RETURN
      END

      REAL*8 FUNCTION XLI2(X)
************************************************************************
*
* Real part of the function Li_2(x) with real x togerther with YLI2
*
* "The Dilogarithm Function of a Real Argument" by R. Morris
*    MATHEMATICS OF COMPUTATION, V33, N146 (1979 Apr.) pp. 778-787
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      PI=2.D0*DASIN(1.D0)
*
      IF(X.GT.2.D0) THEN                     ! 2 <    X < infinity
       XLI2=PI**2/3.D0-(DLOG(X))**2/2.D0-YLI2(1.D0/X)
      ELSEIF(X.GT.1.D0.AND.X.LE.2.D0) THEN   ! 1 <    X <= 2 [typo corrected]
       XLI2=PI**2/6.D0-DLOG(X)*(DLOG(X-1.D0)-DLOG(X)/2.D0)
     .     +YLI2(1.D0-1.D0/X)
      ELSEIF(X.EQ.1.D0) THEN                 !        X = 1
       XLI2=PI**2/6.D0
      ELSEIF(X.GT.0.5D0.AND.X.LT.1.D0) THEN  ! 1/2 <  X <  1
       XLI2=PI**2/6.D0-DLOG(X)*DLOG(1.0D0-X)-YLI2(1.0D0-X)
      ELSEIF(X.GT.0.D0.AND.X.LE.0.5D0) THEN  ! 0   <  X <= 1/2
       XLI2=YLI2(X)
      ELSEIF(ABS(X).EQ.0.D0) THEN            !        X = 0
       XLI2=0.D0
      ELSEIF(X.GT.-1.D0.AND.X.LT.0.D0) THEN  !-1   <  X < 0
       XLI2=-(DLOG(1.0D0-X))**2/2.D0-YLI2(X/(X-1.D0))
      ELSEIF(X.EQ.-1.D0) THEN                !        X = -1
       XLI2=-PI**2/12.D0
      ELSEIF(X.LT.-1.D0) THEN                !        X < -1
       XLI2=-PI**2/6.D0
     .      -DLOG(1.D0-X)/2.D0*(2.D0*DLOG(-X)-DLOG(1.D0-X))
     .      +YLI2(1.D0/(1.D0-X))
      ENDIF
*
      RETURN
      END

      REAL*8 FUNCTION YLI2(X)
************************************************************************
*
* Li_2(x) = sum[k=1..infinity] x^k/k**2 with real x and |x| <= 1/2
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      PI=2.D0*DASIN(1.D0)
*
      IF(ABS(X).LE.0.5D0) THEN
       YLI2=X
       DO K=2,1000
        DYLI2=X**K/DBLE(K)**2
        YLI2=YLI2+DYLI2
        IF(DYLI2.LT.1.D-10) GOTO 88
       ENDDO
      ELSEIF(ABS(X).GT.0.5D0) THEN
       print*,'ERROR: YLI2: |X| should be less than 1/2',X
       RETURN
      ENDIF
*
 88   CONTINUE
*
      RETURN
      END



      REAL*8 FUNCTION B2SG_F22(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /FIJ_BODE/ DGAM_FIJ,Z_FIJ
      COMPLEX*16 XI,G
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      D=DGAM_FIJ
      Z=Z_FIJ
      T=X/Z
      IF(T.LT.4.D0) THEN
       G=-2.D0*(DATAN(DSQRT(T/(4.D0-T))))**2
      ELSEIF(T.GE.4.D0) THEN
       G=2.D0*(DLOG((DSQRT(T)+DSQRT(T-4.D0))/2.D0)-XI*PI/2.D0)**2
      ELSE
       print*,'ERROR: B2SG_F22'
       RETURN
      ENDIF
*
      XD=DMAX1(X,1.D0-D)
*
      B2SG_F22=16.D0/27.D0*(1.D0-X)*(1.D0-XD)*CDABS(Z/X*G+0.5D0)**2
*      B2SG_F22=DLOG(X)
*      print*,X,1.D0-D,XD,Z,T,G,B2SG_F22
*      print*,X,T,G,B2SG_F22
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F27(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /FIJ_BODE/ DGAM_FIJ,Z_FIJ
      COMPLEX*16 XI,G
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      D=DGAM_FIJ
      Z=Z_FIJ
      T=X/Z
      IF(T.LT.4.D0) THEN
       G=-2.D0*(DATAN(DSQRT(T/(4.D0-T))))**2
      ELSEIF(T.GE.4.D0) THEN
       G=2.D0*(DLOG((DSQRT(T)+DSQRT(T-4.D0))/2.D0)-XI*PI/2.D0)**2
      ELSE
       print*,'ERROR: B2SG_F27'
       RETURN
      ENDIF
*
      XD=DMAX1(X,1.D0-D)
*
      B2SG_F27=-8.D0*Z/9.D0*(1.D0-XD)*DREAL(G+X/2.D0/Z)
*      B2SG_F27=DLOG(X)
*      print*,X,1.D0-D,XD,Z,T,G,B2SG_F27
*      print*,X,T,G,B2SG_F27
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F7_1(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F7_1=X*(7.D0-5.D0*X-8.D0*X**2)/24.D0/(X-1.D0)**3
     .    +X**2*(3.D0*X-2.D0)/4.D0/(X-1.D0)**4*DLOG(X)
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F7_2(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F7_2=X*(3.D0-5.D0*X)/12.D0/(X-1.D0)**2
     .    +X*(3.D0*X-2.D0)/6.D0/(X-1.D0)**3*DLOG(X)
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F7_3(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F7_3=(5.D0-7.D0*X)/6.D0/(X-1.D0)**2
     .    +X*(3.D0*X-2.D0)/3.D0/(X-1.D0)**3*DLOG(X)
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F8_1(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F8_1=X*(2.D0+5.D0*X-X**2)/8.D0/(X-1.D0)**3
     .    -3.D0*X**2/4.D0/(X-1.D0)**4*DLOG(X)
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F8_2(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F8_2=X*(3.D0-X)/4.D0/(X-1.D0)**2
     .    -X/2.D0/(X-1.D0)**3*DLOG(X)
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_F8_3(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_F8_3=(1.D0+X)/2.D0/(X-1.D0)**2
     .    -X/(X-1.D0)**3*DLOG(X)
*
      RETURN
      END



      REAL*8 FUNCTION B2SG_G7_1(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_G7_1=(-16.D0*X**4-122.D0*X**3+80.D0*X**2-8.D0*X)
     .          /9.D0/(X-1.D0)**4*XLI2(1.D0-1.D0/X)
     .         +(6.D0*X**4+46.D0*X**3-28.D0*X**2)
     .          /3.D0/(X-1.D0)**5*(DLOG(X))**2
     .         +(-102.D0*X**5-588.D0*X**4-2262.D0*X**3+3244.D0*X**2
     .           -1364.D0*X+208.D0)
     .          /81.D0/(X-1.D0)**5*DLOG(X)
     .         +(1646.D0*X**4+12205.D0*X**3-10740.D0*X**2
     .           +2509.D0*X-436.D0)
     .          /486.D0/(X-1.D0)**4
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_G8_1(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_G8_1=(-4.D0*X**4+40.D0*X**3+41.D0*X**2+X)
     .          /6.D0/(X-1.D0)**4*XLI2(1.D0-1.D0/X)
     .         +(-17.D0*X**3-31.D0*X**2)
     .          /2.D0/(X-1.D0)**5*(DLOG(X))**2
     .         +(-210.D0*X**5+1086.D0*X**4+4893.D0*X**3+2857.D0*X**2
     .           -1994.D0*X+280.D0)
     .          /216.D0/(X-1.D0)**5*DLOG(X)
     .         +(737.D0*X**4-14102.D0*X**3-28209.D0*X**2
     .           +610.D0*X-508.D0)
     .          /1296.D0/(X-1.D0)**4
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_E(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      B2SG_E=X*(18.D0-11.D0*X-X**2)/12.D0/(1.D0-X)**3
     .      +X**2*(15.D0-16.D0*X+4.D0*X**2)/6.D0/(1.D0-X)**4*DLOG(X)
     .      -2.D0/3.D0*DLOG(X)
*
      RETURN
      END

      COMPLEX*16 FUNCTION B2SG_G7_H(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*
      B2SG_G7_H=-AD_CH*4.D0/3.D0*X*(
     .          4.D0*(-3.D0+7.D0*X-2.D0*X**2)
     .          /3.D0/(X-1.D0)**3*XLI2(1.D0-1.D0/X)
     .         +(8.D0-14.D0*X-3.D0*X**2)
     .          /3.D0/(X-1.D0)**4*(DLOG(X))**2
     .         +2.D0*(-3.D0-X+12.D0*X**2-2.D0*X**3)
     .          /3.D0/(X-1.D0)**4*DLOG(X)
     .         +(7.D0-13.D0*X+2.D0*X**2)
     .          /(X-1.D0)**3 )
     .         +AU_CH*2.D0/9.D0*X*(
     .          X*(18.D0-37.D0*X+8.D0*X**2)
     .          /(X-1.D0)**4*XLI2(1.D0-1.D0/X)
     .         +X*(-14.D0+23.D0*X+3.D0*X**2)
     .          /(X-1.D0)**5*(DLOG(X))**2
     .         +(-50.D0+251.D0*X-174.D0*X**2-192.D0*X**3+21.D0*X**4)
     .          /9.D0/(X-1.D0)**5*DLOG(X)
     .         +(797.D0-5436.D0*X+7569.D0*X**2-1202.D0*X**3)
     .          /108.D0/(X-1.D0)**4 )
*
*      print*,'>>>>>>>>>>>>G7_H',au_ch,ad_ch
*
      RETURN
      END

      COMPLEX*16 FUNCTION B2SG_G8_H(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*
      B2SG_G8_H=-AD_CH*1.D0/3.D0*X*(
     .          (-36.D0+25.D0*X-17.D0*X**2)
     .          /2.D0/(X-1.D0)**3*XLI2(1.D0-1.D0/X)
     .         +(19.D0+17.D0*X)
     .          /(X-1.D0)**4*(DLOG(X))**2
     .         +(-3.D0-187.D0*X+12.D0*X**2-14.D0*X**3)
     .          /4.D0/(X-1.D0)**4*DLOG(X)
     .         +3.D0*(143.D0-44.D0*X+29.D0*X**2)
     .          /8.D0/(X-1.D0)**3 )
     .         +AU_CH*1.D0/6.D0*X*(
     .          X*(30.D0-17.D0*X+13.D0*X**2)
     .          /(X-1.D0)**4*XLI2(1.D0-1.D0/X)
     .         -X*(31.D0+17.D0*X)
     .          /(X-1.D0)**5*(DLOG(X))**2
     .         +(-226.D0+817.D0*X+1353.D0*X**2+318.D0*X**3+42.D0*X**4)
     .          /36.D0/(X-1.D0)**5*DLOG(X)
     .         +(1130.D0-18153.D0*X+7650.D0*X**2-4451.D0*X**3)
     .          /216.D0/(X-1.D0)**4 )
*
*      print*,'G8_H',au_ch,ad_ch
*
      RETURN
      END

      COMPLEX*16 FUNCTION B2SG_D7_H(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*
      B2SG_D7_H=-AD_CH*2.D0/9.D0*X*(
     .          (21.D0-47.D0*X+8.D0*X**2)
     .          /(X-1.D0)**3
     .         +2.D0*(-8.D0+14.D0*X+3.D0*X**2)
     .          /(X-1.D0)**4*DLOG(X) )
     .         +AU_CH*2.D0/9.D0*X*(
     .          (-31.D0-18.D0*X+135.D0*X**2-14.D0*X**3)
     .          /6.D0/(X-1.D0)**4
     .         +X*(14.D0-23.D0*X-3.D0*X**2)
     .          /(X-1.D0)**5*DLOG(X) )
*
      RETURN
      END


      COMPLEX*16 FUNCTION B2SG_D8_H(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*
      B2SG_D8_H=-AD_CH*1.D0/3.D0*X*(
     .          (81.D0-16.D0*X+7.D0*X**2)
     .          /2.D0/(X-1.D0)**3
     .         -(19.D0+17.D0*X)
     .          /(X-1.D0)**4*DLOG(X) )
     .         +AU_CH*1.D0/6.D0*X*(
     .          (-38.D0-261.D0*X+18.D0*X**2-7.D0*X**3)
     .          /6.D0/(X-1.D0)**4
     .         +X*(31.D0+17.D0*X)
     .          /(X-1.D0)**5*DLOG(X) )
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_EH(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 AU_CH,AD_CH
      COMMON /C78_1_CH/ AU_CH,AD_CH
*
      B2SG_EH=AU_CH*(
     .        X*(16.D0-29.D0*X+7.D0*X**2)/36.D0/(X-1.D0)**3
     .         +X*(-2.D0+3.D0*X)/6.D0/(X-1.D0)**4*DLOG(X) )
*
      RETURN
      END

      REAL*8 FUNCTION B2SG_GG(X)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /FIJ_BODE/ DGAM_FIJ,Z_FIJ
      COMPLEX*16 XI,G
*
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
      D=DGAM_FIJ
      Z=Z_FIJ
      Y=X
      IF(Y.LT.4.D0*Z) THEN
       B2SG_GG=0.D0
       RETURN
      ELSE
       B2SG_GG=(Y**2-4.D0*Y*Z+6.D0*Z**2)
     .        *DLOG(DSQRT(Y/4.D0/Z)+DSQRT(Y/4.D0/Z-1.D0))
     .        -3.D0*Y*(Y-2.D0*Z)/4.D0*DSQRT(1.D0-4.D0*Z/Y)
       RETURN
      ENDIF
*
      RETURN
      END

      SUBROUTINE MAGIC_MISIAK(ISM,AS_MW,AS_160
     .         ,C7_0_MW_SM,C8_0_MW_SM,C7_1_MW_SM,C8_1_MW_SM
     .         ,C7_0_MW,C8_0_MW,C7_1_MW,C8_1_MW)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 C7_0_MW,C8_0_MW,C7_1_MW,C8_1_MW
      COMPLEX*16 C70NP,C80NP,C71NP,C81NP
*
*      print*,'>> MISIAK <<'
*      print*,AS_MW,AS_160
*      print*,'C7(MW)^0',C7_0_MW,C7_0_MW_SM
*      print*,'C8(MW)^0',C8_0_MW,C8_0_MW_SM
*      print*,'C7(MW)^1',C7_1_MW,C7_1_MW_SM
*      print*,'C8(MW)^1',C8_1_MW,C8_1_MW_SM
*
*NEW-PHYSICS Wilson coefficients at 160 GeV
      ETA=AS_MW/AS_160
*
      C70NP=ETA**(16.D0/3.D0/7.D0)*(C7_0_MW-C7_0_MW_SM)
     .            +8.D0/3.D0*(ETA**(14.D0/3.D0/7.D0)
     .                       -ETA**(16.D0/3.D0/7.D0))
     .                      *(C8_0_MW-C8_0_MW_SM)
      C80NP=ETA**(14.D0/3.D0/7.D0)*(C8_0_MW-C8_0_MW_SM)
*
      C71NP=ETA**(16.D0/3.D0/7.D0)*(C7_1_MW-C7_1_MW_SM)
     .            +8.D0/3.D0*(ETA**(14.D0/3.D0/7.D0)
     .                       -ETA**(16.D0/3.D0/7.D0))
     .                      *(C8_1_MW-C8_1_MW_SM)
      C81NP=ETA**(14.D0/3.D0/7.D0)*(C8_1_MW-C8_1_MW_SM)
*
*BR
      BR=( 3.150055981433902D0
     .   - 7.824786698647555D0*DREAL(c70np)
     .   - 2.0945629943122532D0*DREAL(c80np)
     .   + 5.133484906701248D0*CDABS(c70np)**2
     .   + 0.518633676234877D0*CDABS(c80np)**2
     .   - 0.7501776763621714D0*DIMAG(c70np)
     .   + 0.09209321186259296D0*DIMAG(c70np)*DIMAG(c71np)
     .   + 0.5977860427960031D0*DIMAG(c80np)
     .   + 2.5928603736700957D0*DIMAG(c70np)*DIMAG(c80np) 
     .   + 0.027487426566116562D0*DIMAG(c71np)*DIMAG(c80np)
     .   + 0.027487426566116562D0*DIMAG(c70np)*DIMAG(c81np)
     .   + 0.008204281335686028D0*DIMAG(c80np)*DIMAG(c81np)
     .   - 0.9686837773848204D0*DIMAG(c80np)*DREAL(c70np)
     .   - 0.08142794331649322D0*DREAL(c71np)
     .   + 0.09209321186259296D0*DREAL(c70np)*DREAL(c71np)
     .   + 0.9686837773848204D0*DIMAG(c70np)*DREAL(c80np)
     .   + 2.5928603736700957D0*DREAL(c70np)*DREAL(c80np)
     .   + 0.027487426566116562D0*DREAL(c71np)*DREAL(c80np)
     .   - 0.02430412152072149D0*DREAL(c81np)
     .   + 0.027487426566116562D0*DREAL(c70np)*DREAL(c81np)
     .   + 0.008204281335686028D0*DREAL(c80np)*DREAL(c81np) )
*
*      print*,'>> MISIAK <<',ism,br
*
      RETURN
      END


