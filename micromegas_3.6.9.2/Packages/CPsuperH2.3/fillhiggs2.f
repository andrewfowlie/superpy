      SUBROUTINE FILLHIGGS2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .                     ,MCH,HMASS,OMIX_OUT)
************************************************************************
* ---
* |1| RAUX_H(10) is reserved for the charged Higgs boson POLE mass
* --- when IFLAG_H(11)=0 but an effective potential mass when
*     IFLAG_H(11)=1.

* ---
* |2| New flag IFLAG_H(12)
* ---
*    1 : same as the old one but without CP phase regulator
*    2 : 1+Including threshold corrections to CPI_Q
*    3 : 1+Including threhsold corrections to X1L1-X1L4
*    4 : 1+Pole mass improvement
*    5 : full improvement 1+2+3+4+5
*
* ---
* |3| New flag IFLAG_H(60)
* ---
*    IFLAG_H(60)=1 returned when iteration for the on-shell Higgs masses
*    fails.
* ---
* |4| The plan of calculating the Higgs masses depending on IFLAG_H(12) 
* --- is as follows:
*    
*    <i>   Call GET_MASQ for MA^2(s=MCH^2)
*
*    <ii>  Call GET_CMNH for CMNH[4,4](s=0), the neutral Higgs-boson 
*          mass squared matrix at s=0. As a result we have EP3[3] and
*          OMIX_0[3,3]
*
*    Then we treat two cases seperately depending on IFLAG_H(12)
*
*    <iii> -+=>> When IFLAG_H(12)=1,2,3:
*          Call GET_CMNH for CMNH[4,4](s=EP3[3]^2) and, as a result, we
*          have HP3[3] together with OMIX_0[3,3] taking the way used in
*          the previous FILLHIGGS routine
*
*          -+=>> When IFLAG_H(12)=4,5:
*          Call GET_CMNH for CMNH[4,4](s=HP3[3]^2) and, as a result, we
*          have HP3[3] together with OMIX_{1,2,3}[3,3]. The numerical
*          iteration is needed
*
*    Then finally,
*
*    <iv>  -+=>> If IFLAG_H(11)=1:
*          HMASS[3]=EP[3], MCH_OUT=MCH^eff., OMIX_OUT[3,3]=OMIX_0[3,3]
*
*
*          -+=>> If IFLAG_H(11)=0:
*          HMASS[3]=HP[3], MCH_OUT=MCH^pole, OMIX_OUT[3,3]=OMIX_0[3,3]
*          In this case, the reordering may need.
*
*    N.B. We are returning the mixing matrix at s=0 for, specifically, 
*         the couplings of the Higgs bosons
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8    SMPARA_H(NSMIN),SSPARA_H(NSSIN)
      INTEGER*8 IFLAG_H(NFLAG)
      REAL*8    HMASS(3),OMIX_OUT(3,3)
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
*Local
      COMPLEX*16 CDELHB,CDELHB1,CDELHB2
      COMPLEX*16 CDELHT,CDELHT1,CDELHT2
      COMPLEX*16 HT_CP,HB_CP
      COMPLEX*16 CMNH(4,4),CTMP3(3,3)
      COMPLEX*16 CTR23,CDET23,CD23
*
      REAL*8     NH3(3,3),EV3(3),AUX3(3)
      REAL*8     EP3(3),HP3(3)
      REAL*8     HP3_TMP(3),OMIX_TMP(3,3)
      REAL*8     OMIX_0(3,3),OMIX_1(3,3),OMIX_2(3,3),OMIX_3(3,3)
*
      REAL*8     DMH3(3,3),DET_DMH3
*
      PI=2.D0*DASIN(1.D0)
*
*-------------------------------------------------------------------------
*Several SFermion mass scales
*-------------------------------------------------------------------------
* 
      QQT2=SSPARA_H(11)**2+MTPOLE_H**2
      QTT2=SSPARA_H(12)**2+MTPOLE_H**2
      QST2=DMAX1(QQT2,QTT2)
      QQB2=SSPARA_H(11)**2+MBMT_H**2
      QBB2=SSPARA_H(13)**2+MBMT_H**2
      QSB2=DMAX1(QQB2,QBB2)
      QSF2=DMAX1(QST2,QSB2)
*      print*,'>> Check 1 (NEW) : ',QST2,QSB2,QSF2
      RAUX_H(11)=QST2
      RAUX_H(12)=QSB2
      RAUX_H(13)=QSF2
*      print*,'Scales^2:',qst2,qsb2,qsf2
*
*-------------------------------------------------------------------------
* VEVs and Yukawa Couplings at SFermion scales without including the
* Threshold corrections
*-------------------------------------------------------------------------
*
*.....At Mt_pole :
      HTSM=DSQRT(2.D0)*MTMT_H/V_H
      HT  =DSQRT(2.D0)*MTMT_H/V_H/SB_H
      HB  =DSQRT(2.D0)*MBMT_H/V_H/CB_H
      GS2 =4.D0*PI*ASMT_H
      BTSM=1.D0/16.D0/PI**2*(9.D0*HTSM**2/2.D0-8.D0*GS2)
      BBSM=1.D0/16.D0/PI**2*(HTSM**2/2.D0-8.D0*GS2)
      BT=1.D0/16.D0/PI**2*(9.D0*HT**2/2.D0+HB**2/2.D0-8.D0*GS2)
      BB=1.D0/16.D0/PI**2*(9.D0*HB**2/2.D0+HT**2/2.D0-8.D0*GS2)
      V1=V_H*CB_H
      V2=V_H*SB_H
*      print*,'>> Check 2 (NEW) : ',HTSM,HT,HB,GS2

      IF(MCH.GT.MTPOLE_H) THEN ! MTpole < MCH < MSfermion
*                                two-step running from MTpole to Sfermion scale
*
*......VEVs at MCH
       TB_MCH=TB_H*(1.D0
     .       -3.D0*(HT**2-HB**2)/32.D0/PI**2*DLOG(MCH**2/MTPOLE_H**2))
       CB_MCH= 1.D0/DSQRT(1.D0+TB_MCH**2)
       SB_MCH= TB_MCH/DSQRT(1.D0+TB_MCH**2)
       XISM=1.D0+3.D0*HTSM**2/32.D0/PI**2*DLOG(MCH**2/MTPOLE_H**2)
       V_MCH=V_H/XISM
       V1_MCH=CB_MCH*V_MCH
       V2_MCH=SB_MCH*V_MCH
*       print*,'>> Check 3 (NEW) : ',TB_MCH,XISM,V_H,V1_MCH,V2_MCH
*......Yukawa Couplings at MCH
       HT_MCH=HT*(1.D0+2.D0*BTSM*DLOG(MCH**2/MTPOLE_H**2))**0.25D0
       HB_MCH=HB*(1.D0+2.D0*BBSM*DLOG(MCH**2/MTPOLE_H**2))**0.25D0
*......VEVs and anomalous-dim factors at Sfermion scale
       XI1H_ST=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QST2/MCH**2)
       XI1H_SB=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSB2/MCH**2)
       XI1H_SF=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSF2/MCH**2)
       XI2H_ST=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QST2/MCH**2)
       XI2H_SB=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSB2/MCH**2)
       XI2H_SF=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSF2/MCH**2)
       V1_ST=V1_MCH/XI1H_ST
       V1_SB=V1_MCH/XI1H_SB
       V1_SF=V1_MCH/XI1H_SF
       V2_ST=V2_MCH/XI2H_ST
       V2_SB=V2_MCH/XI2H_SB
       V2_SF=V2_MCH/XI2H_SF
       XI1_ST=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QST2/MTPOLE_H**2)
       XI1_SB=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSB2/MTPOLE_H**2)
       XI1_SF=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSF2/MTPOLE_H**2)
       XI2_ST=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QST2/MTPOLE_H**2)
       XI2_SB=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSB2/MTPOLE_H**2)
       XI2_SF=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSF2/MTPOLE_H**2)
*       print*,'>> Check 4 (NEW) : ',V1_ST,V1_SB,V1_SF,V2_ST,V2_SB,V2_SF
*......Yukawa Couplings and its CP phases at fermion and Sfermion scales
       HT_MT=HT
       HT_ST=HT_MCH*(1.D0+2.D0*BT*DLOG(QST2/MCH**2))**0.25D0
       HT_SF=HT_MCH*(1.D0+2.D0*BT*DLOG(QSF2/MCH**2))**0.25D0
       HT_CP=DCMPLX(1.D0,0.D0)
       HB_MT=HB
       HB_SB=HB_MCH*(1.D0+2.D0*BB*DLOG(QSB2/MCH**2))**0.25D0
       HB_SF=HB_MCH*(1.D0+2.D0*BB*DLOG(QSF2/MCH**2))**0.25D0
       HB_CP=DCMPLX(1.D0,0.D0)
*       print*,'>> Check 5 (NEW) : ',HT_ST,HT_SF,HB_SB,HB_SF
*
      ELSE      ! MCH < MTpole < M_Sfermion : One-step from Mtpole to Sfermion scales
*
*......VEVs and anomalous-dim factors at Sfermion scale
       XI1_ST=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QST2/MTPOLE_H**2)
       XI1_SB=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSB2/MTPOLE_H**2)
       XI1_SF=1.D0+3.D0*HB**2/32.D0/PI**2*DLOG(QSF2/MTPOLE_H**2)
       XI2_ST=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QST2/MTPOLE_H**2)
       XI2_SB=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSB2/MTPOLE_H**2)
       XI2_SF=1.D0+3.D0*HT**2/32.D0/PI**2*DLOG(QSF2/MTPOLE_H**2)
       V1_ST=V_H*CB_H/XI1_ST
       V1_SB=V_H*CB_H/XI1_SB
       V1_SF=V_H*CB_H/XI1_SF
       V2_ST=V_H*SB_H/XI2_ST
       V2_SB=V_H*SB_H/XI2_SB
       V2_SF=V_H*SB_H/XI2_SF
*......Yukawa Couplings and its CP phase at fermion and Sfermion scale
       HT_MT=HT
       HT_ST=HT*(1.D0+2.D0*BT*DLOG(QST2/MTPOLE_H**2))**0.25D0
       HT_SF=HT*(1.D0+2.D0*BT*DLOG(QSF2/MTPOLE_H**2))**0.25D0
       HT_CP=DCMPLX(1.D0,0.D0)
       HB_MT=HB
       HB_SB=HB*(1.D0+2.D0*BB*DLOG(QSB2/MTPOLE_H**2))**0.25D0
       HB_SF=HB*(1.D0+2.D0*BB*DLOG(QSF2/MTPOLE_H**2))**0.25D0
       HB_CP=DCMPLX(1.D0,0.D0)
*       print*,'>> Check 4 (NEW) : ',V1_ST,V1_SB,V1_SF,V2_ST,V2_SB,V2_SF
*       print*,'>> Check 5 (NEW) : ',HT_ST,HT_SF,HB_SB,HB_SF
*       print*,HB,BB,2.D0*BB*DLOG(QSF2/MTPOLE_H**2)
*
      ENDIF           ! IF( MCH > MTpole )
*
*      RAUX_H(14)=V1
*      RAUX_H(15)=V1_ST
*      RAUX_H(16)=V1_SB
*      RAUX_H(17)=V1_SF
*      print*,'V1: Before corrections',v1,v1_st,v1_sb,v1_sf
*      RAUX_H(18)=V2
*      RAUX_H(19)=V2_ST
*      RAUX_H(20)=V2_SB
*      RAUX_H(21)=V2_SF
*      print*,'V2: Before corrections',v2,v2_st,v2_sb,v2_sf
      RAUX_H(22)=HT_MT
      RAUX_H(23)=HB_MT
*      print*,'HT,HB (Before Corrections):',ht_mt,hb_mt
*-------------------------------------------------------------------------
* Threshold corrections when IFLAG_H(10)=0: Resummed by iteration with
* N_iteration <= 100. If N_iteration > 100, IFLAG_H(54)=1 returned
* The Yukawa couplings at Sfermions mass scales including the threshold
* corrections has been calculated. The results calculated without corrections
* in the above are used as inputs.
*-------------------------------------------------------------------------
*
*.....Some gauge couplings
      GWY2=GW_H**2/8.D0/CW_H**2
      GXT2=(GW_H**2-5.D0*GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0  
      GXB2=(GW_H**2-GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0  
*.....Alpha_S(Stop, Sbottom)
      B6=( 11.D0 - 2.D0*6.D0/3.D0 )/(4.D0*PI)
      AS_ST = ASMT_H/( 1.D0 + B6*ASMT_H*DLOG(QST2/MTPOLE_H**2) )
      AS_SB = ASMT_H/( 1.D0 + B6*ASMT_H*DLOG(QSB2/MTPOLE_H**2) )

*.....Not Corrected Yukawa couplings
      HB_NC=HB_SB
      HT_NC=HT_ST
*
*      print*,'HT (Before Corrections):',ht_mt,ht_st,ht_sf
*      print*,'HB (Before Corrections):',hb_mt,hb_sb,hb_sf
*
      IF(IFLAG_H(10).EQ.0) THEN
*       print*,'Threshold corrections ON:'
*.....ITeRations
      ITR=0
      ITRMAX=100
      EPS_TR=1.D-3
*      EPS_TR=1.D-9 ! to test precision
*.....An approximated maximal value is used for the initial b-quark Yukawa coupling.
*.....The results are independent of the initial values of the Yukawa couplings.
*JSL;03/MAR/2006, we add 1 GeV to AB to avoid pole in initial HB_R.
      IF(CDABS(CB_H*DCONJG(AB_H)-MU_H*SB_H).EQ.0.D0) THEN
       HB_R=DSQRT((SSPARA_H(11)**2+MBMT_H**2)*QSB2)
     .     /MTPOLE_H/CDABS(CB_H*DCONJG(AB_H+1.D0)-MU_H*SB_H) 
      ELSE
       HB_R=DSQRT((SSPARA_H(11)**2+MBMT_H**2)*QSB2)
     .     /MTPOLE_H/CDABS(CB_H*DCONJG(AB_H)-MU_H*SB_H) 
      ENDIF
*      HT_R=0.D0
      HT_R=HT_NC
*
 111  CONTINUE
      ITR=ITR+1
C.....stop and sbottom masses**2 (Note it's squared!)
 110  CONTINUE
      SUM_ST=SSPARA_H(11)**2+SSPARA_H(12)**2+HT_R**2*V2_ST**2
     .      +GWY2*(V1_ST**2-V2_ST**2)
      DIF_ST2=(SSPARA_H(11)**2-SSPARA_H(12)**2
     .            +GXT2*(V1_ST**2-V2_ST**2)/2.D0)**2
     .       +2.D0*HT_R**2*CDABS(DCONJG(AT_H)*V2_ST-MU_H*V1_ST)**2
      MSQ_ST2=(SUM_ST+DSQRT(DIF_ST2))/2.D0
      MSQ_ST1=(SUM_ST-DSQRT(DIF_ST2))/2.D0
      SUM_SB=SSPARA_H(11)**2+SSPARA_H(13)**2+HB_R**2*V1_SB**2
     .      +GWY2*(V2_SB**2-V1_SB**2)
      DIF_SB2=(SSPARA_H(11)**2-SSPARA_H(13)**2
     .            +GXB2*(V2_SB**2-V1_SB**2)/2.D0)**2
     .       +2.D0*HB_R**2*CDABS(DCONJG(AB_H)*V1_SB-MU_H*V2_SB)**2
      MSQ_SB2=(SUM_SB+DSQRT(DIF_SB2))/2.D0
      MSQ_SB1=(SUM_SB-DSQRT(DIF_SB2))/2.D0
*
      IF(MSQ_SB1.LT.0.D0) THEN ! Tachyonic sbottom_1
*       print*,'WARNING : SB1 becomes tachyonic! HB_R -> 0.9*HB_R',hb_r
       HB_R=0.9D0*HB_R
       GOTO 110
      ENDIF
*JSL:09/MAR/2006 Is it needed? I think so
      IF(MSQ_ST1.LT.0.D0) THEN ! Tachyonic stop_1
*       print*,'WARNING : ST1 becomes tachyonic! HT_R -> 0.9*HT_R',ht_r
       HT_R=0.9D0*HT_R
       GOTO 110
      ENDIF
*      print*,'>> Check 6 (NEW) : ',ITR,MSQ_ST1,MSQ_ST2,MSQ_SB1,
*     .                             MSQ_SB2,HT_R,HB_R
*      print*,'>> Check 6 (NEW) : ',ITR,HT_R,HB_R ! simplified check

*.....Threshold Corrections
      CDELHB1 = -2.D0*AS_SB/(3.D0*PI)*DCONJG(M3_H)*AB_H
     .                              *F_I(MSQ_SB1,MSQ_SB2,CDABS(M3_H)**2)
     .         -HT_R**2/(16.D0*PI**2)*CDABS(MU_H)**2
     .                              *F_I(MSQ_ST1,MSQ_ST2,CDABS(MU_H)**2)
      CDELHB2 =  2.D0*AS_SB/(3.D0*PI)*DCONJG(M3_H*MU_H)
     .                              *F_I(MSQ_SB1,MSQ_SB2,CDABS(M3_H)**2)
     .         +HT_R**2/(16.D0*PI**2)*DCONJG(AT_H*MU_H)
     .                              *F_I(MSQ_ST1,MSQ_ST2,CDABS(MU_H)**2)
      CDELHB = CDELHB1 + CDELHB2*V2_SB/V1_SB
      CDELHT1 =  2.D0*AS_ST/(3.D0*PI)*DCONJG(M3_H*MU_H)
     .                              *F_I(MSQ_ST1,MSQ_ST2,CDABS(M3_H)**2)
     .         +HB_R**2/(16.D0*PI**2)*DCONJG(AB_H*MU_H)
     .                              *F_I(MSQ_SB1,MSQ_SB2,CDABS(MU_H)**2)
      CDELHT2 = -2.D0*AS_ST/(3.D0*PI)*DCONJG(M3_H)*AT_H
     .                              *F_I(MSQ_ST1,MSQ_ST2,CDABS(M3_H)**2)
     .         -HB_R**2/(16.D0*PI**2)*CDABS(MU_H)**2
     .                              *F_I(MSQ_SB1,MSQ_SB2,CDABS(MU_H)**2)
      CDELHT = CDELHT1*V1_ST/V2_ST + CDELHT2
*      print*,'Iterating ...',itr,cdelhb1,cdelhb2
*      print*,'Iterating ...',itr,cdelht1,cdelht2

*
      HB_OLD=HB_R
      HT_OLD=HT_R
      HB_R=HB_NC/CDABS(1.D0+CDELHB)
      HT_R=HT_NC/CDABS(1.D0+CDELHT)

       IF(ITR.GE.ITRMAX) THEN
          IFLAG_H(54)=1
          RETURN
       ENDIF
       IF(DABS(HB_OLD-HB_R)/DABS(HB_OLD+HB_R).GT.EPS_TR .OR.
     .    DABS(HT_OLD-HT_R)/DABS(HT_OLD+HT_R).GT.EPS_TR ) THEN
          GOTO 111
       ENDIF
*
*      print*,'... Iteration ends :',itr,cdabs(cdelhb1),cdabs(cdelhb2)
*      print*,'... Iteration ends :',itr,cdabs(cdelht1),cdabs(cdelht2)
*      print*,'... Iteration ends : CDELHB',itr,cdelhb
*      print*,'... Iteration ends : CDELHT',itr,cdelht
*Killing threshold corrections BY HAND:
*      print*,'Killing threshold corrections BY HAND:',ht_nc,hb_nc
*      HB_R=HB_NC
*      HT_R=HT_NC
*.....Absoulte values of Yukawa couplings and its CP phases 
*.....including the threshold corrections
*JSL 09/Sep/2009, ht(mt^pole) and hb(mt^pole) improved
      HT_ST=HT_R
*      HT_MT=HT_ST*(1.D0+2.D0*BT*DLOG(MTPOLE_H**2/QST2))**0.25D0
       IF(MCH.GT.MTPOLE_H) THEN ! MTpole < MCH < MSfermion
        HT_MCH=HT_ST/(1.D0+2.D0*BT*DLOG(QST2/MCH**2))**0.25D0
        HT_MT =HT_MCH/(1.D0+2.D0*BTSM*DLOG(MCH**2/MTPOLE_H**2))**0.25D0
        HT_SF =HT_MCH*(1.D0+2.D0*BT*DLOG(QSF2/MCH**2))**0.25D0
       ELSE                     ! MCH < MTpole < M_Sfermion 
        HT_MT=HT_ST/(1.D0+2.D0*BT*DLOG(QST2/MTPOLE_H**2))**0.25D0
        HT_SF=HT_MT*(1.D0+2.D0*BT*DLOG(QSF2/MTPOLE_H**2))**0.25D0
       ENDIF                    ! IF(MCH.GT.MTPOLE_H) THEN
      HT_CP=CDABS(1.D0+CDELHT)/(1.D0+CDELHT)
*
      HB_SB=HB_R
*      HB_MT=HB_SB*(1.D0+2.D0*BB*DLOG(MTPOLE_H**2/QSB2))**0.25D0
       IF(MCH.GT.MTPOLE_H) THEN ! MTpole < MCH < MSfermion
        HB_MCH=HB_SB/(1.D0+2.D0*BB*DLOG(QSB2/MCH**2))**0.25D0
        HB_MT =HB_MCH/(1.D0+2.D0*BBSM*DLOG(MCH**2/MTPOLE_H**2))**0.25D0
        HB_SF =HB_MCH*(1.D0+2.D0*BB*DLOG(QSF2/MCH**2))**0.25D0
       ELSE                     ! MCH < MTpole < M_Sfermion 
        HB_MT=HB_SB/(1.D0+2.D0*BB*DLOG(QSB2/MTPOLE_H**2))**0.25D0
        HB_SF=HB_MT*(1.D0+2.D0*BB*DLOG(QSF2/MTPOLE_H**2))**0.25D0
       ENDIF                    ! IF(MCH.GT.MTPOLE_H) THEN
      HB_CP=CDABS(1.D0+CDELHB)/(1.D0+CDELHB)
*      print*,'>> FILLHIGGS ',hb_r,ht_r
*      print*,'>> FILLHIGGS ',hb_mt,ht_mt
*      print*,'>> FILLHIGGS ',bb,bt
*
      ENDIF      ! IF ( IFLAG_H(10)=0 )

*.....Check of perturbativity for h_t(M_SUSY) and h_b(M_SUSY)
      IF(HT_ST.GT.2.D0.OR.HB_SB.GT.2.D0) THEN
       IFLAG_H(55) = 1
       RETURN
      ENDIF
*
      RAUX_H(14)=V1
      RAUX_H(15)=V1_ST
      RAUX_H(16)=V1_SB
      RAUX_H(17)=V1_SF
*      print*,'V1: After corrections (no changes)',v1,v1_st,v1_sb,v1_sf
      RAUX_H(18)=V2
      RAUX_H(19)=V2_ST
      RAUX_H(20)=V2_SB
      RAUX_H(21)=V2_SF
*      print*,'V2: After corrections (no changes)',v2,v2_st,v2_sb,v2_sf
*
*      print*,'>> Check 7 (NEW) : ',HT_ST,HT_SF,HB_SB,HB_SF
*      print*,'>> FILLHIGGS ',cdelhb1,cdelhb2,cdelht1,cdelht2
*
      RAUX_H(24)=HT_MT
      RAUX_H(25)=HT_ST
      RAUX_H(26)=HT_SF
*      print*,'HT (After Corrections):',ht_mt,ht_st,ht_sf
      RAUX_H(27)=HB_MT
      RAUX_H(28)=HB_SB
      RAUX_H(29)=HB_SF
*      print*,'HB (After Corrections):',hb_mt,hb_sb,hb_sf
*
      CAUX_H(1)=HT_CP
      CAUX_H(2)=HB_CP
*      print*,'HT and HB Phases (After Corrections):',ht_cp,hb_cp
*-----------------------------------------------------------------------
*RG-improved Higgs-boson mass matrix in the Weak basis (phi_1,phi_2,a,G)
*at S=0
*-----------------------------------------------------------------------
*.....The squared off-shell momentum of the Higgs-boson propagator
*.....or the c.o.m. energy squared for Higgs production at colliders

*.....RG-improved MASQ: mass squared of the would-be CP-odd scalar
      CALL GET_MASQ(MCH**2,MCH,NFLAG,IFLAG_H                  
     .             ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .             ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .             ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF 
     .             ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .             ,HT_CP,HB_CP                               
     .             ,MASQ_P,REPI22_P,BL4VSQ_P,BARL4_P)
      MASQ  =MASQ_P
      RAUX_H(30)=MASQ
      RAUX_H(31)=REPI22_P
      RAUX_H(32)=BL4VSQ_P
      RAUX_H(33)=BARL4_P
*      print*,'MA^2',masq
*      print*,'Re(Pi22),L4/2*V^2,L4',REPI22_P,BL4VSQ_P,BARL4_P
*
*.....The complex 4X4 neutral Higgs mass matrix CMNH
      CALL GET_CMNH(0.D0,MCH,MASQ,NFLAG,IFLAG_H               
     .             ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .             ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .             ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF 
     .             ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .             ,HT_CP,HB_CP,BARL1,BARL2,BARL34                               
     .             ,CMNH)                                     
*      print*,'L1,L2,L34:',BARL1,BARL2,BARL34
      RAUX_H(34)=BARL1
      RAUX_H(35)=BARL2
      RAUX_H(36)=BARL34
*
      DO I=1,3
       DO J=1,3
        NH3(I,J)=DREAL(CMNH(I,J))
*        if (j.ge.i) print*,'MASS^2 matrix at s=0 ',i,j,nh3(i,j)
       ENDDO
      ENDDO
*
      CALL DIAGRS(3,3,NH3,EV3,AUX3,IERR_DIAGRS)
*
      Z12=NH3(1,1)*NH3(1,2)+NH3(2,1)*NH3(2,2)+NH3(3,1)*NH3(3,2)
      Z13=NH3(1,1)*NH3(1,3)+NH3(2,1)*NH3(2,3)+NH3(3,1)*NH3(3,3)
      Z23=NH3(1,2)*NH3(1,3)+NH3(2,2)*NH3(2,3)+NH3(3,2)*NH3(3,3)
      Z21=NH3(2,1)*NH3(1,1)+NH3(2,2)*NH3(1,2)+NH3(2,3)*NH3(1,3)
      Z31=NH3(3,1)*NH3(1,1)+NH3(3,2)*NH3(1,2)+NH3(3,3)*NH3(1,3)
      Z32=NH3(3,1)*NH3(2,1)+NH3(3,2)*NH3(2,2)+NH3(3,3)*NH3(2,3)
*      print*,'ZEROs? ',z12,z21
*      print*,'ZEROs? ',z13,z31
*      print*,'ZEROs? ',z23,z32
      IF( (DABS(Z12).GT.1.D-14) .OR. (DABS(Z21).GT.1.D-14) .OR.
     .    (DABS(Z13).GT.1.D-14) .OR. (DABS(Z31).GT.1.D-14) .OR.
     .    (DABS(Z23).GT.1.D-14) .OR. (DABS(Z32).GT.1.D-14) ) THEN
*        print*,'Check Orthognality of OMIX_0 at s=0'
        IFLAG_H(52)=1
        RETURN
      ENDIF
*
      IF(EV3(1).LE.0.D0.OR.EV3(2).LE.0.D0.OR.EV3(3).LE.0.D0) THEN
*         print*,'MHSQ(EFF) =',ev3(1),ev3(2),ev3(3)
        IFLAG_H(51)=1
        RETURN
      ENDIF
      IF(IERR_DIAGRS.GT.0) THEN
*        print*,'DIAGRS Error at sqrt(s) = 0'
        IFLAG_H(52)=1
        RETURN
      ENDIF
*The effective potential masses and the mixing matrix OMIX_0 at s=0
      EP3(1)=DSQRT(EV3(1))
      EP3(2)=DSQRT(EV3(2))
      EP3(3)=DSQRT(EV3(3))
      DO I=1,3
       DO J=1,3
        OMIX_0(I,J)=NH3(I,J)
*        print*,'OMIX_0',omix_0(i,j)
       ENDDO
      ENDDO
*      print*,'>> S=0 : Effective potential masses'
*      write(*,3) ep3(1),ep3(2),ep3(3)
*      print*,'>> Mixing matrix O at S=0'
*      write(*,3) oh3(1,1),oh3(1,2),oh3(1,3)
*      write(*,3) oh3(2,1),oh3(2,2),oh3(2,3)
*      write(*,3) oh3(3,1),oh3(3,2),oh3(3,3)
*
*-----------------------------------------------------------------------
*The Pole and On-Shell Masses : Iterative method is needed 
*Try On-Shell mass which is the same as the Pole Masses at one-loop level,
*by solving det{s^OS - Re[M^2(s=s^OS)]}=0 interatively where s^OS is real.
*-----------------------------------------------------------------------
      DO I=1,3
       DO J=1,3
        OMIX_1(I,J)=0.D0
        OMIX_2(I,J)=0.D0
        OMIX_3(I,J)=0.D0
       ENDDO
      ENDDO
*>>>>>
      IF(IFLAG_H(12).EQ.1 .OR.
     .   IFLAG_H(12).EQ.2 .OR. IFLAG_H(12).EQ.3) THEN
*>>>>>
*
      DO IH=1,3
       MPOLE=EP3(IH) ! The effective potential masses 
*                      and mixing matrix at zero momentum transfer are used
       CALL GET_CMNH(MPOLE**2,MCH,MASQ,NFLAG,IFLAG_H           
     .              ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .              ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .              ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF 
     .              ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .              ,HT_CP,HB_CP,BARL1,BARL2,BARL34                               
     .              ,CMNH)                                     
*      print*,'L1,L2,L34:',BARL1,BARL2,BARL34
*
       DO JH=1,3
       CTMP3(IH,JH)=DCMPLX(0.D0,0.D0)
         DO IA1=1,3
         DO IA2=1,3
         CTMP3(IH,JH)=CTMP3(IH,JH)
     .               +OMIX_0(IA1,IH)*CMNH(IA1,IA2)*OMIX_0(IA2,JH)
         ENDDO
         ENDDO
*         print*,'CTMP3(',IH,JH,') = ',ctmp3(ih,jh)
       ENDDO ! JH
*
       IF(DREAL(CTMP3(IH,IH)).LE.0.D0) THEN
*         print*,'MHSQ(OLD) =',ctmp3(1,1),ctmp3(2,2),ctmp3(3,3)
         IFLAG_H(51)=1
         RETURN
       ENDIF
*......M^pole
       MPOLE=DSQRT( DREAL(CTMP3(IH,IH)) )
       HP3(IH)=MPOLE
      ENDDO ! IH
*.....Improving Pole masses
      DH23=DABS( EP3(3)**2 - EP3(2)**2 )/10.D0
      IF(DH23.LT.CDABS(CTMP3(2,3)+CTMP3(3,2))) THEN
*
      HSQ23 = (EP3(2)**2 + EP3(3)**2 )/2.D0
      CALL GET_CMNH(HSQ23,MCH,MASQ,NFLAG,IFLAG_H              
     .             ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .             ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .             ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF 
     .             ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .             ,HT_CP,HB_CP,BARL1,BARL2,BARL34                               
     .             ,CMNH)                                     
*      print*,'L1,L2,L34:',BARL1,BARL2,BARL34
      DO IH=1,3
      DO JH=1,3
      CTMP3(IH,JH)=DCMPLX(0.D0,0.D0)
        DO IA1=1,3
        DO IA2=1,3
        CTMP3(IH,JH)=CTMP3(IH,JH)
     .              +OMIX_0(IA1,IH)*CMNH(IA1,IA2)*OMIX_0(IA2,JH)
        ENDDO
        ENDDO
      ENDDO ! JH
      ENDDO ! IH
      CTR23  = CTMP3(2,2)+CTMP3(3,3)
      CDET23 = CTMP3(2,2)*CTMP3(3,3) - CTMP3(2,3)*CTMP3(3,2)
      CD23   = CTR23**2 - 4.D0*CDET23
      HP3(2)  = DSQRT( 0.5D0*DREAL( CTR23 - CDSQRT(CD23) ) )
      HP3(3)  = DSQRT( 0.5D0*DREAL( CTR23 + CDSQRT(CD23) ) )
      IFLAG_H(53)=1
*
*      print*,'>> Improved Pole masses :'
*      print*,hp3(1)
*      print*,hp3(2)
*      print*,hp3(3)
*
      ENDIF ! Improving Pole masses
*
*>>>>>
      ELSEIF (IFLAG_H(12).EQ.4 .OR. IFLAG_H(12).EQ.5) THEN
*>>>>>
      ITRHP_MAX=200
      DO IH=1,3
       MPOLE=EP3(IH) ! Initially, the effective potential masses 
*                      and mixing matrix at zero momentum transfer are used
*       print*,'%-------- MPOLE_START = ',MPOLE
       DO ITRHP=1,ITRHP_MAX
       CALL GET_CMNH(MPOLE**2,MCH,MASQ,NFLAG,IFLAG_H           
     .              ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .              ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .              ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF 
     .              ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .              ,HT_CP,HB_CP,BARL1,BARL2,BARL34                              
     .              ,CMNH)                                     
*      print*,'L1,L2,L34:',BARL1,BARL2,BARL34
       DO I=1,3
        DO J=1,3
         NH3(I,J)=DREAL(CMNH(I,J))
         DMH3(I,J)=-DREAL(CMNH(I,J))
         IF(I.EQ.J) DMH3(I,J)=MPOLE**2+DMH3(I,J) ! s^OS - Re[M^2(s=s^OS)]}
        ENDDO
       ENDDO

       CALL DIAGRS(3,3,NH3,EV3,AUX3,IERR_DIAGRS)
       IF(EV3(IH).LE.0.D0) THEN
*         print*,'MHSQ(OS) =',ev3(1),ev3(2),ev3(3)
         IFLAG_H(51)=1
         RETURN
       ENDIF
       IF(IERR_DIAGRS.GT.0) THEN
*         print*,'DIAGRS Error at sqrt(s) = ', MPOLE
         IFLAG_H(52)=1
         RETURN
       ENDIF

       DET_DMH3=DMH3(1,1)*DMH3(2,2)*DMH3(3,3)
     .         -DMH3(1,1)*DMH3(2,3)*DMH3(3,2)
     .         -DMH3(2,1)*DMH3(1,2)*DMH3(3,3)
     .         +DMH3(2,1)*DMH3(1,3)*DMH3(3,2)
     .         +DMH3(3,1)*DMH3(1,2)*DMH3(2,3)
     .         -DMH3(3,1)*DMH3(1,3)*DMH3(2,2)
       DMH_IH=MPOLE-DSQRT(EV3(IH))

*Monitoring Iterations...
*       print*,'IH = ',ih,' ITR = ',itrhp
*     .       ,mpole,dsqrt(ev3(ih)),det_dmh3/mch**6,dmh_ih

*Mpole and Omix(M^pole)
       if(itrhp.le.itrhp_max-150) then
       MPOLE=DSQRT(EV3(IH))
       else  ! to treat some oscillatory case
       mpole=(mpole+DSQRT(EV3(IH)))/2.d0
       endif
*
       DO I=1,3
        DO J=1,3
         IF(IH.EQ.1) OMIX_1(I,J)=NH3(I,J)
         IF(IH.EQ.2) OMIX_2(I,J)=NH3(I,J)
         IF(IH.EQ.3) OMIX_3(I,J)=NH3(I,J)
        ENDDO
       ENDDO
*Exit iteration...
       IF(DABS(DET_DMH3/MCH**6).LT.1.D-6 .AND. 
     .    DABS(DMH_IH).LT.1.D-3) GOTO 88
       IF(ITRHP.EQ.ITRHP_MAX) THEN
         IFLAG_H(60)=1
*         print*,'ITERATION for the on-shell Higgs masses FAILS !'
         RETURN
       ENDIF
*
       ENDDO ! ITRHP 
 88    HP3(IH)=MPOLE
      ENDDO ! IH
*
*>>>>>
      ELSE
*       print*,'Invalid IFLAG_H(12)'
      RETURN
*>>>>>
      ENDIF !  IFLAG_H(12)
*>>>>>
*-----------------------------------------------------------------------
*Effective-Potential Mass or Pole mass??
*-----------------------------------------------------------------------
*>>>>>
      IF(IFLAG_H(11).EQ.1) THEN 
*>>>>>
        HMASS(1)=EP3(1)
        HMASS(2)=EP3(2)
        HMASS(3)=EP3(3)
        CALL GET_MASQ(0.D0,MCH,NFLAG,IFLAG_H                    
     .               ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2        
     .               ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF 
     .               ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF
     .               ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF 
     .               ,HT_CP,HB_CP                               
     .               ,MASQ_0,REPI22_0,BL4VSQ_0,BARL4_0)
*-The inverse of propagator    s-M_0^2+Pi^hat(s) with the tree-leve mass M_0
*-The effective potential mass M_eff^2  = M_0^2-Real[Pi^hat(0)]
*-The pole mass                M_pole^2 = M_0^2-Real[Pi^hat(M_pole^2)]
*-Threfore we have             M_eff^2  = M_pole^2+Real[Pi^hat(M_pole^2)-Pi-hat(0)]
        MCHSQ_EP=MCH**2+REPI22_P-REPI22_0 
        MCH_OUT =DSQRT(MCHSQ_EP)
        DO IH=1,3
         DO JH=1,3
          OMIX_OUT(IH,JH)=OMIX_0(IH,JH)
         ENDDO
        ENDDO
*>>>>>
      ELSEIF(IFLAG_H(11).EQ.0) THEN
*Pole Masses for the neutral and charged Higgs bosons
*>>>>>
*.....Reordering of HP(3) and OMIX_0 ..............................
        IF(HP3(1)-HP3(2).GT.0.D0 .OR.
     .     HP3(2)-HP3(3).GT.0.D0 .OR.
     .     HP3(1)-HP3(3).GT.0.D0) THEN
*
        IH_MIN=1
         DO IH=1,3
           IF(HP3(IH).LT.HP3(IH_MIN)) IH_MIN=IH
         ENDDO
*
        IH_MAX=3
         DO IH=1,3
           IF(HP3(IH).GT.HP3(IH_MAX)) IH_MAX=IH
         ENDDO
*
        IH_MID=6-IH_MIN-IH_MAX
*        print*,'>>>> NEW Reordering...',HP3(1),HP3(2),HP3(3),' --> '
*     .        ,IH_MIN,IH_MID,IH_MAX
*
        DO IH=1,3
         HP3_TMP(IH)=HP3(IH)
         DO IA=1,3
          OMIX_TMP(IA,IH)=OMIX_0(IA,IH)
         ENDDO
        ENDDO
*
        HP3(1)=HP3_TMP(IH_MIN)
        HP3(2)=HP3_TMP(IH_MID)
        HP3(3)=HP3_TMP(IH_MAX)
        DO IA=1,3
          OMIX_0(IA,1)=OMIX_TMP(IA,IH_MIN)
          OMIX_0(IA,2)=OMIX_TMP(IA,IH_MID)
          OMIX_0(IA,3)=OMIX_TMP(IA,IH_MAX)
        ENDDO
*
        ENDIF 
*.....End of Reordering.............................................
        HMASS(1)=HP3(1)
        HMASS(2)=HP3(2)
        HMASS(3)=HP3(3)
        MCH_OUT =MCH
        DO IH=1,3
         DO JH=1,3
          OMIX_OUT(IH,JH)=OMIX_0(IH,JH)
         ENDDO
        ENDDO
      ELSE
        WRITE(6,*) 'INVALID OPTION OF IFLAG(11)'
        RETURN
      ENDIF
*Fill some HC_AUX COMMON Blocks with MCH_OUT 
      RAUX_H(10)=MCH_OUT
*      print*,'RAUX_H[10]',mch_out
*-----------------------------------------------------------------------
*Print results
*-----------------------------------------------------------------------
*      IF(IFLAG_H(2).EQ.1) 
*     .            CALL DUMP_HIGGS(NFLAG,IFLAG_H,MCH_OUT,HMASS,OMIX_OUT)
************************************************************************
      RETURN
      END

      SUBROUTINE GET_MASQ(SPRO,MCH,NFLAG,IFLAG_H
     .                   ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2
     .                   ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF
     .                   ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF
     .                   ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF
     .                   ,HT_CP,HB_CP 
     .                   ,MASQ,REPI22,BL4VSQ,BARL4)
************************************************************************
*
* RG-improved MASQ: mass squared of the would-be CP-odd scalar at m_t-pole
*
* M_A^2 = MCH^pole - lambda4/2 v^2 + Re [Pi-hat(s=MCH^pole^2,mtpole)]
*                                          [Eq.(3.8) of NPB625(2002)345]
*
* MASQ = MCH^2 - BARL4/2*V^2 + Re[CPICH(H^+,H^-)]
*      = MCH^2 - BARL4/2*V^2 - Re[CMCH_1L(H^+,H^-)]
* where
* CMCH_1L(H^+,H^-)=CMCH_1L(1,1)*SB^2-[CMCH_1L(1,2)+CMCH_1L(2,1)]*CB*SB
*                 +CMCH_1L(2,2)*CB^2 at mtpole and 
* CMCH_1L(I,J) = -CM0_1L(I,J)-CPI_SQ(I,J)-CPI_Q(I,J) is a 2X2 charged 
* Higgs-boson mass matrix in (phi_1^pm,phi_2^pm0) basis dropping the two-loop 
* Born-improved term at the scale m_t-pole
*
* CM0_1L : One-loop part of tree-level-form mass matrix at the scale Q_tb
*          after multiplying anomalous dim. factors XI_I and XI_J
*          given by Eq. (2.9)
* CPI_SQ : Squark contributions to the self energy at Q_tb
*          after multiplying anomalous dim. factors XI_I and XI_J
* CPI_Q  : Quark contributions to the self energy at m_t
*          
* For explict forms, see (B.12), (B.13), (B.15), and (B.16)(d) of 
* NPB625(2002)345
*
*---> SPRO,MCH                             ! S and Charged Higgs-boson pole mass
*---> QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2   ! sfermion scales 
*---> XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF ! anomalous dimensions
*       : ST=Stop mass scale
*       : SB=Sbottom mass scale
*       : SF=Sfermion mass scale
*---> V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF ! vevs
*       : V1, V2 at Mtpole
*---> HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF ! Yukawa couplings
*       : HT   , HB    at Mtpole without threshold corrections
*       : HT_MT, HB_MT at Mtpole with    threshold corrections
*---> HT_CP,HB_CP                               ! CP phases of Yukawa couplings
*       : HT_X(complex)=HT_X*HT_CP with X=none, MT, ST, SF
*       : HB_X(complex)=HB_X*HB_CP with X=none, MT, SB, SF
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 IFLAG_H(NFLAG)
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
*-----------------------------------------------------------------------
*Local
      COMPLEX*16 HT_CP,HB_CP
      COMPLEX*16 DELTA(2,2),UT3(2,2),UB3(2,2)
      COMPLEX*16 CF1TT(2,2),CF2TT(2,2),CF1BB(2,2),CF2BB(2,2),CA1TT(2,2),
     .           CA2TT(2,2),CA1BB(2,2),CA2BB(2,2),CP1TB(2,2),CP2TB(2,2),
     .           CPP1TT(2,2),CPP2TT(2,2),CPP1BB(2,2),CPP2BB(2,2)
      COMPLEX*16 CHTR,CHBR,CXI
      COMPLEX*16 C0U,C1UB,C1UT,C2U
      COMPLEX*16 CTB11A,CTB11B,CTB11C,CTAD11
      COMPLEX*16 CTB12A,CTB12B,CTB12C,CTAD12
      COMPLEX*16 CTB21A,CTB21B,CTB21C,CTAD21
      COMPLEX*16 CTB22A,CTB22B,CTB22C,CTAD22
      COMPLEX*16 CB0_H
      COMPLEX*16 CM0_1L(2,2),CPI_SQ(2,2),CPI_Q(2,2),CMCH_1L(2,2)
      COMPLEX*16 CMCH22,CPICH22
*
      PI=2.D0*DASIN(1.D0)
      CXI=DCMPLX(0.D0,1.D0)
*
************************************************************************
      S=SPRO
*at top-quark pole mass scale
      GS2 =4.D0*PI*ASMT_H
*-------------------------------------------------------------------------
*BARL4 at mtpole
*     Eq.(3.6) NPB586(2000)92
      IF(IFLAG_H(12).EQ.3 .OR. IFLAG_H(12).EQ.5) THEN
       HBX=HB_MT
       HTX=HT_MT
      ELSE
       HBX=HB
       HTX=HT
      ENDIF
*
      X1L4=3.D0/(16.D0*PI**2)*( HTX**2*HBX**2*(DLOG(QQT2/MTPOLE_H**2)
     .    +DLOG(DMAX1(QTT2,QBB2)/MTPOLE_H**2))
     .    -GW_H**2/2.D0*(HTX**2-GW_H**2/4.D0)
     .     *DLOG(QQT2/MTPOLE_H**2)-GW_H**2/2.D0
     .     *(HBX**2-GW_H**2/4.D0)*DLOG(QQB2/MTPOLE_H**2) )

*Threshold corrections have NOT included in two-loop couplings
*     Eqs.(3.11) and (3.12) NPB586(2000)92
      X2L4=3.D0*HB**2*HT**2/(16.D0*PI**2)**2*(HT**2+HB**2-8.D0*GS2)
     .*( (DLOG(QQT2/MTPOLE_H**2))**2 
     .  +(DLOG(DMAX1(QTT2,QBB2)/MTPOLE_H**2))**2 )
*
*.....Computation of the chargino and neutrlino contributions to the
*.....couplings \lambda_1, \lambda_2, \lambda_{34} = \lambda_3+\lambda_4,
*.....and \lambda_4: XINO1, XINO2, XINO34, XINO4
*.....In this determination the formulas in Appendix C of H.E. Haber and
*.....R. Hempfling, Phys. Rev. D48 (1993) 4280, are used.
*
      TWINO=0.D0
      THINO=0.D0
      TCHI1=0.D0
      TCHI2=0.D0
      TCHI12=0.D0
*
      IF( QSF2.GT.CDABS(M2_H)**2 ) THEN
       TWINO=-DLOG( QSF2/CDABS(M2_H)**2 )
       IF(MTPOLE_H.GT.CDABS(M2_H)) TWINO=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.CDABS(MU_H)**2 ) THEN
       THINO=-DLOG( QSF2/CDABS(MU_H)**2 )
       IF(MTPOLE_H.GT.CDABS(MU_H)) THINO=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.DMAX1(CDABS(MU_H)**2,CDABS(M1_H)**2) ) THEN
       TCHI1=-DLOG( QSF2/DMAX1(CDABS(MU_H)**2,CDABS(M1_H)**2) )
       IF(MTPOLE_H.GT.DMAX1(CDABS(MU_H),CDABS(M1_H)))
     .  TCHI1=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.DMAX1(CDABS(MU_H)**2,CDABS(M2_H)**2) ) THEN
       TCHI2=-DLOG( QSF2/ DMAX1(CDABS(MU_H)**2,CDABS(M2_H)**2) )
       IF(MTPOLE_H.GT.DMAX1(CDABS(MU_H),CDABS(M1_H)))
     .  TCHI2=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( DSQRT(QSF2).GT.
     .    DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) ) THEN
       TCHI12=-2.D0*DLOG(DSQRT(QSF2)
     .     /DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) )
       IF( MTPOLE_H.GT.
     .     DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) )
     . TCHI12=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      XINO4=GW_H**4/(192.D0*PI**2)*(
     . 6.D0*SW_H**2/CW_H**2*TCHI1 + 24.D0*SW_H**2/CW_H**2*TCHI12
     .-6.D0*TCHI2 - 4.D0*THINO - 8.D0*TWINO )
*
      BARL4 = GW_H**2/2.D0 + X1L4 + X2L4 + XINO4
*
*-------------------------------------------------------------------------
*CM0_1L at Q_tb and including running from Q_tb to mtpole
*
*     Eq.(3.6) NPB586(2000)92
      X1L4Q=3.D0/(16.D0*PI**2)*( HT_SF**2*HB_SF**2*(DLOG(QQT2/QSF2)
     .     +DLOG(DMAX1(QTT2,QBB2)/QSF2))
     .     -GW_H**2/2.D0*(HT_SF**2-GW_H**2/4.D0)
     .      *DLOG(QQT2/QSF2)-GW_H**2/2.D0
     .      *(HB_SF**2-GW_H**2/4.D0)*DLOG(QQB2/QSF2) )
*     Eq.(3.7) NPB586(2000)92
      XRM12Q=3.D0/(16.D0*PI**2)
     .*(HT_SF**2*DREAL(MU_H*AT_H)*DLOG(QST2/QSF2) 
     . +HB_SF**2*DREAL(MU_H*AB_H)*DLOG(QSB2/QSF2) )
*     Eq.(2.8) NPB586(2000)92
      TANB_SF=V2_SF/V1_SF
      FAC_CM0=X1L4Q*V1_SF*V2_SF/2.D0+XRM12Q
      CM0_1L(1,1)=DCMPLX(FAC_CM0*TANB_SF,0.D0)/XI1_SF/XI1_SF
      CM0_1L(1,2)=DCMPLX(-FAC_CM0       ,0.D0)/XI1_SF/XI2_SF
      CM0_1L(2,1)=DCMPLX(-FAC_CM0       ,0.D0)/XI2_SF/XI1_SF
      CM0_1L(2,2)=DCMPLX(FAC_CM0/TANB_SF,0.D0)/XI2_SF/XI2_SF
*      print*,DCMPLX(FAC_CM0*TANB_SF,0.D0),XI1_SF,XI1_SF
*      print*,DCMPLX(-FAC_CM0       ,0.D0),XI1_SF,XI2_SF
*      print*,DCMPLX(-FAC_CM0       ,0.D0),XI2_SF,XI1_SF
*      print*,DCMPLX(FAC_CM0/TANB_SF,0.D0),XI2_SF,XI2_SF
*      print*,'>> GET_MASQ 2 <<',XRM12Q
*-------------------------------------------------------------------------
*CPI_Q at mtpole One-loop t- and b- contributions to the charged Higgs-boson
*mass matrix in the weak basis (\phi+_1,\phi+_2):
*  [Eq.(B.15) + diag(T_\phi_i(e)/v_i) of Eq.(B.16)] in NPB625(2002)345 
*
*The relation A0(m**2)=m**2[1+B0_H(0,m**2,m**2)] has been used.
*
*For the quark masses inside the loop, we use top-quark pole mass for 
*the top-quark case and m_b(m_t) for the b-quark case are used. 
*The mixed uses of the t'Hooft scale [m_t-pole or m_t(m_t^pole)] will be 
*clarified later <-- ERROR ?? 
*
*See also CPI_Q in GET_CMNH
*
      IF(IFLAG_H(12).EQ.2 .OR. IFLAG_H(12).EQ.5) THEN
       HBX=HB_MT
       HTX=HT_MT
      ELSE
       HBX=HB
       HTX=HT
      ENDIF
*
      CPI_Q(1,1) =-3.D0/(16.D0*PI**2)*HBX**2
     .*( MTMT_H**2
     .  +MTMT_H**2*B0_H(0.D0,MTPOLE_H**2,MTPOLE_H**2,MTPOLE_H**2)
     .  -MBMT_H**2
     .  -MBMT_H**2*B0_H(0.D0,MBMT_H**2,MBMT_H**2,MTPOLE_H**2)
     .  -(S-MTMT_H**2-MBMT_H**2)
     .                      *CB0_H(S,MBMT_H**2,MTPOLE_H**2,MTPOLE_H**2))
      CPI_Q(1,2) = 3.D0/(8.D0*PI**2)*HBX*HTX*MTMT_H*MBMT_H
     .            *CB0_H(S,MBMT_H**2,MTMT_H**2,MTMT_H**2)
      CPI_Q(2,1) = CPI_Q(1,2)
      CPI_Q(2,2) =-3.D0/(16.D0*PI**2)*HTX**2
     .*(-MTMT_H**2
     .  -MTMT_H**2*B0_H(0.D0,MTPOLE_H**2,MTPOLE_H**2,MTPOLE_H**2)
     .  +MBMT_H**2
     .  +MBMT_H**2*B0_H(0.D0,MBMT_H**2,MBMT_H**2,MTPOLE_H**2)
     .  -(S-MTMT_H**2-MBMT_H**2)
     .                       *CB0_H(S,MBMT_H**2,MTMT_H**2,MTMT_H**2))
*-------------------------------------------------------------------------
*CPI_SQ at Q_tb and running from Q_tb to mtpole
*.....One-loop ~t-, ~b- contributions to the charged Higgs-boson
*.....mass matrix in the weak basis (\phi+_1,\phi+_2):

*SQURK masses squared and the mixing matrices UT3 and UB3 at SF scale:
*     Eq.(B.10) NPB625(2002)345 
*UT3 = Ut \tau_3 (Ut)dagger ;    UB3 = Ub \tau_3 (Ub)dagger
      GWY2=GW_H**2/8.D0/CW_H**2
      GXT2=(GW_H**2-5.D0*GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0
      GXB2=(GW_H**2-GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0

      SUM_ST=MQ3_H**2+MU3_H**2+HT_SF**2*V2_SF**2
     .      +GWY2*(V1_SF**2-V2_SF**2)
      DIF_ST2=(MQ3_H**2-MU3_H**2
     .            +GXT2*(V1_SF**2-V2_SF**2)/2.D0)**2
     .       +2.D0*HT_SF**2*CDABS(DCONJG(AT_H)*V2_SF-MU_H*V1_SF)**2
      MSQ_ST2=(SUM_ST+DSQRT(DIF_ST2))/2.D0
      MSQ_ST1=(SUM_ST-DSQRT(DIF_ST2))/2.D0
      SUM_SB=MQ3_H**2+MD3_H**2+HB_SF**2*V1_SF**2
     .      +GWY2*(V2_SF**2-V1_SF**2)
      DIF_SB2=(MQ3_H**2-MD3_H**2
     .            +GXB2*(V2_SF**2-V1_SF**2)/2.D0)**2
     .       +2.D0*HB_SF**2*CDABS(DCONJG(AB_H)*V1_SF-MU_H*V2_SF)**2
      MSQ_SB2=(SUM_SB+DSQRT(DIF_SB2))/2.D0
      MSQ_SB1=(SUM_SB-DSQRT(DIF_SB2))/2.D0

      IF(MSQ_ST1.LT.0.D0.OR.MSQ_ST2.LT.0.D0.OR.
     .   MSQ_SB1.LT.0.D0.OR.MSQ_SB2.LT.0.D0) THEN
       IFLAG_H(50)=1
       RETURN
      ENDIF

      UT3(1,1)=(MQ3_H**2-MU3_H**2
     .          +GXT2*(V1_SF**2-V2_SF**2)/2.D0)/(MSQ_ST2-MSQ_ST1)
      UT3(2,1)=DSQRT(2.D0)*HT_SF*HT_CP*(AT_H*V2_SF-DCONJG(MU_H)*V1_SF)
     .                          /(MSQ_ST2-MSQ_ST1)
      UT3(1,2)=DCONJG( UT3(2,1) )
      UT3(2,2)=-UT3(1,1)

      UB3(1,1)=(MQ3_H**2-MD3_H**2
     .           -GXB2*(V1_SF**2-V2_SF**2)/2.D0)/(MSQ_SB2-MSQ_SB1)
      UB3(2,1)=DSQRT(2.D0)*HB_SF*HB_CP*(AB_H*V1_SF-DCONJG(MU_H)*V2_SF)
     .                             /(MSQ_SB2-MSQ_SB1)
      UB3(1,2)=DCONJG( UB3(2,1) )
      UB3(2,2)=-UB3(1,1)
*
*      print*,'>> GET_MASQ 4 <<',MSQ_ST2,MSQ_ST1,MSQ_SB2,MSQ_SB1
*      print*,'>> GET_MASQ 5 <<',UT3(1,1),UT3(1,2),UT3(2,2)
*      print*,'>> GET_MASQ 6 <<',UB3(1,1),UB3(1,2),UB3(2,2)

*Here charged Higgs-SQURK-SQUARK couplings : 
      CHTR = HT_SF*HT_CP
      CHBR = HB_SF*HB_CP
      V1RT = V1_SF
      V2RT = V2_SF
      V1RB = V1_SF
      V2RB = V2_SF
* The complex Yukawa couplings with threhsold corrections should be used ?:
*=>One can show that the phases of the resummed h_t and h_b in
*=>the squark sector can be entirely absorbed to the right-handed stops and
*=>sbottoms. This property applies to (A.1)[or (B.10) : See UT3(2,1) and
*=>UB3(2,1)], (A.5) and (A.7). Notice that in (A.1) the quark masses are 
*=>not the physical masses, but the resummed h_tand h_b \times v_2/sqrt(2) 
*=>and v_1/sqrt(2), respectively.
*
*=>We have checked this rephasing invariance by comparing the old version of
*=>CPsuperH in which only the absolute values of the Yukawa couplings are used
*=>and the current version considering the full complex Yukawa couplings.
*--->Start
* Gamma{\phi_1 ~t* ~t}:       CF1TT(2,2)                              *
      CF1TT(1,1)=-V1RT*(GW_H**2-GP_H**2/3.D0)/4.D0
      CF1TT(1,2)=DCONJG(CHTR)*MU_H/DSQRT(2.D0)
      CF1TT(2,1)=DCONJG( CF1TT(1,2) )
      CF1TT(2,2)=-V1RT*GP_H**2/3.D0
* Gamma{\phi_2 ~t* ~t}:       CF2TT(2,2)                              *
      CF2TT(1,1)=-CDABS(CHTR)**2*V2RT+V2RT*(GW_H**2-GP_H**2/3.D0)/4.D0
      CF2TT(1,2)=-DCONJG(CHTR*AT_H)/DSQRT(2.D0)
      CF2TT(2,1)=DCONJG( CF2TT(1,2) )
      CF2TT(2,2)=-CDABS(CHTR)**2*V2RT+V2RT*GP_H**2/3.D0
* Gamma{\phi_1 ~b* ~b}:       CF1BB(2,2)                              *
      CF1BB(1,1)=-CDABS(CHBR)**2*V1RB+V1RB*(GW_H**2+GP_H**2/3.D0)/4.D0
      CF1BB(1,2)=-DCONJG(CHBR*AB_H)/DSQRT(2.D0)
      CF1BB(2,1)=DCONJG( CF1BB(1,2) )
      CF1BB(2,2)=-CDABS(CHBR)**2*V1RB+V1RB*GP_H**2/6.D0
* Gamma{\phi_2 ~b* ~b}:       CF2BB(2,2)                              *
      CF2BB(1,1)=-V2RB*(GW_H**2+GP_H**2/3.D0)/4.D0
      CF2BB(1,2)=DCONJG(CHBR)*MU_H/DSQRT(2.D0)
      CF2BB(2,1)=DCONJG( CF2BB(1,2) )
      CF2BB(2,2)=-V2RB*GP_H**2/6.D0
* Gamma{a_1 ~t* ~t}:          CA1TT(2,2)                              *
      CA1TT(1,1)=0.D0
      CA1TT(1,2)=-CXI*DCONJG(CHTR)/DSQRT(2.D0) * MU_H
      CA1TT(2,1)=DCONJG( CA1TT(1,2) )
      CA1TT(2,2)=0.D0
* Gamma{a_2 ~t* ~t}:          CA2TT(2,2)                              *
      CA2TT(1,1)=0.D0
      CA2TT(1,2)=CXI*DCONJG(CHTR)/DSQRT(2.D0) * DCONJG(AT_H)
      CA2TT(2,1)=DCONJG( CA2TT(1,2) )
      CA2TT(2,2)=0.D0
* Gamma{a_1 ~b* ~b}:          CA1BB(2,2)                              *
      CA1BB(1,1)=0.D0
      CA1BB(1,2)=-CXI*DCONJG(CHBR)/DSQRT(2.D0) * DCONJG(AB_H)
      CA1BB(2,1)=DCONJG( CA1BB(1,2) )
      CA1BB(2,2)=0.D0
* Gamma{a_2 ~b* ~b}:          CA2BB(2,2)                              *
      CA2BB(1,1)=0.D0
      CA2BB(1,2)=CXI*DCONJG(CHBR)/DSQRT(2.D0) * MU_H
      CA2BB(2,1)=DCONJG( CA2BB(1,2) )
      CA2BB(2,2)=0.D0
* Gamma{\phi+_1 ~t* ~b}:     CP1TB(2,2)                               *
      CP1TB(1,1)=-V1_SF*(HB_SF**2-GW_H**2/2.D0)/DSQRT(2.D0)
      CP1TB(1,2)=-DCONJG(HB_SF*HB_CP*AB_H)
      CP1TB(2,1)=-HT_SF*HT_CP*DCONJG(MU_H)
      CP1TB(2,2)=-HT_SF*HT_CP*DCONJG(HB_SF*HB_CP)*V2_SF/DSQRT(2.D0)
* Gamma{\phi+_2 ~t* ~b}:     CP2TB(2,2)                               *
      CP2TB(1,1)=V2_SF*(HT_SF**2-GW_H**2/2.D0)/DSQRT(2.D0)
      CP2TB(1,2)=DCONJG(HB_SF*HB_CP)*MU_H
      CP2TB(2,1)=HT_SF*HT_CP*AT_H
      CP2TB(2,2)=HT_SF*HT_CP*DCONJG(HB_SF*HB_CP)*V1_SF/DSQRT(2.D0)
* Gamma{\phi+_1 \phi+_1 ~t* ~t}:  CPP1TT(2,2)                         *
      CPP1TT(1,1)=-HB_SF**2+(GW_H**2+GP_H**2/3.D0)/4.D0
      CPP1TT(1,2)=0.D0
      CPP1TT(2,1)=0.D0
      CPP1TT(2,2)=-GP_H**2/3.D0
* Gamma{\phi+_2 \phi+_2 ~t* ~t}:  CPP2TT(2,2)                         *
      CPP2TT(1,1)=-(GW_H**2+GP_H**2/3.D0)/4.D0
      CPP2TT(1,2)=0.D0
      CPP2TT(2,1)=0.D0
      CPP2TT(2,2)=-HT_SF**2+GP_H**2/3.D0
* Gamma{\phi+_1 \phi+_1 ~b* ~b}:  CPP1BB(2,2)                         *
      CPP1BB(1,1)=-(GW_H**2-GP_H**2/3.D0)/4.D0
      CPP1BB(1,2)=0.D0
      CPP1BB(2,1)=0.D0
      CPP1BB(2,2)=-HB_SF**2+GP_H**2/6.D0
* Gamma{\phi+_2 \phi+_2 ~b* ~b}:  CPP2BB(2,2)                         *
      CPP2BB(1,1)=-HT_SF**2+(GW_H**2-GP_H**2/3.D0)/4.D0
      CPP2BB(1,2)=0.D0
      CPP2BB(2,1)=0.D0
      CPP2BB(2,2)=-GP_H**2/6.D0
*-->EOF

*Some constatnt for CTB..(A,B,C)
      C0U  =CB0_H(S,MSQ_ST1,MSQ_SB1,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB2,QSF2)
     .     +CB0_H(S,MSQ_ST1,MSQ_SB2,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB1,QSF2)
      C1UB=-CB0_H(S,MSQ_ST1,MSQ_SB1,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB2,QSF2)
     .     +CB0_H(S,MSQ_ST1,MSQ_SB2,QSF2)-CB0_H(S,MSQ_ST2,MSQ_SB1,QSF2)
      C1UT=-CB0_H(S,MSQ_ST1,MSQ_SB1,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB2,QSF2)
     .     -CB0_H(S,MSQ_ST1,MSQ_SB2,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB1,QSF2)
      C2U  =CB0_H(S,MSQ_ST1,MSQ_SB1,QSF2)+CB0_H(S,MSQ_ST2,MSQ_SB2,QSF2)
     .     -CB0_H(S,MSQ_ST1,MSQ_SB2,QSF2)-CB0_H(S,MSQ_ST2,MSQ_SB1,QSF2)
*      print*,'B functions : ',c0u,c1ub,c1ut,c2u

*CTB..A : The 1st part of Eq.(B.12) NPB625(2002)345
      CTB11A=DCMPLX(0.D0,0D0)
      CTB12A=DCMPLX(0.D0,0D0)
      CTB21A=DCMPLX(0.D0,0D0)
      CTB22A=DCMPLX(0.D0,0D0)
      DO I=1,2
       DO J=1,2
        CTB11A= 3.D0/(64.D0*PI**2)*C0U
     .        *CP1TB(I,J)*DCONJG(CP1TB(I,J))+CTB11A
        CTB12A= 3.D0/(64.D0*PI**2)*C0U
     .        *CP1TB(I,J)*DCONJG(CP2TB(I,J))+CTB12A
        CTB21A= 3.D0/(64.D0*PI**2)*C0U
     .        *CP2TB(I,J)*DCONJG(CP1TB(I,J))+CTB21A
        CTB22A= 3.D0/(64.D0*PI**2)*C0U
     .        *CP2TB(I,J)*DCONJG(CP2TB(I,J))+CTB22A
       ENDDO
      ENDDO
*      print*,'>> GET_MASQ 7 <<',CTB11A
*      print*,'>> GET_MASQ 7 <<',CTB12A
*      print*,'>> GET_MASQ 7 <<',CTB21A
*      print*,'>> GET_MASQ 7 <<',CTB22A
*CTB..B : The 2nd and 3rd parts of Eq.(B.12) NPB625(2002)345
      CTB11B=DCMPLX(0.D0,0D0)
      CTB12B=DCMPLX(0.D0,0D0)
      CTB21B=DCMPLX(0.D0,0D0)
      CTB22B=DCMPLX(0.D0,0D0)
      DO I=1,2
       DO J=1,2
        DO K=1,2
         CTB11B= 3.D0/(64.D0*PI**2)
     .         *(C1UB*CP1TB(I,J)*UB3(J,K)*DCONJG(CP1TB(I,K))
     .          +C1UT*UT3(I,J)*CP1TB(J,K)*DCONJG(CP1TB(I,K)))+CTB11B
         CTB12B= 3.D0/(64.D0*PI**2)
     .         *(C1UB*CP1TB(I,J)*UB3(J,K)*DCONJG(CP2TB(I,K))
     .          +C1UT*UT3(I,J)*CP1TB(J,K)*DCONJG(CP2TB(I,K)))+CTB12B
         CTB21B= 3.D0/(64.D0*PI**2)
     .         *(C1UB*CP2TB(I,J)*UB3(J,K)*DCONJG(CP1TB(I,K))
     .          +C1UT*UT3(I,J)*CP2TB(J,K)*DCONJG(CP1TB(I,K)))+CTB21B
         CTB22B= 3.D0/(64.D0*PI**2)
     .         *(C1UB*CP2TB(I,J)*UB3(J,K)*DCONJG(CP2TB(I,K))
     .          +C1UT*UT3(I,J)*CP2TB(J,K)*DCONJG(CP2TB(I,K)))+CTB22B
        ENDDO
       ENDDO
      ENDDO
*      print*,'>> GET_MASQ 8 <<',CTB11B
*      print*,'>> GET_MASQ 8 <<',CTB12B
*      print*,'>> GET_MASQ 8 <<',CTB21B
*      print*,'>> GET_MASQ 8 <<',CTB22B
*CTB..C : The 4th part of Eq.(B.12) NPB625(2002)345
      CTB11C=DCMPLX(0.D0,0D0)
      CTB12C=DCMPLX(0.D0,0D0)
      CTB21C=DCMPLX(0.D0,0D0)
      CTB22C=DCMPLX(0.D0,0D0)
      DO I=1,2
       DO J=1,2
        DO K=1,2
         DO L=1,2
          CTB11C= 3.D0/(64.D0*PI**2)*C2U
     .          *UT3(I,J)*CP1TB(J,K)*UB3(K,L)*DCONJG(CP1TB(I,L))+CTB11C
          CTB12C= 3.D0/(64.D0*PI**2)*C2U
     .          *UT3(I,J)*CP1TB(J,K)*UB3(K,L)*DCONJG(CP2TB(I,L))+CTB12C
          CTB21C= 3.D0/(64.D0*PI**2)*C2U
     .          *UT3(I,J)*CP2TB(J,K)*UB3(K,L)*DCONJG(CP1TB(I,L))+CTB21C
          CTB22C= 3.D0/(64.D0*PI**2)*C2U
     .          *UT3(I,J)*CP2TB(J,K)*UB3(K,L)*DCONJG(CP2TB(I,L))+CTB22C
         ENDDO
        ENDDO
       ENDDO
      ENDDO
*      print*,'>> GET_MASQ 9 <<',CTB11C
*      print*,'>> GET_MASQ 9 <<',CTB12C
*      print*,'>> GET_MASQ 9 <<',CTB21C
*      print*,'>> GET_MASQ 9 <<',CTB22C
*CTAD..  : tadpoles including seagull graphs 
*Eq.(B.13) and T^(d) in (B.16) NPB625(2002)345. Also see Eq.(2.5)
      DELTA(1,1)=1.D0
      DELTA(1,2)=0.D0
      DELTA(2,1)=0.D0
      DELTA(2,2)=1.D0
*
      CTAD11=DCMPLX(0.D0,0D0)
      CTAD12=DCMPLX(0.D0,0D0)
      CTAD21=DCMPLX(0.D0,0D0)
      CTAD22=DCMPLX(0.D0,0D0)
      DO I=1,2
       DO J=1,2
      CTAD11 = 3.D0/(32.D0*PI**2*V1_SF)
     .       *(A0_H(MSQ_ST2,QSF2)+A0_H(MSQ_ST1,QSF2))
     .       *(CF1TT(I,J)-V1_SF*CPP1TT(I,J))*DELTA(J,I)
     .        +3.D0/(32.D0*PI**2*V1_SF)
     .       *(A0_H(MSQ_ST2,QSF2)-A0_H(MSQ_ST1,QSF2))
     .       *UT3(I,J)*(CF1TT(J,I)-V1_SF*CPP1TT(J,I))
     .        +3.D0/(32.D0*PI**2*V1_SF)
     .       *(A0_H(MSQ_SB2,QSF2)+A0_H(MSQ_SB1,QSF2))
     .       *(CF1BB(I,J)-V1_SF*CPP1BB(I,J))*DELTA(J,I)
     .        +3.D0/(32.D0*PI**2*V1_SF)
     .       *(A0_H(MSQ_SB2,QSF2)-A0_H(MSQ_SB1,QSF2))
     .       *UB3(I,J)*(CF1BB(J,I)-V1_SF*CPP1BB(J,I))
     .       + CTAD11
      CTAD22 = 3.D0/(32.D0*PI**2*V2_SF)
     .       *(A0_H(MSQ_ST2,QSF2)+A0_H(MSQ_ST1,QSF2))
     .       *(CF2TT(I,J)-V2_SF*CPP2TT(I,J))*DELTA(J,I)
     .        +3.D0/(32.D0*PI**2*V2_SF)
     .       *(A0_H(MSQ_ST2,QSF2)-A0_H(MSQ_ST1,QSF2))
     .       *UT3(I,J)*(CF2TT(J,I)-V2_SF*CPP2TT(J,I))
     .        +3.D0/(32.D0*PI**2*V2_SF)
     .       *(A0_H(MSQ_SB2,QSF2)+A0_H(MSQ_SB1,QSF2))
     .       *(CF2BB(I,J)-V2_SF*CPP2BB(I,J))*DELTA(J,I)
     .        +3.D0/(32.D0*PI**2*V2_SF)
     .       *(A0_H(MSQ_SB2,QSF2)-A0_H(MSQ_SB1,QSF2))
     .       *UB3(I,J)*(CF2BB(J,I)-V2_SF*CPP2BB(J,I))
     .       + CTAD22
      CTAD21 =  DCMPLX(0.D0,3.D0)/(32.D0*PI**2*V2_SF) ! iT_a1/v2
     .       *(A0_H(MSQ_ST2,QSF2)-A0_H(MSQ_ST1,QSF2))
     .       *UT3(I,J)*CA1TT(J,I)
     .         +DCMPLX(0.D0,3.D0)/(32.D0*PI**2*V2_SF)
     .       *(A0_H(MSQ_SB2,QSF2)-A0_H(MSQ_SB1,QSF2))
     .       *UB3(I,J)*CA1BB(J,I)
     .       + CTAD21
      CTAD12 =  DCMPLX(0.D0,3.D0)/(32.D0*PI**2*V1_SF) ! iT_a2/v1
     .       *(A0_H(MSQ_ST2,QSF2)-A0_H(MSQ_ST1,QSF2))
     .       *UT3(I,J)*CA2TT(J,I)
     .         +DCMPLX(0.D0,3.D0)/(32.D0*PI**2*V1_SF)
     .       *(A0_H(MSQ_SB2,QSF2)-A0_H(MSQ_SB1,QSF2))
     .       *UB3(I,J)*CA2BB(J,I)
     .       + CTAD12
       ENDDO
      ENDDO
*      print*,'>> GET_MASQ 10 <<',CTAD11
*      print*,'>> GET_MASQ 10 <<',CTAD12
*      print*,'>> GET_MASQ 10 <<',CTAD21
*      print*,'>> GET_MASQ 10 <<',CTAD22

      CPI_SQ(1,1)=(CTB11A+CTB11B+CTB11C+CTAD11)/XI1_SF/XI1_SF
      CPI_SQ(1,2)=(CTB12A+CTB12B+CTB12C+CTAD12)/XI1_SF/XI2_SF
      CPI_SQ(2,1)=(CTB21A+CTB21B+CTB21C+CTAD21)/XI2_SF/XI1_SF
      CPI_SQ(2,2)=(CTB22A+CTB22B+CTB22C+CTAD22)/XI2_SF/XI2_SF
*      print*,CTB11A,CTB11B,CTB11C,CTAD11,XI1_SF,XI1_SF
*      print*,CTB12A,CTB12B,CTB12C,CTAD12,XI1_SF,XI2_SF
*      print*,CTB21A,CTB21B,CTB21C,CTAD21,XI2_SF,XI1_SF
*      print*,CTB22A,CTB22B,CTB22C,CTAD22,XI2_SF,XI2_SF
*
*-------------------------------------------------------------------------
* CMCH_1L(I,J) = -CM0_1L(I,J)-CPI_SQ(I,J)-CPI_Q(I,J) is a 2X2 charged 
* Higgs-boson mass matrix in (phi_1^pm,phi_2^pm) basis dropping the two-loop 
* Born-improved term at the scale m_t-pole
*
      CMCH_1L(1,1)=-CM0_1L(1,1)-CPI_SQ(1,1)-CPI_Q(1,1)
      CMCH_1L(1,2)=-CM0_1L(1,2)-CPI_SQ(1,2)-CPI_Q(1,2)
      CMCH_1L(2,1)=-CM0_1L(2,1)-CPI_SQ(2,1)-CPI_Q(2,1)
      CMCH_1L(2,2)=-CM0_1L(2,2)-CPI_SQ(2,2)-CPI_Q(2,2)
*      print*,CM0_1L(1,1),CPI_SQ(1,1),CPI_Q(1,1)
*      print*,CM0_1L(1,2),CPI_SQ(1,2),CPI_Q(1,2)
*      print*,CM0_1L(2,1),CPI_SQ(2,1),CPI_Q(2,1)
*      print*,CM0_1L(2,2),CPI_SQ(2,2),CPI_Q(2,2)

*.....CMCH_1L(H^+,H^-)
      CMCH22=CMCH_1L(1,1)*SB_H**2-(CMCH_1L(1,2)+CMCH_1L(2,1))*SB_H*CB_H
     .      +CMCH_1L(2,2)*CB_H**2
*      print*,CMCH_1L(1,1)*SB_H**2,(CMCH_1L(1,2)+CMCH_1L(2,1))*SB_H*CB_H
*     .      ,CMCH_1L(2,2)*CB_H**2
*      print*,'>> GET_MASQ 13 <<',CMCH22
*-------------------------------------------------------------------------
* MASQ = MCHpole^2 - BARL4/2*V^2 + Re[CPICH(H^+,H^-)]

      CPICH22 =-CMCH22
      MASQ    = MCH**2-BARL4*(V1**2+V2**2)/2.D0+DREAL(CPICH22)
      REPI22  = DREAL(CPICH22) 
      BL4VSQ  = BARL4*(V1**2+V2**2)/2.D0
*      print*,MCH**2,BARL4*(V1**2+V2**2)/2.D0,DREAL(CPICH22)
*      print*,'>> GET_MASQ 15 <<',s,MASQ
*-------------------------------------------------------------------------
************************************************************************
      RETURN
      END

*      REAL*8 FUNCTION F_I(A,B,C)
************************************************************************
*JSL 18/APR/2005, A=B, B=C, C=A, A=B=C cases considered.
*JSL:28/AUG/2006 : improved treatment for the degenerate cases
************************************************************************
*      IMPLICIT REAL*8(A-H,M,O-Z)
**
*      EPS=1.D-6
*      IF(DABS((B-A)/A).LT.EPS) B=A
*      IF(DABS((C-A)/A).LT.EPS) C=A
*      IF(DABS((C-B)/B).LT.EPS) B=C
**
*      IF(A.EQ.B .AND. B.EQ.C .AND. C.EQ.A) THEN
*       F_I=1.D0/2.D0/A
*      ELSEIF(A.EQ.B) THEN
*       F_I=(B-C+C*DLOG(C/B))/(B-C)**2
*      ELSEIF(B.EQ.C) THEN
*       F_I=(C-A+A*DLOG(A/C))/(C-A)**2
*      ELSEIF(C.EQ.A) THEN
*       F_I=(A-B+B*DLOG(B/A))/(A-B)**2
*      ELSE
*       F_I=(A*B*DLOG(A/B)+B*C*DLOG(B/C)+C*A*DLOG(C/A))
*     .    /((A-B)*(B-C)*(A-C))
*      ENDIF
**
*      RETURN
*      END



      REAL*8 FUNCTION A0_H(A1,Q2)
************************************************************************
* A_0(a) at the 't-Hooft scale Q**2
************************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      A0_H=A1*(1.D0+B0_H(0.D0,A1,A1,Q2))
      RETURN
      END

      REAL*8 FUNCTION B0_H(S,A1,A2,Q2)
************************************************************************
* Dispersive part of B_0(s,a,b) at the 't-Hooft scale Q**2
************************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16 (C)
      IF (S.GT.0.D0) THEN
      B0_H=DREAL(CB0_H(S,A1,A2,Q2))
c      B0_H=DIMAG(CB0_H(S,A1,A2,Q2))
      ELSEIF (A1.EQ.A2) THEN
      B0_H = -DLOG(DSQRT(A1*A2)/Q2)
      ELSE
      B0_H=1.D0-DLOG(DSQRT(A1*A2)/Q2)+0.5D0*(A1+A2)/(A1-A2)*DLOG(A2/A1)
      ENDIF
      RETURN
      END

      COMPLEX*16 FUNCTION CB0_H(S,A1,A2,Q2)
************************************************************************
*  B_0(q**2, m_1**2, m_2**2) at the 't-Hooft scale Q**2
************************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16 (C)
      DATA PI/3.141592653589793238462643D0/
      ARCH(X)=DLOG(X+DSQRT(X**2-1.D0))
      X=Q2
      B1=DSQRT(A1)
      B2=DSQRT(A2)
      AI=0.D0
ccc
      IF (A1.EQ.0.D0.AND.A2.EQ.0.D0) THEN
       CB0_H=-CLN_H(DCMPLX(-S/X),DCMPLX(-1.D0))+2.D0
      ELSE IF (A1.EQ.0.D0.OR.A2.EQ.0.D0) THEN
       AM=DMAX1(A1,A2)
       IF (S.EQ.0.D0) THEN
        H=-1.D0
       ELSE IF (S.LT.AM) THEN
        H=(AM/S-1.D0)*DLOG(1.D0-S/AM)
       ELSE IF (S.EQ.AM) THEN
        H=0.D0
       ELSE
        H=(AM/S-1.D0)*DLOG(S/AM-1.D0)
        AI=-PI*(AM/S-1.D0)
       END IF
       R=DLOG(X/AM)+2.D0+H
       CB0_H=DCMPLX(R,AI)
      ELSE IF (S.EQ.0.D0.AND.A1.EQ.A2) THEN
       CB0_H=DCMPLX(DLOG(X/A1))
      ELSE IF (S.EQ.0.D0) THEN
       CB0_H=DCMPLX(DLOG(X/(B1*B2))+1.D0-(A1+A2)/(A1-A2)*DLOG(B1/B2))
      ELSE
       A=(S-A1-A2)**2-4.D0*A1*A2
       B=0.5D0*(A1+A2-S)/(B1*B2)
       IF (S.LE.(B1-B2)**2) THEN
        H=DSQRT(A)*ARCH(B)
       ELSE IF (S.GT.(B1+B2)**2) THEN
        H=-DSQRT(A)*ARCH(-B)
        AI=PI*DSQRT(A)/S
       ELSE
        H=-DSQRT(-A)*DACOS(B)
       END IF
       R=DLOG(X/(B1*B2))+2.D0+(-(A1-A2)*DLOG(B1/B2)+H)/S
       CB0_H=DCMPLX(R,AI)
      END IF
      RETURN
      END

      COMPLEX*16 FUNCTION CLN_H(CX,CE)
************************************************************************
* ln{x + i*epsilon*sign[Re(e)]}
************************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16 (C)
      DATA PI/3.141592653589793238462643D0/
      IF (DREAL(CX).LT.0.D0.AND.DIMAG(CX).EQ.0.D0) THEN
      CLN_H=CDLOG(-CX)+DSIGN(1.D0,DREAL(CE))*(0.D0,1.D0)*PI
      ELSE
      CLN_H=CDLOG(CX)
      END IF
      RETURN
      END

      COMPLEX*16 FUNCTION CB1_H(S,A1,A2,Q2)
************************************************************************
*  B_1(q**2, m_1**2, m_2**2) at the 't-Hooft scale Q**2
************************************************************************
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16 (C)
      DATA PI/3.141592653589793238462643D0/
      IF (S.GT.0.d0) THEN
      CB1_H=(A2-A1)*(CB0_H(S,A1,A2,Q2)-CB0_H(0.D0,A1,A2,Q2))/(2.D0*S)
     #-0.5D0*CB0_H(S,A1,A2,Q2)
      ELSE
      CB1_H=-0.25D0+0.5D0*DLOG(A2/Q2)+0.5D0*A1**2/(A1-A2)**2*DLOG(A1/A2)
     #-0.5D0*A1/(A1-A2)
      ENDIF
      RETURN
      END

      SUBROUTINE GET_CMNH(SPRO,MCH,MASQ,NFLAG,IFLAG_H
     .                   ,QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2
     .                   ,XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF
     .                   ,V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF
     .                   ,HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF
     .                   ,HT_CP,HB_CP,BARL1,BARL2,BARL34
     .                   ,CMNH)
************************************************************************
*
* RG-improved 4x4 Complex Neurtal Higgs Mass matrix 
* in (phi_1, phi_2, a, G) basis
*
* M^2_ij = M^2(0)_ij - M^2(0)[1 loop]_ij/(xi_i xi_j)
*         - Pi-hat_ij(Squarks)/(xi_i xi_j) - Pi-hat_ij(Quarks)
*
* CMNH(I,J) = CM0(I,J)-CM0_1L(I,J)-CPI_SQ(I,J)-CPI_Q(I,J)
*
* CM0    : Two-loop Born-improved mass matrix at the scale m_t
*          Eqs. (2.21-24) (3.29) of NPB586(2000)92
* CM0_1L : One-loop part of tree-level-form mass matrix at the scale Q
*          after multiplying anomalous dim. factors XI_I and XI_J
* CPI_SQ : Squark contributions to the self energy at Q
*          after multiplying anomalous dim. factors XI_I and XI_J
*          Eqs. (2.11-14), (B.5)(a), (B.11)(a), (B.6)(b), (B.16)(d)
*          of NPB625(2002)345
* CPI_Q  : Quark contributions to the self energy at m_t
*          Eqs. (B.14)(c), (B.16)(e) of NPB625(2002)345
*
*---> SPRO,MCH,MASQ          ! S and Charged Higgs-boson pole mass
*                              and the would-be CP-odd scalar
*---> QQT2,QTT2,QST2,QQB2,QBB2,QSB2,QSF2   ! S and sfermion scales 
*---> XI1_ST,XI1_SB,XI1_SF,XI2_ST,XI2_SB,XI2_SF ! anomalous dimensions
*       : ST=Stop mass scale
*       : SB=Sbottom mass scale
*       : SF=Sfermion mass scale
*---> V1,V1_ST,V1_SB,V1_SF,V2,V2_ST,V2_SB,V2_SF ! vevs
*       : V1, V2 at Mtpole
*---> HT,HT_MT,HT_ST,HT_SF,HB,HB_MT,HB_SB,HB_SF ! Yukawa couplings
*       : HT   , HB    at Mtpole without threshold corrections
*       : HT_MT, HB_MT at Mtpole with    threshold corrections
*---> HT_CP,HB_CP                               ! CP phases of Yukawa couplings
*       : HT_X(complex)=HT_X*HT_CP with X=none, MT, ST, SF
*       : HB_X(complex)=HB_X*HB_CP with X=none, MT, SB, SF
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 IFLAG_H(NFLAG)
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
*-----------------------------------------------------------------------
*Local
      COMPLEX*16 HT_CP,HB_CP
      COMPLEX*16 CF1TT(2,2),CF2TT(2,2),CF1BB(2,2),CF2BB(2,2),CA1TT(2,2),
     .           CA2TT(2,2),CA1BB(2,2),CA2BB(2,2),CP1TB(2,2),CP2TB(2,2),
     .           CPP1TT(2,2),CPP2TT(2,2),CPP1BB(2,2),CPP2BB(2,2)
      COMPLEX*16 CHTR,CHBR,CXI
      COMPLEX*16 UT3(2,2),UB3(2,2)
      COMPLEX*16 C0U,C1UB,C1UT,C2U
      COMPLEX*16 CB0_H
      COMPLEX*16 CM0(4,4),CM0_1L(4,4),CPI_SQ(4,4),CPI_Q(4,4)
      COMPLEX*16 CSTA(4,4),CSTB(4,4),CSTC(4,4),CSTS(4,4),CSTT(4,4)
      COMPLEX*16 CSBA(4,4),CSBB(4,4),CSBC(4,4),CSBS(4,4),CSBT(4,4)
      COMPLEX*16 CST_F,CSB_F
      COMPLEX*16 CST_FP,CST_FM,CSB_FP,CSB_FM
      COMPLEX*16 CMNH(4,4)
*
      REAL*8     XIST(4),XISB(4),UROT(4,4)
*
      PI=2.D0*DASIN(1.D0)
      CXI=DCMPLX(0.D0,1.D0)
*
************************************************************************
      S=SPRO
*at top-quark pole mass scale
      GS2 =4.D0*PI*ASMT_H
*-------------------------------------------------------------------------
* CM0    : Two-loop Born-improved mass matrix at the scale m_t
*           in (phi_1, phi_2, a, G) basis
*          Eqs. (2.21-24) (3.29) of NPB586(2000)92
*
*.....1-loop couplings at M_t^pole : 
*For M_{ij}^2(0) [Eq.(3.15) NPB586]
      IF(IFLAG_H(12).EQ.3 .OR. IFLAG_H(12).EQ.5) THEN
       HBX=HB_MT
       HTX=HT_MT
      ELSE
       HBX=HB
       HTX=HT
      ENDIF
*
      X1L1=-3.D0/(32.D0*PI**2)*( 
     . (GW_H**2/4.D0-GP_H**2/12.D0)**2*DLOG(QQT2/MTPOLE_H**2)+
     .  GP_H**4/9.D0*DLOG(QTT2/MTPOLE_H**2)+
     . (HBX**2-GW_H**2/4.D0-GP_H**2/12.D0)**2*DLOG(QQB2/MTPOLE_H**2)+
     . (HBX**2-GP_H**2/6.D0)**2*DLOG(QBB2/MTPOLE_H**2) )
      X1L2=-3.D0/(32.D0*PI**2)*( 
     . (HTX**2-GW_H**2/4.D0+GP_H**2/12.D0)**2*DLOG(QQT2/MTPOLE_H**2)+
     . (HTX**2-GP_H**2/3.D0)**2*DLOG(QTT2/MTPOLE_H**2)+
     . (GW_H**2/4.D0+GP_H**2/12.D0)**2*DLOG(QQB2/MTPOLE_H**2)+
     . GP_H**4/36.D0*DLOG(QBB2/MTPOLE_H**2) )
      X1L3=-3.D0/(16.D0*PI**2)*(HTX**2*HBX**2*(DLOG(QQT2/MTPOLE_H**2)
     .                         +DLOG(DMAX1(QTT2,QBB2)/MTPOLE_H**2))
     .-((GW_H**2/4.D0+GP_H**2/12.D0)*(HTX**2-GW_H**2/4.D0)-
     .   GP_H**2/12.D0*(GW_H**2/4.D0-GP_H**2/12.D0))*
     .                                         DLOG(QQT2/MTPOLE_H**2)
     .-((GW_H**2/4.D0-GP_H**2/12.D0)*(HBX**2-GW_H**2/4.D0)+
     .   GP_H**2/12.D0*(GW_H**2/4.D0+GP_H**2/12.D0) )*
     .                                         DLOG(QQB2/MTPOLE_H**2)
     .+GP_H**2/3.D0*(HTX**2-GP_H**2/3.D0)*DLOG(QTT2/MTPOLE_H**2)
     .+GP_H**2/6.D0*(HBX**2-GP_H**2/6.D0)*DLOG(QBB2/MTPOLE_H**2)  )
      X1L4= 3.D0/(16.D0*PI**2)*(HTX**2*HBX**2*(DLOG(QQT2/MTPOLE_H**2)
     .                         +DLOG(DMAX1(QTT2,QBB2)/MTPOLE_H**2))
     ,-GW_H**2/2.D0*(HTX**2-GW_H**2/4.D0)*DLOG(QQT2/MTPOLE_H**2)
     .-GW_H**2/2.D0*(HBX**2-GW_H**2/4.D0)*DLOG(QQB2/MTPOLE_H**2) )
      X1L34 = X1L3 + X1L4

*Threshold corrections have NOT included in two-loop couplings
*.....2-loop couplings at M_t^pole
      X2L1=-6.D0*HB**4/(32.D0*PI**2)**2
     .     *(3.D0*HB**2/2.D0+HT**2/2.D0-8.D0*GS2)
     .     *((DLOG(QQB2/MTPOLE_H**2))**2+(DLOG(QBB2/MTPOLE_H**2))**2)
      X2L2=-6.D0*HT**4/(32.D0*PI**2)**2
     .     *(3.D0*HT**2/2.D0+HB**2/2.D0-8.D0*GS2)
     .     *((DLOG(QQT2/MTPOLE_H**2))**2+(DLOG(QTT2/MTPOLE_H**2))**2)
      X2L34 = 0.D0
*.....Computation of the chargino and neutrlino contributions to the
*.....couplings \lambda_1, \lambda_2, \lambda_{34} = \lambda_3+\lambda_4,
*.....and \lambda_4: XINO1, XINO2, XINO34, XINO4
*.....In this determination the formulas in Appendix C of H.E. Haber and
*.....R. Hempfling, Phys. Rev. D48 (1993) 4280, are used.
      TWINO=0.D0
      THINO=0.D0
      TCHI1=0.D0
      TCHI2=0.D0
      TCHI12=0.D0
*
      IF( QSF2.GT.CDABS(M2_H)**2 ) THEN
       TWINO=-DLOG( QSF2/CDABS(M2_H)**2 )
       IF(MTPOLE_H.GT.CDABS(M2_H)) TWINO=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.CDABS(MU_H)**2 ) THEN
       THINO=-DLOG( QSF2/CDABS(MU_H)**2 )
       IF(MTPOLE_H.GT.CDABS(MU_H)) THINO=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.DMAX1(CDABS(MU_H)**2,CDABS(M1_H)**2) ) THEN
       TCHI1=-DLOG( QSF2/DMAX1(CDABS(MU_H)**2,CDABS(M1_H)**2) )
       IF(MTPOLE_H.GT.DMAX1(CDABS(MU_H),CDABS(M1_H)))
     .  TCHI1=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( QSF2.GT.DMAX1(CDABS(MU_H)**2,CDABS(M2_H)**2) ) THEN
       TCHI2=-DLOG( QSF2/ DMAX1(CDABS(MU_H)**2,CDABS(M2_H)**2) )
       IF(MTPOLE_H.GT.DMAX1(CDABS(MU_H),CDABS(M1_H)))
     .  TCHI2=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      IF( DSQRT(QSF2).GT.
     .    DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) ) THEN
       TCHI12=-2.D0*DLOG(DSQRT(QSF2)
     .     /DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) )
       IF( MTPOLE_H.GT.
     .     DMAX1(CDABS(MU_H),DMAX1(CDABS(M1_H),CDABS(M2_H))) )
     . TCHI12=-DLOG( QSF2/MTPOLE_H**2 )
      ENDIF
*
      XINO1=-GW_H**2**2/(768.D0*PI**2*CW_H**2**2)*( 
     . 6.D0*SW_H**2*(1.D0-2.D0*SW_H**2)*TCHI1
     .-24.D0*SW_H**2*CW_H**2*TCHI12
     .-(42.D0-102.D0*SW_H**2+60.D0*SW_H**2**2)*TCHI2
     .-4.D0*(SW_H**2**2+CW_H**2**2)*THINO-8.D0*CW_H**2**2*TWINO )
      XINO2 = XINO1
      XINO34=GW_H**2**2/(384.D0*PI**2*CW_H**2**2)*( 
     . 6.D0*SW_H**2*(1.D0+2.D0*SW_H**2)*TCHI1
     .+24.D0*SW_H**2*CW_H**2*TCHI12
     .+(30.D0-42.D0*SW_H**2+12.D0*SW_H**2**2)*TCHI2
     .-4.D0*(SW_H**2**2+CW_H**2**2)*THINO-8.D0*CW_H**2**2*TWINO )
      XINO4=GW_H**4/(192.D0*PI**2)*(
     . 6.D0*SW_H**2/CW_H**2*TCHI1 + 24.D0*SW_H**2/CW_H**2*TCHI12
     .-6.D0*TCHI2 - 4.D0*THINO - 8.D0*TWINO )
*.....Two-loop improved couplings
      BARL1 = -(GW_H**2+GP_H**2)/8.D0 + X1L1  + X2L1  + XINO1
      BARL2 = -(GW_H**2+GP_H**2)/8.D0 + X1L2  + X2L2  + XINO2
      BARL34 = (GW_H**2+GP_H**2)/4.D0 + X1L34 + X2L34 + XINO34
*.....Finally, in the (phi_1, phi_2, a, G) basis
      DO I=1,4
        DO J=1,4
         CM0(I,J)=DCMPLX(0.D0,0.D0)
        ENDDO
      ENDDO
      CM0(1,1)=DCMPLX( MASQ*SB_H**2   - 2.D0*BARL1*V1**2,0.D0)
      CM0(2,2)=DCMPLX( MASQ*CB_H**2   - 2.D0*BARL2*V2**2,0.D0)
      CM0(1,2)=DCMPLX(-MASQ*SB_H*CB_H - BARL34*V1*V2    ,0.D0)
      CM0(2,1)=CM0(1,2)
      CM0(3,3)=DCMPLX( MASQ,0.D0)
*      print*,'> CM0(1,1) <',cm0(1,1)
*      print*,'> CM0(1,2) <',cm0(1,2)
*      print*,'> CM0(2,2) <',cm0(2,2)
*      print*,'> CM0(3,3) <',cm0(3,3)
*-------------------------------------------------------------------------
* CM0_1L : One-loop part of tree-level-form mass matrix at the scale Q
*          in the (phi_1, phi_2, a_1, a_2) basis. The anomalous dim. factors 
*          XI_I and XI_J are multiplied
*
*.....1-loop couplings at stop and sbottom scales:
*..... : Note that  CM0_1L(I,J)=0 identically when M_Q=M_U=M_D. 
*..... : Yukawa couplings WITH threshold corrections used
      X1L1T=-3.D0/(32.D0*PI**2)*(
     . (GW_H**2/4.D0-GP_H**2/12.D0)**2*DLOG(QQT2/QST2)+
     .  GP_H**4/9.D0*DLOG(QTT2/QST2) )
      X1L1B=-3.D0/(32.D0*PI**2)*(
     . (HB_SB**2-GW_H**2/4.D0-GP_H**2/12.D0)**2*DLOG(QQB2/QSB2)+
     . (HB_SB**2-GP_H**2/6.D0)**2*DLOG(QBB2/QSB2) )

      X1L2T=-3.D0/(32.D0*PI**2)*(
     . (HT_ST**2-GW_H**2/4.D0+GP_H**2/12.D0)**2*DLOG(QQT2/QST2)+
     . (HT_ST**2-GP_H**2/3.D0)**2*DLOG(QTT2/QST2) )
      X1L2B=-3.D0/(32.D0*PI**2)*(
     . (GW_H**2/4.D0+GP_H**2/12.D0)**2*DLOG(QQB2/QSB2)+
     . GP_H**4/36.D0*DLOG(QBB2/QSB2) )

      X1L3T=-3.D0/(16.D0*PI**2)*(
     .-((GW_H**2/4.D0+GP_H**2/12.D0)*(HT_ST**2-GW_H**2/4.D0)-
     .   GP_H**2/12.D0*(GW_H**2/4.D0-GP_H**2/12.D0))*
     .                                         DLOG(QQT2/QST2)
     .+GP_H**2/3.D0*(HT_ST**2-GP_H**2/3.D0)*DLOG(QTT2/QST2) )
      X1L3B=-3.D0/(16.D0*PI**2)*(
     .-((GW_H**2/4.D0-GP_H**2/12.D0)*(HB_SB**2-GW_H**2/4.D0)+
     .   GP_H**2/12.D0*(GW_H**2/4.D0+GP_H**2/12.D0) )*
     .                                         DLOG(QQB2/QSB2)
     .+GP_H**2/6.D0*(HB_SB**2-GP_H**2/6.D0)*DLOG(QBB2/QSB2)  )

      X1L4T=-3.D0/(16.D0*PI**2)*
     , GW_H**2/2.D0*(HT_ST**2-GW_H**2/4.D0)*DLOG(QQT2/QST2)
      X1L4B=-3.D0/(16.D0*PI**2)*
     . GW_H**2/2.D0*(HB_SB**2-GW_H**2/4.D0)*DLOG(QQB2/QSB2) 

      X1L34T = X1L3T + X1L4T
      X1L34B = X1L3B + X1L4B
*.....Note that XRM12T and XRM12B are always vanishing! 
      XRM12T=3.D0/(16.D0*PI**2)*HT_ST**2*DREAL(MU_H*AT_H)
     .      *DLOG(DMAX1(QQT2,QTT2)/QST2)
      XRM12B=3.D0/(16.D0*PI**2)*HB_SB**2*DREAL(MU_H*AB_H)
     .      *DLOG(DMAX1(QQB2,QBB2)/QSB2)
*.....Finally, in the (phi_1, phi_2, a_1, a_2) basis
      DO I=1,4
        DO J=1,4
         CM0_1L(I,J)=DCMPLX(0.D0,0.D0)
        ENDDO
      ENDDO
      CM0_1L(1,1)= DCMPLX( XRM12T*V2_ST/V1_ST-2.D0*X1L1T*V1_ST**2,0.D0)
     .            /XI1_ST/XI1_ST
     .            +DCMPLX( XRM12B*V2_SB/V1_SB-2.D0*X1L1B*V1_SB**2,0.D0)
     .            /XI1_SB/XI1_SB
      CM0_1L(2,2)= DCMPLX( XRM12T*V1_ST/V2_ST-2.D0*X1L2T*V2_ST**2,0.D0)
     .            /XI2_ST/XI2_ST
     .            +DCMPLX( XRM12B*V1_SB/V2_SB-2.D0*X1L2B*V2_SB**2,0.D0)
     .            /XI2_SB/XI2_SB
      CM0_1L(1,2)= DCMPLX(-XRM12T            -X1L34T*V1_ST*V2_ST ,0.D0)
     .            /XI1_ST/XI2_ST
     .            +DCMPLX(-XRM12B            -X1L34B*V1_SB*V2_SB ,0.D0)
     .            /XI1_SB/XI2_SB
      CM0_1L(2,1)=CM0_1L(1,2)
      CM0_1L(3,3)= DCMPLX( XRM12T*V2_ST/V1_ST,0.D0)/XI1_ST/XI1_ST
     .            +DCMPLX( XRM12B*V2_SB/V1_SB,0.D0)/XI1_SB/XI1_SB
      CM0_1L(4,4)= DCMPLX( XRM12T*V1_ST/V2_ST,0.D0)/XI2_ST/XI2_ST
     .            +DCMPLX( XRM12B*V1_SB/V2_SB,0.D0)/XI2_SB/XI2_SB
      CM0_1L(3,4)= DCMPLX(-XRM12T            ,0.D0)/XI1_ST/XI2_ST
     .            +DCMPLX(-XRM12B            ,0.D0)/XI1_SB/XI2_SB
      CM0_1L(4,3)=CM0_1L(3,4)
*      print*,'> CM0_1L(1,1) <',cm0_1l(1,1)
*      print*,'> CM0_1L(1,2) <',cm0_1l(1,2)
*      print*,'> CM0_1L(2,2) <',cm0_1l(2,2)
*      print*,'> CM0_1L(3,3) <',cm0_1l(3,3)
*      print*,'> CM0_1L(3,4) <',cm0_1l(3,4)
*      print*,'> CM0_1L(4,4) <',cm0_1l(4,4)
*-------------------------------------------------------------------------
* CPI_SQ : Squark contributions to the self energy at Q in the 
*          (phi_1, phi_2, a_1, a_2) basis. The anomalous dim. factors 
*          XI_I and XI_J are multiplied
*          Eqs. (2.11-14), (B.5)(a), (B.11)(a), (B.6)(b), (B.16)(d)
*          of NPB625(2002)345
*
*.....Squark sector:
      GWY2=GW_H**2/8.D0/CW_H**2
      GXT2=(GW_H**2-5.D0*GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0
      GXB2=(GW_H**2-GW_H**2*SW_H**2/CW_H**2/3.D0)/4.D0

      SUM_ST=MQ3_H**2+MU3_H**2+HT_ST**2*V2_ST**2
     .      +GWY2*(V1_ST**2-V2_ST**2)
      DIF_ST2=(MQ3_H**2-MU3_H**2
     .            +GXT2*(V1_ST**2-V2_ST**2)/2.D0)**2
     .       +2.D0*HT_ST**2*CDABS(DCONJG(AT_H)*V2_ST-MU_H*V1_ST)**2
      MSQ_ST2=(SUM_ST+DSQRT(DIF_ST2))/2.D0
      MSQ_ST1=(SUM_ST-DSQRT(DIF_ST2))/2.D0
      SUM_SB=MQ3_H**2+MD3_H**2+HB_SB**2*V1_SB**2
     .      +GWY2*(V2_SB**2-V1_SB**2)
      DIF_SB2=(MQ3_H**2-MD3_H**2
     .            +GXB2*(V2_SB**2-V1_SB**2)/2.D0)**2
     .       +2.D0*HB_SB**2*CDABS(DCONJG(AB_H)*V1_SB-MU_H*V2_SB)**2
      MSQ_SB2=(SUM_SB+DSQRT(DIF_SB2))/2.D0
      MSQ_SB1=(SUM_SB-DSQRT(DIF_SB2))/2.D0

      IF(MSQ_ST1.LT.0.D0.OR.MSQ_ST2.LT.0.D0.OR.
     .   MSQ_SB1.LT.0.D0.OR.MSQ_SB2.LT.0.D0) THEN
       IFLAG_H(50)=1
       RETURN
      ENDIF

      UT3(1,1)=(MQ3_H**2-MU3_H**2
     .          +GXT2*(V1_ST**2-V2_ST**2)/2.D0)/(MSQ_ST2-MSQ_ST1)
      UT3(2,1)=DSQRT(2.D0)*HT_ST*HT_CP*(AT_H*V2_ST-DCONJG(MU_H)*V1_ST)
     .                          /(MSQ_ST2-MSQ_ST1)
      UT3(1,2)=DCONJG( UT3(2,1) )
      UT3(2,2)=-UT3(1,1)

      UB3(1,1)=(MQ3_H**2-MD3_H**2
     .           -GXB2*(V1_SB**2-V2_SB**2)/2.D0)/(MSQ_SB2-MSQ_SB1)
      UB3(2,1)=DSQRT(2.D0)*HB_SB*HB_CP*(AB_H*V1_SB-DCONJG(MU_H)*V2_SB)
     .                             /(MSQ_SB2-MSQ_SB1)
      UB3(1,2)=DCONJG( UB3(2,1) )
      UB3(2,2)=-UB3(1,1)

*Here charged Higgs-SQURK-SQUARK couplings : 
      CHTR = HT_ST*HT_CP
      CHBR = HB_SB*HB_CP
      V1RT = V1_ST
      V2RT = V2_ST
      V1RB = V1_SB
      V2RB = V2_SB
* The complex Yukawa couplings with threhsold corrections should be used ?:
*=> See the above arguments.
*--->Start
* Gamma{\phi_1 ~t* ~t}:       CF1TT(2,2)                              *
      CF1TT(1,1)=-V1RT*(GW_H**2-GP_H**2/3.D0)/4.D0
      CF1TT(1,2)=DCONJG(CHTR)*MU_H/DSQRT(2.D0)
      CF1TT(2,1)=DCONJG( CF1TT(1,2) )
      CF1TT(2,2)=-V1RT*GP_H**2/3.D0
* Gamma{\phi_2 ~t* ~t}:       CF2TT(2,2)                              *
      CF2TT(1,1)=-CDABS(CHTR)**2*V2RT+V2RT*(GW_H**2-GP_H**2/3.D0)/4.D0
      CF2TT(1,2)=-DCONJG(CHTR*AT_H)/DSQRT(2.D0)
      CF2TT(2,1)=DCONJG( CF2TT(1,2) )
      CF2TT(2,2)=-CDABS(CHTR)**2*V2RT+V2RT*GP_H**2/3.D0
* Gamma{\phi_1 ~b* ~b}:       CF1BB(2,2)                              *
      CF1BB(1,1)=-CDABS(CHBR)**2*V1RB+V1RB*(GW_H**2+GP_H**2/3.D0)/4.D0
      CF1BB(1,2)=-DCONJG(CHBR*AB_H)/DSQRT(2.D0)
      CF1BB(2,1)=DCONJG( CF1BB(1,2) )
      CF1BB(2,2)=-CDABS(CHBR)**2*V1RB+V1RB*GP_H**2/6.D0
* Gamma{\phi_2 ~b* ~b}:       CF2BB(2,2)                              *
      CF2BB(1,1)=-V2RB*(GW_H**2+GP_H**2/3.D0)/4.D0
      CF2BB(1,2)=DCONJG(CHBR)*MU_H/DSQRT(2.D0)
      CF2BB(2,1)=DCONJG( CF2BB(1,2) )
      CF2BB(2,2)=-V2RB*GP_H**2/6.D0
* Gamma{a_1 ~t* ~t}:          CA1TT(2,2)                              *
      CA1TT(1,1)=0.D0
      CA1TT(1,2)=-CXI*DCONJG(CHTR)/DSQRT(2.D0) * MU_H
      CA1TT(2,1)=DCONJG( CA1TT(1,2) )
      CA1TT(2,2)=0.D0
* Gamma{a_2 ~t* ~t}:          CA2TT(2,2)                              *
      CA2TT(1,1)=0.D0
      CA2TT(1,2)=CXI*DCONJG(CHTR)/DSQRT(2.D0) * DCONJG(AT_H)
      CA2TT(2,1)=DCONJG( CA2TT(1,2) )
      CA2TT(2,2)=0.D0
* Gamma{a_1 ~b* ~b}:          CA1BB(2,2)                              *
      CA1BB(1,1)=0.D0
      CA1BB(1,2)=-CXI*DCONJG(CHBR)/DSQRT(2.D0) * DCONJG(AB_H)
      CA1BB(2,1)=DCONJG( CA1BB(1,2) )
      CA1BB(2,2)=0.D0
* Gamma{a_2 ~b* ~b}:          CA2BB(2,2)                              *
      CA2BB(1,1)=0.D0
      CA2BB(1,2)=CXI*DCONJG(CHBR)/DSQRT(2.D0) * MU_H
      CA2BB(2,1)=DCONJG( CA2BB(1,2) )
      CA2BB(2,2)=0.D0
* Gamma{\phi+_1 ~t* ~b}:     CP1TB(2,2)                               *
      CP1TB(1,1)=-V1_SF*(HB_SF**2-GW_H**2/2.D0)/DSQRT(2.D0)
      CP1TB(1,2)=-DCONJG(HB_SF*HB_CP*AB_H)
      CP1TB(2,1)=-HT_SF*HT_CP*DCONJG(MU_H)
      CP1TB(2,2)=-HT_SF*HT_CP*DCONJG(HB_SF*HB_CP)*V2_SF/DSQRT(2.D0)
* Gamma{\phi+_2 ~t* ~b}:     CP2TB(2,2)                               *
      CP2TB(1,1)=V2_SF*(HT_SF**2-GW_H**2/2.D0)/DSQRT(2.D0)
      CP2TB(1,2)=DCONJG(HB_SF*HB_CP)*MU_H
      CP2TB(2,1)=HT_SF*HT_CP*AT_H
      CP2TB(2,2)=HT_SF*HT_CP*DCONJG(HB_SF*HB_CP)*V1_SF/DSQRT(2.D0)
* Gamma{\phi+_1 \phi+_1 ~t* ~t}:  CPP1TT(2,2)                         *
      CPP1TT(1,1)=-HB_SF**2+(GW_H**2+GP_H**2/3.D0)/4.D0
      CPP1TT(1,2)=0.D0
      CPP1TT(2,1)=0.D0
      CPP1TT(2,2)=-GP_H**2/3.D0
* Gamma{\phi+_2 \phi+_2 ~t* ~t}:  CPP2TT(2,2)                         *
      CPP2TT(1,1)=-(GW_H**2+GP_H**2/3.D0)/4.D0
      CPP2TT(1,2)=0.D0
      CPP2TT(2,1)=0.D0
      CPP2TT(2,2)=-HT_SF**2+GP_H**2/3.D0
* Gamma{\phi+_1 \phi+_1 ~b* ~b}:  CPP1BB(2,2)                         *
      CPP1BB(1,1)=-(GW_H**2-GP_H**2/3.D0)/4.D0
      CPP1BB(1,2)=0.D0
      CPP1BB(2,1)=0.D0
      CPP1BB(2,2)=-HB_SF**2+GP_H**2/6.D0
* Gamma{\phi+_2 \phi+_2 ~b* ~b}:  CPP2BB(2,2)                         *
      CPP2BB(1,1)=-HT_SF**2+(GW_H**2-GP_H**2/3.D0)/4.D0
      CPP2BB(1,2)=0.D0
      CPP2BB(2,1)=0.D0
      CPP2BB(2,2)=-GP_H**2/6.D0
*-->EOF

*.....The 1st term of Eqs.(B.5) and (B.11) of NPB625(2002)345 : CSTA and CSBA
      CST_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_ST2,MSQ_ST2,QST2)+
     .                          CB0_H(S,MSQ_ST1,MSQ_ST1,QST2)+
     .                     2.D0*CB0_H(S,MSQ_ST1,MSQ_ST2,QST2))
      CSB_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_SB2,MSQ_SB2,QSB2)+
     .                          CB0_H(S,MSQ_SB1,MSQ_SB1,QSB2)+
     .                     2.D0*CB0_H(S,MSQ_SB1,MSQ_SB2,QSB2))
      DO I=1,4
       DO J=I,4
        CSTA(I,J)=DCMPLX(0.D0,0.D0)
        CSBA(I,J)=DCMPLX(0.D0,0.D0)
         DO I1=1,2
         DO J1=1,2
*
          IF(I.EQ.1.AND.J.EQ.1) 
     .      CSTA(I,J)=CST_F*CF1TT(I1,J1)*CF1TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.1.AND.J.EQ.2) 
     .      CSTA(I,J)=CST_F*CF1TT(I1,J1)*CF2TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.1.AND.J.EQ.3) 
     .      CSTA(I,J)=CST_F*CF1TT(I1,J1)*CA1TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.1.AND.J.EQ.4) 
     .      CSTA(I,J)=CST_F*CF1TT(I1,J1)*CA2TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.2.AND.J.EQ.2) 
     .      CSTA(I,J)=CST_F*CF2TT(I1,J1)*CF2TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.2.AND.J.EQ.3) 
     .      CSTA(I,J)=CST_F*CF2TT(I1,J1)*CA1TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.2.AND.J.EQ.4) 
     .      CSTA(I,J)=CST_F*CF2TT(I1,J1)*CA2TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.3.AND.J.EQ.3) 
     .      CSTA(I,J)=CST_F*CA1TT(I1,J1)*CA1TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.3.AND.J.EQ.4) 
     .      CSTA(I,J)=CST_F*CA1TT(I1,J1)*CA2TT(J1,I1)+CSTA(I,J)
          IF(I.EQ.4.AND.J.EQ.4) 
     .      CSTA(I,J)=CST_F*CA2TT(I1,J1)*CA2TT(J1,I1)+CSTA(I,J)
*
          IF(I.EQ.1.AND.J.EQ.1)
     .      CSBA(I,J)=CSB_F*CF1BB(I1,J1)*CF1BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.1.AND.J.EQ.2)
     .      CSBA(I,J)=CSB_F*CF1BB(I1,J1)*CF2BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.1.AND.J.EQ.3)
     .      CSBA(I,J)=CSB_F*CF1BB(I1,J1)*CA1BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.1.AND.J.EQ.4)
     .      CSBA(I,J)=CSB_F*CF1BB(I1,J1)*CA2BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.2.AND.J.EQ.2)
     .      CSBA(I,J)=CSB_F*CF2BB(I1,J1)*CF2BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.2.AND.J.EQ.3)
     .      CSBA(I,J)=CSB_F*CF2BB(I1,J1)*CA1BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.2.AND.J.EQ.4)
     .      CSBA(I,J)=CSB_F*CF2BB(I1,J1)*CA2BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.3.AND.J.EQ.3)
     .      CSBA(I,J)=CSB_F*CA1BB(I1,J1)*CA1BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.3.AND.J.EQ.4)
     .      CSBA(I,J)=CSB_F*CA1BB(I1,J1)*CA2BB(J1,I1)+CSBA(I,J)
          IF(I.EQ.4.AND.J.EQ.4)
     .      CSBA(I,J)=CSB_F*CA2BB(I1,J1)*CA2BB(J1,I1)+CSBA(I,J)
*
         ENDDO ! J1
         ENDDO ! I1
*         print*,i,j,s,CSTA(I,J)
*         print*,i,j,s,CSBA(I,J)
       ENDDO
      ENDDO
      CSTA(2,1)=CSTA(1,2)
      CSTA(3,1)=CSTA(1,3)
      CSTA(3,2)=CSTA(2,3)
      CSTA(4,1)=CSTA(1,4)
      CSTA(4,2)=CSTA(2,4)
      CSTA(4,3)=CSTA(3,4)
      CSBA(2,1)=CSBA(1,2)
      CSBA(3,1)=CSBA(1,3)
      CSBA(3,2)=CSBA(2,3)
      CSBA(4,1)=CSBA(1,4)
      CSBA(4,2)=CSBA(2,4)
      CSBA(4,3)=CSBA(3,4)
*.....The 2nd term of Eq.(B.5) and (B.11) of NPB625(2002)345 : CSTB and CSBB
      CST_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_ST2,MSQ_ST2,QST2)-
     .                          CB0_H(S,MSQ_ST1,MSQ_ST1,QST2))
      CSB_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_SB2,MSQ_SB2,QSB2)-
     .                          CB0_H(S,MSQ_SB1,MSQ_SB1,QSB2))
      DO I=1,4
       DO J=I,4
        CSTB(I,J)=DCMPLX(0.D0,0.D0)
        CSBB(I,J)=DCMPLX(0.D0,0.D0)
         DO I1=1,2
         DO J1=1,2
         DO K1=1,2
*
          IF(I.EQ.1.AND.J.EQ.1)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF1TT(J1,K1)*CF1TT(K1,I1)+CF1TT(J1,K1)*CF1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.1.AND.J.EQ.2)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF1TT(J1,K1)*CF2TT(K1,I1)+CF2TT(J1,K1)*CF1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.1.AND.J.EQ.3)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF1TT(J1,K1)*CA1TT(K1,I1)+CA1TT(J1,K1)*CF1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.1.AND.J.EQ.4)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF1TT(J1,K1)*CA2TT(K1,I1)+CA2TT(J1,K1)*CF1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.2.AND.J.EQ.2)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF2TT(J1,K1)*CF2TT(K1,I1)+CF2TT(J1,K1)*CF2TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.2.AND.J.EQ.3)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF2TT(J1,K1)*CA1TT(K1,I1)+CA1TT(J1,K1)*CF2TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.2.AND.J.EQ.4)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CF2TT(J1,K1)*CA2TT(K1,I1)+CA2TT(J1,K1)*CF2TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.3.AND.J.EQ.3)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CA1TT(J1,K1)*CA1TT(K1,I1)+CA1TT(J1,K1)*CA1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.3.AND.J.EQ.4)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CA1TT(J1,K1)*CA2TT(K1,I1)+CA2TT(J1,K1)*CA1TT(K1,I1))
     .               +CSTB(I,J)
          IF(I.EQ.4.AND.J.EQ.4)
     .      CSTB(I,J)=CST_F*UT3(I1,J1)
     .      *(CA2TT(J1,K1)*CA2TT(K1,I1)+CA2TT(J1,K1)*CA2TT(K1,I1))
     .               +CSTB(I,J)
*
          IF(I.EQ.1.AND.J.EQ.1)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF1BB(J1,K1)*CF1BB(K1,I1)+CF1BB(J1,K1)*CF1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.1.AND.J.EQ.2)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF1BB(J1,K1)*CF2BB(K1,I1)+CF2BB(J1,K1)*CF1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.1.AND.J.EQ.3)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF1BB(J1,K1)*CA1BB(K1,I1)+CA1BB(J1,K1)*CF1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.1.AND.J.EQ.4)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF1BB(J1,K1)*CA2BB(K1,I1)+CA2BB(J1,K1)*CF1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.2.AND.J.EQ.2)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF2BB(J1,K1)*CF2BB(K1,I1)+CF2BB(J1,K1)*CF2BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.2.AND.J.EQ.3)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF2BB(J1,K1)*CA1BB(K1,I1)+CA1BB(J1,K1)*CF2BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.2.AND.J.EQ.4)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CF2BB(J1,K1)*CA2BB(K1,I1)+CA2BB(J1,K1)*CF2BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.3.AND.J.EQ.3)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CA1BB(J1,K1)*CA1BB(K1,I1)+CA1BB(J1,K1)*CA1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.3.AND.J.EQ.4)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CA1BB(J1,K1)*CA2BB(K1,I1)+CA2BB(J1,K1)*CA1BB(K1,I1))
     .               +CSBB(I,J)
          IF(I.EQ.4.AND.J.EQ.4)
     .      CSBB(I,J)=CSB_F*UB3(I1,J1)
     .      *(CA2BB(J1,K1)*CA2BB(K1,I1)+CA2BB(J1,K1)*CA2BB(K1,I1))
     .               +CSBB(I,J)

*
         ENDDO ! K1
         ENDDO ! J1
         ENDDO ! I1
*         print*,i,j,s,CSTB(I,J)
*         print*,i,j,s,CSBB(I,J)
       ENDDO
      ENDDO
      CSTB(2,1)=CSTB(1,2)
      CSTB(3,1)=CSTB(1,3)
      CSTB(3,2)=CSTB(2,3)
      CSTB(4,1)=CSTB(1,4)
      CSTB(4,2)=CSTB(2,4)
      CSTB(4,3)=CSTB(3,4)
      CSBB(2,1)=CSBB(1,2)
      CSBB(3,1)=CSBB(1,3)
      CSBB(3,2)=CSBB(2,3)
      CSBB(4,1)=CSBB(1,4)
      CSBB(4,2)=CSBB(2,4)
      CSBB(4,3)=CSBB(3,4)
*.....The 3nd term of Eq.(B.5) and (B.11) of NPB625(2002)345 : CSTC and CSBC
      CST_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_ST2,MSQ_ST2,QST2)+
     .                          CB0_H(S,MSQ_ST1,MSQ_ST1,QST2)-
     .                     2.D0*CB0_H(S,MSQ_ST1,MSQ_ST2,QST2))
      CSB_F=3.D0/(64.D0*PI**2)*(CB0_H(S,MSQ_SB2,MSQ_SB2,QSB2)+
     .                          CB0_H(S,MSQ_SB1,MSQ_SB1,QSB2)-
     .                     2.D0*CB0_H(S,MSQ_SB1,MSQ_SB2,QSB2))
      DO I=1,4
       DO J=I,4
        CSTC(I,J)=DCMPLX(0.D0,0.D0)
        CSBC(I,J)=DCMPLX(0.D0,0.D0)
         DO I1=1,2
         DO J1=1,2
         DO K1=1,2
         DO L1=1,2
*
          IF(I.EQ.1.AND.J.EQ.1)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF1TT(J1,K1)*UT3(K1,L1)*CF1TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.1.AND.J.EQ.2)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF1TT(J1,K1)*UT3(K1,L1)*CF2TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.1.AND.J.EQ.3)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF1TT(J1,K1)*UT3(K1,L1)*CA1TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.1.AND.J.EQ.4)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF1TT(J1,K1)*UT3(K1,L1)*CA2TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.2.AND.J.EQ.2)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF2TT(J1,K1)*UT3(K1,L1)*CF2TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.2.AND.J.EQ.3)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF2TT(J1,K1)*UT3(K1,L1)*CA1TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.2.AND.J.EQ.4)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CF2TT(J1,K1)*UT3(K1,L1)*CA2TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.3.AND.J.EQ.3)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CA1TT(J1,K1)*UT3(K1,L1)*CA1TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.3.AND.J.EQ.4)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CA1TT(J1,K1)*UT3(K1,L1)*CA2TT(L1,I1)+CSTC(I,J)
          IF(I.EQ.4.AND.J.EQ.4)
     .      CSTC(I,J)=CST_F
     .       *UT3(I1,J1)*CA2TT(J1,K1)*UT3(K1,L1)*CA2TT(L1,I1)+CSTC(I,J)
*
          IF(I.EQ.1.AND.J.EQ.1)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF1BB(J1,K1)*UB3(K1,L1)*CF1BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.1.AND.J.EQ.2)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF1BB(J1,K1)*UB3(K1,L1)*CF2BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.1.AND.J.EQ.3)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF1BB(J1,K1)*UB3(K1,L1)*CA1BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.1.AND.J.EQ.4)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF1BB(J1,K1)*UB3(K1,L1)*CA2BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.2.AND.J.EQ.2)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF2BB(J1,K1)*UB3(K1,L1)*CF2BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.2.AND.J.EQ.3)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF2BB(J1,K1)*UB3(K1,L1)*CA1BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.2.AND.J.EQ.4)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CF2BB(J1,K1)*UB3(K1,L1)*CA2BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.3.AND.J.EQ.3)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CA1BB(J1,K1)*UB3(K1,L1)*CA1BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.3.AND.J.EQ.4)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CA1BB(J1,K1)*UB3(K1,L1)*CA2BB(L1,I1)+CSBC(I,J)
          IF(I.EQ.4.AND.J.EQ.4)
     .      CSBC(I,J)=CSB_F
     .       *UB3(I1,J1)*CA2BB(J1,K1)*UB3(K1,L1)*CA2BB(L1,I1)+CSBC(I,J)
*
         ENDDO ! L1
         ENDDO ! K1
         ENDDO ! J1
         ENDDO ! I1
*         print*,i,j,s,CSTC(I,J)
*         print*,i,j,s,CSBC(I,J)
       ENDDO
      ENDDO
      CSTC(2,1)=CSTC(1,2)
      CSTC(3,1)=CSTC(1,3)
      CSTC(3,2)=CSTC(2,3)
      CSTC(4,1)=CSTC(1,4)
      CSTC(4,2)=CSTC(2,4)
      CSTC(4,3)=CSTC(3,4)
      CSBC(2,1)=CSBC(1,2)
      CSBC(3,1)=CSBC(1,3)
      CSBC(3,2)=CSBC(2,3)
      CSBC(4,1)=CSBC(1,4)
      CSBC(4,2)=CSBC(2,4)
      CSBC(4,3)=CSBC(3,4)
*.....Seagull graphs : Eq.(B.6) and (B.11) of NPB625(2002)345 : CSTS and CSBS
      CST_FP=3.D0/32.D0/PI**2*(A0_H(MSQ_ST2,QST2)+A0_H(MSQ_ST1,QST2))
      CST_FM=3.D0/32.D0/PI**2*(A0_H(MSQ_ST2,QST2)-A0_H(MSQ_ST1,QST2))
      CSB_FP=3.D0/32.D0/PI**2*(A0_H(MSQ_SB2,QSB2)+A0_H(MSQ_SB1,QSB2))
      CSB_FM=3.D0/32.D0/PI**2*(A0_H(MSQ_SB2,QSB2)-A0_H(MSQ_SB1,QSB2))
      DO I=1,4
       DO J=1,4
        CSTS(I,J)=DCMPLX(0.D0,0.D0)
        CSBS(I,J)=DCMPLX(0.D0,0.D0)
       ENDDO
      ENDDO
      CSTS(1,1)=-CST_FP*(CF1TT(1,1)+CF1TT(2,2))/V1_ST
     .          -CST_FM*(UT3(1,1)*CF1TT(1,1)+UT3(2,2)*CF1TT(2,2))/V1_ST
      CSTS(2,2)=-CST_FP*(CF2TT(1,1)+CF2TT(2,2))/V2_ST
     .          -CST_FM*(UT3(1,1)*CF2TT(1,1)+UT3(2,2)*CF2TT(2,2))/V2_ST
      CSTS(3,3)=CSTS(1,1)
      CSTS(4,4)=CSTS(2,2)
*
      CSBS(1,1)=-CSB_FP*(CF1BB(1,1)+CF1BB(2,2))/V1_SB
     .          -CSB_FM*(UB3(1,1)*CF1BB(1,1)+UB3(2,2)*CF1BB(2,2))/V1_SB
      CSBS(2,2)=-CSB_FP*(CF2BB(1,1)+CF2BB(2,2))/V2_SB
     .          -CSB_FM*(UB3(1,1)*CF2BB(1,1)+UB3(2,2)*CF2BB(2,2))/V2_SB
      CSBS(3,3)=CSBS(1,1)
      CSBS(4,4)=CSBS(2,2)
*.....Tadpole graphs: Eq.(B.6) and (B.11) of NPB625(2002)345 : CSTT and CSBT
      CST_FP=3.D0/32.D0/PI**2*(A0_H(MSQ_ST2,QST2)+A0_H(MSQ_ST1,QST2))
      CST_FM=3.D0/32.D0/PI**2*(A0_H(MSQ_ST2,QST2)-A0_H(MSQ_ST1,QST2))
      CSB_FP=3.D0/32.D0/PI**2*(A0_H(MSQ_SB2,QSB2)+A0_H(MSQ_SB1,QSB2))
      CSB_FM=3.D0/32.D0/PI**2*(A0_H(MSQ_SB2,QSB2)-A0_H(MSQ_SB1,QSB2))
      DO I=1,4
       DO J=1,4
        CSTT(I,J)=DCMPLX(0.D0,0.D0)
        CSBT(I,J)=DCMPLX(0.D0,0.D0)
       ENDDO
      ENDDO
      CSTT(1,1)=( CST_FP*(CF1TT(1,1)+CF1TT(2,2))
     .           +CST_FM*(UT3(1,1)*CF1TT(1,1)+UT3(1,2)*CF1TT(2,1)
     .                   +UT3(2,1)*CF1TT(1,2)+UT3(2,2)*CF1TT(2,2)) )
     .          /V1_ST
      CSTT(2,2)=( CST_FP*(CF2TT(1,1)+CF2TT(2,2))
     .           +CST_FM*(UT3(1,1)*CF2TT(1,1)+UT3(1,2)*CF2TT(2,1)
     .                   +UT3(2,1)*CF2TT(1,2)+UT3(2,2)*CF2TT(2,2)) )
     .          /V2_ST
      CSTT(3,3)=CSTT(1,1)
      CSTT(4,4)=CSTT(2,2)
*
      CSTT(1,4)=( CST_FM*(UT3(1,1)*CA2TT(1,1)+UT3(1,2)*CA2TT(2,1)
     .                   +UT3(2,1)*CA2TT(1,2)+UT3(2,2)*CA2TT(2,2)) )
     .          /V1_ST
      CSTT(2,3)=-CSTT(1,4)
      CSTT(3,2)=CSTT(2,3)
      CSTT(4,1)=CSTT(1,4)
*
      CSBT(1,1)=( CSB_FP*(CF1BB(1,1)+CF1BB(2,2))
     .           +CSB_FM*(UB3(1,1)*CF1BB(1,1)+UB3(1,2)*CF1BB(2,1)
     .                   +UB3(2,1)*CF1BB(1,2)+UB3(2,2)*CF1BB(2,2)) )
     .          /V1_SB
      CSBT(2,2)=( CSB_FP*(CF2BB(1,1)+CF2BB(2,2))
     .           +CSB_FM*(UB3(1,1)*CF2BB(1,1)+UB3(1,2)*CF2BB(2,1)
     .                   +UB3(2,1)*CF2BB(1,2)+UB3(2,2)*CF2BB(2,2)) )
     .          /V2_SB
      CSBT(3,3)=CSBT(1,1)
      CSBT(4,4)=CSBT(2,2)
*
      CSBT(1,4)=( CSB_FM*(UB3(1,1)*CA2BB(1,1)+UB3(1,2)*CA2BB(2,1)
     .                   +UB3(2,1)*CA2BB(1,2)+UB3(2,2)*CA2BB(2,2)) )
     .          /V1_SB
      CSBT(2,3)=-CSBT(1,4)
      CSBT(3,2)=CSBT(2,3)
      CSBT(4,1)=CSBT(1,4)
*
*      print*,CSTS(1,1)+CSTT(1,1)
*      print*,CSTS(1,4)+CSTT(1,4)
*      print*,CSTS(2,2)+CSTT(2,2)
*      print*,CSTS(2,3)+CSTT(2,3)
*      print*,CSTS(3,3)+CSTT(3,3)
*      print*,CSTS(4,4)+CSTT(4,4)
*      print*,CSBS(1,1)+CSBT(1,1)
*      print*,CSBS(1,4)+CSBT(1,4)
*      print*,CSBS(2,2)+CSBT(2,2)
*      print*,CSBS(2,3)+CSBT(2,3)
*      print*,CSBS(3,3)+CSBT(3,3)
*      print*,CSBS(4,4)+CSBT(4,4)
*.....Finally, in the (phi_1, phi_2, a_1, a_2) basis
      XIST(1)=XI1_ST
      XIST(2)=XI2_ST
      XIST(3)=XI1_ST
      XIST(4)=XI2_ST
      XISB(1)=XI1_SB
      XISB(2)=XI2_SB
      XISB(3)=XI1_SB
      XISB(4)=XI2_SB
      DO I=1,4
        DO J=1,4
         CPI_SQ(I,J)=
     .     (CSTA(I,J)+CSTB(I,J)+CSTC(I,J)+CSTS(I,J)+CSTT(I,J))
     .     /XIST(I)/XIST(J)
     .    +(CSBA(I,J)+CSBB(I,J)+CSBC(I,J)+CSBS(I,J)+CSBT(I,J))
     .     /XISB(I)/XISB(J)
*        print*,i,j,s,cm0_1l(i,j)+cpi_sq(i,j)
        ENDDO
      ENDDO
*-------------------------------------------------------------------------
*CPI_Q at mtpole : Quark-loop contributions to the Neutral Higgs-boson 
*self-energies. See  Eqs. (B.14)(c)+(B.16)(e) of NPB625(2002)345
*
*NOTE: The overall minus signs are missing in PI^P in (B.14) 
*
*For the quark masses inside the loop and the t'Hooft scale, we should use
*the same conventions for CPI_Q in GET_MASQ. The mixed uses of m_t-pole 
*and m_t(m_t^pole) will be clarified later <-- ERROR ??
*
      IF(IFLAG_H(12).EQ.2 .OR. IFLAG_H(12).EQ.5) THEN
       HBX=HB_MT
       HTX=HT_MT
      ELSE
       HBX=HB
       HTX=HT
      ENDIF
*
      DO I=1,4
       DO J=1,4
        CPI_Q(I,J)=DCMPLX(0.D0,0.D0)
        CPI_Q(I,J)=DCMPLX(0.D0,0.D0)
       ENDDO
      ENDDO
*
      CPI_Q(1,1)=3.D0*HBX**2/16.D0/PI**2
     .    *(S-4.D0*MBMT_H**2)*CB0_H(S,MBMT_H**2,MBMT_H**2,MTPOLE_H**2)
      CPI_Q(2,2)=3.D0*HTX**2/16.D0/PI**2
     .    *(S-4.D0*MTMT_H**2)*CB0_H(S,MTMT_H**2,MTMT_H**2,MTMT_H**2)
      CPI_Q(3,3)=3.D0*HBX**2/16.D0/PI**2
     .    *S*CB0_H(S,MBMT_H**2,MBMT_H**2,MTPOLE_H**2)
      CPI_Q(4,4)=3.D0*HTX**2/16.D0/PI**2
     .    *S*CB0_H(S,MTPOLE_H**2,MTPOLE_H**2,MTPOLE_H**2)
*      print*,cpi_q(1,1)
*      print*,cpi_q(2,2)
*      print*,cpi_q(3,3)
*      print*,cpi_q(4,4)
*
*-------------------------------------------------------------------------
* CMNH(I,J) = CM0(I,J)-CM0_1L(I,J)-CPI_SQ(I,J)-CPI_Q(I,J)
*          in (phi_1, phi_2, a, G) basis:
* Rotation from (phi_1, phi_2, a_1, a_2) to (phi_1, phi_2, a, G) needed
*          for CM0_IL, CPI_SQ, and CPI_Q
* (phi_1 phi_2 a_1 a_2)^T = U (phi_1 phi_2 a G)^T

* where        / 1  0  0  0 \   
*          U = | 0  1  0  0 |  and U^-1 = U
*              | 0  0 -s  c |
*              \ 0  0  c  s /
      DO I=1,4
       DO J=1,4
        UROT(I,J) = 0.D0
       ENDDO
      ENDDO
      UROT(1,1) = 1.D0
      UROT(2,2) = 1.D0
      UROT(3,3) =-SB_H
      UROT(4,4) = SB_H
      UROT(3,4) = CB_H
      UROT(4,3) = CB_H
*
      DO I=1,4
       DO J=1,4
        CMNH(I,J) = CM0(I,J)
         DO I1=1,4
         DO J1=1,4
          CMNH(I,J) =-UROT(I,I1)
     .               *(CM0_1L(I1,J1)+CPI_SQ(I1,J1)+CPI_Q(I1,J1))
     .               *UROT(J1,J)+CMNH(I,J)
         ENDDO ! J1
         ENDDO ! I1
*        print*,i,j,s,cmnh(i,j)
       ENDDO
      ENDDO
*-------------------------------------------------------------------------
************************************************************************
      RETURN
      END


      SUBROUTINE DUMP_HIGGS(NFLAG,IFLAG_H,MCH,HMASS,OMIX)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8 HMASS(3),OMIX(3,3)
      INTEGER*8 IFLAG_H(NFLAG)
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
*-----------------------------------------------------------------------
      print*,'---------------------------------------------------------'
      print*,' Masses and Mixing Matrix of Higgs bosons :'
      print*,'                               HMASS_H(I) and OMIX_H(A,I)'
      print*,'---------------------------------------------------------'
      DO IH=1,3
      IF(IFLAG_H(11).EQ.0) WRITE(*,1) IH,HMASS(IH)
      IF(IFLAG_H(11).EQ.1) WRITE(*,6) IH,HMASS(IH)
      ENDDO
      IF(IFLAG_H(11).EQ.0) WRITE(*,2) MCH
      IF(IFLAG_H(11).EQ.1) WRITE(*,7) MCH
      print*,'                         [H1]        [H2]        [H3]'
      WRITE(*,3) OMIX(1,1),OMIX(1,2),OMIX(1,3)
      WRITE(*,4) OMIX(2,1),OMIX(2,2),OMIX(2,3)
      WRITE(*,5) OMIX(3,1),OMIX(3,2),OMIX(3,3)
      print*,'---------------------------------------------------------'
      print*,' '
*-----------------------------------------------------------------------
 1    FORMAT(2X,'H',I1,'  Pole Mass           = ',E10.4,' GeV')
 2    FORMAT(2X,'Charged Higgs Pole Mass = ',E10.4,' GeV [SSPARA_H(2)]')
 3    FORMAT(2X,'         [phi_1] [',3(1X,E10.4,1X),' ]')
 4    FORMAT(2X,'O(IA,IH)=[phi_2] [',3(1X,E10.4,1X),' ]')
 5    FORMAT(2X,'         [  a  ] [',3(1X,E10.4,1X),' ]')
 6    FORMAT(2X,'H',I1,'  Eff. Pot. Mass      = ',E10.4,' GeV')
 7    FORMAT(2X,'C. Higgs Eff. Pot. Mass = ',E10.4,' GeV ')
*-----------------------------------------------------------------------
      RETURN
      END


      SUBROUTINE INCL_STAU(NFLAG,IFLAG_H,HMASS,OMIX)
************************************************************************
*
* Include stau contribution to the Higgs masses and mixing by modifying
* the sbottom contributions as appeared as in:
*
* Phys.Lett. B481 (2000) 57-66, e-Print: hep-ph/0002287
* by Choi, Drees, Lee
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
*Input/Output Arrays
      INTEGER*8 IFLAG_H(NFLAG)
      REAL*8    HMASS(3),OMIX(3,3)
*-----------------------------------------------------------------------
*Local
      REAL*8    HMASS_IN(3),OMIX_IN(3,3)
      REAL*8    HMASS_OUT(3),OMIX_OUT(3,3)
      REAL*8    DMH3_STAU(3,3)
      REAL*8    NH3(3,3),EV3(3),AUX3(3)
*
* slepton mass and mixining matrices:
      REAL*8     STAUMASS(2),SNU3MASS
      COMPLEX*16 STAUMIX(2,2)
* Radiative corrections to Htautau Yukawa couplings
      COMPLEX*16 HB_H,HT_H,HTAU_H
      COMPLEX*16 CKTAU_H
*
      COMPLEX*16 XRL,XLR
*-----------------------------------------------------------------------
      PI=2.D0*DASIN(1.D0)
*
      DO I=1,3
       HMASS_IN(I)=HMASS(I)
       DO J=1,3
       OMIX_IN(I,J)=OMIX(I,J)
       ENDDO
      ENDDO
*
*      print*,(hmass_in(i),i=1,3)
*      print*,(omix_in(1,i),i=1,3)
*      print*,(omix_in(2,i),i=1,3)
*      print*,(omix_in(3,i),i=1,3)
*
      RAUX_H(511)=HMASS_IN(1)
      RAUX_H(512)=HMASS_IN(2)
      RAUX_H(513)=HMASS_IN(3)
      RAUX_H(520)=OMIX_IN(1,1)
      RAUX_H(521)=OMIX_IN(1,2)
      RAUX_H(522)=OMIX_IN(1,3)
      RAUX_H(523)=OMIX_IN(2,1)
      RAUX_H(524)=OMIX_IN(2,2)
      RAUX_H(525)=OMIX_IN(2,3)
      RAUX_H(526)=OMIX_IN(3,1)
      RAUX_H(527)=OMIX_IN(3,2)
      RAUX_H(528)=OMIX_IN(3,3)
***********************************************************************
*masses of stau leptons:
*
      HTAU_H=DCMPLX(DSQRT(2.D0)*MTAU_H/V_H/CB_H,0.D0)
      HB_H  =DCMPLX(DSQRT(2.D0)*MBMT_H/V_H/CB_H,0.D0)
      HT_H  =DCMPLX(DSQRT(2.D0)*MTMT_H/V_H/SB_H,0.D0)
      CALL SLMIX(HB_H,HT_H,HTAU_H,STAUMASS,SNU3MASS,STAUMIX)
*----
*      print*,'>> INCL_STAU:',htau_h,staumass(1)
*-------------------------------------------------------------------------*
*---- JSL[21/JAN/2004] Radiative correction  to H-tau-tau           ------*
*----      Eq.(6) of hep-ph/0106027 by Guasch, Hollik, and Penaranda.-----*
*----      CKTAU_H = (Delta m_tau/tb) with M -> M^* and mu -> mu^*  ------*
*-------------------------------------------------------------------------*
      DO ITER_ITAU=1,10
*
      IF(IFLAG_H(10).EQ.0) THEN
       CKTAU_H =DCONJG(MU_H)*AEM_H/4.D0/PI*
     . (
     .  -DCONJG(M2_H)/SW_H**2*
     .   (F_I(SNU3MASS**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .   +F_I(STAUMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(1,1)*DCONJG(STAUMIX(1,1))/2.D0
     .   +F_I(STAUMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(1,2)*DCONJG(STAUMIX(1,2))/2.D0)  ! M2_H
     .  +DCONJG(M1_H)/CW_H**2*
     .   (F_I(STAUMASS(1)**2,STAUMASS(2)**2,CDABS(M1_H)**2)
     .   +F_I(STAUMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(1,1)*DCONJG(STAUMIX(1,1))/2.D0
     .   +F_I(STAUMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(1,2)*DCONJG(STAUMIX(1,2))/2.D0
     .   -F_I(STAUMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(2,1)*DCONJG(STAUMIX(2,1))
     .   -F_I(STAUMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *STAUMIX(2,2)*DCONJG(STAUMIX(2,2)))       ! M1_H
     . )
      ELSE
       CKTAU_H =DCMPLX(0.D0,0.D0)
      ENDIF
*
       HTAU_H=DSQRT(2.D0)*MTAU_H/V_H/CB_H/(1.D0+CKTAU_H*TB_H)
       CALL SLMIX(HB_H,HT_H,HTAU_H,STAUMASS,SNU3MASS,STAUMIX)
       IF(SNU3MASS.LE.0.D0 .OR. STAUMASS(1).LE.0.D0) THEN
         IFLAG_H(56)=1
         RETURN
       ENDIF
*      print*,'>> INCL_STAU:',iter_itau,htau_h,staumass(1),staumass(2)
      ENDDO ! ITER_ITAU
*
*Slepton scale
      QL2  = DMAX1(ML3_H**2,ME3_H**2)
*
      HTAU=CDABS(HTAU_H)
      HB  =CDABS(HB_H)
      HT  =CDABS(HT_H)
*----
      XI1 = 1.D0+3.D0*HB**2/(32.D0*PI**2)*DLOG(QL2/MTPOLE_H**2)
      XI2 = 1.D0+3.D0*HT**2/(32.D0*PI**2)*DLOG(QL2/MTPOLE_H**2)
      V1  = V_H*CB_H
      V2  = V_H*SB_H
*Running V1 and V2 at slepton scale
      V1R = V1/XI1
      V2R = V2/XI2
*
      V_R    = DSQRT(V1R**2+V2R**2)
      SB_R   = DSQRT(V2R**2/(V1R**2+V2R**2))
      CB_R   = DSQRT(V1R**2/(V1R**2+V2R**2))
      MTAU_R = HTAU*V_R*CB_R/DSQRT(2.D0)
*      print*,mtau_h,mtau_r
*
      XLL = ML3_H**2+MTAU_R**2
     .     +MZ_H**2*(CB_R**2-SB_R**2)*(SW_H**2-1.D0/2.D0)
      XRR = ME3_H**2+MTAU_R**2
     .     -MZ_H**2*SW_H**2*(CB_R**2-SB_R**2)
      XRL = MTAU_R*(ATAU_H-DCONJG(MU_H)*SB_R/CB_R)
      XLR = DCONJG(XRL)
*
      DELTA = DSQRT((XLL-XRR)**2+4.D0*CDABS(XRL)**2)
      MAVG  = (XLL+XRR)/2.D0
      IF(MAVG-DELTA/2.D0.LE.0.D0) THEN
*       print*,'Negative Stau mass squared'
       IFLAG_H(56)=1
       RETURN
      ENDIF
*
      MSTAU1=DSQRT(DABS(MAVG-DELTA/2.D0))
      MSTAU2=DSQRT(DABS(MAVG+DELTA/2.D0))
*      print*,mstau1,mstau2
***********************************************************************
      Q_0=MTPOLE_H
      GHSQ=(GW_H**2+GP_H**2)/4.D0
      XTAU=(3.D0*GP_H**2-GW_H**2)/4.D0
     .    *(ML3_H**2-ME3_H**2)/(MSTAU2**2-MSTAU1**2)
      DTAU=-DIMAG(ATAU_H*MU_H)/(MSTAU2**2-MSTAU1**2)
      GSTAU_12=2.D0-(MSTAU1**2+MSTAU2**2)/(MSTAU1**2-MSTAU2**2)
     .             *DLOG(MSTAU1**2/MSTAU2**2)
      RTAU =(CDABS(MU_H)**2*TB_H-DREAL(ATAU_H*MU_H))
     .     /(MSTAU2**2-MSTAU1**2)
      RTAUP=(CDABS(ATAU_H)**2-DREAL(ATAU_H*MU_H)*TB_H)
     .     /(MSTAU2**2-MSTAU1**2)
*      print*,q_0,ghsq,xtau,dtau
*      print*,gstau_12,rtau,rtaup
      F_STAU1=2.D0*MSTAU1**2*(DLOG(MSTAU1**2/Q_0**2)-1.D0)
      F_STAU2=2.D0*MSTAU2**2*(DLOG(MSTAU2**2/Q_0**2)-1.D0)
      FSTAU_12=1.D0/32.D0/PI**2*(F_STAU1-F_STAU2)/(MSTAU2**2-MSTAU1**2)
*      print*,f_stau1,f_stau2,fstau_12
*
*(phi_1,phi_1) component
*
      DMH3_STAU(1,1)=MTAU_R**2/8.D0/PI**2*(
     . HTAU**2*DLOG(MSTAU1**2*MSTAU2**2/MTAU_R**4)
     .-GHSQ*DLOG(MSTAU1**2*MSTAU2**2/Q_0**4)
     .+GSTAU_12*RTAUP*(HTAU**2*RTAUP+XTAU)
     .+DLOG(MSTAU2**2/MSTAU1**2)*(XTAU+(2.D0*HTAU**2-GHSQ)*RTAUP)
     . )
*      print*,'(1,1)',dmh3_stau(1,1),dsqrt(dabs(dmh3_stau(1,1)))
*
*(phi_1,phi_2) component
*
      DMH3_STAU(1,2)=MTAU_R**2/8.D0/PI**2*(
     . GSTAU_12*(HTAU**2*RTAU*RTAUP
     .          +XTAU/2.D0*(RTAU-RTAUP*TB_H))
     .+GHSQ/2.D0*TB_H*DLOG(MSTAU1**2*MSTAU2**2/Q_0**4)
     .+DLOG(MSTAU2**2/MSTAU1**2)*(HTAU**2*RTAU-XTAU/2.D0*TB_H
     .                           +GHSQ/2.D0*(RTAUP*TB_H-RTAU))
     . )
      DMH3_STAU(2,1)=DMH3_STAU(1,2)
*      print*,'(1,2)',dmh3_stau(1,2),dsqrt(dabs(dmh3_stau(1,2)))
*
*(phi_2,phi_2) component
*
      DMH3_STAU(2,2)=MTAU_R**2/8.D0/PI**2*(
     . GSTAU_12*RTAU*(HTAU**2*RTAU-TB_H*XTAU)
     .+GHSQ*TB_H*RTAU*DLOG(MSTAU2**2/MSTAU1**2)
     . )
*      print*,'(2,2)',dmh3_stau(2,2),dsqrt(dabs(dmh3_stau(2,2)))
*
*(a,a) component
*
      DMH3_STAU(3,3)=
     .-1.D0/CB_H/SB_H*HTAU**2*DREAL(ATAU_H*MU_H)*FSTAU_12
     .+1.D0/8.D0/PI**2*HTAU**2*MTAU_R**2/CB_H**2*GSTAU_12*DTAU**2
*      print*,'(3,3)',dmh3_stau(3,3)
*NEGLECTING the contribution to M_A
      DMH3_STAU(3,3)=0.D0
*
*(phi_1,a) component
      DMH3_STAU(1,3)=1.D0/16.D0/PI**2*MTAU_R**2*DTAU/CB_H
     .*(-GSTAU_12*(XTAU+2.D0*HTAU**2*RTAUP)
     .  +(GHSQ-2.D0*HTAU**2)*DLOG(MSTAU2**2/MSTAU1**2))
      DMH3_STAU(3,1)=DMH3_STAU(1,3)
*      print*,dmh3_stau(1,3),dmh3_stau(3,1)
*
*(phi_2,a) component
      DMH3_STAU(2,3)=1.D0/16.D0/PI**2*MTAU_R**2*DTAU/CB_H
     .*( GSTAU_12*(XTAU*TB_H-2.D0*HTAU**2*RTAU)
     .  -GHSQ*TB_H*DLOG(MSTAU2**2/MSTAU1**2))
      DMH3_STAU(3,2)=DMH3_STAU(2,3)
*      print*,dmh3_stau(2,3),dmh3_stau(3,2)
*
*      print*,(dmh3_stau(1,i),i=1,3)
*      print*,(dmh3_stau(2,i),i=1,3)
*      print*,(dmh3_stau(3,i),i=1,3)
*-------------------------------------------------------------------------
      DO I=1,3
       DO J=1,3
        NH3(I,J)=OMIX_IN(I,1)*HMASS_IN(1)**2*OMIX_IN(J,1)
     .          +OMIX_IN(I,2)*HMASS_IN(2)**2*OMIX_IN(J,2)
     .          +OMIX_IN(I,3)*HMASS_IN(3)**2*OMIX_IN(J,3)
     .          +DMH3_STAU(I,J)
       ENDDO
      ENDDO
**
*      print*,((hmass_in(i)),i=1,3)
*      print*,(nh3(1,i),i=1,3)
*      print*,(nh3(2,i),i=1,3)
*      print*,(nh3(3,i),i=1,3)
      CALL DIAGRS(3,3,NH3,EV3,AUX3,IERR_DIAGRS)
      IF(EV3(1).LE.0.D0.OR.EV3(2).LE.0.D0.OR.EV3(3).LE.0.D0) THEN
*         print*,'MHSQ(EFF) =',ev3(1),ev3(2),ev3(3)
        IFLAG_H(51)=1
        RETURN
      ENDIF
      IF(IERR_DIAGRS.GT.0) THEN
*        print*,'DIAGRS Error at sqrt(s) = 0'
        IFLAG_H(52)=1
        RETURN
      ENDIF
**
*      print*,'--------------------------------------------------'
*      print*,(dsqrt(ev3(i)),i=1,3)
*      print*,'--------------------------------------------------'
*      print*,(nh3(1,i),i=1,3)
*      print*,(nh3(2,i),i=1,3)
*      print*,(nh3(3,i),i=1,3)
*
*      print*,'--------------------------------------------------'
*      ddd=-htau**4*v_h**2*cdabs(mu_h)**4/96.d0/pi**2
*     .    /mstau1**2/mstau2**2/2.d0/hmass_in(1)
*      print*,tb_h,dsqrt(ev3(1))-hmass_in(1),ddd
*-------------------------------------------------------------------------
      DO I=1,3
       HMASS_OUT(I)=DSQRT(EV3(I))
       DO J=1,3
        OMIX_OUT(I,J)=NH3(I,J)
       ENDDO
      ENDDO
*
*Print results
      MCH_OUT=RAUX_H(10)
      IF(IFLAG_H(2).EQ.1)
     .            CALL DUMP_HIGGS(NFLAG,IFLAG_H
     .                           ,MCH_OUT,HMASS_OUT,OMIX_OUT)
*Return results
      DO I=1,3
       HMASS(I)=HMASS_OUT(I)
       DO J=1,3
        OMIX(I,J)=OMIX_OUT(I,J)
       ENDDO
      ENDDO
***********************************************************************
*
      RETURN
      END
