      SUBROUTINE FILLCOUPL2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .                     ,MCH,HMASS,OMIX,STMASS,STMIX,SBMASS,SBMIX
     .                     ,STAUMASS,STAUMIX,SNU3MASS,M_C,U_L,U_R
     .                     ,M_N,N_N,NCMAX,NHC_H,SHC_H,CHC_H)
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
*
* Higgs mass and mixining matrices:
      REAL*8     MCH,HMASS(3),OMIX(3,3)
* squark mass and mixining matrices:
      REAL*8     STMASS(2),SBMASS(2)
      COMPLEX*16 STMIX(2,2),SBMIX(2,2)
* slepton mass and mixining matrices:
      REAL*8     STAUMASS(2),SNU3MASS
      COMPLEX*16 STAUMIX(2,2)
* chargino/neutralino mass and mixing matrices:
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN
      REAL*8     M_C(2),M_N(4)
      COMPLEX*16 U_L(2,2),U_R(2,2),N_N(4,4)
* Higgs potential couplings
      REAL*8     LR_H(4)
      COMPLEX*16 LC_H(3)
* H-tau-tau Yukawa couplings
      COMPLEX*16 HTAU_H
* Radiative corrections to Htt Hbb Yukawa couplings
      COMPLEX*16 HB_H,HT_H,CKB_H,CKT_H
      COMPLEX*16 RB_H,RT_H,CKBB_H,CKBT_H
*JSL 10/Jun/2009: Including threshold corrections
      COMPLEX*16 CKD_H,CKS_H
*H^+-Sneutrino3-Stau couplings
      COMPLEX*16 CHN3TUL,CHN3TUR
* Radiative corrections to Htautau Yukawa couplings
      COMPLEX*16 CKTAU_H,CKBEW_H,CKB_SUM
*-----------------------------------------------------------------------
* Input Array:
      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
      INTEGER*8 IFLAG_H(NFLAG)
*-----------------------------------------------------------------------
* Ouput Array:
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*-----------------------------------------------------------------------
* Local Variables:
      COMPLEX*16 V_C(2,2),XI,GP1,GP2,GSF,GPF
      COMPLEX*16 SGG(7),PGG(3)
      COMPLEX*16 SPP(15),PPP(7)
      COMPLEX*16 SGGSM(3),PGGSM(3)
      COMPLEX*16 SPPSM(9),PPPSM(9)
*-----------------------------------------------------------------------
      PI = 2.D0*DASIN(1.D0)
      XI = DCMPLX(0.D0,1.D0)
*-----------------------------------------------------------------------
*--> Finite radiatve corrections to H-t-t and H-b-b Yukawa couplings.
      CALL RADNHTB(NFLAG,IFLAG_H,SBMASS,STMASS,HB_H,HT_H,CKB_H,CKT_H)
*       print*,'FILLCOUPL2:CKB_H',ckb_h
* JSL[14/AUG/2011]: CAUX_H(10)=CKB_H -> CAUX_H(10)=CKB_SUM 
      CAUX_H(212)=CKT_H/TB_H/(1.D0+CKT_H/TB_H)
*--> Stop and sbottom masses and mixing
      CALL SQMIX(HB_H,HT_H,STMASS,SBMASS,STMIX,SBMIX)
*JSL[10/JUN/2009] Tachyonic stop or sbottom
      IF(STMASS(1).LE.0.D0 .OR. SBMASS(1).LE.0.D0) THEN
        IFLAG_H(56)=2
        RETURN
      ENDIF
      IF(IFLAG_H(3).EQ.1) CALL DUMP_SQ(STMASS,SBMASS,STMIX,SBMIX)
*--> Stau and sneutrino3 masses and mixing
      HTAU_H=DCMPLX(DSQRT(2.D0)*MTAU_H/V_H/CB_H,0.D0)  ! initially
      CALL SLMIX(HB_H,HT_H,HTAU_H,STAUMASS,SNU3MASS,STAUMIX)
      IF(SNU3MASS.LE.0.D0 .OR. STAUMASS(1).LE.0.D0) THEN
        IFLAG_H(56)=1
        RETURN
      ENDIF
*
*      print*,'INITIA:STAU',htau_h,staumass(1),staumix(1,2)
*      print*,'INITIA:STAU',htau_h,staumass(1)
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
*      print*,'>> FILLCOUPL2:'
*     .       ,iter_itau,htau_h,staumass(1),staumix(1,2)
*      print*,'>> FILLCOUPL2:',iter_itau,htau_h,staumass(1)
      ENDDO ! ITER_ITAU
*
      CAUX_H(13)=HTAU_H
      IF(IFLAG_H(3).EQ.1) CALL DUMP_SL(STAUMASS,SNU3MASS,STAUMIX)
*
*--> Finite radiatve corrections to H^pm-t-b Yukawa couplings.
      CALL RADCHTB(NFLAG,IFLAG_H,SBMIX,SBMASS,STMIX,STMASS
     .            ,RB_H,RT_H,CKBB_H,CKBT_H)
*
*      print*,'>> FILLCOUPL'
*      print*,DSQRT(2.D0)*MBMT_H/V_H/CB_H,DSQRT(2.D0)*MTMT_H/V_H/SB_H
*      print*,ckb_h,ckt_h
*      print*,rb_h,rt_h
*      print*,ckbb_h,ckbt_h
*
*The differences between the following quantities are due to, for example,
*        (1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))
*       *(1.D0-2.D0*BHB*DLOG(QB2/MTPOLE_H**2)) =\= 1 .
*In other words, it's becasue of the uncertainity in the higher-order 
*corrections. For details, look into the subroutines RADNHTB and
*FILLHIGGS2.
*      print*,cdabs(hb_h),raux_h(27)
*      print*,cdabs(ht_h),raux_h(24)
*      print*,hb_h/cdabs(hb_h),caux_h(2)
*      print*,ht_h/cdabs(ht_h),caux_h(1)
*
*JSL 10/Jun/2009: Including threshold corrections in first two
*generations
       B6    = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
       R_123Q=SSPARA_H(22)
       R_123D=SSPARA_H(24)
       MQ12  =R_123Q*MQ3_H
       MD12  =R_123D*MD3_H
       Q12SQ =DMAX1(MQ12**2,MD12**2)
       AS_M12=ASMT_H/(1.D0+B6*ASMT_H*DLOG(Q12SQ/MTPOLE_H**2))
       CKS_H=2.D0*AS_M12/3.D0/PI*DCONJG(MU_H*M3_H)
     .           *F_I(MD12**2,MQ12**2,CDABS(M3_H)**2)
       CKD_H=CKS_H
*       print*,'>>> FILLCOUPL2: cks_h',cks_h
       CAUX_H(11)=CKS_H
*starng-quark Yukawa couping
       CAUX_H(12)=DSQRT(2.D0)*MSMT_H/V_H/CB_H/(1.D0+CKS_H*TB_H)
*tachyonic scalar strange
       ADmBC=MQ12**2*MD12**2
     .      -(CDABS(CAUX_H(12))*V_H*CB_H/DSQRT(2.D0))**2
     .      *(CDABS(    SSPARA_H(37)*DCOS(SSPARA_H(38)/180.D0*PI)
     .              -XI*SSPARA_H(37)*DSIN(SSPARA_H(38)/180.D0*PI)
     .              -MU_H*TB_H
     .             )
     .       )**2
*       print*,'FILLCOUPL2',SSPARA_H(37)*DCOS(SSPARA_H(38)/180.D0*PI)
*     .                 -XI*SSPARA_H(37)*DSIN(SSPARA_H(38)/180.D0*PI)
*     .       ,CDABS(CAUX_H(12))*V_H*CB_H/DSQRT(2.D0),ADmBC
       IF(ADmBC.LE.0.D0) THEN
        IFLAG_H(56)=3
        RETURN
       ENDIF
*-----------------------------------------------------------------------
*--> Chargino and neutralino masses and mixing
* For /WEINBERG/ S2W_CN,MW_CN,MZ_CN
      S2W_CN = SW_H**2
      MW_CN  = MW_H
      MZ_CN  = MZ_H
*
      TH1  = SSPARA_H(6)*PI/180.D0
      TH2  = SSPARA_H(8)*PI/180.D0
      THMU = SSPARA_H(4)*PI/180.D0
*
      CALL CHARDIAG(CDABS(M2_H),CDABS(MU_H),TB_H
     .             ,TH2,THMU,M_C,U_L,V_C)
        U_R(1,1)=DCONJG(V_C(1,1))
        U_R(1,2)=DCONJG(V_C(1,2))
        U_R(2,1)=DCONJG(V_C(2,1))
        U_R(2,2)=DCONJG(V_C(2,2))
*
      CALL NEUTDIAG(CDABS(M1_H),CDABS(M2_H),CDABS(MU_H),TB_H
     .             ,TH1,TH2,THMU,M_N,N_N)
*
      IF(IFLAG_H(4).EQ.1) 
     .CALL DUMP_CN(TB_H,CDABS(M1_H),CDABS(M2_H),CDABS(MU_H)
     .            ,TH1,TH2,THMU,M_C,M_N,U_L,U_R,N_N)
*-----------------------------------------------------------------------
*--> Higgs Potential Couplings
      CALL HPLAMBDA(HB_H,HT_H,STMASS,LR_H,LC_H)
*      print*,LR_H(1),LR_H(2),LR_H(3),LR_H(4)
*      print*,LC_H(1),LC_H(2),LC_H(3)
      CAUX_H(201)=DCMPLX(LR_H(1),0.D0)
      CAUX_H(202)=DCMPLX(LR_H(2),0.D0)
      CAUX_H(203)=DCMPLX(LR_H(3),0.D0)
      CAUX_H(204)=DCMPLX(LR_H(4),0.D0)
      CAUX_H(205)=LC_H(1)
      CAUX_H(206)=LC_H(2)
      CAUX_H(207)=LC_H(3)
*-----------------------------------------------------------------------
*===> Neutral HiggsIH coupling to two particles : NHC_H(NC,IH)    <====*
*-----------------------------------------------------------------------
      DO IH=1,3
*>>>>> Quarks masses at MH
      CALL QMASS_MH(HMASS(IH),ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH)
      CALL QMASS_MH(HMASS(IH)/2.D0,ASMH2,MTMH2,MBMH2,MCMH2
     .                                  ,MSMH2,MUMH2,MDMH2)
*|--->Higgs-electron-electron                             [NC= 1, 2, 3]
      NHC_H( 1,IH)=DCMPLX(GW_H*ME_H/2.D0/MW_H,0.D0)
      NHC_H( 2,IH)=DCMPLX(OMIX(1,IH)/CB_H,0.D0)
      NHC_H( 3,IH)=DCMPLX(-OMIX(3,IH)*SB_H/CB_H,0.D0)
*|--->Higgs-muon-muon                                     [NC= 4, 5, 6]
      NHC_H( 4,IH)=DCMPLX(GW_H*MMU_H/2.D0/MW_H,0.D0)
      NHC_H( 5,IH)=DCMPLX(OMIX(1,IH)/CB_H,0.D0)
      NHC_H( 6,IH)=DCMPLX(-OMIX(3,IH)*SB_H/CB_H,0.D0)
*|--->Higgs-tau-tau                                       [NC= 7, 8, 9]
      NHC_H( 7,IH)=DCMPLX(GW_H*MTAU_H/2.D0/MW_H,0.D0)
*----
*      NHC_H( 8,IH)=DCMPLX(OMIX(1,IH)/CB_H,0.D0)
*      NHC_H( 9,IH)=DCMPLX(-OMIX(3,IH)*SB_H/CB_H,0.D0)
*      print*,'G^S_tau, G^P_tau : Before including correction',
*     .        NHC_H(8,IH),NHC_H(9,IH)
*----
      NHC_H( 8,IH)=DCMPLX(
     . DREAL(1.D0/(1.D0+CKTAU_H*TB_H))*OMIX(1,IH)/CB_H
     .+DREAL(CKTAU_H/(1.D0+CKTAU_H*TB_H))*OMIX(2,IH)/CB_H
     .+DIMAG(CKTAU_H*(TB_H**2+1.D0)/(1.D0+CKTAU_H*TB_H))*OMIX(3,IH)
     .            ,0.D0)
      NHC_H( 9,IH)=DCMPLX(
     .-DREAL((TB_H-CKTAU_H)/(1.D0+CKTAU_H*TB_H))*OMIX(3,IH)
     .+DIMAG(CKTAU_H*TB_H/(1.D0+CKTAU_H*TB_H))*OMIX(1,IH)/CB_H
     .-DIMAG(CKTAU_H/(1.D0+CKTAU_H*TB_H))*OMIX(2,IH)/CB_H
     .            ,0.D0)
*----
*      print*,'G^S_tau, G^P_tau : After including correction',
*     .        NHC_H(8,IH),NHC_H(9,IH)
*----
*|--->Higgs-down-down                                     [NC=10,11,12]
      NHC_H(10,IH)=DCMPLX(GW_H*MDMT_H/2.D0/MW_H,0.D0)
*      NHC_H(11,IH)=DCMPLX(OMIX(1,IH)/CB_H,0.D0)
*      NHC_H(12,IH)=DCMPLX(-OMIX(3,IH)*SB_H/CB_H,0.D0)
*JSL 10/Jun/2009: Including threshold corrections
      NHC_H(11,IH)=DCMPLX(
     . DREAL(1.D0/(1.D0+CKD_H*TB_H))*OMIX(1,IH)/CB_H
     .+DREAL(CKD_H/(1.D0+CKD_H*TB_H))*OMIX(2,IH)/CB_H
     .+DIMAG(CKD_H*(TB_H**2+1.D0)/(1.D0+CKD_H*TB_H))*OMIX(3,IH)
     .            ,0.D0)
      NHC_H(12,IH)=DCMPLX(
     .-DREAL((TB_H-CKD_H)/(1.D0+CKD_H*TB_H))*OMIX(3,IH)
     .+DIMAG(CKD_H*TB_H/(1.D0+CKD_H*TB_H))*OMIX(1,IH)/CB_H
     .-DIMAG(CKD_H/(1.D0+CKD_H*TB_H))*OMIX(2,IH)/CB_H
     .            ,0.D0)
*|--->Higgs-strange-strange                               [NC=13,14,15]
      NHC_H(13,IH)=DCMPLX(GW_H*MSMT_H/2.D0/MW_H,0.D0)
*      NHC_H(14,IH)=DCMPLX(OMIX(1,IH)/CB_H,0.D0)
*      NHC_H(15,IH)=DCMPLX(-OMIX(3,IH)*SB_H/CB_H,0.D0)
*       print*,'>>> FILLCOUPL2: cks_h',cks_h
*JSL 10/Jun/2009: Including threshold corrections
      NHC_H(14,IH)=DCMPLX(
     . DREAL(1.D0/(1.D0+CKS_H*TB_H))*OMIX(1,IH)/CB_H
     .+DREAL(CKS_H/(1.D0+CKS_H*TB_H))*OMIX(2,IH)/CB_H
     .+DIMAG(CKS_H*(TB_H**2+1.D0)/(1.D0+CKS_H*TB_H))*OMIX(3,IH)
     .            ,0.D0)
      NHC_H(15,IH)=DCMPLX(
     .-DREAL((TB_H-CKS_H)/(1.D0+CKS_H*TB_H))*OMIX(3,IH)
     .+DIMAG(CKS_H*TB_H/(1.D0+CKS_H*TB_H))*OMIX(1,IH)/CB_H
     .-DIMAG(CKS_H/(1.D0+CKS_H*TB_H))*OMIX(2,IH)/CB_H
     .            ,0.D0)
*|--->Higgs-bottom-bottom                                 [NC=16,17,18]
      NHC_H(16,IH)=DCMPLX(GW_H*MBMT_H/2.D0/MW_H,0.D0)
*-------------------------------------------------------------------------*
*---- JSL[21/JAN/2004] EW Radiative correction  to H-b-b            ------*
*----      Eq.(5) of hep-ph/0106027 by Guasch, Hollik, and Penaranda.-----*
*----      CKB_H = CKB_H (from SUBROUTINE RADNHTB) + CKBEW_H        ------*
*----              where CKBEW_H is the term propotional to alpha   ------*
*----              with CKB_H=Delta/tb, M -> M^* and mu -> mu^*     ------*
*-------------------------------------------------------------------------*
      IF(IFLAG_H(10).EQ.0) THEN
       CKBEW_H =DCONJG(MU_H)*AEM_H/4.D0/PI*
     . (
     .  -DCONJG(M2_H)/SW_H**2*
     .   (F_I(STMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *STMIX(1,1)*DCONJG(STMIX(1,1))
     .   +F_I(STMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *STMIX(1,2)*DCONJG(STMIX(1,2))
     .   +F_I(SBMASS(1)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(1,1)*DCONJG(SBMIX(1,1))/2.D0
     .   +F_I(SBMASS(2)**2,CDABS(M2_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(1,2)*DCONJG(SBMIX(1,2))/2.D0)  ! M2_H
     .  -DCONJG(M1_H)/CW_H**2/3.D0*
     .   (F_I(SBMASS(1)**2,SBMASS(2)**2,CDABS(M1_H)**2)/3.D0
     .   +F_I(SBMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(1,1)*DCONJG(SBMIX(1,1))/2.D0
     .   +F_I(SBMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(1,2)*DCONJG(SBMIX(1,2))/2.D0
     .   +F_I(SBMASS(1)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(2,1)*DCONJG(SBMIX(2,1))
     .   +F_I(SBMASS(2)**2,CDABS(M1_H)**2,CDABS(MU_H)**2)
     .    *SBMIX(2,2)*DCONJG(SBMIX(2,2)))       ! M1_H
     . )
      ELSE
       CKBEW_H =DCMPLX(0.D0,0.D0)
      ENDIF
*      print*,'CKB_H and CKBEW_H',CKB_H,CKBEW_H
      CKB_SUM=CKB_H+CKBEW_H
      CAUX_H(10)=CKB_SUM 
      CAUX_H(211)=CKB_SUM*TB_H/(1+CKB_SUM*TB_H)
      NHC_H(17,IH)=DCMPLX(
     . DREAL(1.D0/(1.D0+CKB_SUM*TB_H))*OMIX(1,IH)/CB_H
     .+DREAL(CKB_SUM/(1.D0+CKB_SUM*TB_H))*OMIX(2,IH)/CB_H
     .+DIMAG(CKB_SUM*(TB_H**2+1.D0)/(1.D0+CKB_SUM*TB_H))*OMIX(3,IH)
     .            ,0.D0)
      NHC_H(18,IH)=DCMPLX(
     .-DREAL((TB_H-CKB_SUM)/(1.D0+CKB_SUM*TB_H))*OMIX(3,IH)
     .+DIMAG(CKB_SUM*TB_H/(1.D0+CKB_SUM*TB_H))*OMIX(1,IH)/CB_H
     .-DIMAG(CKB_SUM/(1.D0+CKB_SUM*TB_H))*OMIX(2,IH)/CB_H
     .            ,0.D0)
*|--->Higgs-up-up                                         [NC=19,20,21]
      NHC_H(19,IH)=DCMPLX(GW_H*MUMT_H/2.D0/MW_H,0.D0)
      NHC_H(20,IH)=DCMPLX(OMIX(2,IH)/SB_H,0.D0)
      NHC_H(21,IH)=DCMPLX(-OMIX(3,IH)*CB_H/SB_H,0.D0)
*|--->Higgs-charm-charm                                   [NC=22,23,24]
      NHC_H(22,IH)=DCMPLX(GW_H*MCMT_H/2.D0/MW_H,0.D0)
      NHC_H(23,IH)=DCMPLX(OMIX(2,IH)/SB_H,0.D0)
      NHC_H(24,IH)=DCMPLX(-OMIX(3,IH)*CB_H/SB_H,0.D0)
*|--->Higgs-top-top                                       [NC=25,26,27]
      NHC_H(25,IH)=DCMPLX(GW_H*MTMT_H/2.D0/MW_H,0.D0)
      NHC_H(26,IH)=DCMPLX(
     . DREAL(1.D0/(1.D0+CKT_H/TB_H))*OMIX(2,IH)/SB_H
     .+DREAL(CKT_H/(1.D0+CKT_H/TB_H))*OMIX(1,IH)/SB_H
     .+DIMAG(CKT_H*(1.D0/TB_H**2+1.D0)/(1.D0+CKT_H/TB_H))*OMIX(3,IH)
     .            ,0.D0)
      NHC_H(27,IH)=DCMPLX(
     .-DREAL((1.D0/TB_H-CKT_H)/(1.D0+CKT_H/TB_H))*OMIX(3,IH)
     .+DIMAG(CKT_H/TB_H/(1.D0+CKT_H/TB_H))*OMIX(2,IH)/SB_H
     .-DIMAG(CKT_H/(1.D0+CKT_H/TB_H))*OMIX(1,IH)/SB_H
     .            ,0.D0)
*      DO NXX=0,8
*      NX=3*NXX+1
*      print*,nx,ih,ih,ih,nhc_h(nx,ih),nhc_h(nx+1,ih),nhc_h(nx+2,ih)
*      ENDDO
*|--->Higgs-neutralino-neutrino
      GP1 = OMIX(1,IH)-XI*SB_H*OMIX(3,IH)
      GP2 = OMIX(2,IH)-XI*CB_H*OMIX(3,IH)
*                I J
      CALL HNINJ(1,1,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(28,IH),NHC_H(29,IH),NHC_H(30,IH)) ! [NC=28,29,30]
      CALL HNINJ(2,2,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(31,IH),NHC_H(32,IH),NHC_H(33,IH)) ! [NC=31,32,33]
      CALL HNINJ(3,3,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(34,IH),NHC_H(35,IH),NHC_H(36,IH)) ! [NC=34,35,36]
      CALL HNINJ(4,4,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(37,IH),NHC_H(38,IH),NHC_H(39,IH)) ! [NC=37,38,39]
      CALL HNINJ(1,2,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(40,IH),NHC_H(41,IH),NHC_H(42,IH)) ! [NC=40,41,42]
      CALL HNINJ(1,3,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(43,IH),NHC_H(44,IH),NHC_H(45,IH)) ! [NC=43,44,45]
      CALL HNINJ(1,4,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(46,IH),NHC_H(47,IH),NHC_H(48,IH)) ! [NC=46,47,48]
      CALL HNINJ(2,3,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(49,IH),NHC_H(50,IH),NHC_H(51,IH)) ! [NC=49,50,51]
      CALL HNINJ(2,4,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(52,IH),NHC_H(53,IH),NHC_H(54,IH)) ! [NC=52,53,54]
      CALL HNINJ(3,4,TW_H,GP1,GP2,N_N,GW_H
     .          ,NHC_H(55,IH),NHC_H(56,IH),NHC_H(57,IH)) ! [NC=55,56,57]
*      do ni=0,9
*      nn=3*ni+28
*      print*,ih,nhc_h(nn,ih),nhc_h(nn+1,ih),nhc_h(nn+2,ih)
*      enddo
*|--->Higgs-chargino-chargino
*                I J
      CALL HCICJ(1,1,GP1,GP2,U_L,U_R,GW_H
     .          ,NHC_H(58,IH),NHC_H(59,IH),NHC_H(60,IH)) ! [NC=58,59,60]
      CALL HCICJ(1,2,GP1,GP2,U_L,U_R,GW_H
     .          ,NHC_H(61,IH),NHC_H(62,IH),NHC_H(63,IH)) ! [NC=61,62,63]
      CALL HCICJ(2,1,GP1,GP2,U_L,U_R,GW_H
     .          ,NHC_H(64,IH),NHC_H(65,IH),NHC_H(66,IH)) ! [NC=64,65,66]
      CALL HCICJ(2,2,GP1,GP2,U_L,U_R,GW_H
     .          ,NHC_H(67,IH),NHC_H(68,IH),NHC_H(69,IH)) ! [NC=67,68,69]
*      do ni=0,3
*      nc=3*ni+58
*      print*,ih,nhc_h(nc,ih),nhc_h(nc+1,ih),nhc_h(nc+2,ih)
*      enddo
*|--->Higgs-Vector boson-Vector boson                      [NC=70]
      NHC_H(70,IH)=DCMPLX(CB_H*OMIX(1,IH)+SB_H*OMIX(2,IH),0.D0)
*      print*,ih,nhc_h(70,ih)
*|--->Higgs-stop-stop 
      CALL HSTST(IH,1,1,OMIX,HT_H,STMIX,NHC_H(71,IH))   !  [NC=71]
      CALL HSTST(IH,1,2,OMIX,HT_H,STMIX,NHC_H(72,IH))   !  [NC=72]
      CALL HSTST(IH,2,1,OMIX,HT_H,STMIX,NHC_H(73,IH))   !  [NC=73]
      CALL HSTST(IH,2,2,OMIX,HT_H,STMIX,NHC_H(74,IH))   !  [NC=74]
*      print*,ih,nhc_h(71,ih),nhc_h(72,ih),nhc_h(73,ih),nhc_h(74,ih)
*|--->Higgs-sbottom-sbottom 
      CALL HSBSB(IH,1,1,OMIX,HB_H,SBMIX,NHC_H(75,IH))   !  [NC=75]
      CALL HSBSB(IH,1,2,OMIX,HB_H,SBMIX,NHC_H(76,IH))   !  [NC=76]
      CALL HSBSB(IH,2,1,OMIX,HB_H,SBMIX,NHC_H(77,IH))   !  [NC=77]
      CALL HSBSB(IH,2,2,OMIX,HB_H,SBMIX,NHC_H(78,IH))   !  [NC=78]
*      print*,ih,nhc_h(75,ih),nhc_h(76,ih),nhc_h(77,ih),nhc_h(78,ih)
*|--->Higgs-stau-stau 
      CALL HSTUSTU(IH,1,1,OMIX,HTAU_H,STAUMIX,NHC_H(79,IH)) !  [NC=79]
      CALL HSTUSTU(IH,1,2,OMIX,HTAU_H,STAUMIX,NHC_H(80,IH)) !  [NC=80]
      CALL HSTUSTU(IH,2,1,OMIX,HTAU_H,STAUMIX,NHC_H(81,IH)) !  [NC=81]
      CALL HSTUSTU(IH,2,2,OMIX,HTAU_H,STAUMIX,NHC_H(82,IH)) !  [NC=82]
*|--->Higgs-snu3-snu3 
      NHC_H(83,IH)=GW_H*MZ_H/2.D0/CW_H/V_H                  !  [NC=83]
     .            *(-OMIX(1,IH)*CB_H+OMIX(2,IH)*SB_H)
*|--->Higgs-glue-glue                                      [NC=84,85]
      CALL HGG(IH,HMASS,STMASS,SBMASS,NCMAX,NHC_H
     .        ,NHC_H(84,IH),NHC_H(85,IH),SGG,PGG)
*JSL 10/Jun/2009: Save quark running masses in RAUX_H(501-509)
      IF(IH.EQ.1) THEN
       RAUX_H(501)=MB_MH
       RAUX_H(502)=MT_MH
       RAUX_H(503)=MC_MH
      ENDIF
      IF(IH.EQ.2) THEN
       RAUX_H(504)=MB_MH
       RAUX_H(505)=MT_MH
       RAUX_H(506)=MC_MH
      ENDIF
      IF(IH.EQ.3) THEN
       RAUX_H(507)=MB_MH
       RAUX_H(508)=MT_MH
       RAUX_H(509)=MC_MH
      ENDIF
*
      CALL HGGSM(HMASS(IH),SGGSM,PGGSM)
      IF(IH.EQ.1) THEN
       CAUX_H(221)=SGGSM(3)
      ELSEIF(IH.EQ.2) THEN
       CAUX_H(222)=SGGSM(3)
      ELSEIF(IH.EQ.3) THEN
       CAUX_H(223)=SGGSM(3)
      ELSE
       print*,'IH is out of range, IH=',IH
       RETURN
      ENDIF
*      print*,'---------------------------------------------------------'
*      print*,'H-G-G with IH = ',ih
*      print*,'---------------------------------------------------------'
*      print*,'|--> MSSM contributions :'
*      print*,'bottom,top    ',sgg(1),sgg(2)
*      print*,'stop11+stop22 ',sgg(3)+sgg(4)
*      print*,'sbot11+sbot22 ',sgg(5)+sgg(6)
*      print*,'S(H-g-g)[MSSM]',nhc_h(84,ih)
*      print*,'zero?         ',nhc_h(84,ih)-sgg(7)
*      print*,'bottom,top    ',pgg(1),pgg(2)
*      print*,'P(H-g-g)[MSSM]',nhc_h(85,ih)
*      print*,'zero?         ',nhc_h(85,ih)-pgg(3)
*      print*,'|--> SM   contributions :'
*      print*,'bottom,top    ',sggsm(1),sggsm(2)
*      print*,'S(H-g-g)  [SM]',sggsm(3)
**      print*,'bottom,top    ',pgg(1),pgg(2)
**      print*,'zeros?    [SM]',pggsm(3)
*      print*,'---------------------------------------------------------'
*|--->Higgs-ChargedHiggs-ChargedHiggs                      [NC=86]
      CALL HCHCH(IH,CB_H,SB_H,LR_H,LC_H,OMIX,NHC_H(86,IH))
*      print*,ih,nhc_h(86,ih)
*|--->Higgs-ChargedHiggs^(+)-W^(-)                         [NC=87]
      NHC_H(87,IH)=DCMPLX(CB_H*OMIX(2,IH)-SB_H*OMIX(1,IH),-OMIX(3,IH))
*      print*,ih,abs(nhc_h(70,ih))**2+abs(nhc_h(87,ih))**2
*      print*,ih,nhc_h(87,ih)
*|--->Higgs-photon-photon                                  [NC=88,89]
      CALL HPP(IH,ASMH2,MTMH2,MBMH2,MCMH2
     .        ,M_C,MCH,HMASS,STMASS,SBMASS,STAUMASS,NCMAX,NHC_H
     .        ,NHC_H(88,IH),NHC_H(89,IH),SPP,PPP)
      CALL HPPSM(HMASS(IH),MBMH2,MTMH2,MCMH2,SPPSM,PPPSM)
      IF(IH.EQ.1) THEN
       CAUX_H(231)=SPPSM(9)
      ELSEIF(IH.EQ.2) THEN
       CAUX_H(232)=SPPSM(9)
      ELSEIF(IH.EQ.3) THEN
       CAUX_H(233)=SPPSM(9)
      ELSE
       print*,'IH is out of range, IH=',IH
       RETURN
      ENDIF
*      print*,'---> FILLCOUPL '
*      print*,'H-P-P with IH = ',ih
*      print*,'bottom,top     ',spp(1),spp(2)
*      print*,'charm,tau      ',spp(3),spp(4)
*      print*,'c.ino11,c.ino22',spp(5),spp(6)
*      print*,'stop11,stop22  ',spp(7),spp(8)
*      print*,'sbot11,sbot22  ',spp(9),spp(10)
*      print*,'ww,chch        ',spp(11),spp(12)
*      print*,'stau11,stau22  ',spp(13),spp(14)
*      print*,'S(H-p-p)       ',nhc_h(88,ih)
*      print*,'zero?          ',nhc_h(88,ih)-spp(15)
*      print*,'bottom,top     ',ppp(1),ppp(2)
*      print*,'charm,tau      ',ppp(3),ppp(4)
*      print*,'c.ino11,c.ino22',ppp(5),ppp(6)
*      print*,'P(H-p-p)       ',nhc_h(89,ih)
*      print*,'zero?',nhc_h(89,ih)-ppp(7)
*|--->Higgs-glue-glue in the limit of vanishing HMASS(I)   [NC=90,91]
       NHC_H(90,IH)=NHC_H(17,IH)*2.D0/3.D0
     .             +NHC_H(26,IH)*2.D0/3.D0
     .             -NHC_H(71,IH)*V_H**2/4.D0/STMASS(1)**2/3.D0
     .             -NHC_H(74,IH)*V_H**2/4.D0/STMASS(2)**2/3.D0
     .             -NHC_H(75,IH)*V_H**2/4.D0/SBMASS(1)**2/3.D0
     .             -NHC_H(78,IH)*V_H**2/4.D0/SBMASS(2)**2/3.D0
       NHC_H(91,IH)=NHC_H(18,IH)+NHC_H(27,IH)
*      print*,nhc_h(90,ih),nhc_h(91,ih)
*|--->Higgs-gamma-gamma in the limit of vanishing HMASS(I) [NC=92,93]
       NHC_H(92,IH)=2.D0/3.D0*NHC_H(17,IH)*2.D0/3.D0
     .             +8.D0/3.D0*NHC_H(26,IH)*2.D0/3.D0
     .             +2.D0*NHC_H(58,IH)*NHC_H(59,IH)*V_H/M_C(1)*2.D0/3.D0
     .             +2.D0*NHC_H(67,IH)*NHC_H(68,IH)*V_H/M_C(2)*2.D0/3.D0
     .             -8.D0/3.D0*NHC_H(71,IH)*V_H**2/4.D0/STMASS(1)**2/3.D0
     .             -8.D0/3.D0*NHC_H(74,IH)*V_H**2/4.D0/STMASS(2)**2/3.D0
     .             -2.D0/3.D0*NHC_H(75,IH)*V_H**2/4.D0/SBMASS(1)**2/3.D0
     .             -2.D0/3.D0*NHC_H(78,IH)*V_H**2/4.D0/SBMASS(2)**2/3.D0
     .             -NHC_H(70,IH)*7.D0
     .             -NHC_H(86,IH)*V_H**2/2.D0/MCH**2/3.D0
     .             -2.D0*NHC_H(79,IH)*V_H**2/4.D0/STAUMASS(1)**2/3.D0
     .             -2.D0*NHC_H(82,IH)*V_H**2/4.D0/STAUMASS(2)**2/3.D0
       NHC_H(93,IH)=2.D0/3.D0*NHC_H(18,IH)+8.D0/3.D0*NHC_H(27,IH)
     .             +2.D0*NHC_H(58,IH)*NHC_H(60,IH)*V_H/M_C(1)
     .             +2.D0*NHC_H(67,IH)*NHC_H(69,IH)*V_H/M_C(2)
*      print*,nhc_h(92,ih),nhc_h(93,ih)
      ENDDO ! DO IH=1,3
*-----------------------------------------------------------------------
*===> Neutral Higgs self coupling             : SHC_H(NC)         <====*
*-----------------------------------------------------------------------
*                  I J K : HI-HJ-HK Coupling with I>=J>=K
      CALL NHSELF3(3,3,3,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 1))   ! NC= 1
      CALL NHSELF3(3,3,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 2))   ! NC= 2
      CALL NHSELF3(3,3,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 3))   ! NC= 3
      CALL NHSELF3(3,2,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 4))   ! NC= 4
      CALL NHSELF3(3,2,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 5))   ! NC= 5
      CALL NHSELF3(3,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 6))   ! NC= 6
      CALL NHSELF3(2,2,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 7))   ! NC= 7
      CALL NHSELF3(2,2,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 8))   ! NC= 8
      CALL NHSELF3(2,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H( 9))   ! NC= 9
      CALL NHSELF3(1,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(10))   ! NC=10
      SHC_H(11)=DREAL(NHC_H(86,1))                             ! NC=11
      SHC_H(12)=DREAL(NHC_H(86,2))                             ! NC=12
      SHC_H(13)=DREAL(NHC_H(86,3))                             ! NC=13
*To compare with old symmetrized conventions
*      print*,'H333   ',shc_h(1)
*      print*,'H332/3 ',shc_h(2)/3.d0
*      print*,'H331/3 ',shc_h(3)/3.d0
*      print*,'H322/3 ',shc_h(4)/3.d0
*      print*,'H321/6 ',shc_h(5)/6.d0
*      print*,'H311/3 ',shc_h(6)/3.d0
*      print*,'H222   ',shc_h(7)
*      print*,'H221/3 ',shc_h(8)/3.d0
*      print*,'H211/3 ',shc_h(9)/3.d0
*-----
*      print*,'H111   ',shc_h(10)
*      print*,'H333   ',shc_h(1)
*      print*,'H332   ',shc_h(2)
*      print*,'H331   ',shc_h(3)
*      print*,'H322   ',shc_h(4)
*      print*,'H321   ',shc_h(5)
*      print*,'H311   ',shc_h(6)
*      print*,'H222   ',shc_h(7)
*      print*,'H221   ',shc_h(8)
*      print*,'H211   ',shc_h(9)
*      print*,'H111   ',shc_h(10)
*                  I J K L: HI-HJ-HK-HL Coupling with I>=J>=K>=L
      CALL NHSELF4(3,3,3,3,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(14)) ! NC=14
      CALL NHSELF4(3,3,3,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(15)) ! NC=15
      CALL NHSELF4(3,3,3,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(16)) ! NC=16
      CALL NHSELF4(3,3,2,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(17)) ! NC=17
      CALL NHSELF4(3,3,2,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(18)) ! NC=18
      CALL NHSELF4(3,3,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(19)) ! NC=19
      CALL NHSELF4(3,2,2,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(20)) ! NC=10
      CALL NHSELF4(3,2,2,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(21)) ! NC=21
      CALL NHSELF4(3,2,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(22)) ! NC=22
      CALL NHSELF4(3,1,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(23)) ! NC=23
      CALL NHSELF4(2,2,2,2,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(24)) ! NC=24
      CALL NHSELF4(2,2,2,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(25)) ! NC=25
      CALL NHSELF4(2,2,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(26)) ! NC=26
      CALL NHSELF4(2,1,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(27)) ! NC=27
      CALL NHSELF4(1,1,1,1,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(28)) ! NC=28
      CALL NHSELF4(3,3,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(29)) ! NC=29
      CALL NHSELF4(3,2,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(30)) ! NC=30
      CALL NHSELF4(3,1,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(31)) ! NC=31
      CALL NHSELF4(2,2,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(32)) ! NC=32
      CALL NHSELF4(2,1,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(33)) ! NC=33
      CALL NHSELF4(1,1,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(34)) ! NC=34
      CALL NHSELF4(5,5,5,5,SB_H,CB_H,OMIX,LR_H,LC_H,SHC_H(35)) ! NC=35
*-----------------------------------------------------------------------
*==> Positively Charged Higgs coupling to two particles : CHC_H(NC) <==*
*-----------------------------------------------------------------------
*Charged Higgs(+)-electron-neutrino                       NC=[ 1, 2, 3]
      CHC_H(1)=DCMPLX(-GW_H*ME_H/DSQRT(2.D0)/MW_H,0.D0)
      CHC_H(2)=DCMPLX(TB_H/2.D0,0.D0)
      CHC_H(3)=DCMPLX(0.D0,-TB_H/2.D0)
*Charged Higgs(+)-muon-neutrino                           NC=[ 4, 5, 6]
      CHC_H(4)=DCMPLX(-GW_H*MMU_H/DSQRT(2.D0)/MW_H,0.D0)
      CHC_H(5)=DCMPLX(TB_H/2.D0,0.D0)
      CHC_H(6)=DCMPLX(0.D0,-TB_H/2.D0)
*Charged Higgs(+)-tau-neutrino                            NC=[ 7, 8, 9]
      CHC_H(7)=DCMPLX(-GW_H*MTAU_H/DSQRT(2.D0)/MW_H,0.D0)
      CHC_H(8)=DCMPLX(TB_H/2.D0,0.D0)
      CHC_H(9)=DCMPLX(0.D0,-TB_H/2.D0)
*Charged Higgs(+)-u quark-d quark                         NC=[10,11,12]
      CHC_H(10)=DCMPLX(-GW_H*MUMT_H/DSQRT(2.D0)/MW_H,0.D0)
*      CHC_H(11)=DCMPLX((1.D0/TB_H+MDMT_H/MUMT_H*TB_H)/2.D0,0.D0)
*      CHC_H(12)=DCMPLX(0.D0,(1.D0/TB_H-MDMT_H/MUMT_H*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections 
      CHC_H(11)=(1.D0/TB_H
     .          +TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMT_H/MUMT_H)/2.D0
      CHC_H(12)=(1.D0/TB_H
     .          -TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMT_H/MUMT_H)*XI/2.D0
*Charged Higgs(+)-c quark-s quark                         NC=[13,14,15]
      CHC_H(13)=DCMPLX(-GW_H*MCMT_H/DSQRT(2.D0)/MW_H,0.D0)
*      CHC_H(14)=DCMPLX((1.D0/TB_H+MSMT_H/MCMT_H*TB_H)/2.D0,0.D0)
*      CHC_H(15)=DCMPLX(0.D0,(1.D0/TB_H-MSMT_H/MCMT_H*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections 
*      print*,'>>> FILLCOUPL2: cks_h',cks_h
      CHC_H(14)=(1.D0/TB_H
     .          +TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMT_H/MCMT_H)/2.D0
      CHC_H(15)=(1.D0/TB_H
     .          -TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMT_H/MCMT_H)*XI/2.D0
*       print*,'FILLCOUPL2',msmt_h,mcmt_h,msmt_h/mcmt_h
*Charged Higgs(+)-t quark-b quark                         NC=[16,17,18]
      CHC_H(16)=DCMPLX(-GW_H*MTMT_H/DSQRT(2.D0)/MW_H,0.D0)
      CHC_H(17)=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMT_H/MTMT_H)/2.D0
      CHC_H(18)=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMT_H/MTMT_H)*XI/2.D0
*      print*,chc_h(1),chc_h(2),chc_h(3)
*      print*,chc_h(4),chc_h(5),chc_h(6)
*      print*,chc_h(7),chc_h(8),chc_h(9)
*      print*,chc_h(10),chc_h(11),chc_h(12)
*      print*,chc_h(13),chc_h(14),chc_h(15)
*      print*,chc_h(16),chc_h(17),chc_h(18)
*
*Charged Higgs(+)-neutralinoI-charginoJ   
*                I J                 
      CALL HNICJ(1,1,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(19),CHC_H(20),CHC_H(21))         ! [NC=19,20,21]
      CALL HNICJ(1,2,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(22),CHC_H(23),CHC_H(24))         ! [NC=22,23,24]
      CALL HNICJ(2,1,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(25),CHC_H(26),CHC_H(27))         ! [NC=25,26,27]
      CALL HNICJ(2,2,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(28),CHC_H(29),CHC_H(30))         ! [NC=28,29,30]
      CALL HNICJ(3,1,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(31),CHC_H(32),CHC_H(33))         ! [NC=31,32,33]
      CALL HNICJ(3,2,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(34),CHC_H(35),CHC_H(36))         ! [NC=34,35,36]
      CALL HNICJ(4,1,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(37),CHC_H(38),CHC_H(39))         ! [NC=37,38,39]
      CALL HNICJ(4,2,N_N,U_L,U_R,SB_H,CB_H,TW_H,GW_H
     .          ,CHC_H(40),CHC_H(41),CHC_H(42))         ! [NC=40,41,42]
*      do ii=0,7
*      nc=3*ii+19
*      print*,chc_h(nc),chc_h(nc+1),chc_h(nc+2)
*      enddo
*
*Charged Higgs(+)-Stop*-Sbottom   
*
      CALL HSTSB(1,1,HT_H,HB_H,STMIX,SBMIX,CHC_H(43))   !  [NC=43]
      CALL HSTSB(1,2,HT_H,HB_H,STMIX,SBMIX,CHC_H(44))   !  [NC=44]
      CALL HSTSB(2,1,HT_H,HB_H,STMIX,SBMIX,CHC_H(45))   !  [NC=45]
      CALL HSTSB(2,2,HT_H,HB_H,STMIX,SBMIX,CHC_H(46))   !  [NC=46]
*      do nc=43,46
*      print*,chc_h(nc)
*      enddo
*Charged Higgs(+)-Snu3*-Stau                               [NC=47,48]
      CHN3TUL=(CDABS(HTAU_H)**2-GW_H**2)*V_H*SB_H*CB_H/DSQRT(2.D0)
      CHN3TUR=DCONJG(HTAU_H)*(SB_H*DCONJG(ATAU_H)+CB_H*MU_H)
      CHC_H(47)=(CHN3TUL*STAUMIX(1,1)+CHN3TUR*STAUMIX(2,1))/V_H
      CHC_H(48)=(CHN3TUL*STAUMIX(1,2)+CHN3TUR*STAUMIX(2,2))/V_H
   
*-----------------------------------------------------------------------
      IF(IFLAG_H(5).EQ.1) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,1)
      IF(IFLAG_H(5).EQ.2) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,2)
      IF(IFLAG_H(5).EQ.3) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,3)
      IF(IFLAG_H(5).EQ.4) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,4)
      IF(IFLAG_H(5).EQ.5) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,5)
      IF(IFLAG_H(5).EQ.6) CALL DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,6)
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE RADCHTB(NFLAG,IFLAG_H,SBMIX,SBMASS,STMIX,STMASS
     .                  ,RB,RT,CKBB,CKBT)
************************************************************************
*
* Calculation of the finite radiatve corrections to H-t-t and
* H-t-b Yukawa couplings
*
* Ref: 'Collider Probes of the MSSM Higgs Sector with Explicit
*       CP Violation', hep-ph/0211467 by
*       Carena, Ellis, Mrenna, Pilaftsis, and Wagner
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
*-----------------------------------------------------------------------
*Input Array
      INTEGER*8 IFLAG_H(NFLAG)
      COMPLEX*16 SBMIX(2,2),STMIX(2,2)
      REAL*8 SBMASS(2),STMASS(2)
*Output
      COMPLEX*16 RB,RT,CKBB,CKBT
*Local
      COMPLEX*16 DSB,DLB,DST,DLT
      COMPLEX*16 DSBN,DSTN
*Stop scale
      QT2  = DMAX1(MQ3_H**2+MTPOLE_H**2,MU3_H**2+MTPOLE_H**2)
*Sbottom scale
      QB2  = DMAX1(MQ3_H**2+MBMT_H**2,MD3_H**2+MBMT_H**2)
*Stop-Sbottom scale
      QTB2=DMAX1(QT2,QB2)
*alpha_S at sferimon mass scale
      PI=2.D0*DASIN(1.D0)
      B6=(11.D0-2.D0*6.D0/3.D0)/4.D0/PI
      AS=ASMT_H/(1.D0+B6*ASMT_H*DLOG(QTB2/MTPOLE_H**2))
      ASMST=ASMT_H/(1.D0+B6*ASMT_H*DLOG(QT2/MTPOLE_H**2))
      ASMSB=ASMT_H/(1.D0+B6*ASMT_H*DLOG(QB2/MTPOLE_H**2))
*tree-level Yukawa couplings at sferimon mass scale
      HBT=DSQRT(2.D0)*MBMT_H/V_H/CB_H
      HTT=DSQRT(2.D0)*MTMT_H/V_H/SB_H
      BHB=(9.D0*HBT**2/2.D0+HTT**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
      BHT=(9.D0*HTT**2/2.D0+HBT**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
      HBR=HBT*(1.D0+2.D0*BHB*DLOG(QTB2/MTPOLE_H**2))**0.25D0
      HTR=HTT*(1.D0+2.D0*BHT*DLOG(QTB2/MTPOLE_H**2))**0.25D0
*
      HBRN=HBT*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
      HTRN=HTT*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
*      print*,asmt_h,as
*      print*,hbt,hbr
*      print*,htt,htr
      CT2=CDABS(STMIX(1,2))**2
      ST2=CDABS(STMIX(1,1))**2
*      print*,ct2,st2
      CB2=CDABS(SBMIX(1,2))**2
      SB2=CDABS(SBMIX(1,1))**2
*      print*,cb2,sb2
*
      IF(IFLAG_H(10).EQ.0) THEN
       FI1M3=CT2*SB2*F_I(STMASS(2)**2,SBMASS(2)**2,CDABS(M3_H)**2)
     .      +ST2*CB2*F_I(STMASS(1)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .      +CT2*CB2*F_I(STMASS(2)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .      +ST2*SB2*F_I(STMASS(1)**2,SBMASS(2)**2,CDABS(M3_H)**2)
       FI2M3=ST2*CB2*F_I(STMASS(2)**2,SBMASS(2)**2,CDABS(M3_H)**2)
     .      +CT2*SB2*F_I(STMASS(1)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .      +ST2*SB2*F_I(STMASS(2)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .      +CT2*CB2*F_I(STMASS(1)**2,SBMASS(2)**2,CDABS(M3_H)**2)
       FI1MU=CT2*SB2*F_I(STMASS(2)**2,SBMASS(2)**2,CDABS(MU_H)**2)
     .      +ST2*CB2*F_I(STMASS(1)**2,SBMASS(1)**2,CDABS(MU_H)**2)
     .      +CT2*CB2*F_I(STMASS(2)**2,SBMASS(1)**2,CDABS(MU_H)**2)
     .      +ST2*SB2*F_I(STMASS(1)**2,SBMASS(2)**2,CDABS(MU_H)**2)
       FI2MU=ST2*CB2*F_I(STMASS(2)**2,SBMASS(2)**2,CDABS(MU_H)**2)
     .      +CT2*SB2*F_I(STMASS(1)**2,SBMASS(1)**2,CDABS(MU_H)**2)
     .      +ST2*SB2*F_I(STMASS(2)**2,SBMASS(1)**2,CDABS(MU_H)**2)
     .      +CT2*CB2*F_I(STMASS(1)**2,SBMASS(2)**2,CDABS(MU_H)**2)
*       print*,fi1m3,fi2m3
*       print*,fi1mu,fi2mu
       DSB=-2.D0*AS/3.D0/PI*DCONJG(M3_H)*AB_H*FI1M3
     .     -HTR**2/16.D0/PI**2*CDABS(MU_H)**2*FI2MU
       DLB= 2.D0*AS/3.D0/PI*DCONJG(M3_H*MU_H)*FI1M3
     .     +HTR**2/16.D0/PI**2*DCONJG(AT_H*MU_H)*FI2MU
       DST=-2.D0*AS/3.D0/PI*DCONJG(M3_H)*AT_H*FI2M3
     .     -HBR**2/16.D0/PI**2*CDABS(MU_H)**2*FI1MU
       DLT= 2.D0*AS/3.D0/PI*DCONJG(M3_H*MU_H)*FI2M3
     .     +HBR**2/16.D0/PI**2*DCONJG(AB_H*MU_H)*FI1MU
       DSBN=-2.D0*ASMSB/3.D0/PI*DCONJG(M3_H)*AB_H
     .      *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .      -HTRN**2/16.D0/PI**2*CDABS(MU_H)**2
     .      *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(MU_H)**2)
       DSTN=-2.D0*ASMST/3.D0/PI*DCONJG(M3_H)*AT_H
     .      *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(M3_H)**2)
     .      -HBRN**2/16.D0/PI**2*CDABS(MU_H)**2
     .      *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(MU_H)**2)
       RB=(DSB-DSBN)/(1.D0+DSBN)
       RT=(DST-DSTN)/(1.D0+DSTN)
       CKBB=DLB/(1.D0+DSBN)
       CKBT=DLT/(1.D0+DSTN)
      ELSE
       RB=DCMPLX(0.D0,0.D0)
       RT=DCMPLX(0.D0,0.D0)
       CKBB=DCMPLX(0.D0,0.D0)
       CKBT=DCMPLX(0.D0,0.D0)
      ENDIF
*
      RETURN
      END



      SUBROUTINE RADNHTB(NFLAG,IFLAG_H,SBMASS,STMASS,HB,HT,CKB,CKT)
************************************************************************
*
* Calculation of the finite radiatve corrections to H-t-t and 
* H-b-b Yukawa couplings
*
* Ref: 'Collider Probes of the MSSM Higgs Sector with Explicit
*       CP Violation', hep-ph/0211467 by
*       Carena, Ellis, Mrenna, Pilaftsis, and Wagner
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
*-----------------------------------------------------------------------
*Input Array
      INTEGER*8 IFLAG_H(NFLAG) 
*Output 
      REAL*8 SBMASS(2),STMASS(2)
      COMPLEX*16 HB,HT,CKB,CKT
*Local 
      COMPLEX*16 HB0,HT0
      COMPLEX*16 DSB,DLB,DST,DLT
      COMPLEX*16 CKBOLD
      COMPLEX*16 SBMIX(2,2),STMIX(2,2) ! SQMIX is to be called
*Stop scale
      QT2  = DMAX1(MQ3_H**2+MTPOLE_H**2,MU3_H**2+MTPOLE_H**2)
*Sbottom scale
      QB2  = DMAX1(MQ3_H**2+MBMT_H**2,MD3_H**2+MBMT_H**2)
*alpha_S at sferimon mass scale
      PI=2.D0*DASIN(1.D0)
      B6=(11.D0-2.D0*6.D0/3.D0)/4.D0/PI
      ASMST=ASMT_H/(1.D0+B6*ASMT_H*DLOG(QT2/MTPOLE_H**2))
      ASMSB=ASMT_H/(1.D0+B6*ASMT_H*DLOG(QB2/MTPOLE_H**2))
*tree-level Yukawa couplings at sferimon mass scale
      HBT=DSQRT(2.D0)*MBMT_H/V_H/CB_H
      HTT=DSQRT(2.D0)*MTMT_H/V_H/SB_H
      BHB=(9.D0*HBT**2/2.D0+HTT**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
      BHT=(9.D0*HTT**2/2.D0+HBT**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
      HBR=HBT*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
      HTR=HTT*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
*      print*,hbt,htt
*      print*,hbr,htr
*

*Sfermion masses with tree level coupling
      HB0=DCMPLX(HBT,0.D0)
      HT0=DCMPLX(HTT,0.D0)
      CALL SQMIX(HB0,HT0,STMASS,SBMASS,STMIX,SBMIX)
*
      IF(IFLAG_H(10).EQ.0) THEN
*before entering iteration, if the tree level HB gives the tachyonic sbottom,
       IF(SBMASS(1).LT.0.D0) THEN
 997    HB0=0.9D0*HB0
        CALL SQMIX(HB0,HT0,STMASS,SBMASS,STMIX,SBMIX)
        IF(SBMASS(1).LT.0.D0) GOTO 997
       ENDIF
*
       CKBOLD=DCMPLX(0.D0,0.D0)
       MSB1OLD=SBMASS(1)
       ICOUNT=0
 998   CONTINUE
       DSB=-2.D0*ASMSB/3.D0/PI*DCONJG(M3_H)*AB_H
     .     *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .     -HTR**2/16.D0/PI**2*CDABS(MU_H)**2
     .     *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(MU_H)**2)
       DLB= 2.D0*ASMSB/3.D0/PI*DCONJG(M3_H*MU_H)
     .     *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(M3_H)**2)
     .     +HTR**2/16.D0/PI**2*DCONJG(AT_H*MU_H)
     .     *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(MU_H)**2)
       DST=-2.D0*ASMST/3.D0/PI*DCONJG(M3_H)*AT_H
     .     *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(M3_H)**2)
     .     -HBR**2/16.D0/PI**2*CDABS(MU_H)**2
     .     *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(MU_H)**2)
       DLT= 2.D0*ASMST/3.D0/PI*DCONJG(M3_H*MU_H)
     .     *F_I(STMASS(2)**2,STMASS(1)**2,CDABS(M3_H)**2)
     .     +HBR**2/16.D0/PI**2*DCONJG(AB_H*MU_H)
     .     *F_I(SBMASS(2)**2,SBMASS(1)**2,CDABS(MU_H)**2)
       HB=DSQRT(2.D0)*MBMT_H/V_H/CB_H/(1.D0+DSB+DLB*TB_H)
       HT=DSQRT(2.D0)*MTMT_H/V_H/SB_H/(1.D0+DST+DLT/TB_H)
       CKB=DLB/(1.D0+DSB)
       CKT=DLT/(1.D0+DST)
       CALL SQMIX(HB,HT,STMASS,SBMASS,STMIX,SBMIX)
        IF(SBMASS(1).LT.1.D0) THEN
          SBMASS(1)=0.1D0
        ENDIF
*--iteration for sfermion masses
       IF(CDABS(1.D0-CKBOLD/CKB).GT.0.0001D0) THEN
*       IF(DABS(SBMASS(1)-MSB1OLD).GT.0.1D0) THEN
**
*        print*,'RADNHTB : MSB(1)_OLD, MSB(1), |KB| : ',
*     .  ICOUNT,MSB1OLD,SBMASS(1),CDABS(CKB)
**
        CKBOLD=CKB
        MSB1OLD=SBMASS(1)
*        BHB=(9.D0*CDABS(HB)**2/2.D0+CDABS(HT)**2/2.D0
*     .       -32.D0*PI*ASMT_H)/16.D0/PI**2
*        BHT=(9.D0*CDABS(HT)**2/2.D0+CDABS(HB)**2/2.D0
*     .       -32.D0*PI*ASMT_H)/16.D0/PI**2
        BHB=(9.D0*DABS(HBT)**2/2.D0+DABS(HTT)**2/2.D0
     .       -32.D0*PI*ASMT_H)/16.D0/PI**2
        BHT=(9.D0*DABS(HTT)**2/2.D0+DABS(HBT)**2/2.D0
     .       -32.D0*PI*ASMT_H)/16.D0/PI**2
        HBR=CDABS(HB)*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
        HTR=CDABS(HT)*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.100) THEN
*          print*,'ERROR : <fillcoupl.f:RADNHTB> Iteration Failed '
          IFLAG_H(54)=1
          RETURN
        ENDIF
        GOTO 998 ! iterate
       ELSE
**
*        print*,'RADNHTB : OUT OF ITERATION'
*        CALL SQMIX(HB,HT,STMASS,SBMASS,STMIX,SBMIX)
*        print*,'RADNHTB : MSB(1)_OLD, MSB(1), |KB| : ',
*     .  ICOUNT,MSB1OLD,SBMASS(1),CDABS(CKB)
**
        GOTO 999 ! out of iteration
       ENDIF
 999   CONTINUE
*      print*,'>> FILLCOUPL:',dsb,dlb,dst,dlt
*      print*,'>> FILLCOUPL:',hbr,htr
*      print*,'>> FILLCOUPL:',cdabs(hb),cdabs(ht)
*      print*,'>> FILLCOUPL:',bhb,bht
*      print*,'>> FILLCOUPL:',qb2,qt2
*--iteration
      ELSE
       HB=DCMPLX(DSQRT(2.D0)*MBMT_H/V_H/CB_H,0.D0)
       HT=DCMPLX(DSQRT(2.D0)*MTMT_H/V_H/SB_H,0.D0)
       CKB=DCMPLX(0.D0,0.D0)
       CKT=DCMPLX(0.D0,0.D0)
      ENDIF
*
      RETURN
      END

      REAL*8 FUNCTION F_I(A,B,C)
************************************************************************
*JSL 18/APR/2005, A=B, B=C, C=A, A=B=C cases considered.
*JSL:28/AUG/2006 : improved treatment for the degenerate cases 
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      EPS=1.D-6
      IF(DABS((B-A)/A).LT.EPS) B=A
      IF(DABS((C-A)/A).LT.EPS) C=A
      IF(DABS((C-B)/B).LT.EPS) B=C
*
      IF(A.EQ.B .AND. B.EQ.C .AND. C.EQ.A) THEN
       F_I=1.D0/2.D0/A
      ELSEIF(A.EQ.B) THEN
       F_I=(B-C+C*DLOG(C/B))/(B-C)**2
      ELSEIF(B.EQ.C) THEN
       F_I=(C-A+A*DLOG(A/C))/(C-A)**2
      ELSEIF(C.EQ.A) THEN
       F_I=(A-B+B*DLOG(B/A))/(A-B)**2
      ELSE
       F_I=(A*B*DLOG(A/B)+B*C*DLOG(B/C)+C*A*DLOG(C/A))
     .    /((A-B)*(B-C)*(A-C))
      ENDIF
*
      RETURN
      END

      SUBROUTINE DUMP_COUPLING(NCMAX,NHC_H,SHC_H,CHC_H,ICPRI)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*Input Array
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 NHC_H(NCMAX,3),CHC_H(NCMAX)
      INTEGER    ICPRI
*
      IF(ICPRI.EQ.6.OR.ICPRI.LE.3) THEN
      IF(ICPRI.LE.3) IH=ICPRI
      IC = 0
 999  CONTINUE
      IC = IC+1
      IF(ICPRI.EQ.6) IH=IC
      print*,'---------------------------------------------------------'
      IF(IH.EQ.1)
     .print*,'The Lightest Higgs H_1 Couplings : NHC_H(NC,1)'
      IF(IH.EQ.2)
     .print*,'The Second Lightest Higgs H_2 Couplings : NHC_H(NC,2)'
      IF(IH.EQ.3)
     .print*,'The Heaviest Higgs H_3 Couplings : NHC_H(NC,3)'
      print*,'---------------------------------------------------------'
      NC=1 ! ee
      WRITE(*,  1) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=4 ! muon muon
      WRITE(*,  2) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=7 ! tau tau
      WRITE(*,  3) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=10 ! d d
      WRITE(*,  4) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=13 ! s s
      WRITE(*,  5) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=16 ! b b
      WRITE(*,  6) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=19 ! u u
      WRITE(*,  7) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=22 ! c c
      WRITE(*,  8) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=25 ! t t
      WRITE(*,  9) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=28 ! n1 n1
      WRITE(*, 10) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=31 ! n2 n2
      WRITE(*, 11) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=34 ! n3 n3
      WRITE(*, 12) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=37 ! n4 n4
      WRITE(*, 13) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=40 ! n1 n2
      WRITE(*, 14) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=43 ! n1 n3
      WRITE(*, 15) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=46 ! n1 n4
      WRITE(*, 16) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=49 ! n2 n3
      WRITE(*, 17) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=52 ! n2 n4
      WRITE(*, 18) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=55 ! n3 n4
      WRITE(*, 19) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=58 ! c1 c1
      WRITE(*, 20) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=61 ! c1+ c2-
      WRITE(*, 21) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=64 ! c2+ c1-
      WRITE(*, 22) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=67 ! c2 c2
      WRITE(*, 23) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,101)    NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      WRITE(*,102)    NC+2,DREAL(NHC_H(NC+2,IH)),DIMAG(NHC_H(NC+2,IH))
      NC=70 ! V V
      WRITE(*, 24) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=71 ! st1 st1
      WRITE(*, 25) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=72 ! st1* st2
      WRITE(*, 26) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=73 ! st2* st1
      WRITE(*, 27) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=74 ! st2 st2
      WRITE(*, 28) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=75 ! sb1 sb1
      WRITE(*, 29) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=76 ! sb1* sb2
      WRITE(*, 30) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=77 ! sb2* sb1
      WRITE(*, 31) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=78 ! sb2 sb2
      WRITE(*, 32) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=79 ! stau1* stau1
      WRITE(*,279) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=80 ! stau1* stau2
      WRITE(*,280) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=81 ! stau2* stau1
      WRITE(*,281) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=82 ! stau2* stau2
      WRITE(*,282) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=83 ! snu3* snu3
      WRITE(*,283) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=84 ! g g
      WRITE(*, 33) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,103) NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      NC=86 ! H+ H-
      WRITE(*, 34) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=87 ! H+ W-
      WRITE(*, 35) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      NC=88 ! p p
      WRITE(*, 36) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,103) NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      NC=90 ! g g (0)
      WRITE(*,137) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,103) NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      NC=92 ! p p (0)
      WRITE(*,138) IH,NC,DREAL(NHC_H(NC,IH)),DIMAG(NHC_H(NC,IH))
      WRITE(*,103) NC+1,DREAL(NHC_H(NC+1,IH)),DIMAG(NHC_H(NC+1,IH))
      IF(ICPRI.EQ.6) THEN
       IF(IC.EQ.1.OR.IC.EQ.2) GOTO 999
      ENDIF
      print*,'---------------------------------------------------------'
      ENDIF
*
      IF(ICPRI.EQ.6.OR.ICPRI.EQ.4) THEN
      print*,'---------------------------------------------------------'
      print*,'Charged Higgs Couplings     : CHC_H(NC) '
      print*,'---------------------------------------------------------'
      NC=1 ! e nu
      WRITE(*, 47) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=4 ! mu nu
      WRITE(*, 48) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=7 ! tau nu
      WRITE(*, 49) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=10 ! u d
      WRITE(*, 50) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=13 ! c s
      WRITE(*, 51) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=16 ! t b
      WRITE(*, 52) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=19 ! n1 c1
      WRITE(*, 53) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=22 ! n1 c2
      WRITE(*, 54) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=25 ! n2 c1
      WRITE(*, 55) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=28 ! n2 c2
      WRITE(*, 56) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=31 ! n3 c1
      WRITE(*, 57) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=34 ! n3 c2
      WRITE(*, 58) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=37 ! n4 c1
      WRITE(*, 59) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=40 ! n4 c2
      WRITE(*, 60) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      WRITE(*,101) NC+1,DREAL(CHC_H(NC+1)),DIMAG(CHC_H(NC+1))
      WRITE(*,102) NC+2,DREAL(CHC_H(NC+2)),DIMAG(CHC_H(NC+2))
      NC=43 ! stop1* sbottom1
      WRITE(*,104) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      NC=44 ! stop1* sbottom2
      WRITE(*,105) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      NC=45 ! stop2* sbottom1
      WRITE(*,106) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      NC=46 ! stop2* sbottom2
      WRITE(*,107) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      NC=47 ! snu3* stau1
      WRITE(*,108) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      NC=48 ! snu3* stau2
      WRITE(*,109) NC,DREAL(CHC_H(NC)),DIMAG(CHC_H(NC))
      print*,'---------------------------------------------------------'
      ENDIF
*
      IF(ICPRI.EQ.6.OR.ICPRI.EQ.5) THEN 
      print*,'---------------------------------------------------------'
      print*,'Higgs Boson Self Couplings  : SHC_H(NC) '
      print*,'---------------------------------------------------------'
      NC=1
      WRITE(*, 37) NC,SHC_H(NC)
      NC=2
      WRITE(*, 38) NC,SHC_H(NC)
      NC=3
      WRITE(*, 39) NC,SHC_H(NC)
      NC=4
      WRITE(*, 40) NC,SHC_H(NC)
      NC=5
      WRITE(*, 41) NC,SHC_H(NC)
      NC=6
      WRITE(*, 42) NC,SHC_H(NC)
      NC=7
      WRITE(*, 43) NC,SHC_H(NC)
      NC=8
      WRITE(*, 44) NC,SHC_H(NC)
      NC=9
      WRITE(*, 45) NC,SHC_H(NC)
      NC=10
      WRITE(*, 46) NC,SHC_H(NC)
      NC=11
      WRITE(*,151) NC,SHC_H(NC)
      NC=12
      WRITE(*,152) NC,SHC_H(NC)
      NC=13
      WRITE(*,153) NC,SHC_H(NC)
      NC=14
      WRITE(*,154) NC,SHC_H(NC)
      NC=15
      WRITE(*,155) NC,SHC_H(NC)
      NC=16
      WRITE(*,156) NC,SHC_H(NC)
      NC=17
      WRITE(*,157) NC,SHC_H(NC)
      NC=18
      WRITE(*,158) NC,SHC_H(NC)
      NC=19
      WRITE(*,159) NC,SHC_H(NC)
      NC=20
      WRITE(*,160) NC,SHC_H(NC)
      NC=21
      WRITE(*,161) NC,SHC_H(NC)
      NC=22
      WRITE(*,162) NC,SHC_H(NC)
      NC=23
      WRITE(*,163) NC,SHC_H(NC)
      NC=24
      WRITE(*,164) NC,SHC_H(NC)
      NC=25
      WRITE(*,165) NC,SHC_H(NC)
      NC=26
      WRITE(*,166) NC,SHC_H(NC)
      NC=27
      WRITE(*,167) NC,SHC_H(NC)
      NC=28
      WRITE(*,168) NC,SHC_H(NC)
      NC=29
      WRITE(*,169) NC,SHC_H(NC)
      NC=30
      WRITE(*,170) NC,SHC_H(NC)
      NC=31
      WRITE(*,171) NC,SHC_H(NC)
      NC=32
      WRITE(*,172) NC,SHC_H(NC)
      NC=33
      WRITE(*,173) NC,SHC_H(NC)
      NC=34
      WRITE(*,174) NC,SHC_H(NC)
      NC=35
      WRITE(*,175) NC,SHC_H(NC)
      print*,'---------------------------------------------------------'
      ENDIF
*
 1    FORMAT(1X,'H',I1,' e e               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 2    FORMAT(1X,'H',I1,' mu mu             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 3    FORMAT(1X,'H',I1,' tau tau           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 4    FORMAT(1X,'H',I1,' d d               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 5    FORMAT(1X,'H',I1,' s s               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 6    FORMAT(1X,'H',I1,' b b               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 7    FORMAT(1X,'H',I1,' u u               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 8    FORMAT(1X,'H',I1,' c c               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 9    FORMAT(1X,'H',I1,' t t               [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 10   FORMAT(1X,'H',I1,' N1 N1             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 11   FORMAT(1X,'H',I1,' N2 N2             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 12   FORMAT(1X,'H',I1,' N3 N3             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 13   FORMAT(1X,'H',I1,' N4 N4             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 14   FORMAT(1X,'H',I1,' N1 N2             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 15   FORMAT(1X,'H',I1,' N1 N3             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 16   FORMAT(1X,'H',I1,' N1 N4             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 17   FORMAT(1X,'H',I1,' N2 N3             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 18   FORMAT(1X,'H',I1,' N2 N4             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 19   FORMAT(1X,'H',I1,' N3 N4             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 20   FORMAT(1X,'H',I1,' C1+ C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 21   FORMAT(1X,'H',I1,' C1+ C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 22   FORMAT(1X,'H',I1,' C2+ C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 23   FORMAT(1X,'H',I1,' C2+ C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 24   FORMAT(1X,'H',I1,' V V               [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 25   FORMAT(1X,'H',I1,' ST1* ST1          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 26   FORMAT(1X,'H',I1,' ST1* ST2          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 27   FORMAT(1X,'H',I1,' ST2* ST1          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 28   FORMAT(1X,'H',I1,' ST2* ST2          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 29   FORMAT(1X,'H',I1,' SB1* SB1          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 30   FORMAT(1X,'H',I1,' SB1* SB2          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',          E10.4,')')
 31   FORMAT(1X,'H',I1,' SB2* SB1          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 32   FORMAT(1X,'H',I1,' SB2* SB2          [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
279   FORMAT(1X,'H',I1,' STA1* STA1        [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
280   FORMAT(1X,'H',I1,' STA1* STA2        [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
281   FORMAT(1X,'H',I1,' STA2* STA1        [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
282   FORMAT(1X,'H',I1,' STA2* STA2        [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
283   FORMAT(1X,'H',I1,' SNU3* SNU3        [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 33   FORMAT(1X,'H',I1,' glue glue         [NC=',I2,']:'
     .         ,' S =(',E10.4,',',E10.4,')')
 34   FORMAT(1X,'H',I1,' CH+ CH-           [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 35   FORMAT(1X,'H',I1,' CH+ W-            [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 36   FORMAT(1X,'H',I1,' photon photon     [NC=',I2,']:'
     .         ,' S =(',E10.4,',',E10.4,')')
137   FORMAT(1X,'H',I1,' glue glue    (M=0)[NC=',I2,']:'
     .         ,' S =(',E10.4,',',E10.4,')')
138   FORMAT(1X,'H',I1,' photon photon(M=0)[NC=',I2,']:'
     .         ,' S =(',E10.4,',',E10.4,')')
 37   FORMAT(1X,'H3 H3 H3             [NC=',I2,']:',' G = ',E10.4)
 38   FORMAT(1X,'H3 H3 H2             [NC=',I2,']:',' G = ',E10.4)
 39   FORMAT(1X,'H3 H3 H1             [NC=',I2,']:',' G = ',E10.4)
 40   FORMAT(1X,'H3 H2 H2             [NC=',I2,']:',' G = ',E10.4)
 41   FORMAT(1X,'H3 H2 H1             [NC=',I2,']:',' G = ',E10.4)
 42   FORMAT(1X,'H3 H1 H1             [NC=',I2,']:',' G = ',E10.4)
 43   FORMAT(1X,'H2 H2 H2             [NC=',I2,']:',' G = ',E10.4)
 44   FORMAT(1X,'H2 H2 H1             [NC=',I2,']:',' G = ',E10.4)
 45   FORMAT(1X,'H2 H1 H1             [NC=',I2,']:',' G = ',E10.4)
 46   FORMAT(1X,'H1 H1 H1             [NC=',I2,']:',' G = ',E10.4)
151   FORMAT(1X,'H1 CH+ CH-           [NC=',I2,']:',' G = ',E10.4)
152   FORMAT(1X,'H2 CH+ CH-           [NC=',I2,']:',' G = ',E10.4)
153   FORMAT(1X,'H3 CH+ CH-           [NC=',I2,']:',' G = ',E10.4)
154   FORMAT(1X,'H3 H3 H3 H3          [NC=',I2,']:',' G = ',E10.4)
155   FORMAT(1X,'H3 H3 H3 H2          [NC=',I2,']:',' G = ',E10.4)
156   FORMAT(1X,'H3 H3 H3 H1          [NC=',I2,']:',' G = ',E10.4)
157   FORMAT(1X,'H3 H3 H2 H2          [NC=',I2,']:',' G = ',E10.4)
158   FORMAT(1X,'H3 H3 H2 H1          [NC=',I2,']:',' G = ',E10.4)
159   FORMAT(1X,'H3 H3 H1 H1          [NC=',I2,']:',' G = ',E10.4)
160   FORMAT(1X,'H3 H2 H2 H2          [NC=',I2,']:',' G = ',E10.4)
161   FORMAT(1X,'H3 H2 H2 H1          [NC=',I2,']:',' G = ',E10.4)
162   FORMAT(1X,'H3 H2 H1 H1          [NC=',I2,']:',' G = ',E10.4)
163   FORMAT(1X,'H3 H1 H1 H1          [NC=',I2,']:',' G = ',E10.4)
164   FORMAT(1X,'H2 H2 H2 H2          [NC=',I2,']:',' G = ',E10.4)
165   FORMAT(1X,'H2 H2 H2 H1          [NC=',I2,']:',' G = ',E10.4)
166   FORMAT(1X,'H2 H2 H1 H1          [NC=',I2,']:',' G = ',E10.4)
167   FORMAT(1X,'H2 H1 H1 H1          [NC=',I2,']:',' G = ',E10.4)
168   FORMAT(1X,'H1 H1 H1 H1          [NC=',I2,']:',' G = ',E10.4)
169   FORMAT(1X,'H3 H3 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
170   FORMAT(1X,'H3 H2 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
171   FORMAT(1X,'H3 H1 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
172   FORMAT(1X,'H2 H2 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
173   FORMAT(1X,'H2 H1 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
174   FORMAT(1X,'H1 H1 CH+ CH-        [NC=',I2,']:',' G = ',E10.4)
175   FORMAT(1X,'CH+ CH- CH+ CH-      [NC=',I2,']:',' G = ',E10.4)
 47   FORMAT(1X,'CH+ e nu             [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 48   FORMAT(1X,'CH+ mu nu            [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 49   FORMAT(1X,'CH+ tau nu           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 50   FORMAT(1X,'CH+ u d              [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 51   FORMAT(1X,'CH+ c s              [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 52   FORMAT(1X,'CH+ t b              [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 53   FORMAT(1X,'CH+ N1 C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 54   FORMAT(1X,'CH+ N1 C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 55   FORMAT(1X,'CH+ N2 C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 56   FORMAT(1X,'CH+ N2 C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 57   FORMAT(1X,'CH+ N3 C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 58   FORMAT(1X,'CH+ N3 C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 59   FORMAT(1X,'CH+ N4 C1-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 60   FORMAT(1X,'CH+ N4 C2-           [NC=',I2,']:'
     .         ,' GF=(',E10.4,',',E10.4,')')
 104  FORMAT(1X,'CH+ ST1* SB1         [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 105  FORMAT(1X,'CH+ ST1* SB2         [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 106  FORMAT(1X,'CH+ ST2* SB1         [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 107  FORMAT(1X,'CH+ ST2* SB2         [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 108  FORMAT(1X,'CH+ SNU3* STA1       [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 109  FORMAT(1X,'CH+ SNU3* STA2       [NC=',I2,']:'
     .         ,' G =(',E10.4,',',E10.4,')')
 101  FORMAT (22X,'[NC=',I2,']: GS=(',E10.4,',',E10.4,')')
 102  FORMAT (22X,'[NC=',I2,']: GP=(',E10.4,',',E10.4,')')
 103  FORMAT (22X,'[NC=',I2,']: P =(',E10.4,',',E10.4,')')
*
      RETURN
      END

      SUBROUTINE HNICJ(IN,JC,N_N,U_L,U_R,SB,CB,TW,GW,GF,GS,GP)
************************************************************************
*
* Charged Higgs(+)-NeutralinoI-CharginoJ(-) Coupling
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*Input Array
      COMPLEX*16 N_N(4,4),U_L(2,2),U_R(2,2)
*Output
      COMPLEX*16 GF,GS,GP
*Local
      COMPLEX*16 XI
*
      XI=DCMPLX(0.D0,1.D0)
*
      GF=DCMPLX(GW/DSQRT(2.D0),0.D0)
      GS=1.D0/2.D0*(
     .   SB*(DSQRT(2.D0)*DCONJG(N_N(IN,3))*DCONJG(U_L(JC,1))
     .      -DCONJG(N_N(IN,2)+TW*N_N(IN,1))*DCONJG(U_L(JC,2)))
     .  +CB*(DSQRT(2.D0)*N_N(IN,4)*DCONJG(U_R(JC,1))
     .      +(N_N(IN,2)+TW*N_N(IN,1))*DCONJG(U_R(JC,2))))
      GP=XI/2.D0*(
     .   SB*(DSQRT(2.D0)*DCONJG(N_N(IN,3))*DCONJG(U_L(JC,1))
     .      -DCONJG(N_N(IN,2)+TW*N_N(IN,1))*DCONJG(U_L(JC,2)))
     .  -CB*(DSQRT(2.D0)*N_N(IN,4)*DCONJG(U_R(JC,1))
     .      +(N_N(IN,2)+TW*N_N(IN,1))*DCONJG(U_R(JC,2))))
*
      RETURN
      END

      SUBROUTINE HSTSB(IST,ISB,HT,HB,STMIX,SBMIX,COUPLING)   
************************************************************************
*
* Charged Higgs-stopIJ*-sbottomIK Coupling
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
*-----------------------------------------------------------------------
      COMPLEX*16 STMIX(2,2),SBMIX(2,2),COUPLING
      COMPLEX*16 XI,HT,HB,GABC(2,2)
*-----------------------------------------------------------------------
      V1   = V_H*CB_H
      V2   = V_H*SB_H
*GABC(STOP*_L,SBOTTOM_L)
      GABC(1,1)=-SB_H*(-1.D0/DSQRT(2.D0)*(CDABS(HB)**2-GW_H**2/2.D0)*V1)
     .          +CB_H*( 1.D0/DSQRT(2.D0)*(CDABS(HT)**2-GW_H**2/2.D0)*V2)
*GABC(STOP*_L,SBOTTOM_R)
      GABC(1,2)=DCONJG(HB*AB_H)*SB_H+DCONJG(HB)*MU_H*CB_H
*GABC(STOP*_R,SBOTTOM_L)
      GABC(2,1)=HT*DCONJG(MU_H)*SB_H+HT*AT_H*CB_H
*GABC(STOP*_R,SBOTTOM_R)
      GABC(2,2)=HT*DCONJG(HB)*(V2*SB_H+V1*CB_H)/DSQRT(2.D0)
*-----------------------------------------------------------------------
*                           A IST         B ISB       A B
      COUPLING=DCONJG(STMIX(1,IST))*SBMIX(1,ISB)*GABC(1,1)
     .        +DCONJG(STMIX(1,IST))*SBMIX(2,ISB)*GABC(1,2)
     .        +DCONJG(STMIX(2,IST))*SBMIX(1,ISB)*GABC(2,1)
     .        +DCONJG(STMIX(2,IST))*SBMIX(2,ISB)*GABC(2,2)
      COUPLING=COUPLING/V_H
*-----------------------------------------------------------------------
      RETURN
      END


      SUBROUTINE NHSELF3(I,J,K,SB,CB,OMIX,LR_H,LC_H,G)
************************************************************************
*
* Triple Higgs self-coupling
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*Input Array
      REAL*8     OMIX(3,3),LR_H(4)
      COMPLEX*16 LC_H(3)
*Local
      REAL*8      L1,L2,L3,L4
      COMPLEX*16  L5,L6,L7
      REAL*8      P(10)
*
      L1=LR_H(1)
      L2=LR_H(2)
      L3=LR_H(3)
      L4=LR_H(4)
      L5=LC_H(1)
      L6=LC_H(2)
      L7=LC_H(3)
*-----------------------------------------------------------------------
*G( a- a- a)
      G333=SB*CB*DIMAG(L5)-SB**2*DIMAG(L6)/2.D0-CB**2*DIMAG(L7)/2.D0
*G( a- a-p2)
      G332=CB**2*SB*L2+SB**3*(L3+L4)/2.D0-SB*(1.D0+CB**2)*DREAL(L5)
     .    +SB**2*CB*DREAL(L6)/2.D0+CB*(CB**2-2.D0*SB**2)*DREAL(L7)/2.D0
*G( a- a-p1)
      G331=SB**2*CB*L1+CB**3*(L3+L4)/2.D0-CB*(1.D0+SB**2)*DREAL(L5)
     .    +CB**2*SB*DREAL(L7)/2.D0+SB*(SB**2-2.D0*CB**2)*DREAL(L6)/2.D0
*G( a-p2-p2)
      G322=-SB*CB*DIMAG(L5)-(1.D0+2.D0*SB**2)*DIMAG(L7)/2.D0
*G( a-p2-p1)
      G321=-2.D0*DIMAG(L5)-CB*SB*DIMAG(L6+L7)
*G( a-p1-p1)
      G311=-SB*CB*DIMAG(L5)-(1.D0+2.D0*CB**2)*DIMAG(L6)/2.D0
*G(p2-p2-p2)
      G222=SB*L2+CB*DREAL(L7)/2.D0
*G(p2-p2-p1)
      G221=CB*(L3+L4)/2.D0+CB*DREAL(L5)+3.D0/2.D0*SB*DREAL(L7)
*G(p2-p1-p1)
      G211=SB*(L3+L4)/2.D0+SB*DREAL(L5)+3.D0/2.D0*CB*DREAL(L6)
*G(p1-p1-p1)
      G111=CB*L1+SB*DREAL(L6)/2.D0
*To compare with old conventions
*      print*,'G333   ',g333
*      print*,'G332/3 ',g332/3.d0
*      print*,'G331/3 ',g331/3.d0
*      print*,'G322/3 ',g322/3.d0
*      print*,'G321/6 ',g321/6.d0
*      print*,'G311/3 ',g311/3.d0
*      print*,'G222   ',g222
*      print*,'G221/3 ',g221/3.d0
*      print*,'G211/3 ',g211/3.d0
*      print*,'G111   ',g111
*-----------------------------------------------------------------------
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.3)) SYMFAC=6.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.2)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.2).AND.(K.EQ.2)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.2).AND.(K.EQ.1)) SYMFAC=1.D0
      IF((I.EQ.3).AND.(J.EQ.1).AND.(K.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.2).AND.(J.EQ.2).AND.(K.EQ.2)) SYMFAC=6.D0
      IF((I.EQ.2).AND.(J.EQ.2).AND.(K.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.2).AND.(J.EQ.1).AND.(K.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1)) SYMFAC=6.D0
*
      DO IP=1,10
        IF(IP.EQ.1) THEN
         IA=3
         IB=3
         IC=3
        ENDIF
        IF(IP.EQ.2) THEN
         IA=3
         IB=3
         IC=2
        ENDIF
        IF(IP.EQ.3) THEN
         IA=3
         IB=3
         IC=1
        ENDIF
        IF(IP.EQ.4) THEN
         IA=3
         IB=2
         IC=2
        ENDIF
        IF(IP.EQ.5) THEN
         IA=3
         IB=2
         IC=1
        ENDIF
        IF(IP.EQ.6) THEN
         IA=3
         IB=1
         IC=1
        ENDIF
        IF(IP.EQ.7) THEN
         IA=2
         IB=2
         IC=2
        ENDIF
        IF(IP.EQ.8) THEN
         IA=2
         IB=2
         IC=1
        ENDIF
        IF(IP.EQ.9) THEN
         IA=2
         IB=1
         IC=1
        ENDIF
        IF(IP.EQ.10) THEN
         IA=1
         IB=1
         IC=1
        ENDIF
*
      P(IP)=(OMIX(IA,I)*OMIX(IB,J)*OMIX(IC,K)
     .      +OMIX(IA,I)*OMIX(IB,K)*OMIX(IC,J)
     .      +OMIX(IA,J)*OMIX(IB,I)*OMIX(IC,K)
     .      +OMIX(IA,J)*OMIX(IB,K)*OMIX(IC,I)
     .      +OMIX(IA,K)*OMIX(IB,I)*OMIX(IC,J)
     .      +OMIX(IA,K)*OMIX(IB,J)*OMIX(IC,I))/SYMFAC
      ENDDO
      G=P( 1)*G333
     . +P( 2)*G332
     . +P( 3)*G331
     . +P( 4)*G322
     . +P( 5)*G321
     . +P( 6)*G311
     . +P( 7)*G222
     . +P( 8)*G221
     . +P( 9)*G211
     . +P(10)*G111
*
      RETURN
      END

      SUBROUTINE NHSELF4(I,J,K,L,SB,CB,OMIX,LR_H,LC_H,G)
************************************************************************
*
* Quartic Higgs self-coupling
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*Input Array
      REAL*8     OMIX(3,3),LR_H(4)
      COMPLEX*16 LC_H(3)
*Local
      REAL*8      L1,L2,L3,L4,L34
      COMPLEX*16  L5,L6,L7
      REAL*8      P(15)
*
      L1=LR_H(1)
      L2=LR_H(2)
      L3=LR_H(3)
      L4=LR_H(4)
      L34=(L3+L4)/2.D0
      L5=LC_H(1)
      L6=LC_H(2)
      L7=LC_H(3)
*-----------------------------------------------------------------------
*G4H+
      G4HP=SB**4*L1+CB**4*L2+SB**2*CB**2*(L3+L4)
     .    +2.D0*SB**2*CB**2*DREAL(L5)-2.D0*SB**3*CB*DREAL(L6)
     .    -2.D0*SB*CB**3*DREAL(L7)
*G( a- a- a- a)
      G3333=G4HP/4.D0
*G( a- a- a-p2)
      G3332=SB**2*CB*DIMAG(L5)
     .     -SB**3*DIMAG(L6)/2.D0-SB*CB**2*DIMAG(L7)/2.D0
*G( a- a- a-p1)
      G3331=SB*CB**2*DIMAG(L5)
     .     -SB**2*CB*DIMAG(L6)/2.D0-CB**3*DIMAG(L7)/2.D0
*G( a- a-p2-p2)
      G3322=(CB**2*L2+SB**2*L34-SB**2*DREAL(L5)-SB*CB*DREAL(L7))/2.D0
*G( a- a-p2-p1)
      G3321=-2.D0*SB*CB*DREAL(L5)
     .     +SB**2*DREAL(L6)/2.D0+CB**2*DREAL(L7)/2.D0
*G( a- a-p1-p1)
      G3311=(SB**2*L1+CB**2*L34-CB**2*DREAL(L5)-SB*CB*DREAL(L6))/2.D0
*G( a-p2-p2-p2)
      G3222=-SB*DIMAG(L7)/2.D0
*G( a-p2-p2-p1)
      G3221=-SB*DIMAG(L5)-CB*DIMAG(L7)/2.D0
*G( a-p2-p1-p1)
      G3211=-CB*DIMAG(L5)-SB*DIMAG(L6)/2.D0
*G( a-p1-p1-p1)
      G3111=-CB*DIMAG(L6)/2.D0
*G(p2-p2-p2-p2)
      G2222=L2/4.D0
*G(p2-p2-p2-p1)
      G2221=DREAL(L7)/2.D0
*G(p2-p2-p1-p1)
      G2211=L34/2.D0+DREAL(L5)/2.D0
*G(p2-p1-p1-p1)
      G2111=DREAL(L6)/2.D0
*G(p1-p1-p1-p1)
      G1111=L1/4.D0
*G( a- a-ch-ch)
      G33=G4HP
*G( a-p2-ch-ch)
      G32=2.D0*SB**2*CB*DIMAG(L5)-SB**3*DIMAG(L6)-SB*CB**2*DIMAG(L7)
*G( a-p1-ch-ch)
      G31=2.D0*CB**2*SB*DIMAG(L5)-CB**3*DIMAG(L7)-CB*SB**2*DIMAG(L6)
*G(p2-p2-ch-ch)
      G22=CB**2*L2+SB**2*L3/2.D0-SB*CB*DREAL(L7)
*G(p2-p1-ch-ch)
      G21=-SB*CB*L4-2.D0*CB*SB*DREAL(L5)+SB**2*DREAL(L6)+CB**2*DREAL(L7)
*G(p1-p1-ch-ch)
      G11=SB**2*L1+CB**2*L3/2.D0-SB*CB*DREAL(L6)
*-----------------------------------------------------------------------
      IF(K.EQ.5) THEN
       IF(I.EQ.5) THEN
        G=G4HP
        RETURN
       ENDIF
        G=(OMIX(3,I)*OMIX(3,J)+OMIX(3,J)*OMIX(3,I))*G33
     .   +(OMIX(3,I)*OMIX(2,J)+OMIX(3,J)*OMIX(2,I))*G32
     .   +(OMIX(3,I)*OMIX(1,J)+OMIX(3,J)*OMIX(1,I))*G31
     .   +(OMIX(2,I)*OMIX(2,J)+OMIX(2,J)*OMIX(2,I))*G22
     .   +(OMIX(2,I)*OMIX(1,J)+OMIX(2,J)*OMIX(1,I))*G21
     .   +(OMIX(1,I)*OMIX(1,J)+OMIX(1,J)*OMIX(1,I))*G11
        IF(I.EQ.J) G=G/2.D0
        RETURN
      ENDIF
    
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.3).AND.(L.EQ.3)) SYMFAC=24.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.3).AND.(L.EQ.2)) SYMFAC=6.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.3).AND.(L.EQ.1)) SYMFAC=6.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.2).AND.(L.EQ.2)) SYMFAC=4.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.2).AND.(L.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.3).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=4.D0 ! 33xxx
      IF((I.EQ.3).AND.(J.EQ.2).AND.(K.EQ.2).AND.(L.EQ.2)) SYMFAC=6.D0
      IF((I.EQ.3).AND.(J.EQ.2).AND.(K.EQ.2).AND.(L.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.2).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=2.D0
      IF((I.EQ.3).AND.(J.EQ.1).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=6.D0 ! 3xxxx
      IF((I.EQ.2).AND.(J.EQ.2).AND.(K.EQ.2).AND.(L.EQ.2)) SYMFAC=24.D0
      IF((I.EQ.2).AND.(J.EQ.2).AND.(K.EQ.2).AND.(L.EQ.1)) SYMFAC=6.D0
      IF((I.EQ.2).AND.(J.EQ.2).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=4.D0
      IF((I.EQ.2).AND.(J.EQ.1).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=6.D0
      IF((I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1).AND.(L.EQ.1)) SYMFAC=24.D0
*
      DO IP=1,15
        IF(IP.EQ.1) THEN
         IA=3
         IB=3
         IC=3
         ID=3
        ENDIF
        IF(IP.EQ.2) THEN
         IA=3
         IB=3
         IC=3
         ID=2
        ENDIF
        IF(IP.EQ.3) THEN
         IA=3
         IB=3
         IC=3
         ID=1
        ENDIF
        IF(IP.EQ.4) THEN
         IA=3
         IB=3
         IC=2
         ID=2
        ENDIF
        IF(IP.EQ.5) THEN
         IA=3
         IB=3
         IC=2
         ID=1
        ENDIF
        IF(IP.EQ.6) THEN
         IA=3
         IB=3
         IC=1
         ID=1
        ENDIF
        IF(IP.EQ.7) THEN
         IA=3
         IB=2
         IC=2
         ID=2
        ENDIF
        IF(IP.EQ.8) THEN
         IA=3
         IB=2
         IC=2
         ID=1
        ENDIF
        IF(IP.EQ.9) THEN
         IA=3
         IB=2
         IC=1
         ID=1
        ENDIF
        IF(IP.EQ.10) THEN
         IA=3
         IB=1
         IC=1
         ID=1
        ENDIF
        IF(IP.EQ.11) THEN
         IA=2
         IB=2
         IC=2
         ID=2
        ENDIF
        IF(IP.EQ.12) THEN
         IA=2
         IB=2
         IC=2
         ID=1
        ENDIF
        IF(IP.EQ.13) THEN
         IA=2
         IB=2
         IC=1
         ID=1
        ENDIF
        IF(IP.EQ.14) THEN
         IA=2
         IB=1
         IC=1
         ID=1
        ENDIF
        IF(IP.EQ.15) THEN
         IA=1
         IB=1
         IC=1
         ID=1
        ENDIF
*
      P(IP)=(0.D0
     .      +OMIX(IA,I)*OMIX(IB,J)*OMIX(IC,K)*OMIX(ID,L)
     .      +OMIX(IA,I)*OMIX(IB,J)*OMIX(IC,L)*OMIX(ID,K)
     .      +OMIX(IA,I)*OMIX(IB,K)*OMIX(IC,J)*OMIX(ID,L)
     .      +OMIX(IA,I)*OMIX(IB,K)*OMIX(IC,L)*OMIX(ID,J)
     .      +OMIX(IA,I)*OMIX(IB,L)*OMIX(IC,J)*OMIX(ID,K)
     .      +OMIX(IA,I)*OMIX(IB,L)*OMIX(IC,K)*OMIX(ID,J)   !  I
     .      +OMIX(IA,J)*OMIX(IB,K)*OMIX(IC,L)*OMIX(ID,I)
     .      +OMIX(IA,J)*OMIX(IB,K)*OMIX(IC,I)*OMIX(ID,L)
     .      +OMIX(IA,J)*OMIX(IB,L)*OMIX(IC,K)*OMIX(ID,I)
     .      +OMIX(IA,J)*OMIX(IB,L)*OMIX(IC,I)*OMIX(ID,K)
     .      +OMIX(IA,J)*OMIX(IB,I)*OMIX(IC,K)*OMIX(ID,L)
     .      +OMIX(IA,J)*OMIX(IB,I)*OMIX(IC,L)*OMIX(ID,K)   !  J
     .      +OMIX(IA,K)*OMIX(IB,L)*OMIX(IC,I)*OMIX(ID,J)
     .      +OMIX(IA,K)*OMIX(IB,L)*OMIX(IC,J)*OMIX(ID,I)
     .      +OMIX(IA,K)*OMIX(IB,I)*OMIX(IC,L)*OMIX(ID,J)
     .      +OMIX(IA,K)*OMIX(IB,I)*OMIX(IC,J)*OMIX(ID,L)
     .      +OMIX(IA,K)*OMIX(IB,J)*OMIX(IC,L)*OMIX(ID,I)
     .      +OMIX(IA,K)*OMIX(IB,J)*OMIX(IC,I)*OMIX(ID,L)   !  K
     .      +OMIX(IA,L)*OMIX(IB,I)*OMIX(IC,J)*OMIX(ID,K)
     .      +OMIX(IA,L)*OMIX(IB,I)*OMIX(IC,K)*OMIX(ID,J)
     .      +OMIX(IA,L)*OMIX(IB,J)*OMIX(IC,I)*OMIX(ID,K)
     .      +OMIX(IA,L)*OMIX(IB,J)*OMIX(IC,K)*OMIX(ID,I)
     .      +OMIX(IA,L)*OMIX(IB,K)*OMIX(IC,I)*OMIX(ID,J)
     .      +OMIX(IA,L)*OMIX(IB,K)*OMIX(IC,J)*OMIX(ID,I)   !  L
     .      )/SYMFAC
      ENDDO
*
      G=P( 1)*G3333          ! 33xx = 6
     . +P( 2)*G3332
     . +P( 3)*G3331
     . +P( 4)*G3322
     . +P( 5)*G3321
     . +P( 6)*G3311
     . +P( 7)*G3222          ! 32xx = 3
     . +P( 8)*G3221
     . +P( 9)*G3211
     . +P(10)*G3111          ! 3111    
     . +P(11)*G2222          ! 22xx = 3
     . +P(12)*G2221
     . +P(13)*G2211
     . +P(14)*G2111          ! 2111
     . +P(15)*G1111          ! 1111
*
*       print*,P( 1)*G3333          ! 33xx = 6
*       print*,P( 2)*G3332
*       print*,P( 3)*G3331
*       print*,P( 4)*G3322
*       print*,P( 5)*G3321
*       print*,P( 6)*G3311
*       print*,P( 7)*G3222          ! 32xx = 3
*       print*,P( 8)*G3221
*       print*,P( 9)*G3211
*       print*,P(10)*G3111          ! 3111    
*       print*,P(11)*G2222          ! 22xx = 3
*       print*,P(12)*G2221
*       print*,P(13)*G2211
*       print*,P(14)*G2111          ! 2111
*       print*,P(15)*G1111          ! 1111
*
      RETURN
      END





      SUBROUTINE HPP(IH,ASMH2,MTMH2,MBMH2,MCMH2
     .              ,M_C,MCH,HMASS,STMASS,SBMASS,STAUMASS
     .              ,NCMAX,NHC_H,SPHO,PPHO,SPP,PPP)
************************************************************************
*
* HiggsIH-photon-photon Coupling
*
* For each IH,
*
* SPP( 1) PPP(1) : bottom
* SPP( 2) PPP(2) : top
* SPP( 3) PPP(3) : charm
* SPP( 4) PPP(4) : tau lepton
* SPP( 5) PPP(5) : chargino1 chargino1
* SPP( 6) PPP(6) : chargino2 chargino2
* SPP( 7)        : stop1 stop1
* SPP( 8)        : stop2 stop2
* SPP( 9)        : sbottom1 sbottom1
* SPP(10)        : sbottom2 sbottom2
* SPP(11)        : W+W-
* SPP(12)        : charged Higgs charged Higgs
* SPP(13)        : stau1 stau1
* SPP(14)        : stau2 stau2
* SPP(15) PPP(7) : total
*
*
*JSL:26/OCT/2006 : (1) The running quark mass used in the loop functions
*                  HC_FSF and HC_FPF
*                  (2) The contributions from tau-lepton and charm-quark
*                  loops have included
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
*Input Array
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     M_C(2),HMASS(3),STMASS(2),SBMASS(2),STAUMASS(2)
*Output
      COMPLEX*16 SPHO,PPHO
      COMPLEX*16 SPP(15),PPP(7)
*Local
      COMPLEX*16 HC_FSF,HC_FPF,HC_F0,HC_F1
      REAL*8     NC
*
*-----------------------------------------------------------------------
      PI      = 2.D0*DASIN(1.D0)
*
      NC   = 3.D0
      QB   =-1.D0/3.D0
      QT   = 2.D0/3.D0
      QCHA = 2.D0/3.D0
      ITAU = 7
      IB   = 16
      ICHA = 22
      IT   = 25
      IC1  = 58
      IC2  = 67
      IHV  = 70
      IT1  = 71
      IT2  = 74
      IB1  = 75
      IB2  = 78
      ITU1 = 79
      ITU2 = 82
      ICH  = 86
*
      SPHO=0.D0
     .+2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+1,IH)*V_H/MBMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MBMH2**2)
     .+2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+1,IH)*V_H/MTMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MTMH2**2)
     .+2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+1,IH)*V_H/MCMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MCMH2**2)
     .+         2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+1,IH)*V_H/MTAU_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MTAU_H**2)
     .+         2.D0*NHC_H(IC1,IH)*NHC_H(IC1+1,IH)*V_H/M_C(1)
     .              *HC_FSF(HMASS(IH)**2/4.D0/M_C(1)**2)
     .+         2.D0*NHC_H(IC2,IH)*NHC_H(IC2+1,IH)*V_H/M_C(2)
     .              *HC_FSF(HMASS(IH)**2/4.D0/M_C(2)**2)
     .-2.D0*NC*QT**2*NHC_H(IT1,IH)*V_H**2/4.D0/STMASS(1)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STMASS(1)**2)
     .-2.D0*NC*QT**2*NHC_H(IT2,IH)*V_H**2/4.D0/STMASS(2)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STMASS(2)**2)
     .-2.D0*NC*QB**2*NHC_H(IB1,IH)*V_H**2/4.D0/SBMASS(1)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/SBMASS(1)**2)
     .-2.D0*NC*QB**2*NHC_H(IB2,IH)*V_H**2/4.D0/SBMASS(2)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/SBMASS(2)**2)
     .-NHC_H(IHV,IH)*HC_F1(HMASS(IH)**2/4.D0/MW_H**2)
     .-NHC_H(ICH,IH)*V_H**2/2.D0/MCH**2
     .              *HC_F0(HMASS(IH)**2/4.D0/MCH**2)
     .-2.D0*NHC_H(ITU1,IH)*V_H**2/4.D0/STAUMASS(1)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STAUMASS(1)**2)
     .-2.D0*NHC_H(ITU2,IH)*V_H**2/4.D0/STAUMASS(2)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STAUMASS(2)**2)

*      print*,'fillcoupl2.f:HPP',ih,mbmt_h,mb_s,mb_pole

      SPP(1)= 2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+1,IH)*V_H/MBMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MBMH2**2)
      SPP(2)= 2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+1,IH)*V_H/MTMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MTMH2**2)
      SPP(3)= 2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+1,IH)*V_H/MCMT_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MCMH2**2)
      SPP(4)= 2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+1,IH)*V_H/MTAU_H
     .              *HC_FSF(HMASS(IH)**2/4.D0/MTAU_H**2)
      SPP(5)= 2.D0*NHC_H(IC1,IH)*NHC_H(IC1+1,IH)*V_H/M_C(1)
     .              *HC_FSF(HMASS(IH)**2/4.D0/M_C(1)**2)
      SPP(6)= 2.D0*NHC_H(IC2,IH)*NHC_H(IC2+1,IH)*V_H/M_C(2)
     .              *HC_FSF(HMASS(IH)**2/4.D0/M_C(2)**2)
      SPP(7)=-2.D0*NC*QT**2*NHC_H(IT1,IH)*V_H**2/4.D0/STMASS(1)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STMASS(1)**2)
      SPP(8)=-2.D0*NC*QT**2*NHC_H(IT2,IH)*V_H**2/4.D0/STMASS(2)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/STMASS(2)**2)
      SPP(9)=-2.D0*NC*QB**2*NHC_H(IB1,IH)*V_H**2/4.D0/SBMASS(1)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/SBMASS(1)**2)
      SPP(10)=-2.D0*NC*QB**2*NHC_H(IB2,IH)*V_H**2/4.D0/SBMASS(2)**2
     .              *HC_F0(HMASS(IH)**2/4.D0/SBMASS(2)**2)
      SPP(11)=-NHC_H(IHV,IH)*HC_F1(HMASS(IH)**2/4.D0/MW_H**2)
      SPP(12)=-NHC_H(ICH,IH)*V_H**2/2.D0/MCH**2
     .              *HC_F0(HMASS(IH)**2/4.D0/MCH**2)
      SPP(13)=-2.D0*NHC_H(ITU1,IH)*V_H**2/4.D0/STAUMASS(1)**2
     .             *HC_F0(HMASS(IH)**2/4.D0/STAUMASS(1)**2)
      SPP(14)=-2.D0*NHC_H(ITU2,IH)*V_H**2/4.D0/STAUMASS(2)**2
     .             *HC_F0(HMASS(IH)**2/4.D0/STAUMASS(2)**2)
      SPP(15)=SPP(1)+SPP(2)+SPP(3)+SPP(4)+SPP(5)+SPP(6)+SPP(7)+SPP(8)
     .       +SPP(9)+SPP(10)+SPP(11)+SPP(12)+SPP(13)+SPP(14)

*JSL:28/AUG/2006 : the typo for the chargino contributions corrected
      PPHO=0.D0
     .+2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+2,IH)*V_H/MBMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MBMH2**2)
     .+2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+2,IH)*V_H/MTMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MTMH2**2)
     .+2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+2,IH)*V_H/MCMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MCMH2**2)
     .+         2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+2,IH)*V_H/MTAU_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MTAU_H**2)
     .+         2.D0*NHC_H(IC1,IH)*NHC_H(IC1+2,IH)*V_H/M_C(1)
     .              *HC_FPF(HMASS(IH)**2/4.D0/M_C(1)**2)
     .+         2.D0*NHC_H(IC2,IH)*NHC_H(IC2+2,IH)*V_H/M_C(2)
     .              *HC_FPF(HMASS(IH)**2/4.D0/M_C(2)**2)

      PPP(1)= 2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+2,IH)*V_H/MBMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MBMH2**2)
      PPP(2)= 2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+2,IH)*V_H/MTMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MTMH2**2)
      PPP(3)= 2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+2,IH)*V_H/MCMT_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MCMH2**2)
      PPP(4)= 2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+2,IH)*V_H/MTAU_H
     .              *HC_FPF(HMASS(IH)**2/4.D0/MTAU_H**2)
      PPP(5)= 2.D0*NHC_H(IC1,IH)*NHC_H(IC1+2,IH)*V_H/M_C(1)
     .              *HC_FPF(HMASS(IH)**2/4.D0/M_C(1)**2)
      PPP(6)= 2.D0*NHC_H(IC2,IH)*NHC_H(IC2+2,IH)*V_H/M_C(2)
     .              *HC_FPF(HMASS(IH)**2/4.D0/M_C(2)**2)
      PPP(7)=PPP(1)+PPP(2)+PPP(3)+PPP(4)+PPP(5)+PPP(6)
*
      RETURN
      END


      SUBROUTINE HPPSM(MH,MB_S,MT_S,MC_S,SPP,PPP)
************************************************************************
*
* SM Higgs-photon-photon Coupling
*
* For each IH,
*
* SPP(1) PPP(1) : bottom
* SPP(2) PPP(2) : top
* SPP(3) PPP(3) : charm
* SPP(4) PPP(4) : tau lepton
* SPP(5)        : W+W-
* SPP(9) PPP(9) : total
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
*Output
      COMPLEX*16 SPP(9),PPP(9)
*Local
      COMPLEX*16 HC_FSF,HC_FPF,HC_F0,HC_F1
      REAL*8     NC
*
      NC   = 3.D0
*=======================================================================
*
      DO I=1,9
       SPP(I)=DCMPLX(0.D0,0.D0)
      ENDDO
*
*      print*,'fillcoupl2.f:HPPSM',MH,MB_S,MT_S,MC_S
*
      QB=-1.D0/3.D0
      QT= 2.D0/3.D0
      QC= 2.D0/3.D0

      SPP(1)= 2.D0*NC*QB**2*HC_FSF(MH**2/4.D0/MB_S**2)
      SPP(2)= 2.D0*NC*QT**2*HC_FSF(MH**2/4.D0/MT_S**2)
      SPP(3)= 2.D0*NC*QC**2*HC_FSF(MH**2/4.D0/MC_S**2)
      SPP(4)= 2.D0*HC_FSF(MH**2/4.D0/MTAU_H**2)
      SPP(5)=-HC_F1(MH**2/4.D0/MW_H**2)
*
      SPP(9)=SPP(1)+SPP(2)+SPP(3)+SPP(4)+SPP(5)
*
*      print*,'fillcoupl2.f:HPPSM'
*     .      ,cdabs(SPP(1))
*     .      ,cdabs(SPP(2))
*     .      ,cdabs(SPP(3))
*     .      ,cdabs(SPP(4))
*     .      ,cdabs(SPP(5))
*     .      ,cdabs(SPP(9))
*
      RETURN
      END


      SUBROUTINE HCHCH(IH,CB,SB,LR_H,LC_H,OMIX,G)
************************************************************************
*
* HiggsIH-charged Higgs-charged Higgs coupling 
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*Input Array
      REAL*8     OMIX(3,3)
      REAL*8     LR_H(4)
      COMPLEX*16 LC_H(3)
*Output Coupling
      COMPLEX*16 G
*Local Parameters
      REAL*8      L1,L2,L3,L4
      COMPLEX*16  L5,L6,L7
*
      L1=LR_H(1)
      L2=LR_H(2)
      L3=LR_H(3)
      L4=LR_H(4)
      L5=LC_H(1)
      L6=LC_H(2)
      L7=LC_H(3)
*
      GA=2.D0*SB*CB*DIMAG(L5)-SB**2*DIMAG(L6)-CB**2*DIMAG(L7)
      GP1=2.D0*SB**2*CB*L1+CB**3*L3-SB**2*CB*L4-2.D0*SB**2*CB*DREAL(L5)
     .   +SB*(SB**2-2.D0*CB**2)*DREAL(L6)+SB*CB**2*DREAL(L7)
      GP2=2.D0*CB**2*SB*L2+SB**3*L3-CB**2*SB*L4-2.D0*CB**2*SB*DREAL(L5)
     .   +CB*(CB**2-2.D0*SB**2)*DREAL(L7)+CB*SB**2*DREAL(L6)
*
      G=DCMPLX(OMIX(1,IH)*GP1+OMIX(2,IH)*GP2+OMIX(3,IH)*GA,0.D0)
*
      RETURN
      END

      SUBROUTINE HPLAMBDA(HB_H,HT_H,STMASS,LR_H,LC_H)
************************************************************************
*
* Higgs Potential Coupling: Pilaftsis and Wagner NPB553(1999)3
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
*-----------------------------------------------------------------------
*Input Array 
      COMPLEX*16 HB_H,HT_H
      REAL*8 STMASS(2)
*-----------------------------------------------------------------------
*Output Array : LR_H(I)=LAMBDA_I, LC_H(I)=LAMBDA_(I+4)
      REAL*8     LR_H(4)
      COMPLEX*16 LC_H(3)
*-----------------------------------------------------------------------
      PI=2.D0*DASIN(1.D0)
      MSUSYSQ=1.D0/2.D0*(STMASS(1)**2+STMASS(2)**2)
      MSUSY=DSQRT(MSUSYSQ)
      HT=CDABS(HT_H)
      HB=CDABS(HB_H)
      XT=2.D0*CDABS(AT_H)**2/MSUSYSQ*(1.D0-CDABS(AT_H)**2/12.D0/MSUSYSQ)
      XB=2.D0*CDABS(AB_H)**2/MSUSYSQ*(1.D0-CDABS(AB_H)**2/12.D0/MSUSYSQ)
      XTB=(CDABS(AT_H)**2+CDABS(AB_H)**2+2.D0*DREAL(DCONJG(AB_H)*AT_H))
     .   /2.D0/MSUSYSQ-CDABS(MU_H)**2/MSUSYSQ
     .   -CDABS(CDABS(MU_H)**2-DCONJG(AB_H)*AT_H)**2/6.D0/MSUSYSQ**2
      GS2=4.D0*PI*ASMT_H
      T=DLOG(MSUSYSQ/MTPOLE_H**2)
*
      LR_H(1)=-(GW_H**2+GP_H**2)/8.D0*(1.D0-3.D0*HB**2/8.D0/PI**2*T)
     .        -3.D0/16.D0/PI**2*HB**4*(T+XB/2.D0+1.D0/16.D0/PI**2
     .         *(3.D0/2.D0*HB**2+HT**2/2.D0-8.D0*GS2)*(XB*T+T**2))
     .        +3.D0/192.D0/PI**2*HT**4*CDABS(MU_H)**4/MSUSYSQ**2
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(9.D0*HT**2-5.D0*HB**2-16.D0*GS2)*T)
      LR_H(2)=-(GW_H**2+GP_H**2)/8.D0*(1.D0-3.D0*HT**2/8.D0/PI**2*T)
     .        -3.D0/16.D0/PI**2*HT**4*(T+XT/2.D0+1.D0/16.D0/PI**2
     .         *(3.D0/2.D0*HT**2+HB**2/2.D0-8.D0*GS2)*(XT*T+T**2))
     .        +3.D0/192.D0/PI**2*HB**4*CDABS(MU_H)**4/MSUSYSQ**2
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(9.D0*HB**2-5.D0*HT**2-16.D0*GS2)*T)
      LR_H(3)=-(GW_H**2-GP_H**2)/4.D0
     .         *(1.D0-3.D0*(HT**2+HB**2)/16.D0/PI**2*T)
     .        -3.D0/8.D0/PI**2*HT**2*HB**2*(T+XTB/2.D0+1.D0/16.D0/PI**2
     .         *(HT**2+HB**2-8.D0*GS2)*(XTB*T+T**2))
     .        -3.D0/96.D0/PI**2*HT**4
     .    *(3.D0*CDABS(MU_H)**2/MSUSYSQ-CDABS(MU_H*AT_H)**2/MSUSYSQ**2)
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(6.D0*HT**2-2.D0*HB**2-16.D0*GS2)*T)
     .        -3.D0/96.D0/PI**2*HB**4
     .    *(3.D0*CDABS(MU_H)**2/MSUSYSQ-CDABS(MU_H*AB_H)**2/MSUSYSQ**2)
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(6.D0*HB**2-2.D0*HT**2-16.D0*GS2)*T)
      LR_H(4)=GW_H**2/2.D0
     .         *(1.D0-3.D0*(HT**2+HB**2)/16.D0/PI**2*T)
     .        +3.D0/8.D0/PI**2*HT**2*HB**2*(T+XTB/2.D0+1.D0/16.D0/PI**2
     .         *(HT**2+HB**2-8.D0*GS2)*(XTB*T+T**2))
     .        -3.D0/96.D0/PI**2*HT**4
     .    *(3.D0*CDABS(MU_H)**2/MSUSYSQ-CDABS(MU_H*AT_H)**2/MSUSYSQ**2)
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(6.D0*HT**2-2.D0*HB**2-16.D0*GS2)*T)
     .        -3.D0/96.D0/PI**2*HB**4
     .    *(3.D0*CDABS(MU_H)**2/MSUSYSQ-CDABS(MU_H*AB_H)**2/MSUSYSQ**2)
     .         *(1.D0+1.D0/16.D0/PI**2
     .           *(6.D0*HB**2-2.D0*HT**2-16.D0*GS2)*T)
*
      LC_H(1)=3.D0/192.D0/PI**2*HT**4*(MU_H*AT_H)**2/MSUSYSQ**2
     .         *(1.D0-1.D0/16.D0/PI**2
     .           *(2.D0*HB**2-6.D0*HT**2+16.D0*GS2)*T)
     .       +3.D0/192.D0/PI**2*HB**4*(MU_H*AB_H)**2/MSUSYSQ**2
     .         *(1.D0-1.D0/16.D0/PI**2
     .           *(2.D0*HT**2-6.D0*HB**2+16.D0*GS2)*T)
      LC_H(2)=-3.D0/96.D0/PI**2*HT**4*CDABS(MU_H)**2*MU_H*AT_H
     .         /MSUSYSQ**2
     .         *(1.D0-1.D0/16.D0/PI**2
     .           *(7.D0/2.D0*HB**2-15.D0/2.D0*HT**2+16.D0*GS2)*T)
     .        +3.D0/96.D0/PI**2*HB**4*MU_H/MSUSY
     .         *(6.D0*AB_H/MSUSY-CDABS(AB_H)**2*AB_H/MSUSY**3)
     .          *(1.D0-1.D0/16.D0/PI**2
     .            *(HT**2/2.D0-9.D0/2.D0*HB**2+16.D0*GS2)*T)
      LC_H(3)=-3.D0/96.D0/PI**2*HB**4*CDABS(MU_H)**2*MU_H*AB_H
     .         /MSUSYSQ**2
     .         *(1.D0-1.D0/16.D0/PI**2
     .           *(7.D0/2.D0*HT**2-15.D0/2.D0*HB**2+16.D0*GS2)*T)
     .        +3.D0/96.D0/PI**2*HT**4*MU_H/MSUSY
     .         *(6.D0*AT_H/MSUSY-CDABS(AT_H)**2*AT_H/MSUSY**3)
     .          *(1.D0-1.D0/16.D0/PI**2
     .            *(HB**2/2.D0-9.D0/2.D0*HT**2+16.D0*GS2)*T)
*
      RETURN
      END

      SUBROUTINE HGG(IH,HMASS,STMASS,SBMASS,NCMAX,NHC_H
     .              ,SGLUE,PGLUE,SGG,PGG)
************************************************************************
*
* HiggsIH-glue-glue Coupling
*
* For each IH,
*
* SGG(1) PGG(1) : bottom
* SGG(2) PGG(2) : top
* SGG(3)        : stop1 stop1
* SGG(4)        : stop2 stop2
* SGG(5)        : sbottom1 sbottom1
* SGG(6)        : sbottom2 sbottom2
* SGG(7) PGG(3) : total
*
*JSL:30/OCT/2006 : The running quark mass used in the loop functions
*                  HC_FSF and HC_FPF
*
*JSL:16/Mar/2012 : The pole quark mass used in the loop functions
*                  HC_FSF and HC_FPF
*
*                  Ref. M.~Spira, A.~Djouadi, D.~Graudenz and P.~M.~Zerwas,
*                  ``Higgs boson production at the LHC,''
*                  Nucl.Phys.B453 (1995) 17 [hep-ph/9504378].
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
*-----------------------------------------------------------------------
      COMPLEX*16 NHC_H(NCMAX,3)
      COMPLEX*16 SGLUE,PGLUE
      COMPLEX*16 HC_FSF,HC_FPF,HC_F0,HC_F1
      REAL*8     HMASS(3),STMASS(2),SBMASS(2)
      COMPLEX*16 SGG(7),PGG(3)
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*=======================================================================
*
      IBOT = 16
      ITOP = 25
      IT11 = 71
      IT22 = 74
      IB11 = 75
      IB22 = 78
*
      PI=2.D0*DASIN(1.D0)
*
*      print*,'fillcoupl2.f:HGG',1.d0/NHC_H(IBOT,IH),V_H/MBMT_H
*      print*,'fillcoupl2.f:HGG',1.d0/NHC_H(ITOP,IH),V_H/MTMT_H
*
      SGLUE=0.D0
     .+NHC_H(IBOT,IH)*NHC_H(IBOT+1,IH)*V_H/MBMT_H
     .               *HC_FSF(HMASS(IH)**2/4.D0/RAUX_H(1)**2)
     .+NHC_H(ITOP,IH)*NHC_H(ITOP+1,IH)*V_H/MTMT_H
     .               *HC_FSF(HMASS(IH)**2/4.D0/MTPOLE_H**2)
     .-NHC_H(IT11,IH)*V_H**2/4.D0/STMASS(1)**2
     .               *HC_F0(HMASS(IH)**2/4.D0/STMASS(1)**2)
     .-NHC_H(IT22,IH)*V_H**2/4.D0/STMASS(2)**2
     .               *HC_F0(HMASS(IH)**2/4.D0/STMASS(2)**2)
     .-NHC_H(IB11,IH)*V_H**2/4.D0/SBMASS(1)**2
     .               *HC_F0(HMASS(IH)**2/4.D0/SBMASS(1)**2)
     .-NHC_H(IB22,IH)*V_H**2/4.D0/SBMASS(2)**2
     .               *HC_F0(HMASS(IH)**2/4.D0/SBMASS(2)**2)

      SGG(1)= NHC_H(IBOT,IH)*NHC_H(IBOT+1,IH)*V_H/MBMT_H
     .                *HC_FSF(HMASS(IH)**2/4.D0/RAUX_H(1)**2)
      SGG(2)= NHC_H(ITOP,IH)*NHC_H(ITOP+1,IH)*V_H/MTMT_H
     .                *HC_FSF(HMASS(IH)**2/4.D0/MTPOLE_H**2)
      SGG(3)=-NHC_H(IT11,IH)*V_H**2/4.D0/STMASS(1)**2
     .                *HC_F0(HMASS(IH)**2/4.D0/STMASS(1)**2)
      SGG(4)=-NHC_H(IT22,IH)*V_H**2/4.D0/STMASS(2)**2
     .                *HC_F0(HMASS(IH)**2/4.D0/STMASS(2)**2)
      SGG(5)=-NHC_H(IB11,IH)*V_H**2/4.D0/SBMASS(1)**2
     .                *HC_F0(HMASS(IH)**2/4.D0/SBMASS(1)**2)
      SGG(6)=-NHC_H(IB22,IH)*V_H**2/4.D0/SBMASS(2)**2
     .                *HC_F0(HMASS(IH)**2/4.D0/SBMASS(2)**2)
      SGG(7)=SGG(1)+SGG(2)+SGG(3)+SGG(4)+SGG(5)+SGG(6)

      PGLUE=0.D0
     .+NHC_H(IBOT,IH)*NHC_H(IBOT+2,IH)*V_H/MBMT_H
     .                *HC_FPF(HMASS(IH)**2/4.D0/RAUX_H(1)**2)
     .+NHC_H(ITOP,IH)*NHC_H(ITOP+2,IH)*V_H/MTMT_H
     .                *HC_FPF(HMASS(IH)**2/4.D0/MTPOLE_H**2)

      PGG(1)= NHC_H(IBOT,IH)*NHC_H(IBOT+2,IH)*V_H/MBMT_H
     .                *HC_FPF(HMASS(IH)**2/4.D0/RAUX_H(1)**2)
      PGG(2)= NHC_H(ITOP,IH)*NHC_H(ITOP+2,IH)*V_H/MTMT_H
     .                *HC_FPF(HMASS(IH)**2/4.D0/MTPOLE_H**2)
      PGG(3)=PGG(1)+PGG(2)
*
      RETURN
      END

      SUBROUTINE HGGSM(HMASS,SGGSM,PGGSM)
************************************************************************
*
* The Standard Model Higgs-glue-glue Coupling
*
* For given Higgs Mass, HMASS=HMASS_SM
*
* SGGSM(1) PGGSM(1) : bottom
* SGGSM(2) PGGSM(2) : top
* SGGSM(3) PGGSM(3) : total
*
* USE THE SAME QUARK MASSES AS IN HGG SUBROUTINES !!!
* See message sent to Sasha on  21 Aug 2011
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
      COMPLEX*16 HC_FSF,HC_FPF
      COMPLEX*16 SGGSM(3),PGGSM(3)
*
      PI=2.D0*DASIN(1.D0)
*
      GB  = GW_H*MBMT_H/2.D0/MW_H
      GSB = 1.D0
      GPB = 0.D0
*
      GT  = GW_H*MTMT_H/2.D0/MW_H
      GST = 1.D0
      GPT = 0.D0
*
      MB_POLE=RAUX_H(1)
*
*Dependence on quark masses in the loop
*
*      print*,'fillcoupl2.f:HGGSM',mbmt_h,mb_s,mb_pole
*      print*,'fillcoupl2.f:HGGSM'
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MBMT_H**2))
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MB_S**2))
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MB_POLE**2))
*      print*,'fillcoupl2.f:HGGSM',mtmt_h,mt_s,mtpole_h
*      print*,'fillcoupl2.f:HGGSM'
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MTMT_H**2))
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MT_S**2))
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MTPOLE_H**2))
*
      SGGSM(1)= GB*GSB*V_H/MBMT_H*HC_FSF(HMASS**2/4.D0/RAUX_H(1)**2)
      SGGSM(2)= GT*GST*V_H/MTMT_H*HC_FSF(HMASS**2/4.D0/MTPOLE_H**2)
      SGGSM(3)=SGGSM(1)+SGGSM(2)
*      print*,'fillcoupl2.f:HGGSM'
*     .      ,cdabs(SGGSM(1)),cdabs(HC_FSF(HMASS**2/4.D0/MB_S**2))
*     .      ,cdabs(SGGSM(2)),cdabs(HC_FSF(HMASS**2/4.D0/MT_S**2))
*     .      ,cdabs(SGGSM(3))
*     .      ,cdabs(HC_FSF(HMASS**2/4.D0/MB_S**2)
*     .            +HC_FSF(HMASS**2/4.D0/MT_S**2))

      PGGSM(1)= GB*GPB*V_H/MBMT_H*HC_FPF(HMASS**2/4.D0/RAUX_H(1)**2)
      PGGSM(2)= GT*GPT*V_H/MTMT_H*HC_FPF(HMASS**2/4.D0/MTPOLE_H**2)
      PGGSM(3)=PGGSM(1)+PGGSM(2)
*      print*,'fillcoupl2.f:HGGSM'
*     .      ,cdabs(PGGSM(3))
*
      RETURN
      END

      COMPLEX*16 FUNCTION HC_FSF(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 HC_FTAU
*
      IF(X.EQ.0.D0) THEN
       HC_FSF=DCMPLX(2.D0/3.D0,0.D0)
       RETURN
      ENDIF
      HC_FSF=1.D0/X*(1.D0+(1.D0-1.D0/X)*HC_FTAU(X))
*
      RETURN
      END

      COMPLEX*16 FUNCTION HC_FPF(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 HC_FTAU
*
      IF(X.EQ.0.D0) THEN
       HC_FPF=DCMPLX(1.D0,0.D0)
       RETURN
      ENDIF
      HC_FPF=HC_FTAU(X)/X
*
      RETURN
      END

      COMPLEX*16 FUNCTION HC_F0(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 HC_FTAU
*
      IF(X.EQ.0.D0) THEN
       HC_F0=DCMPLX(1.D0/3.D0,0.D0)
       RETURN
      ENDIF
      HC_F0=1.D0/X*(-1.D0+HC_FTAU(X)/X)
*
      RETURN
      END

      COMPLEX*16 FUNCTION HC_F1(X)
***********************************************************************
*
***********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMPLEX*16 HC_FTAU
*
      IF(X.EQ.0.D0) THEN
       HC_F1=DCMPLX(7.D0,0.D0)
       RETURN
      ENDIF
      HC_F1=2.D0+3.D0/X+3.D0/X*(2.D0-1.D0/X)*HC_FTAU(X)
*
      RETURN
      END

      COMPLEX*16 FUNCTION HC_FTAU(X)
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
       HC_FTAU=DASIN(DSQRT(X))**2
      ELSEIF(X.GT.1.D0) THEN
       HC_FTAU=-0.25D0*(DLOG( (1.D0+SQRT(1.D0-1.D0/X))
     .                       /(1.D0-SQRT(1.D0-1.D0/X)) ) -XI*PI)**2
      ELSE
*JSL 05/FEB/2004, print removed
*       PRINT*,'INVALID INPUT TO FUNCTION HC_FTAU(X) X = ',X
        RETURN
      ENDIF
*
      RETURN
      END


      SUBROUTINE HSTST(IH,IJ,IK,OMIX,HT,STMIX,COUPLING)
************************************************************************
*
* HiggsIH-stopIJ*-stopIK Coupling
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
*-----------------------------------------------------------------------
      COMPLEX*16 STMIX(2,2),COUPLING
      REAL*8     OMIX(3,3)
      COMPLEX*16 XI,HT,GABC(3,2,2)
*-----------------------------------------------------------------------
*GABC(IA,IB,IC) IA=(phi_1,phi_2,a) IB,IC=(L,R)
      XI=DCMPLX(0.D0,1.D0)
*phi_1
      GABC(1,1,1)=-1.D0/4.D0*(GW_H**2-GP_H**2/3.D0)*V_H*CB_H
      GABC(1,1,2)=1.D0/DSQRT(2.D0)*DCONJG(HT)*MU_H
      GABC(1,2,1)=DCONJG(GABC(1,1,2))
      GABC(1,2,2)=-1.D0/3.D0*GP_H**2*V_H*CB_H
*      print*,'Gamma^phi1(stop):LL,LR,RL,RR'
*      print*,gabc(1,1,1),gabc(1,1,2),gabc(1,2,1),gabc(1,2,2)
*phi_2
      GABC(2,1,1)=-CDABS(HT)**2*V_H*SB_H
     .           +1.D0/4.D0*(GW_H**2-GP_H**2/3.D0)*V_H*SB_H
      GABC(2,1,2)=-1.D0/DSQRT(2.D0)*DCONJG(HT*AT_H)
      GABC(2,2,1)=DCONJG(GABC(2,1,2))
      GABC(2,2,2)=-CDABS(HT)**2*V_H*SB_H+GP_H**2/3.D0*V_H*SB_H
*      print*,'Gamma^phi2(stop):LL,LR,RL,RR'
*      print*,gabc(2,1,1),gabc(2,1,2),gabc(2,2,1),gabc(2,2,2)
*a
      GABC(3,1,1)=0.D0
      GABC(3,1,2)=XI/DSQRT(2.D0)*DCONJG(HT)
     .           *(CB_H*DCONJG(AT_H)+SB_H*MU_H)
      GABC(3,2,1)=DCONJG(GABC(3,1,2))
      GABC(3,2,2)=0.D0
*      print*,'Gamma^a(stop):LL,LR,RL,RR'
*      print*,gabc(3,1,1),gabc(3,1,2),gabc(3,2,1),gabc(3,2,2)
*-----------------------------------------------------------------------
*                  (A IH)             (B,IJ)       (C,IK)     (A,B,C)
      COUPLING=OMIX(1,IH)*DCONJG(STMIX(1,IJ))*STMIX(1,IK)*GABC(1,1,1)
     .        +OMIX(1,IH)*DCONJG(STMIX(1,IJ))*STMIX(2,IK)*GABC(1,1,2)
     .        +OMIX(1,IH)*DCONJG(STMIX(2,IJ))*STMIX(1,IK)*GABC(1,2,1)
     .        +OMIX(1,IH)*DCONJG(STMIX(2,IJ))*STMIX(2,IK)*GABC(1,2,2)
     .        +OMIX(2,IH)*DCONJG(STMIX(1,IJ))*STMIX(1,IK)*GABC(2,1,1)
     .        +OMIX(2,IH)*DCONJG(STMIX(1,IJ))*STMIX(2,IK)*GABC(2,1,2)
     .        +OMIX(2,IH)*DCONJG(STMIX(2,IJ))*STMIX(1,IK)*GABC(2,2,1)
     .        +OMIX(2,IH)*DCONJG(STMIX(2,IJ))*STMIX(2,IK)*GABC(2,2,2)
     .        +OMIX(3,IH)*DCONJG(STMIX(1,IJ))*STMIX(1,IK)*GABC(3,1,1)
     .        +OMIX(3,IH)*DCONJG(STMIX(1,IJ))*STMIX(2,IK)*GABC(3,1,2)
     .        +OMIX(3,IH)*DCONJG(STMIX(2,IJ))*STMIX(1,IK)*GABC(3,2,1)
     .        +OMIX(3,IH)*DCONJG(STMIX(2,IJ))*STMIX(2,IK)*GABC(3,2,2)
      COUPLING=COUPLING/V_H
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE HSBSB(IH,IJ,IK,OMIX,HB,SBMIX,COUPLING)
************************************************************************
*
* HiggsIH-sbottomIJ*-sbottomIK Coupling
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
*-----------------------------------------------------------------------
      COMPLEX*16 SBMIX(2,2),COUPLING
      REAL*8     OMIX(3,3)
      COMPLEX*16 XI,HB,GABC(3,2,2)
*-----------------------------------------------------------------------
*GABC(IA,IB,IC) IA=(phi_1,phi_2,a) IB,IC=(L,R)
      XI=DCMPLX(0.D0,1.D0)
*phi_1
      GABC(1,1,1)=-CDABS(HB)**2*V_H*CB_H
     .           +1.D0/4.D0*(GW_H**2+GP_H**2/3.D0)*V_H*CB_H
      GABC(1,1,2)=-1.D0/DSQRT(2.D0)*DCONJG(HB*AB_H)
      GABC(1,2,1)=DCONJG(GABC(1,1,2))
      GABC(1,2,2)=-CDABS(HB)**2*V_H*CB_H+GP_H**2/6.D0*V_H*CB_H
*      print*,'Gamma^phi1(sbottom):LL,LR,RL,RR'
*      print*,gabc(1,1,1),gabc(1,1,2),gabc(1,2,1),gabc(1,2,2)
*phi_2
      GABC(2,1,1)=-1.D0/4.D0*(GW_H**2+GP_H**2/3.D0)*V_H*SB_H
      GABC(2,1,2)=1.D0/DSQRT(2.D0)*DCONJG(HB)*MU_H
      GABC(2,2,1)=DCONJG(GABC(2,1,2))
      GABC(2,2,2)=-1.D0/6.D0*GP_H**2*V_H*SB_H
*      print*,'Gamma^phi2(sbottom):LL,LR,RL,RR'
*      print*,gabc(2,1,1),gabc(2,1,2),gabc(2,2,1),gabc(2,2,2)
*a
      GABC(3,1,1)=0.D0
      GABC(3,1,2)=XI/DSQRT(2.D0)*DCONJG(HB)
     .           *(SB_H*DCONJG(AB_H)+CB_H*MU_H)
      GABC(3,2,1)=DCONJG(GABC(3,1,2))
      GABC(3,2,2)=0.D0
*      print*,'Gamma^a(sbottom):LL,LR,RL,RR'
*      print*,gabc(3,1,1),gabc(3,1,2),gabc(3,2,1),gabc(3,2,2)
*-----------------------------------------------------------------------
*                  (A IH)             (B,IJ)       (C,IK)     (A,B,C)
      COUPLING=OMIX(1,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(1,IK)*GABC(1,1,1)
     .        +OMIX(1,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(2,IK)*GABC(1,1,2)
     .        +OMIX(1,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(1,IK)*GABC(1,2,1)
     .        +OMIX(1,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(2,IK)*GABC(1,2,2)
     .        +OMIX(2,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(1,IK)*GABC(2,1,1)
     .        +OMIX(2,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(2,IK)*GABC(2,1,2)
     .        +OMIX(2,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(1,IK)*GABC(2,2,1)
     .        +OMIX(2,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(2,IK)*GABC(2,2,2)
     .        +OMIX(3,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(1,IK)*GABC(3,1,1)
     .        +OMIX(3,IH)*DCONJG(SBMIX(1,IJ))*SBMIX(2,IK)*GABC(3,1,2)
     .        +OMIX(3,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(1,IK)*GABC(3,2,1)
     .        +OMIX(3,IH)*DCONJG(SBMIX(2,IJ))*SBMIX(2,IK)*GABC(3,2,2)
      COUPLING=COUPLING/V_H
*-----------------------------------------------------------------------
      RETURN
      END


      SUBROUTINE HSTUSTU(IH,IJ,IK,OMIX,HTAU,STAUMIX,GC)
************************************************************************
*
* HiggsIH-stauIJ*-stauIK Coupling
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
*-----------------------------------------------------------------------
      COMPLEX*16 STAUMIX(2,2),GC
      REAL*8     OMIX(3,3)
      COMPLEX*16 XI,HTAU,GABC(3,2,2)
*-----------------------------------------------------------------------
*GABC(IA,IB,IC) IA=(phi_1,phi_2,a) IB,IC=(L,R)
      XI=DCMPLX(0.D0,1.D0)
*phi_1
      GABC(1,1,1)=-CDABS(HTAU)**2*V_H*CB_H
     .           +1.D0/4.D0*(GW_H**2-GP_H**2)*V_H*CB_H
      GABC(1,1,2)=-1.D0/DSQRT(2.D0)*DCONJG(HTAU*ATAU_H)
      GABC(1,2,1)=DCONJG(GABC(1,1,2))
      GABC(1,2,2)=-CDABS(HTAU)**2*V_H*CB_H+GP_H**2/2.D0*V_H*CB_H
*      print*,'Gamma^phi1(stau):LL,LR,RL,RR'
*      print*,gabc(1,1,1),gabc(1,1,2),gabc(1,2,1),gabc(1,2,2)
*phi_2
      GABC(2,1,1)=-1.D0/4.D0*(GW_H**2-GP_H**2)*V_H*SB_H
      GABC(2,1,2)=1.D0/DSQRT(2.D0)*DCONJG(HTAU)*MU_H
      GABC(2,2,1)=DCONJG(GABC(2,1,2))
      GABC(2,2,2)=-1.D0/2.D0*GP_H**2*V_H*SB_H
*      print*,'Gamma^phi2(stau):LL,LR,RL,RR'
*      print*,gabc(2,1,1),gabc(2,1,2),gabc(2,2,1),gabc(2,2,2)
*a
      GABC(3,1,1)=0.D0
      GABC(3,1,2)=XI/DSQRT(2.D0)*DCONJG(HTAU)
     .           *(SB_H*DCONJG(ATAU_H)+CB_H*MU_H)
      GABC(3,2,1)=DCONJG(GABC(3,1,2))
      GABC(3,2,2)=0.D0
*      print*,'Gamma^a(stau):LL,LR,RL,RR'
*      print*,gabc(3,1,1),gabc(3,1,2),gabc(3,2,1),gabc(3,2,2)
*-----------------------------------------------------------------------
*            (A IH)               (B,IJ)         (C,IK)     (A,B,C)
      GC=OMIX(1,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(1,IK)*GABC(1,1,1)
     .  +OMIX(1,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(2,IK)*GABC(1,1,2)
     .  +OMIX(1,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(1,IK)*GABC(1,2,1)
     .  +OMIX(1,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(2,IK)*GABC(1,2,2)
     .  +OMIX(2,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(1,IK)*GABC(2,1,1)
     .  +OMIX(2,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(2,IK)*GABC(2,1,2)
     .  +OMIX(2,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(1,IK)*GABC(2,2,1)
     .  +OMIX(2,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(2,IK)*GABC(2,2,2)
     .  +OMIX(3,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(1,IK)*GABC(3,1,1)
     .  +OMIX(3,IH)*DCONJG(STAUMIX(1,IJ))*STAUMIX(2,IK)*GABC(3,1,2)
     .  +OMIX(3,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(1,IK)*GABC(3,2,1)
     .  +OMIX(3,IH)*DCONJG(STAUMIX(2,IJ))*STAUMIX(2,IK)*GABC(3,2,2)
      GC=GC/V_H
*-----------------------------------------------------------------------
      RETURN
      END



      SUBROUTINE HNINJ(IN,JN,TW,GP1,GP2,NMIX,GW,GF,GS,GP)
************************************************************************
*
* Higgs-NeutralinoI-NeutralinoJ Coupling
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GP1,GP2,NMIX(4,4),GF,GS,GP
*
      GF=DCMPLX(GW/2.D0,0.D0)
      GS=DCMPLX(1.D0/2.D0
     .  *DREAL( DCONJG(NMIX(JN,2)-TW*NMIX(JN,1))
     .        *(DCONJG(NMIX(IN,3))*GP1-DCONJG(NMIX(IN,4))*GP2)
     .        + DCONJG(NMIX(IN,2)-TW*NMIX(IN,1))
     .        *(DCONJG(NMIX(JN,3))*GP1-DCONJG(NMIX(JN,4))*GP2))
     .  ,0.D0)
      GP=DCMPLX(-1.D0/2.D0
     .  *DIMAG( DCONJG(NMIX(JN,2)-TW*NMIX(JN,1))
     .        *(DCONJG(NMIX(IN,3))*GP1-DCONJG(NMIX(IN,4))*GP2)
     .        + DCONJG(NMIX(IN,2)-TW*NMIX(IN,1))
     .        *(DCONJG(NMIX(JN,3))*GP1-DCONJG(NMIX(JN,4))*GP2))
     .  ,0.D0)
*
      RETURN
      END

      SUBROUTINE HCICJ(IC,JC,GP1,GP2,UL,UR,GW,GF,GS,GP)
************************************************************************
*
* Higgs-CharginoI(+)-CharginoJ(-) Coupling
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 XI,GP1,GP2,UL(2,2),UR(2,2),GF,GS,GP
*
      XI = DCMPLX(0.D0,1.D0)
*
      GF=DCMPLX(GW/DSQRT(2.D0),0.D0)
      GS=1.D0/2.D0*(
     .   UR(IC,1)*DCONJG(UL(JC,2))*GP1+UR(IC,2)*DCONJG(UL(JC,1))*GP2
     .  +DCONJG(
     .   UR(JC,1)*DCONJG(UL(IC,2))*GP1+UR(JC,2)*DCONJG(UL(IC,1))*GP2
     .   ))
      GP=XI/2.D0*(
     .   UR(IC,1)*DCONJG(UL(JC,2))*GP1+UR(IC,2)*DCONJG(UL(JC,1))*GP2
     .  -DCONJG(
     .   UR(JC,1)*DCONJG(UL(IC,2))*GP1+UR(JC,2)*DCONJG(UL(IC,1))*GP2
     .   ))
*
      RETURN
      END

      SUBROUTINE SQMIX(HB_H,HT_H,STMASS,SBMASS,STMIX,SBMIX)
************************************************************************
*
* This subroutine calculates the squark masses and mixing matirices.
* The mixing matrices are parametrized as:
*
*               [1]                        [2]
*  U = [L] / cos_theta                 -sin_theta e^{-i phi} \
*      [R] \ sin_theta e^{+i phi}       cos_theta            / 
*
* where m_2 > m_1 and -pi/2 <= theta,phi <= +pi/2 implying cos(theta)
* and cos(phi) are taken always positive or zeros.
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
*-----------------------------------------------------------------------
* Input
      COMPLEX*16 HB_H,HT_H
*-----------------------------------------------------------------------
* OUTPUT: squark mass and mixining matrices
      REAL*8  STMASS(2),SBMASS(2)
      COMPLEX*16 STMIX(2,2),SBMIX(2,2)
*-----------------------------------------------------------------------
* LOCAL VARIABLES:
      COMPLEX*16 XRL,XLR
      COMPLEX*16 MTMT_R,MBMT_R
*-----------------------------------------------------------------------
* Running Parameters at sfermion mass scales
      PI   = 2.D0*DASIN(1.D0)
      HT   = CDABS(HT_H)
      HB   = CDABS(HB_H)
      BHT  = (9.D0*HT**2/2.D0+HB**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
      BHB  = (9.D0*HB**2/2.D0+HT**2/2.D0-32.D0*PI*ASMT_H)/16.D0/PI**2
*Stop scale
      QT2  = DMAX1(MQ3_H**2+MTPOLE_H**2,MU3_H**2+MTPOLE_H**2)
*Sbottom scale
      QB2  = DMAX1(MQ3_H**2+MBMT_H**2,MD3_H**2+MBMT_H**2)
*Running Top-Yukawa coupling at stop scale
      HTR  = HT*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
*Running Bottom-Yukawa coupling at sbottom scale
      HBR  = HB*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
*----
      XI1T = 1.D0+3.D0*HB**2/(32.D0*PI**2)*DLOG(QT2/MTPOLE_H**2)
      XI1B = 1.D0+3.D0*HB**2/(32.D0*PI**2)*DLOG(QB2/MTPOLE_H**2)
      XI2T = 1.D0+3.D0*HT**2/(32.D0*PI**2)*DLOG(QT2/MTPOLE_H**2)
      XI2B = 1.D0+3.D0*HT**2/(32.D0*PI**2)*DLOG(QB2/MTPOLE_H**2)
      V1   = V_H*CB_H
      V2   = V_H*SB_H
*Running V1 at stop and sbottom scales
      V1RT = V1/XI1T
      V1RB = V1/XI1B
*Running V2 at stop and sbottom scales
      V2RT = V2/XI2T
      V2RB = V2/XI2B
*----
*-----------------------------------------------------------------------
*STOP sector
      V_R    = DSQRT(V1RT**2+V2RT**2)
      SB_R   = DSQRT(V2RT**2/(V1RT**2+V2RT**2))
      CB_R   = DSQRT(V1RT**2/(V1RT**2+V2RT**2))
*     MTMT_R = HTR*V_R*SB_R/DSQRT(2.D0)
      MTMT_R = HT_H*(1.D0+2.D0*BHT*DLOG(QT2/MTPOLE_H**2))**0.25D0
     .             *V_R*SB_R/DSQRT(2.D0)

      XLL = MQ3_H**2+CDABS(MTMT_R)**2
     .     +(4.D0*MW_H**2-MZ_H**2)*(CB_R**2-SB_R**2)/6.D0
      XRR = MU3_H**2+CDABS(MTMT_R)**2
     .     +2.D0*MZ_H**2*SW_H**2*(CB_R**2-SB_R**2)/3.D0
      XRL = MTMT_R*(AT_H-DCONJG(MU_H)*CB_R/SB_R)
      XLR = DCONJG(XRL)

      DELTA = SQRT((XLL-XRR)**2+4.D0*CDABS(XRL)**2)
      XMAVG = (XLL+XRR)/2.D0

      STMASS(1)=SQRT(DABS(XMAVG-DELTA/2.D0))
      IF((XMAVG-DELTA/2.D0).LT.0.D0) STMASS(1)=-STMASS(1)
      STMASS(2)=SQRT(XMAVG+DELTA/2.D0)

      IF(DREAL(XRL).EQ.0.D0) XRL=DCMPLX(1.D-10,DIMAG(XRL))
      PHI=ATAN(DIMAG(XRL)/DREAL(XRL))
      THT_ABS = DATAN(-(STMASS(1)**2-XLL)/CDABS(XRL))
      IF (DREAL(XRL).LT.0.D0) THT= THT_ABS
      IF (DREAL(XRL).GT.0.D0) THT=-THT_ABS
      IF(STMASS(1).GT.0.D0) THEN
      IF(COS(PHI).LT.0.D0 .OR. COS(THT).LT.0.D0 .OR. THT_ABS.LT.0.D0)
     .   THEN
         print*,'ERROR in <SQMIX:stop mixing> !!!'
     .         ,COS(PHI),COS(THT),THT_ABS
      ENDIF
      ENDIF

*     STMIX(alpha,i)
      STMIX(1,1)=DCMPLX(COS(THT),0.D0)
      STMIX(1,2)=DCMPLX(-SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      STMIX(2,1)=DCMPLX(SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      STMIX(2,2)=DCMPLX(COS(THT),0.D0)
*      print*,'stop mass squared in TeV^2:',XLL/1.D6,XLR/1.D6
*     .      ,XRL/1.D6,XRR/1.D6
*Check                  A      AB       B
      M1SQ=DCONJG(STMIX(1,1))*XLL*STMIX(1,1)
     .    +DCONJG(STMIX(1,1))*XLR*STMIX(2,1)
     .    +DCONJG(STMIX(2,1))*XRL*STMIX(1,1)
     .    +DCONJG(STMIX(2,1))*XRR*STMIX(2,1)
      M2SQ=DCONJG(STMIX(1,2))*XLL*STMIX(1,2)
     .    +DCONJG(STMIX(1,2))*XLR*STMIX(2,2)
     .    +DCONJG(STMIX(2,2))*XRL*STMIX(1,2)
     .    +DCONJG(STMIX(2,2))*XRR*STMIX(2,2)
      D121=DCONJG(STMIX(1,1))*XLL*STMIX(1,2)
     .    +DCONJG(STMIX(1,1))*XLR*STMIX(2,2)
     .    +DCONJG(STMIX(2,1))*XRL*STMIX(1,2)
     .    +DCONJG(STMIX(2,1))*XRR*STMIX(2,2)
      D122=DCONJG(STMIX(1,2))*XLL*STMIX(1,1)
     .    +DCONJG(STMIX(1,2))*XLR*STMIX(2,1)
     .    +DCONJG(STMIX(2,2))*XRL*STMIX(1,1)
     .    +DCONJG(STMIX(2,2))*XRR*STMIX(2,1)
*      print*,'Stop Mix : All zer0s?'
*      write(*,4) M1SQ-STMASS(1)**2,M2SQ-STMASS(2)**2,D121,D122
*-----------------------------------------------------------------------
* SBOTTOM sector
      V_R    = DSQRT(V1RB**2+V2RB**2)
      SB_R   = DSQRT(V2RB**2/(V1RB**2+V2RB**2))
      CB_R   = DSQRT(V1RB**2/(V1RB**2+V2RB**2))
*     MBMT_R = HBR*V_R*CB_R/DSQRT(2.D0)
      MBMT_R = HB_H*(1.D0+2.D0*BHB*DLOG(QB2/MTPOLE_H**2))**0.25D0
     .             *V_R*CB_R/DSQRT(2.D0)
*      print*,'SQMIX',hb,bhb,hbr,mbmt_r

      XLL = MQ3_H**2+CDABS(MBMT_R)**2
     .     -(2.D0*MW_H**2+MZ_H**2)*(CB_R**2-SB_R**2)/6.D0
      XRR = MD3_H**2+CDABS(MBMT_R)**2
     .     -MZ_H**2*SW_H**2*(CB_R**2-SB_R**2)/3.D0
      XRL = MBMT_R*(AB_H-DCONJG(MU_H)*SB_R/CB_R)
      XLR = DCONJG(XRL)

      DELTA = SQRT((XLL-XRR)**2+4.D0*CDABS(XRL)**2)
      XMAVG = (XLL+XRR)/2.D0

      SBMASS(1)=SQRT(DABS(XMAVG-DELTA/2.D0))
      IF((XMAVG-DELTA/2.D0).LT.0.D0) SBMASS(1)=-SBMASS(1)
      SBMASS(2)=SQRT(XMAVG+DELTA/2.D0)

      IF(DREAL(XRL).EQ.0.D0) XRL=DCMPLX(1.D-10,DIMAG(XRL))
      PHI=ATAN(DIMAG(XRL)/DREAL(XRL))
      THT_ABS = DATAN(-(SBMASS(1)**2-XLL)/CDABS(XRL))
      IF (DREAL(XRL).LT.0.D0) THT= THT_ABS
      IF (DREAL(XRL).GT.0.D0) THT=-THT_ABS
      IF(SBMASS(1).GT.0.D0) THEN
      IF(COS(PHI).LT.0.D0 .OR. COS(THT).LT.0.D0 .OR. THT_ABS.LT.0.D0)
     .   THEN
         print*,'ERROR in <SQMIX:sbottom mixing> !!!'
     .         ,COS(PHI),COS(THT),THT_ABS
      ENDIF
      ENDIF

      SBMIX(1,1)=DCMPLX(COS(THT),0.D0)
      SBMIX(1,2)=DCMPLX(-SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      SBMIX(2,1)=DCMPLX(SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      SBMIX(2,2)=DCMPLX(COS(THT),0.D0)
*      print*,'sbottom mass squared in TeV^2:',XLL/1.D6,XLR/1.D6
*     .      ,XRL/1.D6,XRR/1.D6
*Check                  A      AB       B
      M1SQ=DCONJG(SBMIX(1,1))*XLL*SBMIX(1,1)
     .    +DCONJG(SBMIX(1,1))*XLR*SBMIX(2,1)
     .    +DCONJG(SBMIX(2,1))*XRL*SBMIX(1,1)
     .    +DCONJG(SBMIX(2,1))*XRR*SBMIX(2,1)
      M2SQ=DCONJG(SBMIX(1,2))*XLL*SBMIX(1,2)
     .    +DCONJG(SBMIX(1,2))*XLR*SBMIX(2,2)
     .    +DCONJG(SBMIX(2,2))*XRL*SBMIX(1,2)
     .    +DCONJG(SBMIX(2,2))*XRR*SBMIX(2,2)
      D121=DCONJG(SBMIX(1,1))*XLL*SBMIX(1,2)
     .    +DCONJG(SBMIX(1,1))*XLR*SBMIX(2,2)
     .    +DCONJG(SBMIX(2,1))*XRL*SBMIX(1,2)
     .    +DCONJG(SBMIX(2,1))*XRR*SBMIX(2,2)
      D122=DCONJG(SBMIX(1,2))*XLL*SBMIX(1,1)
     .    +DCONJG(SBMIX(1,2))*XLR*SBMIX(2,1)
     .    +DCONJG(SBMIX(2,2))*XRL*SBMIX(1,1)
     .    +DCONJG(SBMIX(2,2))*XRR*SBMIX(2,1)
*      print*,'Sbottom Mix : All zer0s?'
*      write(*,4) M1SQ-SBMASS(1)**2,M2SQ-SBMASS(2)**2,D121,D122
*----------------------------------------------------------------
 4    FORMAT(2X,4(1X,E10.4,1X))
*
      RETURN
      END

      SUBROUTINE DUMP_SQ(STMASS,SBMASS,STMIX,SBMIX)
C**************************************************************
C
C**************************************************************
      IMPLICIT REAL*8 (A-H,M,O-Z)
*
      REAL*8  STMASS(2),SBMASS(2)
      COMPLEX*16 STMIX(2,2),SBMIX(2,2)
*-----------------------------------------------------------------------
      PI = 2.D0*DASIN(1.D0)
      PHIT=ATAN(DIMAG(STMIX(2,1))/DREAL(STMIX(2,1)))
      THTT=ASIN(DREAL(STMIX(2,1))/COS(PHIT))
*
      PHIB=ATAN(DIMAG(SBMIX(2,1))/DREAL(SBMIX(2,1)))
      THTB=ASIN(DREAL(SBMIX(2,1))/COS(PHIB))
*
      print*,'---------------------------------------------------------'
      print*,' Masses and Mixing Matrix of Stop and Sbottom : '
      print*,'     STMASS_H(I), STMIX_H(A,I), SBMASS_H(I), SBMIX_H(A,I)'
      print*,'---------------------------------------------------------'
      DO IS=1,2
      WRITE(*,1) IS,STMASS(IS)
      ENDDO
      print*,' U[Stop] = '
      print*,'                [1]                     [2]'
      WRITE(*,3) DREAL(STMIX(1,1)),DIMAG(STMIX(1,1))
     .          ,DREAL(STMIX(1,2)),DIMAG(STMIX(1,2))
      WRITE(*,4) DREAL(STMIX(2,1)),DIMAG(STMIX(2,1))
     .          ,DREAL(STMIX(2,2)),DIMAG(STMIX(2,2))
      WRITE(*,5),THTT/PI*180.D0,PHIT/PI*180.D0
      print*,' '
      DO IS=1,2
      WRITE(*,2) IS,SBMASS(IS)
      ENDDO
      print*,' U[Sbottom] = '
      print*,'                [1]                     [2]'
      WRITE(*,3) DREAL(SBMIX(1,1)),DIMAG(SBMIX(1,1))
     .          ,DREAL(SBMIX(1,2)),DIMAG(SBMIX(1,2))
      WRITE(*,4) DREAL(SBMIX(2,1)),DIMAG(SBMIX(2,1))
     .          ,DREAL(SBMIX(2,2)),DIMAG(SBMIX(2,2))
      WRITE(*,5),THTB/PI*180.D0,PHIB/PI*180.D0
      print*,'---------------------------------------------------------'
*-----------------------------------------------------------------------
 1    FORMAT(2X,'Mass of Stop(',I1,')    = ',E10.4,' GeV')
 2    FORMAT(2X,'Mass of Sbottom(',I1,') = ',E10.4,' GeV')
 3    FORMAT(2X,'[L] [','(',E10.4,1X,E10.4,') '
     .                 ,'(',E10.4,1X,E10.4,')',' ]')
 4    FORMAT(2X,'[R] [','(',E10.4,1X,E10.4,') '
     .                  ,'(',E10.4,1X,E10.4,')',' ]')
 5    FORMAT(2X,'Theta = ',E10.4,' Deg.  :  Phi = ',E10.4,' Deg. ')
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SLMIX(HB_H,HT_H,HTAU_H,STAUMASS,SNU3MASS,STAUMIX)
************************************************************************
*
* This subroutine calculates the slepton masses and mixing matirix.
* The mixing matrice is parametrized as:
*
*               [1]                        [2]
*  U = [L] / cos_theta                 -sin_theta e^{-i phi} \
*      [R] \ sin_theta e^{+i phi}       cos_theta            / 
*
* where m_2 > m_1 and -pi/2 <= theta,phi <= +pi/2 implying cos(theta)
* and cos(phi) are taken always positive or zeros.
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
*-----------------------------------------------------------------------
* Input
      COMPLEX*16 HB_H,HT_H,HTAU_H
*-----------------------------------------------------------------------
* OUTPUT: squark mass and mixining matrices
      REAL*8  STAUMASS(2),SNU3MASS
      COMPLEX*16 STAUMIX(2,2)
*-----------------------------------------------------------------------
* LOCAL VARIABLES:
      COMPLEX*16 XRL,XLR
      COMPLEX*16 MTAU_R
*-----------------------------------------------------------------------
* Running Parameters at sfermion mass scales
      PI   = 2.D0*DASIN(1.D0)
      HT   = CDABS(HT_H)
      HB   = CDABS(HB_H)
      HTAU = CDABS(HTAU_H)
*Slepton scale
      QL2  = DMAX1(ML3_H**2,ME3_H**2)
*----
      XI1 = 1.D0+3.D0*HB**2/(32.D0*PI**2)*DLOG(QL2/MTPOLE_H**2)
      XI2 = 1.D0+3.D0*HT**2/(32.D0*PI**2)*DLOG(QL2/MTPOLE_H**2)
      V1  = V_H*CB_H
      V2  = V_H*SB_H
*Running V1 and V2 at slepton scale
      V1R = V1/XI1
      V2R = V2/XI2
*----
*-----------------------------------------------------------------------
* SNU3MASS
      SNU3MSQ=ML3_H**2+(CB_H**2-SB_H**2)*MZ_H**2/2.D0
      SNU3MASS=DSQRT(DABS(SNU3MSQ))
      IF (SNU3MSQ.LT.0.D0) SNU3MASS=-SNU3MASS
*-----------------------------------------------------------------------
* STAU sector
      V_R    = DSQRT(V1R**2+V2R**2)
      SB_R   = DSQRT(V2R**2/(V1R**2+V2R**2))
      CB_R   = DSQRT(V1R**2/(V1R**2+V2R**2))
      MTAU_R = HTAU_H*V_R*CB_R/DSQRT(2.D0)
*      print*,'SLMIX',htau,mtau_r

      XLL = ML3_H**2+CDABS(MTAU_R)**2
     .     +MZ_H**2*(CB_R**2-SB_R**2)*(SW_H**2-1.D0/2.D0)
      XRR = ME3_H**2+CDABS(MTAU_R)**2
     .     -MZ_H**2*SW_H**2*(CB_R**2-SB_R**2)
      XRL = MTAU_R*(ATAU_H-DCONJG(MU_H)*SB_R/CB_R)
      XLR = DCONJG(XRL)

      DELTA = SQRT((XLL-XRR)**2+4.D0*CDABS(XRL)**2)
      XMAVG = (XLL+XRR)/2.D0

      STAUMASS(1)=SQRT(DABS(XMAVG-DELTA/2.D0))
      IF((XMAVG-DELTA/2.D0).LT.0.D0) STAUMASS(1)=-STAUMASS(1)
      STAUMASS(2)=SQRT(XMAVG+DELTA/2.D0)

      IF(DREAL(XRL).EQ.0.D0) XRL=DCMPLX(1.D-10,DIMAG(XRL))
      PHI=ATAN(DIMAG(XRL)/DREAL(XRL))
      THT_ABS = DATAN(-(STAUMASS(1)**2-XLL)/CDABS(XRL))
      IF (DREAL(XRL).LT.0.D0) THT= THT_ABS
      IF (DREAL(XRL).GT.0.D0) THT=-THT_ABS
      IF(STAUMASS(1).GT.0.D0) THEN
      IF(COS(PHI).LT.0.D0 .OR. COS(THT).LT.0.D0 .OR. THT_ABS.LT.0.D0)
     .   THEN
         print*,'ERROR in <SLMIX:stau mixing> !!!'
     .         ,COS(PHI),COS(THT),THT_ABS
      ENDIF
      ENDIF

      STAUMIX(1,1)=DCMPLX(COS(THT),0.D0)
      STAUMIX(1,2)=DCMPLX(-SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      STAUMIX(2,1)=DCMPLX(SIN(THT)*COS(PHI),SIN(THT)*SIN(PHI))
      STAUMIX(2,2)=DCMPLX(COS(THT),0.D0)
*      print*,'stau mass squared in TeV^2:',XLL/1.D6,XLR/1.D6
*     .      ,XRL/1.D6,XRR/1.D6
*Check                    A      AB         B
      M1SL=DCONJG(STAUMIX(1,1))*XLL*STAUMIX(1,1)
     .    +DCONJG(STAUMIX(1,1))*XLR*STAUMIX(2,1)
     .    +DCONJG(STAUMIX(2,1))*XRL*STAUMIX(1,1)
     .    +DCONJG(STAUMIX(2,1))*XRR*STAUMIX(2,1)
      M2SL=DCONJG(STAUMIX(1,2))*XLL*STAUMIX(1,2)
     .    +DCONJG(STAUMIX(1,2))*XLR*STAUMIX(2,2)
     .    +DCONJG(STAUMIX(2,2))*XRL*STAUMIX(1,2)
     .    +DCONJG(STAUMIX(2,2))*XRR*STAUMIX(2,2)
      D121=DCONJG(STAUMIX(1,1))*XLL*STAUMIX(1,2)
     .    +DCONJG(STAUMIX(1,1))*XLR*STAUMIX(2,2)
     .    +DCONJG(STAUMIX(2,1))*XRL*STAUMIX(1,2)
     .    +DCONJG(STAUMIX(2,1))*XRR*STAUMIX(2,2)
      D122=DCONJG(STAUMIX(1,2))*XLL*STAUMIX(1,1)
     .    +DCONJG(STAUMIX(1,2))*XLR*STAUMIX(2,1)
     .    +DCONJG(STAUMIX(2,2))*XRL*STAUMIX(1,1)
     .    +DCONJG(STAUMIX(2,2))*XRR*STAUMIX(2,1)
*      print*,'Stau Mix : All zer0s?'
*      write(*,4) M1SL-STAUMASS(1)**2,M2SL-STAUMASS(2)**2,D121,D122
*----------------------------------------------------------------
 4    FORMAT(2X,4(1X,E10.4,1X))
*
      RETURN
      END

      SUBROUTINE DUMP_SL(STAUMASS,SNU3MASS,STAUMIX)
C**************************************************************
C
C**************************************************************
      IMPLICIT REAL*8 (A-H,M,O-Z)
*
      REAL*8  STAUMASS(2),SNU3MASS
      COMPLEX*16 STAUMIX(2,2)
*-----------------------------------------------------------------------
      PI = 2.D0*DASIN(1.D0)
*
      PHI=ATAN(DIMAG(STAUMIX(2,1))/DREAL(STAUMIX(2,1)))
      THT=ASIN(DREAL(STAUMIX(2,1))/COS(PHI))
*
*      print*,'---------------------------------------------------------'
      print*,' Masses and Mixing Matrix of Sneutrino and Stau : '
      print*,'                SNU3MASS_H, STAUMASS_H(I), STAUMIX_H(A,I)'
      print*,'---------------------------------------------------------'
      WRITE(*,1) SNU3MASS
      print*,' '
      DO IS=1,2
      WRITE(*,2) IS,STAUMASS(IS)
      ENDDO
      print*,' U[Stau] = '
      print*,'                [1]                     [2]'
      WRITE(*,3) DREAL(STAUMIX(1,1)),DIMAG(STAUMIX(1,1))
     .          ,DREAL(STAUMIX(1,2)),DIMAG(STAUMIX(1,2))
      WRITE(*,4) DREAL(STAUMIX(2,1)),DIMAG(STAUMIX(2,1))
     .          ,DREAL(STAUMIX(2,2)),DIMAG(STAUMIX(2,2))
      WRITE(*,5),THT/PI*180.D0,PHI/PI*180.D0
      print*,'---------------------------------------------------------'
*-----------------------------------------------------------------------
 1    FORMAT(2X,'Mass of Sneutrino3 = ',E10.4,' GeV')
 2    FORMAT(2X,'Mass of Stau(',I1,')    = ',E10.4,' GeV')
 3    FORMAT(2X,'[L] [','(',E10.4,1X,E10.4,') '
     .                 ,'(',E10.4,1X,E10.4,')',' ]')
 4    FORMAT(2X,'[R] [','(',E10.4,1X,E10.4,') '
     .                  ,'(',E10.4,1X,E10.4,')',' ]')
 5    FORMAT(2X,'Theta = ',E10.4,' Deg.  :  Phi = ',E10.4,' Deg. ')
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE DUMP_CN(TANB,M1,M2,MU,TH1,TH2,THMU
     .                  ,M_C,M_N,U_L,U_R,N_N)
C**************************************************************
C
C**************************************************************
      IMPLICIT REAL*8 (A-H,M,O-Z)
*
      REAL*8  M_C(2),M_N(4)
      COMPLEX*16 U_L(2,2),U_R(2,2),N_N(4,4)
      COMPLEX*16 V_C(2,2)
*
      PI = 2.D0*DASIN(1.D0)
      V_C(1,1)=DCONJG(U_R(1,1))
      V_C(1,2)=DCONJG(U_R(1,2))
      V_C(2,1)=DCONJG(U_R(2,1))
      V_C(2,2)=DCONJG(U_R(2,2))
*
*      print*,'======================================================='
*      print*,'TANB         M1           M2           MU'
*      WRITE(*,4) TANB,M1,M2,MU
*      print*,' '
*      print*,'TH1          TH2          THMU   : in Degree'
*      WRITE(*,3) TH1*180.D0/PI,TH2*180.D0/PI,THMU*180.D0/PI
      print*,'---------------------------------------------------------'
      print*,' Chargino Masses and Mixing Matrices : '
      print*,'                        MC_H(I), UL_H(I,A), and UR_H(I,A)'
      print*,'---------------------------------------------------------'
      WRITE(*,1) M_C(1),M_C(2)
      print*,' '
      print*,' UL_H ='
      WRITE(*,8) DREAL(U_L(1,1)),DIMAG(U_L(1,1))
     .          ,DREAL(U_L(1,2)),DIMAG(U_L(1,2))
      WRITE(*,9) DREAL(U_L(2,1)),DIMAG(U_L(2,1))
     .          ,DREAL(U_L(2,2)),DIMAG(U_L(2,2))
      print*,' '
      print*,' UR_H ='
      WRITE(*,8) DREAL(U_R(1,1)),DIMAG(U_R(1,1))
     .          ,DREAL(U_R(1,2)),DIMAG(U_R(1,2))
      WRITE(*,9) DREAL(U_R(2,1)),DIMAG(U_R(2,1))
     .          ,DREAL(U_R(2,2)),DIMAG(U_R(2,2))
      print*,'---------------------------------------------------------'
      print*,' Neutralino Masses MN_H(I) and Mixing Matrix N_H(I,A)'
      print*,'---------------------------------------------------------'
      WRITE(*,2) M_N(1),M_N(2)
      WRITE(*,3) M_N(3),M_N(4)
      print*,' '
      WRITE(*,11) DREAL(N_N(1,1)),DIMAG(N_N(1,1))
      WRITE(*,12) DREAL(N_N(1,2)),DIMAG(N_N(1,2))
      WRITE(*,13) DREAL(N_N(1,3)),DIMAG(N_N(1,3))
      WRITE(*,14) DREAL(N_N(1,4)),DIMAG(N_N(1,4))
      print*,' '
      WRITE(*,21) DREAL(N_N(2,1)),DIMAG(N_N(2,1))
      WRITE(*,22) DREAL(N_N(2,2)),DIMAG(N_N(2,2))
      WRITE(*,23) DREAL(N_N(2,3)),DIMAG(N_N(2,3))
      WRITE(*,24) DREAL(N_N(2,4)),DIMAG(N_N(2,4))
      print*,' '
      WRITE(*,31) DREAL(N_N(3,1)),DIMAG(N_N(3,1))
      WRITE(*,32) DREAL(N_N(3,2)),DIMAG(N_N(3,2))
      WRITE(*,33) DREAL(N_N(3,3)),DIMAG(N_N(3,3))
      WRITE(*,34) DREAL(N_N(3,4)),DIMAG(N_N(3,4))
      print*,' '
      WRITE(*,41) DREAL(N_N(4,1)),DIMAG(N_N(4,1))
      WRITE(*,42) DREAL(N_N(4,2)),DIMAG(N_N(4,2))
      WRITE(*,43) DREAL(N_N(4,3)),DIMAG(N_N(4,3))
      WRITE(*,44) DREAL(N_N(4,4)),DIMAG(N_N(4,4))
      print*,'---------------------------------------------------------'
*      print*,' ======================================================='
*
  1   FORMAT(2X,'MC1 = ',E10.4,' GeV',6X,'MC2 = ',E10.4,' GeV')
  2   FORMAT(2X,'MN1 = ',E10.4,' GeV',6X,'MN2 = ',E10.4,' GeV')
  3   FORMAT(2X,'MN3 = ',E10.4,' GeV',6X,'MN4 = ',E10.4,' GeV')
  8   FORMAT(2X,' [','(',E10.4,1X,E10.4,') '
     .              ,'(',E10.4,1X,E10.4,')',' ]')
  9   FORMAT(2X,' [','(',E10.4,1X,E10.4,') '
     .               ,'(',E10.4,1X,E10.4,')',' ]')
 11   FORMAT(2X,'N_H(1,1) = ','(',E10.4,1X,E10.4,') ')
 12   FORMAT(2X,'N_H(1,2) = ','(',E10.4,1X,E10.4,') ')
 13   FORMAT(2X,'N_H(1,3) = ','(',E10.4,1X,E10.4,') ')
 14   FORMAT(2X,'N_H(1,4) = ','(',E10.4,1X,E10.4,') ')
 21   FORMAT(2X,'N_H(2,1) = ','(',E10.4,1X,E10.4,') ')
 22   FORMAT(2X,'N_H(2,2) = ','(',E10.4,1X,E10.4,') ')
 23   FORMAT(2X,'N_H(2,3) = ','(',E10.4,1X,E10.4,') ')
 24   FORMAT(2X,'N_H(2,4) = ','(',E10.4,1X,E10.4,') ')
 31   FORMAT(2X,'N_H(3,1) = ','(',E10.4,1X,E10.4,') ')
 32   FORMAT(2X,'N_H(3,2) = ','(',E10.4,1X,E10.4,') ')
 33   FORMAT(2X,'N_H(3,3) = ','(',E10.4,1X,E10.4,') ')
 34   FORMAT(2X,'N_H(3,4) = ','(',E10.4,1X,E10.4,') ')
 41   FORMAT(2X,'N_H(4,1) = ','(',E10.4,1X,E10.4,') ')
 42   FORMAT(2X,'N_H(4,2) = ','(',E10.4,1X,E10.4,') ')
 43   FORMAT(2X,'N_H(4,3) = ','(',E10.4,1X,E10.4,') ')
 44   FORMAT(2X,'N_H(4,4) = ','(',E10.4,1X,E10.4,') ')
*
      RETURN
      END

      SUBROUTINE CHARDIAG(M2,MU,TANB,TH2,THMU,MC,U,V)
C******************************************************************
C* Diagonalizes the chargino mass matrix of the MSSM in Gunion    *
C* and Haber notation, allowing the gaugino mass M_2 and higgsino *
C* mass mu to be complex. M2 and MU are the absolute values of    *
C* these quantities, and TH2 and THMU their phases. TANB is the   *
C* ratio of vevs tan(beta). U and V are the diagonalization ma-   *
C* trices, and MC(1) and MC(2) the two eigenvalues, 1 standing    *
C* for the lighter one.                                           *
C******************************************************************


      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 U(2,2), V(2,2), CHECK3(2,2), X(2,2), XI
      DIMENSION MC(2)
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN

C  *** Eigenvalues ***

      MW   = MW_CN
      MZ   = MZ_CN

      M2SQ = M2*M2
      MUSQ = MU*MU
      MWSQ = MW*MW

      THETA = TH2 + THMU
      BETA = DATAN(TANB)
      CB = DCOS(BETA)
      SB = DSIN(BETA)

      TERM1 = M2SQ + MUSQ + 2.D0*MWSQ
      TERM2 = DSQRT( (M2SQ - MUSQ - 2.D0*MWSQ*DCOS(2.D0*BETA))**2
     &             + 8.D0*MWSQ*( (M2*CB)**2 + (MU*SB)**2
     &                         + M2*MU*DSIN(2.D0*BETA)*DCOS(THETA) ) )

      MC1SQ = .5D0*(TERM1 - TERM2)
      MC2SQ = .5D0*(TERM1 + TERM2)

      IF(MC1SQ.LT.0.D0) THEN
         WRITE(*,*) ' Squared chargino mass is negative!!'
         STOP 77
      ENDIF

      MC(1) = DSQRT(MC1SQ)
      MC(2) = DSQRT(MC2SQ)

C  *** Phases Gamma_L, Gamma_R  ***

      C2 = DCOS(TH2)
      S2 = DSIN(TH2)
      CMU = DCOS(THMU)
      SMU = DSIN(THMU)

      XR = CB*C2*M2 + SB*CMU*MU
      IF(DABS(XR).LT.1.D-5) XR = 1.D-5
      GAMR = DATAN( (CB*S2*M2 - SB*SMU*MU) / XR )

      XL = SB*C2*M2 + CB*CMU*MU
      IF(DABS(XL).LT.1.D-5) XL = 1.D-5
      GAML = DATAN( (SB*S2*M2 - CB*SMU*MU) / XL )

C  *** Mixing angle theta_L, theta_R  ***

      RT2 = DSQRT(2.D0)
      YR = M2*CB*DCOS(GAMR-TH2) + MU*SB*DCOS(GAMR+THMU)
      IF(DABS(YR).LT.1.D-5) YR = 1.D-5
      THR = DATAN( ( M2SQ + 2.D0*(MW*SB)**2 - MC1SQ )
     &           / (RT2*MW*YR)                       )

      YL = M2*SB*DCOS(GAML-TH2) + MU*CB*DCOS(GAML+THMU)
      IF(DABS(YL).LT.1.D-5) YL = 1.D-5
      THL = DATAN( ( M2SQ + 2.D0*(MW*CB)**2 - MC1SQ )
     &           / (RT2*MW*YL)                       )


C  *** Phases gamma_1, gamma_2     ***

      CTHL = DCOS(THL)
      CTHR = DCOS(THR)
      STHL = DSIN(THL)
      STHR = DSIN(THR)
      IF(DABS(CTHL).LT.1.D-7) CTHL = 1.D-7
      CGR = DCOS(GAMR)
      CGL = DCOS(GAML)
      SGR = DSIN(GAMR)
      SGL = DSIN(GAML)

*JSL:06/Jun/06: CG1 and CG2 are in one line (Thanks to Pukhov)
      CG1 = ( M2*C2*CTHR - RT2*MW*CB*CGR*STHR) / (MC(1)*CTHL)
      IF(DABS(CG1).GT.1.0001D0) THEN
        WRITE(*,*) ' |cos(theta_1)| > 1!!'
        GAM1 = 0.D0
      ELSEIF(CG1.GT.1.D0) THEN
        GAM1 = 0.D0
      ELSEIF(CG1.LT.-1.D0) THEN
        GAM1 = 3.141592654D0
      ELSEIF(DABS(1.D0-CG1).LT.1.D-12) THEN
        GAM1 = 0.D0
      ELSE
        GAM1 = DACOS(CG1)
      ENDIF
      SG1 = ( M2*S2*CTHR - RT2*MW*CB*SGR*STHR ) / (MC(1)*CTHL)
      IF(SG1.LT.0.D0) GAM1 = -GAM1

      CG2 = ( MU*CMU*CTHR + RT2*MW*SB*CGR*STHR) / (MC(2)*CTHL)
      IF(DABS(CG2).GT.1.0001D0) THEN
        WRITE(*,*) ' |cos(gamma_2)| > 1!!'
        GAM2 = 0.D0
      ELSEIF(CG2.GT.1.D0) THEN
        GAM2 = 0.D0
      ELSEIF(CG2.LT.-1.D0) THEN
        GAM2 = 3.141592654D0
      ELSEIF(DABS(1.D0-CG2).LT.1.D-12) THEN
        GAM2 = 0.D0
      ELSE
        GAM2 = DACOS(CG2)
      ENDIF
      SG2 = ( MU*SMU*CTHR - RT2*MW*SB*SGR*STHR ) / (MC(2)*CTHL)
      IF(SG2.LT.0.D0) GAM2 = -GAM2

C  *** Define diagonalizing matrices U, V  ***

      XI = (0.D0,1.D0)
      
      U(1,1) = CTHR
      U(1,2) = -STHR*( CGR - XI*SGR )
      U(2,1) = STHR*( CGR + XI*SGR )
      U(2,2) = CTHR

      V(1,1) = CTHL*( CG1 + XI*SG1 )
      V(1,2) = -STHL*( DCOS(GAM1-GAML) + XI*DSIN(GAM1-GAML) )
      V(2,1) = STHL*( DCOS(GAML+GAM2) + XI*DSIN(GAML+GAM2) )
      V(2,2) = CTHL*( CG2 + XI*SG2 )
*      print*,'JSLEE:1-CG1,1-CG2',1.D0-CG1,1.D0-CG2
*      print*,'JSLEE:GAMR,GAML,GAM1,GAM2',GAMR,GAML,GAM1,GAM2
*      print*,'JSLEE:V(1,2)',V(1,2)

C  *** Checks  ***

      CHECK1 = CG1*CG1 + SG1*SG1
      IF(DABS(CHECK1-1.D0).GT.1.D-5) THEN
        WRITE(*,*) ' gamma_1 not computed correctly!'
        WRITE(1,22) GAML,GAMR,THL,THR,GAM1,GAM2
      ENDIF

      CHECK2 = CG2*CG2 + SG2*SG2
      IF(DABS(CHECK2-1.D0).GT.1.D-5) THEN
        WRITE(*,*) ' gamma_2 not computed correctly!'
        WRITE(1,22) GAML,GAMR,THL,THR,GAM1,GAM2
      ENDIF
 22   FORMAT(' Gamma_L, Gamma_R = ',2(e11.4,2x),/,
     &       ' theta_L, theta_R = ',2(e11.4,2x),/,
     &       ' gamma_1, gamma_2 = ',2(e11.4,2x))

      X(1,1) = M2*( C2 + XI*S2 )
      X(2,2) = MU*( CMU + XI*SMU )
      X(1,2) = RT2*MW*SB
      X(2,1) = RT2*MW*CB

      DO 10 I = 1,2
      DO 10 L = 1,2
      CHECK3(I,L) = (0.D0,0.D0)

      DO 20 J = 1,2
      DO 20 K = 1,2
 20   CHECK3(I,L) = CHECK3(I,L) 
     &            + DCONJG(U(I,J))*X(J,K)*DCONJG(V(L,K))
 10   CONTINUE

      IF(ABS(CHECK3(1,2)).GT.1.D-3.OR.
     &   ABS(CHECK3(2,1)).GT.1.D-3.OR.
     &   DABS(DIMAG(CHECK3(1,1))).GT.1.D-3.OR.
     &   DABS(DIMAG(CHECK3(2,2))).GT.1.D-3)
     & WRITE(*,*) ' Diagonalization of chargino mass matrix failed!!'

c      WRITE(1,*) '       Diagonalized chargino mass matrix:'
c      DO 30 I = 1,2
c 30   WRITE(1,31) CHECK3(I,1), CHECK3(I,2)
c 31   FORMAT(2(2X,E11.4,' + i*',e11.4))

      RETURN
      END


      SUBROUTINE NEUTDIAG(M1,M2,MU,TANB,TH1,TH2,THMU,MN,N)
C********************************************************************
C* Diagonalizes the complex, symmetric MSSM neutralino mass matrix. *
C* Gaugino mass unification is assumed.                             *
C********************************************************************

      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 Y(4,4),N(4,4),MDIAG(4,4),XI
      DIMENSION AUX(8,8),EV(8),H(8),MN(4)
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN

      MW   = MW_CN
      MZ   = MZ_CN
      CW   = DSQRT(1.D0-S2W_CN)
      SW   = DSQRT(S2W_CN)

c     M1 = 5.D0*S2W_CN*M2/(3.D0*CW*CW)     !U(1) gaugino mass

*      print*,m1,m2,mu,tanb,th1,th2,thmu

C  *** Define complex neutralino mass matrix  ***

      XI = (0.D0,1.D0)

      C1 = DCOS(TH1)
      S1 = DSIN(TH1)

      C2 = DCOS(TH2)
      S2 = DSIN(TH2)

      CMU = DCOS(THMU)
      SMU = DSIN(THMU)

      BETA = DATAN(TANB)
      CB   = DCOS(BETA)
      SB   = DSIN(BETA)

      Y(1,1) = M1*(C1+XI*S1)
      Y(1,2) = (0.D0,0.D0)
      Y(1,3) =-MZ*SW*CB
      Y(1,4) = MZ*SW*SB
      Y(2,2) = M2*(C2+XI*S2)
      Y(2,3) = MW*CB
      Y(2,4) =-MW*SB
      Y(3,3) = 0.D0
      Y(3,4) =-MU*(CMU+XI*SMU)
      Y(4,4) = 0.D0

      DO 1 I = 2, 4
      DO 1 J = 1, I-1
 1    Y(I,J) = Y(J,I)


C  *** Define auxiliary real, symmetric 8x8 matrix AUX ***

      DO 10 I = 1,4
      DO 10 J = 1,4
      AUX(I,J)     = DREAL(Y(I,J))
      AUX(I,J+4)   = DIMAG(Y(I,J))
      AUX(I+4,J)   = DIMAG(Y(I,J))
 10   AUX(I+4,J+4) =-DREAL(Y(I,J))


C  *** Diagonalize AUX; eigenvalues in EV, eigenvectors in AUX  ***

      CALL DIAGRS(8,8,AUX,EV,H,IERR)

c      WRITE(*,11) (EV(K),K=1,8)
c 11   FORMAT(' Eigenvalues of AUX:',/2(4(2x,e11.4),/))
c      WRITE(*,*) '                   Eigenvectors of AUX:'
c      DO 12 I = 1,8
c 12   WRITE(*,13) (AUX(K,I),K=1,8)
c 13   FORMAT(8(1X,F7.4))

      DO 20 I = 1,4
 20   MN(I) = EV(I+4)        !First 4 eigenvalues are negative!


C  *** Define N. Recall that the eigenvectors are the COLUMNS of AUX! ***
c
c Re(N)[i,a]=AUX^T[i,a] and Im(N)[i,a]=AUX^T[i,4+a] with i=5-8 for
c the positive mass eigenstates and a=1-4 for the electroweak states

      DO 30 I = 1,4
      DO 30 J = 1,4

 30   N(I,J) = AUX(J,4+I) + XI*AUX(J+4,4+I)


C  *** Check diagonalization  ***

      DO 40 I = 1,4
      DO 40 L = 1,4
      MDIAG(I,L) = (0.D0,0.D0)

      DO 41 J = 1,4
      DO 41 K = 1,4
 41   MDIAG(I,L) = MDIAG(I,L)
     .           + DCONJG(N(I,J))*Y(J,K)*DCONJG(N(L,K))
c      print*,'>> NEUTDIAG << MN(',i,l,')=',mdiag(i,l)
 40   CONTINUE

      RETURN
      END


      SUBROUTINE DIAGRS(NM,N,Z,D,E,IERR)
C*************************************************************
C*  NM: MAXIMAL DIMENSION OF THE MATRIX Z                    *
C*  N : ACTUAL DIMENSION IN THE CALLING PROGRAM              *
C*  E : AUXILIARY VECTOR                                     *
C*  D : VECTOR CONTAINING THE EIGENVALUES                    *
C*  AFTER THE DIAGONALIZATION, THE COLUMNS -                 *
C*  NOT THE ROWS!! - OF THE MATRIX                           *
C*  Z ARE THE EIGENVECTORS                                   *
C*  IERR: ERROR PARAMETER (IERR=0: EVERYTHING OK!)           *
C*************************************************************

C**** EISPACK TRED2
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(NM,N),E(N),D(N)
      IF(N.EQ.1) GOTO 320
      DO 300 II=2,N
      I=N+2-II
      L=I-1
      H=0.D0
      SCALE=0.D0
      IF(L.LT.2) GOTO 130
      DO 120 K=1,L
  120 SCALE=SCALE+DABS(Z(I,K))
      IF(SCALE.NE.0.D0) GOTO 140
  130 E(I)=Z(I,L)
      GOTO 290
  140 DO 150 K=1,L
      Z(I,K)=Z(I,K)/SCALE
      H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
      F=Z(I,L)
      G=-DSIGN(DSQRT(H),F)
      E(I)=SCALE*G
      H=H-F*G
      Z(I,L)=F-G
      F=0.D0
      DO 240 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.D0
      DO 180 K=1,J
  180 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF(L.LT.JP1) GOTO 220
      DO 200 K=JP1,L
  200 G=G+Z(K,J)*Z(I,K)
  220 E(J)=G/H
      F=F+E(J)*Z(I,J)
  240 CONTINUE
      HH=F/(H+H)
      DO 260 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
      DO 260 K=1,J
      Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
  290 D(I)=H
  300 CONTINUE
  320 D(1)=0.D0
      E(1)=0.D0
      DO 500 I=1,N
      L=I-1
      IF(D(I).EQ.0.D0) GOTO 380
      DO 360 J=1,L
      G=0.D0
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
      Z(I,I)=1.D0
      IF(L.LT.1) GOTO 500
      DO 400 J=1,L
      Z(I,J)=0.D0
      Z(J,I)=0.D0
  400 CONTINUE
  500 CONTINUE

C**** EISPACK IMTQL2

      GENAU = 2.D0  **(-40.D0)
      IERR = 0
      IF (N .EQ. 1) GO TO 5001
      DO 5100 I = 2, N
 5100 E(I-1) = E(I)
      E(N) = 0.D0
      DO 5240 L = 1, N
         J = 0
 5105    DO 5110 M = L, N
            IF (M .EQ. N) GO TO 5120
            IF ( DABS(E(M)) .LE. GENAU * ( DABS(D(M)) +  DABS(D(M+1))))
     X         GO TO 5120
 5110    CONTINUE
 5120    P = D(L)
         IF (M .EQ. L) GO TO 5240
         IF (J .EQ. 30) GO TO 5000
         J = J + 1
         G = (D(L+1) - P) / (2.D0   * E(L))
         R =  DSQRT(G*G+1.D0  )
         G = D(M) - P + E(L) / (G +  DSIGN(R,G))
         S = 1.D0
         C = 1.D0
         P = 0.D0
         MML = M - L
         DO 5200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF ( DABS(F) .LT.  DABS(G)) GO TO 5150
            C = G / F
            R =  DSQRT(C*C+1.D0  )
            E(I+1) = F * R
            S = 1.D0   / R
            C = C * S
            GO TO 5160
 5150       S = F / G
            R =  DSQRT(S*S+1.D0  )
            E(I+1) = G * R
            C = 1.D0   / R
            S = S * C
 5160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.D0   * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
            DO 5180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
 5180       CONTINUE
 5200    CONTINUE
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.D0
         GO TO 5105
 5240 CONTINUE
      DO 5300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 5260 J = II, N
            IF (D(J) .GE. P) GO TO 5260
            K = J
            P = D(J)
 5260    CONTINUE
         IF (K .EQ. I) GO TO 5300
         D(K) = D(I)
         D(I) = P
         DO 5280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
 5280    CONTINUE
 5300 CONTINUE
      GO TO 5001
 5000 IERR = L
 5001 RETURN
      END


      SUBROUTINE QMASS_MH(SQRTS,AS_S,MT_S,MB_S,MC_S,MS_S,MU_S,MD_S)
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
*---> running alpha_s and b-quark mass at SQRTS
*     : SQRTS > MS^pole assumed
*-----------------------------------------------------------------------
      PI      = 2.D0*DASIN(1.D0)
      B3      = (11.D0-2.D0/3.D0*3.D0)/4.D0/PI
      B4      = (11.D0-2.D0/3.D0*4.D0)/4.D0/PI
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
*
      AS_MT   = ASMT_H
      AS_MZ   = ASMZ_H
      AS_MB   = RAUX_H(3)
      AS_MC   = RAUX_H(6)
      MB_POLE = RAUX_H(1)
      MC_POLE = RAUX_H(4)
*Quark masses at mb^pole and mc^pole
      MT_MB = MTMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MB_MB = MBMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MC_MB = MCMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MS_MB = MSMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MU_MB = MUMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MD_MB = MDMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
*      print*,'at mb^pole:',mt_mb,mb_mb,mc_mb,ms_mb,mu_mb,md_mb
*      print*,'at mb^pole:',mb_mb,' ?= ',raux_h(2)
      MT_MC = MT_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MB_MC = MB_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MC_MC = MC_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MS_MC = MS_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MU_MC = MU_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MD_MC = MD_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
*      print*,'at mc^pole:',mt_mc,mb_mc,mc_mc,ms_mc,mu_mc,md_mc
*      print*,'at mc^pole:',mc_mc,' ?= ',raux_h(5)
*-----
*AS(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MTPOLE_H) THEN
       AS_S = AS_MT/(1.D0+B6*AS_MT*DLOG(SQRTS**2/MTPOLE_H**2))
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MTPOLE_H .AND. SQRTS.GT.MB_POLE ) THEN
       AS_S = AS_MZ/(1.D0+B5*AS_MZ*DLOG(SQRTS**2/MZ_H**2))
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
*MQ(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MTPOLE_H) THEN
       MT_S = MTMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MB_S = MBMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MC_S = MCMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MS_S = MSMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MU_S = MUMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MD_S = MDMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MTPOLE_H .AND. SQRTS.GT.MB_POLE ) THEN
       MT_S = MTMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MB_S = MBMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MC_S = MCMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MS_S = MSMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MU_S = MUMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MD_S = MDMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
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
*      print*,' > MQ(SQRTS)    :',MT_S,MB_S,MC_S,MS_S
*-----------------------------------------------------------------------
*      print*,'>>> QMASS_MH :',sqrts,as_s,mt_s,mb_s,mc_s,ms_s,mu_s,md_s
*=======================================================================
      RETURN
      END
