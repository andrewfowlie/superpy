      SUBROUTINE FILLCPsuperH2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
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
*Local assays:
*
      REAL*8     HZP_WIDTH(3),HZP_WIDTH_SM(3)
*=======================================================================
      CALL FILLPARA2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
      MCH=SSPARA_H(2)
*
      IF(IFLAG_H(57).GT.0) GOTO 99
*=======================================================================
      CALL FILLHIGGS2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     .               ,HMASS_H,OMIX_H)
*      print*,'>0 ',(hmass_h(i),i=1,3)
*
      CALL INCL_STAU(NFLAG,IFLAG_H,HMASS_H,OMIX_H)
*      print*,'>1 ',(hmass_h(i),i=1,3)
*
      MCH=RAUX_H(10) ! Charged Higgs pole mass or effetive-pot. mass
*
      IERR1=IFLAG_H(50)+IFLAG_H(51)+IFLAG_H(52)+IFLAG_H(54)
     .     +IFLAG_H(56)
     .     +IFLAG_H(55)+IFLAG_H(60)
      IF(IERR1.GT.0) GOTO 99
*
      IF(IFLAG_H(53).EQ.1) THEN
       print*,'WARNING! IFLAG_H(53) = ',IFLAG_H(53)
       IFLAG_H(53)=0
      ENDIF
*=======================================================================
      CALL FILLCOUPL2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     . ,HMASS_H,OMIX_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H)
*
      IERR2=IFLAG_H(54)+IFLAG_H(56)
      IF(IERR2.GT.0) GOTO 99
*=======================================================================
      CALL FILLGAMBR2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,NCMAX,NHC_H,SHC_H,CHC_H,STMASS_H
     . ,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     . ,MC_H,MN_H,NMNH,GAMBRN,NMCH,GAMBRC)
*
      CALL H2ZGAMMA_V0(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC,HZP_WIDTH,HZP_WIDTH_SM)
*
      CALL UPDATEGAMBR_0(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC,HZP_WIDTH,HZP_WIDTH_SM)
*
      MHSM=SMPARA_H(20)
      CALL FILLSMBRS(MHSM,GAMBRN,NMNH)
*=======================================================================
*For the subroutine FILLDHPG, we need \sqrt{s} value as an input. 
      SQRTS=RAUX_H(101)
*
      CALL FILLDHPG(SQRTS,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
      CALL FILLBOBS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .             ,HMASS_H,OMIX_H,MCH,MC_H,UL_H,UR_H,STMASS_H,STMIX_H)
*=======================================================================
      IF(ISKIP_EDM.EQ.0) THEN
      CALL FILLEDMS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)

      CALL HGRA_EDMS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
      ENDIF ! ISKIP_EDM
*=======================================================================
      IF(ISKIP_EDM.EQ.0) THEN
      CALL FILLMUON(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
      ENDIF ! ISKIP_EDM
*=======================================================================
      CALL FILLCOLL(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
      IF(IFLAG_H(30).EQ.1) THEN
      CALL FILLSLHA2(
     . NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
      ENDIF
*=======================================================================
 99   CONTINUE
*
      RETURN
      END
