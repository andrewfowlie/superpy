      SUBROUTINE AURUN(ISKIP_EDM
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
*Local Vairables:
*
*=======================================================================
      print*,'    >>>>> AURUN: START <<<<<'
*=======================================================================
*
*=======================================================================
      print*,'    >>>>> AURUN: E N D <<<<<'
*=======================================================================
 99   CONTINUE
*
      RETURN
      END
