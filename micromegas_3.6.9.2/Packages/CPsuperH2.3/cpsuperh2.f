      PROGRAM CPsuperH2
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
      DATA NAUX_H/999/
*-----------------------------------------------------------------------
*ARRAYS:
      REAL*8 SMPARA_H(20),SSPARA_H(38)
      DATA NSMIN/20/
      DATA NSSIN/38/
*
      INTEGER*8 IFLAG_H(100)
      DATA NFLAG/100/
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(100,3)  ! 100 = NCMAX
      REAL*8     SHC_H(100)
      COMPLEX*16 CHC_H(100)
      DATA NCMAX/100/
*
      REAL*8 GAMBRN(101,3,3)   ! 101 = IFLAG_H(20)+IFLAG_H(21)+1 = NMNH
*                                      ISMN       =ISUSYN        = 50
      REAL*8 GAMBRC(51,3)      !  51 = IFLAG_H(22)+IFLAG_H(23)+1 = NMCH
*                                      ISMC       =ISUSYC        = 25
      DATA NMNH/101/
      DATA NMCH/51/
*=======================================================================
      CALL FILLINIT2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,NCMAX,NHC_H,SHC_H,CHC_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
*
* To use other values for the input parameters than those in the "run"
* file, one can specifiy them here. To scan the phase of, for example,
* the gluino mass parameter M_3, one can do
*
*      DO IVAR=0,72
*       SSPARA_H(10)=5.D0*DBLE(IVAR) ! Phi_3
*       print*,'Phi_3 = ',SSPARA_H(10)
*
* Don't forget commenting in "ENDDO ! IVAR" at the end of this block.
*
*-----------------------------------------------------------------------
* For the \sqrt{s}-dependent propagators and the Higgs couplings to the
* gluons and photons, the following should be specified. If not, the
* value in FILLINIT2 is to be used:
*      RAUX_H(101)= ... !  \sqrt{s} for the subroutine FILLDHPG
*-----------------------------------------------------------------------
* One may skip the time-consuming EDM calculations by commenting in the
* follwing line:
*      ISKIP_EDM=1
*-----------------------------------------------------------------------
*
      CALL FILLCPsuperH2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
*
*Error messages:
*--a stop or sbottom squared mass is negative
       IF(IFLAG_H(50).EQ.1) THEN
        print*,'ERROR! IFLAG_H(50) = ',IFLAG_H(50)
        IFLAG_H(50)=0
        GOTO 99
       ENDIF
*
*--the Higgs--boson mass matrix contains a complex or negative eigenvalue
       IF(IFLAG_H(51).EQ.1) THEN
        print*,'ERROR! IFLAG_H(51) = ',IFLAG_H(51)
        IFLAG_H(51)=0
        GOTO 99
       ENDIF
*
*--the diagonalization of the Higgs mass matrix is not successful
       IF(IFLAG_H(52).EQ.1) THEN
        print*,'ERROR! IFLAG_H(52) = ',IFLAG_H(52)
        IFLAG_H(52)=0
        GOTO 99
       ENDIF
*
*--the iteration resumming the threshold corrections is not convergent
       IF(IFLAG_H(54).EQ.1) THEN
        print*,'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
        IFLAG_H(54)=0
        GOTO 99
       ENDIF
*
*--Yukawa coupling has a non--perturbative value: |h_{t,b}| > 2
       IF(IFLAG_H(55).EQ.1) THEN
        print*,'ERROR! IFLAG_H(55) = ',IFLAG_H(55)
        IFLAG_H(55)=0
        GOTO 99
       ENDIF
*
*-- 1 = a tau sneutrino or a stau squared mass is negative
*-- 2 = tachyonic stop or sbottom
*-- 3 = tachyonic scalar strange
       IF(IFLAG_H(56).GT.0) THEN
        IF(IFLAG_H(56).EQ.1) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IF(IFLAG_H(56).EQ.2) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IF(IFLAG_H(56).EQ.3) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
        IFLAG_H(56)=0
        GOTO 99
       ENDIF
*
*--one of the magnitudes of the complex input parameters is negative
       IF(IFLAG_H(57).EQ.1) THEN
        print*,'ERROR! IFLAG_H(57) = ',IFLAG_H(57)
        IFLAG_H(57)=0
        GOTO 99
       ENDIF
*
*--the iterative method for the neutral Higgs-boson pole masses fails
       IF(IFLAG_H(60).EQ.1) THEN
        print*,'ERROR! IFLAG_H(60) = ',IFLAG_H(60)
        IFLAG_H(60)=0
        GOTO 99
       ENDIF
*
*-----------------------------------------------------------------------
* Users may use the following subroutine for further analysis:
*
*      CALL AURUN(ISKIP_EDM
*     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
*     .,MCH,HMASS_H,OMIX_H
*     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
*     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
*     .,NMNH,GAMBRN,NMCH,GAMBRC)
*-----------------------------------------------------------------------
 99   CONTINUE
*
*      ENDDO ! IVAR
*=======================================================================
*
      STOP
      END
