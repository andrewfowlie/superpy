      SUBROUTINE FILLINIT2(ISKIP_EDM
     .,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,NCMAX,NHC_H,SHC_H,CHC_H,NMNH,GAMBRN,NMCH,GAMBRC)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------

      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
*
      INTEGER*8 IFLAG_H(NFLAG)
*
      COMPLEX*16 NHC_H(NCMAX,3) 
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)  ! 101 = IFLAG_H(20)+IFLAG_H(21)+1 = NMNH
*                                      ISMN       =ISUSYN        = 50
      REAL*8 GAMBRC(NMCH,3)    !  51 = IFLAG_H(22)+IFLAG_H(23)+1 = NMCH
*                                      ISMC       =ISUSYC        = 25
*=======================================================================
*initialize
*
      DO IFLAG=1,NFLAG
       IFLAG_H(IFLAG)=0
      ENDDO
*
      DO INC=1,NCMAX
        DO IH=1,3
         NHC_H(INC,IH)=DCMPLX(0.D0,0.D0)
        ENDDO
         CHC_H(INC)=DCMPLX(0.D0,0.D0)
         SHC_H(INC)=0.D0
      ENDDO
*
      DO IM=1,NMNH
        DO IWB=1,3
          DO IH=1,3
            GAMBRN(IM,IWB,IH)=0.D0
          ENDDO
        ENDDO
      ENDDO
      DO IM=1,NMCH
        DO IWB=1,3
          GAMBRC(IM,IWB)=0.D0
        ENDDO
      ENDDO
*
      DO IAUX=1,NAUX
        RAUX_H(IAUX)=0.D0
        CAUX_H(IAUX)=DCMPLX(0.D0,0.D0)
      ENDDO
*-----------------------------------------------------------------------
*read in input data from the file 'run'
*
      DO IP=1,NSMIN
       READ(*,*) SMPARA_H(IP)
      ENDDO
      DO IP=1,NSSIN
       READ(*,*) SSPARA_H(IP)
      ENDDO
*
*Set some flag's
*
      READ(*,*) IFLAG_H( 1) ! '1' will print input parameters
      READ(*,*) IFLAG_H( 2) ! '1' will print Higgs sector
      READ(*,*) IFLAG_H( 3) ! '1' will print masses and mixings of
*                             stop and sbottom sectors
      READ(*,*) IFLAG_H( 4) ! '1' will print masses and mixings of chargino and
*                             neutralino sectors
      READ(*,*) IFLAG_H( 5) ! '1-6' will print Higgs boson couplings
      READ(*,*) IFLAG_H( 6) ! '1-5' will print Higgs boson decays
      READ(*,*) IFLAG_H(10) ! if 0, include radiative corrections to
*                             H-top-top and H-bot-bot Yukawa couplings
      READ(*,*) IFLAG_H(11) ! if 0, use pole Higgs masses
*                             if 1, use effective potential Higgs mass
      READ(*,*) IFLAG_H(12) ! 5 For full improvemnt
       IF(IFLAG_H(12).EQ.0) IFLAG_H(12)=5
      READ(*,*) IFLAG_H(13) ! 1 Not to include the off-diagonal absorptive parts
      READ(*,*) IFLAG_H(14) ! 1 to print FILLDHPG results
      READ(*,*) IFLAG_H(16) ! 1 to print FILLBOBS results
      READ(*,*) IFLAG_H(17) ! 1 to print BTSGAM   results
      READ(*,*) IFLAG_H(18) ! 1 or 2 or 3 to print FILLEDMS results
      READ(*,*) IFLAG_H(19) ! 1 or 2 to print FILLMUON results
      READ(*,*) IFLAG_H(30) ! 1 to print the SLHA2 output
*
      IFLAG_H(20) = (NMNH-1)/2      ! ISMN
      IFLAG_H(21) = (NMNH-1)/2      ! ISUSYN = ISMN
      IFLAG_H(22) = (NMCH-1)/2      ! ISMC
      IFLAG_H(23) = (NMCH-1)/2      ! ISUSYC = ISMC
*
*\sqrt{s} for the subroutine FILLDHPG. It is HMASS_H(2) in "run.dist"
*
      RAUX_H(101)=0.27182383938642528D+03
*
*  
      ISKIP_EDM=0 ! to include EDM calculations by default
*=======================================================================
*
      RETURN
      END
