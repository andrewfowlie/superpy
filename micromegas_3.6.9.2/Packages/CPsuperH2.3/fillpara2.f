      SUBROUTINE FILLPARA2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
************************************************************************
*
* This subroutine fills the three common blocks /HC_SMPARA/, /HC_RSUSYPARA/,
* and /HC_CSUSYPARA/ from the input array SMPARA_H(NSMIN),SSPARA_H(NSSIN).
*
* RAUX_H(1) = m_b^pole
* RAUX_H(2) = m_b(m_b^pole)
* RAUX_H(3) = a_s(m_b^pole)
* RAUX_H(4) = m_c^pole
* RAUX_H(5) = m_c(m_c^pole)
* RAUX_H(6) = a_s(m_c^pole)
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
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
      COMPLEX*16 CKM(3,3)
*-----------------------------------------------------------------------
      PI     = 2.D0*DASIN(1.D0)
*
      AEM     = 1.D0/SMPARA_H( 1) 
      ASMZ    = SMPARA_H( 2)        
      MZ      = SMPARA_H( 3) 
      SWSQ    = SMPARA_H( 4) 
      ME      = SMPARA_H( 5)
      MMU     = SMPARA_H( 6) 
      MTAU    = SMPARA_H( 7) 
      MDMT    = SMPARA_H( 8) 
      MSMT    = SMPARA_H( 9) 
      MBMT    = SMPARA_H(10) 
      MUMT    = SMPARA_H(11) 
      MCMT    = SMPARA_H(12) 
      MTPOLE  = SMPARA_H(13) 
      GAMW    = SMPARA_H(14) 
      GAMZ    = SMPARA_H(15) 
      XL_CKM  = SMPARA_H(16)
      A_CKM   = SMPARA_H(17)
      RB_CKM  = SMPARA_H(18)
      EB_CKM  = SMPARA_H(19)
      MHSM    = SMPARA_H(20)
      TB      = SSPARA_H( 1) 
      MCHPOLE = SSPARA_H( 2) 
      MUMAG   = SSPARA_H( 3) 
      AMUD    = SSPARA_H( 4) 
      M1MAG   = SSPARA_H( 5) 
      AM1D    = SSPARA_H( 6) 
      M2MAG   = SSPARA_H( 7) 
      AM2D    = SSPARA_H( 8) 
      M3MAG   = SSPARA_H( 9) 
      AM3D    = SSPARA_H(10) 
      MQ3     = SSPARA_H(11) 
      MU3     = SSPARA_H(12) 
      MD3     = SSPARA_H(13) 
      ML3     = SSPARA_H(14) 
      ME3     = SSPARA_H(15) 
      ATMAG   = SSPARA_H(16) 
      AATD    = SSPARA_H(17) 
      ABMAG   = SSPARA_H(18) 
      AABD    = SSPARA_H(19) 
      ATUMAG  = SSPARA_H(20) 
      AATUD   = SSPARA_H(21) 
*JSL:05/APR/2005 Set IFLAG_H(57)=1 when the one of the mangitudes
*                of complex parameter is negative
      IF(MUMAG.LT.0.D0 .OR. M1MAG.LT.0.D0 .OR. M2MAG.LT.0.D0 .OR.
     .   M3MAG.LT.0.D0 .OR. ATMAG.LT.0.D0 .OR. ABMAG.LT.0.D0 .OR.
     .   ATUMAG.LT.0.D0) THEN
       IFLAG_H(57)=1
       RETURN
      ENDIF
*
************************************************************************
* /HC_SMPARA/
************************************************************************
* ---> Directly from SMPARA_H Array
*
      AEM_H    = AEM
      ASMZ_H   = ASMZ
      MZ_H     = MZ
      SW_H     = DSQRT(SWSQ)
      ME_H     = ME
      MMU_H    = MMU
      MTAU_H   = MTAU
      MDMT_H   = MDMT
      MSMT_H   = MSMT
      MBMT_H   = MBMT
      MUMT_H   = MUMT
      MCMT_H   = MCMT
      MTPOLE_H = MTPOLE
      GAMW_H   = GAMW
      GAMZ_H   = GAMZ
*
* ---> Induced SM Parameters
*
      NF     = 5
      B5     = (11.D0-2.D0/3.D0*DBLE(NF))/4.D0/PI
*
      EEM_H  = DSQRT(4.D0*PI*AEM_H)
      ASMT_H = ASMZ_H/(1.D0+B5*ASMZ_H*DLOG(MTPOLE_H**2/MZ_H**2)) ! AS(MT^POLE)
      CW_H   = DSQRT(1.D0-SW_H**2)
      TW_H   = SW_H/CW_H
      MW_H   = MZ_H*DSQRT(1.D0-SW_H**2)
      GW_H   = EEM_H/SW_H
      GP_H   = EEM_H/CW_H
      V_H    = 2.D0*MW_H/GW_H
*-----
*2009/Mar/06 JSL: GF as an in put
*      GF_H   = DSQRT(2.D0)*GW_H**2/8.D0/MW_H**2
      GF_H   = 1.16637D-5 
      MW_H   = DSQRT(GW_H**2/4.D0/DSQRT(2.D0)/GF_H) ! Recalculate MW
      V_H    = 2.D0*MW_H/GW_H                       ! Recaluclate V_H
*-----
      MTMT_H = MTPOLE_H/(1.D0+4.D0*ASMT_H/3.D0/PI) ! MT(MT^POLE)
*
      CALL GET_CKM(XL_CKM,A_CKM,RB_CKM,EB_CKM,CKM)
************************************************************************
* /HC_RSUSYPARA/
************************************************************************
* 
      TB_H  = TB
      CB_H  = 1.D0/DSQRT(1.D0+TB**2)
      SB_H  = TB/DSQRT(1.D0+TB**2)
      MQ3_H = MQ3
      MU3_H = MU3
      MD3_H = MD3
      ML3_H = ML3
      ME3_H = ME3
*
************************************************************************
* /HC_CSUSYPARA/
************************************************************************
* 
      MU_H   = MUMAG*DCMPLX(DCOS(AMUD/180.D0*PI),DSIN(AMUD/180.D0*PI))
      M1_H   = M1MAG*DCMPLX(DCOS(AM1D/180.D0*PI),DSIN(AM1D/180.D0*PI))
      M2_H   = M2MAG*DCMPLX(DCOS(AM2D/180.D0*PI),DSIN(AM2D/180.D0*PI))
      M3_H   = M3MAG*DCMPLX(DCOS(AM3D/180.D0*PI),DSIN(AM3D/180.D0*PI))
      AT_H   = ATMAG*DCMPLX(DCOS(AATD/180.D0*PI),DSIN(AATD/180.D0*PI))
      AB_H   = ABMAG*DCMPLX(DCOS(AABD/180.D0*PI),DSIN(AABD/180.D0*PI))
      ATAU_H = ATUMAG
     .              *DCMPLX(DCOS(AATUD/180.D0*PI),DSIN(AATUD/180.D0*PI))
*
************************************************************************
* /HC_RAUX/
************************************************************************
      B3      = (11.D0-2.D0/3.D0*3.D0)/4.D0/PI
      B4      = (11.D0-2.D0/3.D0*4.D0)/4.D0/PI
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      AS_MZ   = ASMZ_H
      AS_MT   = ASMT_H
*      print*,b3,b4,b5,b6,as_mz,as_mt

      NQMAX   = 100 ! the maximum number of iteration for the quark-pole masses
*b-quark pole mass
      MB_POLE=MBMT_H
      DO I=1,NQMAX
        AS_MB     = AS_MZ/(1.D0+B5*AS_MZ*DLOG(MB_POLE**2/MZ_H**2))
        MB_MB     = MBMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
        MB_POLE_N = MB_MB*(1.D0+4.D0*AS_MB/3.D0/PI)
*        print*,'b-quark pole mass',I,AS_MB,MB_POLE,MB_POLE_N
        IF(DABS((MB_POLE_N-MB_POLE)/MB_POLE).LT.1.D-6) THEN
         MB_POLE  = MB_POLE_N
         GOTO 199
        ENDIF
        MB_POLE   = MB_POLE_N
      ENDDO
 199  CONTINUE
*      print*,'mb^pole, mb(mb^pole), As(mb^pole) :',mb_pole,mb_mb,as_mb
      RAUX_H(1)=MB_POLE
      RAUX_H(2)=MB_MB
      RAUX_H(3)=AS_MB
*c-quark pole mass
      MC_POLE=MCMT_H
      DO I=1,NQMAX
        AS_MC     = AS_MB/(1.D0+B4*AS_MB*DLOG(MC_POLE**2/MB_POLE**2))
        MC_MB     = MCMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
        MC_MC     = MC_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
        MC_POLE_N = MC_MC*(1.D0+4.D0*AS_MC/3.D0/PI)
*        print*,'c-quark pole mass',I,AS_MC,MC_POLE,MC_POLE_N
        IF(DABS((MC_POLE_N-MC_POLE)/MC_POLE).LT.1.D-6) THEN
         MC_POLE  = MC_POLE_N
         GOTO 198
        ENDIF
        MC_POLE   = MC_POLE_N
      ENDDO
 198  CONTINUE
*      print*,'mc^pole, mc(mc^pole), As(mc^pole) :',mc_pole,mc_mc,as_mc
      RAUX_H(4)=MC_POLE
      RAUX_H(5)=MC_MC
      RAUX_H(6)=AS_MC
************************************************************************
      IF(IFLAG_H(1).EQ.1) CALL DUMP_INPUT(MCHPOLE,CKM
     .   ,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
************************************************************************
      RETURN
      END

      SUBROUTINE DUMP_INPUT(MCHPOLE,CKM
     .   ,NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
************************************************************************
*
* This subroutine dumps the contents of the three common blocks /HC_SMPARA/,
* /HC_RSUSYPARA/, and /HC_CSUSYPARA/ filled by the subroutine FILLPARA.
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      REAL*8 SMPARA_H(NSMIN),SSPARA_H(NSSIN)
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
      COMPLEX*16 CKM(3,3)
*-----------------------------------------------------------------------
*
* /HC_SMPARA/
*
      print*,'---------------------------------------------------------'
      print*,'Standard Model Parameters  in /HC_SMPARA/'
      print*,'---------------------------------------------------------'
      WRITE(*, 1) AEM_H
      WRITE(*, 2) ASMZ_H
      WRITE(*, 3) MZ_H
      WRITE(*, 4) SW_H
      WRITE(*, 5) ME_H
      WRITE(*, 6) MMU_H
      WRITE(*, 7) MTAU_H
      WRITE(*, 8) MDMT_H
      WRITE(*, 9) MSMT_H
      WRITE(*,10) MBMT_H
      WRITE(*,11) MUMT_H
      WRITE(*,12) MCMT_H
      WRITE(*,13) MTPOLE_H
      WRITE(*,314) SMPARA_H(20)
      WRITE(*,94) GAMW_H
      WRITE(*,95) GAMZ_H
      WRITE(*,51) EEM_H
      WRITE(*,52) ASMT_H
      WRITE(*,53) CW_H
      WRITE(*,54) TW_H
      WRITE(*,55) MW_H
      WRITE(*,56) GW_H
      WRITE(*,57) GP_H
      WRITE(*,58) V_H
      WRITE(*,59) GF_H
      WRITE(*,60) MTMT_H
      print*,'---------------------------------------------------------'
      print*,'CKM Matrix : '
*      WRITE(*,111) (DREAL(CKM(1,J)),DIMAG(CKM(1,J)),J=1,3)
*      WRITE(*,112) (DREAL(CKM(2,J)),DIMAG(CKM(2,J)),J=1,3)
*      WRITE(*,113) (DREAL(CKM(3,J)),DIMAG(CKM(3,J)),J=1,3)
      WRITE(*,211) DREAL(CKM(1,1)),DIMAG(CKM(1,1)),CDABS(CKM(1,1))
      WRITE(*,212) DREAL(CKM(1,2)),DIMAG(CKM(1,2)),CDABS(CKM(1,2))
      WRITE(*,213) DREAL(CKM(1,3)),DIMAG(CKM(1,3)),CDABS(CKM(1,3))
      WRITE(*,214) DREAL(CKM(2,1)),DIMAG(CKM(2,1)),CDABS(CKM(2,1))
      WRITE(*,215) DREAL(CKM(2,2)),DIMAG(CKM(2,2)),CDABS(CKM(2,2))
      WRITE(*,216) DREAL(CKM(2,3)),DIMAG(CKM(2,3)),CDABS(CKM(2,3))
      WRITE(*,217) DREAL(CKM(3,1)),DIMAG(CKM(3,1)),CDABS(CKM(3,1))
      WRITE(*,218) DREAL(CKM(3,2)),DIMAG(CKM(3,2)),CDABS(CKM(3,2))
      WRITE(*,219) DREAL(CKM(3,3)),DIMAG(CKM(3,3)),CDABS(CKM(3,3))
      print*,'---------------------------------------------------------'
*
* /HC_RSUSYPARA/
* 
      print*,'Real SUSY Parameters  in /HC_RSUSYPARA/'
      print*,'---------------------------------------------------------'
      WRITE(*,14) TB_H
      WRITE(*,61) CB_H
      WRITE(*,62) SB_H
      WRITE(*,24) MQ3_H
      WRITE(*,25) MU3_H
      WRITE(*,26) MD3_H
      WRITE(*,96) ML3_H
      WRITE(*,97) ME3_H
      print*,'---------------------------------------------------------'
* 
* /HC_CSUSYPARA/
* 
      print*,'Complex SUSY Parameters  in /HC_CSUSYPARA/'
      print*,'---------------------------------------------------------'
      WRITE(*,16) CDABS(MU_H)
      WRITE(*,18) CDABS(M1_H)
      WRITE(*,20) CDABS(M2_H)
      WRITE(*,22) CDABS(M3_H)
      WRITE(*,27) CDABS(AT_H)
      WRITE(*,29) CDABS(AB_H)
      WRITE(*,98) CDABS(ATAU_H)
      WRITE(*,17) ADEG(MU_H)
      WRITE(*,19) ADEG(M1_H)
      WRITE(*,21) ADEG(M2_H)
      WRITE(*,23) ADEG(M3_H)
      WRITE(*,28) ADEG(AT_H)
      WRITE(*,30) ADEG(AB_H)
      WRITE(*,99) ADEG(ATAU_H)
      print*,'---------------------------------------------------------'
      print*,'Diagonal Sfermion Mass Matrices [GeV] (Not squared) :'
      WRITE(*,114) SSPARA_H(11),SSPARA_H(22),SSPARA_H(22),1.D0
      WRITE(*,115) SSPARA_H(12),SSPARA_H(23),SSPARA_H(23),1.D0
      WRITE(*,116) SSPARA_H(13),SSPARA_H(24),SSPARA_H(24),1.D0
      WRITE(*,117) SSPARA_H(14),SSPARA_H(25),SSPARA_H(25),1.D0
      WRITE(*,118) SSPARA_H(15),SSPARA_H(26),SSPARA_H(26),1.D0
      print*,'---------------------------------------------------------'
      WRITE(*,15) MCHPOLE
      print*,'---------------------------------------------------------'
      print*,' '
*
************************************************************************
*
* FORMAT's for /HC_SMPARA/
*
************************************************************************
  1   FORMAT(' AEM_H    = ',E10.4,' : alpha_em(MZ)                    ')
  2   FORMAT(' ASMZ_H   = ',E10.4,' : alpha_s(MZ)                     ')
  3   FORMAT(' MZ_H     = ',E10.4,' : Z boson mass in GeV             ')
  4   FORMAT(' SW_H     = ',E10.4,' : sinTheta_W                      ')
  5   FORMAT(' ME_H     = ',E10.4,' : electron mass in GeV            ')
  6   FORMAT(' MMU_H    = ',E10.4,' : muon mass in GeV                ')
  7   FORMAT(' MTAU_H   = ',E10.4,' : tau mass in GeV                 ')
  8   FORMAT(' MDMT_H   = ',E10.4,' : d-quark mass at M_t^pole in GeV ')
  9   FORMAT(' MSMT_H   = ',E10.4,' : s-quark mass at M_t^pole in GeV ')
 10   FORMAT(' MBMT_H   = ',E10.4,' : b-quark mass at M_t^pole in GeV ')
 11   FORMAT(' MUMT_H   = ',E10.4,' : u-quark mass at M_t^pole in GeV ')
 12   FORMAT(' MCMT_H   = ',E10.4,' : c-quark mass at M_t^pole in GeV ')
 13   FORMAT(' MTPOLE_H = ',E10.4,' : t-quark pole mass in GeV        ')
 314  FORMAT(' MHSM     = ',E10.4,' : mass of the SM Higgs in GeV     ')
 94   FORMAT(' GAMW_H   = ',E10.4,' : Gam_W in GeV                    ')
 95   FORMAT(' GAMZ_H   = ',E10.4,' : Gam_Z in GeV                    ')
 111  format(1x,'/',3('(',e10.4,1x,e10.4,')',1x),'\\')
 112  format(1x,'|',3('(',e10.4,1x,e10.4,')',1x),'|')
 113  format(1x,'\\',3('(',e10.4,1x,e10.4,')',1x),'/')
 211  format(' |V_ud|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 212  format(' |V_us|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 213  format(' |V_ub|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 214  format(' |V_cd|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 215  format(' |V_cs|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 216  format(' |V_cb|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 217  format(' |V_td|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 218  format(' |V_ts|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
 219  format(' |V_tb|   = |(',e10.4,1x,e10.4,')| = ',e10.4)
************************************************************************
 51   FORMAT(' EEM_H    = ',E10.4,' : e = (4*pi*alpha_em)^1/2         ')
 52   FORMAT(' ASMT_H   = ',E10.4,' : alpha_s(M_t^pole)               ')
 53   FORMAT(' CW_H     = ',E10.4,' : cosTheta_W                      ')
 54   FORMAT(' TW_H     = ',E10.4,' : tanTheta_W                      ')
*2009/Mar/06 JSL
* 55   FORMAT(' MW_H     = ',E10.4,' : W boson mass MW = MZ*CW         ')
 55   FORMAT(' MW_H     = ',E10.4,' : MW ={gw^2 / [4 sqrt(2) G_F]}^1/2')
 56   FORMAT(' GW_H     = ',E10.4,' : SU(2) gauge coupling  gw=e/s_W  ')
 57   FORMAT(' GP_H     = ',E10.4,' : U(1)_Y gauge coupling gp=e/c_W  ')
 58   FORMAT(' V_H      = ',E10.4,' : V = 2 MW / gw                   ')
*2009/Mar/06 JSL
* 59   FORMAT(' GF_H     = ',E10.4,' : GF=sqrt(2)*gw^2/8 MW^2 in GeV^-2')
 59   FORMAT(' GF_H     = ',E10.4,' : GF                              ')
 60   FORMAT(' MTMT_H   = ',E10.4,' : t-quark mass at M_t^pole in GeV ')
************************************************************************
*
* FORMAT's for /HC_RSUSYPARA/
*
************************************************************************
*
 14   FORMAT(' TB_H     = ',E10.4,' : tan(beta)                       ')
 15   FORMAT(' Charged Higgs boson pole mass : ',E10.4,' GeV          ')
 61   FORMAT(' CB_H     = ',E10.4,' : cos(beta)                       ')
 62   FORMAT(' SB_H     = ',E10.4,' : sin(beta)                       ')
 24   FORMAT(' MQ3_H    = ',E10.4,' : M_tilde{Q_3} in GeV             ')
 25   FORMAT(' MU3_H    = ',E10.4,' : M_tilde{U_3} in GeV             ')
 26   FORMAT(' MD3_H    = ',E10.4,' : M_tilde{D_3} in GeV             ')
 96   FORMAT(' ML3_H    = ',E10.4,' : M_tilde{L_3} in GeV             ')
 97   FORMAT(' ME3_H    = ',E10.4,' : M_tilde{E_3} in GeV             ')
************************************************************************
*
* FORMAT's for /HC_CSUSYPARA/
*
************************************************************************
 16   FORMAT(' |MU_H|     = ',E10.4,':Mag. of MU parameter in GeV     ')
 18   FORMAT(' |M1_H|     = ',E10.4,':Mag. of M1 parameter in GeV     ')
 20   FORMAT(' |M2_H|     = ',E10.4,':Mag. of M2 parameter in GeV     ')
 22   FORMAT(' |M3_H|     = ',E10.4,':Mag. of M3 parameter in GeV     ')
 27   FORMAT(' |AT_H|     = ',E10.4,':Mag. of AT parameter in GeV     ')
 29   FORMAT(' |AB_H|     = ',E10.4,':Mag. of AB parameter in GeV     ')
 98   FORMAT(' |ATAU_H|   = ',E10.4,':Mag. of ATAU parameter in GeV   ')
 17   FORMAT(' ARG(MU_H)  = ',E10.4,':Arg. of MU parameter in Degree  ')
 19   FORMAT(' ARG(M1_H)  = ',E10.4,':Arg. of M1 parameter in Degree  ')
 21   FORMAT(' ARG(M2_H)  = ',E10.4,':Arg. of M2 parameter in Degree  ')
 23   FORMAT(' ARG(M3_H)  = ',E10.4,':Arg. of M3 parameter in Degree  ')
 28   FORMAT(' ARG(AT_H)  = ',E10.4,':Arg. of AT parameter in Degree  ')
 30   FORMAT(' ARG(AB_H)  = ',E10.4,':Arg. of AB parameter in Degree  ')
 99   FORMAT(' ARG(ATAU_H)= ',E10.4,':Arg. of ATAU parameter in Degree')
************************************************************************
 114  FORMAT(1x,'M_Q = ',e10.4,' x Diag(',e10.4,1x,e10.4,1x,e10.4,')')
 115  FORMAT(1x,'M_U = ',e10.4,' x Diag(',e10.4,1x,e10.4,1x,e10.4,')')
 116  FORMAT(1x,'M_D = ',e10.4,' x Diag(',e10.4,1x,e10.4,1x,e10.4,')')
 117  FORMAT(1x,'M_L = ',e10.4,' x Diag(',e10.4,1x,e10.4,1x,e10.4,')')
 118  FORMAT(1x,'M_E = ',e10.4,' x Diag(',e10.4,1x,e10.4,1x,e10.4,')')
************************************************************************
*
      RETURN
      END

      REAL*8 FUNCTION ADEG(CX)
************************************************************************
*
************************************************************************
*
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 CX
      PI     = 2.D0*DASIN(1.D0)
*
      ARG_L=DABS(DATAN(DIMAG(CX)/DREAL(CX)))
      DEG_L=DABS(DATAN(DIMAG(CX)/DREAL(CX)))*180.D0/PI
*
      IF(DIMAG(CX).GT.0.D0.AND.DREAL(CX).GT.0.D0) ADEG=DEG_L
      IF(DIMAG(CX).GT.0.D0.AND.DREAL(CX).LT.0.D0) ADEG=180.D0-DEG_L
      IF(DIMAG(CX).LT.0.D0.AND.DREAL(CX).LT.0.D0) ADEG=DEG_L+180.D0
      IF(DIMAG(CX).LT.0.D0.AND.DREAL(CX).GT.0.D0) ADEG=360.D0-DEG_L
      IF(DIMAG(CX).GT.0.D0.AND.DREAL(CX).EQ.0.D0) ADEG=90.D0
      IF(DIMAG(CX).LT.0.D0.AND.DREAL(CX).EQ.0.D0) ADEG=270.D0
      IF(DIMAG(CX).EQ.0.D0.AND.DREAL(CX).LT.0.D0) ADEG=180.D0
      IF(DIMAG(CX).EQ.0.D0.AND.DREAL(CX).GT.0.D0) ADEG=0.D0
*
      RETURN
      END

      SUBROUTINE GET_CKM(L,A,RB,EB,V)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 V(3,3),VDAG(3,3)
      COMPLEX*16 XI,S13C
      REAL*8     L,A,RB,EB
*
      XI=DCMPLX(0.D0,1.D0)
*
      S12 =L
      S23 =A*L**2
      S13C=A*L**3*(RB+XI*EB)*DSQRT(1.D0-A**2*L**4)/DSQRT(1.D0-L**2)
     .    /(1.D0-A**2*L**4*(RB+XI*EB))
      C12=DSQRT(1.D0-S12**2)
      C23=DSQRT(1.D0-S23**2)
      C13=DSQRT(1.D0-CDABS(S13C)**2)
*
       V(1,1)= C12*C13
       V(1,2)= S12*C13
       V(1,3)= DCONJG(S13C)
       V(2,1)=-S12*C23-C12*S23*S13C
       V(2,2)= C12*C23-S12*S23*S13C
       V(2,3)= S23*C13
       V(3,1)= S12*S23-C12*C23*S13C
       V(3,2)=-C12*S23-S12*C23*S13C
       V(3,3)= C23*C13
*
      RETURN
      END


