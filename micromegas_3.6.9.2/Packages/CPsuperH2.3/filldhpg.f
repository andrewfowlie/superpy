      SUBROUTINE FILLDHPG(SQRTS,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
************************************************************************
*
* To call subroutines HPROP, HGG_S, and HPP_S which depend on sqrt{s}.
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
*Local:
      INTEGER*8 IFLAG_H(NFLAG)
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)  
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)   
      REAL*8 GAMBRC(NMCH,3)
*=======================================================================
*RAUX_H(101) is reserved for sqrt{s}
*      RAUX_H(101)=SQRTS
*
*                      H1          H2          H3          G0
*              H1 / CAUX_H(100) CAUX_H(101) CAUX_H(102) CAUX_U(103) \
*  DNH4[4,4] = H2 | CAUX_H(104) CAUX_H(105) CAUX_H(106) CAUX_H(107) |
*              H3 | CAUX_H(108) CAUX_H(109) CAUX_H(110) CAUX_H(111) |
*              G0 \ CAUX_H(112) CAUX_H(113) CAUX_H(114) CAUX_H(115) /
*
      IKILL=IFLAG_H(13)
      CALL HPROP(SQRTS**2,IKILL,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
*
      IF(IFLAG_H(14).EQ.1) THEN
      print*,'---------------------------------------------------------'
      write(*,10) raux_h(101)
      print*,'---------------------------------------------------------'
      write(*,11) dreal(caux_h(100)),dimag(caux_h(100))
     .           ,cdabs(caux_h(100))
      write(*,12) dreal(caux_h(105)),dimag(caux_h(105))
     .           ,cdabs(caux_h(105))
      write(*,13) dreal(caux_h(110)),dimag(caux_h(110))
     .           ,cdabs(caux_h(110))
      write(*,14) dreal(caux_h(101)),dimag(caux_h(101))
     .           ,cdabs(caux_h(101))
      write(*,15) dreal(caux_h(102)),dimag(caux_h(102))
     .           ,cdabs(caux_h(102))
      write(*,16) dreal(caux_h(106)),dimag(caux_h(106))
     .           ,cdabs(caux_h(106))
      write(*,17) dreal(caux_h(112)),dimag(caux_h(112))
     .           ,cdabs(caux_h(112))
      write(*,18) dreal(caux_h(113)),dimag(caux_h(113))
     .           ,cdabs(caux_h(113))
      write(*,19) dreal(caux_h(114)),dimag(caux_h(114))
     .           ,cdabs(caux_h(114))
      write(*,20) dreal(caux_h(115)),dimag(caux_h(115))
     .           ,cdabs(caux_h(115))
      ENDIF
*=======================================================================
*
*                          H+          G+
*     DCH2[2,2] =  H+ / CAUX_H(116) CAUX_U(117) \
*                  G+ \ CAUX_H(118) CAUX_H(119) /
*
      IF(IFLAG_H(14).EQ.1) THEN
      print*,'---------------------------------------------------------'
      write(*,21) raux_h(101)
      print*,'---------------------------------------------------------'
      write(*,22) dreal(caux_h(116)),dimag(caux_h(116))
     .           ,cdabs(caux_h(116))
      write(*,23) dreal(caux_h(117)),dimag(caux_h(117))
     .           ,cdabs(caux_h(117))
      write(*,24) dreal(caux_h(118)),dimag(caux_h(118))
     .           ,cdabs(caux_h(118))
      write(*,25) dreal(caux_h(119)),dimag(caux_h(119))
     .           ,cdabs(caux_h(119))
      ENDIF
*=======================================================================
* Higgs-Photon-Photon couplings at sqrt(s)=SQRTS:
*
* CP-even couplings:       | CP-odd couplings:
* CAUX_H(130) = SPHO(IH=1) | CAUX_H(131) = PPHO(IH=1)
* CAUX_H(132) = SPHO(IH=2) | CAUX_H(133) = PPHO(IH=2)
* CAUX_H(134) = SPHO(IH=3) | CAUX_H(135) = PPHO(IH=3)
*
      CALL HPP_S(SQRTS,MC_H,MCH,STMASS_H,SBMASS_H,STAUMASS_H
     .          ,NCMAX,NHC_H)
*
      IF(IFLAG_H(14).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,' Comparisons of the H-photon-photon coupings at MH^pole' 
      write(*,30) raux_h(101)
      print*,'---------------------------------------------------------'
      print*,'                 S couplings             P couplings'
      write(*,31) dreal(nhc_h(88,1)),dimag(nhc_h(88,1))
     .           ,dreal(nhc_h(89,1)),dimag(nhc_h(89,1))
      write(*,32) dreal(caux_h(130)),dimag(caux_h(130))
     .           ,dreal(caux_h(131)),dimag(caux_h(131))
      write(*,33) dreal(nhc_h(88,2)),dimag(nhc_h(88,2))
     .           ,dreal(nhc_h(89,2)),dimag(nhc_h(89,2))
      write(*,34) dreal(caux_h(132)),dimag(caux_h(132))
     .           ,dreal(caux_h(133)),dimag(caux_h(133))
      write(*,35) dreal(nhc_h(88,3)),dimag(nhc_h(88,3))
     .           ,dreal(nhc_h(89,3)),dimag(nhc_h(89,3))
      write(*,36) dreal(caux_h(134)),dimag(caux_h(134))
     .           ,dreal(caux_h(135)),dimag(caux_h(135))
      ENDIF
*=======================================================================
* Higgs-Gluon-Gluon couplings at sqrt(s)=SQRTS:
*
* CP-even couplings:       | CP-odd couplings:
* CAUX_H(140) = SGLUE(IH=1) | CAUX_H(141) = PGLUE(IH=1)
* CAUX_H(142) = SGLUE(IH=2) | CAUX_H(143) = PGLUE(IH=2)
* CAUX_H(144) = SGLUE(IH=3) | CAUX_H(145) = PGLUE(IH=3)
*
      CALL HGG_S(SQRTS,STMASS_H,SBMASS_H,NCMAX,NHC_H)
*
      IF(IFLAG_H(14).EQ.1) THEN
      print*,'---------------------------------------------------------'
      print*,' Comparisons of the H-glue-glue coupings at MH^pole'
      write(*,30) raux_h(101)
      print*,'---------------------------------------------------------'
      print*,'                 S couplings             P couplings'
      write(*,41) dreal(nhc_h(84,1)),dimag(nhc_h(84,1))
     .           ,dreal(nhc_h(85,1)),dimag(nhc_h(85,1))
      write(*,42) dreal(caux_h(140)),dimag(caux_h(140))
     .           ,dreal(caux_h(141)),dimag(caux_h(141))
      write(*,43) dreal(nhc_h(84,2)),dimag(nhc_h(84,2))
     .           ,dreal(nhc_h(85,2)),dimag(nhc_h(85,2))
      write(*,44) dreal(caux_h(142)),dimag(caux_h(142))
     .           ,dreal(caux_h(143)),dimag(caux_h(143))
      write(*,45) dreal(nhc_h(84,3)),dimag(nhc_h(84,3))
     .           ,dreal(nhc_h(85,3)),dimag(nhc_h(85,3))
      write(*,46) dreal(caux_h(144)),dimag(caux_h(144))
     .           ,dreal(caux_h(145)),dimag(caux_h(145))
      print*,'---------------------------------------------------------'
      print*,' '
*
      ENDIF
*=======================================================================
*HAHA
  10  FORMAT(2X,'DNH4 at sqrt{s} = ',E10.4,' GeV')
  11  FORMAT(2X,'DNH4[H1,H1]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  12  FORMAT(2X,'DNH4[H2,H2]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  13  FORMAT(2X,'DNH4[H3,H3]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  14  FORMAT(2X,'DNH4[H1,H2]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  15  FORMAT(2X,'DNH4[H1,H3]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  16  FORMAT(2X,'DNH4[H2,H3]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  17  FORMAT(2X,'DNH4[G0,H1]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  18  FORMAT(2X,'DNH4[G0,H2]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  19  FORMAT(2X,'DNH4[G0,H3]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  20  FORMAT(2X,'DNH4[G0,G0]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  21  FORMAT(2X,'DCH2 at sqrt{s} = ',E10.4,' GeV')
  22  FORMAT(2X,'DCH2[H+,H+]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  23  FORMAT(2X,'DCH2[H+,G+]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  24  FORMAT(2X,'DCH2[G+,H+]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  25  FORMAT(2X,'DCH2[G+,G+]: |(',E10.4,1X,E10.4,')| = ',E10.4)
  30  FORMAT(2X,'and those at sqrt{s} = ',E10.4,' GeV')
  31  FORMAT(2X,'H1PP(M): ',2('(',E10.4,1X,E10.4,')',1X))
  32  FORMAT(2X,'H1PP(S): ',2('(',E10.4,1X,E10.4,')',1X))
  33  FORMAT(2X,'H2PP(M): ',2('(',E10.4,1X,E10.4,')',1X))
  34  FORMAT(2X,'H2PP(S): ',2('(',E10.4,1X,E10.4,')',1X))
  35  FORMAT(2X,'H3PP(M): ',2('(',E10.4,1X,E10.4,')',1X))
  36  FORMAT(2X,'H3PP(S): ',2('(',E10.4,1X,E10.4,')',1X))
  41  FORMAT(2X,'H1GG(M): ',2('(',E10.4,1X,E10.4,')',1X))
  42  FORMAT(2X,'H1GG(S): ',2('(',E10.4,1X,E10.4,')',1X))
  43  FORMAT(2X,'H2GG(M): ',2('(',E10.4,1X,E10.4,')',1X))
  44  FORMAT(2X,'H2GG(S): ',2('(',E10.4,1X,E10.4,')',1X))
  45  FORMAT(2X,'H3GG(M): ',2('(',E10.4,1X,E10.4,')',1X))
  46  FORMAT(2X,'H3GG(S): ',2('(',E10.4,1X,E10.4,')',1X))
*
      RETURN
      END


      SUBROUTINE HPROP(S,IKILL,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
************************************************************************
* To calculate the propagate matrix:
*
*     DNH4(I,J)=s[(s-M_I^2) delta_IJ + i*PI^hat_IJ]^-1 
*     DCH2(I,J)=s[(s-MCH^2) delta_IJ + i*PI^hat_IJ]^-1 
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
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3) 
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)   
      REAL*8 GAMBRC(NMCH,3)
*-----------------------------------------------------------------------
* Local Variables:

      COMPLEX*16 XI
*
      REAL*8     HATPI(3,3),G0_HATPI(4),CH_HATPI(2,2)
      COMPLEX*16 DNH3(3,3),DNH4(4,4),DCH2(2,2)
      COMPLEX*16 X11,X12,X13,X21,X22,X23,X31,X32,X33,DETX
      COMPLEX*16 Y11,Y12,Y13,Y14,Y21,Y22,Y23,Y24
      COMPLEX*16 Y31,Y32,Y33,Y34,Y41,Y42,Y43,Y44,DETY
      COMPLEX*16 Z11,Z12,Z21,Z22,DETZ
*
      XI=DCMPLX(0.D0,1.D0)
*-----------------------------------------------------------------------
*
      DO I=1,3
       DO J=1,3
        CALL HSELFIJ(S,I,J,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,HATPI(I,J))
*        print*,i,j,hatpi(I,J),hmass_h(i)*gambrn(101,1,i)
       ENDDO
      ENDDO
*
      DO I=1,4
        CALL NGSELFI(S,I,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,G0_HATPI(I))
*       print*,'>> FILLDHPG: G0_HATPI(I) ',i,g0_hatpi(I)
      ENDDO
*
      DO I=1,2
       DO J=1,2
        CALL CHSELFIJ(S,I,J,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMCH,GAMBRC,CH_HATPI(I,J))
*        print*,i,j,ch_hatpi(I,J)
       ENDDO
      ENDDO
*-----------------------------------------------------------------------
*--> Calculating DNH3(I,J)=s[(s-M_I^2) d_IJ + i*PI^hat_IJ]^-1 
*
      X11=1.D0-HMASS_H(1)**2/S+XI*HATPI(1,1)/S
      X12=XI*HATPI(1,2)/S
      X13=XI*HATPI(1,3)/S
*
      X21=XI*HATPI(2,1)/S
      X22=1.D0-HMASS_H(2)**2/S+XI*HATPI(2,2)/S
      X23=XI*HATPI(2,3)/S
*
      X31=XI*HATPI(3,1)/S
      X32=XI*HATPI(3,2)/S
      X33=1.D0-HMASS_H(3)**2/S+XI*HATPI(3,3)/S

*To kill off-diagonal elements
      IF (IKILL.EQ.1) THEN
       X12=DCMPLX(0.D0,0.D0)
       X21=DCMPLX(0.D0,0.D0)
       X13=DCMPLX(0.D0,0.D0)
       X31=DCMPLX(0.D0,0.D0)
       X23=DCMPLX(0.D0,0.D0)
       X32=DCMPLX(0.D0,0.D0)
      ENDIF
*
      DETX=X11*X22*X33-X11*X23*X32-X21*X12*X33+X21*X13*X32
     .    +X31*X12*X23-X31*X13*X22
      IF(CDABS(DETX).EQ.0.D0) THEN
       PRINT*,'ERROR!: (s-m^2)^-1'
       STOP
      ENDIF
*
      DNH3(1,1)= (X22*X33-X23*X32)/DETX
      DNH3(1,2)=-(X12*X33-X13*X32)/DETX
      DNH3(1,3)= (X12*X23-X13*X22)/DETX
      DNH3(2,1)=-(X21*X33-X23*X31)/DETX
      DNH3(2,2)= (X11*X33-X13*X31)/DETX
      DNH3(2,3)=-(X11*X23-X13*X21)/DETX
      DNH3(3,1)= (X21*X32-X22*X31)/DETX
      DNH3(3,2)=-(X11*X32-X12*X31)/DETX
      DNH3(3,3)= (X11*X22-X12*X21)/DETX
*
* checking the diagonalization
*
*      print*,'(1,0,0) (0,1,0) (0,0,1) ?'
*      print*,DNH3(1,1)*X11+DNH3(1,2)*X21+DNH3(1,3)*X31  
*      print*,DNH3(1,1)*X12+DNH3(1,2)*X22+DNH3(1,3)*X32  
*      print*,DNH3(1,1)*X13+DNH3(1,2)*X23+DNH3(1,3)*X33  
*      print*,DNH3(2,1)*X11+DNH3(2,2)*X21+DNH3(2,3)*X31  
*      print*,DNH3(2,1)*X12+DNH3(2,2)*X22+DNH3(2,3)*X32  
*      print*,DNH3(2,1)*X13+DNH3(2,2)*X23+DNH3(2,3)*X33  
*      print*,DNH3(3,1)*X11+DNH3(3,2)*X21+DNH3(3,3)*X31  
*      print*,DNH3(3,1)*X12+DNH3(3,2)*X22+DNH3(3,3)*X32  
*      print*,DNH3(3,1)*X13+DNH3(3,2)*X23+DNH3(3,3)*X33  
*
*-----------------------------------------------------------------------
*--> Calculating DNH4(I,J)=s[(s-M_I^2) d_IJ + i*PI^hat_IJ]^-1 
      Y11=X11
      Y12=X12
      Y13=X13
      Y14=XI*G0_HATPI(1)/S
       IF (IKILL.EQ.1) Y14=DCMPLX(0.D0,0.D0)
      Y21=X21
      Y22=X22
      Y23=X23
      Y24=XI*G0_HATPI(2)/S
       IF (IKILL.EQ.1) Y24=DCMPLX(0.D0,0.D0)
      Y31=X31
      Y32=X32
      Y33=X33
      Y34=XI*G0_HATPI(3)/S
       IF (IKILL.EQ.1) Y34=DCMPLX(0.D0,0.D0)
      Y41=Y14
      Y42=Y24
      Y43=Y34
      Y44=1.D0+XI*G0_HATPI(4)/S
* 
      DETY= Y11*Y22*Y33*Y44 - Y11*Y22*Y34*Y43 - Y11*Y32*Y23*Y44 
     .    + Y11*Y32*Y24*Y43 + Y11*Y42*Y23*Y34 - Y11*Y42*Y24*Y33
     .    - Y21*Y12*Y33*Y44 + Y21*Y12*Y34*Y43 + Y21*Y32*Y13*Y44 
     .    - Y21*Y32*Y14*Y43 - Y21*Y42*Y13*Y34 + Y21*Y42*Y14*Y33 
     .    + Y31*Y12*Y23*Y44 - Y31*Y12*Y24*Y43 - Y31*Y22*Y13*Y44 
     .    + Y31*Y22*Y14*Y43 + Y31*Y42*Y13*Y24 - Y31*Y42*Y14*Y23 
     .    - Y41*Y12*Y23*Y34 + Y41*Y12*Y24*Y33 + Y41*Y22*Y13*Y34
     .    - Y41*Y22*Y14*Y33 - Y41*Y32*Y13*Y24 + Y41*Y32*Y14*Y23
*      print*,'DETX,DETY',DETX,DETY

      DNH4(1,1)=( Y22*Y33*Y44 - Y22*Y34*Y43 - Y32*Y23*Y44 
     .          +Y32*Y24*Y43 + Y42*Y23*Y34 - Y42*Y24*Y33)/DETY
      DNH4(1,2)=(-Y12*Y33*Y44 + Y12*Y34*Y43 + Y32*Y13*Y44 
     .          -Y32*Y14*Y43 - Y42*Y13*Y34 + Y42*Y14*Y33)/DETY
      DNH4(1,3)=( Y12*Y23*Y44 - Y12*Y24*Y43 - Y22*Y13*Y44 
     .          +Y22*Y14*Y43 + Y42*Y13*Y24 - Y42*Y14*Y23)/DETY
      DNH4(1,4)=(-Y12*Y23*Y34 + Y12*Y24*Y33 + Y22*Y13*Y34 
     .          -Y22*Y14*Y33 - Y32*Y13*Y24 + Y32*Y14*Y23)/DETY
      DNH4(2,1)=(-Y21*Y33*Y44 + Y21*Y34*Y43 + Y31*Y23*Y44 
     .          -Y31*Y24*Y43 - Y41*Y23*Y34 + Y41*Y24*Y33)/DETY
      DNH4(2,2)=( Y11*Y33*Y44 - Y11*Y34*Y43 - Y31*Y13*Y44 
     .          +Y31*Y14*Y43 + Y41*Y13*Y34 - Y41*Y14*Y33)/DETY
      DNH4(2,3)=(-Y11*Y23*Y44 + Y11*Y24*Y43 + Y21*Y13*Y44 
     .          -Y21*Y14*Y43 - Y41*Y13*Y24 + Y41*Y14*Y23)/DETY
      DNH4(2,4)=( Y11*Y23*Y34 - Y11*Y24*Y33 - Y21*Y13*Y34 
     .          +Y21*Y14*Y33 + Y31*Y13*Y24 - Y31*Y14*Y23)/DETY
      DNH4(3,1)=( Y21*Y32*Y44 - Y21*Y34*Y42 - Y31*Y22*Y44 
     .          +Y31*Y24*Y42 + Y41*Y22*Y34 - Y41*Y24*Y32)/DETY
      DNH4(3,2)=(-Y11*Y32*Y44 + Y11*Y34*Y42 + Y31*Y12*Y44 
     .          -Y31*Y14*Y42 - Y41*Y12*Y34 + Y41*Y14*Y32)/DETY
      DNH4(3,3)=( Y11*Y22*Y44 - Y11*Y24*Y42 - Y21*Y12*Y44 
     .          +Y21*Y14*Y42 + Y41*Y12*Y24 - Y41*Y14*Y22)/DETY
      DNH4(3,4)=(-Y11*Y22*Y34 + Y11*Y24*Y32 + Y21*Y12*Y34 
     .          -Y21*Y14*Y32 - Y31*Y12*Y24 + Y31*Y14*Y22)/DETY
      DNH4(4,1)=(-Y21*Y32*Y43 + Y21*Y33*Y42 + Y31*Y22*Y43 
     .          -Y31*Y23*Y42 - Y41*Y22*Y33 + Y41*Y23*Y32)/DETY
      DNH4(4,2)=( Y11*Y32*Y43 - Y11*Y33*Y42 - Y31*Y12*Y43 
     .          +Y31*Y13*Y42 + Y41*Y12*Y33 - Y41*Y13*Y32)/DETY
      DNH4(4,3)=(-Y11*Y22*Y43 + Y11*Y23*Y42 + Y21*Y12*Y43 
     .          -Y21*Y13*Y42 - Y41*Y12*Y23 + Y41*Y13*Y22)/DETY
      DNH4(4,4)=( Y11*Y22*Y33 - Y11*Y23*Y32 - Y21*Y12*Y33 
     .          +Y21*Y13*Y32 + Y31*Y12*Y23 - Y31*Y13*Y22)/DETY
*      do i=1,3
*       do j=1,3
*        print*,dnh3(i,j),dnh4(i,j)
*     .        ,cdabs(dnh3(i,j)-dnh4(i,j))/cdabs(dnh3(i,j))
*       enddo
*      enddo
*      do i=1,4
**       print*,dnh4(i,4),dnh4(4,i),cdabs(dnh4(4,i)/dnh4(i,i))
*       print*,dnh4(i,4),cdabs(dnh4(4,i)),cdabs(dnh4(4,i)/dnh4(i,i))
*      enddo
*      do i=1,4
*       do j=i,4
*        print*,'DNH4[',i,j,']=',dnh4(i,j),cdabs(dnh4(i,j))
*       enddo
*      enddo
* Filling the COMMON block /HC_CAUX/ CAUX_H
      CAUX_H(100)=DNH4(1,1)
      CAUX_H(101)=DNH4(1,2)
      CAUX_H(102)=DNH4(1,3)
      CAUX_H(103)=DNH4(1,4)
      CAUX_H(104)=DNH4(2,1)
      CAUX_H(105)=DNH4(2,2)
      CAUX_H(106)=DNH4(2,3)
      CAUX_H(107)=DNH4(2,4)
      CAUX_H(108)=DNH4(3,1)
      CAUX_H(109)=DNH4(3,2)
      CAUX_H(110)=DNH4(3,3)
      CAUX_H(111)=DNH4(3,4)
      CAUX_H(112)=DNH4(4,1)
      CAUX_H(113)=DNH4(4,2)
      CAUX_H(114)=DNH4(4,3)
      CAUX_H(115)=DNH4(4,4)
*-----------------------------------------------------------------------
*--> Calculating DCH2(I,J)=s[(s-M_I^2) d_IJ + i*PI^hat_IJ]^-1
      Z11=1.D0-MCH**2/S+XI*CH_HATPI(1,1)/S
      Z12=XI*CH_HATPI(1,2)/S
       IF (IKILL.EQ.1) Z12=DCMPLX(0.D0,0.D0)
      Z21=XI*CH_HATPI(2,1)/S
       IF (IKILL.EQ.1) Z21=DCMPLX(0.D0,0.D0)
      Z22=1.D0+XI*HATPI(2,2)/S
*      print*,z11,z12,z21,z22
*
      DETZ=Z11*Z22 - Z12*Z21
*      print*,detz
*
      DCH2(1,1)= Z22/DETZ
      DCH2(1,2)=-Z12/DETZ
      DCH2(2,1)=-Z21/DETZ
      DCH2(2,2)= Z11/DETZ
*
*      do i=1,2
*       do j=1,2
*        print*,'DCH2[',i,j,']=',dch2(i,j),cdabs(dch2(i,j))
*       enddo
*      enddo
* Filling the COMMON block /HC_CAUX/ CAUX_H
      CAUX_H(116)=DCH2(1,1)
      CAUX_H(117)=DCH2(1,2)
      CAUX_H(118)=DCH2(2,1)
      CAUX_H(119)=DCH2(2,2)
*-----------------------------------------------------------------------
*  3   FORMAT(1X,3(2X,E12.6,2X))
*
      RETURN
      END

      SUBROUTINE HSELFIJ(S,I,J,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,HATPI)
**********************************************************************
*
* OUTPUT : HATPI(I,J) = Im[\hat{\Pi}_{ij}(s)] absorptic part of the 
*          propagator. Note the diagonal element satisfies:
*          Im[\hat{\Pi}_{ii}] = M_Hi Gamma_Hi
*
**********************************************************************
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
*Input Arrays
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)  
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRN(NMNH,3,3)
*-----------------------------------------------------------------------
*Local
      INTEGER NSHC_INDEX
      REAL*8 K1,K2,NC
      REAL*8 N11,N12,N13,N14,N22,N23,N24,N33,N34,N44   
      COMPLEX*16 GI_S,GI_P,GJ_S,GJ_P,GIC,GJC
*-----------------------------------------------------------------------
      PI=2.D0*DASIN(1.D0)
*---> running alpha_s and b- and t- quark mass at COM energy s
*     : S > MB^pole assumed
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      IF(DSQRT(S).LE.MTPOLE_H) THEN                               ! AS(S)
       AS_S   = ASMZ_H/(1.D0+B5*ASMZ_H*DLOG(S/MZ_H**2))
      ELSE
       AS_S   = ASMT_H/(1.D0+B6*ASMT_H*DLOG(S/MTPOLE_H**2))
      ENDIF
      IF(DSQRT(S).LE.MTPOLE_H) THEN                               ! MQ(S)
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
      ELSE
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
      ENDIF
*      print*,'sqrt(S),AS(S),MT(S),MB(S)=',DSQRT(S),AS_S,MT_S,MB_S

*--->b-quark
      SF=1.D0
      NC=3.D0
      K1=MB_S**2/S
      K2=MB_S**2/S
      GI  =MB_S/MBMT_H*DREAL(NHC_H(16,I))
      GJ  =MB_S/MBMT_H*DREAL(NHC_H(16,J))
      GI_S=NHC_H(17,I)
      GI_P=NHC_H(18,I)
      GJ_S=NHC_H(17,J)
      GJ_P=NHC_H(18,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,BBB)
      BB=(1.D0+5.67D0*AS_S/PI)*BBB
*      if(i.eq.j) print*,'BB',i,bb,hmass_h(i)*gambrn(6,1,i)
*--->t-quark
      SF=1.D0
      NC=3.D0
      K1=MT_S**2/S
      K2=MT_S**2/S
      GI  =MT_S/MTMT_H*DREAL(NHC_H(25,I))
      GJ  =MT_S/MTMT_H*DREAL(NHC_H(25,J))
      GI_S=NHC_H(26,I)
      GI_P=NHC_H(27,I)
      GJ_S=NHC_H(26,J)
      GJ_P=NHC_H(27,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TTT)
      TT=(1.D0+5.67D0*AS_S/PI)*TTT
*      if(i.eq.j) print*,'TT',i,tt,hmass_h(i)*gambrn(9,1,i)
*--->tau
      SF=1.D0
      NC=1.D0
      K1=MTAU_H**2/S
      K2=MTAU_H**2/S
      GI  =DREAL(NHC_H(7,I))
      GJ  =DREAL(NHC_H(7,J))
      GI_S=NHC_H(8,I)
      GI_P=NHC_H(9,I)
      GJ_S=NHC_H(8,J)
      GJ_P=NHC_H(9,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TATA)
*      if(i.eq.j) print*,'TATA',i,tata,hmass_h(i)*gambrn(3,1,i)
*--->N1 N1
      SF=2.D0
      NC=1.D0
      K1=MN_H(1)**2/S
      K2=MN_H(1)**2/S
      GI  =DREAL(NHC_H(28,I))
      GJ  =DREAL(NHC_H(28,J))
      GI_S=NHC_H(29,I)
      GI_P=NHC_H(30,I)
      GJ_S=NHC_H(29,J)
      GJ_P=NHC_H(30,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N11)
*      if(i.eq.j) print*,'N11',i,n11,hmass_h(i)*gambrn(51,1,i)
*--->N1 N2
      SF=4.D0
      NC=1.D0
      K1=MN_H(1)**2/S
      K2=MN_H(2)**2/S
      GI  =DREAL(NHC_H(40,I))
      GJ  =DREAL(NHC_H(40,J))
      GI_S=NHC_H(41,I)
      GI_P=NHC_H(42,I)
      GJ_S=NHC_H(41,J)
      GJ_P=NHC_H(42,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N12)
*      if(i.eq.j) print*,'N12',i,n12,hmass_h(i)*gambrn(52,1,i)
*--->N1 N3
      SF=4.D0
      NC=1.D0
      K1=MN_H(1)**2/S
      K2=MN_H(3)**2/S
      GI  =DREAL(NHC_H(43,I))
      GJ  =DREAL(NHC_H(43,J))
      GI_S=NHC_H(44,I)
      GI_P=NHC_H(45,I)
      GJ_S=NHC_H(44,J)
      GJ_P=NHC_H(45,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N13)
*      if(i.eq.j) print*,'N13',i,n13,hmass_h(i)*gambrn(53,1,i)
*--->N1 N4
      SF=4.D0
      NC=1.D0
      K1=MN_H(1)**2/S
      K2=MN_H(4)**2/S
      GI  =DREAL(NHC_H(46,I))
      GJ  =DREAL(NHC_H(46,J))
      GI_S=NHC_H(47,I)
      GI_P=NHC_H(48,I)
      GJ_S=NHC_H(47,J)
      GJ_P=NHC_H(48,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N14)
*      if(i.eq.j) print*,'N14',i,n14,hmass_h(i)*gambrn(54,1,i)
*--->N2 N2
      SF=2.D0
      NC=1.D0
      K1=MN_H(2)**2/S
      K2=MN_H(2)**2/S
      GI  =DREAL(NHC_H(31,I))
      GJ  =DREAL(NHC_H(31,J))
      GI_S=NHC_H(32,I)
      GI_P=NHC_H(33,I)
      GJ_S=NHC_H(32,J)
      GJ_P=NHC_H(33,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N22)
*      if(i.eq.j) print*,'N22',i,n22,hmass_h(i)*gambrn(55,1,i)
*--->N2 N3
      SF=4.D0
      NC=1.D0
      K1=MN_H(2)**2/S
      K2=MN_H(3)**2/S
      GI  =DREAL(NHC_H(49,I))
      GJ  =DREAL(NHC_H(49,J))
      GI_S=NHC_H(50,I)
      GI_P=NHC_H(51,I)
      GJ_S=NHC_H(50,J)
      GJ_P=NHC_H(51,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N23)
*      if(i.eq.j) print*,'N23',i,n23,hmass_h(i)*gambrn(56,1,i)
*--->N2 N4
      SF=4.D0
      NC=1.D0
      K1=MN_H(2)**2/S
      K2=MN_H(4)**2/S
      GI  =DREAL(NHC_H(52,I))
      GJ  =DREAL(NHC_H(52,J))
      GI_S=NHC_H(53,I)
      GI_P=NHC_H(54,I)
      GJ_S=NHC_H(53,J)
      GJ_P=NHC_H(54,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N24)
*      if(i.eq.j) print*,'N24',i,n24,hmass_h(i)*gambrn(57,1,i)
*--->N3 N3
      SF=2.D0
      NC=1.D0
      K1=MN_H(3)**2/S
      K2=MN_H(3)**2/S
      GI  =DREAL(NHC_H(34,I))
      GJ  =DREAL(NHC_H(34,J))
      GI_S=NHC_H(35,I)
      GI_P=NHC_H(36,I)
      GJ_S=NHC_H(35,J)
      GJ_P=NHC_H(36,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N33)
*      if(i.eq.j) print*,'N33',i,n33,hmass_h(i)*gambrn(58,1,i)
*--->N3 N4
      SF=4.D0
      NC=1.D0
      K1=MN_H(3)**2/S
      K2=MN_H(4)**2/S
      GI  =DREAL(NHC_H(55,I))
      GJ  =DREAL(NHC_H(55,J))
      GI_S=NHC_H(56,I)
      GI_P=NHC_H(57,I)
      GJ_S=NHC_H(56,J)
      GJ_P=NHC_H(57,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N34)
*      if(i.eq.j) print*,'N34',i,n34,hmass_h(i)*gambrn(59,1,i)
*--->N4 N4
      SF=2.D0
      NC=1.D0
      K1=MN_H(4)**2/S
      K2=MN_H(4)**2/S
      GI  =DREAL(NHC_H(37,I))
      GJ  =DREAL(NHC_H(37,J))
      GI_S=NHC_H(38,I)
      GI_P=NHC_H(39,I)
      GJ_S=NHC_H(38,J)
      GJ_P=NHC_H(39,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,N44)
*      if(i.eq.j) print*,'N44',i,n44,hmass_h(i)*gambrn(60,1,i)
*--->C1+ C1-
      SF=1.D0
      NC=1.D0
      K1=MC_H(1)**2/S
      K2=MC_H(1)**2/S
      GI  =DREAL(NHC_H(58,I))
      GJ  =DREAL(NHC_H(58,J))
      GI_S=NHC_H(59,I)
      GI_P=NHC_H(60,I)
      GJ_S=NHC_H(59,J)
      GJ_P=NHC_H(60,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,C11)
*      if(i.eq.j) print*,'C11',i,c11,hmass_h(i)*gambrn(61,1,i)
*--->C1+ C2-
      SF=1.D0
      NC=1.D0
      K1=MC_H(1)**2/S
      K2=MC_H(2)**2/S
      GI  =DREAL(NHC_H(61,I))
      GJ  =DREAL(NHC_H(61,J))
      GI_S=NHC_H(62,I)
      GI_P=NHC_H(63,I)
      GJ_S=NHC_H(62,J)
      GJ_P=NHC_H(63,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,C12)
*      if(i.eq.j) print*,'C12',i,c12,hmass_h(i)*gambrn(62,1,i)
*--->C2+ C1-
      SF=1.D0
      NC=1.D0
      K1=MC_H(2)**2/S
      K2=MC_H(1)**2/S
      GI  =DREAL(NHC_H(64,I))
      GJ  =DREAL(NHC_H(64,J))
      GI_S=NHC_H(65,I)
      GI_P=NHC_H(66,I)
      GJ_S=NHC_H(65,J)
      GJ_P=NHC_H(66,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,C21)
*      if(i.eq.j) print*,'C21',i,c21,hmass_h(i)*gambrn(63,1,i)
*--->C2+ C2-
      SF=1.D0
      NC=1.D0
      K1=MC_H(2)**2/S
      K2=MC_H(2)**2/S
      GI  =DREAL(NHC_H(67,I))
      GJ  =DREAL(NHC_H(67,J))
      GI_S=NHC_H(68,I)
      GI_P=NHC_H(69,I)
      GJ_S=NHC_H(68,J)
      GJ_P=NHC_H(69,J)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,C22)
*      if(i.eq.j) print*,'C22',i,c22,hmass_h(i)*gambrn(64,1,i)
*--->Z Z
      SF=1.D0
      K1=MZ_H**2/S
      MI=HMASS_H(I)
      MJ=HMASS_H(J)
      GI=DREAL(NHC_H(70,I))
      GJ=DREAL(NHC_H(70,J))
      CALL HATPI_VV(S,SF,K1,GW_H,MW_H,MI,MJ,GI,GJ,MZ_H,ZZ)
*      if(i.eq.j) print*,'ZZ',i,zz,hmass_h(i)*gambrn(11,1,i)
*--->W W
      SF=2.D0
      K1=MW_H**2/S
      MI=HMASS_H(I)
      MJ=HMASS_H(J)
      GI=DREAL(NHC_H(70,I))
      GJ=DREAL(NHC_H(70,J))
      CALL HATPI_VV(S,SF,K1,GW_H,MW_H,MI,MJ,GI,GJ,MW_H,WW)
*      if(i.eq.j) print*,'WW',i,ww,hmass_h(i)*gambrn(10,1,i)
*--->H1 Z
      KK=1
      K1=MZ_H**2/S
      K2=HMASS_H(KK)**2/S
      MI=HMASS_H(I)
      MJ=HMASS_H(J)
      GI=AEPS(I,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(I,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(I,KK,3)*DREAL(NHC_H(70,3))
      GJ=AEPS(J,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(J,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(J,KK,3)*DREAL(NHC_H(70,3))
      MK=HMASS_H(KK)
      CALL HATPI_HZ(S,K1,K2,GW_H,MW_H,MI,MJ,GI,GJ,MK,MZ_H,H1Z)
*      if(i.eq.j) print*,'H1Z',i,h1z,hmass_h(i)*gambrn(12,1,i)
*--->H2 Z
      KK=2
      K1=MZ_H**2/S
      K2=HMASS_H(KK)**2/S
      MI=HMASS_H(I)
      MJ=HMASS_H(J)
      GI=AEPS(I,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(I,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(I,KK,3)*DREAL(NHC_H(70,3))
      GJ=AEPS(J,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(J,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(J,KK,3)*DREAL(NHC_H(70,3))
      MK=HMASS_H(KK)
      CALL HATPI_HZ(S,K1,K2,GW_H,MW_H,MI,MJ,GI,GJ,MK,MZ_H,H2Z)
*      if(i.eq.j) print*,'H2Z',i,h2z,hmass_h(i)*gambrn(13,1,i)
*--->H3 Z
      KK=3
      K1=MZ_H**2/S
      K2=HMASS_H(KK)**2/S
      MI=HMASS_H(I)
      MJ=HMASS_H(J)
      GI=AEPS(I,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(I,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(I,KK,3)*DREAL(NHC_H(70,3))
      GJ=AEPS(J,KK,1)*DREAL(NHC_H(70,1))
     .  +AEPS(J,KK,2)*DREAL(NHC_H(70,2))
     .  +AEPS(J,KK,3)*DREAL(NHC_H(70,3))
      MK=HMASS_H(KK)
      CALL HATPI_HZ(S,K1,K2,GW_H,MW_H,MI,MJ,GI,GJ,MK,MZ_H,H3Z)
*      if(i.eq.j) print*,'H3Z',i,h3z
*--->H W
      K1=MW_H**2/S
      K2=MCH**2/S
      IF((DSQRT(K1)+DSQRT(K2)).GT.1.D0) THEN
       CHW=0.D0
      ELSE
       XLAM=(1.D0-K1-K2)**2-4.D0*K1*K2
       CHW=GW_H**2/32.D0/PI/MW_H**2
     .    *DREAL( NHC_H(87,I)*DCONJG(NHC_H(87,J)) )
     .    *DSQRT(XLAM)*(-4.D0*S*MW_H**2+(MW_H**2-MCH**2)**2
     .    +(MW_H**2-MCH**2)*(HMASS_H(I)**2+HMASS_H(J)**2)
     .    +HMASS_H(I)**2*HMASS_H(J)**2)
      ENDIF
*      if(i.eq.j) print*,'CHW',i,chw
*--->H1 H1
      DKL=1.D0
      NC=1.D0
      IK=1
      IL=1
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H11)
*      if(i.eq.j) print*,'H11',i,h11,hmass_h(i)*gambrn(14,1,i)
*--->H1 H2
      DKL=0.D0
      NC=1.D0
      IK=1
      IL=2
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H12)
*      if(i.eq.j) print*,'H12',i,h12,hmass_h(i)*gambrn(15,1,i)
*--->H1 H3
      DKL=0.D0
      NC=1.D0
      IK=1
      IL=3
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H13)
*      if(i.eq.j) print*,'H13',i,h13
*--->H2 H2
      DKL=1.D0
      NC=1.D0
      IK=2
      IL=2
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H22)
*      if(i.eq.j) print*,'H22',i,h22
*--->H2 H3
      DKL=0.D0
      NC=1.D0
      IK=2
      IL=3
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H23)
*      if(i.eq.j) print*,'H23',i,h23
*--->H3 H3
      DKL=1.D0
      NC=1.D0
      IK=3
      IL=3
      K1=HMASS_H(IK)**2/S
      K2=HMASS_H(IL)**2/S
      II=NSHC_INDEX(I,IK,IL)
       IF( (II.EQ.1) .OR. (II.EQ.7) .OR. (II.EQ. 10) ) THEN
        SFI=3.D0
       ELSEIF(II.EQ.5) THEN
        SFI=1.D0
       ELSE
        SFI=2.D0
       ENDIF
      JJ=NSHC_INDEX(J,IK,IL)
       IF( (JJ.EQ.1) .OR. (JJ.EQ.7) .OR. (JJ.EQ. 10) ) THEN
        SFJ=3.D0
       ELSEIF(JJ.EQ.5) THEN
        SFJ=1.D0
       ELSE
        SFJ=2.D0
       ENDIF
      SF=SFI*SFJ/(1.D0+DKL)
      GIC=DCMPLX(SHC_H(II),0.D0)
      GJC=DCMPLX(SHC_H(JJ),0.D0)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,H33)
*      if(i.eq.j) print*,'H33',i,h33
*--->ST1 ST1*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(71,I)
      GJC=NHC_H(71,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST11)
*      if(i.eq.j) print*,'ST11',i,st11,hmass_h(i)*gambrn(65,1,i)
*--->ST1 ST2*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(72,I)
      GJC=NHC_H(72,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST12)
*      if(i.eq.j) print*,'ST12',i,st12,hmass_h(i)*gambrn(66,1,i)
*--->ST2 ST1*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(73,I)
      GJC=NHC_H(73,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST21)
*      if(i.eq.j) print*,'ST21',i,st21,hmass_h(i)*gambrn(67,1,i)
*--->ST2 ST2*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(74,I)
      GJC=NHC_H(74,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST22)
*      if(i.eq.j) print*,'ST22',i,st22,hmass_h(i)*gambrn(68,1,i)
*--->SB1 SB1*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=1
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(75,I)
      GJC=NHC_H(75,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB11)
*      if(i.eq.j) print*,'SB11',i,sb11,hmass_h(i)*gambrn(69,1,i)
*--->SB1 SB2*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=2
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(76,I)
      GJC=NHC_H(76,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB12)
*      if(i.eq.j) print*,'SB12',i,sb12,hmass_h(i)*gambrn(70,1,i)
*--->SB2 SB1*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=1
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(77,I)
      GJC=NHC_H(77,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB21)
*      if(i.eq.j) print*,'SB21',i,sb21,hmass_h(i)*gambrn(71,1,i)
*--->SB2 SB2*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=2
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(78,I)
      GJC=NHC_H(78,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB22)
*      if(i.eq.j) print*,'SB22',i,sb22,hmass_h(i)*gambrn(72,1,i)
*--->STA1 STA1*
      SF=1.D0
      NC=1.D0
      IK=1
      IL=1
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(79,I)
      GJC=NHC_H(79,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA11)
*      if(i.eq.j) print*,'STA11',i,sta11,hmass_h(i)*gambrn(73,1,i)
*--->STA1 STA2*
      SF=1.D0
      NC=1.D0
      IK=1
      IL=2
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(80,I)
      GJC=NHC_H(80,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA12)
*      if(i.eq.j) print*,'STA12',i,sta12,hmass_h(i)*gambrn(74,1,i)
*--->STA2 STA1*
      SF=1.D0
      NC=1.D0
      IK=2
      IL=1
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(81,I)
      GJC=NHC_H(81,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA21)
*      if(i.eq.j) print*,'STA21',i,sta21,hmass_h(i)*gambrn(75,1,i)
*--->STA2 STA2*
      SF=1.D0
      NC=1.D0
      IK=2
      IL=2
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(82,I)
      GJC=NHC_H(82,J)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA22)
*      if(i.eq.j) print*,'STA22',i,sta22,hmass_h(i)*gambrn(76,1,i)
*
*---> Collecting all
*============================================================================
      HATPI=BB+TT+TATA                                 ! quarks and leptons
     .     +N11+N12+N13+N14+N22+N23+N24+N33+N34+N44    ! neutralinos
     .     +C11+C12+C21+C22                            ! charginos
     .     +ZZ+WW                                      ! vector bosons
     .     +H1Z+H2Z+H3Z+CHW                            ! Higgs and Z/W bosons
     .     +H11+H12+H13+H22+H23+H33                    ! Higgs bosons
     .     +ST11+ST12+ST21+ST22                        ! stops
     .     +SB11+SB12+SB21+SB22                        ! sbottoms
     .     +STA11+STA12+STA21+STA22                    ! staus

*      HATPI=TATA
*============================================================================
*
      RETURN
      END
  
      INTEGER FUNCTION NSHC_INDEX(I,J,K)
************************************************************************
*
************************************************************************
      IF((I.EQ.1) .AND. (J.EQ.1) .AND. (K.EQ.1)) NSHC_INDEX=10
      IF((I.EQ.1) .AND. (J.EQ.1) .AND. (K.EQ.2)) NSHC_INDEX=9
      IF((I.EQ.1) .AND. (J.EQ.1) .AND. (K.EQ.3)) NSHC_INDEX=6
      IF((I.EQ.1) .AND. (J.EQ.2) .AND. (K.EQ.1)) NSHC_INDEX=9
      IF((I.EQ.1) .AND. (J.EQ.2) .AND. (K.EQ.2)) NSHC_INDEX=8
      IF((I.EQ.1) .AND. (J.EQ.2) .AND. (K.EQ.3)) NSHC_INDEX=5
      IF((I.EQ.1) .AND. (J.EQ.3) .AND. (K.EQ.1)) NSHC_INDEX=6
      IF((I.EQ.1) .AND. (J.EQ.3) .AND. (K.EQ.2)) NSHC_INDEX=5
      IF((I.EQ.1) .AND. (J.EQ.3) .AND. (K.EQ.3)) NSHC_INDEX=3
*
      IF((I.EQ.2) .AND. (J.EQ.1) .AND. (K.EQ.1)) NSHC_INDEX=9
      IF((I.EQ.2) .AND. (J.EQ.1) .AND. (K.EQ.2)) NSHC_INDEX=8
      IF((I.EQ.2) .AND. (J.EQ.1) .AND. (K.EQ.3)) NSHC_INDEX=5
      IF((I.EQ.2) .AND. (J.EQ.2) .AND. (K.EQ.1)) NSHC_INDEX=8
      IF((I.EQ.2) .AND. (J.EQ.2) .AND. (K.EQ.2)) NSHC_INDEX=7
      IF((I.EQ.2) .AND. (J.EQ.2) .AND. (K.EQ.3)) NSHC_INDEX=4
      IF((I.EQ.2) .AND. (J.EQ.3) .AND. (K.EQ.1)) NSHC_INDEX=5
      IF((I.EQ.2) .AND. (J.EQ.3) .AND. (K.EQ.2)) NSHC_INDEX=4
      IF((I.EQ.2) .AND. (J.EQ.3) .AND. (K.EQ.3)) NSHC_INDEX=2
*
      IF((I.EQ.3) .AND. (J.EQ.1) .AND. (K.EQ.1)) NSHC_INDEX=6
      IF((I.EQ.3) .AND. (J.EQ.1) .AND. (K.EQ.2)) NSHC_INDEX=5
      IF((I.EQ.3) .AND. (J.EQ.1) .AND. (K.EQ.3)) NSHC_INDEX=3
      IF((I.EQ.3) .AND. (J.EQ.2) .AND. (K.EQ.1)) NSHC_INDEX=5
      IF((I.EQ.3) .AND. (J.EQ.2) .AND. (K.EQ.2)) NSHC_INDEX=4
      IF((I.EQ.3) .AND. (J.EQ.2) .AND. (K.EQ.3)) NSHC_INDEX=2
      IF((I.EQ.3) .AND. (J.EQ.3) .AND. (K.EQ.1)) NSHC_INDEX=3
      IF((I.EQ.3) .AND. (J.EQ.3) .AND. (K.EQ.2)) NSHC_INDEX=2
      IF((I.EQ.3) .AND. (J.EQ.3) .AND. (K.EQ.3)) NSHC_INDEX=1
*
      RETURN
      END

      REAL*8 FUNCTION AEPS(I,J,K)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      IF    ( (I.EQ.1) .AND. (J.EQ.2) .AND. (K.EQ.3) ) THEN
       AEPS=1.D0
      ELSEIF( (I.EQ.1) .AND. (J.EQ.3) .AND. (K.EQ.2) ) THEN
       AEPS=-1.D0
      ELSEIF( (I.EQ.2) .AND. (J.EQ.3) .AND. (K.EQ.1) ) THEN
       AEPS=1.D0
      ELSEIF( (I.EQ.2) .AND. (J.EQ.1) .AND. (K.EQ.3) ) THEN
       AEPS=-1.D0
      ELSEIF( (I.EQ.3) .AND. (J.EQ.1) .AND. (K.EQ.2) ) THEN
       AEPS=1.D0
      ELSEIF( (I.EQ.3) .AND. (J.EQ.2) .AND. (K.EQ.1) ) THEN
       AEPS=-1.D0
      ELSE
       AEPS=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,PFF)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*Local
      REAL*8 K1,K2,NC
      COMPLEX*16 GI_S,GI_P,GJ_S,GJ_P,XFF
*-----------------------------------------------------------------------
      IF((DSQRT(K1)+DSQRT(K2)).GT.1.D0) THEN
       PFF=0.D0
       RETURN
      ENDIF
*
      PI=2.D0*DASIN(1.D0)
      XLAM=(1.D0-K1-K2)**2-4.D0*K1*K2
      XFF=S/8.D0/PI*SF*NC*GI*GJ*(
     .(1.D0-K1-K2)*(GI_S*DCONJG(GJ_S)+GI_P*DCONJG(GJ_P))
     .-2.D0*DSQRT(K1*K2)*(GI_S*DCONJG(GJ_S)-GI_P*DCONJG(GJ_P))
     .)*DSQRT(XLAM)
      PFF=DREAL(XFF) ! XFF could be complex in offdiaginal chargino case.
*                      Our purpoose is to calculate HATPI(I,J) where all the
*                      contributions are to be summed after all. Therefore
*                      it's enough to take real part here. Also see HATPI_SS
*
      RETURN
      END

      SUBROUTINE HATPI_VV(S,SF,KV,GW,MW,MI,MJ,GI,GJ,MV,PVV)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*Local
      REAL*8 KV
*-----------------------------------------------------------------------
      IF(2.D0*(DSQRT(KV)).GT.1.D0) THEN
       PVV=0.D0
       RETURN
      ENDIF
*
      PI=2.D0*DASIN(1.D0)
      BETAV=DSQRT(1.D0-4.D0*KV)
      PVV=GW**2*GI*GJ*SF/128.D0/PI/MW**2*BETAV*(
     .-4.D0*MV**2*(2.D0*S-3.D0*MV**2)+2.D0*MV**2*(MI**2+MJ**2)
     .+MI**2*MJ**2)
*   
      RETURN
      END

      SUBROUTINE HATPI_HZ(S,K1,K2,GW,MW,MI,MJ,GI,GJ,MK,MZ,PHZ)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*Local
      REAL*8 K1,K2
*-----------------------------------------------------------------------
      IF((DSQRT(K1)+DSQRT(K2)).GT.1.D0) THEN
       PHZ=0.D0
       RETURN
      ENDIF
*
      PI=2.D0*DASIN(1.D0)
      XLAM=(1.D0-K1-K2)**2-4.D0*K1*K2
      PHZ=GW**2/64.D0/PI/MW**2*GI*GJ*DSQRT(XLAM)*(-4.D0*S*MZ**2
     .+(MZ**2-MK**2)**2+(MZ**2-MK**2)*(MI**2+MJ**2)+MI**2*MJ**2)
*
      RETURN
      END

      SUBROUTINE HATPI_SS(S,SF,NC,V,K1,K2,GI,GJ,PSS)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*Local
      REAL*8 K1,K2,NC
      COMPLEX*16 GI,GJ,XSS
*-----------------------------------------------------------------------
      IF((DSQRT(K1)+DSQRT(K2)).GT.1.D0) THEN
       PSS=0.D0
       RETURN
      ENDIF
*
      PI=2.D0*DASIN(1.D0)
      XLAM=(1.D0-K1-K2)**2-4.D0*K1*K2
      XSS=V**2/16.D0/PI*SF*NC*GI*DCONJG(GJ)*DSQRT(XLAM)
      PSS=DREAL(XSS) ! XSS can be complex for offdiagonal sfermions
*
      RETURN
      END


      SUBROUTINE HPP_S(SQRTS,M_C,MCH,STMASS,SBMASS,STAUMASS
     .                ,NCMAX,NHC_H)
************************************************************************
*
* HiggsIH-photon-photon Coupling at sqrt(s)=SQRTS
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
      REAL*8     M_C(2),STMASS(2),SBMASS(2),STAUMASS(2)
*Local
      COMPLEX*16 SPHO,PPHO
      COMPLEX*16 SPP(15),PPP(7)
      COMPLEX*16 HC_FSF,HC_FPF,HC_F0,HC_F1
      REAL*8     NC
*
*=======================================================================
*---> running alpha_s and b-quark mass at SQRTS scale
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
       print*,'HPP_S>> SQRTS = ',sqrts,' is out of range !!!'
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
       print*,'HPP_S>> SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*
*      print*,'SQRTS,AS(SQRTS) =',SQRTS,AS_S
*      print*,' > MQ(SQRTS)    :',MT_S,MB_S,MC_S,MS_S
*-----------------------------------------------------------------------
*=======================================================================
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

      DO IH=1,3
*
* Note that, for the quarks (and leptons) NHC_H(IQ,IH)*V_H/MQMT_H=1 and
* the running quarks masses are only used in the loop functions HC_FSF 
* and HC_FPF.  Check it as:
*      print*,' 1? ',NHC_H(IB,IH)*V_H/MBMT_H
*      print*,' 1? ',NHC_H(IT,IH)*V_H/MTMT_H
*      print*,' 1? ',NHC_H(ICHA,IH)*V_H/MCMT_H
*      print*,' 1? ',NHC_H(ITAU,IH)*V_H/MTAU_H
*
      SPP(1)= 2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+1,IH)*V_H/MBMT_H
     .              *HC_FSF(SQRTS**2/4.D0/MB_S**2)
      SPP(2)= 2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+1,IH)*V_H/MTMT_H
     .              *HC_FSF(SQRTS**2/4.D0/MT_S**2)
      SPP(3)= 2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+1,IH)*V_H/MCMT_H
     .              *HC_FSF(SQRTS**2/4.D0/MC_S**2)
      SPP(4)= 2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+1,IH)*V_H/MTAU_H
     .              *HC_FSF(SQRTS**2/4.D0/MTAU_H**2)
      SPP(5)= 2.D0*NHC_H(IC1,IH)*NHC_H(IC1+1,IH)*V_H/M_C(1)
     .              *HC_FSF(SQRTS**2/4.D0/M_C(1)**2)
      SPP(6)= 2.D0*NHC_H(IC2,IH)*NHC_H(IC2+1,IH)*V_H/M_C(2)
     .              *HC_FSF(SQRTS**2/4.D0/M_C(2)**2)
      SPP(7)=-2.D0*NC*QT**2*NHC_H(IT1,IH)*V_H**2/4.D0/STMASS(1)**2
     .              *HC_F0(SQRTS**2/4.D0/STMASS(1)**2)
      SPP(8)=-2.D0*NC*QT**2*NHC_H(IT2,IH)*V_H**2/4.D0/STMASS(2)**2
     .              *HC_F0(SQRTS**2/4.D0/STMASS(2)**2)
      SPP(9)=-2.D0*NC*QB**2*NHC_H(IB1,IH)*V_H**2/4.D0/SBMASS(1)**2
     .              *HC_F0(SQRTS**2/4.D0/SBMASS(1)**2)
      SPP(10)=-2.D0*NC*QB**2*NHC_H(IB2,IH)*V_H**2/4.D0/SBMASS(2)**2
     .              *HC_F0(SQRTS**2/4.D0/SBMASS(2)**2)
      SPP(11)=-NHC_H(IHV,IH)*HC_F1(SQRTS**2/4.D0/MW_H**2)
      SPP(12)=-NHC_H(ICH,IH)*V_H**2/2.D0/MCH**2
     .              *HC_F0(SQRTS**2/4.D0/MCH**2)
      SPP(13)=-2.D0*NHC_H(ITU1,IH)*V_H**2/4.D0/STAUMASS(1)**2
     .             *HC_F0(SQRTS**2/4.D0/STAUMASS(1)**2)
      SPP(14)=-2.D0*NHC_H(ITU2,IH)*V_H**2/4.D0/STAUMASS(2)**2
     .             *HC_F0(SQRTS**2/4.D0/STAUMASS(2)**2)
      SPP(15)=SPP(1)+SPP(2)+SPP(3)+SPP(4)+SPP(5)+SPP(6)+SPP(7)+SPP(8)
     .       +SPP(9)+SPP(10)+SPP(11)+SPP(12)+SPP(13)+SPP(14)
      SPHO=SPP(15)
*      SPHO=SPP(15)-SPP(3)-SPP(4) ! to compare the results w/o c and tau
*
      PPP(1)= 2.D0*NC*QB**2*NHC_H(IB,IH)*NHC_H(IB+2,IH)*V_H/MBMT_H
     .              *HC_FPF(SQRTS**2/4.D0/MB_S**2)
      PPP(2)= 2.D0*NC*QT**2*NHC_H(IT,IH)*NHC_H(IT+2,IH)*V_H/MTMT_H
     .              *HC_FPF(SQRTS**2/4.D0/MT_S**2)
      PPP(3)= 2.D0*NC*QCHA**2*NHC_H(ICHA,IH)*NHC_H(ICHA+2,IH)*V_H/MCMT_H
     .              *HC_FPF(SQRTS**2/4.D0/MC_S**2)
      PPP(4)= 2.D0*NHC_H(ITAU,IH)*NHC_H(ITAU+2,IH)*V_H/MTAU_H
     .              *HC_FPF(SQRTS**2/4.D0/MTAU_H**2)
      PPP(5)= 2.D0*NHC_H(IC1,IH)*NHC_H(IC1+2,IH)*V_H/M_C(1)
     .              *HC_FPF(SQRTS**2/4.D0/M_C(1)**2)
      PPP(6)= 2.D0*NHC_H(IC2,IH)*NHC_H(IC2+2,IH)*V_H/M_C(2)
     .              *HC_FPF(SQRTS**2/4.D0/M_C(2)**2)
      PPP(7)=PPP(1)+PPP(2)+PPP(3)+PPP(4)+PPP(5)+PPP(6)
      PPHO=PPP(7)
*      PPHO=PPP(7)-PPP(3)-PPP(4) ! to compare the results w/o c and tau
*
*Filling the COMMON block /HC_CAUX/ CAUX_H
      CAUX_H(128+2*IH)=SPHO
      CAUX_H(129+2*IH)=PPHO
*
      ENDDO ! IH
*
 29   format(6(1x,e14.6))
*
      RETURN
      END



      SUBROUTINE HGG_S(SQRTS,STMASS,SBMASS,NCMAX,NHC_H)
************************************************************************
*
* HiggsIH-glue-glue Coupling at sqrt(S)
*
* SGG(1) PGG(1) : bottom
* SGG(2) PGG(2) : top
* SGG(3)        : stop1 stop1
* SGG(4)        : stop2 stop2
* SGG(5)        : sbottom1 sbottom1
* SGG(6)        : sbottom2 sbottom2
* SGG(7) PGG(3) : total
*
*JSL:16/Mar/2012 : The pole quark mass used in the loop functions
*                  HC_FSF and HC_FPF
*
*                  Ref. M.~Spira, A.~Djouadi, D.~Graudenz and P.~M.~Zerwas,
*                  ``Higgs boson production at the LHC,''
*                  Nucl.\ Phys.\ B {\bf 453} (1995) 17 [hep-ph/9504378].
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
      COMPLEX*16 NHC_H(NCMAX,3)
      COMPLEX*16 SGLUE,PGLUE
      COMPLEX*16 HC_FSF,HC_FPF,HC_F0,HC_F1
      REAL*8     STMASS(2),SBMASS(2)
      COMPLEX*16 SGG(7),PGG(3)
*=======================================================================
*---> running alpha_s and b-quark mass at SQRTS scale
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
       print*,'HGG_S>> SQRTS = ',sqrts,' is out of range !!!'
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
       print*,'HGG_S>> SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*
*      print*,'SQRTS,AS(SQRTS) =',SQRTS,AS_S
*      print*,' > MQ(SQRTS)    :',MT_S,MB_S,MC_S,MS_S
*-----------------------------------------------------------------------
*=======================================================================
*
      IBOT = 16
      ITOP = 25
      IT11 = 71
      IT22 = 74
      IB11 = 75
      IB22 = 78
*
*      print*,'filldhpg.f:HGG',mtpole_h,mb_pole,mc_pole
*     .      ,' : ',raux_h(1),raux_h(4)
*
      DO IH=1,3
*
      SGG(1)= NHC_H(IBOT,IH)*NHC_H(IBOT+1,IH)*V_H/MBMT_H
     .                *HC_FSF(SQRTS**2/4.D0/RAUX_H(1)**2)
      SGG(2)= NHC_H(ITOP,IH)*NHC_H(ITOP+1,IH)*V_H/MTMT_H
     .                *HC_FSF(SQRTS**2/4.D0/MTPOLE_H**2)
      SGG(3)=-NHC_H(IT11,IH)*V_H**2/4.D0/STMASS(1)**2
     .                *HC_F0(SQRTS**2/4.D0/STMASS(1)**2)
      SGG(4)=-NHC_H(IT22,IH)*V_H**2/4.D0/STMASS(2)**2
     .                *HC_F0(SQRTS**2/4.D0/STMASS(2)**2)
      SGG(5)=-NHC_H(IB11,IH)*V_H**2/4.D0/SBMASS(1)**2
     .                *HC_F0(SQRTS**2/4.D0/SBMASS(1)**2)
      SGG(6)=-NHC_H(IB22,IH)*V_H**2/4.D0/SBMASS(2)**2
     .                *HC_F0(SQRTS**2/4.D0/SBMASS(2)**2)
      SGG(7)=SGG(1)+SGG(2)+SGG(3)+SGG(4)+SGG(5)+SGG(6)
      SGLUE=SGG(7)

      PGG(1)= NHC_H(IBOT,IH)*NHC_H(IBOT+2,IH)*V_H/MBMT_H
     .                *HC_FPF(SQRTS**2/4.D0/RAUX_H(1)**2)
      PGG(2)= NHC_H(ITOP,IH)*NHC_H(ITOP+2,IH)*V_H/MTMT_H
     .                *HC_FPF(SQRTS**2/4.D0/MTPOLE_H**2)
      PGG(3)=PGG(1)+PGG(2)
      PGLUE=PGG(3)
*
*Filling the COMMON block /HC_CAUX/ CAUX_H
      CAUX_H(138+2*IH)=SGLUE
      CAUX_H(139+2*IH)=PGLUE
*
      ENDDO ! IH
*
      RETURN
      END


      SUBROUTINE NGSELFI(S,I,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H      
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,G0_HATPI)
**********************************************************************
*
* OUTPUT : G0_HATPI(I) = Im[\hat{\Pi}_{G^0 H_I)}(s)] absorptic part 
*          of the propagator. 
*          Note G0_HATPI(I=4) = Im[\hat{\Pi}_{G^0 G^0}(s)]
*
* Convention of the G^0 couplings:
*      G0_COUPL( 1) = g^S_{G0-tau-tau}
*      G0_COUPL( 2) = g^P_{G0-tau-tau}
*      G0_COUPL( 3) = g^S_{G0-b-b}
*      G0_COUPL( 4) = g^P_{G0-b-b}
*      G0_COUPL( 5) = g^S_{G0-t-t}
*      G0_COUPL( 6) = g^P_{G0-t-t}
*      G0_COUPL( 7) = g_{G0-stau_1^*-stau_1}
*      G0_COUPL( 8) = g_{G0-stau_1^*-stau_2}
*      G0_COUPL( 9) = g_{G0-stau_2^*-stau_1}
*      G0_COUPL(10) = g_{G0-stau_2^*-stau_2}
*      G0_COUPL(11) = g_{G0-sbot_1^*-sbot_1}
*      G0_COUPL(12) = g_{G0-sbot_1^*-sbot_2}
*      G0_COUPL(13) = g_{G0-sbot_2^*-sbot_1}
*      G0_COUPL(14) = g_{G0-sbot_2^*-sbot_2}
*      G0_COUPL(15) = g_{G0-stop_1^*-stop_1}
*      G0_COUPL(16) = g_{G0-stop_1^*-stop_2}
*      G0_COUPL(17) = g_{G0-stop_2^*-stop_1}
*      G0_COUPL(18) = g_{G0-stop_2^*-stop_2}
*
**********************************************************************
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
*Input Arrays
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*-----------------------------------------------------------------------
*Local
*
      COMPLEX*16 XI,HTAU,HB_MT,HT_MT,HB_S,HT_S,SMIX(2,2)
*
      REAL*8 K1,K2,NC
      COMPLEX*16 GI_S,GI_P,GJ_S,GJ_P,GIC,GJC
*G^0 Couplings
      COMPLEX*16 GAM_LR,GAM_RL,G0_COUPL(18)
*-----------------------------------------------------------------------
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
*Yukawa couplings at mt^pole including Threshold corrections
      HTAU  = DCMPLX(DSQRT(2.D0)*MTAU_H/V_H/CB_H,0.D0) ! use just tree-level coupling
      HB_MT = RAUX_H(27)*CAUX_H(2)
      HT_MT = RAUX_H(24)*CAUX_H(1)
*      print*,'hb@mt^pole ',hb_mt,' @',mtpole_h,' GeV'
*      print*,'ht@mt^pole ',ht_mt,' @',mtpole_h,' GeV'
*
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      IF(DSQRT(S).LE.MTPOLE_H) THEN         ! AS(S) : S > MB^pole assumed
       AS_S   = ASMT_H/(1.D0+B5*ASMT_H*DLOG(S/MTPOLE_H**2))
      ELSE
       AS_S   = ASMT_H/(1.D0+B6*ASMT_H*DLOG(S/MTPOLE_H**2))
      ENDIF
      GS2=4.D0*PI*AS_S
*      print*,'As@sqrts(s)',as_s,dsqrt(s)
*
      HB0=DSQRT(2.D0)*MBMT_H/V_H/CB_H
      HT0=DSQRT(2.D0)*MTMT_H/V_H/SB_H
*      print*,'b-q Yukawa (tree)',hb0,raux_h(23)
*      print*,'t-q Yukawa (tree)',ht0,raux_h(22)
      BB_HB=1.D0/16.D0/PI**2*(9.D0*HB0**2/2.D0+HT0**2/2.D0-8.D0*GS2)
      BT_HT=1.D0/16.D0/PI**2*(9.D0*HT0**2/2.D0+HB0**2/2.D0-8.D0*GS2)
*
      HB_S=HB_MT*(1.D0+2.D0*BB_HB*DLOG(S/MTPOLE_H**2))**0.25D0
      HT_S=HT_MT*(1.D0+2.D0*BT_HT*DLOG(S/MTPOLE_H**2))**0.25D0
*      print*,'hb@sqrts(s)',hb_s,' @',dsqrt(s),' GeV'
*      print*,'ht@sqrts(s)',ht_s,' @',dsqrt(s),' GeV'
*
*G^0 Couplings
      G0_COUPL(1) = DCMPLX(0.D0,0.D0)  ! g^S_{G0-tau-tau}
      G0_COUPL(2) = DCMPLX(1.D0,0.D0)  ! g^P_{G0-tau-tau}
      G0_COUPL(3) = DCMPLX(0.D0,0.D0)  ! g^S_{G0-b-b}
      G0_COUPL(4) = DCMPLX(1.D0,0.D0)  ! g^P_{G0-b-b}
      G0_COUPL(5) = DCMPLX(0.D0,0.D0)  ! g^S_{G0-t-t}
      G0_COUPL(6) =-DCMPLX(1.D0,0.D0)  ! g^P_{G0-t-t}
      DO ISF_TYPE=1,3 ! (1=Stau, 2=Sbottom, 3=Stop)
       IF(ISF_TYPE.EQ.1) THEN
        GAM_LR    =-XI*DCONJG(HTAU)*(CB_H*DCONJG(ATAU_H)-SB_H*MU_H)
     .             /DSQRT(2.D0)
        GAM_RL    = DCONJG(GAM_LR)
        SMIX(1,1) = STAUMIX_H(1,1)
        SMIX(1,2) = STAUMIX_H(1,2)
        SMIX(2,1) = STAUMIX_H(2,1)
        SMIX(2,2) = STAUMIX_H(2,2)
        NSF       = 7
       ENDIF ! ISF_TYPE.EQ.1
       IF(ISF_TYPE.EQ.2) THEN
        GAM_LR    =-XI*DCONJG(HB_S)*(CB_H*DCONJG(AB_H)  -SB_H*MU_H)
     .             /DSQRT(2.D0)
        GAM_RL    = DCONJG(GAM_LR)
        SMIX(1,1) = SBMIX_H(1,1)
        SMIX(1,2) = SBMIX_H(1,2)
        SMIX(2,1) = SBMIX_H(2,1)
        SMIX(2,2) = SBMIX_H(2,2)
        NSF       = 11
       ENDIF ! ISF_TYPE.EQ.2
       IF(ISF_TYPE.EQ.3) THEN
        GAM_LR    = XI*DCONJG(HT_S)*(SB_H*DCONJG(AT_H)  -CB_H*MU_H)
     .             /DSQRT(2.D0)
        GAM_RL    = DCONJG(GAM_LR)
        SMIX(1,1) = STMIX_H(1,1)
        SMIX(1,2) = STMIX_H(1,2)
        SMIX(2,1) = STMIX_H(2,1)
        SMIX(2,2) = STMIX_H(2,2)
        NSF       = 15
       ENDIF ! ISF_TYPE.EQ.3
*g_[G^0 - (f^tilde_ISF)^* - f^tilde_JSF]
       DO ISF=1,2
        DO JSF=1,2
         IF(ISF.EQ.1 .AND. JSF.EQ.1) IJ=0
         IF(ISF.EQ.1 .AND. JSF.EQ.2) IJ=1
         IF(ISF.EQ.2 .AND. JSF.EQ.1) IJ=2
         IF(ISF.EQ.2 .AND. JSF.EQ.2) IJ=3
*         print*,nsf+ij
         G0_COUPL(NSF+IJ)
     .    = GAM_LR*DCONJG(SMIX(1,ISF))*SMIX(2,JSF)/V_H
     .     +GAM_RL*DCONJG(SMIX(2,ISF))*SMIX(1,JSF)/V_H
        ENDDO ! JSF
       ENDDO ! ISF
*
      ENDDO ! ISF_TYPE
*
*      if (i.eq.1) then
*       do icoupl=1,18
*        print*,icoupl,g0_coupl(icoupl)
*       enddo
*      endif
*-----------------------------------------------------------------------
*---> running alpha_s and b- and t- quark mass at COM energy s
      IF(DSQRT(S).LE.MTPOLE_H) THEN  ! MQ(S): S > MB^pole assumed
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
      ELSE
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
      ENDIF
*      print*,'sqrt(S),AS(S),MT(S),MB(S)=',DSQRT(S),AS_S,MT_S,MB_S
*      print*,mt_s,cdabs(ht_s/ht_mt)*mtmt_h
*      print*,mb_s,cdabs(hb_s/hb_mt)*mbmt_h
*--->b-quark
      SF=1.D0
      NC=3.D0
      K1=MB_S**2/S
      K2=MB_S**2/S
      GI  =MB_S/MBMT_H*DREAL(NHC_H(16,1))
      GJ  =MB_S/MBMT_H*DREAL(NHC_H(16,1))
      GI_S=NHC_H(17,I)
      GI_P=NHC_H(18,I)
      IF(I.EQ.4) THEN
       GI_S=G0_COUPL(3)
       GI_P=G0_COUPL(4)
      ENDIF
      GJ_S=G0_COUPL(3)
      GJ_P=G0_COUPL(4)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,BBB)
      BB=(1.D0+5.67D0*AS_S/PI)*BBB
*      print*,'BotBot',i,bb
*--->t-quark
      SF=1.D0
      NC=3.D0
      K1=MT_S**2/S
      K2=MT_S**2/S
      GI  =MT_S/MTMT_H*DREAL(NHC_H(25,1))
      GJ  =MT_S/MTMT_H*DREAL(NHC_H(25,1))
      GI_S=NHC_H(26,I)
      GI_P=NHC_H(27,I)
      IF(I.EQ.4) THEN
       GI_S=G0_COUPL(5)
       GI_P=G0_COUPL(6)
      ENDIF
      GJ_S=G0_COUPL(5)
      GJ_P=G0_COUPL(6)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TTT)
      TT=(1.D0+5.67D0*AS_S/PI)*TTT
*      print*,'TopTop',i,tt
*--->tau
      SF=1.D0
      NC=1.D0
      K1=MTAU_H**2/S
      K2=MTAU_H**2/S
      GI  =DREAL(NHC_H(7,1))
      GJ  =DREAL(NHC_H(7,1))
      GI_S=NHC_H(8,I)
      GI_P=NHC_H(9,I)
      IF(I.EQ.4) THEN
       GI_S=G0_COUPL(1)
       GI_P=G0_COUPL(2)
      ENDIF
      GJ_S=G0_COUPL(1)
      GJ_P=G0_COUPL(2)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TATA)
*      print*,'TauTau',i,tata
*--->ST1 ST1*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(71,I)
       IF(I.EQ.4) GIC=G0_COUPL(15)
      GJC=G0_COUPL(15)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST11)
*      print*,'STop11',i,st11
*--->ST1 ST2*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(72,I)
       IF(I.EQ.4) GIC=G0_COUPL(16)
      GJC=G0_COUPL(16)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST12)
*      print*,'STop12',i,st12
*--->ST2 ST1*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(73,I)
       IF(I.EQ.4) GIC=G0_COUPL(17)
      GJC=G0_COUPL(17)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST21)
*      print*,'STop21',i,st21
*--->ST2 ST2*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=STMASS_H(IL)**2/S
      GIC=NHC_H(74,I)
       IF(I.EQ.4) GIC=G0_COUPL(18)
      GJC=G0_COUPL(18)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST22)
*      print*,'STop22',i,st22
*--->SB1 SB1*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=1
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(75,I)
       IF(I.EQ.4) GIC=G0_COUPL(11)
      GJC=G0_COUPL(11)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB11)
*      print*,'SBot11',i,sb11
*--->SB1 SB2*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=2
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(76,I)
       IF(I.EQ.4) GIC=G0_COUPL(12)
      GJC=G0_COUPL(12)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB12)
*      print*,'SBot12',i,sb12
*--->SB2 SB1*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=1
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(77,I)
       IF(I.EQ.4) GIC=G0_COUPL(13)
      GJC=G0_COUPL(13)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB21)
*      print*,'SBot21',i,sb21
*--->SB2 SB2*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=2
      K1=SBMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=NHC_H(78,I)
       IF(I.EQ.4) GIC=G0_COUPL(14)
      GJC=G0_COUPL(14)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SB22)
*      print*,'SBot22',i,sb22
*--->STA1 STA1*
      SF=1.D0
      NC=1.D0
      IK=1
      IL=1
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(79,I)
       IF(I.EQ.4) GIC=G0_COUPL(7)
      GJC=G0_COUPL(7)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA11)
*      print*,'STau11',i,sta11
*--->STA1 STA2*
      SF=1.D0
      NC=1.D0
      IK=1
      IL=2
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(80,I)
       IF(I.EQ.4) GIC=G0_COUPL(8)
      GJC=G0_COUPL(8)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA12)
*      print*,'STau12',i,sta12
*--->STA2 STA1*
      SF=1.D0
      NC=1.D0
      IK=2
      IL=1
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(81,I)
       IF(I.EQ.4) GIC=G0_COUPL(9)
      GJC=G0_COUPL(9)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA21)
*      print*,'STau21',i,sta21
*--->STA2 STA2*
      SF=1.D0
      NC=1.D0
      IK=2
      IL=2
      K1=STAUMASS_H(IK)**2/S
      K2=STAUMASS_H(IL)**2/S
      GIC=NHC_H(82,I)
       IF(I.EQ.4) GIC=G0_COUPL(10)
      GJC=G0_COUPL(10)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,STA22)
*      print*,'STau22',i,sta22
*
*---> Collecting all
*============================================================================
      G0_HATPI=BB+TT+TATA              ! quarks and leptons
     .        +ST11+ST12+ST21+ST22     ! stops
     .        +SB11+SB12+SB21+SB22     ! sbottoms
     .        +STA11+STA12+STA21+STA22 ! staus

*============================================================================
*
      RETURN
      END


      SUBROUTINE CHSELFIJ(S,I,J,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMCH,GAMBRC,CH_HATPI)
**********************************************************************
*
* OUTPUT : CH_HATPI(1,1) = Im[\hat{\Pi}_{H^+ G^+)}(s)] 
*          CH_HATPI(1,2) = Im[\hat{\Pi}_{H^+ G^+)}(s)] 
*          CH_HATPI(2,1) = Im[\hat{\Pi}_{G^+ H^+)}(s)] 
*          CH_HATPI(2,2) = Im[\hat{\Pi}_{G^+ H^+)}(s)] 
*
* Convention of the G^+ couplings:
*      GP_COUPL( 1) = g^S_{G+-nu-tau}
*      GP_COUPL( 2) = g^P_{G+-nu-tau}
*      GP_COUPL( 3) = g^S_{G+-t-b}
*      GP_COUPL( 4) = g^P_{G+-t-b}
*      GP_COUPL( 5) = g_{G^+-stop_1^*-sbot_1}
*      GP_COUPL( 6) = g_{G^+-stop_1^*-sbot_2}
*      GP_COUPL( 7) = g_{G^+-stop_2^*-sbot_1}
*      GP_COUPL( 8) = g_{G^+-stop_2^*-sbot_2}
*      GP_COUPL( 9) = g_{G^+-snu^*-stau_1}
*      GP_COUPL(10) = g_{G^+-snu^*-stau_2}
*
* Convention of the H^+ couplings:
*      HP_COUPL( 1) = g^S_{H+-nu-tau}         : CHC_H( 8)
*      HP_COUPL( 2) = g^P_{H+-nu-tau}         : CHC_H( 9)
*      HP_COUPL( 3) = g^S_{H+-t-b}            : CHC_H(17)
*      HP_COUPL( 4) = g^P_{H+-t-b}            : CHC_H(18)
*      HP_COUPL( 5) = g_{H^+-stop_1^*-sbot_1} : CHC_H(43)
*      HP_COUPL( 6) = g_{H^+-stop_1^*-sbot_2} : CHC_H(44)
*      HP_COUPL( 7) = g_{H^+-stop_2^*-sbot_1} : CHC_H(45)
*      HP_COUPL( 8) = g_{H^+-stop_2^*-sbot_2} : CHC_H(46)
*      HP_COUPL( 9) = g_{H^+-snu^*-stau_1}    : CHC_H(47)
*      HP_COUPL(10) = g_{H^+-snu^*-stau_2}    : CHC_H(48)
*
**********************************************************************
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
*Input Arrays
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*
      REAL*8 GAMBRC(NMCH,3)
*-----------------------------------------------------------------------
*Local
*
      COMPLEX*16 XI,HTAU,HB_MT,HT_MT,HB_S,HT_S
*
      REAL*8 K1,K2,NC
      COMPLEX*16 GI_S,GI_P,GJ_S,GJ_P,GIC,GJC
*G^+ and H^+ Couplings
      COMPLEX*16 GAM_GP(2,2),GAM_HP(2,2)
      COMPLEX*16 GP_COUPL(10),HP_COUPL(10),HG_COUPL(2,10)
*-----------------------------------------------------------------------
      XI=DCMPLX(0.D0,1.D0)
      PI=2.D0*DASIN(1.D0)
*
*Yukawa couplings at mt^pole including Threshold corrections
      HTAU  = DCMPLX(DSQRT(2.D0)*MTAU_H/V_H/CB_H,0.D0) 
*                                        ! use just tree-level coupling
      HB_MT = RAUX_H(27)*CAUX_H(2)
      HT_MT = RAUX_H(24)*CAUX_H(1)
*      print*,'hb@mt^pole ',hb_mt,' @',mtpole_h,' GeV'
*      print*,'ht@mt^pole ',ht_mt,' @',mtpole_h,' GeV'
*
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      IF(DSQRT(S).LE.MTPOLE_H) THEN         ! AS(S) : S > MB^pole assumed
       AS_S   = ASMT_H/(1.D0+B5*ASMT_H*DLOG(S/MTPOLE_H**2))
      ELSE
       AS_S   = ASMT_H/(1.D0+B6*ASMT_H*DLOG(S/MTPOLE_H**2))
      ENDIF
      GS2=4.D0*PI*AS_S
*      print*,'As@sqrts(s)',as_s,dsqrt(s)
*
      HB0=DSQRT(2.D0)*MBMT_H/V_H/CB_H
      HT0=DSQRT(2.D0)*MTMT_H/V_H/SB_H
*      print*,'b-q Yukawa (tree)',hb0,raux_h(23)
*      print*,'t-q Yukawa (tree)',ht0,raux_h(22)
      BB_HB=1.D0/16.D0/PI**2*(9.D0*HB0**2/2.D0+HT0**2/2.D0-8.D0*GS2)
      BT_HT=1.D0/16.D0/PI**2*(9.D0*HT0**2/2.D0+HB0**2/2.D0-8.D0*GS2)
*
      HB_S=HB_MT*(1.D0+2.D0*BB_HB*DLOG(S/MTPOLE_H**2))**0.25D0
      HT_S=HT_MT*(1.D0+2.D0*BT_HT*DLOG(S/MTPOLE_H**2))**0.25D0
*      print*,'hb@sqrts(s)',hb_s,' @',dsqrt(s),' GeV'
*      print*,'ht@sqrts(s)',ht_s,' @',dsqrt(s),' GeV'
*
*---> running alpha_s and b- and t- quark mass at COM energy s
      IF(DSQRT(S).LE.MTPOLE_H) THEN  ! MQ(S): S > MB^pole assumed
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B5/PI)
      ELSE
       MT_S   = MTMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
       MB_S   = MBMT_H*(AS_S/ASMT_H)**(1.D0/B6/PI)
      ENDIF
*      print*,'sqrt(S),AS(S),MT(S),MB(S)=',DSQRT(S),AS_S,MT_S,MB_S
*      print*,mt_s,cdabs(ht_s/ht_mt)*mtmt_h
*      print*,mb_s,cdabs(hb_s/hb_mt)*mbmt_h
*
*-----------------------------------------------------------------------
*G^+ Couplings
      GAM_GP(1,1)=(CDABS(HT_S)**2*SB_H**2
     .            -CDABS(HB_S)**2*CB_H**2)*V_H/DSQRT(2.D0)
     .           +GW_H**2*(CB_H**2-SB_H**2)*V_H/2.D0/DSQRT(2.D0)
      GAM_GP(1,2)=-DCONJG(HB_S)*(CB_H*DCONJG(AB_H)-SB_H*MU_H)
      GAM_GP(2,1)=HT_S*(SB_H*AT_H-CB_H*DCONJG(MU_H))
      GAM_GP(2,2)=DCMPLX(0.D0,0.D0)
*GP_COUPL( 1) = g^S_{G+-nu-tau}
      GP_COUPL( 1) =-DCMPLX(1.D0,0.D0)/2.D0
*GP_COUPL( 2) = g^P_{G+-nu-tau}
      GP_COUPL( 2) = DCMPLX(0.D0,1.D0)/2.D0
*GP_COUPL( 3) = g^S_{G+-t-b}
      GP_COUPL( 3) = DCMPLX(1.D0,0.D0)*(1.D0-MB_S/MT_S)/2.D0
*GP_COUPL( 4) = g^P_{G+-t-b}
      GP_COUPL( 4) = DCMPLX(0.D0,1.D0)*(1.D0+MB_S/MT_S)/2.D0
*GP_COUPL( 5) = g_{G^+-stop_1^*-sbot_1}
*GP_COUPL( 6) = g_{G^+-stop_1^*-sbot_2}
*GP_COUPL( 7) = g_{G^+-stop_2^*-sbot_1}
*GP_COUPL( 8) = g_{G^+-stop_2^*-sbot_2}
      DO ISF=1,2
       DO JSF=1,2
        IF(ISF.EQ.1 .AND. JSF.EQ.1) NSF=5
        IF(ISF.EQ.1 .AND. JSF.EQ.2) NSF=6
        IF(ISF.EQ.2 .AND. JSF.EQ.1) NSF=7
        IF(ISF.EQ.2 .AND. JSF.EQ.2) NSF=8
        GP_COUPL(NSF)=GAM_GP(1,1)*DCONJG(STMIX_H(1,ISF))*SBMIX_H(1,JSF)
     .               +GAM_GP(1,2)*DCONJG(STMIX_H(1,ISF))*SBMIX_H(2,JSF)
     .               +GAM_GP(2,1)*DCONJG(STMIX_H(2,ISF))*SBMIX_H(1,JSF)
     .               +GAM_GP(2,2)*DCONJG(STMIX_H(2,ISF))*SBMIX_H(2,JSF)
        GP_COUPL(NSF)=GP_COUPL(NSF)/V_H
       ENDDO
      ENDDO
*GP_COUPL( 9) = g_{G^+-snu^*-stau_1}
      GP_COUPL( 9)=(-CDABS(HTAU)**2*CB_H**2*V_H/DSQRT(2.D0)
     .              +GW_H**2*(CB_H**2-SB_H**2)*V_H/2.D0/DSQRT(2.D0))
     .            *STAUMIX_H(1,1)/V_H
     .            -DCONJG(HTAU)*(CB_H*DCONJG(ATAU_H)-SB_H*MU_H)
     .            *STAUMIX_H(2,1)/V_H
*GP_COUPL(10) = g_{G^+-snu^*-stau_2}
      GP_COUPL(10)=(-CDABS(HTAU)**2*CB_H**2*V_H/DSQRT(2.D0)
     .              +GW_H**2*(CB_H**2-SB_H**2)*V_H/2.D0/DSQRT(2.D0))
     .            *STAUMIX_H(1,2)/V_H
     .            -DCONJG(HTAU)*(CB_H*DCONJG(ATAU_H)-SB_H*MU_H)
     .            *STAUMIX_H(2,2)/V_H
*
*      if (i.eq.1.and.j.eq.1) then
*       do icoupl=1,10
*        print*,icoupl,gp_coupl(icoupl)
*       enddo
*      endif
*-----------------------------------------------------------------------
*H^+ Couplings
      GAM_HP(1,1)=(CDABS(HT_S)**2+CDABS(HB_S)**2-GW_H**2)
     .           *V_H*SB_H*CB_H/DSQRT(2.D0)
      GAM_HP(1,2)=DCONJG(HB_S)*(SB_H*DCONJG(AB_H)+CB_H*MU_H)
      GAM_HP(2,1)=HT_S*(CB_H*AT_H+SB_H*DCONJG(MU_H))
      GAM_HP(2,2)=HT_S*DCONJG(HB_S)*V_H/DSQRT(2.D0)
*HP_COUPL( 1) = g^S_{H+-nu-tau}         : CHC_H( 8)
*HP_COUPL( 2) = g^P_{H+-nu-tau}         : CHC_H( 9)
      HP_COUPL(1) = DCMPLX(1.D0,0.D0)*TB_H/2.D0
      HP_COUPL(2) =-DCMPLX(0.D0,1.D0)*TB_H/2.D0
*HP_COUPL( 3) = g^S_{H+-t-b}            : CHC_H(17)
*HP_COUPL( 4) = g^P_{H+-t-b}            : CHC_H(18)
      HP_COUPL(3) = DCMPLX((1.D0/TB_H+MB_S/MT_S*TB_H)/2.D0,0.D0)
      HP_COUPL(4) = DCMPLX(0.D0,(1.D0/TB_H-MB_S/MT_S*TB_H)/2.D0)
* ...To include threshold corrections, it's better to use CHC_H
      HP_COUPL(3) = CHC_H(17)
      HP_COUPL(4) = CHC_H(18)
*HP_COUPL( 5) = g_{H^+-stop_1^*-sbot_1} : CHC_H(43)
*HP_COUPL( 6) = g_{H^+-stop_1^*-sbot_2} : CHC_H(44)
*HP_COUPL( 7) = g_{H^+-stop_2^*-sbot_1} : CHC_H(45)
*HP_COUPL( 8) = g_{H^+-stop_2^*-sbot_2} : CHC_H(46)
      DO ISF=1,2
       DO JSF=1,2
        IF(ISF.EQ.1 .AND. JSF.EQ.1) NSF=5
        IF(ISF.EQ.1 .AND. JSF.EQ.2) NSF=6
        IF(ISF.EQ.2 .AND. JSF.EQ.1) NSF=7
        IF(ISF.EQ.2 .AND. JSF.EQ.2) NSF=8
        HP_COUPL(NSF)=GAM_HP(1,1)*DCONJG(STMIX_H(1,ISF))*SBMIX_H(1,JSF)
     .               +GAM_HP(1,2)*DCONJG(STMIX_H(1,ISF))*SBMIX_H(2,JSF)
     .               +GAM_HP(2,1)*DCONJG(STMIX_H(2,ISF))*SBMIX_H(1,JSF)
     .               +GAM_HP(2,2)*DCONJG(STMIX_H(2,ISF))*SBMIX_H(2,JSF)
        HP_COUPL(NSF)=HP_COUPL(NSF)/V_H
       ENDDO
      ENDDO
*HP_COUPL( 9) = g_{H^+-snu^*-stau_1}    : CHC_H(47)
      HP_COUPL( 9)=(CDABS(HTAU)**2-GW_H**2)*V_H*CB_H*SB_H/DSQRT(2.D0)
     .            *STAUMIX_H(1,1)/V_H
     .            +DCONJG(HTAU)*(SB_H*DCONJG(ATAU_H)+CB_H*MU_H)
     .            *STAUMIX_H(2,1)/V_H
*HP_COUPL(10) = g_{H^+-snu^*-stau_2}    : CHC_H(48)
      HP_COUPL(10)=(CDABS(HTAU)**2-GW_H**2)*V_H*CB_H*SB_H/DSQRT(2.D0)
     .            *STAUMIX_H(1,2)/V_H
     .            +DCONJG(HTAU)*(SB_H*DCONJG(ATAU_H)+CB_H*MU_H)
     .            *STAUMIX_H(2,2)/V_H
*
*      if (i.eq.1.and.j.eq.1) then
*       print*,hp_coupl(1),chc_h(8)
*       print*,hp_coupl(2),chc_h(9)
*       print*,hp_coupl(3),chc_h(17)
*       print*,hp_coupl(4),chc_h(18)
*       print*,hp_coupl(5),chc_h(43)
*       print*,hp_coupl(6),chc_h(44)
*       print*,hp_coupl(7),chc_h(45)
*       print*,hp_coupl(8),chc_h(46)
*       print*,hp_coupl(9),chc_h(47)
*       print*,hp_coupl(10),chc_h(48)
*      endif
*-----------------------------------------------------------------------
       DO JJ=1,10
        HG_COUPL(1,JJ)=HP_COUPL(JJ)  ! Charged Higgs
        HG_COUPL(2,JJ)=GP_COUPL(JJ)  ! Charged Goldstone
*        if(i.eq.1.and.j.eq.1) print*,hg_coupl(1,jj),hg_coupl(2,jj)
       ENDDO
*-----------------------------------------------------------------------
*--->t-b-quark
      SF=1.D0
      NC=3.D0
      K1=MT_S**2/S
      K2=MB_S**2/S
      GI  =DCMPLX(-GW_H*MT_S/DSQRT(2.D0)/MW_H,0.D0)
      GJ  =DCMPLX(-GW_H*MT_S/DSQRT(2.D0)/MW_H,0.D0)
*      print*,gi,gj,chc_h(16)*mt_s/mtmt_h
      GI_S=HG_COUPL(I,3)
      GI_P=HG_COUPL(I,4)
      GJ_S=HG_COUPL(J,3)
      GJ_P=HG_COUPL(J,4)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TOPB)
      TOPB=(1.D0+5.67D0*AS_S/PI)*TOPB
*      print*,'TopB',i,j,TOPB,topb*mch/s,GAMBRC(6,1)
*--->tau-nu
      SF=1.D0
      NC=1.D0
      K1=MTAU_H**2/S
      K2=0.D0
      GI  =-GW_H*MTAU_H/DSQRT(2.D0)/MW_H
      GJ  =-GW_H*MTAU_H/DSQRT(2.D0)/MW_H
      GI_S=HG_COUPL(I,1) 
      GI_P=HG_COUPL(I,2)
      GJ_S=HG_COUPL(J,1)
      GJ_P=HG_COUPL(J,2)
      CALL HATPI_FF(S,SF,NC,K1,K2,GI,GJ,GI_S,GI_P,GJ_S,GJ_P,TAUNU)
*      print*,'TauNu',i,j,TAUNU,taunu*mch/s,GAMBRC(3,1)
*     . ,gi**2*mch/8.d0/pi*(cdabs(GI_S)**2+cdabs(GI_P)**2)
*--->ST1 SB1*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=HG_COUPL(I,5)
      GJC=HG_COUPL(J,5)
*      print*,i,j,gic*dconjg(gjc)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST1SB1)
*      print*,'St1Sb1',i,j,ST1SB1,st1sb1/mch
*     .      ,GAMBRC(25+9,1)*cdabs(hp_coupl(5)/chc_h(43))**2
*--->ST1 SB2*
      SF=1.D0
      NC=3.D0
      IK=1
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=HG_COUPL(I,6)
      GJC=HG_COUPL(J,6)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST1SB2)
*      print*,'St1Sb2',i,j,ST1SB2,st1sb2/mch
*     .      ,GAMBRC(25+10,1)*cdabs(hp_coupl(6)/chc_h(44))**2
*--->ST2 SB1*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=1
      K1=STMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=HG_COUPL(I,7)
      GJC=HG_COUPL(J,7)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST2SB1)
*      print*,'St2Sb1',i,j,ST2SB1,st2sb1/mch
*     .      ,GAMBRC(25+11,1)*cdabs(hp_coupl(7)/chc_h(45))**2
*--->ST2 SB2*
      SF=1.D0
      NC=3.D0
      IK=2
      IL=2
      K1=STMASS_H(IK)**2/S
      K2=SBMASS_H(IL)**2/S
      GIC=HG_COUPL(I,8)
      GJC=HG_COUPL(J,8)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,ST2SB2)
*      print*,'St2Sb2',i,j,ST2SB2,st2sb2/mch
*     .      ,GAMBRC(25+12,1)*cdabs(hp_coupl(8)/chc_h(46))**2
*--->SNU3 STA1*
      SF=1.D0
      NC=1.D0
      IK=1
      K1=SNU3MASS_H**2/S
      K2=STAUMASS_H(IK)**2/S
      GIC=HG_COUPL(I,9)
      GJC=HG_COUPL(J,9)
*      print*,i,j,gic*dconjg(gjc)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SNUSTAU1)
*      print*,'SnuStau1',i,j,SNUSTAU1,snustau1/mch,GAMBRC(38,1)
*      print*,'STau11',i,sta11
*--->SNU3 STA2*
      SF=1.D0
      NC=1.D0
      IK=2
      K1=SNU3MASS_H**2/S
      K2=STAUMASS_H(IK)**2/S
      GIC=HG_COUPL(I,10)
      GJC=HG_COUPL(J,10)
      CALL HATPI_SS(S,SF,NC,V_H,K1,K2,GIC,GJC,SNUSTAU2)
*      print*,'SnuStau2',i,j,SNUSTAU2,snustau2/mch,GAMBRC(39,1)
*---> Collecting all
*============================================================================
*      if(i.eq.1 .and. j.eq.1) then
*      print*,'TopB    ',TOPB,topb*mch/s,GAMBRC(6,1)
*      print*,'TauNu   ',TAUNU,taunu*mch/s,GAMBRC(3,1)
*      print*,'St1Sb1  ',ST1SB1,st1sb1/mch
*     .      ,GAMBRC(25+9,1)*cdabs(hp_coupl(5)/chc_h(43))**2
*      print*,'St1Sb2  ',ST1SB2,st1sb2/mch
*     .      ,GAMBRC(25+10,1)*cdabs(hp_coupl(6)/chc_h(44))**2
*      print*,'St2Sb1  ',ST2SB1,st2sb1/mch
*     .      ,GAMBRC(25+11,1)*cdabs(hp_coupl(7)/chc_h(45))**2
*      print*,'St2Sb2',ST2SB2,st2sb2/mch
*     .      ,GAMBRC(25+12,1)*cdabs(hp_coupl(8)/chc_h(46))**2
*      print*,'SnuStau1',SNUSTAU1,snustau1/mch,GAMBRC(38,1)
*      print*,'SnuStau2',SNUSTAU2,snustau2/mch,GAMBRC(39,1)
*      endif
      CH_HATPI=TopB+TauNu                  ! quarks and leptons
     .        +St1Sb1+St1Sb2+St2Sb1+St2Sb2 ! stops and sbottoms
     .        +SnuStau1+SnuStau2           ! snus and staus
*      print*,i,j,ch_hatpi
*============================================================================
*
      RETURN
      END
