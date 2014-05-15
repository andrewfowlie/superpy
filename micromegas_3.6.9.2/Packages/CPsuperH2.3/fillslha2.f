      SUBROUTINE FILLSLHA2(
     . NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .,MCH,HMASS_H,OMIX_H
     .,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     .,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     .,NMNH,GAMBRN,NMCH,GAMBRC)
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
* Input Array:
      REAL*8     SMPARA_H(NSMIN),SSPARA_H(NSSIN)
      INTEGER*8  IFLAG_H(NFLAG)
      REAL*8     MCH,HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
      REAL*8     GAMBRN(NMNH,3,3)   
      REAL*8     GAMBRC(NMCH,3)      
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
*-----------------------------------------------------------------------
* Local :
      COMPLEX*16 VMIX(2,2),FMIX(2,2)
      REAL*8     SMIX(3,4),RTMP,RTMP34(3,4),RTMP44(4,4)
      COMPLEX*16 CTMP
*-----------------------------------------------------------------------
      NOUT=80 
      OPEN(NOUT,FILE='cpsuperh2_slha.out',STATUS='UNKNOWN')
*-----------------------------------------------------------------------
      PI     = 2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
*Head
      WRITE(NOUT,11) '--------------------------------------------------
     .----------------------------'
      WRITE(NOUT,11) 'SUSY Les Houches Accord 2.0 - MSSM Spectrum + Deca
     .ys'
      WRITE(NOUT,11) 'CPsuperH2.3'
      WRITE(NOUT,11) 'J.S. Lee, A. Pilaftsis, M. Carena, S.Y. Choi, M. D
     .rees, J. Ellis, C.E.M.Wagner'
      WRITE(NOUT,11) 'arXiv:1208.2212 [hep-ph]'
      WRITE(NOUT,11) 'Comput.Phys.Commun.180:312-331,2009, arXiv:0712.23
     .60 [hep-ph]'
      WRITE(NOUT,11) 'Comput.Phys.Commun.156:283-317,2004, e-Print: hep-
     .ph/0307377'
      WRITE(NOUT,11) 'In case of problem(s): send email to jslee@jnu.ac.
     .kr'
      WRITE(NOUT,11) '-------------------------------------------------
     .----------------------------'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'SPINFO','Spectrum Program INFOrmation'
      WRITE(NOUT,30) 1,'CPsuperH    # Spectrum Calculator'
      WRITE(NOUT,30) 2,'2.3         # Version Number'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'DCINFO','DeCay program INFOrmation'
      WRITE(NOUT,30) 1,'CPsuperH    # Decay Package'
      WRITE(NOUT,30) 2,'2.3         # Version Number'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'MODSEL','MODel SELection'
      WRITE(NOUT,40) 5,2,'CP violation: Completely general CP phases'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
*Quark masses at arbitrary scale RSCALE:
*
      MTPOLE=MTPOLE_H
      MBPOLE=RAUX_H(1)
*       print*,'MBPOLE',mbpole
      MCPOLE=RAUX_H(4)
*    
      AS_MT =ASMT_H
      AS_MZ =ASMZ_H
      AS_MB =RAUX_H(3)
      AS_MC =RAUX_H(6)
*   
      MB_MBPOLE=RAUX_H(2)
      MC_MCPOLE=RAUX_H(5)
*   
      RSCALE=2.D0
      CALL MQ_RUN(RSCALE,MTPOLE,MZ_H,MBPOLE,MCPOLE
     .           ,MTMT_H,MBMT_H,MCMT_H,MSMT_H,MUMT_H,MDMT_H
     .           ,AS_MT,AS_MZ,AS_MB,AS_MC
     .           ,MT_RSCALE,MB_RSCALE,MC_RSCALE
     .                               ,MS_RSCALE,MU_RSCALE,MD_RSCALE)
*      print*,MT_RSCALE,MB_RSCALE,MC_RSCALE,MS_RSCALE,MU_RSCALE,MD_RSCALE
*
      WRITE(NOUT,20) 'SMINPUTS','Standard Model INPUTS'
      WRITE(NOUT,50) 1,SMPARA_H(1),'alpha_em^-1(M_Z)^MSbar'
      WRITE(NOUT,50) 2,GF_H,'G_F [GeV^-2]'
      WRITE(NOUT,50) 3,ASMZ_H,'alpha_S(M_Z)^MSbar'
      WRITE(NOUT,50) 4,MZ_H,'M_Z pole mass'
      WRITE(NOUT,50) 5,MB_MBPOLE,'mb(mb)^MSbar'
      WRITE(NOUT,50) 6,MTPOLE_H,'mt pole mass'
      WRITE(NOUT,50) 7,MTAU_H,'mtau pole mass'
      WRITE(NOUT,50) 11,ME_H,'me pole mass'
      WRITE(NOUT,50) 13,MMU_H,'mmu pole mass'
      WRITE(NOUT,50) 21,MD_RSCALE,'md(2 GeV)^MSbar'
      WRITE(NOUT,50) 22,MU_RSCALE,'mu(2 GeV)^MSbar'
      WRITE(NOUT,50) 23,MS_RSCALE,'ms(2 GeV)^MSbar'
      WRITE(NOUT,50) 24,MC_MCPOLE,'mc(mc)^MSbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'VCKMIN','CKM Matrix'
      WRITE(NOUT,50) 1,SMPARA_H(16),'lambda(M_Z)^MSbar [assumed]'
      WRITE(NOUT,50) 2,SMPARA_H(17),'A(M_Z)^MSbar [assumed]'
      WRITE(NOUT,50) 3,SMPARA_H(18),'rho^bar(M_Z)^MSbar [assumed]'
      WRITE(NOUT,50) 4,SMPARA_H(19),'eta^bar(M_Z)^MSbar [assumed]'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
*      WRITE(NOUT,20) 'MINPAR','MINimal model input PARameters'
*      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      S_INPUT=DSQRT(RAUX_H(13)) ! max(Q_stop,Q_sbottom)
      WRITE(NOUT,20) 'EXTPAR','non-minimal(EXT) model input PARameters [
     .Real Part]'
      WRITE(NOUT,50) 0,S_INPUT,'M_input: Input scale for EXTPAR'
      WRITE(NOUT,50) 1,DREAL(M1_H),'U(1)_Y gaugino (Bino) mass'
      WRITE(NOUT,50) 2,DREAL(M2_H),'SU(2)_L gaugino (Wino) mass'
      WRITE(NOUT,50) 3,DREAL(M3_H),'SU(3)_C gaugino (Gluino) mass'
      WRITE(NOUT,50) 11,DREAL(AT_H),'Top trilinear coupling'
      WRITE(NOUT,50) 12,DREAL(AB_H),'Bottom trilinear coupling'
      WRITE(NOUT,50) 13,DREAL(ATAU_H),'Tau trilinear coupling'
      WRITE(NOUT,50) 23,DREAL(MU_H),'mu parameter'
      WRITE(NOUT,50) 25,RAUX_H(21)/RAUX_H(17),'tan_beta(M_input)'
      WRITE(NOUT,50) 27,SSPARA_H(2),'Charged Higgs pole mass'
      WRITE(NOUT,50) 31,ML3_H*SSPARA_H(25),'Left 1st-gen. scalar lepton 
     .mass'
      WRITE(NOUT,50) 32,ML3_H*SSPARA_H(25),'Left 2nd-gen. scalar lepton 
     .mass'
      WRITE(NOUT,50) 33,ML3_H,'Left 3rd-gen. scalar lepton mass'
      WRITE(NOUT,50) 34,ME3_H*SSPARA_H(26),'Right scalar electron mass'
      WRITE(NOUT,50) 35,ME3_H*SSPARA_H(26),'Right scalar muon mass'
      WRITE(NOUT,50) 36,ME3_H,'Right scalar tau mass'
      WRITE(NOUT,50) 41,MQ3_H*SSPARA_H(22),'Left 1st-gen. scalar quark m
     .ass'
      WRITE(NOUT,50) 42,MQ3_H*SSPARA_H(22),'Left 2nd-gen. scalar quark m
     .ass'
      WRITE(NOUT,50) 43,MQ3_H,'Left 3rd-gen. scalar quark mass'
      WRITE(NOUT,50) 44,MU3_H*SSPARA_H(23),'Right scalar up mass'
      WRITE(NOUT,50) 45,MU3_H*SSPARA_H(23),'Right scalar charm mass'
      WRITE(NOUT,50) 46,MU3_H,'Right scalar top mass'
      WRITE(NOUT,50) 47,MD3_H*SSPARA_H(24),'Right scalar down mass'
      WRITE(NOUT,50) 48,MD3_H*SSPARA_H(24),'Right scalar strange mass'
      WRITE(NOUT,50) 49,MD3_H,'Right scalar bottom mass'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMEXTPAR','non-minimal(EXT) model input PARameters
     . [Imaginary Part]'
      WRITE(NOUT,50) 1,DIMAG(M1_H),'U(1)_Y gaugino (Bino) mass'
      WRITE(NOUT,50) 2,DIMAG(M2_H),'SU(2)_L gaugino (Wino) mass'
      WRITE(NOUT,50) 3,DIMAG(M3_H),'SU(3)_C gaugino (Gluino) mass'
      WRITE(NOUT,50) 11,DIMAG(AT_H),'Top trilinear coupling'
      WRITE(NOUT,50) 12,DIMAG(AB_H),'Bottom trilinear coupling'
      WRITE(NOUT,50) 13,DIMAG(ATAU_H),'Tau trilinear coupling'
      WRITE(NOUT,50) 23,DIMAG(MU_H),'mu parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'MASS','Mass Spectrum [Pole masses]'
      write(NOUT,11) 'PDG code           Mass       Particle'
      WRITE(NOUT,50) 24,MW_H,'MW'
      WRITE(NOUT,50) 25,HMASS_H(1),'H1'
      WRITE(NOUT,50) 35,HMASS_H(2),'H2'
      WRITE(NOUT,50) 36,HMASS_H(3),'H3'
      WRITE(NOUT,50) 37,SSPARA_H(2),'H+'
      WRITE(NOUT,50) 1000005,SBMASS_H(1),'~b_1'
      WRITE(NOUT,50) 2000005,SBMASS_H(2),'~b_2'
      WRITE(NOUT,50) 1000006,STMASS_H(1),'~t_1'
      WRITE(NOUT,50) 2000006,STMASS_H(2),'~t_2'
      WRITE(NOUT,50) 1000015,STAUMASS_H(1),'~tau_1'
      WRITE(NOUT,50) 2000015,STAUMASS_H(2),'~tau_2'
      WRITE(NOUT,50) 1000016,SNU3MASS_H,'~nu_tauL'
      WRITE(NOUT,50) 1000022,MN_H(1),'~chi_10'
      WRITE(NOUT,50) 1000023,MN_H(2),'~chi_20'
      WRITE(NOUT,50) 1000025,MN_H(3),'~chi_30'
      WRITE(NOUT,50) 1000035,MN_H(4),'~chi_40'
      WRITE(NOUT,50) 1000024,MC_H(1),'~chi_1+'
      WRITE(NOUT,50) 1000037,MC_H(2),'~chi_2+'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'HCOUPLINGS','Couplings of Higgs self interaction [
     .Real]'
      WRITE(NOUT,50) 1,DREAL(CAUX_H(201)),' RE lambda_1'
      WRITE(NOUT,50) 2,DREAL(CAUX_H(202)),' RE lambda_2'
      WRITE(NOUT,50) 3,DREAL(CAUX_H(203)),' RE lambda_3'
      WRITE(NOUT,50) 4,DREAL(CAUX_H(204)),' RE lambda_4'
      WRITE(NOUT,50) 5,DREAL(CAUX_H(205)),' RE lambda_5'
      WRITE(NOUT,50) 6,DREAL(CAUX_H(206)),' RE lambda_6'
      WRITE(NOUT,50) 7,DREAL(CAUX_H(207)),' RE lambda_7'
      WRITE(NOUT,20) 'IMHCOUPLINGS','Couplings of Higgs self interaction
     . [Imag.]'
      WRITE(NOUT,50) 1,DIMAG(CAUX_H(201)),' IM lambda_1'
      WRITE(NOUT,50) 2,DIMAG(CAUX_H(202)),' IM lambda_2'
      WRITE(NOUT,50) 3,DIMAG(CAUX_H(203)),' IM lambda_3'
      WRITE(NOUT,50) 4,DIMAG(CAUX_H(204)),' IM lambda_4'
      WRITE(NOUT,50) 5,DIMAG(CAUX_H(205)),' IM lambda_5'
      WRITE(NOUT,50) 6,DIMAG(CAUX_H(206)),' IM lambda_6'
      WRITE(NOUT,50) 7,DIMAG(CAUX_H(207)),' IM lambda_7'
      WRITE(NOUT,10) 
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'THRESHOLD','Td3,Tu3  threshold corrections'
*      Td3=CKB_SUM*TB_H/(1+CKB_SUM*TB_H)
*      Tu3=CKT_H/TB_H/(1.D0+CKT_H/TB_H)
*      print*,CAUX_H(10),' =? ',CAUX_H(211)/(1.D0-CAUX_H(211))/TB_H
      WRITE(NOUT,60) 5,1,DREAL(CAUX_H(211)),'Re(Td3)'
      WRITE(NOUT,60) 5,2,DIMAG(CAUX_H(211)),'Im(Td3)'
      WRITE(NOUT,60) 6,1,DREAL(CAUX_H(212)),'Re(Tu3)'
      WRITE(NOUT,60) 6,2,DIMAG(CAUX_H(212)),'Im(Tu3)'
      WRITE(NOUT,10) 
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'NMIX','Neutralino Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(N_H(1,1)),'Real[N_11] parameter'
      WRITE(NOUT,60) 1,2,DREAL(N_H(1,2)),'Real[N_12] parameter'
      WRITE(NOUT,60) 1,3,DREAL(N_H(1,3)),'Real[N_13] parameter'
      WRITE(NOUT,60) 1,4,DREAL(N_H(1,4)),'Real[N_14] parameter'
      WRITE(NOUT,60) 2,1,DREAL(N_H(2,1)),'Real[N_21] parameter'
      WRITE(NOUT,60) 2,2,DREAL(N_H(2,2)),'Real[N_22] parameter'
      WRITE(NOUT,60) 2,3,DREAL(N_H(2,3)),'Real[N_23] parameter'
      WRITE(NOUT,60) 2,4,DREAL(N_H(2,4)),'Real[N_24] parameter'
      WRITE(NOUT,60) 3,1,DREAL(N_H(3,1)),'Real[N_31] parameter'
      WRITE(NOUT,60) 3,2,DREAL(N_H(3,2)),'Real[N_32] parameter'
      WRITE(NOUT,60) 3,3,DREAL(N_H(3,3)),'Real[N_33] parameter'
      WRITE(NOUT,60) 3,4,DREAL(N_H(3,4)),'Real[N_34] parameter'
      WRITE(NOUT,60) 4,1,DREAL(N_H(4,1)),'Real[N_41] parameter'
      WRITE(NOUT,60) 4,2,DREAL(N_H(4,2)),'Real[N_42] parameter'
      WRITE(NOUT,60) 4,3,DREAL(N_H(4,3)),'Real[N_43] parameter'
      WRITE(NOUT,60) 4,4,DREAL(N_H(4,4)),'Real[N_44] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMNMIX','Neutralino Mixing Matrix [Imaginary Part]
     .'
      WRITE(NOUT,60) 1,1,DIMAG(N_H(1,1)),'Imag[N_11] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(N_H(1,2)),'Imag[N_12] parameter'
      WRITE(NOUT,60) 1,3,DIMAG(N_H(1,3)),'Imag[N_13] parameter'
      WRITE(NOUT,60) 1,4,DIMAG(N_H(1,4)),'Imag[N_14] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(N_H(2,1)),'Imag[N_21] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(N_H(2,2)),'Imag[N_22] parameter'
      WRITE(NOUT,60) 2,3,DIMAG(N_H(2,3)),'Imag[N_23] parameter'
      WRITE(NOUT,60) 2,4,DIMAG(N_H(2,4)),'Imag[N_24] parameter'
      WRITE(NOUT,60) 3,1,DIMAG(N_H(3,1)),'Imag[N_31] parameter'
      WRITE(NOUT,60) 3,2,DIMAG(N_H(3,2)),'Imag[N_32] parameter'
      WRITE(NOUT,60) 3,3,DIMAG(N_H(3,3)),'Imag[N_33] parameter'
      WRITE(NOUT,60) 3,4,DIMAG(N_H(3,4)),'Imag[N_34] parameter'
      WRITE(NOUT,60) 4,1,DIMAG(N_H(4,1)),'Imag[N_41] parameter'
      WRITE(NOUT,60) 4,2,DIMAG(N_H(4,2)),'Imag[N_42] parameter'
      WRITE(NOUT,60) 4,3,DIMAG(N_H(4,3)),'Imag[N_43] parameter'
      WRITE(NOUT,60) 4,4,DIMAG(N_H(4,4)),'Imag[N_44] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'UMIX','Chargino U Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(UL_H(1,1)),'Real[U_11] parameter'
      WRITE(NOUT,60) 1,2,DREAL(UL_H(1,2)),'Real[U_12] parameter'
      WRITE(NOUT,60) 2,1,DREAL(UL_H(2,1)),'Real[U_21] parameter'
      WRITE(NOUT,60) 2,2,DREAL(UL_H(2,2)),'Real[U_22] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMUMIX','Chargino U Mixing Matrix [Imaginary Part]
     .'
      WRITE(NOUT,60) 1,1,DIMAG(UL_H(1,1)),'Imag[U_11] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(UL_H(1,2)),'Imag[U_12] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(UL_H(2,1)),'Imag[U_21] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(UL_H(2,2)),'Imag[U_22] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      VMIX(1,1)=DCONJG(UR_H(1,1))  ! V = UR^*
      VMIX(1,2)=DCONJG(UR_H(1,2))
      VMIX(2,1)=DCONJG(UR_H(2,1))
      VMIX(2,2)=DCONJG(UR_H(2,2))
      WRITE(NOUT,20) 'VMIX','Chargino V Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(VMIX(1,1)),'Real[V_11] parameter'
      WRITE(NOUT,60) 1,2,DREAL(VMIX(1,2)),'Real[V_12] parameter'
      WRITE(NOUT,60) 2,1,DREAL(VMIX(2,1)),'Real[V_21] parameter'
      WRITE(NOUT,60) 2,2,DREAL(VMIX(2,2)),'Real[V_22] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMVMIX','Chargino V Mixing Matrix [Imaginary Part]
     .'
      WRITE(NOUT,60) 1,1,DIMAG(VMIX(1,1)),'Imag[V_11] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(VMIX(1,2)),'Imag[V_12] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(VMIX(2,1)),'Imag[V_21] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(VMIX(2,2)),'Imag[V_22] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      FMIX(1,1)=DCONJG(STMIX_H(1,1))  ! F = U~f^dagger
      FMIX(1,2)=DCONJG(STMIX_H(2,1))
      FMIX(2,1)=DCONJG(STMIX_H(1,2))
      FMIX(2,2)=DCONJG(STMIX_H(2,2))
      WRITE(NOUT,20) 'STOPMIX','Stop Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(FMIX(1,1)),'Real[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DREAL(FMIX(1,2)),'Real[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DREAL(FMIX(2,1)),'Real[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DREAL(FMIX(2,2)),'Real[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMSTOPMIX','Stop Mixing Matrix [Imaginary Part]'
      WRITE(NOUT,60) 1,1,DIMAG(FMIX(1,1)),'Imag[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(FMIX(1,2)),'Imag[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(FMIX(2,1)),'Imag[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(FMIX(2,2)),'Imag[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      FMIX(1,1)=DCONJG(SBMIX_H(1,1))  ! F = U~f^dagger
      FMIX(1,2)=DCONJG(SBMIX_H(2,1))
      FMIX(2,1)=DCONJG(SBMIX_H(1,2))
      FMIX(2,2)=DCONJG(SBMIX_H(2,2))
      WRITE(NOUT,20) 'SBOTMIX','Sbottom Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(FMIX(1,1)),'Real[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DREAL(FMIX(1,2)),'Real[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DREAL(FMIX(2,1)),'Real[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DREAL(FMIX(2,2)),'Real[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMSBOTMIX','Sbottom Mixing Matrix [Imaginary Part]
     .'
      WRITE(NOUT,60) 1,1,DIMAG(FMIX(1,1)),'Imag[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(FMIX(1,2)),'Imag[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(FMIX(2,1)),'Imag[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(FMIX(2,2)),'Imag[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      FMIX(1,1)=DCONJG(STAUMIX_H(1,1))  ! F = U~f^dagger
      FMIX(1,2)=DCONJG(STAUMIX_H(2,1))
      FMIX(2,1)=DCONJG(STAUMIX_H(1,2))
      FMIX(2,2)=DCONJG(STAUMIX_H(2,2))
      WRITE(NOUT,20) 'STAUMIX','Stau Mixing Matrix [Real Part]'
      WRITE(NOUT,60) 1,1,DREAL(FMIX(1,1)),'Real[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DREAL(FMIX(1,2)),'Real[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DREAL(FMIX(2,1)),'Real[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DREAL(FMIX(2,2)),'Real[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,20) 'IMSTAUMIX','Stau Mixing Matrix [Imaginary Part]'
      WRITE(NOUT,60) 1,1,DIMAG(FMIX(1,1)),'Imag[F_1L] parameter'
      WRITE(NOUT,60) 1,2,DIMAG(FMIX(1,2)),'Imag[F_1R] parameter'
      WRITE(NOUT,60) 2,1,DIMAG(FMIX(2,1)),'Imag[F_2L] parameter'
      WRITE(NOUT,60) 2,2,DIMAG(FMIX(2,2)),'Imag[F_2R] parameter'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      DO I=1,3
       DO J=1,4
        IF(J.EQ.4) THEN
         RTMP34(I,J)=0.D0
        ELSE
         RTMP34(I,J)=OMIX_H(J,I) ! Transpose
        ENDIF
       ENDDO
*       print*,(rtmp34(i,j),j=1,4)
      ENDDO
      DO I=1,4
       DO J=1,4
        RTMP44(I,J)=0.D0
        IF(I.EQ.1 .AND. J.EQ.1) RTMP44(I,J)=1.D0
        IF(I.EQ.2 .AND. J.EQ.2) RTMP44(I,J)=1.D0
        IF(I.EQ.3 .AND. J.EQ.3) RTMP44(I,J)=-SB_H
        IF(I.EQ.3 .AND. J.EQ.4) RTMP44(I,J)= CB_H
        IF(I.EQ.4 .AND. J.EQ.3) RTMP44(I,J)= CB_H
        IF(I.EQ.4 .AND. J.EQ.4) RTMP44(I,J)= SB_H
       ENDDO
*       print*,(rtmp44(i,j),j=1,4)
      ENDDO
      DO I=1,3
       DO J=1,4
        SMIX(I,J)=RTMP34(I,1)*RTMP44(1,J)
     .           +RTMP34(I,2)*RTMP44(2,J)
     .           +RTMP34(I,3)*RTMP44(3,J)
     .           +RTMP34(I,4)*RTMP44(4,J)
       ENDDO
*       print*,(smix(i,j),j=1,4)
      ENDDO
      WRITE(NOUT,20) 'CVHMIX','3 x 4 Higgs Mixing Matrix'
      WRITE(NOUT,60) 1,1,SMIX(1,1),'CVHMIX_11'
      WRITE(NOUT,60) 1,2,SMIX(1,2),'CVHMIX_12'
      WRITE(NOUT,60) 1,3,SMIX(1,3),'CVHMIX_13'
      WRITE(NOUT,60) 1,4,SMIX(1,4),'CVHMIX_14'
      WRITE(NOUT,60) 2,1,SMIX(2,1),'CVHMIX_21'
      WRITE(NOUT,60) 2,2,SMIX(2,2),'CVHMIX_22'
      WRITE(NOUT,60) 2,3,SMIX(2,3),'CVHMIX_23'
      WRITE(NOUT,60) 2,4,SMIX(2,4),'CVHMIX_24'
      WRITE(NOUT,60) 3,1,SMIX(3,1),'CVHMIX_31'
      WRITE(NOUT,60) 3,2,SMIX(3,2),'CVHMIX_32'
      WRITE(NOUT,60) 3,3,SMIX(3,3),'CVHMIX_33'
      WRITE(NOUT,60) 3,4,SMIX(3,4),'CVHMIX_34'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'AU Q=',QSCALE,'Up-type quark trilinear coupling at
     . Q [Real]'
      WRITE(NOUT,60) 3,3,DREAL(AT_H),'At(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(33)*DCOS(SSPARA_H(34)/180.D0*PI),'Ac(Q
     .)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(31)*DCOS(SSPARA_H(32)/180.D0*PI),'Au(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'IMAU Q=',QSCALE,'Up-type quark trilinear coupling 
     .at Q [Imag.]'
      WRITE(NOUT,60) 3,3,DIMAG(AT_H),'At(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(33)*DSIN(SSPARA_H(34)/180.D0*PI),'Ac(Q
     .)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(31)*DSIN(SSPARA_H(32)/180.D0*PI),'Au(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'AD Q=',QSCALE,'Down-type quark trilinear coupling 
     .at Q [Real]'
      WRITE(NOUT,60) 3,3,DREAL(AB_H),'Ab(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(37)*DCOS(SSPARA_H(38)/180.D0*PI),'As(Q
     .)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(35)*DCOS(SSPARA_H(36)/180.D0*PI),'Ad(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'IMAD Q=',QSCALE,'Down-type quark trilinear couplin
     .g at Q [Imag.]'
      WRITE(NOUT,60) 3,3,DIMAG(AB_H),'Ab(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(37)*DSIN(SSPARA_H(38)/180.D0*PI),'As(Q
     .)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(35)*DSIN(SSPARA_H(36)/180.D0*PI),'Ad(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'AE Q=',QSCALE,'Charged lepton trilinear coupling a
     .t Q [Real]'
      WRITE(NOUT,60) 3,3,DREAL(ATAU_H),'Atau(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(29)*DCOS(SSPARA_H(30)/180.D0*PI),'Amu(
     .Q)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(27)*DCOS(SSPARA_H(28)/180.D0*PI),'Ae(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=S_INPUT
      WRITE(NOUT,21) 'IMAE Q=',QSCALE,'Charged lepton trilinear coupling
     . at Q [Imag.]'
      WRITE(NOUT,60) 3,3,DIMAG(ATAU_H),'Atau(Q)MSSM DRbar'
      WRITE(NOUT,60) 2,2,SSPARA_H(29)*DSIN(SSPARA_H(30)/180.D0*PI),'Amu(
     .Q)MSSM DRbar'
      WRITE(NOUT,60) 1,1,SSPARA_H(27)*DSIN(SSPARA_H(28)/180.D0*PI),'Ae(Q
     .)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'YU Q=',QSCALE,'Up-type quark Yukawa coupling at Q 
     .[Real]'
      WRITE(NOUT,60) 3,3,RAUX_H(24)*DREAL(CAUX_H(1)),'y_t(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'IMYU Q=',QSCALE,'Up-type quark Yukawa coupling at 
     .Q [Imag.]'
      WRITE(NOUT,60) 3,3,RAUX_H(24)*DIMAG(CAUX_H(1)),'y_t(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'YD Q=',QSCALE,'Down-type quark Yukawa coupling at 
     .Q [Real]'
      WRITE(NOUT,60) 3,3,RAUX_H(27)*DREAL(CAUX_H(2)),'y_b(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'IMYD Q=',QSCALE,'Down-type quark Yukawa coupling a
     .t Q [Imag.]'
      WRITE(NOUT,60) 3,3,RAUX_H(27)*DIMAG(CAUX_H(2)),'y_b(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'YE Q=',QSCALE,'Charged lepton Yukawa coupling at Q
     . [Real]'
      WRITE(NOUT,60) 3,3,DREAL(CAUX_H(13)),'y_tau(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      QSCALE=MTPOLE_H
      WRITE(NOUT,21) 'IMYE Q=',QSCALE,'Charged lepton Yukawa coupling at
     . Q [Imag.]'
      WRITE(NOUT,60) 3,3,DIMAG(CAUX_H(13)),'y_tau(Q)MSSM DRbar'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
*Mid-Line
      WRITE(NOUT,11) '--------------------------------------------------
     .----------------------------'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      IH=1
      IPDG=25
      WRITE(NOUT,11) '     PDG code          Width'
      WRITE(NOUT,100) IPDG,GAMBRN(101,1,IH),'H1 Decays'
      WRITE(NOUT,11) '         BR         NDA        ID1       ID2'
      IF(GAMBRN(1,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(1,3,IH),2,11,-11,'H1 -> e- e+'
      IF(GAMBRN(2,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(2,3,IH),2,13,-13,'H1 -> mu- mu+'
      IF(GAMBRN(3,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(3,3,IH),2,15,-15,'H1 -> tau- tau+'
      IF(GAMBRN(4,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(4,3,IH),2,1,-1,'H1 -> d dbar'
      IF(GAMBRN(5,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(5,3,IH),2,3,-3,'H1 -> s sbar'
      IF(GAMBRN(6,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(6,3,IH),2,5,-5,'H1 -> b bbar'
      IF(GAMBRN(7,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(7,3,IH),2,2,-2,'H1 -> u ubar'
      IF(GAMBRN(8,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(8,3,IH),2,4,-4,'H1 -> c cbar'
      IF(GAMBRN(9,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(9,3,IH),2,6,-6,'H1 -> t tbar'
      IF(GAMBRN(10,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(10,3,IH),2,-24,24,'H1 -> W- W+'
      IF(GAMBRN(11,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(11,3,IH),2,23,23,'H1 -> Z Z'
      IF(GAMBRN(12,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(12,3,IH),2,25,23,'H1 -> H1 Z'
      IF(GAMBRN(13,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(13,3,IH),2,35,23,'H1 -> H2 Z'
      IF(GAMBRN(14,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(14,3,IH),2,25,25,'H1 -> H1 H1'
      IF(GAMBRN(15,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(15,3,IH),2,25,35,'H1 -> H1 H2'
      IF(GAMBRN(16,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(16,3,IH),2,35,35,'H1 -> H2 H2'
      IF(GAMBRN(17,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(17,3,IH),2,22,22,'H1 -> gam gam'
      IF(GAMBRN(18,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(18,3,IH),2,21,21,'H1 -> g g'
      IF(GAMBRN(19,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(19,3,IH),2,23,22,'H1 -> Z gam'
      IF(GAMBRN(51,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(51,3,IH),2,1000022,1000022,'H1 -> ~chi_10 ~
     .chi_10'
      IF(GAMBRN(52,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(52,3,IH),2,1000022,1000023,'H1 -> ~chi_10 ~
     .chi_20'
      IF(GAMBRN(53,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(53,3,IH),2,1000022,1000025,'H1 -> ~chi_10 ~
     .chi_30'
      IF(GAMBRN(54,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(54,3,IH),2,1000022,1000035,'H1 -> ~chi_10 ~
     .chi_40'
      IF(GAMBRN(55,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(55,3,IH),2,1000023,1000023,'H1 -> ~chi_20 ~
     .chi_20'
      IF(GAMBRN(56,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(56,3,IH),2,1000023,1000025,'H1 -> ~chi_20 ~
     .chi_30'
      IF(GAMBRN(57,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(57,3,IH),2,1000023,1000035,'H1 -> ~chi_20 ~
     .chi_40'
      IF(GAMBRN(58,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(58,3,IH),2,1000025,1000025,'H1 -> ~chi_30 ~
     .chi_30'
      IF(GAMBRN(59,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(59,3,IH),2,1000025,1000035,'H1 -> ~chi_30 ~
     .chi_40'
      IF(GAMBRN(60,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(60,3,IH),2,1000035,1000035,'H1 -> ~chi_40 ~
     .chi_40'
      IF(GAMBRN(61,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(61,3,IH),2,1000024,-1000024,'H1 -> ~chi_1+ 
     .~chi_1-'
      IF(GAMBRN(62,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(62,3,IH),2,1000024,-1000037,'H1 -> ~chi_1+ 
     .~chi_2-'
      IF(GAMBRN(63,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(63,3,IH),2,1000037,-1000024,'H1 -> ~chi_2+ 
     .~chi_1-'
      IF(GAMBRN(64,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(64,3,IH),2,1000037,-1000037,'H1 -> ~chi_2+ 
     .~chi_2-'
      IF(GAMBRN(65,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(65,3,IH),2,1000006,-1000006,'H1 -> ~t_1 ~t_
     .1*'
      IF(GAMBRN(66,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(66,3,IH),2,2000006,-1000006,'H1 -> ~t_2 ~t_
     .1*'
      IF(GAMBRN(67,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(67,3,IH),2,1000006,-2000006,'H1 -> ~t_1 ~t_
     .2*'
      IF(GAMBRN(68,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(68,3,IH),2,2000006,-2000006,'H1 -> ~t_2 ~t_
     .2*'
      IF(GAMBRN(69,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(69,3,IH),2,1000005,-1000005,'H1 -> ~b_1 ~b_
     .1*'
      IF(GAMBRN(70,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(70,3,IH),2,2000005,-1000005,'H1 -> ~b_2 ~b_
     .1*'
      IF(GAMBRN(71,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(71,3,IH),2,1000005,-2000005,'H1 -> ~b_1 ~b_
     .2*'
      IF(GAMBRN(72,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(72,3,IH),2,2000005,-2000005,'H1 -> ~b_2 ~b_
     .2*'
      IF(GAMBRN(73,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(73,3,IH),2,1000015,-1000015,'H1 -> ~tau_1- 
     .~tau_1+'
      IF(GAMBRN(74,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(74,3,IH),2,2000015,-1000015,'H1 -> ~tau_2- 
     .~tau_1+'
      IF(GAMBRN(75,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(75,3,IH),2,1000015,-2000015,'H1 -> ~tau_1- 
     .~tau_2+'
      IF(GAMBRN(76,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(76,3,IH),2,2000015,-2000015,'H1 -> ~tau_2- 
     .~tau_2+'
      IF(GAMBRN(77,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(77,3,IH),2,1000016,-1000016,'H1 -> ~nu_tauL
     . ~nu_tauL*'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      IH=2
      IPDG=35
      WRITE(NOUT,11) '     PDG code          Width'
      WRITE(NOUT,100) IPDG,GAMBRN(101,1,IH),'H2 Decays'
      WRITE(NOUT,11) '         BR         NDA        ID1       ID2'
      IF(GAMBRN(1,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(1,3,IH),2,11,-11,'H2 -> e- e+'
      IF(GAMBRN(2,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(2,3,IH),2,13,-13,'H2 -> mu- mu+'
      IF(GAMBRN(3,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(3,3,IH),2,15,-15,'H2 -> tau- tau+'
      IF(GAMBRN(4,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(4,3,IH),2,1,-1,'H2 -> d dbar'
      IF(GAMBRN(5,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(5,3,IH),2,3,-3,'H2 -> s sbar'
      IF(GAMBRN(6,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(6,3,IH),2,5,-5,'H2 -> b bbar'
      IF(GAMBRN(7,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(7,3,IH),2,2,-2,'H2 -> u ubar'
      IF(GAMBRN(8,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(8,3,IH),2,4,-4,'H2 -> c cbar'
      IF(GAMBRN(9,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(9,3,IH),2,6,-6,'H2 -> t tbar'
      IF(GAMBRN(10,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(10,3,IH),2,-24,24,'H2 -> W- W+'
      IF(GAMBRN(11,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(11,3,IH),2,23,23,'H2 -> Z Z'
      IF(GAMBRN(12,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(12,3,IH),2,25,23,'H2 -> H1 Z'
      IF(GAMBRN(13,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(13,3,IH),2,35,23,'H2 -> H2 Z'
      IF(GAMBRN(14,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(14,3,IH),2,25,25,'H2 -> H1 H1'
      IF(GAMBRN(15,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(15,3,IH),2,25,35,'H2 -> H1 H2'
      IF(GAMBRN(16,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(16,3,IH),2,35,35,'H2 -> H2 H2'
      IF(GAMBRN(17,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(17,3,IH),2,22,22,'H2 -> gam gam'
      IF(GAMBRN(18,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(18,3,IH),2,21,21,'H2 -> g g'
      IF(GAMBRN(19,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(19,3,IH),2,23,22,'H2 -> Z gam'
      IF(GAMBRN(51,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(51,3,IH),2,1000022,1000022,'H2 -> ~chi_10 ~
     .chi_10'
      IF(GAMBRN(52,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(52,3,IH),2,1000022,1000023,'H2 -> ~chi_10 ~
     .chi_20'
      IF(GAMBRN(53,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(53,3,IH),2,1000022,1000025,'H2 -> ~chi_10 ~
     .chi_30'
      IF(GAMBRN(54,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(54,3,IH),2,1000022,1000035,'H2 -> ~chi_10 ~
     .chi_40'
      IF(GAMBRN(55,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(55,3,IH),2,1000023,1000023,'H2 -> ~chi_20 ~
     .chi_20'
      IF(GAMBRN(56,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(56,3,IH),2,1000023,1000025,'H2 -> ~chi_20 ~
     .chi_30'
      IF(GAMBRN(57,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(57,3,IH),2,1000023,1000035,'H2 -> ~chi_20 ~
     .chi_40'
      IF(GAMBRN(58,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(58,3,IH),2,1000025,1000025,'H2 -> ~chi_30 ~
     .chi_30'
      IF(GAMBRN(59,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(59,3,IH),2,1000025,1000035,'H2 -> ~chi_30 ~
     .chi_40'
      IF(GAMBRN(60,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(60,3,IH),2,1000035,1000035,'H2 -> ~chi_40 ~
     .chi_40'
      IF(GAMBRN(61,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(61,3,IH),2,1000024,-1000024,'H2 -> ~chi_1+ 
     .~chi_1-'
      IF(GAMBRN(62,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(62,3,IH),2,1000024,-1000037,'H2 -> ~chi_1+ 
     .~chi_2-'
      IF(GAMBRN(63,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(63,3,IH),2,1000037,-1000024,'H2 -> ~chi_2+ 
     .~chi_1-'
      IF(GAMBRN(64,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(64,3,IH),2,1000037,-1000037,'H2 -> ~chi_2+ 
     .~chi_2-'
      IF(GAMBRN(65,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(65,3,IH),2,1000006,-1000006,'H2 -> ~t_1 ~t_
     .1*'
      IF(GAMBRN(66,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(66,3,IH),2,2000006,-1000006,'H2 -> ~t_2 ~t_
     .1*'
      IF(GAMBRN(67,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(67,3,IH),2,1000006,-2000006,'H2 -> ~t_1 ~t_
     .2*'
      IF(GAMBRN(68,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(68,3,IH),2,2000006,-2000006,'H2 -> ~t_2 ~t_
     .2*'
      IF(GAMBRN(69,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(69,3,IH),2,1000005,-1000005,'H2 -> ~b_1 ~b_
     .1*'
      IF(GAMBRN(70,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(70,3,IH),2,2000005,-1000005,'H2 -> ~b_2 ~b_
     .1*'
      IF(GAMBRN(71,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(71,3,IH),2,1000005,-2000005,'H2 -> ~b_1 ~b_
     .2*'
      IF(GAMBRN(72,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(72,3,IH),2,2000005,-2000005,'H2 -> ~b_2 ~b_
     .2*'
      IF(GAMBRN(73,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(73,3,IH),2,1000015,-1000015,'H2 -> ~tau_1- 
     .~tau_1+'
      IF(GAMBRN(74,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(74,3,IH),2,2000015,-1000015,'H2 -> ~tau_2- 
     .~tau_1+'
      IF(GAMBRN(75,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(75,3,IH),2,1000015,-2000015,'H2 -> ~tau_1- 
     .~tau_2+'
      IF(GAMBRN(76,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(76,3,IH),2,2000015,-2000015,'H2 -> ~tau_2- 
     .~tau_2+'
      IF(GAMBRN(77,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(77,3,IH),2,1000016,-1000016,'H2 -> ~nu_tauL
     . ~nu_tauL*'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      IH=3
      IPDG=36
      WRITE(NOUT,11) '     PDG code          Width'
      WRITE(NOUT,100) IPDG,GAMBRN(101,1,IH),'H3 Decays'
      WRITE(NOUT,11) '         BR         NDA        ID1       ID2'
      IF(GAMBRN(1,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(1,3,IH),2,11,-11,'H3 -> e- e+'
      IF(GAMBRN(2,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(2,3,IH),2,13,-13,'H3 -> mu- mu+'
      IF(GAMBRN(3,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(3,3,IH),2,15,-15,'H3 -> tau- tau+'
      IF(GAMBRN(4,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(4,3,IH),2,1,-1,'H3 -> d dbar'
      IF(GAMBRN(5,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(5,3,IH),2,3,-3,'H3 -> s sbar'
      IF(GAMBRN(6,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(6,3,IH),2,5,-5,'H3 -> b bbar'
      IF(GAMBRN(7,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(7,3,IH),2,2,-2,'H3 -> u ubar'
      IF(GAMBRN(8,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(8,3,IH),2,4,-4,'H3 -> c cbar'
      IF(GAMBRN(9,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(9,3,IH),2,6,-6,'H3 -> t tbar'
      IF(GAMBRN(10,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(10,3,IH),2,-24,24,'H3 -> W- W+'
      IF(GAMBRN(11,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(11,3,IH),2,23,23,'H3 -> Z Z'
      IF(GAMBRN(12,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(12,3,IH),2,25,23,'H3 -> H1 Z'
      IF(GAMBRN(13,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(13,3,IH),2,35,23,'H3 -> H2 Z'
      IF(GAMBRN(14,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(14,3,IH),2,25,25,'H3 -> H1 H1'
      IF(GAMBRN(15,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(15,3,IH),2,25,35,'H3 -> H1 H2'
      IF(GAMBRN(16,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(16,3,IH),2,35,35,'H3 -> H2 H2'
      IF(GAMBRN(17,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(17,3,IH),2,22,22,'H3 -> gam gam'
      IF(GAMBRN(18,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(18,3,IH),2,21,21,'H3 -> g g'
      IF(GAMBRN(19,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(19,3,IH),2,23,22,'H3 -> Z gam'
      IF(GAMBRN(51,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(51,3,IH),2,1000022,1000022,'H3 -> ~chi_10 ~
     .chi_10'
      IF(GAMBRN(52,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(52,3,IH),2,1000022,1000023,'H3 -> ~chi_10 ~
     .chi_20'
      IF(GAMBRN(53,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(53,3,IH),2,1000022,1000025,'H3 -> ~chi_10 ~
     .chi_30'
      IF(GAMBRN(54,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(54,3,IH),2,1000022,1000035,'H3 -> ~chi_10 ~
     .chi_40'
      IF(GAMBRN(55,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(55,3,IH),2,1000023,1000023,'H3 -> ~chi_20 ~
     .chi_20'
      IF(GAMBRN(56,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(56,3,IH),2,1000023,1000025,'H3 -> ~chi_20 ~
     .chi_30'
      IF(GAMBRN(57,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(57,3,IH),2,1000023,1000035,'H3 -> ~chi_20 ~
     .chi_40'
      IF(GAMBRN(58,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(58,3,IH),2,1000025,1000025,'H3 -> ~chi_30 ~
     .chi_30'
      IF(GAMBRN(59,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(59,3,IH),2,1000025,1000035,'H3 -> ~chi_30 ~
     .chi_40'
      IF(GAMBRN(60,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(60,3,IH),2,1000035,1000035,'H3 -> ~chi_40 ~
     .chi_40'
      IF(GAMBRN(61,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(61,3,IH),2,1000024,-1000024,'H3 -> ~chi_1+ 
     .~chi_1-'
      IF(GAMBRN(62,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(62,3,IH),2,1000024,-1000037,'H3 -> ~chi_1+ 
     .~chi_2-'
      IF(GAMBRN(63,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(63,3,IH),2,1000037,-1000024,'H3 -> ~chi_2+ 
     .~chi_1-'
      IF(GAMBRN(64,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(64,3,IH),2,1000037,-1000037,'H3 -> ~chi_2+ 
     .~chi_2-'
      IF(GAMBRN(65,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(65,3,IH),2,1000006,-1000006,'H3 -> ~t_1 ~t_
     .1*'
      IF(GAMBRN(66,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(66,3,IH),2,2000006,-1000006,'H3 -> ~t_2 ~t_
     .1*'
      IF(GAMBRN(67,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(67,3,IH),2,1000006,-2000006,'H3 -> ~t_1 ~t_
     .2*'
      IF(GAMBRN(68,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(68,3,IH),2,2000006,-2000006,'H3 -> ~t_2 ~t_
     .2*'
      IF(GAMBRN(69,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(69,3,IH),2,1000005,-1000005,'H3 -> ~b_1 ~b_
     .1*'
      IF(GAMBRN(70,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(70,3,IH),2,2000005,-1000005,'H3 -> ~b_2 ~b_
     .1*'
      IF(GAMBRN(71,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(71,3,IH),2,1000005,-2000005,'H3 -> ~b_1 ~b_
     .2*'
      IF(GAMBRN(72,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(72,3,IH),2,2000005,-2000005,'H3 -> ~b_2 ~b_
     .2*'
      IF(GAMBRN(73,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(73,3,IH),2,1000015,-1000015,'H3 -> ~tau_1- 
     .~tau_1+'
      IF(GAMBRN(74,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(74,3,IH),2,2000015,-1000015,'H3 -> ~tau_2- 
     .~tau_1+'
      IF(GAMBRN(75,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(75,3,IH),2,1000015,-2000015,'H3 -> ~tau_1- 
     .~tau_2+'
      IF(GAMBRN(76,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(76,3,IH),2,2000015,-2000015,'H3 -> ~tau_2- 
     .~tau_2+'
      IF(GAMBRN(77,3,IH).GT.0.D0)
     .WRITE(NOUT,110) GAMBRN(77,3,IH),2,1000016,-1000016,'H3 -> ~nu_tauL
     . ~nu_tauL*'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      ISMC   = IFLAG_H(22)
      ISUSYC = IFLAG_H(23)
      WRITE(NOUT,11) '     PDG code          Width'
      WRITE(NOUT,100) 37,GAMBRC(ISMC+ISUSYC+1,1),'H+ Decays'
      WRITE(NOUT,11) '         BR         NDA        ID1       ID2'
      IF(GAMBRC(1,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(1,3),2,-11,12,'H+ -> e+ nu_e'
      IF(GAMBRC(2,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(2,3),2,-13,14,'H+ -> mu+ nu_mu'
      IF(GAMBRC(3,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(3,3),2,-15,16,'H+ -> tau+ nu_tau'
      IF(GAMBRC(4,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(4,3),2,2,-1,'H+ -> u dbar'
      IF(GAMBRC(5,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(5,3),2,4,-3,'H+ -> c sbar'
      IF(GAMBRC(6,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(6,3),2,6,-5,'H+ -> t bbar'
      IF(GAMBRC(7,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(7,3),2,24,25,'H+ -> W+ H1'
      IF(GAMBRC(8,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(8,3),2,24,35,'H+ -> W+ H2'
      IF(GAMBRC(9,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(9,3),2,24,36,'H+ -> W+ H3'
      IF(GAMBRC(ISMC+1,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+1,3),2,1000024,1000022,'H+ -> ~chi_1+ 
     .~chi_10'
      IF(GAMBRC(ISMC+2,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+2,3),2,1000024,1000023,'H+ -> ~chi_1+ 
     .~chi_20'
      IF(GAMBRC(ISMC+3,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+3,3),2,1000024,1000025,'H+ -> ~chi_1+ 
     .~chi_30'
      IF(GAMBRC(ISMC+4,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+4,3),2,1000024,1000035,'H+ -> ~chi_1+ 
     .~chi_40'
      IF(GAMBRC(ISMC+5,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+5,3),2,1000037,1000022,'H+ -> ~chi_2+ 
     .~chi_10'
      IF(GAMBRC(ISMC+6,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+6,3),2,1000037,1000023,'H+ -> ~chi_2+ 
     .~chi_20'
      IF(GAMBRC(ISMC+7,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+7,3),2,1000037,1000025,'H+ -> ~chi_2+ 
     .~chi_30'
      IF(GAMBRC(ISMC+8,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+8,3),2,1000037,1000035,'H+ -> ~chi_2+ 
     .~chi_40'
      IF(GAMBRC(ISMC+9,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+9,3),2,1000006,-1000005,'H+ -> ~t_1 ~b
     ._1*'
      IF(GAMBRC(ISMC+10,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+10,3),2,1000006,-2000005,'H+ -> ~t_1 ~
     .b_2*'
      IF(GAMBRC(ISMC+11,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+11,3),2,2000006,-1000005,'H+ -> ~t_2 ~
     .b_1*'
      IF(GAMBRC(ISMC+12,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+12,3),2,2000006,-2000005,'H+ -> ~t_2 ~
     .b_2*'
      IF(GAMBRC(ISMC+13,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+13,3),2,-1000015,1000016,'H+ -> ~tau_1
     .+ ~nu_tauL'
      IF(GAMBRC(ISMC+14,3).GT.0.D0)
     .WRITE(NOUT,110) GAMBRC(ISMC+14,3),2,-2000015,1000016,'H+ -> ~tau_2
     .+ ~nu_tauL'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,11) '     PDG code          Width'
      WRITE(NOUT,100) 6,RAUX_H(50)+RAUX_H(51),'Top Decays'
      WRITE(NOUT,11) '         BR         NDA        ID1       ID2'
      WRITE(NOUT,110) RAUX_H(52),2,24,5,'t -> W+ b' 
      WRITE(NOUT,110) RAUX_H(53),2,37,5,'t -> H+ b' 
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,FMT='("Block FOBS  # Flavour observables")')
      WRITE(NOUT,FMT='("# ParentPDG type   M     value            q  NDA
     .        ID1         ID2")')
      WRITE(NOUT,120)5,1,1,RAUX_H(135)*1.D-4,0,2,3,22,'BR(b->s gamma)'
      WRITE(NOUT,120)5,3,1,RAUX_H(136)/100.D0,0,2,3,22,'ACP[BR(b->s gamm
     .a)]'
      WRITE(NOUT,120)531,1,1,RAUX_H(130)*1.D-7,0,2,13,-13,'BR(B_s->mu+ m
     .u-)'
      WRITE(NOUT,120)521,2,1,RAUX_H(134),0,2,-15,16,'BR(B_u->tau nu)/BR(
     .SM)'
      WRITE(NOUT,120)511,1,1,RAUX_H(131)*1.D-7,0,2,-15,15,'BR(B_d->tau+ 
     .tau-)'
      WRITE(NOUT,121)511,7,1,RAUX_H(132)/2.D0,0,1,-511,'Delta M [B_d] (S
     .USY) in ps^-1'
      WRITE(NOUT,121)531,7,1,RAUX_H(133)/2.D0,0,1,-533,'Delta M [B_s] (S
     .USY) in ps^-1'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
      WRITE(NOUT,FMT='("Block FDIPOLE  # Electric and Magnetic dipole mo
     .ments")')
      WRITE(NOUT,FMT='("# PDG_code     type  M     value            comm
     .ent")')
      WRITE(NOUT,122)1000812050,1,1,RAUX_H(300),'EDM of Thallium in [e c
     .m]'
      WRITE(NOUT,122)2112,1,1,RAUX_H(310),'EDM of neutron in [e cm] (Chi
     .ral Quark Model)'
      WRITE(NOUT,122)2112,1,1,RAUX_H(320),'EDM of neutron in [e cm] (Par
     .ton Quark Model)'
      WRITE(NOUT,122)2112,1,1,RAUX_H(330),'EDM of neutron in [e cm] (QCD
     . Sum-rule Approach)'
*      WRITE(NOUT,122)1000801990,1,1,RAUX_H(340),'EDM of Mercury in [e cm
*     .]'
*      WRITE(NOUT,122)1000801990,1,1,RAUX_H(410),'EDM of Mercury in [e cm
*     .]'
      WRITE(NOUT,122)1000801990,1,1,RAUX_H(411),'EDM of Mercury in [e cm
     .] (I)'
      WRITE(NOUT,122)1000801990,1,1,RAUX_H(412),'EDM of Mercury in [e cm
     .] (II)'
      WRITE(NOUT,122)1000801990,1,1,RAUX_H(413),'EDM of Mercury in [e cm
     .] (III)'
      WRITE(NOUT,122)1000801990,1,1,RAUX_H(414),'EDM of Mercury in [e cm
     .] (IV)'
      WRITE(NOUT,122)1000882250,1,1,RAUX_H(420),'EDM of Radium in [e cm]
     .'
      WRITE(NOUT,122)1000010020,1,1,RAUX_H(350),'EDM of Deuteron in [e c
     .m]'
      WRITE(NOUT,122)13,1,1,RAUX_H(360),'EDM of muon in [e cm]'
      WRITE(NOUT,122)13,2,1,RAUX_H(380),'(a_mu)_SUSY of muon'
      WRITE(NOUT,10)
*-----------------------------------------------------------------------
c HiggsBounds
      WRITE(NOUT,FMT='("Block HiggsBoundsInputHiggsCouplingsBosons")')
      WRITE(NOUT,130),CDABS(NHC_H(70,1))**2,3,25,24,24,'# H1-W-W couplin
     .g^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,2))**2,3,35,24,24,'# H2-W-W couplin
     .g^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,3))**2,3,36,24,24,'# H3-W-W couplin
     .g^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,1))**2,3,25,23,23,'# H1-Z-Z couplin
     .g^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,2))**2,3,35,23,23,'# H2-Z-Z couplin
     .g^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,3))**2,3,36,23,23,'# H3-Z-Z couplin
     .g^2, normalized to SM'
      RTMP=(CDABS(NHC_H(84,1))**2+CDABS(NHC_H(85,1))**2)
     .    /CDABS(CAUX_H(221))**2
      WRITE(NOUT,130),RTMP,3,25,21,21,'# H1-gluon-gluon coupling^2, norm
     .alized to SM'
      RTMP=(CDABS(NHC_H(84,2))**2+CDABS(NHC_H(85,2))**2)
     .    /CDABS(CAUX_H(222))**2
      WRITE(NOUT,130),RTMP,3,35,21,21,'# H2-gluon-gluon coupling^2, norm
     .alized to SM'
      RTMP=(CDABS(NHC_H(84,3))**2+CDABS(NHC_H(85,3))**2)
     .    /CDABS(CAUX_H(223))**2
      WRITE(NOUT,130),RTMP,3,36,21,21,'# H3-gluon-gluon coupling^2, norm
     .alized to SM'
      WRITE(NOUT,130),0.D0,3,25,25,23,'# H1-H1-Z coupling^2 = 0'
      WRITE(NOUT,130),CDABS(NHC_H(70,3))**2,3,25,35,23,'# H1-H2-Z coupli
     .ng^2 = H3-W-W coupling^2, normalized to SM'
      WRITE(NOUT,130),CDABS(NHC_H(70,2))**2,3,25,36,23,'# H1-H3-Z coupli
     .ng^2 = H2-W-W coupling^2, normalized to SM'
      WRITE(NOUT,130),0.D0,3,35,35,23,'# H2-H2-Z coupling^2 = 0'
      WRITE(NOUT,130),CDABS(NHC_H(70,1))**2,3,35,36,23,'# H2-H3-Z coupli
     .ng^2 = H1-W-W coupling^2, normalized to SM'
      WRITE(NOUT,130),0.D0,3,36,36,23,'# H3-H3-Z coupling^2 = 0'
      WRITE(NOUT,10)
*
      WRITE(NOUT,FMT='("Block HiggsBoundsInputHiggsCouplingsFermions")')
*      do ih=1,3
*       print*,cdabs(nhc_h( 8,ih))**2,cdabs(nhc_h( 9,ih))**2   ! tau
*       print*,cdabs(nhc_h(17,ih))**2,cdabs(nhc_h(18,ih))**2   ! b
*       print*,cdabs(nhc_h(26,ih))**2,cdabs(nhc_h(27,ih))**2   ! top
*      enddo
      WRITE(NOUT,140),CDABS(NHC_H( 8,1))**2,CDABS(NHC_H( 9,1))**2,3,25,1
     .5,15,'# H1-tau-tau eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H( 8,2))**2,CDABS(NHC_H( 9,2))**2,3,35,1
     .5,15,'# H2-tau-tau eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H( 8,3))**2,CDABS(NHC_H( 9,3))**2,3,36,1
     .5,15,'# H3-tau-tau eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(17,1))**2,CDABS(NHC_H(18,1))**2,3,25,5
     .,5,'# H1-b-b eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(17,2))**2,CDABS(NHC_H(18,2))**2,3,35,5
     .,5,'# H2-b-b eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(17,3))**2,CDABS(NHC_H(18,3))**2,3,36,5
     .,5,'# H3-b-b eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(26,1))**2,CDABS(NHC_H(27,1))**2,3,25,6
     .,6,'# H1-top-top eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(26,2))**2,CDABS(NHC_H(27,2))**2,3,35,6
     .,6,'# H2-top-top eff. coupling^2, normmalized to SM'
      WRITE(NOUT,140),CDABS(NHC_H(26,3))**2,CDABS(NHC_H(27,3))**2,3,36,6
     .,6,'# H3-top-top eff. coupling^2, normmalized to SM'
*-----------------------------------------------------------------------
*Tail
      WRITE(NOUT,11) '--------------------------------------------------
     .----------------------------'
*-----------------------------------------------------------------------
*Formats
 10   FORMAT('#')
 11   FORMAT('#',1X,A)
 20   FORMAT('BLOCK',1X,A,2X,'#',1X,A)
 21   FORMAT('BLOCK',1X,A,1P,E16.8,2x,'#',1X,A)
 30   FORMAT(1X,I5,3X,A)
 40   FORMAT(1X,I5,1X,I5,3X,'#',1X,A)
 50   FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 60   FORMAT(1X,I2,1X,I2,3X,1P,E16.7,0P,3X,'#',1X,A)
 100  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 110  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,(I9,1X),(I9,1X),2X,'#',1X,A)
* 120  FORMAT(1X,I9,3X,I2,3X,I2,3X,1P,E16.8,0P,3X,1P,E16.8,0P,3X,I1,3X,I9
*     .,3X,I9,3x,'#',1X,A)
* 121  FORMAT(1X,I9,3X,I2,3X,I2,3X,1P,E16.8,0P,3X,1P,E16.8,0P,3X,I1,3X,I9
*     .,3X,9X,3x,'#',1X,A)
 120  FORMAT(1X,I9,3X,I2,3X,I2,3X,1P,E16.8,3X,I1,3X,I1,3X,I9
     .,3X,I9,3x,'#',1X,A)
 121  FORMAT(1X,I9,3X,I2,3X,I2,3X,1P,E16.8,3X,I1,3X,I1,3X,I9
     .,3X,9X,3x,'#',1X,A)
 122  FORMAT(1X,I12,3X,I2,3X,I1,3X,1P,E16.8,3x,'#',1X,A)
 130  FORMAT(1X,1P,E16.8,3X,I2,3X,I9,3X,I9,3X,I9,3X,A)
 140  FORMAT(1X,1P,E16.8,0P,3X,1P,E16.8,0P,3X,I2,3X,I9,3X,I9,3X,I9,3X,A)
*-----------------------------------------------------------------------
      CLOSE(NOUT)
*
      RETURN
      END
