      IMPLICIT REAL*8(A-H,L,M,O-Z)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HSELF_HDEC/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7
      real*8 initQCD5
      character*200 inFile
      character*200 outFile

      if(iargc().lt.2) then
         call getarg(0,inFile)
         write(*,fmt='(A16,1x,A16,1x,A32)')'Correct usage:',inFile,
     .   'SLHA_file.input SLHA_file.output'
         stop
      endif

      call getarg(1,inFile)
      call getarg(2,outFile)

      err=slhaRead(inFile,0)

      GF=         1.16639E-05
      ALPH=       0.00729735254
      AMTAU=      1.777
      AMMUON=     0.105658367
      AMZ=        91.1884003
      AMW=        79.9467442
      PI=4*DATAN(1D0)
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)

      zero=0
      MTOP=     slhaVal1('SMINPUTS',zero,6)
      MA=       slhaVal1('EXTPAR',zero,26)
      MGLU=     slhaVal1('MASS',zero,1000021)  
      MCHI0=    slhaVal1('MASS',zero,1000037)
      MGLU=     slhaVal1('MASS',zero,1000021)  

      qcdL=initQCD5(slhaVal1('SMINPUTS',zero,3), 
     . 1.4D0 , slhaVal1('SMINPUTS',zero,5), MTOP)
    
      susy=dsqrt(slhaVal1('MASS',zero,1000006)
     . *slhaVal1('MASS',zero,2000006))
      AMSQ= slhaVal1('MSOFT',susy,43)
      AMUR= slhaVal1('MSOFT',susy,46)
      AMDR= slhaVal1('MSOFT',susy,49)
c     susy = DSQRT(2*AMSQ**2+AMUR**2+AMDR**2)/2

      TANB=     slhaVal1('HMIX',susy,2)
      MU=       slhaVal1('HMIX',susy,1) 
      AD=slhaVal2('AD',susy,3,3)
      AU=slhaVal2('AU',susy,3,3)
 

      ft=1.017363287
        BET=DATAN(TANB)
        SB = DSIN(BET)
        CB = DCOS(BET)

C  ============ TRANSFORMATION OF INPUT FOR SUBH ========== 
      CF = 4/3.D0
      Q0 = susy
      ALSP = alphaQCD(Q0)/PI
      ALTP = MtRun(Q0)**2/2/PI/V**2/SB**2 / PI *ft
      ALBP = MbRun(Q0)**2/2/PI/V**2/CB**2 / PI
      RMT = MtRun(MTOP)*ft
      RMB = MbRun(MTOP)

      QT = DSQRT(DMAX1(AMSQ**2+RMT**2,AMUR**2+RMT**2))
      QB = DSQRT(DMAX1(AMSQ**2+RMB**2,AMDR**2+RMB**2))
      AMH12 = MA**2*SB**2 - AMZ**2/2*(CB**2-SB**2) - AMU**2
      AMH22 = MA**2*CB**2 + AMZ**2/2*(CB**2-SB**2) - AMU**2
      XB = AMSQ**2 + AMDR**2 + AMH12 + AD**2
      XT = AMSQ**2 + AMUR**2 + AMH22 + AU**2
      AD0 = AD + (CF*ALSP*MGLU + 3*ALTP/2*AU + ALBP/4*AD)
     .         * DLOG(QB**2/Q0**2)
      AMDL0 = DSQRT(AMSQ**2 + (-CF*ALSP*MGLU**2 + (ALTP*XT+ALBP*XB)/4)
     .                      * DLOG(QB**2/Q0**2))
      AMDR0 = DSQRT(AMDR**2 + (-CF*ALSP*MGLU**2 + ALBP*XB/4)
     .                      * DLOG(QB**2/Q0**2))
      AU0 = AU + (CF*ALSP*MGLU + ALTP/4*AU + 3*ALBP/2*AD)
     .         * DLOG(QT**2/Q0**2)
      AMUL0 = DSQRT(AMSQ**2 + (-CF*ALSP*MGLU**2 + (ALTP*XT+ALBP*XB)/4)
     .                      * DLOG(QT**2/Q0**2))
      AMUR0 = DSQRT(AMUR**2 + (-CF*ALSP*MGLU**2 + ALTP*XT/4)
     .                      * DLOG(QT**2/Q0**2))

#ifdef NNN

        write(*,*) 'MA', MA,         1410.
        write(*,*) 'TANB=',TANB,     45.8286811
        write(*,*) 'AMUL0=', AMUL0,  2279.80856
        write(*,*) 'AMDL0=', AMDL0,  2278.37708
        write(*,*) 'AMUR0=', AMUR0,  1059.57155
        write(*,*) 'AMDR0=', AMDR0,  2701.35256
        write(*,*) 'MTOP=',  MTOP,   175.
        write(*,*) 'AU0=',   AU0,    2931.98355
        write(*,*) 'AD0=',   AD0,    144.30904
        write(*,*) 'MU=',    MU,     1844.55714
        write(*,*) 'MCHI0=', MCHI0,  1842.36806
        write(*,*) 'MGLU=',  MGLU,   1639.63265

#endif
            
      call SUBH1_HDEC(MA,TANB,AMUL0,AMDL0,AMUR0,AMDR0,MTOP,AU0,AD0,
     * MU,MCHI0, MHP,HMP,MCH,SA,CA,TANBA,MGLU,dmb)
              
c      write(*,fmt=
c     .'("HdecayLambdas",1PE12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4)') 
c     . Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7

      open(22,FILE=outFile,STATUS='UNKNOWN')

      write(22,fmt='("BLOCK LAMBDAS")')
      write(22,fmt='(1x,I3,1PE12.4)') 0,DMBAPP_HDEC()
      write(22,fmt='(1x,I3,1PE12.4)') 1,Lambda1
      write(22,fmt='(1x,I3,1PE12.4)') 2,Lambda2
      write(22,fmt='(1x,I3,1PE12.4)') 3,Lambda3
      write(22,fmt='(1x,I3,1PE12.4)') 4,Lambda4
      write(22,fmt='(1x,I3,1PE12.4)') 5,Lambda5
      write(22,fmt='(1x,I3,1PE12.4)') 6,Lambda6
      write(22,fmt='(1x,I3,1PE12.4)') 7,Lambda7
      close(22)

      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", Nucl. Phys.
c       B461 (1996) 407. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SUBH1_HDEC(MA,TANB,MUL,MDL,MUR,MD,MTOP,AU,AD,MU,MCHI0,
     *                 MHP,HMP,MCH,SA,CA,TANBA,MGLU,DELTAMB)

      IMPLICIT REAL*8(A-H,L,M,O-Z)
      real*8 MbRun,MtRun
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HSELF_HDEC/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

      MCHI = MCHI0
      TANBA = TANB
      TANBT = TANB
      
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = alphaQCD(AMZ)
      ALPHA3  = alphaQCD(MTOP)
      MB      = MbRun(MTOP)
      RMTOP   = MtRun(MTOP)
      TUL = LOG((MUL**2+MTOP**2)/MTOP**2)
      TDL = LOG((MDL**2+MTOP**2)/MTOP**2)
      TU = LOG((MUR**2 + MTOP**2)/MTOP**2)
      TD = LOG((MD**2 + MTOP**2)/MTOP**2)
      SINB = TANB/DSQRT(1.D0 + TANB**2)
      COSB = SINB/TANB

      IF(MA.GT.MTOP)
     *       TANBA = TANB*(1.D0-3.D0/32.D0/PI**2*
     *       (RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *       DLOG(MA**2/MTOP**2))
      IF(MA.LT.MTOP.OR.MA.EQ.MTOP) TANBT = TANBA

      SINB = TANBT/DSQRT(1.D0 + TANBT**2)
      COSB = 1.D0/DSQRT(1.D0 + TANBT**2)
      COS2B = (TANBT**2 - 1.D0)/(TANBT**2 + 1.D0)
      G1 = DSQRT(ALPHA1*4.D0*PI)
      G2 = DSQRT(ALPHA2*4.D0*PI)
      G3 = DSQRT(ALPHA3*4.D0*PI)
      HU = RMTOP/V/SINB
      HD =  MB/V/COSB
C

      IF(MUL.GT.MUR) TP = TUL - TU
      IF(MUL.LT.MUR.OR.MUL.EQ.MUR) TP = TU - TUL
      IF(MUL.GT.MUR) TDP = TU
      IF(MUL.LT.MUR.OR.MUL.EQ.MUR) TDP = TUL
      IF(MDL.GT.MD) TPD = TDL - TD
      IF(MDL.LT.MD.OR.MDL.EQ.MD) TPD = TD - TDL
      IF(MDL.GT.MD) TDPD = TD
      IF(MDL.LT.MD.OR.MDL.EQ.MD) TDPD = TDL

      IF(MDL.GT.MD) DLAMBDA1 = 6./96./PI**2*G1**2*HD**2*TPD
      IF(MDL.LT.MD.OR.MDL.EQ.MD) DLAMBDA1 = 3./32./PI**2*
     * HD**2*(G1**2/3.+G2**2)*TPD

      IF(MUL.GT.MUR) DLAMBDA2 =12./96./PI**2*G1**2*HU**2*TP
      IF(MUL.LT.MUR.OR.MUL.EQ.MUR) DLAMBDA2 = 3./32./PI**2*
     * HU**2*(-G1**2/3.+G2**2)*TP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  dlambdap1 and dlambdap2 are the new log corrections due to
c  the presence of the gluino mass. They are in general very small,
c  and only present if there is a hierarchy of masses between the
c  two stops.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        dlambdap2 = 0
        tglu = log(mglu**2/mtop**2)

        if(mglu.lt.mur.or.mglu.lt.mul) then
        if(mul.gt.mur.and.mglu.gt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tglu**2)
        endif

        if(mul.gt.mur.and.mglu.lt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tu**2)
        endif

        if(mul.gt.mur.and.mglu.eq.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tu**2)
        endif

        if(mur.gt.mul.and.mglu.gt.mul) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tglu**2)
        endif

        if(mur.gt.mul.and.mglu.lt.mul) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tul**2)
        endif

        if(mur.gt.mul.and.mglu.eq.mul) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tul**2)
        endif
        endif

      DLAMBDA3 = 0.
      DLAMBDA4 = 0.

      IF(MDL.GT.MD) DLAMBDA3 = -1./32./PI**2*G1**2*HD**2*TPD
      IF(MDL.LT.MD.OR.MDL.EQ.MD) DLAMBDA3 = 3./64./PI**2*HD**2*
     *(G2**2-G1**2/3.)*TPD
      
      IF(MUL.GT.MUR) DLAMBDA3 = DLAMBDA3 - 
     *1./16./PI**2*G1**2*HU**2*TP
      IF(MUL.LT.MUR.OR.MUL.EQ.MUR) DLAMBDA3 = DLAMBDA3 + 
     * 3./64./PI**2*HU**2*(G2**2+G1**2/3.)*TP

      IF(MUL.LT.MUR) DLAMBDA4 = -3./32./PI**2*G2**2*HU**2*TP
      IF(MDL.LT.MD) DLAMBDA4 = DLAMBDA4 - 3./32./PI**2*G2**2*
     *                        HD**2*TPD
C
      LAMBDA1 = ((G1**2 + G2**2)/4.)*
     *(1.-3.*HD**2*(TPD + TDPD)/8./PI**2)
     *+(3.*HD**4./16./PI**2) *TPD*(1.   
     *+ (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * (TPD + 2.*TDPD)/16./PI**2) 
     *+(3.*HD**4./8./PI**2) *TDPD*(1.  + (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * TDPD/16./PI**2) + DLAMBDA1 
C
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*
     *(TP + TDP)/8./PI**2)
     *+(3.*HU**4./16./PI**2) *TP*(1.   
     *+ (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * (TP + 2.*TDP)/16./PI**2) 
     *+(3.*HU**4./8./PI**2) *TDP*(1. + (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * TDP/16./PI**2) + DLAMBDA2  + DLAMBDAP2
C
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2)*(TP + TDP)/16./PI**2 -3.*
     *(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA3 
C
      LAMBDA4 = (- G2**2/2.)*(1.
     *-3.*(HU**2)*(TP + TDP)/16./PI**2
     *-3.*(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA4
C     
	LAMBDA5 = 0.
	LAMBDA6 = 0.
	LAMBDA7 = 0.
C
C     THIS IS THE CONTRIBUTION FROM LIGHT CHARGINOS/NEUTRALINOS
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  	 MSSUSY=DSQRT(0.5D0*(MUL**2+MUR**2)+MTOP**2)
	IF(MCHI.GT.MSSUSY)GOTO 3790
	IF(MCHI.LT.MTOP) MCHI=MTOP
	TCHAR=LOG(MSSUSY**2/MCHI**2)
	DELTAL12=(9./64./PI**2*G2**4+5./192./PI**2*G1**4)*TCHAR
	DELTAL3P4=(3./64./PI**2*G2**4+7./192./PI**2*G1**4
     *       +4./32/PI**2*G1**2*G2**2)*TCHAR
	DELTAM112=2.*DELTAL12*V**2*COSB**2
	DELTAM222=2.*DELTAL12*V**2*SINB**2
	DELTAM122=2.*DELTAL3P4*V**2*SINB*COSB
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
        DLAM1 = DELTAM112/2.D0/V**2/COSB**2
        DLAM2 = DELTAM222/2.D0/V**2/SINB**2
        DLAM3 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(G1**2-G2**2)/(G1**2+G2**2)
        DLAM4 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(2*G2**2)/(G1**2+G2**2)
        LAMBDA1 = LAMBDA1+DLAM1
        LAMBDA2 = LAMBDA2+DLAM2
        LAMBDA3 = LAMBDA3+DLAM3
        LAMBDA4 = LAMBDA4+DLAM4
C--END OF EXTENSION
 3790	CONTINUE
CCCCCCCCCCCCCCC    END OF CHARGINOS AND NEUTRALINOS  CCCCCCCCCCCC 


C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I

      CALL GFUN_HDEC(MA,TANBA,MUL,MDL,MUR,MD,MTOP,AU,AD,MU,MGLU,
     *              DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7,deltamb)

      LAMBDA1 = LAMBDA1+DLAM1
      LAMBDA2 = LAMBDA2+DLAM2
      LAMBDA3 = LAMBDA3+DLAM3
      LAMBDA4 = LAMBDA4+DLAM4
      LAMBDA5 = LAMBDA5+DLAM5
      LAMBDA6 = LAMBDA6+DLAM6
      LAMBDA7 = LAMBDA7+DLAM7
      
      M2(1,1) = 2.*V**2*(LAMBDA1*COSB**2+2.*LAMBDA6*
     *COSB*SINB + LAMBDA5*SINB**2) + MA**2*SINB**2
      M2(2,2) = 2.*V**2*(LAMBDA5*COSB**2+2.*LAMBDA7*
     *COSB*SINB + LAMBDA2*SINB**2) + MA**2*COSB**2
      M2(1,2) = 2.*V**2*(LAMBDA6*COSB**2+(LAMBDA3+LAMBDA4)*
     *COSB*SINB + LAMBDA7*SINB**2) - MA**2*SINB*COSB
      M2(2,1) = M2(1,2)

      M2P(1,1) = M2(1,1)
      M2P(2,2) = M2(2,2)
      M2P(1,2) = M2(1,2)
      M2P(2,1) = M2(2,1)

C--END OF EXTENSION

      TRM2P  = M2P(1,1) + M2P(2,2)
      DETM2P = M2P(1,1)*M2P(2,2) - M2P(1,2)*M2P(2,1)

      MH2P = (TRM2P - DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
      HM2P = (TRM2P + DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
C !!!!!!!!!!!!!!!!!!!
      MCH2=MA**2+(LAMBDA5-LAMBDA4)*V**2
C !!!!!!!!!!!!!!!!!!!
      MCH=DSQRT(MCH2)
      HMP = DSQRT(HM2P) 
      IF(MH2P.LT.0.)GOTO 5555
      MHP = DSQRT(MH2P) 
C
      SIN2ALPHA = 2.*M2P(1,2)/DSQRT(TRM2P**2-4.D0*DETM2P)
      COS2ALPHA = (M2P(1,1)-M2P(2,2))/DSQRT(TRM2P**2-4.D0*DETM2P)
      IF(COS2ALPHA.GT.0.) ALPHA = DASIN(SIN2ALPHA)/2.D0
      IF(COS2ALPHA.LT.0.) ALPHA = -PI/2.D0-DASIN(SIN2ALPHA)/2.D0
      SA = DSIN(ALPHA)
      CA = DCOS(ALPHA)  
      SQBMA = (SINB*CA - COSB*SA)**2

5555  RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCC NON DEGENERATE STOP/SBOTTOM EFFECTS CCCCCCCCC
C
        SUBROUTINE GFUN_HDEC(MA,TANB,MUL,MDL,MUR,MD,MTOP,AT,AB,MU,MGLU,
     *                DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7,deltamb)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        real*8 MbRun,MtRun
        DIMENSION VH(2,2),VH1(2,2),VH2(2,2),
     *            VH3T(2,2),VH3B(2,2),AL(2,2)
        COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
        G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y)

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MUL2  = MUL**2
        MDL2  = MDL**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB
     
       MB = MbRun(MTOP)
  
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = alphaQCD(AMZ)
      ALPHA3  = alphaQCD(MTOP)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MUL.GT.MUR) MST = MUL
        IF(MUR.GT.MUL.OR.MUR.EQ.MUL) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MDL.GT.MD) MSB = MDL
	IF(MD.GT.MDL.OR.MD.EQ.MDL) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        RMTOP   = MtRun(MTOP) *1.017363287

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        AL(1,1) = AL1
        AL(1,2) = (AL2+AL1)/2.
        AL(2,1) = (AL2+AL1)/2.
        AL(2,2) = AL2

	IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
	ELSE
	VI =  V
	H1I = VI*COSB
	H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
	END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

      CALL DELMB_HDEC(MA,TANB,MUL,MDL,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)

        IF(STOP22.LT.0.) GOTO 4237
        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
        MBOT2 = DSQRT(MBOT4)
        mb = mb/(1+deltamb)

        VH1(1,1) = 1./TANBST
        VH1(2,1) = -1.
        VH1(1,2) = -1.
        VH1(2,2) = TANBST
        VH2(1,1) = TANBST
        VH2(1,2) = -1.
        VH2(2,1) = -1.
        VH2(2,2) = 1./TANBST

C CCCCCCCCCCCCCCCCCCCCCCCCCCC  D-terms CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	STW=SW

	F1T=(MUL2-MUR2)/(STOP12-STOP22)*(.5-4./3.*STW)*
     *         LOG(STOP1/STOP2)
     *        +(.5-2./3.*STW)*LOG(STOP1*STOP2/(MUL2+MTOP2))
     *        + 2./3.*STW*LOG(STOP1*STOP2/(MUR2+MTOP2))

	F1B=(MDL2-MD2)/(SBOT12-SBOT22)*(-.5+2./3.*STW)*
     *        LOG(SBOT1/SBOT2)
     *        +(-.5+1./3.*STW)*LOG(SBOT1*SBOT2/(MDL2+MBOT2))
     *        - 1./3.*STW*LOG(SBOT1*SBOT2/(MD2+MBOT2))

	F2T=1/(STOP12-STOP22)*
     *         (-.5*LOG(STOP12/STOP22)
     *        +(4./3.*STW-.5)*(MUL2-MUR2)/(STOP12-STOP22)*
     *         G(STOP12,STOP22))

	F2B=1/(SBOT12-SBOT22)*
     *         (.5*LOG(SBOT12/SBOT22)
     *        +(-2./3.*STW+.5)*(MDL2-MD2)/(SBOT12-SBOT22)*
     *        G(SBOT12,SBOT22))

C*************************************************************
C
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
C
C TRAFOS APPROXIMATE -> EXACT:
C
C (i)  1/M_{SUSY}^2 -> LOG(M1^2/M2^2) / (M1^2-M2^2)
C
C (ii) 1/M_{SUSY}^4 -> -6 G(M1^2,M2^2) / (M1^2-M2^2)^2
C
C Then use results of Phys. Lett. B355 (1995) 209 in order to
C obtain the results for lambda_1 - lambda_7 according to
C Nucl. Phys. B461 (1996) 407. Perform a full evolution from
C M_SUSY -> m_t for lambdas (anomalous dimensions, v_i).
C
C - ht^2*hb^2 terms neglected in lambda_3,4 (according to
C   Nucl. Phys. B461 (1996) 407)
C
C*************************************************************

        DLAM1T = MTOP4/(SINBT**4)*(MU**2/(STOP1**2
     *    -STOP2**2))**2*G(STOP12,STOP22)
     *  - MZ**2*MTOP2*MU**2/TANBST**2*F2T/COSBT**2

        DLAM1B = MBOT4/(COSBB**4)*(LOG(SBOT1**2*SBOT2**2/
     *    (MDL2+MBOT2)/(MD2+MBOT2))
     *    + 2*AB**2/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2))
     *  + MBOT4/(COSBB**4)*(AB**2/
     *    (SBOT1**2-SBOT2**2))**2*G(SBOT12,SBOT22)
     *  + MZ**2*(2*MBOT2*F1B-MBOT2*AB**2*F2B)/COSBB**2

        DLAM2T = MTOP4/(SINBT**4)*(LOG(STOP1**2*STOP2**2/
     *    (MUL2+MTOP2)/(MUR2+MTOP2))
     *  + 2*AT**2/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2))
     *  + MTOP4/(SINBT**4)*(AT**2/
     *    (STOP1**2-STOP2**2))**2*G(STOP12,STOP22)
     *  + MZ**2*(-2*MTOP2*F1T+MTOP2*AT**2*F2T)/SINBT**2
 
        DLAM2B = MBOT4/(COSBB**4)*MU**4/(SBOT1**2
     *    -SBOT2**2)**2*G(SBOT12,SBOT22)
     *    + MZ**2*MBOT2*MU**2*TANBSB**2*F2B/SINBB**2
 
        DLAM3T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  + MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM3B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM4T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  - MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM4B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM5T = MTOP4/(SINBT**4)*
     *    (MU**2*AT**2)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)

        DLAM5B = MBOT4/(COSBB**4)*
     *    (MU**2*AB**2)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)

        DLAM6T = MTOP4/(SINBT**4)*
     *    (-MU**3*AT)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)
     *  + MZ**2*MTOP2*MU*AT/TANBST*F2T/(2*SINBT*COSBT)

        DLAM6B = MBOT4/(COSBB**4)*MU*AB*
     *    (-1.D0/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2)
     *    -AB**2/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22))
     *  - MZ**2*(-MBOT2*AB*MU*TANBSB*F2B)/(2*SINBB*COSBB)

        DLAM7T = MTOP4/(SINBT**4)*MU*AT*
     *    (-1.D0/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2)
     *    -AT**2/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22))
     *  - MZ**2*MTOP2*AT*MU/TANBST*F2T/(2*SINBT*COSBT)

        DLAM7B = MBOT4/(COSBB**4)*
     *    (-MU**3*AB)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)
     *    - MZ**2*MBOT2*MU*AB*TANBSB*F2B/(2*SINBB*COSBB)

       TQ = LOG((MUL2 + MTOP2)/MTOP2)
       TU = LOG((MUR2+MTOP2)/MTOP2)
       TQD = LOG((MDL2 + MB**2)/MB**2)
       TD = LOG((MD2+MB**2)/MB**2)

        FACT = 3.D0/(16.D0*PI**2*(H1T**2+H2T**2)**2)
        FACB = 3.D0/(16.D0*PI**2*(H1B**2+H2B**2)**2)

        DLAM1 = FACT*DLAM1T*(1.-AL1*TT) + FACB*DLAM1B*(1.-AL1*TB)

        DLAM2 = FACT*DLAM2T*(1.-AL2*TT) + FACB*DLAM2B*(1.-AL2*TB)

        DLAM3 = FACT*DLAM3T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM3B*(1.-(AL1+AL2)/2*TB)

        DLAM4 = FACT*DLAM4T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM4B*(1.-(AL1+AL2)/2*TB)

        DLAM5 = FACT*DLAM5T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM5B*(1.-(AL1+AL2)/2*TB)

        DLAM6 = FACT*DLAM6T*(1.-(3*AL1+AL2)/4*TT)
     *        + FACB*DLAM6B*(1.-(3*AL1+AL2)/4*TB)

        DLAM7 = FACT*DLAM7T*(1.-(AL1+3*AL2)/4*TT)
     *        + FACB*DLAM7B*(1.-(AL1+3*AL2)/4*TB)

        FACTOR = 1.D0
        DLAM1 = DLAM1 * FACTOR
        DLAM2 = DLAM2 * FACTOR
        DLAM3 = DLAM3 * FACTOR
        DLAM4 = DLAM4 * FACTOR
        DLAM5 = DLAM5 * FACTOR
        DLAM6 = DLAM6 * FACTOR
        DLAM7 = DLAM7 * FACTOR

C--END OF EXTENSION

        GOTO 4236
 4237   CONTINUE

        DLAM1 = -1.D+15
        DLAM2 = -1.D+15
        DLAM3 = -1.D+15
        DLAM4 = -1.D+15
        DLAM5 = -1.D+15
        DLAM6 = -1.D+15
        DLAM7 = -1.D+15

4236    RETURN
        END


      SUBROUTINE DELMB_HDEC(MA,TANB,MUL,MDL,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        real*8 MbRun,MtRun

        COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MUL2  = MUL**2
        MDL2  = MDL**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      RMTOP = MtRun(MTOP)
      MB = MbRun(MTOP)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = alphaQCD(AMZ)
      ALPHA3  = alphaQCD(MTOP)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MUL.GT.MUR) MST = MUL
        IF(MUR.GT.MUL.OR.MUR.EQ.MUL) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MDL.GT.MD) MSB = MDL
	IF(MD.GT.MDL.OR.MD.EQ.MDL) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
        ELSE
        VI =  V
        H1I = VI*COSB
        H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
        END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

        deltamt = 0
        deltamb = 0

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MUL2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MUL2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MUL2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MUL2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MDL2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MDL2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MDL2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MDL2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Here is the definition of deltamb and deltamt, which
c     are the vertex corrections to the bottom and top quark
c     mass, keeping the dominant QCD and top Yukawa coupling
c     induced corrections.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        deltamb = -2*alpha3/3./pi*mglu*(ab-mu*tanb)*
     *  T_HDEC(sbot1,sbot2,mglu)
     *  + ht**2/(4.*pi)**2*(at-mu/tanb)*mu*tanb*
     *  T_HDEC(stop1,stop2,mu)


        deltamt = -2.*alpha3/3./pi*(at-mu/tanb)*mglu*
     *  T_HDEC(stop1,stop2,mglu)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Here the new values of the top and bottom quark masses at
c   the scale MS are defined, to be used in the effective
c   potential approximation. They are just the old ones, but
c   including the finite corrections deltamt and deltamb.
c   The deltamb corrections can become large and are resummed
c   to all orders, as suggested in the two recent works by M. Carena,
c   S. Mrenna and C.E.M. Wagner, as well as in the work by M. Carena,
c   D. Garcia, U. Nierste and C.E.M. Wagner, to appear. The top
c   quark mass corrections are small and are kept in the perturbative
c   formulation. The function T(X,Y,Z) is necessary for the calculation.
c   the entries are masses and NOT their squares !
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MUL2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MUL2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MUL2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MUL2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MDL2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MDL2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MDL2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MDL2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5
4237    RETURN
        END

      FUNCTION T_HDEC(X,Y,Z)
      implicit real*8(a-h,l,m,o-z)
      if(x.eq.y) x = x - 0.00001
      if(x.eq.z) x = x - 0.00002
      if(y.eq.z) y = y - 0.00003
c       write(*,*) 'xyz',x,y,z
      T_HDEC = (X**2*Y**2*log(X**2/Y**2) + X**2*Z**2*log(Z**2/X**2)
     * + Y**2*Z**2*log(Y**2/Z**2))/((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      return
      end

c      SUBROUTINE DMBAPP_HDEC(ICASE,DGLB,DGHB,DGAB,SCALE,IL)
       function DMBAPP_HDEC()
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT NONE
       real*8 DMBAPP_HDEC
c      INTEGER ICASE,IL
c      real*8 DGLB,DGHB,DGAB,SCALE
      real*8 GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW

c global 
      real*8 amt
      real*8 AMSB1,AMSB2,AMST1,AMST2,AMG
      real*8 A,B 
      real*8 ad,au,amu
      real*8 AMSU(2),AMSD(2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      real*8 QSUSY1,QSUSY2
       
      
#ifdef QQQQ
      DIMENSION AMCHAR(2),AMNEUT(4),XMNEUT(4),
     .          XMST(2),XMSB(2),AMSL(2),
     .          AMSU(2),AMSD(2),AMSE(2),AMSN(2)
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/HMASS_HDEC/AMSM,AMA,AMHL,AMHH,AMCH,AMAR
      COMMON/GLUINO_HDEC/AMG,AMSB1,AMSB2,STH,CTH,
     .              GLBB(2,2),GHBB(2,2),GABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              GLTT(2,2),GHTT(2,2),GATT(2,2)
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/SMASS_HDEC/AMNEUT,XMNEUT,AMCHAR,XMST,XMSB,AMSL,
     .              AMSU,AMSD,AMSE,AMSN 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      COMMON/DAVID/QSUSY1,QSUSY2,LOOP
      COMMON/ALSLAM_HDEC/XLB1(6),XLB2(6),XLB3(6)
#endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c Local 

      double precision 
     . ab,ash,at,celw,cqcd,deltamb,felw,fqcd,ht,pi,rmtop,sb,
     .sbot1,sbot2,scalelw,scalqcd,stop1,stop2,tana,tanb,v,xxt,q

c Functions
      double precision  mtrun,alphaqcd,felw_hdec,fqcd_hdec,t_hdec,
     .slhaVal0,slhaVal1,slhaVal2

      amt=    slhaVal1('SMINPUTS',0.D0,6)
      AMSB1=  slhaVal1('MASS',0.D0,1000005)
      AMSB2=  slhaVal1('MASS',0.D0,2000005)
      AMST1=  slhaVal1('MASS',0.D0,1000006)
      AMST2=  slhaVal1('MASS',0.D0,2000006)
      AMG=    slhaVal1('MASS',0.D0,1000021)

      QSUSY1=1
      QSUSY2=1
      Q=sqrt(AMST1*AMST2)
      A=      slhaVal0('ALPHA',0.D0)
      TANB=      slhaVal1('HMIX',Q  ,2) 
      ad=     slhaVal2('AD',  Q  ,3,3)
      au=     slhaVal2('AU',  Q   ,3,3)
      amu=    slhaVal1('HMIX',Q  ,1)
      AMSU(1)=slhaVal1('MASS',0.D0,1000002)
      AMSU(2)=slhaVal1('MASS',0.D0,2000002)
      AMSD(1)=slhaVal1('MASS',0.D0,1000001)
      AMSD(2)=slhaVal1('MASS',0.D0,2000001)
      
c      write(*,*) 'LMBD:',amt,AMSB1,AMSB2,AMST1,AMST2,AMG,A,TANB, 
c     .ad,au,amu,AMSU,AMSD
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   
      PI = 4*DATAN(1.D0)
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
c      TANB = DTAN(B)
      TANA = DTAN(A)
      SB = TANB/DSQRT(1+TANB**2)
      AT = AU
      AB = AD
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SCALELW = QSUSY2*(AMST1+AMST2+DABS(AMU))/3
      SCALQCD = QSUSY1*(AMSB1+AMSB2+DABS(AMG))/3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      RMTOP   = MtRun(SCALELW)*1.017363287
      HT = RMTOP/V/SB
      STOP1 = AMST1
      STOP2 = AMST2
      SBOT1 = AMSB1
      SBOT2 = AMSB2

      FELW = 1
      FQCD = 1
c      IF(IL.EQ.2)THEN
       ASH = alphaQCD(SCALELW)
       CELW = FELW_HDEC(SCALELW,AMU,AMG,SBOT1,SBOT2,STOP1,STOP2,AMT)
       FELW = 1+ASH/PI*CELW
       ASH = alphaQCD(SCALQCD)
       XXT = AU - AMU/TANB
       CQCD = FQCD_HDEC(SCALQCD,AMT,AMG,SBOT1,SBOT2,STOP1,STOP2,
     .                  AMSU(1),AMSU(2),AMSD(1),AMSD(2),XXT)
       FQCD = 1+ASH/PI*CQCD
c       write(*,*) 'FQCD_HDEC(',SCALQCD,AMT,AMG,SBOT1,SBOT2,STOP1,STOP2,
c     .                  AMSU(1),AMSU(2),AMSD(1),AMSD(2),XXT,')'
c       write(*,*) 'ASH=',ASH,'PI=',PI,'CQCD=',CQCD
c      ENDIF
 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     FQCD = 1
c     FELW = 1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ASH = alphaQCD(SCALQCD)
       DELTAMB = 2*ASH/3/PI*AMG*AMU*TANB*T_HDEC(SBOT1,SBOT2,AMG)*FQCD
     *         /(1-2*ASH/3/PI*AMG*AB*T_HDEC(SBOT1,SBOT2,AMG))
     * + HT**2/(4*PI)**2*AT*AMU*TANB*T_HDEC(STOP1,STOP2,AMU)*FELW

c       write(*,*) 'HT=',HT,'FQCD=',FQCD,'FELW=',FELW
       DMBAPP_HDEC=DELTAMB
      END
       
      double precision function t134p_hdec(am1,am2,mu)
      implicit double precision (a-h,o-z)
      double precision m1,m2,mu,ll1,ll2
      complex*16 li2_hdec
      sp(a) = dreal(li2_hdec(dcmplx(a,0.d0)))
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6
      if(am1.lt.am2)then
       m1 = am1
       m2 = am2
      else
       m1 = am2
       m2 = am1
      endif
      ll1 = dlog(mu**2/m1**2)
      ll2 = dlog(mu**2/m2**2)
      if(m1.eq.m2)then
       dummy = 7*(m1**2+m2**2)/2
     .       + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .       - m1**2/2*dlog(m1**2/m2**2)**2
      else
       dummy = 7*(m1**2+m2**2)/2
     .       + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .       + (m1**2-m2**2)*(dlog(m1**2/m2**2)*dlog(1-m1**2/m2**2)
     .                       + sp(m1**2/m2**2)-zeta2)
     .       - m1**2/2*dlog(m1**2/m2**2)**2
      endif
      t134p_hdec = dummy/mu**2
      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function t134_hdec(am1,am2,am3,mu)
      implicit double precision (a-h,o-z)
      double precision m1,m2,m3,mu,ll1,ll2,ll3
      complex*16 sp,li2_hdec,lam,ll,cx,cy,phi
      sp(cx) = li2_hdec(cx)
      lam(cx,cy) = (1-cx-cy)**2-4*cx*cy
      phi(cx,cy,ll) = (2*cdlog((1+cx-cy-ll)/2)*cdlog((1-cx+cy-ll)/2)
     .              -cdlog(cx)*cdlog(cy)+2*zeta2
     .              -2*sp((1+cx-cy-ll)/2)-2*sp((1-cx+cy-ll)/2))/ll
      eps = 1.d-15
      rim = dcmplx(1.d0,eps)
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6
      m1 = dmin1(am1,am2,am3)
      m3 = dmax1(am1,am2,am3)
      if(m1.eq.am2.and.m3.eq.am3.or.m1.eq.am3.and.m3.eq.am2) m2 = am1
      if(m1.eq.am1.and.m3.eq.am3.or.m1.eq.am3.and.m3.eq.am1) m2 = am2
      if(m1.eq.am1.and.m3.eq.am2.or.m1.eq.am2.and.m3.eq.am1) m2 = am3
      cx = m1**2/m3**2*rim
      cy = m2**2/m3**2*rim
      ll1 = dlog(mu**2/m1**2)
      ll2 = dlog(mu**2/m2**2)
      ll3 = dlog(mu**2/m3**2)
      ll = cdsqrt(lam(cx,cy))
      dummy = 7*(m1**2+m2**2+m3**2)/2
     .      + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .      + m3**2*(ll3**2+3*ll3)
     .      +  (m1**2-m2**2-m3**2)/4*dlog(m2**2/m3**2)**2
     .      + (-m1**2+m2**2-m3**2)/4*dlog(m1**2/m3**2)**2
     .      + (-m1**2-m2**2+m3**2)/4*dlog(m1**2/m2**2)**2
     .      + m3**2/2*lam(cx,cy)*phi(cx,cy,ll)
      t134_hdec = dummy/mu**2
      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function
     .       felw_hdec(scale,amu,amg,amsb1,amsb2,amst1,amst2,amt)
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amg,amsb1,amsb2,amst1,amst2,amu,amt
      double precision mg,mb1,mb2,mt1,mt2,mu,mt
      double precision anomalous
      double precision a
      complex*16 sp,li2_hdec,xgl,xt1,xt2,xt,xb1,xb2
      sp(r) = li2_hdec(r)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      fnorm = 1/amu**2

      cf = 4/3.d0

      rim = dcmplx(1.d0,eps)

      mt  = amt
      mu  = dabs(amu)
      mg  = amg
      mt1 = amst1
      mt2 = amst2
      mb1 = amsb1
      mb2 = amsb2

      xt  = amt**2/amu**2   * rim
      xgl = amg**2/amu**2   * rim
      xt1 = amst1**2/amu**2 * rim
      xt2 = amst2**2/amu**2 * rim
      xb1 = amsb1**2/amu**2 * rim
      xb2 = amsb2**2/amu**2 * rim

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      ans3=t134(mt2,mg,mt,mg)*xgl*(-xb1*xgl*xt1+xb1*xgl-xb1*xt*xt1+
     . xb1*xt+xb1*xt1*xt2-xb1*xt2-xb2*xgl*xt1+xb2*xgl-xb2*xt*xt1+xb2*
     . xt+xb2*xt1*xt2-xb2*xt2+2*xgl**2*xt1-2*xgl**2+2*xgl*xt*xt1-2*
     . xgl*xt-2*xgl*xt1*xt2+2*xgl*xt2)+t134(mu,mb1,mt,mg)*xgl*(-xb1*
     . xb2*xt1+xb1*xb2*xt2+xb1*xgl*xt1-xb1*xgl*xt2-xb2*xt*xt1+xb2*xt*
     . xt2+xb2*xt1-xb2*xt2+xgl*xt*xt1-xgl*xt*xt2-xgl*xt1+xgl*xt2)+
     . t134(mu,mb2,mt,mg)*xgl*(-xb1*xb2*xt1+xb1*xb2*xt2-xb1*xt*xt1+
     . xb1*xt*xt2+xb1*xt1-xb1*xt2+xb2*xgl*xt1-xb2*xgl*xt2+xgl*xt*xt1-
     . xgl*xt*xt2-xgl*xt1+xgl*xt2)+t134(mu,mg,mt,mg)*xgl*(xb1*xgl*xt1
     . -xb1*xgl*xt2+xb1*xt*xt1-xb1*xt*xt2-xb1*xt1+xb1*xt2+xb2*xgl*xt1
     . -xb2*xgl*xt2+xb2*xt*xt1-xb2*xt*xt2-xb2*xt1+xb2*xt2-2*xgl**2*
     . xt1+2*xgl**2*xt2-2*xgl*xt*xt1+2*xgl*xt*xt2+2*xgl*xt1-2*xgl*xt2
     . )
      ans2=t134(mt1,mb1,mt,mg)*xgl*(-xb1*xb2*xt2+xb1*xb2+xb1*xgl*xt2-
     . xb1*xgl-xb2*xt*xt2+xb2*xt+xb2*xt1*xt2-xb2*xt1+xgl*xt*xt2-xgl*
     . xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mb2,mt,mg)*xgl*(-xb1*xb2*xt2+
     . xb1*xb2-xb1*xt*xt2+xb1*xt+xb1*xt1*xt2-xb1*xt1+xb2*xgl*xt2-xb2*
     . xgl+xgl*xt*xt2-xgl*xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mg,mt,mg)*
     . xgl*(xb1*xgl*xt2-xb1*xgl+xb1*xt*xt2-xb1*xt-xb1*xt1*xt2+xb1*xt1
     . +xb2*xgl*xt2-xb2*xgl+xb2*xt*xt2-xb2*xt-xb2*xt1*xt2+xb2*xt1-2*
     . xgl**2*xt2+2*xgl**2-2*xgl*xt*xt2+2*xgl*xt+2*xgl*xt1*xt2-2*xgl*
     . xt1)+t134(mt2,mb1,mt,mg)*xgl*(xb1*xb2*xt1-xb1*xb2-xb1*xgl*xt1+
     . xb1*xgl+xb2*xt*xt1-xb2*xt-xb2*xt1*xt2+xb2*xt2-xgl*xt*xt1+xgl*
     . xt+xgl*xt1*xt2-xgl*xt2)+t134(mt2,mb2,mt,mg)*xgl*(xb1*xb2*xt1-
     . xb1*xb2+xb1*xt*xt1-xb1*xt-xb1*xt1*xt2+xb1*xt2-xb2*xgl*xt1+xb2*
     . xgl-xgl*xt*xt1+xgl*xt+xgl*xt1*xt2-xgl*xt2)+ans3
      ans1=log(xb1)*log(xt1)*xb1*xt1*(-xb2*xt2+xb2+xgl*xt2-xgl)+log(
     . xb1)*log(xt2)*xb1*xt2*(xb2*xt1-xb2-xgl*xt1+xgl)+log(xb2)*log(
     . xt1)*xb2*xt1*(-xb1*xt2+xb1+xgl*xt2-xgl)+log(xb2)*log(xt2)*xb2*
     . xt2*(xb1*xt1-xb1-xgl*xt1+xgl)+log(xgl)*log(xt1)*xt1*(4*xb1*xb2
     . *xt2-4*xb1*xb2-3*xb1*xgl*xt2+3*xb1*xgl-3*xb2*xgl*xt2+3*xb2*xgl
     . +2*xgl**2*xt2-2*xgl**2)+log(xgl)*log(xt2)*xt2*(-4*xb1*xb2*xt1+
     . 4*xb1*xb2+3*xb1*xgl*xt1-3*xb1*xgl+3*xb2*xgl*xt1-3*xb2*xgl-2*
     . xgl**2*xt1+2*xgl**2)+log(xt1)**2*xt1*(-xb1*xb2*xt2+xb1*xb2+xb1
     . *xgl*xt2-xb1*xgl+xb2*xgl*xt2-xb2*xgl-xgl**2*xt2+xgl**2)+4*log(
     . xt1)*xt1*(xb1*xb2*xt2-xb1*xb2-xb1*xgl*xt2+xb1*xgl-xb2*xgl*xt2+
     . xb2*xgl+xgl**2*xt2-xgl**2)+log(xt2)**2*xt2*(xb1*xb2*xt1-xb1*
     . xb2-xb1*xgl*xt1+xb1*xgl-xb2*xgl*xt1+xb2*xgl+xgl**2*xt1-xgl**2)
     . +4*log(xt2)*xt2*(-xb1*xb2*xt1+xb1*xb2+xb1*xgl*xt1-xb1*xgl+xb2*
     . xgl*xt1-xb2*xgl-xgl**2*xt1+xgl**2)+ans2
      r52=ans1/(4*(xb1*xb2*xt1**2*xt2-xb1*xb2*xt1**2-xb1*xb2*xt1*xt2
     . **2+xb1*xb2*xt1+xb1*xb2*xt2**2-xb1*xb2*xt2-xb1*xgl*xt1**2*xt2+
     . xb1*xgl*xt1**2+xb1*xgl*xt1*xt2**2-xb1*xgl*xt1-xb1*xgl*xt2**2+
     . xb1*xgl*xt2-xb2*xgl*xt1**2*xt2+xb2*xgl*xt1**2+xb2*xgl*xt1*xt2
     . **2-xb2*xgl*xt1-xb2*xgl*xt2**2+xb2*xgl*xt2+xgl**2*xt1**2*xt2-
     . xgl**2*xt1**2-xgl**2*xt1*xt2**2+xgl**2*xt1+xgl**2*xt2**2-xgl**
     . 2*xt2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      ans5=2*((2*((14*xt2-15+14*xt1+5*xt)*xt-((xt2-3)*xt2+xt1**2+(4*
     . xt2-3)*xt1))*xgl**3+(4*xt2-3+4*xt1-8*xt)*xgl**4-((2*xt+3)*(xt-
     . xt1)**2*(xt-xt2)**2+2*xgl**5)-(2*(2*((3*xt2-2)*xt2+3*xt1**2)+(
     . 5*xt2-4)*xt1)*xt-((4*(xt2-3)*xt1-3*xt2)*xt2+(4*xt2-3)*xt1**2)-
     . 2*(24*xt2-23+24*xt1)*xt**2-10*xt**3)*xgl**2-2*((2*((3*xt2-2)*
     . xt2+3*xt1**2)+(5*xt2-4)*xt1)*xt**2+(3*((xt2-1)*xt1**2-xt2**2)+
     . (3*xt2+5)*xt1*xt2)*xt-(14*xt2-15+14*xt1)*xt**3+((xt2-3)*xt1-3*
     . xt2)*xt1*xt2+4*xt**4)*xgl)*(xt1-1)*(xt2-1)-(2*(xt+xt1)*xgl-(xt
     . -xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt-1+xgl
     . )*(xt2-2+xt1)*t134(mu,mt,mg,mg)*xgl)*(xt1-xt2)
      ans4=-2*((2*((xt-xt1)**2*(xt-xt2)**2+xgl**4-(6*xt**2+10*xt*xt1+
     . 10*xt*xt2-12*xt-xt1**2-4*xt1*xt2-xt2**2)*xgl**2-2*(xt1+xt2-xt)
     . *xgl**3+2*(((2*xt2-3)*xt2+2*xt1**2+(7*xt2-3)*xt1)*xt-(5*xt2-6+
     . 5*xt1)*xt**2-(xt1+xt2)*xt1*xt2+xt**3)*xgl)*(xt1-1)*(xt2-1)+(2*
     . (xt+xt1)*xgl-(xt-xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-
     . xgl**2)*(xt2-2+xt1)*log(xgl))*(xt1-xt2)+(2*(xt+xt1)*xgl+(xt-
     . xt1)**2*(xt1-2)-xgl**2*xt1)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2
     . )*(log(xgl)-log(xt1))*(xt2-1)**2-(2*(xt+xt1)*xgl-(xt-xt1)**2-
     . xgl**2)*(2*(xt+xt2)*xgl+(xt-xt2)**2*(xt2-2)-xgl**2*xt2)*(log(
     . xgl)-log(xt2))*(xt1-1)**2)*(log(xgl)-log(xt))*xt+ans5
      ans3=-(4*((xt-xt2**2+xt2)*(xt-xt2)**2+xgl**3-((xt2+1)*xt2+xt)*
     . xgl**2-((xt+4*xt2**2)*xt-(2*xt2-1)*xt2**2)*xgl)+(((xt+xt2)*(xt
     . -xt2)*(xt2-2)-4*xt*xt2)*xgl+(xt+2*xt2)*(xt2-2)*xgl**2-(xt-xt2)
     . **2*(xt2-2)*xt-(xt2-2)*xgl**3)*(log(xgl)-log(xt2)))*(2*(xt+xt1
     . )*xgl-(xt-xt1)**2-xgl**2)*(log(xgl)-log(xt2))*(xt1-1)**2+(4*((
     . xt-xt1**2+xt1)*(xt-xt1)**2+xgl**3-((xt1+1)*xt1+xt)*xgl**2-((xt
     . +4*xt1**2)*xt-(2*xt1-1)*xt1**2)*xgl)+(((xt+xt1)*(xt-xt1)*(xt1-
     . 2)-4*xt*xt1)*xgl+(xt+2*xt1)*(xt1-2)*xgl**2-(xt-xt1)**2*(xt1-2)
     . *xt-(xt1-2)*xgl**3)*(log(xgl)-log(xt1)))*(2*(xt+xt2)*xgl-(xt-
     . xt2)**2-xgl**2)*(log(xgl)-log(xt1))*(xt2-1)**2+ans4
      ans2=(2*(((xt+2*xt2)*xt+(xt2-4)*xt2)*xgl+(xt2+2+xt)*xgl**2-(xt+
     . xt2-2)*(xt-xt2)**2-xgl**3)*(xt1-1)**2*t134(mt2,mt,mg,mg)*xgl-(
     . 2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt2-2+xt1)*(log(xgl)+4)*(
     . xgl+xt)*(xt1-xt2)*log(xgl))*(2*(xt+xt1)*xgl-(xt-xt1)**2-xgl**2
     . )-2*(((xt+2*xt1)*xt+(xt1-4)*xt1)*xgl+(xt1+2+xt)*xgl**2-(xt+xt1
     . -2)*(xt-xt1)**2-xgl**3)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(
     . xt2-1)**2*t134(mt1,mt,mg,mg)*xgl+((2*(xt1+xt2)-xgl)*xgl**3-(xt
     . -xt1)**2*(xt-xt2)**2+(2*xt**2+6*xt*xt1+6*xt*xt2-8*xt-xt1**2-4*
     . xt1*xt2-xt2**2)*xgl**2-2*(((xt2-2)*xt2+xt1**2+2*(3*xt2-1)*xt1)
     . *xt-(3*xt2-4+3*xt1)*xt**2-(xt1+xt2)*xt1*xt2)*xgl)*(log(xgl)-
     . log(xt))**2*(xt1-xt2)*(xt1-1)*(xt2-1)*xt+ans3
      ans1=-ans2
      r82=ans1/(4*((xt-xt1)**2+xgl**2-2*(xt+xt1)*xgl)*((xt-xt2)**2+
     . xgl**2-2*(xt+xt2)*xgl)*(xt1-xt2)*(xt1-1)**2*(xt2-1)**2)

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      ans6=-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(((xt1-xt2)*xt1+
     . 4*xt**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+4*log(xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2-4*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-
     . xt1)**2+xgl**2)))*xt1+(xt1-xt2)*log(-(xt-xt1+xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1)*(xt-xt1+xgl)*xt2
      ans5=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(((xt1-xt2)*xt2-
     . 4*xt**2)*(log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . )-4*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*
     . xt**2+4*log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2)))*xt2+(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt2)*(xt-xt2+xgl)*xt1+ans6
      ans4=4*(xgl+xt-xt1)*(xgl-xt-xt1)*log(xt-xt1-xgl+sqrt(-2*(xt+xt1
     . )*xgl+(xt-xt1)**2+xgl**2))*xt**2*xt2+4*(xgl+xt-xt1)*(xgl-xt-
     . xt1)*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl
     . +sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-
     . xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*
     . xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(
     . xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(xt1-xt2)*log(-(xt-
     . xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1*xt2+
     . ans5
      ans3=(2*((xt1+xt2-2*xt-2*xgl)*(xt1-xt2)+(xt1-xt2-8*xt)*log(xt2)
     . *xt2+(xt1-xt2+8*xt)*log(xt1)*xt1-6*(xt1-xt2)*log(xt)*xt-(xt1+
     . xt2+2*xt)*(xt1-xt2)*log(xgl))*xt*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)
     . *(4*xt**2-xt1*xt2+xt2**2)*log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+
     . (xt-xt2)**2+xgl**2))+(xgl+xt-xt2)*(xgl-xt-xt2)*(4*xt**2-xt1*
     . xt2+xt2**2)*log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2)))*xt1-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*
     . xt2)*log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt2-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*xt2)*log(xt+
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*xt2-4*(xgl+
     . xt-xt2)*(xgl-xt-xt2)*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2))*xt**2*xt1-4*(xgl+xt-xt2)*(xgl-xt-xt2)*log(xt-
     . xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*xt**2*xt1+
     . ans4
      ans7=((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)
      ans2=ans3*ans7
      ans1=-ans2
      rat2=ans1/(16*(xt1-xt2)**2*(xt1-1)*(xt2-1)*xt**2*xt1*xt2)

      ans5=((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . )))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-xt2)*log(xgl)
      ans4=sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(-(xt-xt2+xgl+sqrt(-2*
     . (xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))+log(-(xt-xt2+xgl-sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))))*(xt-xt2+xgl)*((xt1-1)*log(
     . xt2)*xt2-(xt2-1)*log(xt1)*xt1)+2*((log(xgl)-log(xt2))**2*(2*xt
     . -xt2)*(xt1-1)*xt2-2*(xt1-xt2)*log(xgl)**2*xt-(log(xgl)-log(xt1
     . ))**2*(2*xt-xt1)*(xt2-1)*xt1)*xt+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)
     . **2+xgl**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+
     . xgl**2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**
     . 2))-log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+log(-(xt-xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))))
     . *(xt-xt1+xgl)*((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)+ans5
      ans3=-ans4
      ans2=-(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)
     . )+log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt2-1)-2*((log(xgl)-log(xt2))*xt2-log(xgl))*(xt1-xt2)*xt)*(log
     . (xgl)-log(xt1))*xt1+ans3
      ans1=(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-1)-2*(xt1-xt2)*log(xgl)*xt)*(log(xgl)-log(xt2))*xt2+ans2
      rlt2=ans1/(8*(xt1-xt2)*(xt1-1)*(xt2-1)*xt**2)

      ans3=(2*(((log(xgl)-log(xt2))*(xt1-1)**2*xt2-(xt1-xt2)**2*log(
     . xgl)+(xt2+1-7*xt)*(xt2-1)*xt1)*xt1-(7*xt1**2-3*xt2-2*(xt2+1)*
     . xt1)*(xt2-1)*xgl)*xt1+(2*((2*(xt2+1)*xt1+3*xt2)*xt-2*xt1**3)*
     . xt1-(log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)*(xt1**2-xt2))*(xt2-1))*(log(xgl)-log(xt1))
      ans2=sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(log(xt+xt1-xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+sqrt(-2*(
     . xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt-xt1-xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*((xt1**2-xt2)*(xt2-1)*log(xt1)-
     . (xt1-1)**2*log(xt2)*xt2-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt1+xgl
     . )+((log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)+2*(5*xt+xt1+5*xgl)*xt1)*((log(xgl)-log(xt2))*(
     . xt1-1)**2*xt2-(xt1-xt2)**2*log(xgl))+ans3
      ans1=2*(((log(xgl)-log(xt2))**2*(xt1-1)**2*xt2-(xt1-xt2)**2*log
     . (xgl)**2)*(xgl+xt)+(xt2+1-5*xt)*(xt2-1)*xt1**2)*xt1-(log(xt-
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-
     . xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl
     . )*(xt1-xt2)*(xt1-1)*(xt2-1)-2*(xt+xt1+xgl)*(log(xgl)-log(xt1))
     . **2*(xt1**2-xt2)*(xt2-1)*xt1-2*(5*xt*xt2+xt1**3+5*(xt1-xt2)*(
     . xt1-1)*xgl-(5*(xt2+1)*xt-xt2)*xt1)*(xt2-1)*xt1-2*((log(xgl)-
     . log(xt1))*(xt1**2-xt2)*(xt2-1)-(log(xgl)-log(xt2))*(xt1-1)**2*
     . xt2+((xt1-xt2)*log(xgl)+(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt1+ans2
      rmt12=ans1/(4*(xt1-xt2)**2*(xt1-1)**2*(xt2-1)*xt1)

      ans3=((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . -log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(log(xgl)-log(xt2))*xt2**2+2*(5*xt+xt2+5*xgl)
     . *xt2)*(log(xgl)-log(xt1))*(xt2-1)**2*xt1+(2*(((2*xt2+3)*xt+xt2
     . **2)*xt1**2-(xt1-xt2)**2*log(xgl)*xt2-((7*xt2**2+3)*xt+2*xt2**
     . 3)*xt1)*xt2+(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2))-log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2))+log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)*(xt1-xt2**2)*(xt1-1)+2*(((7*xt-1)*xt2-2*
     . (xt-xt2**2))*xt2-((7*xt2-2)*xt2-(2*xt2+3)*xt1)*(xt1-1)*xgl)*
     . xt2)*(log(xgl)-log(xt2))
      ans2=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(xt-xt2-xgl+sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt-xt2-xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*((xt1-xt2**2)*(xt1-1)*log(xt2)+
     . (xt2-1)**2*log(xt1)*xt1-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt2+xgl
     . )-((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(5*xt+xt2+5*xgl)*xt2)*(xt1-xt2)**2*log(xgl)+
     . ans3
      ans1=(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)+2*((log(xgl)-log(xt1)
     . )**2*(xt2-1)**2*xt1-(xt1-xt2)**2*log(xgl)**2)*(xgl+xt)*xt2+2*(
     . xt+xt2+xgl)*(log(xgl)-log(xt2))**2*(xt1-xt2**2)*(xt1-1)*xt2+2*
     . (5*xt+xt2+5*xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)*xt2+2*((log(xgl)-
     . log(xt1))*(xt2-1)**2*xt1+(log(xgl)-log(xt2))*(xt1-xt2**2)*(xt1
     . -1)-((xt1-xt2)*log(xgl)-(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt2+ans2
      rmt22=ans1/(4*(xt1-xt2)**2*(xt1-1)*(xt2-1)**2*xt2)

c     relw = r22+r42+r52+r72+r82+r92+rat2+rlt2+rmt12+rmt22
      relw = r22+r42+r52+r72+r82+r92+rlt2+rmt12+rmt22
      bo  = fi(amst1**2,amst2**2,amu**2)
      relw = cf*relw/bo

      anomalous =-cf
      finscale = 2*dlog(scale**2/amg**2)
      felw_hdec = dreal(relw)*fnorm + anomalous + finscale

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     write(6,*)'dat: ',dreal(cf*rat2/bo)*fnorm
c     write(6,*)'tot: ',felw_hdec
c     write(6,*)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     write(6,*)
c     write(6,*)'ind: ',dreal(r22),dreal(r42),dreal(r52),dreal(r72),
c    . dreal(r82),dreal(r92),dreal(rat2),dreal(rlt2),dreal(rmt12),
c    . dreal(rmt22)
c     write(6,*)'     ',-1/4.d0,3/4.d0,-1/2.d0,-1/12.d0,-1/3.d0,1/6.d0,
c    .       (5+3*dlog(amg**2/amt**2))/8,0.d0,1/2.d0,1/2.d0
c     write(6,*)
c     write(6,*)'dat: ',dreal(cf*rat2/bo)*fnorm
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function fqcd_hdec(scale,amt,amg,amsb1,amsb2,
     .                         amst1,amst2,amsu1,amsu2,amsd1,amsd2,xxt)
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amt,amg,amsb1,amsb2,amst1,amst2,
     .                 amsu1,amsu2,amsd1,amsd2
      double precision mt,mg,mb1,mb2,mt1,mt2,ms1,ms2,mu
      complex*16 sp,li2_hdec,xt,xgl,xt1,xt2,xb1,xb2,xs1,xs2,xq
      double precision m,mq
      double precision anomalous
      double precision a
      complex*16 xp,xm
      sp(r) = li2_hdec(r)
      xp(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . +cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      xm(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . -cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      ca = 3
      cf = 4/3.d0
      tr = 1/2.d0
      nu = 2
      nd = 2
      nf = nu+nd+1

      fnorm = 4/amg**2

      rim = dcmplx(1.d0,eps)

      mq  = amt
      mt  = amt
      mg  = amg
      mu  = amg
      mt1 = amsb1
      mt2 = amsb2
      mb1 = amsb1
      mb2 = amsb2

      xq  = amt**2/amg**2   * rim
      xt  = amt**2/amg**2   * rim
      xgl = amg**2/amg**2   * rim
      xt1 = amsb1**2/amg**2 * rim
      xt2 = amsb2**2/amg**2 * rim
      xb1 = amsb1**2/amg**2 * rim
      xb2 = amsb2**2/amg**2 * rim

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c--comparison with luminita
c     amsu1 =  998.540d0
c     amsu2 =  999.385d0
c     amsd1 = 1001.770d0
c     amsd2 = 1001.310d0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      r12=(log(xb1)**2*xb1*(-2*xb1*xb2**2+4*xb1*xb2-2*xb1+xb2**2-2*
     . xb2+1)+4*log(xb1)*xb1*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1
     . )+log(xb2)**2*xb2*(2*xb1**2*xb2-xb1**2-4*xb1*xb2+2*xb1+2*xb2-1
     . )+4*log(xb2)*xb2*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1)+4*
     . t134p(mb1,mb1,mg)*xb1*(xb2**2-2*xb2+1)+4*t134p(mb1,mg,mg)*(-
     . xb1*xb2**2+2*xb1*xb2-xb1-xb2**2+2*xb2-1)+4*t134p(mb2,mb2,mg)*
     . xb2*(-xb1**2+2*xb1-1)+4*t134p(mb2,mg,mg)*(xb1**2*xb2+xb1**2-2*
     . xb1*xb2-2*xb1+xb2+1)+4*t134p(mg,mg,mg)*(-xb1**2+2*xb1+xb2**2-2
     . *xb2))/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*
     . xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2
     . **2-xb2))

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r32=(log(xb1)**2*xb1*(-xb2**2+2*xb2-1)+2*log(xb1)*xb1*(2*xb1*
     . xb2**2-4*xb1*xb2+2*xb1+xb2**2-2*xb2+1)+log(xb2)**2*xb2*(xb1**2
     . -2*xb1+1)+2*log(xb2)*xb2*(-2*xb1**2*xb2-xb1**2+4*xb1*xb2+2*xb1
     . -2*xb2-1)+2*t134p(mb1,mg,mg)*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2
     . **2+2*xb2-1)+2*t134p(mb2,mg,mg)*(-2*xb1**2*xb2+xb1**2+4*xb1*
     . xb2-2*xb1-2*xb2+1)+3*t134p(mg,mg,mg)*(xb1**2*xb2-xb1**2-xb1*
     . xb2**2+xb1+xb2**2-xb2)+14*(xb1**2*xb2-xb1*xb2**2-xb1+xb2))/(8*
     . (xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-
     . 2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r52=(log(xb1)**2*xb1*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2**2+2*xb2
     . -1)+4*log(xb1)*xb1*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+
     . log(xb2)**2*xb2*(-2*xb1**2*xb2+xb1**2+4*xb1*xb2-2*xb1-2*xb2+1)
     . +4*log(xb2)*xb2*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+2*
     . t134p(mb1,mb2,mg)*(xb1**2*xb2-xb1**2-xb1*xb2**2+xb1+xb2**2-xb2
     . )+2*t134p(mb1,mg,mg)*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1
     . )+2*t134p(mb2,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1)
     . )/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**
     . 2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-
     . xb2))

      r62=(log(xb1)**2*xb1*(xb2**2-2*xb2+1)+2*log(xb1)*xb1*(-xb2**2+2
     . *xb2-1)+log(xb2)**2*xb2*(-xb1**2+2*xb1-1)+2*log(xb2)*xb2*(xb1
     . **2-2*xb1+1)+2*t134p(mb1,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1+xb2
     . **2-2*xb2+1)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2-xb1**2+2*xb1*xb2+
     . 2*xb1-xb2-1)+2*t134p(mg,mg,mg)*(xb1**2*xb2+xb1**2-xb1*xb2**2-3
     . *xb1-xb2**2+3*xb2))/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1
     . **2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1
     . -xb2**3+2*xb2**2-xb2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      r82=(log(xb1)**2*(xb1*xb2**2-2*xb1*xb2+xb1-2*xb2**2+4*xb2-2)+4*
     . log(xb1)*(-xb1**2*xb2**2+2*xb1**2*xb2-xb1**2+xb1*xb2**2-2*xb1*
     . xb2+xb1+xb2**2-2*xb2+1)+log(xb2)**2*(-xb1**2*xb2+2*xb1**2+2*
     . xb1*xb2-4*xb1-xb2+2)+4*log(xb2)*(xb1**2*xb2**2-xb1**2*xb2-xb1
     . **2-2*xb1*xb2**2+2*xb1*xb2+2*xb1+xb2**2-xb2-1)+2*t134p(mb1,mg,
     . mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+2*t134p(mb2,mg,
     . mg)*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+10*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2))/(8*(xb1**3*xb2**2-2*xb1**3*
     . xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*
     . xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      mt1 = amst1
      mt2 = amst2
      xt1 = amst1**2/amg**2 * rim
      xt2 = amst2**2/amg**2 * rim

      ralsca2=(-((9*log(xb1)-10)*(xb2-1)*log(xb1)*xb1-(9*log(xb2)-10)
     . *(xb1-1)*log(xb2)*xb2))/(48*(xb1-xb2)*(xb1-1)*(xb2-1))

      ralscf2=(-((xb1-1)*log(xb2)*xb2-(xb2-1)*log(xb1)*xb1))/(8*(xb1-
     . xb2)*(xb1-1)*(xb2-1))

      rmb12=(-((((2*xb1**3-3*xb2-(xb2-6)*xb1**2-2*(xb2+1)*xb1)*xb1+(
     . xb1**2-xb2)*(xb1-1)**2*log(-(xb1-1)))*(xb2-1)+(xb1-1)**2*log(
     . xb2)*xb1**2*xb2)*log(xb1)-(((xb1**2-xb2)*(xb1+1)*(xb2-1)*log(
     . xb1)**2-(xb1-1)**2*log(xb2)**2*xb2)*xb1+((xb1-xb2)*(xb2-1)+(
     . xb1-1)*log(xb2)*xb2)*((xb1+5)*xb1+(xb1-1)**2*log(-(xb1-1)))*(
     . xb1-1))))/(4*(xb1-xb2)**2*(xb1-1)**2*(xb2-1)*xb1)

      rmb22=(-(((xb1-xb2**2)*(xb1-1)*(xb2+1)*log(xb2)**2+(xb2-1)**2*
     . log(xb1)**2*xb1)*xb2+((xb2+5)*xb2+(xb2-1)**2*log(-(xb2-1)))*(
     . xb1-xb2)*(xb1-1)*(xb2-1)-((xb2-1)**2*log(-(xb2-1))-log(xb2)*
     . xb2**2+(xb2+5)*xb2)*(xb2-1)**2*log(xb1)*xb1+((2*(xb2**2+3*xb2-
     . 1)*xb2-(xb2**2+2*xb2+3)*xb1)*xb2-(xb1-xb2**2)*(xb2-1)**2*log(-
     . (xb2-1)))*(xb1-1)*log(xb2)))/(4*(xb1-xb2)**2*(xb1-1)*(xb2-1)**
     . 2*xb2)

      rmgca2=(-((3*log(xb1)-14)*(xb1+1)*(xb2-1)**2*log(xb1)*xb1-(3*
     . log(xb2)-14)*(xb1-1)**2*(xb2+1)*log(xb2)*xb2-28*(xb1-xb2)*(xb1
     . -1)*(xb2-1)))/(16*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsu1
      ms2 = amsu2
      xs1 = amsu1**2/amg**2 * rim
      xs2 = amsu2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      ru102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rualstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rumgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsd1
      ms2 = amsd2
      xs1 = amsd1**2/amg**2 * rim
      xs2 = amsd2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      rd102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rdalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rdmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsb1
      ms2 = amsb2
      xs1 = amsb1**2/amg**2 * rim
      xs2 = amsb2**2/amg**2 * rim

      ans1=log(xb1)**2*xb1*(xb1**2*xb2-xb1**2-4*xb1*xb2**2+6*xb1*xb2-
     . 2*xb1-xb2**3+3*xb2**2-2*xb2)+2*log(xb1)*log(xb2)*xb1*xb2*(xb1
     . **2-2*xb1-xb2**2+2*xb2)+4*log(xb1)*xb1*(-2*xb1**2*xb2+2*xb1**2
     . +3*xb1*xb2**2-2*xb1*xb2-xb1+xb2**3-4*xb2**2+3*xb2)+log(xb2)**2
     . *xb2*(xb1**3+4*xb1**2*xb2-3*xb1**2-xb1*xb2**2-6*xb1*xb2+2*xb1+
     . xb2**2+2*xb2)+4*log(xb2)*xb2*(-xb1**3-3*xb1**2*xb2+4*xb1**2+2*
     . xb1*xb2**2+2*xb1*xb2-3*xb1-2*xb2**2+xb2)+2*t134p(mb1,mb2,mg)*(
     . -xb1**3+xb1**2*xb2+2*xb1**2-xb1*xb2**2-2*xb1+xb2**3-2*xb2**2+2
     . *xb2)+2*t134p(mb1,mg,mg)*(xb1**3-2*xb1**2*xb2-xb1**2+xb1*xb2**
     . 2+2*xb1*xb2-xb2**2)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2+xb1**2+2*
     . xb1*xb2**2-2*xb1*xb2-xb2**3+xb2**2)+12*(xb1**3*xb2-xb1**3-xb1*
     . xb2**3+xb1+xb2**3-xb2)
      rb102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rbalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rbmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amst1
      ms2 = amst2
      xs1 = amst1**2/amg**2 * rim
      xs2 = amst2**2/amg**2 * rim

      ans14=-8*log(xs2)*xb1*xb2**2*xs2-4*log(xs2)*xb1*xq**2*xs2-4*log
     . (xs2)*xb1*xq*xs2**2-4*log(xs2)*xb1*xq*xs2+8*log(xs2)*xb1*xs2**
     . 3-16*log(xs2)*xb1*xs2**2+8*log(xs2)*xb1*xs2-4*log(xs2)*xb2**2*
     . xq**2*xs2-4*log(xs2)*xb2**2*xq*xs2**2-4*log(xs2)*xb2**2*xq*xs2
     . +8*log(xs2)*xb2**2*xs2**3-16*log(xs2)*xb2**2*xs2**2+8*log(xs2)
     . *xb2**2*xs2+4*log(xs2)*xb2*xq**2*xs2+4*log(xs2)*xb2*xq*xs2**2+
     . 4*log(xs2)*xb2*xq*xs2-8*log(xs2)*xb2*xs2**3+16*log(xs2)*xb2*
     . xs2**2-8*log(xs2)*xb2*xs2
      ans13=-log(xs2)**2*xb1*xb2**2*xq**2*xs2+log(xs2)**2*xb1*xb2**2*
     . xs2**3-2*log(xs2)**2*xb1*xb2**2*xs2**2+log(xs2)**2*xb1*xb2**2*
     . xs2+log(xs2)**2*xb1*xq**2*xs2-log(xs2)**2*xb1*xs2**3+2*log(xs2
     . )**2*xb1*xs2**2-log(xs2)**2*xb1*xs2+log(xs2)**2*xb2**2*xq**2*
     . xs2-log(xs2)**2*xb2**2*xs2**3+2*log(xs2)**2*xb2**2*xs2**2-log(
     . xs2)**2*xb2**2*xs2-log(xs2)**2*xb2*xq**2*xs2+log(xs2)**2*xb2*
     . xs2**3-2*log(xs2)**2*xb2*xs2**2+log(xs2)**2*xb2*xs2-4*log(xs2)
     . *xb1**2*xb2*xq**2*xs2-4*log(xs2)*xb1**2*xb2*xq*xs2**2-4*log(
     . xs2)*xb1**2*xb2*xq*xs2+8*log(xs2)*xb1**2*xb2*xs2**3-16*log(xs2
     . )*xb1**2*xb2*xs2**2+8*log(xs2)*xb1**2*xb2*xs2+4*log(xs2)*xb1**
     . 2*xq**2*xs2+4*log(xs2)*xb1**2*xq*xs2**2+4*log(xs2)*xb1**2*xq*
     . xs2-8*log(xs2)*xb1**2*xs2**3+16*log(xs2)*xb1**2*xs2**2-8*log(
     . xs2)*xb1**2*xs2+4*log(xs2)*xb1*xb2**2*xq**2*xs2+4*log(xs2)*xb1
     . *xb2**2*xq*xs2**2+4*log(xs2)*xb1*xb2**2*xq*xs2-8*log(xs2)*xb1*
     . xb2**2*xs2**3+16*log(xs2)*xb1*xb2**2*xs2**2+ans14
      ans12=-4*log(xb2)**2*xb1*xb2*xq*xs1-8*log(xb2)**2*xb1*xb2*xq*
     . xs2**2+4*log(xb2)**2*xb1*xb2*xq*xs2-4*log(xb2)**2*xb1*xb2*xq+2
     . *log(xb2)**2*xb1*xb2*xs1*xs2**2-4*log(xb2)**2*xb1*xb2*xs1*xs2+
     . 2*log(xb2)**2*xb1*xb2*xs1+2*log(xb2)**2*xb1*xb2*xs2**3-4*log(
     . xb2)**2*xb1*xb2*xs2**2+2*log(xb2)**2*xb1*xb2*xs2+2*log(xb2)**2
     . *xb2*xq**3-log(xb2)**2*xb2*xq**2*xs1-5*log(xb2)**2*xb2*xq**2*
     . xs2-4*log(xb2)**2*xb2*xq**2+2*log(xb2)**2*xb2*xq*xs1*xs2+2*log
     . (xb2)**2*xb2*xq*xs1+4*log(xb2)**2*xb2*xq*xs2**2-2*log(xb2)**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb2*xq-log(xb2)**2*xb2*xs1*xs2**2+2*
     . log(xb2)**2*xb2*xs1*xs2-log(xb2)**2*xb2*xs1-log(xb2)**2*xb2*
     . xs2**3+2*log(xb2)**2*xb2*xs2**2-log(xb2)**2*xb2*xs2+log(xs2)**
     . 2*xb1**2*xb2*xq**2*xs2-log(xs2)**2*xb1**2*xb2*xs2**3+2*log(xs2
     . )**2*xb1**2*xb2*xs2**2-log(xs2)**2*xb1**2*xb2*xs2-log(xs2)**2*
     . xb1**2*xq**2*xs2+log(xs2)**2*xb1**2*xs2**3-2*log(xs2)**2*xb1**
     . 2*xs2**2+log(xs2)**2*xb1**2*xs2+ans13
      ans11=4*log(xb1)**2*xb1*xq**2-2*log(xb1)**2*xb1*xq*xs1*xs2-2*
     . log(xb1)**2*xb1*xq*xs1-4*log(xb1)**2*xb1*xq*xs2**2+2*log(xb1)
     . **2*xb1*xq*xs2-2*log(xb1)**2*xb1*xq+log(xb1)**2*xb1*xs1*xs2**2
     . -2*log(xb1)**2*xb1*xs1*xs2+log(xb1)**2*xb1*xs1+log(xb1)**2*xb1
     . *xs2**3-2*log(xb1)**2*xb1*xs2**2+log(xb1)**2*xb1*xs2+2*log(xb2
     . )**2*xb1**2*xb2*xq**3-log(xb2)**2*xb1**2*xb2*xq**2*xs1-5*log(
     . xb2)**2*xb1**2*xb2*xq**2*xs2-4*log(xb2)**2*xb1**2*xb2*xq**2+2*
     . log(xb2)**2*xb1**2*xb2*xq*xs1*xs2+2*log(xb2)**2*xb1**2*xb2*xq*
     . xs1+4*log(xb2)**2*xb1**2*xb2*xq*xs2**2-2*log(xb2)**2*xb1**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb1**2*xb2*xq-log(xb2)**2*xb1**2*xb2*
     . xs1*xs2**2+2*log(xb2)**2*xb1**2*xb2*xs1*xs2-log(xb2)**2*xb1**2
     . *xb2*xs1-log(xb2)**2*xb1**2*xb2*xs2**3+2*log(xb2)**2*xb1**2*
     . xb2*xs2**2-log(xb2)**2*xb1**2*xb2*xs2-4*log(xb2)**2*xb1*xb2*xq
     . **3+2*log(xb2)**2*xb1*xb2*xq**2*xs1+10*log(xb2)**2*xb1*xb2*xq
     . **2*xs2+8*log(xb2)**2*xb1*xb2*xq**2-4*log(xb2)**2*xb1*xb2*xq*
     . xs1*xs2+ans12
      ans10=5*log(xb1)**2*xb1*xb2**2*xq**2*xs2+4*log(xb1)**2*xb1*xb2
     . **2*xq**2-2*log(xb1)**2*xb1*xb2**2*xq*xs1*xs2-2*log(xb1)**2*
     . xb1*xb2**2*xq*xs1-4*log(xb1)**2*xb1*xb2**2*xq*xs2**2+2*log(xb1
     . )**2*xb1*xb2**2*xq*xs2-2*log(xb1)**2*xb1*xb2**2*xq+log(xb1)**2
     . *xb1*xb2**2*xs1*xs2**2-2*log(xb1)**2*xb1*xb2**2*xs1*xs2+log(
     . xb1)**2*xb1*xb2**2*xs1+log(xb1)**2*xb1*xb2**2*xs2**3-2*log(xb1
     . )**2*xb1*xb2**2*xs2**2+log(xb1)**2*xb1*xb2**2*xs2+4*log(xb1)**
     . 2*xb1*xb2*xq**3-2*log(xb1)**2*xb1*xb2*xq**2*xs1-10*log(xb1)**2
     . *xb1*xb2*xq**2*xs2-8*log(xb1)**2*xb1*xb2*xq**2+4*log(xb1)**2*
     . xb1*xb2*xq*xs1*xs2+4*log(xb1)**2*xb1*xb2*xq*xs1+8*log(xb1)**2*
     . xb1*xb2*xq*xs2**2-4*log(xb1)**2*xb1*xb2*xq*xs2+4*log(xb1)**2*
     . xb1*xb2*xq-2*log(xb1)**2*xb1*xb2*xs1*xs2**2+4*log(xb1)**2*xb1*
     . xb2*xs1*xs2-2*log(xb1)**2*xb1*xb2*xs1-2*log(xb1)**2*xb1*xb2*
     . xs2**3+4*log(xb1)**2*xb1*xb2*xs2**2-2*log(xb1)**2*xb1*xb2*xs2-
     . 2*log(xb1)**2*xb1*xq**3+log(xb1)**2*xb1*xq**2*xs1+5*log(xb1)**
     . 2*xb1*xq**2*xs2+ans11
      ans9=-4*t134(mb1,mq,ms2,mg)*xb1*xb2+2*t134(mb1,mq,ms2,mg)*xb1*
     . xq**2-4*t134(mb1,mq,ms2,mg)*xb1*xq*xs2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xq+2*t134(mb1,mq,ms2,mg)*xb1*xs2**2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xs2+2*t134(mb1,mq,ms2,mg)*xb1+2*t134(mb1,mq,ms2,mg)*xb2**2
     . *xq**3-6*t134(mb1,mq,ms2,mg)*xb2**2*xq**2*xs2-4*t134(mb1,mq,
     . ms2,mg)*xb2**2*xq**2+6*t134(mb1,mq,ms2,mg)*xb2**2*xq*xs2**2+2*
     . t134(mb1,mq,ms2,mg)*xb2**2*xq-2*t134(mb1,mq,ms2,mg)*xb2**2*xs2
     . **3+4*t134(mb1,mq,ms2,mg)*xb2**2*xs2**2-2*t134(mb1,mq,ms2,mg)*
     . xb2**2*xs2-4*t134(mb1,mq,ms2,mg)*xb2*xq**3+12*t134(mb1,mq,ms2,
     . mg)*xb2*xq**2*xs2+8*t134(mb1,mq,ms2,mg)*xb2*xq**2-12*t134(mb1,
     . mq,ms2,mg)*xb2*xq*xs2**2-4*t134(mb1,mq,ms2,mg)*xb2*xq+4*t134(
     . mb1,mq,ms2,mg)*xb2*xs2**3-8*t134(mb1,mq,ms2,mg)*xb2*xs2**2+4*
     . t134(mb1,mq,ms2,mg)*xb2*xs2+2*t134(mb1,mq,ms2,mg)*xq**3-6*t134
     . (mb1,mq,ms2,mg)*xq**2*xs2-4*t134(mb1,mq,ms2,mg)*xq**2+6*t134(
     . mb1,mq,ms2,mg)*xq*xs2**2+2*t134(mb1,mq,ms2,mg)*xq-2*t134(mb1,
     . mq,ms2,mg)*xs2**3+4*t134(mb1,mq,ms2,mg)*xs2**2-2*t134(mb1,mq,
     . ms2,mg)*xs2-2*log(xb1)**2*xb1*xb2**2*xq**3+log(xb1)**2*xb1*xb2
     . **2*xq**2*xs1+ans10
      ans8=4*t134(mb1,mq,ms1,mg)*xb2*xq**2*xs1+8*t134(mb1,mq,ms1,mg)*
     . xb2*xq**2*xs2+8*t134(mb1,mq,ms1,mg)*xb2*xq**2-8*t134(mb1,mq,
     . ms1,mg)*xb2*xq*xs1*xs2-8*t134(mb1,mq,ms1,mg)*xb2*xq*xs1-4*t134
     . (mb1,mq,ms1,mg)*xb2*xq*xs2**2+8*t134(mb1,mq,ms1,mg)*xb2*xq*xs2
     . -4*t134(mb1,mq,ms1,mg)*xb2*xq+4*t134(mb1,mq,ms1,mg)*xb2*xs1*
     . xs2**2-8*t134(mb1,mq,ms1,mg)*xb2*xs1*xs2+4*t134(mb1,mq,ms1,mg)
     . *xb2*xs1+2*t134(mb1,mq,ms1,mg)*xq**3-2*t134(mb1,mq,ms1,mg)*xq
     . **2*xs1-4*t134(mb1,mq,ms1,mg)*xq**2*xs2-4*t134(mb1,mq,ms1,mg)*
     . xq**2+4*t134(mb1,mq,ms1,mg)*xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*
     . xq*xs1+2*t134(mb1,mq,ms1,mg)*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xq*xs2+2*t134(mb1,mq,ms1,mg)*xq-2*t134(mb1,mq,ms1,mg)*xs1*xs2
     . **2+4*t134(mb1,mq,ms1,mg)*xs1*xs2-2*t134(mb1,mq,ms1,mg)*xs1+2*
     . t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2**2*xq*xs2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq+2*t134(mb1,
     . mq,ms2,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*
     . xs2+2*t134(mb1,mq,ms2,mg)*xb1*xb2**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2*xq**2+8*t134(mb1,mq,ms2,mg)*xb1*xb2*xq*xs2+8*t134(mb1,mq,
     . ms2,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms2,mg)*xb1*xb2*xs2**2+8*t134
     . (mb1,mq,ms2,mg)*xb1*xb2*xs2+ans9
      ans7=2*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms1,
     . mg)*xb1*xb2**2*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq+2*
     . t134(mb1,mq,ms1,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xb1*xb2**2*xs2+2*t134(mb1,mq,ms1,mg)*xb1*xb2**2-4*t134(mb1,mq,
     . ms1,mg)*xb1*xb2*xq**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xq*xs2+8*
     . t134(mb1,mq,ms1,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms1,mg)*xb1*xb2*
     . xs2**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xs2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xb2+2*t134(mb1,mq,ms1,mg)*xb1*xq**2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xq+2*t134(mb1,mq,ms1,mg)
     . *xb1*xs2**2-4*t134(mb1,mq,ms1,mg)*xb1*xs2+2*t134(mb1,mq,ms1,mg
     . )*xb1+2*t134(mb1,mq,ms1,mg)*xb2**2*xq**3-2*t134(mb1,mq,ms1,mg)
     . *xb2**2*xq**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2**2*xq**2*xs2-4*
     . t134(mb1,mq,ms1,mg)*xb2**2*xq**2+4*t134(mb1,mq,ms1,mg)*xb2**2*
     . xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs1+2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs2+2
     . *t134(mb1,mq,ms1,mg)*xb2**2*xq-2*t134(mb1,mq,ms1,mg)*xb2**2*
     . xs1*xs2**2+4*t134(mb1,mq,ms1,mg)*xb2**2*xs1*xs2-2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2*xq**3+ans8
      ans15=(xq**2-2*xq*xs1-2*xq+xs1**2-2*xs1+1)
      ans6=ans7*ans15
      ans5=-ans6
      ans4=(4*(2*(xs1-1)**2-xq**2-(xs1+1)*xq)-(xs1-1+xq)*(xs1-1-xq)*
     . log(xs1))*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*(xb1-1)*(
     . xb2-1)*log(xs1)*xs1+ans5
      ans3=-ans4
      ans17=4*(3*((xs1+1)*xs1+xs2**2+xs2-1)+2*(3*xq**2+2)*xq+((2*xq**
     . 2-5*xq+8)*xq-3*(xs1+xs2)**2)*(xs1+xs2)-(4*xq-3*xs1*xs2)*(2*xs1
     . **2+3*xs1*xs2+2*xs2**2)+2*(3*xs1**2-8*xs1*xs2+3*xs2**2)*(xs1+
     . xs2)*xq-(13*xq**2+4*xs1*xs2)*xq**2+((3*xs1**2-2*xs1*xs2+3*xs2
     . **2)*xq-2*(xs1+xs2)*xq**2-3*(xs1+xs2)*xs1*xs2+6*xq**3)*(xq-xs1
     . )*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)+4*(2*((xs2-1)**2+xs1**2+
     . 2*(2*xs2-1)*xs1+(5*xq**2-1)*xq)-(xs1**2-4*xs1*xs2+xs2**2)*xq-2
     . *(xq**2+2*xs1*xs2)*(xs1+xs2)+(xs1-1-xq)*(2*(xs2+1)*xq-(xs2-1)
     . **2-xq**2)*log(xs1)*xs1+(2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xs2-1
     . -xq)*log(xs2)*xs2+(xs1+xs2-6*xq)*xq+((xs1+xs2)*xq-4*xq**2+2*
     . xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq)*
     . xq
      ans16=-2*(2*((xs1-1)**2*xs1-xq**3-(3*xs1-1)*xq*xs1+(3*xs1+1)*xq
     . **2)+(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1)*xb2-(2*((xs1+1
     . )*xq-(xs1-1)**2)*xb2-(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1
     . ))*xb1)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*t134(mg,mq,
     . ms1,mg)-2*((xs2-1)**2+xs1**2+2*(2*xs2-1)*xs1+2*(xq+1)*(xq-1)*
     . xq+(xq-2*xs1*xs2)*(xs1+xs2)-(xs1**2-4*xs1*xs2+xs2**2)*xq-(xq**
     . 2-xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq
     . )**2*xq-2*((2*((xs2-1)**2*xs2-xq**3-(3*xs2-1)*xq*xs2+(3*xs2+1)
     . *xq**2)+(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-xs2-1)*xb2-(2*((
     . xs2+1)*xq-(xs2-1)**2)*xb2-(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-
     . xs2-1))*xb1)*(xb1-xb2)*t134(mg,mq,ms2,mg)+((xq-xs1+xb2)*t134(
     . mb2,mq,ms1,mg)+(xq-xs2+xb2)*t134(mb2,mq,ms2,mg))*(2*(xs2+1)*xq
     . -(xs2-1)**2-xq**2)*(xb1-1)**2)*(2*(xs1+1)*xq-(xs1-1)**2-xq**2)
     . +ans17
      ans2=2*(2*(xs1+xs2-2*xq)-log(xs2)*xs2-log(xs1)*xs1+2*log(xq)*xq
     . )*((xb1-1)**2*log(xb2)*xb2-(xb2-1)**2*log(xb1)*xb1)*(2*(xs1+1)
     . *xq-(xs1-1)**2-xq**2)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)+ans3+
     . ans16
      ans1=-ans2
      rt102=ans1/(8*((xs1-1)**2+xq**2-2*(xs1+1)*xq)*((xs2-1)**2+xq**2-
     . 2*(xs2+1)*xq)*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      rtalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rtmgtr2=((xs2-6+xs1-2*xq+(xq+1)*log(xs2)+(xq+1)*log(xs1)-2*log(
     . xq)*xq-(xs2-1-xq)*log(xp(ms2))*xp(ms2)-(xs1-1-xq)*log(xp(ms1))
     . *xp(ms1)-(xs2-1-xq)*log(xm(ms2))*xm(ms2)-(xs1-1-xq)*log(xm(ms1
     . ))*xm(ms1)+(xs2-1-xq)*log(xp(ms2)-1)*xp(ms2)+(xs1-1-xq)*log(xp
     . (ms1)-1)*xp(ms1)+(xs2-1-xq)*log(xm(ms2)-1)*xm(ms2)+(xs1-1-xq)*
     . log(xm(ms1)-1)*xm(ms1))*((xb1+1)*(xb2-1)**2*log(xb1)*xb1-(xb1-
     . 1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs1)+log(xs2))*(xq+1)-2*
     . log(xq)*xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*((log(xm(ms2)
     . -1)-log(xm(ms2)))*xm(ms2)+(log(xp(ms2)-1)-log(xp(ms2)))*xp(ms2
     . ))*(xs2-1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*((log(xm(ms1)-1)-log
     . (xm(ms1)))*xm(ms1)+(log(xp(ms1)-1)-log(xp(ms1)))*xp(ms1))*(xs1
     . -1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*(xs2-6+xs1-2*xq)*(xb1-xb2)*
     . (xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ans2=(((3*xs1-1)*(xs1-1)+3*xq**2-2*(3*xs1+2)*xq-(xq-xs1+1)*(xq-
     . xs1-1)*xb2-(((xs1-1)*(xs1-3)+xq**2-2*(xs1+2)*xq)*xb2+(xq-xs1+1
     . )*(xq-xs1-1))*xb1)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*t134(mg,mq,
     . ms1,mg)-2*(2*(2*xs2-1+2*xq-(xq-xs2)**2)*log(xs1)*xs1-3*((xs2+2
     . +xs1)*xq+(xs1-1)*(xs2-1)-3*xq**2)*(xs1-xs2)-2*(2*xs1-1+2*xq-(
     . xq-xs1)**2)*log(xs2)*xs2)*(xb1-1)*(xb2-1)*log(xq)*xq)*(xb1-xb2
     . )
      ans1=((2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xs2-1+xq)*(log(xs2)-6)*
     . log(xs2)*xs2-(xs1-1+xq)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(log(
     . xs1)-6)*log(xs1)*xs1)*(xb1-xb2)*(xb1-1)*(xb2-1)+((t134(mb1,mq,
     . ms1,mg)-t134(mb1,mq,ms2,mg))*(xb1+1)*(xb2-1)**2-(t134(mb2,mq,
     . ms1,mg)-t134(mb2,mq,ms2,mg))*(xb1-1)**2*(xb2+1))*(2*(xs1+1)*xq
     . -(xs1-1)**2-xq**2)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)-((xs2+2+xs1
     . )*xq+(xs1-1)*(xs2-1)-3*xq**2)*(xb1-xb2)*(xb1-1)*(xb2-1)*(xs1-
     . xs2)*log(xq)**2*xq+14*(xs2-1+xs1-(xq+xs1)*(xq+xs2)+2*(xq**2-
     . xs1*xs2)*xq)*(xb1-xb2)*(xb1-1)*(xb2-1)*(xs1-xs2)-((3*xs2-1)*(
     . xs2-1)+3*xq**2-2*(3*xs2+2)*xq-(xq-xs2+1)*(xq-xs2-1)*xb2-(((xs2
     . -1)*(xs2-3)+xq**2-2*(xs2+2)*xq)*xb2+(xq-xs2+1)*(xq-xs2-1))*xb1
     . )*(2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xb1-xb2)*t134(mg,mq,ms2,mg)
     . +ans2
      rtp102=ans1/(2*((xs1-1)**2+xq**2-2*(xs1+1)*xq)*((xs2-1)**2+xq**2-
     . 2*(xs2+1)*xq)*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2*(xs1-xs2))

      rmgtrt2=(-((log(xp(ms2)-1)-log(xp(ms2)))*xp(ms2)+log(xs1)-log(
     . xs2)-(log(xp(ms1)-1)-log(xp(ms1)))*xp(ms1)+(log(xm(ms2)-1)-log
     . (xm(ms2)))*xm(ms2)-(log(xm(ms1)-1)-log(xm(ms1)))*xm(ms1))*((
     . xb1+1)*(xb2-1)**2*log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2
     . +2*(xb1-xb2)*(xb1-1)*(xb2-1)))/(4*(xb1-xb2)*(xb1-1)**2*(xb2-1)
     . **2)

      fact = amt**2*xxt/amg**3
      rtp102 = fact * rtp102
      fact = 2*amt**2*xxt/(amst1**2-amst2**2)/amg
      rmgtrt2 = fact * rmgtrt2
      rtp102 = rtp102 + rmgtrt2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      rctca = ralsca2 + rmgca2
      rctcf = ralscf2 + rmb12 + rmb22
      ructtr = rualstr2 + rumgtr2
      rdcttr = rdalstr2 + rdmgtr2
      rbcttr = rbalstr2 + rbmgtr2
      rtcttr = rtalstr2 + rtmgtr2

      bo  =-2*fi(amsb1**2,amsb2**2,amg**2)

      rca = r12 + r22/4 + r32 + r52 + r62 + rctca
      rcf = -r22/2 - r42/2 - 2*r52 - r72/2 + r82 - r92/2 + rctcf
      rtr = nu*(ru102 + ructtr) + nd*(rd102 + rdcttr)
      rtrb= rb102 + rbcttr
      rtrt= rt102 + rtp102 + rtcttr + bo/3*dlog(amg**2/amt**2)/fnorm

      rqcd = ca*rca + cf*rcf + tr*(rtr+rtrb+rtrt)
      rqcd = rqcd/bo
      anomalous = - cf/4
      finscale = (11*ca-4*tr*nf)/12*dlog(scale**2/amg**2)
      fqcd_hdec = dreal(rqcd)*fnorm + anomalous + finscale

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     fac = fnorm/bo
c     ysu = tr*dreal(nu*(ru102+ructtr)*fac
c    .              -4*nu/12.d0*dlog(scale**2/amg**2))
c     ysd = tr*dreal(nd*(rd102+rdcttr)*fac
c    .              -4*nd/12.d0*dlog(scale**2/amg**2))
c     ysb = tr*dreal(rtrb*fac-4/12.d0*dlog(scale**2/amg**2))
c     yst = tr*dreal(rtrt-rtp102)*fac
c     ystp= tr*dreal(rtp102/xxt)*fac
c     fac = 0.7705157934627561d0*0.088374016132112168d0/pi
c     ysu  = ysu * fac
c     ysd  = ysd * fac
c     ysb  = ysb * fac
c     yst  = yst * fac
c     ystp = ystp* fac
c     write(6,*)'dmb: '
c     write(6,*)'full ',dreal(rqcd)*fnorm*fac,anomalous*fac,finscale*fac
c     write(6,*)'CA = ',ca*(dreal(rca)*fnorm/bo
c    .                     +11/12.d0*dlog(scale**2/amg**2))*fac
c     write(6,*)'CF = ',cf*(dreal(rcf)*fnorm/bo-1/4.d0)*fac
c     write(6,*)'TR = ',tr*(dreal(rtr+rtrb+rtrt)*fnorm/bo
c    .                     -4*nf/12.d0*dlog(scale**2/amg**2))*fac
c     write(6,*)'su,sd,sb,st: ',ysu,ysd,ysb,yst,ystp*xxt,ystp
c     write(6,*)'xt = ',xxt
c     write(6,*)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      return
      end
 
      COMPLEX*16 FUNCTION LI2_HDEC(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CLI2_HDEC
      COMMON/CONST_HDEC/ZETA2,ZETA3
      ZERO=1.D-16
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      LI2_HDEC=0
      IF(R2.LE.ZERO)THEN
        LI2_HDEC=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI2_HDEC=DCMPLX(ZETA2)
        ELSE
          LI2_HDEC=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
        LI2_HDEC=CLI2_HDEC(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
     .          +0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        LI2_HDEC=-CLI2_HDEC(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        LI2_HDEC=-CLI2_HDEC(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        LI2_HDEC=CLI2_HDEC(Y)
        RETURN
      ENDIF
      END


      COMPLEX*16 FUNCTION CLI2_HDEC(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI_HDEC/B2(18),B12(18),B3(18)
      COMMON/POLY_HDEC/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI2_HDEC=B2(NBER)
      DO 111 I=N,1,-1
        CLI2_HDEC=Z*CLI2_HDEC+B2(I)
111   CONTINUE
      CLI2_HDEC=Z**2*CLI2_HDEC+Z
      RETURN
      END
