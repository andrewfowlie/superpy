      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*  Input data corresponds to GMSB model. In addition:
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*
*      The input file contains lower and upper bounds
*      as well as number of steps for these paramaters
*
*  On output:
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar at the scale
*      Q2 which is either user defined or computed as (2*mQ2+mU2+mD2)/4
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      PMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
*        the convention is 
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(5)=0
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> gluon gluon
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRSS(i)        "       -> ss
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino 
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  ERRORS: IFAIL = 0..19
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         l, tan(beta) or mu = 0
*          10        Violation of phenomenological constraint(s)
*          11,12,13  Problem in integration of RGEs
*          14,15,16  Convergence problem
*          17        No electroweak symmetry breaking
*          18        MSMES=/=MSREF (scenario 2 of ArXiv:0803.2962)
*          19        IFAIL = 10 & 18
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..50: OK
*            
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by t -> bH+ (LHC)
*      PROB(46) =/= 0  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0  chi2gam > chi2max
*      PROB(48) =/= 0  chi2bb > chi2max
*      PROB(49) =/= 0  chi2zz > chi2max
*      PROB(50) =/= 0  excluded by sparticle searches at the LHC
*
************************************************************************

      IMPLICIT NONE

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=19,NPROB=50,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDIM
      INTEGER IMSUSYEFF,IMMESS,ITB,IL,IAL,IMUP,IMSP,IDH
      INTEGER NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,NDH
      INTEGER N1,N2,I1,I2,ITER,ITERMU,Q2FIX
      INTEGER OMGFLAG,MAFLAG,GMUFLAG,HFLAG,ALFLAG

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,MESTEST,DELMB
      DOUBLE PRECISION MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX
      DOUBLE PRECISION TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX
      DOUBLE PRECISION XIFMIN,XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX
      DOUBLE PRECISION MSPMIN,MSPMAX,MSMIN,MSMAX,DHMIN,DHMAX
      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,DHN,DHNN,MUN,MUNN
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DHINP
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION DETM,MUFAIL,SIGMU,Q2,Q2MIN
      DOUBLE PRECISION PI,ALP1,ALP2
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTOPMES
      DOUBLE PRECISION HBOTMES,HTAUMES,MSREF,D,DMIN

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/SCEN2/MSREF,D,DMIN,ALFLAG
      COMMON/SIGMU/SIGMU
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/DELMB/DELMB
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     . HBOTMES,HTAUMES
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DHINP
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/MINMAX/MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX,
     . TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,XIFMIN,
     . XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX,
     . MSMIN,MSMAX,DHMIN,DHMAX
      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,DHN,DHNN,MUN,MUNN
      COMMON/STEPS/NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,
     . NDH,N1,N2
      COMMON/IDIM/IDIM

      PI=4d0*DATAN(1d0)

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      ITOT=0
      IDIM=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      MSUSYEFFN=1d99
      MSUSYEFFNN=-1d99
      MMESSN=1d99
      MMESSNN=-1d99
      TBN=1d99
      TBNN=-1d99
      LN=1d99
      LNN=-1d99
      ALN=1d99
      ALNN=-1d99
      MUN=1d99
      MUNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      KN=1d99
      KNN=-1d99
      MSN=1d99
      MSNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      DHN=1d99
      DHNN=-1d99
      
*   Beginning of the scan

      DO IMSUSYEFF=1,NMSUSYEFF
      IF(NMSUSYEFF.EQ.1)THEN
       MSUSYEFF=MSUSYEFFMIN
      ELSE
       MSUSYEFF=MSUSYEFFMIN+(MSUSYEFFMAX-MSUSYEFFMIN)
     .        *DFLOAT(IMSUSYEFF-1)/DFLOAT(NMSUSYEFF-1)
      ENDIF

      DO IMMESS=1,NMMESS
      IF(NMMESS.EQ.1)THEN
       MMESS=MMESSMIN
      ELSE
       MMESS=MMESSMIN+(MMESSMAX-MMESSMIN)
     .     *DFLOAT(IMMESS-1)/DFLOAT(NMMESS-1)
      ENDIF

      DO ITB=1,NTB
      IF(NTB.EQ.1)THEN
       PAR(3)=TBMIN
      ELSE
       PAR(3)=TBMIN+(TBMAX-TBMIN)*DFLOAT(ITB-1)/DFLOAT(NTB-1)
      ENDIF

      DO IL=1,NL
      IF(NL.EQ.1)THEN
       PAR(1)=LMIN
      ELSE
       PAR(1)=LMIN+(LMAX-LMIN)*DFLOAT(IL-1)/DFLOAT(NL-1)
      ENDIF

      DO I1=1,N1
      IF(N1.EQ.1)THEN
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFMIN
       ELSE
        PAR(2)=KMIN
       ENDIF
      ELSE
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFMIN+(XIFMAX-XIFMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ELSE
        PAR(2)=KMIN+(KMAX-KMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ENDIF
      ENDIF

      DO I2=1,N2
      IF(N2.EQ.1)THEN
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISMIN
       ELSE
        MSINP=MSMIN
       ENDIF
      ELSE
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISMIN+(XISMAX-XISMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ELSE
        MSINP=MSMIN+(MSMAX-MSMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ENDIF
      ENDIF

      DO IMUP=1,NMUP
      IF(NMUP.EQ.1)THEN
       MUPINP=MUPMIN
      ELSE
       MUPINP=MUPMIN+(MUPMAX-MUPMIN)*DFLOAT(IMUP-1)/DFLOAT(NMUP-1)
      ENDIF

      DO IMSP=1,NMSP
      IF(NMSP.EQ.1)THEN
       MSPINP=MSPMIN
      ELSE
       MSPINP=MSPMIN+(MSPMAX-MSPMIN)*DFLOAT(IMSP-1)/DFLOAT(NMSP-1)
      ENDIF

      DO IDH=1,NDH
      IF(NDH.EQ.1)THEN
       DHINP=DHMIN
      ELSE
       DHINP=DHMIN+(DHMAX-DHMIN)*DFLOAT(IDH-1)/DFLOAT(NDH-1)
      ENDIF

      DO IAL=1,NAL
      IF(NAL.EQ.1)THEN
       IF(ALFLAG.EQ.1)THEN
        ALINP=-3d0/(16d0*PI**2)*DHINP*MSUSYEFF
       ELSE
        ALINP=ALMIN
       ENDIF
      ELSE
       ALINP=ALMIN+(ALMAX-ALMIN)*DFLOAT(IAL-1)/DFLOAT(NAL-1)
      ENDIF

      ITOT=ITOT+1

!      WRITE(0,*)""
!      WRITE(0,*)"------------------------------------------------------"
!      WRITE(0,*)""
!      WRITE(0,*)"Point ",ITOT
!      WRITE(0,*)""
!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)"MSUSYEFF =",MSUSYEFF
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"TANB =",PAR(3)
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      WRITE(0,*)"ALAMBDA=",ALINP
!      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)WRITE(0,*)"KAPPA =",PAR(2)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)WRITE(0,*)"XIF =",XIFINP
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)WRITE(0,*)"XIS =",XISINP
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)WRITE(0,*)"MS =",MSINP
!      WRITE(0,*)"MUP =",MUPINP
!      WRITE(0,*)"MSP =",MSPINP
!      WRITE(0,*)"DH =",DHINP
!      WRITE(0,*)""

*   Initialization of PROB

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO      

*   Guess parameters at Q2/MMESS

      MUFAIL=1d0
      ALP1=G1**2/(4d0*PI)
      ALP2=G2**2/(4d0*PI)
      IF(Q2FIX.EQ.0)THEN
       Q2=MAX(N5*MSUSYEFF/(4d0*PI),Q2MIN)
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       PAR(2)=PAR(1)/5d0
      ELSE
       XIFMES=0d0
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       MSMES=0d0
      ELSE
       XISMES=0d0
      ENDIF
      PAR(4)=SIGMU*DSQRT(Q2)
      PAR(5)=ALINP
      IF(PAR(2).NE.0d0)THEN
       PAR(6)=3d0*ALINP
      ELSE
       PAR(6)=0d0
      ENDIF
      PAR(7)=5d0/54d0*ALP1**2+3d0/2d0*ALP2**2+8d0/3d0*ALSMZ**2
      PAR(8)=40d0/27d0*ALP1**2+8d0/3d0*ALSMZ**2
      PAR(9)=10d0/27d0*ALP1**2+8d0/3d0*ALSMZ**2
      PAR(10)=5d0/6d0*ALP1**2+3d0/2d0*ALP2**2
      PAR(11)=10d0/3d0*ALP1**2
      DO I=12,14
       PAR(I)=0d0
      ENDDO
      PAR(15)=5d0/54d0*ALP1**2+3d0/2d0*ALP2**2+8d0/3d0*ALSMZ**2
      PAR(16)=40d0/27d0*ALP1**2+8d0/3d0*ALSMZ**2
      PAR(17)=10d0/27d0*ALP1**2+8d0/3d0*ALSMZ**2
      PAR(18)=5d0/6d0*ALP1**2+3d0/2d0*ALP2**2
      PAR(19)=10d0/3d0*ALP1**2
      DO I=7,19
       PAR(I)=PAR(I)*N5*(MSUSYEFF/(4d0*PI))**2
      ENDDO
      PAR(20)=5d0/3d0*ALP1
      PAR(21)=ALP2
      PAR(22)=ALSMZ
      DO I=20,22
       PAR(I)=PAR(I)*N5*MSUSYEFF/(4d0*PI)
      ENDDO
      PAR(23)=DSQRT(Q2)
      PAR(24)=DSQRT(Q2)
      PAR(25)=0d0
      DELMB=.1d0
      IFAIL=0

*   Check for singular parameters l, tan(beta)

      IF(PAR(1)*PAR(3).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Guess for couplings at MMESS

      CALL RGESGM(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   External loop to compute the soft parameters at Q2

      ITER=0
 21   ITER=ITER+1
      !WRITE(0,*)"ITER =",ITER
      !WRITE(0,*)""
      !WRITE(0,*)""

      CALL RGESINVGM(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   Internal loop to compute mu, k and ms

      ITERMU=0
 22   ITERMU=ITERMU+1
      !WRITE(0,*)"ITERMU =",ITERMU
      !WRITE(0,*)""
      !WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)GOTO 11

      CALL MINIMIZE(PAR,CHECK)
      IF(DETM.LT.0)IFAIL=17
      IF(IFAIL.NE.0)GOTO 11

      IF(CHECK.GT.1.D-12.AND.ITERMU.LE.100) GOTO 22
      IF(CHECK.GT.1.D-8) IFAIL=14
      IF(IFAIL.NE.0) GOTO 11
      
      CALL RGESGM(PAR,IFAIL)
      IF(PROB(27).NE.0.)IFAIL=11
      IF(IFAIL.NE.0) GOTO 11

      CALL RGESUNIGM(PAR,IFAIL,MESTEST)
      IF(IFAIL.NE.0)GOTO 11

      IF(MESTEST.GT.1.D-12.AND.ITER.LE.100) GOTO 21
      IF(MESTEST.GT.1.D-8) IFAIL=15
      IF(IFAIL.NE.0) GOTO 11

*   Computation of sfermion masses   

      CALL MSFERM(PAR,IFAIL,1)
      IF(IFAIL.NE.0)GOTO 11
      
*   Computation of Higgs masses

      CALL MHIGGS(PAR,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=16
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)      
      
*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)
      
*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR)
      CALL TDECAY(PAR)

*   Exp. constraints

      CALL SUBEXP(PAR,PROB)

*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)
      
*   Anom. magn. moment of the Muon

      IF(GMUFLAG.EQ.1)CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*  GUT scale

c      CALL RGESGMGUT(PROB,IFAIL)

*   Computation of the fine-tuning

      CALL FTPAR(PAR,2)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)IFAIL=10
      ENDDO

*   Check MSMES=MSREF if ALFLAG=1  

      IF(ALFLAG.EQ.1)THEN
       MSREF=(7d0/5d0*DHINP**2-(16d0*G3MES+6d0*G2MES+10d0/3d0*G1MES)
     .       *DHINP/5d0-4d0*KMES**2*DHINP)*MSUSYEFF**2/(16d0*PI**2)**2
       D=2d0*DABS(MSMES-MSREF)/DABS(MSMES+MSREF)
       IF(D.GT.DMIN)THEN
        IF(IFAIL.EQ.0)IFAIL=18
        IF(IFAIL.EQ.10)IFAIL=19
       ENDIF
!       WRITE(0,*)MSREF,MSMES,D
      ENDIF

*   Recording of the results

11    CALL OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT)

      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       MSUSYEFFN=MIN(MSUSYEFF,MSUSYEFFN)
       MSUSYEFFNN=MAX(MSUSYEFF,MSUSYEFFNN)
       MMESSN=MIN(MMESS,MMESSN)
       MMESSNN=MAX(MMESS,MMESSNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       ALN=MIN(ALINP,ALN)
       ALNN=MAX(ALINP,ALNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       XIFN=MIN(XIFMES,XIFN)
       XIFNN=MAX(XIFMES,XIFNN)
       XISN=MIN(XISMES,XISN)
       XISNN=MAX(XISMES,XISNN)
       MSN=MIN(MSMES,MSN)
       MSNN=MAX(MSMES,MSNN)
       MUPN=MIN(MUPMES,MUPN)
       MUPNN=MAX(MUPMES,MUPNN)
       MSPN=MIN(MSPMES,MSPN)
       MSPNN=MAX(MSPMES,MSPNN)
       DHN=MIN(DHINP,DHN)
       DHNN=MAX(DHINP,DHNN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

*   Summary of the scan:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

      CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)
      
*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX
      INTEGER N0,NLOOP,NBER,NPAR,ERR
      INTEGER OMGFLAG,MAFLAG,PFLAG,GMUFLAG,HFLAG,ALFLAG
      INTEGER NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,NDH
      INTEGER N1,N2,NK,NXIF,NXIS,NMS

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
      DOUBLE PRECISION MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX
      DOUBLE PRECISION TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX
      DOUBLE PRECISION XIFMIN,XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX
      DOUBLE PRECISION MSPMIN,MSPMAX,MSMIN,MSMAX,DHMIN,DHMAX
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MSREF,D,DMIN

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MINMAX/MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX,
     . TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,XIFMIN,
     . XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX,
     . MSMIN,MSMAX,DHMIN,DHMAX
      COMMON/STEPS/NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,
     . NDH,N1,N2
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/PFLAG/PFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/SCEN2/MSREF,D,DMIN,ALFLAG

*   MAXIMUM DEV (%) FOR MSMES/MSREF (SCENARIO 2)
      DMIN=1d-1

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO      

*   INITIALIZATION OF THE SCANNING PARAMETERS
      N5=1d99
      SIGMU=1d99
      MSUSYEFFMIN=1d99
      MSUSYEFFMAX=1d99
      MMESSMIN=1d99
      MMESSMAX=1d99
      TBMIN=1d99
      TBMAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=1d99
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      XIFMIN=1d99
      XIFMAX=1d99
      XISMIN=1d99
      XISMAX=1d99
      MUPMIN=0d0
      MUPMAX=1d99
      MSPMIN=0d0
      MSPMAX=1d99
      MSMIN=1d99
      MSMAX=1d99
      DHMIN=0d0
      DHMAX=1d99
      NMSUSYEFF=1
      NMMESS=1
      NTB=1
      NL=1
      NK=0
      NAL=1
      NXIF=0
      NXIS=0
      NMUP=1
      NMSP=1
      NMS=0
      NDH=1

*   DEFAULT VALUES FOR FLAGS
      OMGFLAG=0
      PFLAG=0
      GMUFLAG=1
      HFLAG=0
      ALFLAG=0

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(5,'(A120)',END=29,ERR=999) CHINL
      
*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       END DO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.12) ALFLAG=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL
 
*   READ MESS PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.17) MSUSYEFFMIN=VAL
       IF(IX.EQ.18) MSUSYEFFMAX=VAL
       IF(IX.EQ.27) MMESSMIN=VAL
       IF(IX.EQ.28) MMESSMAX=VAL
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL
       IF(IX.EQ.4)  SIGMU=VAL
       IF(IX.EQ.5)  N5=VAL
 
*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.668) XIFMAX=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.678) XISMAX=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.688) MUPMAX=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.698) MSPMAX=VAL
       IF(IX.EQ.707) MSMIN=VAL
       IF(IX.EQ.708) MSMAX=VAL
       IF(IX.EQ.717) DHMIN=VAL
       IF(IX.EQ.718) DHMAX=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.19) NMSUSYEFF=IVAL
       IF(IX.EQ.29) NMMESS=IVAL
       IF(IX.EQ.39) NTB=IVAL
       IF(IX.EQ.619) NL=IVAL
       IF(IX.EQ.629) NK=IVAL
       IF(IX.EQ.639) NAL=IVAL
       IF(IX.EQ.669) NXIF=IVAL
       IF(IX.EQ.679) NXIS=IVAL
       IF(IX.EQ.689) NMUP=IVAL
       IF(IX.EQ.699) NMSP=IVAL
       IF(IX.EQ.709) NMS=IVAL
       IF(IX.EQ.719) NDH=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(ALFLAG.NE.0 .AND. ALFLAG.NE.1)THEN
       WRITE(0,1)"ALFLAG IS EITHER 0 OR 1"
       ERR=1
      ENDIF
      IF(MSUSYEFFMIN.EQ.1d99)THEN
       WRITE(0,1)"MSUSYEFF_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(MMESSMIN.EQ.1d99)THEN
       WRITE(0,1)"MMESS_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANB_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.EQ.1d99)THEN
       WRITE(0,1)"N5 MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.LE.0d0.OR. ANINT(N5).NE.N5)THEN
       WRITE(0,1)"N5 MUST BE A POSITIVE INTEGER IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"LAMBDA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"ALAMBDA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KMIN.EQ.1d99 .AND. KMAX.NE.1d99)THEN
       WRITE(0,1)"KAPPA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.EQ.1d99 .AND. XIFMAX.NE.1d99)THEN
       WRITE(0,1)"XIF_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.EQ.1d99 .AND. XISMAX.NE.1d99)THEN
       WRITE(0,1)"XIS_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSMIN.EQ.1d99 .AND. MSMAX.NE.1d99)THEN
       WRITE(0,1)"MS^2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUPMIN.EQ.1d99 .AND. MUPMAX.NE.1d99)THEN
       WRITE(0,1)"MU'_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSPMIN.EQ.1d99 .AND. MSPMAX.NE.1d99)THEN
       WRITE(0,1)"MS'^2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(DHMIN.EQ.1d99 .AND. DHMAX.NE.1d99)THEN
       WRITE(0,1)"DH_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.NE.1d99 .AND. KMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIF AND KAPPA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.NE.1d99 .AND. MSMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIS AND MS^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALFLAG.EQ.1 .AND. ALMIN.NE.1d99)THEN
       WRITE(0,1)"ALFLAG=1 => ALMIN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALFLAG.EQ.1 .AND. MSMIN.NE.1d99)THEN
       WRITE(0,1)"ALFLAG=1 => MSMIN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

*   Set default values

      IF(KMIN.EQ.1d99 .AND. XIFMIN.EQ.1d99)XIFMIN=0d0
      IF(MSMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)XISMIN=0d0
      IF(ALMIN.EQ.1d99)ALMIN=0d0

      IF(MSUSYEFFMAX.EQ.1d99)MSUSYEFFMAX=MSUSYEFFMIN
      IF(MMESSMAX.EQ.1d99)MMESSMAX=MMESSMIN
      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(XIFMAX.EQ.1d99)XIFMAX=XIFMIN
      IF(XISMAX.EQ.1d99)XISMAX=XISMIN
      IF(MUPMAX.EQ.1d99)MUPMAX=MUPMIN
      IF(MSPMAX.EQ.1d99)MSPMAX=MSPMIN
      IF(MSMAX.EQ.1d99)MSMAX=MSMIN
      IF(DHMAX.EQ.1d99)DHMAX=DHMIN

*   Set MAFLAG

      IF(KMIN.EQ.1d99 .AND. MSMIN.EQ.1d99)MAFLAG=-1
      IF(KMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)MAFLAG=-2
      IF(XIFMIN.EQ.1d99 .AND. MSMIN.EQ.1d99)MAFLAG=-3
      IF(XIFMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)MAFLAG=-4

*   Number of points

      IF(MAFLAG.EQ.-1)THEN
       IF(NK.NE.0)THEN
        WRITE(0,1)"NK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.0)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIF.EQ.0)NXIF=1
       IF(NXIS.EQ.0)NXIS=1
       N1=NXIF
       N2=NXIS
      ENDIF

      IF(MAFLAG.EQ.-2)THEN
       IF(NK.NE.0)THEN
        WRITE(0,1)"NK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIF.EQ.0)NXIF=1
       IF(NMS.EQ.0)NMS=1
       N1=NXIF
       N2=NMS
      ENDIF

      IF(MAFLAG.EQ.-3)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.0)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NK.EQ.0)NK=1
       IF(NXIS.EQ.0)NXIS=1
       N1=NK
       N2=NXIS
      ENDIF

      IF(MAFLAG.EQ.-4)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NK.EQ.0)NK=1
       IF(NMS.EQ.0)NMS=1
       N1=NK
       N2=NMS
      ENDIF

      IF(ALFLAG.EQ.1)THEN
       IF(NAL.NE.1)THEN
        WRITE(0,1)"NAL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.1)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(NMSUSYEFF.LE.0 .OR. NMMESS.LE.0 .OR. NTB.LE.0 .OR.
     . NL.LE.0.OR. NAL.LE.0 .OR. NMUP.LE.0 .OR. NMSP.LE.0
     . .OR. NDH.LE.0.OR. N1.LE.0 .OR. N2.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Warnings

      IF(MSUSYEFFMIN.EQ.MSUSYEFFMAX .AND. NMSUSYEFF.GT.1)THEN
       WRITE(0,1)"WARNING MSUSYEFF_min=MSUSYEFF_max => NMSUSYEFF=1"
       NMSUSYEFF=1
      ENDIF
      IF(MMESSMIN.EQ.MMESSMAX .AND. NMMESS.GT.1)THEN
       WRITE(0,1)"WARNING MMESS_min=MMESS_max => NMMESS=1"
       NMMESS=1
      ENDIF
      IF(TBMIN.EQ.TBMAX .AND. NTB.GT.1)THEN
       WRITE(0,1)"WARNING TANB_min=TANB_max => NTB=1"
       NTB=1
      ENDIF
      IF(LMIN.EQ.LMAX .AND. NL.GT.1)THEN
       WRITE(0,1)"WARNING LAMBDA_min=LAMBDA_max => NL=1"
       NL=1
      ENDIF
      IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
       WRITE(0,1)"WARNING ALAMBDA_min=ALAMBDA_max => NAL=1"
       NAL=1
      ENDIF
      IF(MUPMIN.EQ.MUPMAX .AND. NMUP.GT.1)THEN
       WRITE(0,1)"WARNING MU'_min=MU'_max => NMUP=1"
       NMUP=1
      ENDIF
      IF(MSPMIN.EQ.MSPMAX .AND. NMSP.GT.1)THEN
       WRITE(0,1)"WARNING MS'^2_min=MS'^2_max => NMSP=1"
       NMSP=1
      ENDIF
      IF(DHMIN.EQ.DHMAX .AND. NDH.GT.1)THEN
       WRITE(0,1)"WARNING DH_min=DH_max => NDH=1"
       NDH=1
      ENDIF

      IF(MSUSYEFFMIN.NE.MSUSYEFFMAX .AND. NMSUSYEFF.EQ.1)THEN
       WRITE(0,10)"WARNING NMSUSYEFF=1 => MSUSYEFF_max=MSUSYEFF_min=",
     . MSUSYEFFMIN
       MSUSYEFFMAX=MSUSYEFFMIN
      ENDIF
      IF(MMESSMIN.NE.MMESSMAX .AND. NMMESS.EQ.1)THEN
       WRITE(0,10)"WARNING NMMESS=1 => MMESS_max=MMESS_min=",MMESSMIN
       MMESSMAX=MMESSMIN
      ENDIF
      IF(TBMIN.NE.TBMAX .AND. NTB.EQ.1)THEN
       WRITE(0,10)"WARNING NB=1 => TANB_max=TANB_min=",TBMIN
       TBMAX=TBMIN
      ENDIF
      IF(LMIN.NE.LMAX .AND. NL.EQ.1)THEN
       WRITE(0,10)"WARNING NL=1 => LAMBDA_max=LAMBDA_min=",LMIN
       LMAX=LMIN
      ENDIF
      IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
       WRITE(0,10)"WARNING NAL=1 => ALAMBDA_max=ALAMBDA_min=",ALMIN
       ALMAX=ALMIN
      ENDIF
      IF(MUPMIN.NE.MUPMAX .AND. NMUP.EQ.1)THEN
       WRITE(0,10)"WARNING NMUP=1 => MU'_max=MU'_min=",MUPMIN
       MUPMAX=MUPMIN
      ENDIF
      IF(MSPMIN.NE.MSPMAX .AND. NMSP.EQ.1)THEN
       WRITE(0,10)"WARNING NMSP=1 => MS'^2_max=MS'^2_min=",MSPMIN
       MSPMAX=MSPMIN
      ENDIF
      IF(DHMIN.NE.DHMAX .AND. NDH.EQ.1)THEN
       WRITE(0,10)"WARNING NDH=1 => DH_max=DH_min=",DHMIN
       DHMAX=DHMIN
      ENDIF

      IF(MAFLAG.EQ.-1)THEN
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIF_min=XIF_max => NXIF=1"
        NXIF=1
        N1=NXIF
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XIS_min=XIS_max => NXIS=1"
        NXIS=1
        N2=NXIS
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIF_max=XIF_min=", XIFMIN
        XIFMAX=XIFMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XIS_max=XIS_min=", XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-2)THEN
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIF_min=XIF_max => NXIF=1"
        NXIF=1
        N1=NXIF
       ENDIF
       IF(MSMIN.EQ.MSMAX .AND. NMS.GT.1)THEN
        WRITE(0,1)"WARNING MS^2_min=MS^2_max => NMS=1"
        NMS=1
        N2=NMS
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIF_max=XIF_min=", XIFMIN
        XIFMAX=XIFMIN
       ENDIF
       IF(MSMIN.NE.MSMAX .AND. NMS.EQ.1)THEN
        WRITE(0,10)"WARNING NMS=1 => MS^2_max=MS^2_min=", MSMIN
        MSMAX=MSMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-3)THEN
       IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
        WRITE(0,1)"WARNING KAPPA_min=KAPPA_max => NK=1"
        NK=1
        N1=NK
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XIS_min=XIS_max => NXIS=1"
        NXIS=1
        N2=NXIS
       ENDIF
       IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
        WRITE(0,10)"WARNING NK=1 => KAPPA_max=KAPPA_min=", KMIN
        KMAX=KMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XIS_max=XIS_min=", XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-4)THEN
       IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
        WRITE(0,1)"WARNING KAPPA_min=KAPPA_max => NK=1"
        NK=1
        N1=NK
       ENDIF
       IF(MSMIN.EQ.MSMAX .AND. NMS.GT.1)THEN
        WRITE(0,1)"WARNING MS^2_min=MS^2_max => NMS=1"
        NMS=1
        N2=NMS
       ENDIF
       IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
        WRITE(0,10)"WARNING NK=1 => KAPPA_max=KAPPA_min=", KMIN
        KMAX=KMIN
       ENDIF
       IF(MSMIN.NE.MSMAX .AND. NMS.EQ.1)THEN
        WRITE(0,10)"WARNING NMS=1 => MS^2_max=MS^2_min=", MSMIN
        MSMAX=MSMIN
       ENDIF
      ENDIF

*   total number of points

      NTOT=NMSUSYEFF*NMMESS*NTB*NL*NAL*NMUP*NMSP*NDH*N1*N2

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=2
      NBER=18
      ACC=1.D-8
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of 
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(A)
 10   FORMAT(A,E10.3)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT) 

*********************************************************************      
*   Subroutine writing all the results in the the output file.
*********************************************************************      
 
      IMPLICIT NONE 
 
      INTEGER IFAIL,OMGFLAG,MAFLAG,ITOT,NTOT,IDIM,I,J
      INTEGER NRES,NSUSY,NGUT,NMES,IMAX,JMAX
      PARAMETER (NSUSY=14,NGUT=21,NMES=15,IMAX=1000,JMAX=100)
 
      DOUBLE PRECISION PI,RES(IMAX,JMAX),PAR(*),PROB(*),SIG(3,8),S
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,Q2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     . HBOTMES,HTAUMES      
      DOUBLE PRECISION SIGMU
      DOUBLE PRECISION MHUQ,MHDQ,MSQ,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DHINP
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT,MGUT
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION DELMB,PX,PA(6),PB(2),PL(7),PK(8),MH(3),MMH(3)
      DOUBLE PRECISION DMH(3),MA(2),MMA(2),DMA(2),MHC,MMHC,DMHC
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DHINP
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SOFTMES/M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      COMMON/SIGMU/SIGMU
      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     . HBOTMES,HTAUMES
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT      
      COMMON/MGUT/MGUT
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFHIGM/MH,MMH,DMH,MA,MMA,DMA,MHC,MMHC,DMHC
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/IDIM/IDIM
      COMMON/RES/RES

      PI=4d0*DATAN(1d0)
      NRES=33

      IDIM=IDIM+1
      RES(IDIM,1)=PAR(1)
      RES(IDIM,2)=PAR(2)
      RES(IDIM,3)=PAR(3)
      RES(IDIM,4)=MMESS
      RES(IDIM,5)=MSUSYEFF
      RES(IDIM,6)=ALINP
      RES(IDIM,7)=MSMES
      RES(IDIM,8)=XIFMES
      RES(IDIM,9)=XISMES
      RES(IDIM,10)=MUPINP
      RES(IDIM,11)=MSPINP
      RES(IDIM,12)=DHINP

      IF(IFAIL.EQ.0.OR.IFAIL.EQ.10.OR.IFAIL.EQ.18.OR.IFAIL.EQ.19)THEN

       DO I=1,3
        RES(IDIM,11+2*I)=SMASS(I)
        RES(IDIM,12+2*I)=SCOMP(I,3)**2
       ENDDO
       DO I=1,2
        RES(IDIM,17+2*I)=PMASS(I)
        RES(IDIM,18+2*I)=PCOMP(I,2)**2
       ENDDO
       RES(IDIM,23)=CMASS
       RES(IDIM,24)=DABS(MNEU(1))
       RES(IDIM,25)=NEU(1,1)**2
       RES(IDIM,26)=NEU(1,3)**2+NEU(1,4)**2
       RES(IDIM,27)=NEU(1,5)**2

*   SUSY
       S=0d0
       DO I=1,2
        S=S+PROB(I)
       ENDDO
       DO I=20,25
        S=S+PROB(I)
       ENDDO
       RES(IDIM,28)=S

*   Higgs
       S=0d0
       DO I=3,19
        S=S+PROB(I)
       ENDDO
       S=S+PROB(41)
       RES(IDIM,29)=S

*   B Physics
       S=0d0
       DO I=32,36
        S=S+DABS(PROB(I))
       ENDDO
       DO I=38,40
        S=S+DABS(PROB(I))
       ENDDO
       RES(IDIM,30)=S

*   (g-2)_muon
       RES(IDIM,31)=PROB(37)

*   Theoretical
       S=0d0
       DO I=27,29
        S=S+PROB(I)
       ENDDO
       RES(IDIM,32)=S

*   Charged LSP
       RES(IDIM,33)=PROB(26)

      ELSE

       DO I=13,33
        RES(IDIM,I)=0d0
       ENDDO

      ENDIF

      IF(ITOT.EQ.NTOT.OR.IDIM.EQ.IMAX)THEN
       DO I=1,IDIM
        WRITE(6,10)(RES(I,J),J=1,NRES)
 10     FORMAT(100E16.8)
       ENDDO
       IF(IDIM.EQ.IMAX)IDIM=0
      ENDIF

      END


      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************      
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************      
 
      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG,ALFLAG,I,S,TOT,NTOT,IDIM,NFAIL(*)

      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,DHN,DHNN,MUN,MUNN
      DOUBLE PRECISION MSREF,D,DMIN

      COMMON/IDIM/IDIM
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,DHN,DHNN,MUN,MUNN
      COMMON/SCEN2/MSREF,D,DMIN,ALFLAG

      WRITE(0,*)
      WRITE(0,20)"Number of points:                       "
      WRITE(0,*)
      WRITE(0,20)"  scanned                               ",NTOT
      WRITE(0,20)"  l, tan(beta) or mu=0                  ",NFAIL(9)
      WRITE(0,20)"  no electroweak symmetry breaking      ",NFAIL(17)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,20)"  with mh1^2 or ma1^2 or mhc^2 < 0      ",S
      WRITE(0,20)"  with m_sfermion^2 < 0                 ",NFAIL(8)
      WRITE(0,20)"  violating phenomenological constraints",NFAIL(10)
      IF(ALFLAG.EQ.1)THEN
       WRITE(0,20)"  MSMES=/=MSREF                         ",NFAIL(18)
       WRITE(0,20)"  MSMES=/=MSREF + violating constraints ",NFAIL(19)
      ENDIF
      S=NFAIL(11)+NFAIL(12)+NFAIL(13)
      WRITE(0,20)"  RGE integration problem               ",S
      S=NFAIL(14)+NFAIL(15)+NFAIL(16)
      WRITE(0,20)"  convergence problem                   ",S
      WRITE(0,*)
      WRITE(0,20)"Remaining good points                   ",TOT
      IF(TOT.GT.0)THEN
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:       "
       WRITE(0,*)
       WRITE(0,30)" MSUSYEFF: ",MSUSYEFFN,MSUSYEFFNN
       WRITE(0,30)" MMESS: ",MMESSN,MMESSNN
       WRITE(0,30)" TANB: ",TBN,TBNN
       WRITE(0,30)" LAMBDA: ",LN,LNN
       WRITE(0,30)" ALAMBDA: ",ALN,ALNN
       IF(ALFLAG.EQ.1)
     .  WRITE(0,20)"(ALAMBDA is not an input parameter)"
       WRITE(0,30)"MUEFF: ",MUN,MUNN
       WRITE(0,20)"(MU is not an input parameter)"
       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,30)" XIF: ",XIFN,XIFNN
        WRITE(0,20)"(XIF is not an input parameter)"
       ELSE
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,20)"(KAPPA is not an input parameter)"
        WRITE(0,30)" XIF: ",XIFN,XIFNN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,20)"(MS^2 is not an input parameter)"
       ELSE
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,20)"(XIS is not an input parameter)"
       ENDIF
       WRITE(0,30)" MU': ",MUPN,MUPNN
       WRITE(0,30)" MS'^2: ",MSPN,MSPNN
       WRITE(0,30)" DH: ",DHN,DHNN
       
      ENDIF

 20   FORMAT(A40,I10)
 30   FORMAT(A15,2E15.4)

      END
