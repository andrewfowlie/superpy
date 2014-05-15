      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*   On input:      
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*      m0, M1/2, A0 at the scale MGUT      
*      non-universal parameters are allowed in the Higgs/gaugino sector
*
*      The input file contains lower and upper bounds
*      as well as number of steps for these parameters
*
*   On output:
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
*  ERRORS: IFAIL = 0..17
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
*      PROB(26) =/= 0  lightest neutralino is not LSP or < 511 keV
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(30) =/= 0  excluded by WMAP (checked only if OMGFLAG>=1)
*      PROB(31) =/= 0  excluded by Xenon100 (checked only if OMGFLAG=2 or 4)
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
      PARAMETER (NFL=17,NPROB=50,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDIM,IDUM
      INTEGER ITER,ITERMU,MUFLAG,SFFLAG,Q2FIX,IM
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,OMGFLAG,MAFLAG,GMUFLAG,HFLAG

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,GUTEST,DELMB
      DOUBLE PRECISION DETM,R,K,MUM,MDM,MSM,MUT,MDT,MST,Q2,Q2MIN
      DOUBLE PRECISION M0,M12,A0,SIGMU,MUFAIL,MUSTEP,MUINIT
      DOUBLE PRECISION M0MIN,M0MAX,M12MIN,M12MAX,TBMIN,TBMAX
      DOUBLE PRECISION A0MIN,A0MAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION MHDMIN,MHDMAX,MHUMIN,MHUMAX,MSMIN,MSMAX
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,AKMIN,AKMAX
      DOUBLE PRECISION MUMIN,MUMAX,XIFMIN,XIFMAX,XISMIN,XISMAX
      DOUBLE PRECISION MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX
      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN
      DOUBLE PRECISION MHUN,MHUNN,MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN
      DOUBLE PRECISION AKN,AKNN,MUN,MUNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION MU3GUT,MD3GUT,RAN2
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT

      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,
     . MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DELMB/DELMB
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/MINMAX/M0MIN,M0MAX,M12MIN,M12MAX,TBMIN,TBMAX,
     . A0MIN,A0MAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . MHDMIN,MHDMAX,MHUMIN,MHUMAX,MSMIN,MSMAX,
     . LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,AKMIN,AKMAX,
     . MUMIN,MUMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX
      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN,AKN,AKNN,MUN,MUNN,
     . XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/STEPS/NTOT,IDUM
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/MSAVE/MUM,MDM,MSM,MUT,MDT,MST,IM
      COMMON/IDIM/IDIM
      COMMON/GMUFLAG/GMUFLAG,HFLAG

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      ITOT=0
      IDIM=0
      IFAIL=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      M0N=1d99
      M0NN=-1d99
      M12N=1d99
      M12NN=-1d99
      A0N=1d99
      A0NN=-1d99
      TBN=1d99
      TBNN=-1d99
      M1N=1d99
      M1NN=-1d99
      M2N=1d99
      M2NN=-1d99
      M3N=1d99
      M3NN=-1d99
      LN=1d99
      LNN=-1d99
      KN=1d99
      KNN=-1d99
      ALN=1d99
      ALNN=-1d99
      AKN=1d99
      AKNN=-1d99
      MUN=1d99
      MUNN=-1d99
      MHDN=1d99
      MHDNN=-1d99
      MHUN=1d99
      MHUNN=-1d99
      MSN=1d99
      MSNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      M3HN=1d99
      M3HNN=-1d99

*   Beginning of the scan

      DO ITOT=1,NTOT

      IF(NTOT.EQ.1)THEN

       M0=M0MIN
       M12=M12MIN
       PAR(3)=TBMIN
       A0=A0MIN
       PAR(1)=LMIN
       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSE
        ALINP=ALMIN
       ENDIF
       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSE
        AKINP=AKMIN
       ENDIF
       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSE
        M1INP=M1MIN
       ENDIF
       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSE
        M2INP=M2MIN
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSE
        M3INP=M3MIN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFMIN
       ELSE
        PAR(2)=KMIN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISMIN
       ELSE
        MSINP=MSMIN
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        PAR(4)=MUMIN
        PAR(2)=KMIN
        XIFINP=XIFMIN
        XISINP=XISMIN
       ELSE
        IF(MHDFLAG.EQ.0)THEN
         MHDINP=M0**2
        ELSE
         MHDINP=MHDMIN
        ENDIF
        IF(MHUFLAG.EQ.0)THEN
         MHUINP=M0**2
        ELSE
         MHUINP=MHUMIN
        ENDIF
       ENDIF
       MUPINP=MUPMIN
       MSPINP=MSPMIN
       M3HINP=M3HMIN

      ELSE

       IF(M0MIN.EQ.M0MAX)THEN
        M0=M0MIN
       ELSE
        M0=M0MIN+(M0MAX-M0MIN)*RAN2(IDUM)
       ENDIF
       IF(M12MIN.EQ.M12MAX)THEN
        M12=M12MIN
       ELSE
        M12=M12MIN+(M12MAX-M12MIN)*RAN2(IDUM)
       ENDIF
       IF(TBMIN.EQ.TBMAX)THEN
        PAR(3)=TBMIN
       ELSE
        PAR(3)=TBMIN+(TBMAX-TBMIN)*RAN2(IDUM)
       ENDIF
       IF(A0MIN.EQ.A0MAX)THEN
        A0=A0MIN
       ELSE
        A0=A0MIN+(A0MAX-A0MIN)*RAN2(IDUM)
       ENDIF
       IF(LMIN.EQ.LMAX)THEN
        PAR(1)=LMIN
       ELSE
        PAR(1)=LMIN*(LMAX/LMIN)**RAN2(IDUM)
       ENDIF
       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSEIF(ALMIN.EQ.ALMAX)THEN
        ALINP=ALMIN
       ELSE
        ALINP=ALMIN+(ALMAX-ALMIN)*RAN2(IDUM)
       ENDIF
       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSEIF(AKMIN.EQ.AKMAX)THEN
        AKINP=AKMIN
       ELSE
        AKINP=AKMIN+(AKMAX-AKMIN)*RAN2(IDUM)
       ENDIF
       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSEIF(M1MIN.EQ.M1MAX)THEN
        M1INP=M1MIN
       ELSE
        M1INP=M1MIN+(M1MAX-M1MIN)*RAN2(IDUM)
       ENDIF
       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSEIF(M2MIN.EQ.M2MAX)THEN
        M2INP=M2MIN
       ELSE
        M2INP=M2MIN+(M2MAX-M2MIN)*RAN2(IDUM)
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSEIF(M3MIN.EQ.M3MAX)THEN
        M3INP=M3MIN
       ELSE
        M3INP=M3MIN+(M3MAX-M3MIN)*RAN2(IDUM)
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        IF(XIFMIN.EQ.XIFMAX)THEN
         XIFINP=XIFMIN
        ELSE
         XIFINP=XIFMIN+(XIFMAX-XIFMIN)*RAN2(IDUM)
        ENDIF
       ELSE
        IF(KMIN.EQ.KMAX)THEN
         PAR(2)=KMIN
        ELSE
         PAR(2)=KMIN+(KMAX-KMIN)*RAN2(IDUM)
        ENDIF
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        IF(XISMIN.EQ.XISMAX)THEN
         XISINP=XISMIN
        ELSE
         XISINP=XISMIN+(XISMAX-XISMIN)*RAN2(IDUM)
        ENDIF
       ELSE
        IF(MSFLAG.EQ.0)THEN
         MSINP=M0**2
        ELSEIF(MSMIN.EQ.MSMAX)THEN
         MSINP=MSMIN
        ELSE
         MSINP=MSMIN+(MSMAX-MSMIN)*RAN2(IDUM)
        ENDIF
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        IF(MUMIN.EQ.MUMAX)THEN
         PAR(4)=MUMIN
        ELSE
         PAR(4)=MUMIN+(MUMAX-MUMIN)*RAN2(IDUM)
        ENDIF
        IF(KMIN.EQ.KMAX)THEN
         PAR(2)=KMIN
        ELSE
         PAR(2)=KMIN+(KMAX-KMIN)*RAN2(IDUM)
        ENDIF
        IF(XIFMIN.EQ.XIFMAX)THEN
         XIFINP=XIFMIN
        ELSE
         XIFINP=XIFMIN+(XIFMAX-XIFMIN)*RAN2(IDUM)
        ENDIF
        IF(XISMIN.EQ.XISMAX)THEN
         XISINP=XISMIN
        ELSE
         XISINP=XISMIN+(XISMAX-XISMIN)*RAN2(IDUM)
        ENDIF
       ELSE
        IF(MHDFLAG.EQ.0)THEN
         MHDINP=M0**2
        ELSEIF(MHDMIN.EQ.MHDMAX)THEN
         MHDINP=MHDMIN
        ELSE
         MHDINP=MHDMIN+(MHDMAX-MHDMIN)*RAN2(IDUM)
        ENDIF
        IF(MHUFLAG.EQ.0)THEN
         MHUINP=M0**2
        ELSEIF(MHUMIN.EQ.MHUMAX)THEN
         MHUINP=MHUMIN
        ELSE
         MHUINP=MHUMIN+(MHUMAX-MHUMIN)*RAN2(IDUM)
        ENDIF
       ENDIF
       IF(MUPMIN.EQ.MUPMAX)THEN
        MUPINP=MUPMIN
       ELSE
        MUPINP=MUPMIN+(MUPMAX-MUPMIN)*RAN2(IDUM)
       ENDIF
       IF(MSPMIN.EQ.MSPMAX)THEN
        MSPINP=MSPMIN
       ELSE
        MSPINP=MSPMIN+(MSPMAX-MSPMIN)*RAN2(IDUM)
       ENDIF
       IF(M3HMIN.EQ.M3HMAX)THEN
        M3HINP=M3HMIN
       ELSE
        M3HINP=M3HMIN+(M3HMAX-M3HMIN)*RAN2(IDUM)
       ENDIF
      ENDIF

!      WRITE(0,*)""
!      WRITE(0,*)"------------------------------------------------------"
!      WRITE(0,*)""
!      WRITE(0,*)"Point ",ITOT
!      WRITE(0,*)""
!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)"M0 =",M0
!      IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)WRITE(0,*)"M12 =",M12
!      WRITE(0,*)"TANB =",PAR(3)
!      IF(MAFLAG.NE.-5)WRITE(0,*)"SIGMU =",SIGMU
!      WRITE(0,*)"A0 =",A0
!      IF(M1FLAG.EQ.1)WRITE(0,*)"M1 =",M1INP
!      IF(M2FLAG.EQ.1)WRITE(0,*)"M2 =",M2INP
!      IF(M3FLAG.EQ.1)WRITE(0,*)"M3 =",M3INP
!      IF(MHUFLAG.EQ.1 .AND. MAFLAG.NE.-5)WRITE(0,*)"MHU =",MHUINP
!      IF(MHDFLAG.EQ.1 .AND. MAFLAG.NE.-5)WRITE(0,*)"MHD =",MHDINP
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)WRITE(0,*)"MS =",MSINP
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4 .OR. MAFLAG.EQ.-5)
!     . WRITE(0,*)"KAPPA =",PAR(2)
!      IF(ALFLAG.EQ.1)WRITE(0,*)"AL =",ALINP
!      IF(AKFLAG.EQ.1)WRITE(0,*)"AK =",AKINP
!      IF(MAFLAG.EQ.-5)WRITE(0,*)"MUEFF =",PAR(4)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)WRITE(0,*)"XIF =",XIFINP
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)WRITE(0,*)"XIS =",XISINP
!      WRITE(0,*)"MUP =",MUPINP
!      WRITE(0,*)"MSP =",MSPINP
!      WRITE(0,*)"M3H =",M3HINP
!      WRITE(0,*)""

*   Initialization of PROB

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO      

*   Initialization of algorithm parameters

      SFFLAG=0
      MUINIT=SIGMU*DSQRT(PAR(3)*Q2MIN)
      MUSTEP=SIGMU*10d0
 1    MUINIT=MUINIT+MUSTEP
      MUFAIL=MUINIT-2d0*MUSTEP
      MUFLAG=0
      MUT=0d0
      MDT=0d0
      MST=0d0
      IM=0

*   Guess parameters at Q2/MGUT

      IF(Q2FIX.EQ.0)THEN
       Q2=MAX(M0**2+4d0*M12**2,Q2MIN)
      ENDIF
      R=(1d0+PAR(3)**2)/(1.29d0*PAR(3)**2)
      K=(1d0-R)*(A0-2.24d0*M12)**2+7.84d0*M12**2
      PAR(5)=ALINP-R/2d0*(A0-2.24d0*M12)-.59d0*M12
      PAR(6)=AKINP
      PAR(7)=(1d0-R/2d0)*M0**2+7.02d0*M12**2-R/6d0*K
      PAR(8)=(1d0-R)*M0**2+6.6d0*M12**2-R/3d0*K
      PAR(9)=M0**2+6.55d0*M12**2
      PAR(10)=M0**2+.52d0*M12**2
      PAR(11)=M0**2+.15d0*M12**2
      PAR(12)=A0-R*(A0-2.24d0*M12)-3.97d0*M12
      PAR(13)=A0-R/6d0*(A0-2.24d0*M12)-3.93d0*M12
      PAR(14)=A0-.69d0*M12
      PAR(15)=M0**2+7.02d0*M12**2
      PAR(16)=M0**2+6.6d0*M12**2
      PAR(17)=M0**2+6.55d0*M12**2
      PAR(18)=M0**2+.52d0*M12**2
      PAR(19)=M0**2+.15d0*M12**2
      PAR(20)=.4d0*M1INP
      PAR(21)=.8d0*M2INP
      PAR(22)=2.4d0*M3INP
      PAR(23)=DSQRT(Q2)
      PAR(24)=DSQRT(Q2)
      PAR(25)=A0-.69d0*M12
      DELMB=.1d0
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       PAR(2)=PAR(1)/5d0
      ELSE
       XIFGUT=0d0
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       MSM=M0**2
      ELSE
       XISGUT=0d0
      ENDIF
      IF(MAFLAG.EQ.-5)THEN
       MUM=10d0*M0**2
       MDM=M0**2
       MSM=M0**2
      ELSE
       PAR(4)=MUINIT
      ENDIF
      IFAIL=0

!      WRITE(0,*)"New try"
!      WRITE(0,*)""
!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"AL =",PAR(5)
!      WRITE(0,*)"AK =",PAR(6)
!      WRITE(0,*)"ATOP =",PAR(12)
!      WRITE(0,*)"ABOT =",PAR(13)
!      WRITE(0,*)"ATAU =",PAR(14)
!      WRITE(0,*)"AMUON =",PAR(25)
!      WRITE(0,*)"MQ3 =",PAR(7)
!      WRITE(0,*)"MU3 =",PAR(8)
!      WRITE(0,*)"MD3 =",PAR(9)
!      WRITE(0,*)"MQ =",PAR(15)
!      WRITE(0,*)"MU =",PAR(16)
!      WRITE(0,*)"MD =",PAR(17)
!      WRITE(0,*)"ML3 =",PAR(10)
!      WRITE(0,*)"ME3 =",PAR(11)
!      WRITE(0,*)"ML =",PAR(18)
!      WRITE(0,*)"ME =",PAR(19)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
!       WRITE(0,*)"K =",PAR(2)
!      ELSE
!      WRITE(0,*)"XIF =",XIFGUT
!       ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
!       WRITE(0,*)"MS =",MSM
!      ELSE
!       WRITE(0,*)"XIS =",XISGUT
!      ENDIF
!      IF(MAFLAG.EQ.-5)THEN
!       WRITE(0,*)"MHU =",MUM
!       WRITE(0,*)"MHD =",MDM
!       WRITE(0,*)"MS =",MSM
!      ELSE
!       WRITE(0,*)"MU =",PAR(4)
!      ENDIF
!      WRITE(0,*)""

*   Check for singular parameters l, tan(beta), mu

      IF(PAR(1)*PAR(3)*PAR(4).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Guess for GUT scale and GUT couplings

      CALL RGES(PAR,PROB,IFAIL)
      IF(IFAIL.NE.0.OR.PROB(27).NE.0d0)THEN
!       WRITE(0,*)"RGE integration problem 1"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=11)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=11
        GOTO 11
       ENDIF
      ENDIF

*   External loop to compute the soft parameters at Q2

      ITER=0
 21   ITER=ITER+1
!      WRITE(0,*)"ITER =",ITER
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RGESINV(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 2"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=13)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=13
        GOTO 11
       ENDIF
      ENDIF

*   Internal loop to compute mu, k and ms

      ITERMU=0
 22   ITERMU=ITERMU+1
!      WRITE(0,*)"ITERMU =",ITERMU
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Negative sfermion mass squared"
       IF(MAFLAG.NE.-5)THEN
        IF(DABS(MUINIT).GT.30d0)THEN
         SFFLAG=1
         MUINIT=MUINIT-2d0*MUSTEP
         GOTO 1
        ELSE
!         WRITE(0,*)"Exit (IFAIL=8)"
!         WRITE(0,*)""
!         WRITE(0,*)""
         IFAIL=8
         GOTO 11
        ENDIF
       ENDIF
      ENDIF

      CALL MINIMIZE(PAR,CHECK)
      IF(ITER.GT.10.AND.DETM.GE.0d0)MUFLAG=1

      IF(CHECK.GT.1.D-12.AND.ITERMU.LT.10) GOTO 22
      IF(CHECK.GT.1.D-9.AND.ITERMU.LT.100) GOTO 22
      IF(CHECK.GT.1.D-6)THEN
!       WRITE(0,*)"Exit (IFAIL=14)"
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=14
       GOTO 11
      ENDIF
      
      CALL RGES(PAR,PROB,IFAIL)
      IF(IFAIL.NE.0.OR.PROB(27).NE.0d0)THEN
!       WRITE(0,*)"RGE integration problem 3"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=11)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=11
        GOTO 11
       ENDIF
      ENDIF

      CALL RGESUNI(PAR,IFAIL,GUTEST)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 4"
       IF(DABS(MUINIT).LT.100d0 .AND. DABS(MUSTEP).GT.1.D-1
     .    .AND. MAFLAG.NE.-5)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=12)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=12
        GOTO 11
       ENDIF
      ENDIF

      IF(GUTEST.GT.1.D-12.AND.ITER.LT.10) GOTO 21
      IF(GUTEST.GT.1.D-9.AND.ITER.LT.100) GOTO 21
      IF(GUTEST.GT.1.D-6.AND.ITER.LT.1000) GOTO 21
      IF(GUTEST.GT.1.D-3)THEN
!       WRITE(0,*)"No convergence"
       IF(DABS(MUSTEP).GT.1.D-1 .AND. MAFLAG.NE.-5)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSE
!        WRITE(0,*)"Exit (IFAIL=15)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=15
        GOTO 11
       ENDIF
      ENDIF

*   Check if correct EWSB

      IF(DETM.LE.0d0 .AND. MAFLAG.NE.-5)THEN
!       WRITE(0,*)"Convergence in a false minimum"
       IF(MUFLAG.EQ.1 .AND. DABS(MUSTEP).GT.1.D-1)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSEIF(MUFLAG.EQ.0 .AND. DABS(MUSTEP).LT.DABS(MUFAIL))THEN
        MUFAIL=MUFAIL-MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
        ITER=0
!        WRITE(0,*)"New try"
!        WRITE(0,*)"MUFAIL =",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSEIF(DETM.LT.0d0)THEN
!        WRITE(0,*)"Exit (IFAIL=17)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=17
        GOTO 11
       ENDIF
      ENDIF

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
      CALL SUSYLHCSP(PAR,PROB)

*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)

*   Anom. magn. moment of the Muon

      IF(GMUFLAG.EQ.1)CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*   Relic density

      CALL RELDEN(PAR,PROB)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)IFAIL=10
      ENDDO

*   Computation of the fine-tuning

      CALL FTPAR(PAR,1)

*   Sparticle decays

      CALL NMSDECAY(PAR)

*   Recording of the results

 11   CALL OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT)

      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       M0N=MIN(M0,M0N)
       M0NN=MAX(M0,M0NN)
       M12N=MIN(M12,M12N)
       M12NN=MAX(M12,M12NN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       A0N=MIN(A0,A0N)
       A0NN=MAX(A0,A0NN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       ALN=MIN(ALGUT,ALN)
       ALNN=MAX(ALGUT,ALNN)
       AKN=MIN(AKGUT,AKN)
       AKNN=MAX(AKGUT,AKNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       MHDN=MIN(MHDGUT,MHDN)
       MHDNN=MAX(MHDGUT,MHDNN)
       MHUN=MIN(MHUGUT,MHUN)
       MHUNN=MAX(MHUGUT,MHUNN)
       M1N=MIN(M1GUT,M1N)
       M1NN=MAX(M1GUT,M1NN)
       M2N=MIN(M2GUT,M2N)
       M2NN=MAX(M2GUT,M2NN)
       M3N=MIN(M3GUT,M3N)
       M3NN=MAX(M3GUT,M3NN)
       XIFN=MIN(XIFGUT,XIFN)
       XIFNN=MAX(XIFGUT,XIFNN)
       XISN=MIN(XISGUT,XISN)
       XISNN=MAX(XISGUT,XISNN)
       MSN=MIN(MSGUT,MSN)
       MSNN=MAX(MSGUT,MSNN)
       MUPN=MIN(MUPGUT,MUPN)
       MUPNN=MAX(MUPGUT,MUPNN)
       MSPN=MIN(MSPGUT,MSPN)
       MSPNN=MAX(MSPGUT,MSPNN)
       M3HN=MIN(M3HGUT,M3HN)
       M3HNN=MAX(M3HGUT,M3HNN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      ENDDO

*   Summary of the scanning:
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
      INTEGER N0,NLOOP,NBER,NPAR,ERR,NTOT,ISEED,GMUFLAG,HFLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,PFLAG,NMSFLAG

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
      DOUBLE PRECISION M0MIN,M0MAX,M12MIN,M12MAX,TBMIN,TBMAX
      DOUBLE PRECISION A0MIN,A0MAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION MHDMIN,MHDMAX,MHUMIN,MHUMAX,MSMIN,MSMAX
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,AKMIN,AKMAX
      DOUBLE PRECISION MUMIN,MUMAX,XIFMIN,XIFMAX,XISMIN,XISMAX
      DOUBLE PRECISION MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MINMAX/M0MIN,M0MAX,M12MIN,M12MAX,TBMIN,TBMAX,
     . A0MIN,A0MAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . MHDMIN,MHDMAX,MHUMIN,MHUMAX,MSMIN,MSMAX,
     . LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,AKMIN,AKMAX,
     . MUMIN,MUMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX
      COMMON/STEPS/NTOT,ISEED
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO      

*   INITIALIZATION OF THE SCANNING PARAMETERS
      SIGMU=1d99
      M0MIN=1d99
      M0MAX=1d99
      M12MIN=1d99
      M12MAX=1d99
      TBMIN=1d99
      TBMAX=1d99
      A0MIN=1d99
      A0MAX=1d99
      M1MIN=1d99
      M1MAX=1d99
      M2MIN=1d99
      M2MAX=1d99
      M3MIN=1d99
      M3MAX=1d99
      MHDMIN=1d99
      MHDMAX=1d99
      MHUMIN=1d99
      MHUMAX=1d99
      MSMIN=1d99
      MSMAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=1d99
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      AKMIN=1d99
      AKMAX=1d99
      MUMIN=1d99
      MUMAX=1d99
      XIFMIN=1d99
      XIFMAX=1d99
      XISMIN=1d99
      XISMAX=1d99
      MUPMIN=0d0
      MUPMAX=1d99
      MSPMIN=0d0
      MSPMAX=1d99
      M3HMIN=0d0
      M3HMAX=1d99
      NTOT=0

*   DEFAULT VALUES FOR FLAGS
      OMGFLAG=0
      PFLAG=0
      NMSFLAG=0
      M1FLAG=0
      M2FLAG=0
      M3FLAG=0
      MHDFLAG=0
      MHUFLAG=0
      MSFLAG=1
      ALFLAG=0
      AKFLAG=0
      GMUFLAG=1
      HFLAG=0

*   DEFAULT VALUE FOR THE RANDOM SEED
      ISEED=-1

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
       ENDDO
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
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL
 
*   READ GUT PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.17) M0MIN=VAL
       IF(IX.EQ.18) M0MAX=VAL
       IF(IX.EQ.27) M12MIN=VAL
       IF(IX.EQ.28) M12MAX=VAL
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL
       IF(IX.EQ.4) SIGMU=VAL
       IF(IX.EQ.57) A0MIN=VAL
       IF(IX.EQ.58) A0MAX=VAL
 
*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.107) M1MIN=VAL
       IF(IX.EQ.108) M1MAX=VAL
       IF(IX.EQ.207) M2MIN=VAL
       IF(IX.EQ.208) M2MAX=VAL
       IF(IX.EQ.307) M3MIN=VAL
       IF(IX.EQ.308) M3MAX=VAL
       IF(IX.EQ.217) MHDMIN=VAL
       IF(IX.EQ.218) MHDMAX=VAL
       IF(IX.EQ.227) MHUMIN=VAL
       IF(IX.EQ.228) MHUMAX=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.647) AKMIN=VAL
       IF(IX.EQ.648) AKMAX=VAL
       IF(IX.EQ.657) MUMIN=VAL
       IF(IX.EQ.658) MUMAX=VAL
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
       IF(IX.EQ.727) M3HMIN=VAL
       IF(IX.EQ.728) M3HMAX=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.0) NTOT=IVAL
       IF(IX.EQ.1) ISEED=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(M0MIN.EQ.1d99)THEN
       WRITE(0,1)"M0_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(A0MIN.EQ.1d99)THEN
       WRITE(0,1)"A0_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"LAMBDA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANB_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0 .AND. MUMIN.EQ.1d99)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(SIGMU.NE.1d99 .AND. MUMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUEFF AND SIGMU CANNOT BE GIVEN"
       ERR=1
      ENDIF
      IF(M12MIN.EQ.1d99 .AND. M12MAX.NE.1d99)THEN
       WRITE(0,1)"M12_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(M1MIN.EQ.1d99 .AND. M1MAX.NE.1d99)THEN
       WRITE(0,1)"M1_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2MIN.EQ.1d99 .AND. M2MAX.NE.1d99)THEN
       WRITE(0,1)"M2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3MIN.EQ.1d99 .AND. M3MAX.NE.1d99)THEN
       WRITE(0,1)"M3_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MHDMIN.EQ.1d99 .AND. MHDMAX.NE.1d99)THEN
       WRITE(0,1)"MHD^2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MHUMIN.EQ.1d99 .AND. MHUMAX.NE.1d99)THEN
       WRITE(0,1)"MHU^2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSMIN.EQ.1d99 .AND. MSMAX.NE.1d99)THEN
       WRITE(0,1)"MS^2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KMIN.EQ.1d99 .AND. KMAX.NE.1d99)THEN
       WRITE(0,1)"KAPPA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"ALAMBDA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKMIN.EQ.1d99 .AND. AKMAX.NE.1d99)THEN
       WRITE(0,1)"AKAPPA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUMIN.EQ.1d99 .AND. MUMAX.NE.1d99)THEN
       WRITE(0,1)"MUEFF_min MUST BE GIVEN IN BLOCK EXTPAR"
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
      IF(MUPMIN.EQ.1d99 .AND. MUPMAX.NE.1d99)THEN
       WRITE(0,1)"MUP_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSPMIN.EQ.1d99 .AND. MSPMAX.NE.1d99)THEN
       WRITE(0,1)"MSP_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3HMIN.EQ.1d99 .AND. M3HMAX.NE.1d99)THEN
       WRITE(0,1)"M3H_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUMIN.NE.1d99)THEN
       IF(MHDMIN.NE.1d99)THEN
        WRITE(0,1)"BOTH MUEFF AND MHD^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(MHUMIN.NE.1d99)THEN
        WRITE(0,1)"BOTH MUEFF AND MHU^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(MSMIN.NE.1d99)THEN
        WRITE(0,1)"BOTH MUEFF AND MS^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
      ELSE
       IF(XIFMIN.NE.1d99 .AND. KMIN.NE.1d99)THEN
        WRITE(0,1)"BOTH XIF AND KAPPA CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(XISMIN.NE.1d99 .AND. MSMIN.NE.1d99)THEN
        WRITE(0,1)"BOTH XIS AND MS^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
      ENDIF

*   Set default values

      IF(MUMIN.NE.1d99)THEN
       IF(KMIN.EQ.1d99)KMIN=0d0
       IF(XIFMIN.EQ.1d99)XIFMIN=0d0
       IF(XISMIN.EQ.1d99)XISMIN=0d0
      ELSE
       IF(KMIN.EQ.1d99 .AND. XIFMIN.EQ.1d99)XIFMIN=0d0
       IF(MSMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)XISMIN=0d0
      ENDIF

      IF(M0MAX.EQ.1d99)M0MAX=M0MIN
      IF(M12MAX.EQ.1d99)M12MAX=M12MIN
      IF(A0MAX.EQ.1d99)A0MAX=A0MIN
      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(M1MAX.EQ.1d99)M1MAX=M1MIN
      IF(M2MAX.EQ.1d99)M2MAX=M2MIN
      IF(M3MAX.EQ.1d99)M3MAX=M3MIN
      IF(MHUMAX.EQ.1d99)MHUMAX=MHUMIN
      IF(MHDMAX.EQ.1d99)MHDMAX=MHDMIN
      IF(MSMAX.EQ.1d99)MSMAX=MHDMIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(AKMAX.EQ.1d99)AKMAX=AKMIN
      IF(MUMAX.EQ.1d99)MUMAX=MUMIN
      IF(XIFMAX.EQ.1d99)XIFMAX=XIFMIN
      IF(XISMAX.EQ.1d99)XISMAX=XISMIN
      IF(MUPMAX.EQ.1d99)MUPMAX=MUPMIN
      IF(MSPMAX.EQ.1d99)MSPMAX=MSPMIN
      IF(M3HMAX.EQ.1d99)M3HMAX=M3HMIN

*   Set MAFLAG, SCANFLAGS

      IF(MUMIN.NE.1d99)THEN
       MAFLAG=-5
       MHDFLAG=1
       MHUFLAG=1
      ELSE
       MAFLAG=-1
       IF(KMIN.NE.1d99)MAFLAG=MAFLAG-2
       IF(MSMIN.NE.1d99)MAFLAG=MAFLAG-1
       IF(MHDMIN.NE.1d99)MHDFLAG=1
       IF(MHUMIN.NE.1d99)MHUFLAG=1
       IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        IF(MSMIN.EQ.M0MIN**2 .AND. MSMAX.EQ.M0MAX**2)MSFLAG=0
       ENDIF
      ENDIF

      IF(M1MIN.NE.1d99)M1FLAG=1
      IF(M2MIN.NE.1d99)M2FLAG=1
      IF(M3MIN.NE.1d99)M3FLAG=1
      IF(ALMIN.NE.1d99)ALFLAG=1
      IF(AKMIN.NE.1d99)AKFLAG=1

      IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)THEN
       IF(M12MIN.EQ.1d99)THEN
        WRITE(0,1)"M12_min MUST BE GIVEN IN BLOCK MINPAR"
        ERR=1
       ENDIF
      ELSE
       IF(M12MIN.NE.1d99)WRITE(0,1)"WARNING: M12 IS NOT USED"
      ENDIF

      IF(KMIN.EQ.0d0 .AND. KMAX.EQ.0d0)THEN
       IF((AKMIN.NE.0d0 .AND. AKMIN.NE.1d99) .OR.
     .  (AKMIN.EQ.0d0 .AND. AKMAX.NE.0d0))THEN
        WRITE(0,1)
     .  "WARNING KAPPA_min=KAPPA_max=0 => AKAPPA_min=AKAPPA_max=0"
       ENDIF
       AKMIN=0d0
       AKMAX=0d0
       AKFLAG=1
      ENDIF

*   Total number of points

      IF(NTOT.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

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

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT) 

*********************************************************************      
*   Subroutine writing all the results in the the output file.
*********************************************************************      
 
      IMPLICIT NONE 
 
      INTEGER NBIN,IFAIL,ITOT,NTOT,IDIM,I,J
      INTEGER NRES,NSUSY,NGUT,NMES,IMAX,JMAX
      PARAMETER (NSUSY=14,NGUT=21,NMES=15,IMAX=1,JMAX=100)
 
      DOUBLE PRECISION RES(IMAX,JMAX),PAR(*),PROB(*),SIG(3,8),S
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
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT
      DOUBLE PRECISION MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT,M0,M12,A0,SIGMU
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION MHUQ,MHDQ,MSQ,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION DELMB,PX,PA(6),PB(2),PL(7),PK(8),MH(3),MMH(3)
      DOUBLE PRECISION DMH(3),MA(2),MMA(2),DMA(2),MHC,MMHC,DMHC
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
*
      DOUBLE PRECISION chartot2(2),chartot(2),chartotmulti(2)
      DOUBLE PRECISION brcharst1(2),brcharst2(2),brcharsb1(2),
     .          brcharsb2(2),brcharsupl(2),brcharsupr(2),
     .          brcharsdownl(2),
     .          brcharsdownr(2),brcharsnel(2),brcharsn1(2),
     .          brcharsn2(2),brcharsell(2),brcharstau1(2),brcharselr(2),
     .          brcharstau2(2),brcharhcneut(2,5),brcharwneut(2,5),
     .          brcharzchic,brcharHchic(3),brcharAchic(2)
      DOUBLE PRECISION brntaunut(2,5),brnelnue(2,5),brnmunumu(2,5),
     .          brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .          brglupdb(2),brglchsb(2),brgltopbb(2),
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .          brneutsb2(5),
     .          brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .          brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .          brneutsn2(5),brneutsell(5),brneutselr(5),
     .          brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .          brneuthcchar(5,2),brneutzneut(5,5),
     .          brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .          brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .          brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .          brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .          brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .          brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .          brglch(5),brglst(5),brgltop(5),brglbot(5)
*
      DOUBLE PRECISION selltot,selltot2,selltot3,
     . selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     . sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .          brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .          brsntauneut(5),brsntauchar(2),brsntau1hcstau(2),
     .          brsntau1wstau(2),brstau1neut(5),brstau2neut(5),
     .          brstau1char(2),
     .          brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .          brstau1wsn(2),brstau2wsn(2),brstau2H(3),brstau2A(2),
     .          brstau2ztau
      DOUBLE PRECISION brselrstau,brselrstaustar,brstau2stau1star,
     .    brstau2stau1,brstau2stau1nn,
     .    brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .          brsuprnup(5),brsuprcdow(2),brsuprglui,
     .          brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .          brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION 
     .          stoptot(2),stoptot2(2),stoptotmulti(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .          brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .          brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .          brstsntau(2,2),brstsel(2,2),brstsnel(2),
     .          brstbsbst(2,2),brstbbsbt(2,2),brsttausbnu(2,2),
     .          brstelsbnu(2,2),brstupsbdow(2,2),
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottotmulti(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .          brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .          brsb1glui,brsb2glui,brsb1wst(2),
     .          brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .          brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .          brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
      DOUBLE PRECISION gluitot,gluitot2,gluitotmulti,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),
     .         brgobt(5),brgoud(2),brgocs(2),brgotb(2),brhcst1b,brwst1b
*************************************************************
*
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartotmulti
      COMMON/CHARGINO_BR_2BD/brcharst1,brcharst2,brcharsb1,
     .          brcharsb2,brcharsupl,brcharsupr,brcharsdownl,
     .          brcharsdownr,brcharsnel,brcharsn1,
     .          brcharsn2,brcharsell,brcharstau1,brcharselr,
     .          brcharstau2,brcharhcneut,brcharwneut,
     .          brcharzchic,brcharHchic,brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .          brnupdb,brnchsb,brntopbb,
     .          brglupdb,brglchsb,brgltopbb,
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs, 
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
*
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,
     .selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     .sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .          brselrcharnue,brsnellneut,brsnellchar,
     .          brsntauneut,brsntauchar,brsntau1hcstau,
     .          brsntau1wstau,brstau1neut,brstau2neut,brstau1char,
     .          brstau2char,brstau1hcsn,brstau2hcsn,
     .          brstau1wsn,brstau2wsn,brstau2H,brstau2A,brstau2ztau
      COMMON/SLEPTON_BR_3BD/brselrstau,brselrstaustar,brstau2stau1star,
     .   brstau2stau1,brstau2stau1nn,
     .   brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .   brsellstau1star,brsellstau1,brsellstau1nutau,
     .   brsnestau1star,brsnestau1,brsnestau1nutau
*
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .          brsuprnup,brsuprcdow,brsuprglui,
     .          brsdowlndow,brsdowlchup,brsdowlglui,
     .          brsdowrndow,brsdowrchup,brsdowrglui
*
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptotmulti,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .          brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .          brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .          brstsntau,brstsel,brstsnel,
     .          brstbsbst,brstbbsbt,brsttausbnu,
     .          brstelsbnu,brstupsbdow,
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottotmulti
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .          brsb2chart,brsb1hcst,brsb2hcst,
     .          brsb1glui,brsb2glui,brsb1wst,
     .          brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .          brsbtstsb,brsbtbstb,brsbtaustnu,
     .          brsbelstnu,brsbupstdow,brsbsnel,
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
*
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitotmulti,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*
      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
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
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFHIGM/MH,MMH,DMH,MA,MMA,DMA,MHC,MMHC,DMHC
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/IDIM/IDIM
      COMMON/RES/RES

      NRES=30

      IF(IFAIL.EQ.0)THEN
      
       IDIM=IDIM+1      
       RES(IDIM,1)=M0
       RES(IDIM,2)=M12
       RES(IDIM,3)=A0
       RES(IDIM,4)=ALINP
       RES(IDIM,5)=PAR(4)
       RES(IDIM,6)=PAR(3)
       RES(IDIM,7)=PAR(1)
       RES(IDIM,8)=XIFINP
       RES(IDIM,9)=XISINP

       DO I=1,3
        RES(IDIM,7+3*I)=SMASS(I)
        RES(IDIM,8+3*I)=SCOMP(I,3)**2
        RES(IDIM,9+3*I)=BRNEU(I,1,1)
       ENDDO
       DO I=1,2
        RES(IDIM,17+2*I)=PMASS(I)
        RES(IDIM,18+2*I)=PCOMP(I,2)**2
       ENDDO
       RES(IDIM,23)=CMASS
       RES(IDIM,24)=chi2gam+chi2bb+chi2zz
       RES(IDIM,25)=DABS(MNEU(1))
       RES(IDIM,26)=NEU(1,1)**2
       RES(IDIM,27)=NEU(1,3)**2+NEU(1,4)**2
       RES(IDIM,28)=NEU(1,5)**2
       RES(IDIM,29)=OMG
       RES(IDIM,30)=PROB(31)

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

      INTEGER I,S,TOT,NTOT,NFAIL(*)
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,OMGFLAG,MAFLAG

      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN
      DOUBLE PRECISION MHUN,MHUNN,MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN
      DOUBLE PRECISION AKN,AKNN,MUN,MUNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN

      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN,AKN,AKNN,MUN,MUNN,
     . XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG

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

       WRITE(0,30)" LAMBDA: ",LN,LNN
       WRITE(0,30)" TANB: ",TBN,TBNN
       WRITE(0,30)" M0: ",M0N,M0NN
       IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)WRITE(0,30)" M12: ",M12N,M12NN
       IF(M1FLAG.EQ.1)WRITE(0,30)"M1: ",M1N,M1NN
       IF(M2FLAG.EQ.1)WRITE(0,30)"M2: ",M2N,M2NN
       IF(M3FLAG.EQ.1)WRITE(0,30)"M3: ",M3N,M3NN
       WRITE(0,30)" A0: ",A0N,A0NN
       IF(ALFLAG.EQ.1)WRITE(0,30)"ALAMBDA: ",ALN,ALNN
       IF(AKFLAG.EQ.1)WRITE(0,30)"AKAPPA: ",AKN,AKNN

       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,30)" XIF: ",XIFN,XIFNN
        WRITE(0,20)"(XIF is not an input parameter)"
       ELSEIF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,20)"(KAPPA is not an input parameter)"
        WRITE(0,30)" XIF: ",XIFN,XIFNN
       ENDIF

       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,20)"(MS^2 is not an input parameter)"
       ELSEIF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        IF(MSFLAG.EQ.0)THEN
         WRITE(0,20)"(MS^2=M0^2)                     "
        ELSE
         WRITE(0,30)" MS^2: ",MSN,MSNN
        ENDIF
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,20)"(XIS is not an input parameter)"
       ENDIF

       IF(MAFLAG.EQ.-5)THEN
        WRITE(0,30)"MUEFF: ",MUN,MUNN
        WRITE(0,30)"KAPPA: ",KN,KNN
        WRITE(0,30)" XIF: ",XIFN,XIFNN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,30)"MHD^2: ",MHDN,MHDNN
        WRITE(0,20)"(MHD^2 is not an input parameter)"
        WRITE(0,30)"MHU^2: ",MHUN,MHUNN
        WRITE(0,20)"(MHU^2 is not an input parameter)"
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,20)"(MS^2 is not an input parameter)"
       ELSE
        WRITE(0,30)"MUEFF: ",MUN,MUNN
        WRITE(0,20)"(MU is not an input parameter)"
        IF(MHDFLAG.EQ.1)WRITE(0,30)"MHD^2: ",MHDN,MHDNN
        IF(MHUFLAG.EQ.1)WRITE(0,30)"MHU^2: ",MHUN,MHUNN
       ENDIF

       WRITE(0,30)" MU': ",MUPN,MUPNN
       WRITE(0,30)" MS'^2: ",MSPN,MSPNN
       WRITE(0,30)" M3H^2: ",M3HN,M3HNN

      ENDIF

 20   FORMAT(A40,I10)
 30   FORMAT(A15,2E15.4)

      END
