c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: VEGAS.F
c     External auxiliary numerical routine: Vegas Monte Carlo
c     integration.  

c     CAUTION: standard rand() library function seems to cause problems
c     on some compilers, like ifort (Intel Linux compiler). Thus, below
c     we defined other RNG generator, sufficient for simple Monte Carlo
c     integration in gluon CDM (only place where it is used in the
c     current version of SUSY_FLAVOR). Depending on your preferences,
c     uncomment standard library RNG functions rand and srand in
c     routines veg_rand, veg_srand, or use the ones defined there,
c     result does not seem to change.

      double precision function veg_rand()
c     define here your favorite RNG for VEGAS integration
      implicit double precision(a-h,o-z)
c     rand() works fine for GNU Fortran
ccc   veg_rand = rand()
c     in case of problems, use generator defined below, it is 
      real temp(1)
      len = 1
      call RANMAR(temp, len)
      veg_rand = dble(temp(1))
      return
      end

      subroutine veg_srand(iseed)
c     define here seed base for your favorite RNG for VEGAS integration
      implicit double precision(a-h,o-z)
c     srand() works fine for GNU Fortran
ccc   call srand(iseed)
c     in case of problems, use seed for the self-defined generator
      IJ = 100*abs(iseed)/31328
      IJ = 100*abs(iseed) - 31328*IJ 
      KL = 200*abs(iseed)/30081
      KL = 200*abs(iseed) - 30081*KL
      call rmarin(ij,kl)
      return
      end

      subroutine RMARIN(IJ,KL)
C     This is the initialization routine for the random number generator RANMAR()
C     The seed variables can have values between:    0 <= IJ <= 31328
C                                                    0 <= KL <= 30081
      real U(97), C, CD, CM
      integer I97, J97
      logical TEST
      common /raset1/ U, C, CD, CM, I97, J97, TEST
      TEST = .FALSE.
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.
     *    KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value 
     *between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and         
     *30081'
            stop
      endif
      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169) 
      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         U(ii) = s
2     continue
      C =  362436.0  / 16777216
      CD = 7654321.0 / 16777216
c      CM = 16777213.0 /16777216
      CM = 1.0 - 3.0 / 16777216
      I97 = 97
      J97 = 33
      TEST = .TRUE.
      return
      end

      subroutine ranmar(RVEC, LEN)
C This is the random number generator proposed by George Marsaglia in
C Florida State University Report: FSU-SCRI-87-50 It was slightly
C modified by F. James to produce an array of pseudorandom numbers.
      REAL RVEC(*)
      real U(97), C, CD, CM
      integer I97, J97
      logical TEST
      common /raset1/ U, C, CD, CM, I97, J97, TEST
 
      integer ivec
 
      if( .NOT. TEST ) then
         print '(A)',' Call the init routine (RMARIN) before calling RAN        
     *MAR'  
         stop
      endif

      do 100 ivec = 1, LEN
         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         uni = uni - C
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RVEC(ivec) = uni
100   continue
      return
      end

C=======================================================================
C **********************************************************************
C
C     VEGAS VERSION INCLUDING HBOOK CALLS POSSIBILITY,
C     SIMPLIFIED BUT MORE CLEAR AND EQUALLY EFFICIENT.
C
C                              M.MARTINEZ  DESY   09/09/85
C -------------------------------------------------------------------
      SUBROUTINE VEGAS(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      integer seed
      COMMON/BVEG2/ND,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     . ,D(50,10),DI(50,10)
C           COMMON/VEGN0/NZERO
      DIMENSION X(10),XIN(50),R(50),IA(10)
      COMMON/VEGAS_RESULT/AVGI,SD,TI,TSI,nnew
      COMMON/VEGAS_ININIC/ITMX0,IT0,CHI2A,NDIM0
      COMMON/VEGAS_CONVER/ALPH
      common/vegas_random/seed
      DATA ALPH0/1.5d0/
      ND=50
C===========================================================
C A))  INITIALIZING SOME VARIABLES
C===========================================================
      ITMX0=ITMX
      NDIM0=NDIM
      alph = 0
      IF(ALPH.EQ.0.d0)ALPH=ALPH0
      CALLS=NCALL
      XND=ND
      NDM=ND-1
      call veg_srand(seed)
CC    IF(IGRAPH.EQ.0)CALL INBOOK(0,FUN,WEIGHT)
C.............................................
C  INITIALIZING CUMMULATIVE VARIABLES
C.............................................
      IT=0
      SI=0.d0
      SI2=0.d0
      SWGT=0.d0
      SCHI=0.d0
      SCALLS=0.d0
C.............................................
C  DEFINING THE INITIAL INTERVALS DISTRIBUTION
C.............................................
      RC=1.d0/XND
      DO 7 J=1,NDIM
      XI(ND,J)=1.d0
      DR=0.d0
      DO 7 I=1,NDM
      DR=DR+RC
      XI(I,J)=DR
7     CONTINUE
***      IF(NPRN.NE.0.d0AND.NPRN.NE.10)PRINT 200,NDIM,CALLS,IT,ITMX,
***     . ACC,ND,ALPH
***      IF(NPRN.EQ.10)PRINT 290,NDIM,CALLS,ITMX,ACC,ND
C===========================================================
C B))  ITERATIONS LOOP
C===========================================================
9     IT=IT+1
C.............................................
C  INITIALIZING ITERATION VARIABLES
C.............................................
      NZERO=0
      TI=0.d0
      SFUN2=0.d0
      DO 10 J=1,NDIM
      DO 10 I=1,ND
      D(I,J)=0.d0
      DI(I,J)=0.d0
10    CONTINUE
CC    IF(IGRAPH.EQ.0)CALL REBOOK(0,FUN,WEIGHT)
      DO 11 JJ=1,NCALL
      WGT=1.d0
C.............................................
C  COMPUTING THE POINT POSITION
C.............................................
      DO 15 J=1,NDIM
      XN=veg_rand()*XND+1.d0
      IA(J)=int(XN)
      XIM1=0.d0
      IF(IA(J).GT.1)XIM1=XI(IA(J)-1,J)
      XO=XI(IA(J),J)-XIM1
      X(J)=XIM1+(XN-IA(J))*XO
      WGT=WGT*XO*XND
15    CONTINUE
C.............................................
C  COMPUTING THE FUNCTION VALUE
C.............................................
      FUN=FXN(X)*WGT/CALLS
      IF(FUN.NE.0.d0)NZERO=NZERO+1
      FUN2=FUN*FUN
c      WEIGHT=WGT/CALLS
CC    IF(IGRAPH.EQ.0)CALL XBOOK(0,FUN,WEIGHT)
      TI=TI+FUN
      SFUN2=SFUN2+FUN2
      DO 16 J=1,NDIM
      IAJ=IA(J)
      DI(IAJ,J)=DI(IAJ,J)+FUN
      D(IAJ,J)=D(IAJ,J)+FUN2
16    CONTINUE
11    CONTINUE
c.............................................
c  we modify the original program in order to 
c  pass safely the case of all zero values of
c  the function calls; in such a case the program
c  returns in the common /results/ an additional
c  variable nnew=1
c.............................................
      if (nzero.lt.4) then
        nnew=1
        return
      else
        nnew=0
      endif
C.............................................
C  COMPUTING THE INTEGRAL AND ERROR VALUES
C.............................................
      TI2=TI*TI
      TSI=DSQRT((SFUN2*CALLS-TI2)/(CALLS-1.d0))
      WGT=TI2/TSI**2
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      SCALLS=SCALLS+CALLS
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=0.d0
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      SD=1.d0/DSQRT(SD)
c      ERR=SD*100.d0/AVGI
      IT0=IT
C.............................................
C  PRINTING AND PLOTTING
C.............................................
***      IF(NPRN.NE.10)PRINT 201,IT,TI,TSI,NZERO,AVGI,SD,ERR,CHI2A
***      IF(NPRN.EQ.10)PRINT 203,IT,TI,TSI,AVGI,SD,CHI2A
c      IF(NPRN.GE.0) GO TO 21
c      DO 20 J=1,NDIM
***      PRINT 202,J
c      XIN(1)=XI(1,J)
c      DO 2020 L=2,ND
c2020  XIN(L)=XI(L,J)-XI(L-1,J)
c20    continue
***      PRINT 204,(XI(I,J),XIN(I),DI(I,J),D(I,J),I=1,ND)
c21    CONTINUE
      IF(DABS(SD/AVGI).GT.DABS(ACC).AND.IT.LT.ITMX)GO TO 98
***      PRINT 777,AVGI,SD
CC    IF(IGRAPH.EQ.0)CALL BOOKIT(2,FUN,WEIGHT)
c      call ranget(seed)
      seed = seed + int(10000*(2*veg_rand()-1))
      RETURN
98    CONTINUE
CC    IF(IGRAPH.EQ.0)CALL BOOKIT(0,FUN,WEIGHT)
C===========================================================
C C))  REDEFINING THE GRID
C===========================================================
C.............................................
C  SMOOTHING THE F**2 VALUED STORED FOR EACH INTERVAL
C.............................................
      DO 23 J=1,NDIM
      XO=D(1,J)
      XN=D(2,J)
      D(1,J)=(XO+XN)/2.d0
      X(J)=D(1,J)
      DO 22 I=2,NDM
      D(I,J)=XO+XN
      XO=XN
      XN=D(I+1,J)
      D(I,J)=(D(I,J)+XN)/3.d0
      X(J)=X(J)+D(I,J)
22    CONTINUE
      D(ND,J)=(XN+XO)/2.d0
      X(J)=X(J)+D(ND,J)
23    CONTINUE
C.............................................
C  COMPUTING THE 'IMPORTANCE FUNCTION' OF EACH INTERVAL
C.............................................
      DO 28 J=1,NDIM
      RC=0.d0
      DO 24 I=1,ND
      R(I)=0.d0
      IF(D(I,J).LE.0.d0) GO TO 224
      IF(D(I,J).GT.1.D-100) THEN
         XO=X(J)/D(I,J)
ccc      R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
         R(I)=((XO-1.d0)/XO/DLOG(XO))**ALPH
      ENDIF
224   RC=RC+R(I)
24    CONTINUE
C.............................................
C  REDEFINING THE SIZE OF EACH INTERVAL
C.............................................
      RC=RC/XND
      K=0
      XN=0.d0
      DR=0.d0
      I=0
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
      XI(I,J)=XIN(I)
27    CONTINUE
      XI(ND,J)=1.d0
28    CONTINUE
C
      GO TO 9
C===========================================================
C D))  FORMATS FOR THE PRINTOUTS
C===========================================================
c200   FORMAT('0INPUT INTEGRATION PARAMETERS:  NDIM=',I3
c     .,8H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X
c     .,6H  ACC=,G9.3/28X,6H   ND=,I4/28X
c     .,8H  ALPHA=,G9.3///)
c290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5
c     . ,6H  ACC=,G9.3,6H   ND=,I4)
c201   FORMAT(//' ITERATION NO',I3,
c     .  '.   INTEGRAL  =',G14.8,'+/-',G10.4/
c     .  20X,'NUM(F=/=0)=',I10//
c     .  ' ACCUMULATED RESULTS.   INTEGRAL =',G14.8,
c     .  '+/-',G10.4 / 12X,' % ERROR =',G10.4,
c     .  'CHI**2 PER ITN   =',G10.4)
c202   FORMAT(14H0DATA FOR AXIS,I2 /
c     . 7X,'X',7X,'  DELT X  ',2X,'  DELT I  ',2X,' CONVCE',
c     . 14X,'X',7X,'  DELT X  ',2X,'  DELT I  ',2X,' CONVCE'/)
c204   FORMAT(1X,4G12.4,5X,4G12.4)
c203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
c777   FORMAT(' '//' ',25('+'),' FINAL RESULT ',25('+')//
c     . '        INTEGRAL VALUE       =',G14.6,'  +/-',G14.6//
c     . ' ',64('+'))
      END
