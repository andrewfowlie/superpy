      SUBROUTINE INITIALIZE()
      
*******************************************************************
*   This subroutine serves to
*     a) set default values for the SM parameters
*     b) set default values for the parameters (as limits on
*      sparticle masses) used for actual experimental constraints
*     c) read the tables that are used for actual experimental
*      constraints from the directory EXPCON
*******************************************************************

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER NBIN,I,J
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,ALEM0
      DOUBLE PRECISION GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU
      DOUBLE PRECISION MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2),OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd

      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd

*   The EXPCON_PATH variable is set:

      CALL getenv('EXPCON_PATH',EXPCON_PATH)
      if(EXPCON_PATH.eq.' ')  EXPCON_PATH='../EXPCON'
      
      call  getarg(0, EXPCON_PATH)
      I=len(EXPCON_PATH)+1
11    I=I-1      
      if(EXPCON_PATH(I:I).eq.' ') goto 11 
22     EXPCON_PATH(I:I)=' '
      I=I-1
      if(EXPCON_PATH(I:I).ne.'/') goto 22
      EXPCON_PATH(I:I)=' '
      EXPCON_PATH=catpath(EXPCON_PATH,'../EXPCON')      
*   SM inputs:

*   Alpha_EM(0)
      ALEM0= 1d0/137.036d0
*   Alpha_s at MZ:
      ALSMZ= 0.1172d0
*   Electroweak parameters:
      GF= 1.16639D-5
*   Alpha_em at MZ (used in RGES):
      ALEMMZ= 1d0/127.92d0
*   Z, W pole masses:
      MZ= 91.187d0
      MW= 80.42d0
*   Lepton masses:            
      MTAU= 1.777d0
      MMUON= 0.10566d0
*   Quark pole masses:
      MS= 0.19d0
      MC= 1.4d0
*      MBP= 4.94d0
      MBP=5.0d0
      MT= 175d0
*   Running MS_bar bottom mass MB at the scale MB:
      MB= 4.214d0      
*   Elements of the Kobayashi-Maskawa matrix:
      VUS= 0.22d0
      VCB= 0.04d0
      VUB= 0.004d0

*   Dark matter constraints
*   (Used only if OMGFLAG=1)
      OMGMIN=0.1187d0-5d0*0.0017d0
      OMGMAX=0.1187d0+5d0*0.0017d0
c      OMGMIN=1.d-10
c      OMGMAX=1.d0
      sigmaPiN=34d0
      sigmaS=42d0
      NBIN=10
      EMIN=1.D-3
      DO I=1,100
       x(I)=0d0
       dNdx(I)=0d0
      ENDDO
      DO I=1,NBIN
        x(I)=(NBIN-I)*DLOG10(Emin)/(NBIN-1d0)
      ENDDO

*   Collider constraints on sparticles:
*   Limit on the Z width from Z -> h(i) + a(j):
      GZMAX=5.78D-3
*   Limit on the inv. Z width from Z -> neutralinos:      
      GZINVMAX=1.71D-3
*   Limit on sigma(e+e- -> neutralinos (1,1)):      
      SIGNEU1=1.D-2
*   Limit on sigma(e+e- -> neutralinos (i,j)) (i*j > 1):            
      SIGNEU=1.D-1
*   Lower limit on chargino masses:
      MCHAMIN=103.5d0
*   Lower limit on slepton masses:
      MSLMIN=99.9d0
*   Lower limit on stau masses:      
      MSTMIN=93.2d0
*   Lower limit on squark masses:
      MSQMIN=100d0
*   Lower limit on gluino mass:
      MGLMIN=180d0
*   Lower limit on charged Higgs mass:
      MCMIN=78.6d0

      FILENAME=catpath(EXPCON_PATH,'hZind.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 111  READ(11,*,END=112,ERR=2)(hZind(I,J),J=1,2)
      I=I+1
      GOTO 111
 112  NhZind=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZbb.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 211  READ(11,*,END=212,ERR=2)(hZbb(I,J),J=1,2)
      I=I+1
      GOTO 211
 212  NhZbb=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZll.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 311  READ(11,*,END=312,ERR=2)(hZll(I,J),J=1,2)
      I=I+1
      GOTO 311
 312  NhZll=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZinv.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 411  READ(11,*,END=412,ERR=2)(hZinv(I,J),J=1,2)
      I=I+1
      GOTO 411
 412  NhZinv=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZjj.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 511  READ(11,*,END=512,ERR=2)(hZjj(I,J),J=1,2)
      I=I+1
      GOTO 511
 512  NhZjj=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZgg.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 611  READ(11,*,END=612,ERR=2)(hZgg(I,J),J=1,2)
      I=I+1
      GOTO 611
 612  NhZgg=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 711  READ(11,*,END=712,ERR=2)(hA4b(I,J),J=1,3)
      I=I+1
      GOTO 711
 712  NhA4b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA4tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 811  READ(11,*,END=812,ERR=2)(hA4tau(I,J),J=1,3)
      I=I+1
      GOTO 811
 812  NhA4tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA2b2tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 821  READ(11,*,END=822,ERR=2)(hA2b2tau(I,J),J=1,3)
      I=I+1
      GOTO 821
 822  NhA2b2tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA2tau2b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 831  READ(11,*,END=832,ERR=2)(hA2tau2b(I,J),J=1,3)
      I=I+1
      GOTO 831
 832  NhA2tau2b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAA6b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 911  READ(11,*,END=912,ERR=2)(AAA6b(I,J),J=1,3)
      I=I+1
      GOTO 911
 912  NAAA6b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAA6tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 921  READ(11,*,END=922,ERR=2)(AAA6tau(I,J),J=1,3)
      I=I+1
      GOTO 921
 922  NAAA6tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1101 READ(11,*,END=1102,ERR=2)(AAZ4b(I,J),J=1,3)
      I=I+1
      GOTO 1101
 1102 NAAZ4b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ4tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1111 READ(11,*,END=1112,ERR=2)(AAZ4tau(I,J),J=1,3)
      I=I+1
      GOTO 1111
 1112 NAAZ4tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ2b2tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1121 READ(11,*,END=1122,ERR=2)(AAZ2b2tau(I,J),J=1,3)
      I=I+1
      GOTO 1121
 1122 NAAZ2b2tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1201 READ(11,*,END=1202)cccc02(I,1)
      READ(11,*,END=1202,ERR=2)cccc02(I,2)
      I=I+1
      GOTO 1201
 1202 Ncccc02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1211 READ(11,*,END=1212)cccc04(I,1)
      READ(11,*,END=1212,ERR=2)cccc04(I,2)
      I=I+1
      GOTO 1211
 1212 Ncccc04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1221 READ(11,*,END=1222)cccc05(I,1)
      READ(11,*,END=1222,ERR=2)cccc05(I,2)
      I=I+1
      GOTO 1221
 1222 Ncccc05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1231 READ(11,*,END=1232)cccc06(I,1)
      READ(11,*,END=1232,ERR=2)cccc06(I,2)
      I=I+1
      GOTO 1231
 1232 Ncccc06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1241 READ(11,*,END=1242)cccc08(I,1)
      READ(11,*,END=1242,ERR=2)cccc08(I,2)
      I=I+1
      GOTO 1241
 1242 Ncccc08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1251 READ(11,*,END=1252)cccc1(I,1)
      READ(11,*,END=1252,ERR=2)cccc1(I,2)
      I=I+1
      GOTO 1251
 1252 Ncccc1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1301 READ(11,*,END=1302)ccgg02(I,1)
      READ(11,*,END=1302,ERR=2)ccgg02(I,2)
      I=I+1
      GOTO 1301
 1302 Nccgg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1311 READ(11,*,END=1312)ccgg04(I,1)
      READ(11,*,END=1312,ERR=2)ccgg04(I,2)
      I=I+1
      GOTO 1311
 1312 Nccgg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1321 READ(11,*,END=1322)ccgg05(I,1)
      READ(11,*,END=1322,ERR=2)ccgg05(I,2)
      I=I+1
      GOTO 1321
 1322 Nccgg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1331 READ(11,*,END=1332)ccgg06(I,1)
      READ(11,*,END=1332,ERR=2)ccgg06(I,2)
      I=I+1
      GOTO 1331
 1332 Nccgg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1341 READ(11,*,END=1342)ccgg08(I,1)
      READ(11,*,END=1342,ERR=2)ccgg08(I,2)
      I=I+1
      GOTO 1341
 1342 Nccgg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1351 READ(11,*,END=1352)ccgg1(I,1)
      READ(11,*,END=1352,ERR=2)ccgg1(I,2)
      I=I+1
      GOTO 1351
 1352 Nccgg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1401 READ(11,*,END=1402)cctt02(I,1)
      READ(11,*,END=1402,ERR=2)cctt02(I,2)
      I=I+1
      GOTO 1401
 1402 Ncctt02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1411 READ(11,*,END=1412)cctt04(I,1)
      READ(11,*,END=1412,ERR=2)cctt04(I,2)
      I=I+1
      GOTO 1411
 1412 Ncctt04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1421 READ(11,*,END=1422)cctt05(I,1)
      READ(11,*,END=1422,ERR=2)cctt05(I,2)
      I=I+1
      GOTO 1421
 1422 Ncctt05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1431 READ(11,*,END=1432)cctt06(I,1)
      READ(11,*,END=1432,ERR=2)cctt06(I,2)
      I=I+1
      GOTO 1431
 1432 Ncctt06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1441 READ(11,*,END=1442)cctt08(I,1)
      READ(11,*,END=1442,ERR=2)cctt08(I,2)
      I=I+1
      GOTO 1441
 1442 Ncctt08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1451 READ(11,*,END=1452)cctt1(I,1)
      READ(11,*,END=1452,ERR=2)cctt1(I,2)
      I=I+1
      GOTO 1451
 1452 Ncctt1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1501 READ(11,*,END=1502)gggg02(I,1)
      READ(11,*,END=1502,ERR=2)gggg02(I,2)
      I=I+1
      GOTO 1501
 1502 Ngggg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1511 READ(11,*,END=1512)gggg04(I,1)
      READ(11,*,END=1512,ERR=2)gggg04(I,2)
      I=I+1
      GOTO 1511
 1512 Ngggg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1521 READ(11,*,END=1522)gggg05(I,1)
      READ(11,*,END=1522,ERR=2)gggg05(I,2)
      I=I+1
      GOTO 1521
 1522 Ngggg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1531 READ(11,*,END=1532)gggg06(I,1)
      READ(11,*,END=1532,ERR=2)gggg06(I,2)
      I=I+1
      GOTO 1531
 1532 Ngggg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1541 READ(11,*,END=1542)gggg08(I,1)
      READ(11,*,END=1542,ERR=2)gggg08(I,2)
      I=I+1
      GOTO 1541
 1542 Ngggg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1551 READ(11,*,END=1552)gggg1(I,1)
      READ(11,*,END=1552,ERR=2)gggg1(I,2)
      I=I+1
      GOTO 1551
 1552 Ngggg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1601 READ(11,*,END=1602)ttgg02(I,1)
      READ(11,*,END=1602,ERR=2)ttgg02(I,2)
      I=I+1
      GOTO 1601
 1602 Nttgg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1611 READ(11,*,END=1612)ttgg04(I,1)
      READ(11,*,END=1612,ERR=2)ttgg04(I,2)
      I=I+1
      GOTO 1611
 1612 Nttgg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1621 READ(11,*,END=1622)ttgg05(I,1)
      READ(11,*,END=1622,ERR=2)ttgg05(I,2)
      I=I+1
      GOTO 1621
 1622 Nttgg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1631 READ(11,*,END=1632)ttgg06(I,1)
      READ(11,*,END=1632,ERR=2)ttgg06(I,2)
      I=I+1
      GOTO 1631
 1632 Nttgg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1641 READ(11,*,END=1642)ttgg08(I,1)
      READ(11,*,END=1642,ERR=2)ttgg08(I,2)
      I=I+1
      GOTO 1641
 1642 Nttgg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1651 READ(11,*,END=1652)ttgg1(I,1)
      READ(11,*,END=1652,ERR=2)ttgg1(I,2)
      I=I+1
      GOTO 1651
 1652 Nttgg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1701 READ(11,*,END=1702)tttt02(I,1)
      READ(11,*,END=1702,ERR=2)tttt02(I,2)
      I=I+1
      GOTO 1701
 1702 Ntttt02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1711 READ(11,*,END=1712)tttt04(I,1)
      READ(11,*,END=1712,ERR=2)tttt04(I,2)
      I=I+1
      GOTO 1711
 1712 Ntttt04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1721 READ(11,*,END=1722)tttt05(I,1)
      READ(11,*,END=1722,ERR=2)tttt05(I,2)
      I=I+1
      GOTO 1721
 1722 Ntttt05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1731 READ(11,*,END=1732)tttt06(I,1)
      READ(11,*,END=1732,ERR=2)tttt06(I,2)
      I=I+1
      GOTO 1731
 1732 Ntttt06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1741 READ(11,*,END=1742)tttt08(I,1)
      READ(11,*,END=1742,ERR=2)tttt08(I,2)
      I=I+1
      GOTO 1741
 1742 Ntttt08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1751 READ(11,*,END=1752)tttt1(I,1)
      READ(11,*,END=1752,ERR=2)tttt1(I,2)
      I=I+1
      GOTO 1751
 1752 Ntttt1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'stblsn.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1761 READ(11,*,END=1762,ERR=2)(stblsn(I,J),J=1,2)
      I=I+1
      GOTO 1761
 1762 Nstblsn=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'stnc.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1771 READ(11,*,END=1772,ERR=2)(stnc(I,J),J=1,2)
      I=I+1
      GOTO 1771
 1772 Nstnc=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'sbnb.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1781 READ(11,*,END=1782,ERR=2)(sbnb(I,J),J=1,2)
      I=I+1
      GOTO 1781
 1782 Nsbnb=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'glsq.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1791 READ(11,*,END=1792,ERR=2)(glsq(I,J),J=1,2)
      I=I+1
      GOTO 1791
 1792 Nglsq=I-1
      CLOSE(11)

      RETURN

*   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

      END


      integer function strlen(st)

*     Logical length of a string (omitting blanks)
*
*     | F | O | R | T | R | A | N |  |  |  |  |  |  |  |
*     <------------ full length returned by len() ----->
*     <---- logical length ------> <- trailing blank -->
*
      implicit none
      character      st*(*)
      strlen = len(st)
      do while (st(strlen:strlen) .eq. ' ')
       strlen = strlen - 1
      enddo
      return
      end


      character*(*) function catpath(st1,st2)

*     st1 is a string containing the path.
*     st2 is the file name or an other directory name
*     return value is the string "st1/st2"

      implicit none
      character*(*) st1
      character*(*) st2
      integer strlen
      catpath = st1(1:strlen(st1)) // '/' // st2
      return
      end
