      SUBROUTINE RELDEN(PAR,PROB)

**********************************************************************
*   Subroutine for the computation of the dark matter relic density
*   PROB(26) =/= 0  lightest neutralino is not LSP
*   PROB(30) =/= 0 excluded by WMAP (checked only if OMGFLAG=/=0)
*   PROB(30)  = -1 LSP is not the lightest neutralino in micrOMEGAs
*   PROB(30)  = -2 Problem in micrOMEGAs
*   PROB(31) =/= 0 excluded by LUX (checked if |OMGFLAG|=2 or 4)
*
**********************************************************************

      IMPLICIT NONE

#include "micromegas/sources/micromegas.fh"

      CHARACTER name*10,mess*20

      INTEGER NORD(5),HORD(3),NBIN,OMGFLAG,MAFLAG
      INTEGER err,i,j

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,PI
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),AMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION SST,SSB,SSL,COSB,SINB,TANB
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION tab(250),OMG,OMGMIN,OMGMAX,Xf
      DOUBLE PRECISION sigmaV,x(100),dNdx(100),EMIN,LAM,LUX
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION Nmass,SCcoeff
      DOUBLE PRECISION pA0(2),pA5(2),nA0(2),nA5(2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8),MH(3),MMH(3)
      DOUBLE PRECISION DMH(3),MA(2),MMA(2),DMA(2),MHC,MMHC,DMHC
      DOUBLE PRECISION DELMB,vcsnngg
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)

      COMMON/DELMB/DELMB
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,AMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/LAM/LAM
      COMMON/EFFHIGM/MH,MMH,DMH,MA,MMA,DMA,MHC,MMHC,DMHC
      COMMON/EFFCOUP/PX,PA,PB,PL,PK

      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH

      DATA NORD/1,2,4,3,5/
      DATA HORD/2,1,3/
      real*8 omg_
    
*   Test on LSP

      IF(DABS(MNEU(1)).LT.511d-6)THEN
       PROB(26)=-DDIM(1d0,DABS(MNEU(1))/511d-6)
      ELSE
       PROB(26)=DDIM(1d0,MIN(DABS(MGL),DABS(MCHA(1)),MUR,MUL,MDR,
     .          MDL,MLR,MLL,MNL,MST1,MSB1,MSL1,MSNT)/DABS(MNEU(1)))
      ENDIF

      IF (OMGFLAG.EQ.0) RETURN

*   Input parameters:

      PI=4d0*DATAN(1d0)
      TANB=TANBQ
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      CALL assignValW('alfSMZ',ALSMZ)
      CALL assignValW('MbMb',MB)
      CALL assignValW('Mtp',MT)

      CALL assignValW('At',PAR(12))
      CALL assignValW('Ab',PAR(13))
      CALL assignValW('Al',PAR(14))

      CALL assignValW('Lambda',PAR(1))
      CALL assignValW('Kappa',PAR(2))
      CALL assignValW('aLambda',PAR(5))
      CALL assignValW('aKappa',PAR(6))
      CALL assignValW('tB',TANBQ)
      CALL assignValW('mu',PAR(4))

      CALL assignValW('MG1' ,PAR(20))         
      CALL assignValW('MG2' ,PAR(21))         
      CALL assignValW('MG3' ,PAR(22))         
      CALL assignValW('Ml2' ,DSQRT(PAR(18)))  
      CALL assignValW('Ml3' ,DSQRT(PAR(10)))  
      CALL assignValW('Mr2' ,DSQRT(PAR(19)))  
      CALL assignValW('Mr3' ,DSQRT(PAR(11)))  
      CALL assignValW('Mq2' ,DSQRT(PAR(15)))  
      CALL assignValW('Mq3' ,DSQRT(PAR(7)))   
      CALL assignValW('Mu2' ,DSQRT(PAR(16)))  
      CALL assignValW('Mu3' ,DSQRT(PAR(8)))   
      CALL assignValW('Md2' ,DSQRT(PAR(17)))  
      CALL assignValW('Md3' ,DSQRT(PAR(9)))    

      CALL assignValW('Mha',AMASS(1))
      CALL assignValW('Mhb',AMASS(2))
      CALL assignValW('MHc',CMASS)
      DO i=1,3
       WRITE(name,fmt='(A2,I1)') 'Mh',i
       CALL assignValW(name,SMASS(i))
      ENDDO

      DO i=1,5
       WRITE(name,fmt='(A3,I1)') 'MNE',i
       CALL assignValW(name,MNEU(i))
       DO j=1,5
         WRITE(name,fmt='(A2,I1,I1)') 'Zn',i,j
         CALL assignValW(name,NEU(i,NORD(j)))
        ENDDO
      ENDDO

      CALL assignValW("Zh11",SCOMP(1,2))
      CALL assignValW("Zh12",SCOMP(1,1))
      CALL assignValW("Zh13",SCOMP(1,3))
      CALL assignValW("Zh21",SCOMP(2,2))
      CALL assignValW("Zh22",SCOMP(2,1))
      CALL assignValW("Zh23",SCOMP(2,3))
      CALL assignValW("Zh31",SCOMP(3,2))
      CALL assignValW("Zh32",SCOMP(3,1))
      CALL assignValW("Zh33",SCOMP(3,3))

      CALL assignValW('Za11',SINB*PCOMP(1,1))
      CALL assignValW('Za12',COSB*PCOMP(1,1))
      CALL assignValW('Za13',PCOMP(1,2))
      CALL assignValW('Za21',SINB*PCOMP(2,1))
      CALL assignValW('Za22',COSB*PCOMP(2,1))
      CALL assignValW('Za23',PCOMP(2,2)) 

      CALL assignValW('MSl1',MSL1)
      CALL assignValW('MSl2',MSL2)
      CALL assignValW('Zl11',CSL)
      CALL assignValW('Zl12',SSL)
      CALL assignValW('Zl21',-SSL)
      CALL assignValW('Zl22',CSL)

      CALL assignValW('MSb1',MSB1)
      CALL assignValW('MSb2',MSB2)
      CALL assignValW('Zb11',CSB)
      CALL assignValW('Zb12',SSB)
      CALL assignValW('Zb21',-SSB)
      CALL assignValW('Zb22',CSB)

      CALL assignValW('MSt1',MST1)
      CALL assignValW('MSt2',MST2)
      CALL assignValW('Zt11',CST)
      CALL assignValW('Zt12',SST)
      CALL assignValW('Zt21',-SST)
      CALL assignValW('Zt22',CST)


      CALL assignValW('Zu11',UU(1,1))
      CALL assignValW('Zu12',UU(1,2))
      CALL assignValW('Zu21',UU(2,1))
      CALL assignValW('Zu22',UU(2,2))

      CALL assignValW('Zv11',VV(1,1))
      CALL assignValW('Zv12',VV(1,2))
      CALL assignValW('Zv21',VV(2,1))
      CALL assignValW('Zv22',VV(2,2))

      CALL assignValW('MSeL',MLL)
      CALL assignValW('MSeR',MLR)
      CALL assignValW('MSmL',MLL)
      CALL assignValW('MSmR',MLR)
      CALL assignValW('MSne',MNL)
      CALL assignValW('MSnm',MNL)
      CALL assignValW('MSnl',MSNT)
      CALL assignValW('MSuL',MUL)
      CALL assignValW('MSuR',MUR)
      CALL assignValW('MSdL',MDL)
      CALL assignValW('MSdR',MDR)
      CALL assignValW('MScL',MUL)
      CALL assignValW('MScR',MUR)
      CALL assignValW('MSsL',MDL)
      CALL assignValW('MSsR',MDR)
      CALL assignValW('MSG',MGL)
      CALL assignValW('MC1',MCHA(1))
      CALL assignValW('MC2',MCHA(2))

*   Improved Higgs potential

       CALL assignValW("la1",PL(1))
       CALL assignValW("la2",PL(2))
       CALL assignValW("la3",PL(3))
       CALL assignValW("la4",PL(4))
       CALL assignValW("la5",PL(5))
       CALL assignValW("la6",PL(6))
       CALL assignValW("la7",PL(7))
       CALL assignValW("la1s",PK(1))
       CALL assignValW("la2s",PK(2))
       CALL assignValW("la3s",PK(3))
       CALL assignValW("la4s",PK(4)) 
       CALL assignValW("la5s",PK(5))
       CALL assignValW("la6s",PK(6))
       CALL assignValW("la7s",PK(7))
       CALL assignValW("la8s",PK(8))
       CALL assignValW("aa1",PA(1))
       CALL assignValW("aa2",PA(2))
       CALL assignValW("aa3",PA(3))
       CALL assignValW("aa4",PA(4))
       CALL assignValW("aa5",PA(5)) 
       CALL assignValW("aa6",PA(6))
       CALL assignValW("B1",PB(1))
       CALL assignValW("B2",PB(2))
       CALL assignValW("X",PX)
       CALL assignValW("dMb",DELMB)

*   Sorting sparticles

      err=sortOddParticles(mess)
!      WRITE(*,*) 'micrOMEGAs : err=', err, ' mess=', mess
      IF(mess.ne.'~o1') THEN
        OMG=-1d0
        PROB(30)=-1d0
        RETURN
      ENDIF
      IF(err.ne.0) THEN
        OMG=-2d0
        PROB(30)=-2d0
        RETURN
      ENDIF

      call ModelConstIni(2,WIDTH,err)

*   Computing relic density

      call setVVdecay(0) ! set to 1/0 for vitual decays on/off
      call forceUG(1)    ! set to 1/0 for unitary gauge on/off
      OMG=darkOmega(Xf,1,1.D-3)
      IF(OMG.GT.0d0)THEN
       IF(OMGFLAG.GT.0)THEN
        PROB(30)=DDIM(OMG/OMGMAX,1d0)-DDIM(1d0,OMG/OMGMIN)
       ELSE
        PROB(30)=DDIM(OMG/OMGMAX,1d0)
       ENDIF
      ELSE
       PROB(30)=-1d0
      ENDIF
      IF (IABS(OMGFLAG).EQ.1) RETURN
      IF (IABS(OMGFLAG).EQ.3) GOTO 1

*  Computing WIMP-Nucleon cross sections
*  Muq/Mdq=0.553d0, Msq/Mdq=18.9d0

      CALL calcScalarQuarkFF(0.553d0,18.9d0,sigmaPiN,sigmaS)
      err=nucleonAmplitudes(FeScLoop,pA0,pA5,nA0,nA5)
      Nmass=0.939d0
      SCcoeff=4/PI*3.8937966d8*(Nmass*Mcdm/(Nmass+Mcdm))**2
      csPsi=SCcoeff*pA0(1)**2
      csNsi=SCcoeff*nA0(1)**2
      csPsd=3*SCcoeff*pA5(1)**2
      csNsd=3*SCcoeff*nA5(1)**2
      IF( pA0(1)*nA0(1) .lt. 0) csNsi=-csNsi
      IF( pA5(1)*nA5(1) .lt. 0) csNsd=-csNsd
      PROB(31)=DDIM(OMG/0.1187d0*csPsi/LUX(DABS(MNEU(1))),1d0)

      IF (IABS(OMGFLAG).EQ.2) RETURN

*  Computing indirect detection rate
 1    sigmaV=calcSpectrum(0,tab,NULL,NULL,NULL,NULL,NULL,err)
      IF (err.NE.0) sigmaV=0d0
      IF (sigmaV.NE.0d0) THEN
       DO I=1,NBIN
        dNdx(I)=zInterp(DLOG(10d0)*x(I),tab)*DLOG(10d0)
       ENDDO
      ENDIF

      END


      DOUBLE PRECISION FUNCTION LUX(M)

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=83)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/5.5256e0,5.5715e0,5.6071e0,5.6211e0,5.6681e0,5.7208e0,
     .5.7570e0,5.8309e0,5.9085e0,5.9824e0,6.0702e0,6.1589e0,6.2679e0,
     .6.3670e0,6.4960e0,6.6618e0,6.8335e0,7.0073e0,7.1939e0,7.4151e0,
     .7.6036e0,7.8572e0,8.0718e0,8.3435e0,8.5906e0,8.8332e0,9.1168e0,
     .9.3925e0,9.6694e0,1.0039e1,1.0477d1,
     .1.0822d1,1.1145d1,1.1460d1,1.1874d1,1.2229d1,1.2990d1,1.3665d1,
     .1.4468d1,1.5576d1,1.6764d1,1.8194d1,1.9829d1,2.1741d1,2.4801d1,
     .2.8665d1,3.4190d1,4.1161d1,5.0014d1,6.0198d1,6.9845d1,8.4085d1,
     .9.7513d1,1.1230d2,1.2863d2,1.4424d2,1.6522d2,1.8950d2,2.1265d2,
     .2.4323d2,2.7915d2,3.1981d2,3.6847d2,4.2311d2,4.8925d2,5.6165d2,
     .6.4729d2,7.4703d2,8.4216d2,9.7314d2,1.0884d3,1.2535d3,1.4167d3,
     .1.6330d3,1.8471d3,2.0871d3,2.3954d3,2.7020d3,3.1014d3,3.5339d3,
     .4.0257d3,4.5878d3,5d3/
      DATA Y/1.1198e-3,7.2451e-4,4.5403e-4,2.7959e-4,1.4650e-4,
     .8.4578e-5,4.9318e-5,3.0243e-5,1.7659e-5,1.0986e-5,7.3269e-6,
     .4.5548e-6,3.0539e-6,2.0448e-6,1.3154e-6,8.0084e-7,5.1173e-7,
     .3.5763e-7,2.4993e-7,1.6929e-7,1.1643e-7,8.0997e-8,5.7147e-8,
     .4.3286e-8,3.5039e-8,2.5960e-8,1.9821e-8,1.6403e-8,1.2626e-8,
     .1.0393e-8,9.1692d-9,7.7254d-9,6.5089d-9,5.4840d-9,4.6583d-9,
     .3.9248d-9,3.2799d-9,2.7861d-9,2.3666d-9,1.9778d-9,1.6528d-9,
     .1.4242d-9,1.2122d-9,1.0469d-9,9.1983d-10,8.2610d-10,7.7384d-10,
     .7.8663d-10,8.1065d-10,8.9414d-10,9.9400d-10,1.0816d-9,1.2139d-9,
     .1.3903d-9,1.5745d-9,1.7926d-9,2.0301d-9,2.2948d-9,2.6100d-9,
     .2.9664d-9,3.3496d-9,3.8613d-9,4.3537d-9,5.0093d-9,5.6254d-9,
     .6.4722d-9,7.4165d-9,8.3463d-9,9.4541d-9,1.0621d-8,1.2101d-8,
     .1.3854d-8,1.5953d-8,1.8198d-8,2.0898d-8,2.4013d-8,2.7574d-8,
     .3.1690d-8,3.6445d-8,4.1500d-8,4.7266d-8,5.3807d-8,5.9d-8/

      IF(M.LT.X(1) .OR. M.GE.X(N))THEN
       LUX=1d99
       RETURN
      ENDIF
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        LUX=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      subroutine PrintRelDen(PROB, ch)

      IMPLICIT NONE
      DOUBLE PRECISION PROB(*)
#include "../sources/micromegas/sources/micromegas.fh"

! Common
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION tab(250),OMG,OMGMIN,OMGMAX,Xf
      DOUBLE PRECISION sigmaV,x(100),dNdx(100),EMIN,LAM
      INTEGER  NBIN,OMGFLAG,MAFLAG
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/FLAGS/OMGFLAG,MAFLAG
      DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)
      COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU

! External 
       real*8 vcsnngg,vcsnngz
       external vcsnngg,vcsnngz
!  Local   
      integer ch,I,pdg(4),err
      real*8 weightCh,alph_,omgfo,v,nngg,nngz,ccoeff
      character*50 txt 


      IF(OMGFLAG.EQ.0)THEN

       WRITE(ch,900) 4,"# Omega h^2 not computed (OMGFLAG=0)"
      
      ELSE

       IF(OMG.GT.0d0)THEN
        WRITE(ch,899) "#"
        WRITE(ch,899) "BLOCK ABUNDANCE"
        WRITE(ch,920) 1, Mcdm/Xf,'T_f[GeV]'
        WRITE(ch,920) 4, OMG,"omega h^2"
        WRITE(ch,920) 3, vSigma(Mcdm/Xf,1.D-3,1,alph_)*2.9979d-26,
     .  'vSigma'
        WRITE(ch,899)'# contibutions to vSigma in percents'
        I=1
1111    err=vSigmaTch(I,weightCh, pdg,txt)
        if(err.eq.1 .and. weightCh.gt.1.E-2) then
         write(ch,921)6,I,pdg(1),pdg(2),pdg(3),pdg(4),100*weightCh,txt
         I=I+1  
         goto 1111
        endif
        omgfo=PRINTCHANNELS(Xf,1.D-3,1.D-4,1,0)
        WRITE(ch,899)'# contibutions to 1/Omega in percents'
        I=1
1112    err=omegaCh(I,weightCh, pdg,txt)
        if(err.eq.1 .and. weightCh.gt.1.E-2) then
         write(ch,921)7,I,pdg(1),pdg(2),pdg(3),pdg(4),100*weightCh,txt
         I=I+1  
         goto 1112
        endif
        WRITE(ch,899) "#"
        WRITE(ch,899) "BLOCK LSP"
        WRITE(ch,920) 0, Mcdm,"LSP mass"
        WRITE(ch,920) 1, dabs(NEU(1,1)), "bino"
        WRITE(ch,920) 2, dabs(NEU(1,2)), "wino" 
        WRITE(ch,920) 3, dabs(NEU(1,3)), "higgsino2"
        WRITE(ch,920) 4, dabs(NEU(1,4)), "higgsino1"
        WRITE(ch,920) 5, dabs(NEU(1,5)), "singlino"
        IF(IABS(OMGFLAG).EQ.2 .OR. IABS(OMGFLAG).EQ.4)THEN
         WRITE(ch,899) "#"
         WRITE(ch,899) 'BLOCK NDMCROSSSECT'
         WRITE(ch,920) 1,csPsi*1d-36,"# csPsi [cm^2]"
         WRITE(ch,920) 2,csNsi*1d-36,"# csNsi [cm^2]"
         WRITE(ch,920) 3,csPsd*1d-36,"# csPsd [cm^2]"
         WRITE(ch,920) 4,csNsd*1d-36,"# csNsd [cm^2]"
         WRITE(ch,925)"# Values used for sigma_piN,sigma_S",
     .   " (strange content of the proton)"
         WRITE(ch,924) "#",sigmapiN,"sigma_piN"
         WRITE(ch,924) "#",sigmaS,"sigma_S"
        ENDIF
        IF(IABS(OMGFLAG).EQ.3 .OR. IABS(OMGFLAG).EQ.4)THEN
         WRITE(ch,899) "#"
         WRITE(ch,899) 'BLOCK ANNIHILATION'
         v=0.001d0
         nngg=vcsnngg(v)*2.9979d-26
         nngz=vcsnngz(v)*2.9979d-26 
         ccoeff=sigmaV/(sigmaV+nngg+nngz)
         WRITE(ch,920) 0,sigmaV+nngg+nngz,'sigmaV [cm^2/s]'
         write(ch,923)1,nngg/(sigmaV+nngg+nngz),2,22,22,'~o1,~o1->A,A'
         write(ch,923)2,nngz/(sigmaV+nngg+nngz),2,22,23,'~o1,~o1->A,Z'
         I=3
1113     err=vSigmach(I-2,weightCh, pdg,txt)
         if(err.eq.1 .and. weightCh.gt.1.E-3) then
          write(ch,923) I, weightCh*ccoeff, 2, pdg(3),pdg(4),txt
          I=I+1
          goto 1113
         endif
        ENDIF

       ELSE
        WRITE(ch,900) 4,"# Problem in micrOMEGAs"
       ENDIF

      ENDIF

 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 908  FORMAT(2E16.8)
 910  FORMAT(E16.8,3X,A)
 920  FORMAT(1x,I3,1x, 1PE12.4,3x,'# ',1x,A) 
 921  FORMAT(1x,I3,1x,I3,1x,I9,1x,I9,1x,I9,1x,I9,1x,1PE12.4,3x,'# ',A)
 922  FORMAT(1x,I3,1x,1PE12.4,1x,I9,1x,I9,3x,'# ',A)
 923  FORMAT(1x,I3,1x,1PE12.4,1x,I2,1x,I9,1x,I9,3x,'# ',A)
 924  FORMAT(A,1P,E16.8,1X,A)
 925  FORMAT(A,A)

      end
