
      REAL FUNCTION DTOF(X)
      DOUBLE PRECISION X
      DTOF=X
      END

C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at weak scale
C          MBQ                  = bottom mass at weak scale
C          MLQ                  = tau mass at weak scale
C	   FBMA                 = b-Yukawa at mA scale
C	   VUQ                  = Hu vev at MSUSY
C	   VDQ                  = Hd vev at MSUSY


      subroutine isajet2micromegas()
      IMPLICIT NONE
C      DOUBLE PRECISION M0,MHF,A0,TANB,SGNMU,MT
      REAL DTOF 
      EXTERNAL ALDATA

      character*10 name

      COMMON/SSPAR/
#ifndef lt779
     $ GORGE, 
#endif           
     $ AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      SAVE /SSPAR/

C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/

      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST
      SAVE /SUGPAS/

      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/ 
      integer sing
      DOUBLE PRECISION  findValW
      INTEGER I,J,k, assignVal
      CHARACTER*6 NAMES(24),PNAMES(32)
      REAL*8 DMSS(32),sum, SMC1,SMC2 
      DATA NAMES/'MG3','mu','MH3','tb_Q','Mq1','Md1','Mu1','Ml1','Mr1', 
     &'Mq2','Md2','Mu2','Ml2','Mr2','Mq3','Md3','Mu3','Ml3','Mr3',
     &'At','Ab','Al','MG1','MG2'/ 
      real tb
      DATA PNAMES/'MSG' ,'MSuL','MSuR','MSdL','MSdR','MSsL',
     &            'MSsR','MScL','MScR','MSb1','MSb2','MSt1',
     &            'MSt2','MSne','MSnm','MSnl','MSeL','MSeR',
     &            'MSmL','MSmR','MSl1','MSl2','MNE1','MNE2',
     &            'MNE3','MNE4','MC1' ,'MC2' ,'Mh'   ,'MHH', 
     &            'MH3'  ,'MHc'/
C  MSS( 1)=glss  MSS( 2)=upl   MSS( 3)=upr   MSS( 4)=dnl   MSS( 5)=dnr   MSS( 6)=stl
C  MSS( 7)=str	 MSS( 8)=chl   MSS( 9)=chr   MSS(10)=b1    MSS(11)=b2    MSS(12)=t1
C  MSS(13)=t2    MSS(14)=nuel  MSS(15)=numl  MSS(16)=nutl  MSS(17)=el-   MSS(18)=er-
C  MSS(19)=mul-  MSS(20)=mur-  MSS(21)=tau1  MSS(22)=tau2  MSS(23)=z1ss  MSS(24)=z2ss
C  MSS(25)=z3ss  MSS(26)=z4ss  MSS(27)=w1ss  MSS(28)=w2ss  MSS(29)=hl0   MSS(30)=hh0
C  MSS(31)=ha0   MSS(32)=h+

       
      DO 3  I=1,24
   3     call assignValW(NAMES(I),dble(xisain(I)))

      call assignValW('vev',dsqrt(2*dble(VUQ)**2+dble(VDQ)**2))
      call assignValW('tb_Q',DBLE(VUQ/VDQ))
      call assignValW('QSUSY', DBLE(HIGFRZ))
      call assignValW('g3',DBLE(sqrt(ASMSS*16*ATAN(1.))))
      call assignValW('gY',DBLE( SQRT(.6)*GSS(1)))
      call assignValW('g2',DBLE( GSS(2))) 
      call assignValW('mH1_2',DBLE(GSS(13)))
      call assignValW('mH2_2',DBLE(GSS(14)))
      call assignValW('Yl',DBLE(GSS(4)))
      call assignValW('Yb',DBLE(GSS(5)))
      call assignValW('Yt',DBLE(GSS(6)))
      tb=VUQ/VDQ
      call assignValW('mA_2', DBLE(B*XISAIN(2)*(1+tb*tb)/tb))

      DO 4 I=1,32
   4  DMSS(I)=MSS(I)


      DO 41 I=23,28
41    DMSS(I)=-DMSS(I)    

      DO 5 I=1,32
5     call assignValW(PNAMES(I),DMSS(I))

	do 1  i=1,4
	do 1  j=1,4

	write(name,FMT='(A2,I1,I1)') 'Zn',j,i
        sing=1
        if((i.eq.3).or.(i.eq.4)) sing=-sing 
        if((j.eq.2).or.(j.eq.4)) sing=-sing

1	call assignValW(name,DBLE( sing*ZMIXSS(5-i,j)))

      SMC2=SIGN(1.,TAN(GAMMAL)*TAN(GAMMAR))
 

      call assignValW('Zu11',dble(-sin(GAMMAL)))
      call assignValW('Zu12',dble(-cos(GAMMAL)))
      call assignValW('Zu21',dble(-cos(GAMMAL)*SMC2))
      call assignValW('Zu22',dble(sin(GAMMAL)*SMC2))

      call assignValW('Zv11',dble(-sin(GAMMAR)))
      call assignValW('Zv12',dble(-cos(GAMMAR)))
      call assignValW('Zv21',dble(-cos(GAMMAR)))
      call assignValW('Zv22',dble( sin(GAMMAR)))
     
      call assignValW('alpha',dble(-ALFAH)) 
      
      call assignValW('Zl11',dcos(dble(thetal)))
      call assignValW('Zl21',-dsin(dble(thetal)))
      call assignValW('Zl12',dsin(dble(thetal))) 
      call assignValW('Zl22',dcos(dble(thetal)))

      call assignValW('Zt11',dcos(dble(thetat)))
      call assignValW('Zt21',-dsin(dble(thetat)))
      call assignValW('Zt12',dsin(dble(thetat))) 
      call assignValW('Zt22',dcos(dble(thetat)))

      call assignValW('Zb11',dcos(dble(thetab)))
      call assignValW('Zb21',-dsin(dble(thetab)))
      call assignValW('Zb12',dsin(dble(thetab))) 
      call assignValW('Zb22',dcos(dble(thetab)))

c      call assignValW('MSlth',dble(thetal))   
c      call assignValW('MSbth',dble(thetab))   
c      call assignValW('MStth',dble(thetat))   
       
      end

 
      integer FUNCTION isajetSUGRA(tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,
     >sgn,gMHu,gMHd,gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3)
      IMPLICIT NONE
      real*8 tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,
     >sgn,gMHu,gMHd,gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3

      INTEGER I
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST     
      SAVE /SUGPAS/
      
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/

      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/


      REAL DTOF 
      REAL Mt,Mhlf,A0,M0
      DOUBLE PRECISION  findValW
      INTEGER IALLOW
      XNUSUG(1)   =gMG1
      XNUSUG(2)   =gMG2
      XNUSUG(3)   =gMG3
       Mhlf=findValW('Mhlf')      
      XNUSUG(4)   =gAl
      XNUSUG(5)   =gAb
      XNUSUG(6)   =gAt
       A0=findValW('A0')
      XNUSUG(8)   =gMHu
      XNUSUG(7)   =gMHd
c
      XNUSUG(10)  =gMl2
      XNUSUG(15)  =gMl3
      XNUSUG(9)   =gMr2
      XNUSUG(14)  =gMr3
      XNUSUG(13)  =gMq2
      XNUSUG(18)  =gMq3
      XNUSUG(12)  =gMu2
      XNUSUG(11)  =gMd2
      XNUSUG(17)  =gMu3
      XNUSUG(16)  =gMd3
      M0=findValW('M0')

      isajetSUGRA=0
      ial3un=0
      XNRIN(2)=1.E20
      Mt=findValW('Mtp')  
      call sugra(M0,Mhlf,A0,DTOF(tb),DTOF(sgn),Mt,1)
      if(MHPNEG.EQ.1) NOGOOD=3
      if(nogood.ne.0) then
         isajetSUGRA=-iabs(nogood) 
         return
       endif

        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ Mt,IALLOW,1)
       

       call isajet2micromegas()         
       call assignValW('Am',findValW('Al')) 
       call assignValW('tb',tb)
       return 
      end

      integer FUNCTION isajetBDR(tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,
     >gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,EWmu,EWmA)
c------------------------------------------------------------------------
c     sugra with non-univ Higgs masses and input of mu(Q) and mA(Q) 
c     instead of mHu and mHd; note that some mHu!=mHd!=m0 have to 
c     be given nevertheless as initial guesses for the RGE iteration 
c
c     SK, 2 Oct 2007, last modified on 3/10/2007 

      IMPLICIT NONE
      real*8 tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,
     >gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,EWmu,EWmA

      INTEGER I
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1

      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,
     $MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST     
      SAVE /SUGPAS/
      
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/

      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/


      REAL DTOF 
      REAL Mt,Mhlf,A0,M0,sgn
      DOUBLE PRECISION  findValW
      INTEGER IALLOW

      COMMON/SSPAR/ 
#ifndef lt779
     $ GORGE, 
#endif    
     $ AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/

      if(EWmu.gt.0) then
        sgn=1
      else 
        sgn=-1
      endif
      XNUSUG(1)   =gMG1
      XNUSUG(2)   =gMG2
      XNUSUG(3)   =gMG3
       Mhlf=findValW('Mhlf')      
      XNUSUG(4)   =gAl
      XNUSUG(5)   =gAb
      XNUSUG(6)   =gAt
       A0=findValW('A0')
    
      XNUSUG(10)  =gMl2
      XNUSUG(15)  =gMl3
      XNUSUG(9)   =gMr2
      XNUSUG(14)  =gMr3
      XNUSUG(13)  =gMq2
      XNUSUG(18)  =gMq3
      XNUSUG(12)  =gMu2
      XNUSUG(11)  =gMd2
      XNUSUG(17)  =gMu3
      XNUSUG(16)  =gMd3
      M0=findValW('M0')
      XNUSUG(8)   =M0
      XNUSUG(7)   =M0

      XNUSUG(19)  =EWmu
      XNUSUG(20)  =EWmA
      MU=XNUSUG(19)
      AMHA=XNUSUG(20)
      TWOM1=-MU
      INUHM=1

      isajetBDR=0
      ial3un=0
      XNRIN(2)=1.E20
      Mt=findValW('Mtp')  
      call sugra(M0,Mhlf,A0,DTOF(tb),sgn,Mt,1)

      if(MHPNEG.EQ.1) NOGOOD=3
      if(nogood.ne.0) then
         isajetBDR=-iabs(nogood) 
         return
       endif

        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ Mt,IALLOW,1)
       
       call isajet2micromegas()         
       call assignValW('Am',findValW('Al')) 
       call assignValW('tb',tb)
       return 
      end



      INTEGER  function isajetAMSB(AM0,M32,tanb,sgnmu)
      IMPLICIT NONE
      DOUBLE PRECISION AM0,M32,TANB,SGNMU,MT
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1

      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,
     $MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST

      SAVE /SUGPAS/
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/

      REAL DTOF 
      DOUBLE PRECISION  findValW
      INTEGER I
      
      ial3un=0
      DO 10 I=1,7
10     xamin(i)=1   

      call sugra(DTOF(AM0),DTOF(M32),0.,DTOF(TANB),DTOF(SGNMU),
     &DTOF(findValW('Mtp')),7)
      if(MHPNEG.EQ.1) NOGOOD=3
      isajetAMSB=-iabs(NOGOOD)
              
      if(nogood.ne.0)  return
      call isajet2micromegas() 
      call assignValW('tb',tanb)
      call assignValW('Am',findValW('Al'))
      return 
      end      

      INTEGER  function isajetGMSB(Lambda, Mmess,tb,SGNMU,N5,cGrav,
     * Rsl,dmH_d2,dmH_u2,d_Y,n5_1,n5_2,n5_3)
C  Rsl       factor multiplying gaugino masses at M_mes  
C  dmH_d2   Higgs mass**2 shifts at M_mes                
C  dmH_u2   Higgs mass**2 shifts at M_mes                
C  d_Y      mass**2 shifts proportional to Y at M_mes    
C  n5_1     n5 values for U(1),                         
C  n5_2                   SU(2),                      
C  n5_3                   SU(3)                     


      IMPLICIT NONE
      REAL*8 Lambda, tb, Mmess, cGrav,SGNMU, N5 
      REAL*8 Rsl,dmH_d2,dmH_u2,d_Y,n5_1,n5_2,n5_3

      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,
     $MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST
      SAVE /SUGPAS/

      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/

      COMMON/SSPAR/
#ifndef lt779
     $ GORGE, 
#endif            
     $ AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA, 
     $VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/


      REAL DTOF 
      DOUBLE PRECISION  findValW
      ial3un=0
       AMGVSS=Lambda*Mmess*cGrav/SQRT(3.)/2.4E18
       XGMIN(7)=cGrav
       XGMIN(8)=Rsl 
       XGMIN(9)=dmH_d2 
       XGMIN(10)=dmH_u2
       XGMIN(11)=d_Y
       XGMIN(12)=n5_1  
       XGMIN(13)=n5_2
       XGMIN(14)=n5_3

      call sugra(DTOF(Lambda),DTOF(Mmess),DTOF(N5),DTOF(TB),DTOF(SGNMU),
     &DTOF(findValW('Mtp')),2)
C      if(MHPNEG.EQ.1) NOGOOD=3 
      isajetGMSB=-iabs(NOGOOD)        
      if(nogood.ne.0)  return
      call isajet2micromegas()
      call assignValW('tb',tb)         
      call assignValW('Am',findValW('Al'))
      return 

      end

      INTEGER FUNCTION isajetEwsbTot()

      IMPLICIT NONE
      EXTERNAL ALDATA
    
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,
     $MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1

      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST
      SAVE /SUGPAS/

c      REAL*8 tb,MG1,MG2,MG3,Am,Al,At,Ab,MH3,mu,Ml1,Ml2,Ml3,
c     >Mr1,Mr2,Mr3,Mq1,Mq2,Mq3,Mu1,Mu2,Mu3,Md1,Md2,Md3
       REAL*8 Mu3, Mq3, tb
      CHARACTER*6 NAMES(24)
      DATA NAMES/'gMG3','mu','MH3','tb',
     &'gMq1','gMd1','gMu1','gMl1','gMr1',
     &'gMq2','gMd2','gMu2','gMl2','gMr2',
     &'gMq3','gMd3','gMu3','gMl3','gMr3',
     &'gAt','gAb','gAl','gMG1','gMG2'/
      REAL*8 findValW
      INTEGER I

      REAL QSUSY,ASMB,MBMB,ASMT,MTMT,SUALFS,PI,GG

      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
      COMMON/SSPAR/
#ifndef lt779
     $ GORGE, 
#endif           
     $ AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ
      LOGICAL GORGE
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/

      REAL Mtp 
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      REAL*8 SSMQCD
      REAL AMASS
      SAVE /SUGXIN/

      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/ 

      real*8 ma, MZ
      DO 3  I=1,24
      xisain(I)= findValW(NAMES(I))
3     continue
      Mtp=findValW('Mtp')
      tb=findValW('tb')
      Mu3=findValW('gMu3')
      Mq3=findValW('gMq3') 
      mu=findValW('mu')
      MZ=91.187 
C findValW('MZ')
      ma=findValW('MH3')
      PI=4.*ATAN(1.)
      QSUSY=SQRT(Mq3*Mu3)
      HIGFRZ=QSUSY
      
      ALQCD4=0.177
      AMBT=AMASS(5)
      AMTP=Mtp
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SNGL(SSMQCD(DBLE(MBMB),DBLE(QSUSY)))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)

      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SNGL(SSMQCD(DBLE(MTMT),DBLE(QSUSY)))
      MLQ=1.7463
C     For MSSM solution TANBQ=TANB; for SUGRA, TANBQ=/ TANB
      ALFAEM=1./128.
      SN2THW=.232
      AMW=80.423
      GG=SQRT(4*PI*ALFAEM/SN2THW)
      VUQ=SQRT(2*AMW**2/GG**2/(1.+1./tb**2))
      VDQ=VUQ/tb
      ASMSS=SUALFS(QSUSY**2,.31,AMTP,3)
      ASMSS=ASMSS*(1+ASMSS/4./PI)
      GSS(2)=6.52576625E-01
      GSS(1)=3.57532650E-01/sqrt(0.6)
c    GSS(4,5,6) =Yl,Yb,Yt      
      GSS(4)= MLQ/VDQ
      GSS(5)= MBQ/VDQ
      GSS(6)= MTQ/VUQ
c   mH2^2= ..
      GSS(14)=-(mu**2 +0.5*MZ**2) +(MZ**2+ma**2)/(tb**2+1)
c    mH1^2= ..
      GSS(13)=GSS(14) - (MZ**2+ma**2)*(1-tb**2)/(1+tb**2)    
      B=ma**2*tb/mu/(1+tb**2)
  
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ Mtp, isajetEwsbTot,0)

      IF(NOGOOD.NE.0) isajetEwsbTot=-NOGOOD
      if(isajetEwsbTot.LT.0) return
      call isajet2micromegas()
      call assignValW('Am',findValW('Al'))
      END



c***************************** main program ***************      

      integer model,err,paramExists
      real*8 findValW,Rsl,dmH_d2,dmH_u2,d_Y
      character*100  fInput
      character*100  fOutput 
      character*40 VISAJE
      fInput='slha.in'
      fOutput='slha.out'

      if(iargc().ge. 1) call getarg(1,fInput) 
      if(iargc().ge. 2) call getarg(2,fOutput)

      call readlesh(fInput)
      model=findValW('model')+0.1
      

      if(model.eq.0) then 
       err=isajetEwsbTot()
      else if(model.eq.1) then
        if((paramExists('mu').eq.1).and.(paramExists('MH3').eq.1)) then 
          err=isajetBDR(
     $    findValW('tb'),
     $    findValW('gMG1'),
     $    findValW('gMG2'),
     $    findValW('gMG3'),
     $    findValW('gAl'),
     $    findValW('gAt'),
     $    findValW('gAb'),
     $    findValW('gMl2'),
     $    findValW('gMl3'),
     $    findValW('gMr2'),
     $    findValW('gMr3'),
     $    findValW('gMq2'),
     $    findValW('gMq3'),
     $    findValW('gMu2'),
     $    findValW('gMu3'),
     $    findValW('gMd2'),
     $    findValW('gMd3'),
     $    findValW('mu'),
     $    findValW('MH3')
     $    )
          else 
          err=isajetSUGRA(
     $    findValW('tb'),
     $    findValW('gMG1'),
     $    findValW('gMG2'),
     $    findValW('gMG3'),
     $    findValW('gAl'),
     $    findValW('gAt'),
     $    findValW('gAb'),
     $    findValW('sgn'),
     $    findValW('gMHu'),
     $    findValW('gMHd'),
     $    findValW('gMl2'),
     $    findValW('gMl3'),
     $    findValW('gMr2'),
     $    findValW('gMr3'),
     $    findValW('gMq2'),
     $    findValW('gMq3'),
     $    findValW('gMu2'),
     $    findValW('gMu3'),
     $    findValW('gMd2'),
     $    findValW('gMd3')
     $     )
         endif
       else if(model.eq.2) then 
        Rsl=0
        dmH_d2=0
        dmH_u2=0
        d_Y =0
        err=isajetGMSB(findValW('Lambda'),findValW('Mmess'),
     $  findValW('tb'),findValW('sgn'),findValW('N5'),findValW('Cgrav'),
     $  Rsl,dmH_d2,dmH_u2,d_Y,
     $   findValW('N5_1'),findValW('N5_2'),findValW('N5_3'))
       else if(model.eq.3) then 
          err=isajetAMSB(findValW('M0'),findValW('M32'),findValW('tb'),
     $     findValW('sgn'))
       endif
       call writelesh(-err,VISAJE(),fOutput)
      end
