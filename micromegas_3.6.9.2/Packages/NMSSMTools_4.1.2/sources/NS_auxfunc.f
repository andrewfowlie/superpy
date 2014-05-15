c These routines are generalisations of corresponding routines from
c     SDECAY: A Fortran code for the decays of the supersymmetric 
c         particles in the MSSM
c     by M. Muhlleitner (Karlsruhe, Inst. Technol.),
c        A. Djouadi (Orsay, LPT & CERN, Theory Division),
c        Y. Mambrini (Orsay, LPT),
c     Comput.Phys.Commun.168:46-70 (2005), hep-ph/0311167.
c    SDECAY should be cited whenever NMSDECAY is used.
c
c
c ------------------ Function resum ------------------------------- c
c --- x: tree level, y: rad. corr., returns y or y_modif. --------- c

      double precision function resum(x,y)
      double precision x,y
     
      if(x.ne.0d0) then
        if(y.le.-x) then
          resum=x*y/(x-y)
        else
          resum=y
        endif
      else
        resum=0d0
      endif
      return
      end

c -------------------------------------------------------------------- c
c ------------------ the function lambda ----------------------------- c
c -------------------------------------------------------------------- c

      double precision function NS_lamb(x,y)      

      implicit double precision (a-h,k-z)

      NS_lamb=dsqrt((1.D0-x**2-y**2)**2-4.D0*x**2*y**2) 

      return
      end

c ---------------------------------------------------------------------c
c -------------------------------------------------------------------- c
c ---------------------- the integration limits ---------------------- c
c -------------------------------------------------------------------- c

      double precision function NS_ax(xmu1)
      implicit real*8(a-h,k-z)
      NS_ax=2.d0*dsqrt(xmu1)
      end

c -------------------------------------------------------------------- c

      double precision function NS_BX(xmu1,xmu2,xmu3)
      implicit real*8(a-h,k-z)
      NS_BX=1.D0+xmu1-(dsqrt(xmu2)+dsqrt(xmu3))**2
      end

c -------------------------------------------------------------------- c

      double precision function NS_ay(xmu1,xmu2,xmu3,x1)
      implicit real*8(a-h,k-z)

      a = 1.D0-x1+xmu1
      b = (x1-2.D0)*(x1-1.D0-xmu2+xmu3-xmu1)

      delta = (4.D0*xmu1-x1**2)*
     .        (4.D0*xmu2*xmu3-(x1-1.D0+xmu2+xmu3-xmu1)**2)

      if (delta.lt.0.D0) then
         NS_ay=0.D0
      else
         r1=(b-dsqrt(delta))/(2.D0*a)
         r2=(b+dsqrt(delta))/(2.D0*a)
         NS_ay=r1
      endif

      end

c -------------------------------------------------------------------- c

      double precision function NS_BY(xmu1,xmu2,xmu3,x1)
      implicit real*8(a-h,k-z)

      a = 1.D0-x1+xmu1
      b = (x1-2.D0)*(x1-1.D0-xmu2+xmu3-xmu1)

      delta = (4.D0*xmu1-x1**2)*
     .        (4.D0*xmu2*xmu3-(x1-1.D0+xmu2+xmu3-xmu1)**2)
      
      if (delta.lt.0.d0) then
         NS_BY=0.D0
      else
         r1=(b-dsqrt(delta))/(2.D0*a)
         r2=(b+dsqrt(delta))/(2.D0*a)
         NS_BY=r2
      endif

      end

c -------------------------------------------------------------------- c
c ----------------------- the integration routine -------------------- c
c -------------------------------------------------------------------- c

      SUBROUTINE NS_INTEG2(F,NS_AX,NS_BX,NS_AY,NS_BY,xmu1,xmu2,xmu3,
     .                     NX,NY,SUM)
      IMPLICIT double precision (A-H,O-Z)
      DOUBLE PRECISION NS_AY,NS_BY,NS_AX,NS_BX
      DIMENSION RX(1000),WX(1000),RY(1000),WY(1000)
      EXTERNAL F,NS_AY,NS_BY,NS_AX,NS_BX
      INTEGER NX,NY

      AXX=NS_AX(xmu1)
      BXX=NS_BX(xmu1,xmu2,xmu3)

      CALL NS_GAUS(NX,AXX,BXX,RX,WX)
  
      SX=0.D0
      DO  1 K=1,NX
         X=RX(K)
         AYX=NS_AY(xmu1,xmu2,xmu3,X)
         BYX=NS_BY(xmu1,xmu2,xmu3,X)
         CALL NS_GAUS(NY,AYX,BYX,RY,WY)
         SY=0.D0
         DO 2 J=1,NY
            SY=SY+WY(J)*F(X,RY(J))
 2       CONTINUE
         SX=SX+WX(K)*SY
 1    CONTINUE
      SUM=SX
      END

c -------------------------------------------------------------------- c

      SUBROUTINE NS_GAUS(N,A,B,X,W)
C     GAUSS-LEGENDRE INTEGRATION FROM A TO B (WEIGHT = 1.)
C     CALCULATES GAUSSIAN POINTS X AND WEIGHTS W
C                      FOR N=4,6,8,12,16,24,32 ;
C     IF N IS DIFFERENT FROM THESE VALUES THE PROGRAM DECOMPOSES
C     THE INTERVAL [A,B] AND N INTO CORRESPONDING PEACES
C     N MUST BE EVEN AND >= 4,IF IT IS NOT,IT IS CHANGED !
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XG(16,7),DG(16,7),YG(16),EG(16),X(1),W(1),NI(7),NG(7)
      DATA NBEG,NA,NI,NG /9*0,4,6,8,12,16,24,32/
      DATA XG/ .43056815579702629D 0, .16999052179242813D 0, 14*0.D0,
     X         .46623475710157601D 0, .33060469323313226D 0,
     Y         .11930959304159845D 0,                        13*0.D0,
     Z         .48014492824876812D 0, .39833323870681337D 0,
     1         .26276620495816449D 0, .9171732124782490 D-1, 12*0.D0,
     2         .49078031712335963D 0, .45205862818523743D 0,
     3         .38495133709715234D 0, .29365897714330872D 0,
     4         .18391574949909010D 0, .62616704255734458D-1, 10*0.D0,
     5         .49470046749582497D 0, .47228751153661629D 0,
     6         .43281560119391587D 0, .37770220417750152D 0,
     7         .30893812220132187D 0, .22900838882861369D 0,
     8         .14080177538962946D 0, .47506254918818720D-1,  8*0.D0,
     9         .49759360999851068D 0, .48736427798565475D 0,
     A         .46913727600136638D 0, .44320776350220052D 0,
     B         .41000099298695146D 0, .37006209578927718D 0,
     C         .32404682596848778D 0, .27271073569441977D 0,
     D         .21689675381302257D 0, .15752133984808169D 0,
     E         .95559433736808150D-1, .32028446431302813D-1, 20*0.D0/
      DATA YG/ .49863193092474078D 0, .49280575577263417D 0,
     G         .48238112779375322D 0, .46745303796886984D 0,
     H         .44816057788302606D 0, .42468380686628499D 0,
     I         .39724189798397120D 0, .36609105937014484D 0,
     J         .33152213346510760D 0, .29385787862038116D 0,
     K         .25344995446611470D 0, .21067563806531767D 0,
     L         .16593430114106382D 0, .11964368112606854D 0,
     M         .72235980791398250D-1, .24153832843869158D-1/
      DATA DG/ .17392742256872693D 0, .32607257743127307D 0, 14*0.D0,
     X         .85662246189585173D-1, .18038078652406930D 0,
     Y         .23395696728634552D 0,                        13*0.D0,
     Z         .50614268145188130D-1, .11119051722668724D 0,
     1         .15685332293894364D 0, .18134189168918099D 0, 12*0.D0,
     2         .23587668193255914D-1, .53469662997659215D-1,
     3         .8003916427167311 D-1, .10158371336153296D 0,
     4         .11674626826917740D 0, .12457352290670139D 0, 10*0.D0,
     5         .13576229705877047D-1, .31126761969323946D-1,
     6         .47579255841246392D-1, .62314485627766936D-1,
     7         .74797994408288370D-1, .84578259697501270D-1,
     8         .91301707522461790D-1, .94725305227534250D-1,  8*0.D0,
     9         .61706148999935998D-2, .14265694314466832D-1,
     A         .22138719408709903D-1, .29649292457718390D-1,
     B         .36673240705540153D-1, .43095080765976638D-1,
     C         .48809326052056944D-1, .53722135057982817D-1,
     D         .57752834026862801D-1, .60835236463901696D-1,
     E         .62918728173414148D-1, .63969097673376078D-1, 20*0.D0/
      DATA EG/ .35093050047350483D-2, .8137197365452835 D-2,
     G         .12696032654631030D-1, .17136931456510717D-1,
     H         .21417949011113340D-1, .25499029631188088D-1,
     I         .29342046739267774D-1, .32911111388180923D-1,
     J         .36172897054424253D-1, .39096947893535153D-1,
     K         .41655962113473378D-1, .43826046502201906D-1,
     L         .45586939347881942D-1, .46922199540402283D-1,
     M         .47819360039637430D-1, .48270044257363900D-1/
C
      IF(NBEG.EQ.0) THEN
      NBEG=1
      DO 10 I=1,16
      XG(I,7)=YG(I)
   10 DG(I,7)=EG(I)
      END IF
C
C     N MUST BE EVEN AND >= 4
C
      NN=(N/2)*2
      IF(NN.LT.4) NN=4
      IF(NN.NE.N) THEN
      WRITE (*,*) N,' GAUSS-POINTS WERE NOT POSSIBLE'
      N=NN
      WRITE (*,*) ' INSTEAD NOW ',N,' POINTS ARE USED'
      END IF
      IF(NA.NE.N) THEN
      NA=N
      NR=NA
      DO 11 L=7,1,-1
      NI(L)=NR/NG(L)
      NR=NR-NG(L)*NI(L)
      IF(NR.EQ.2) THEN
      NI(L)=NI(L)-1
      NR=NR+NG(L)
      END IF
 11   CONTINUE
      IF(NR.NE.0) WRITE (*,*) 'ERROR IN GAUSS: NR=',NR
      END IF
C
      DELP=(B-A)/N
      IF(DELP.EQ.0.D0) GO TO 15
      A1=A
      I0=0
      DO 14 L1=7,1,-1
      NIN=NI(L1)
      IF(NIN.EQ.0) GO TO 14
      DEL=DELP*NG(L1)
      M   = NG(L1)/2
      DO 13 K=1,NIN
      ABM =A1+DEL*0.5D0
      DO 12 J=1,M
      I   = M+J
      L   = M+1-J
      J1=J+I0
      I1=I+I0
      X(J1)=ABM-DEL*XG(J,L1)
      X(I1)=ABM+DEL*XG(L,L1)
      W(J1)=    DEL*DG(J,L1)
   12 W(I1)=    DEL*DG(L,L1)
      I0=I0+NG(L1)
   13 A1=A1+DEL
   14 CONTINUE
C
C     TEST
C
      IF(I0.NE.N) WRITE (*,*) 'ERROR IN GAUSS :',I0,'.NE.',N,
     +                        ' A1=',A1,' B=',B
      RETURN
   15 DO 16 L=1,N
      X(L)=A
   16 W(L)=0.D0
      RETURN
      END

c -------------------------------------------------------------------- c
c ---------------- Help functions for radiative decays --------------- c
c -------------------------------------------------------------------- c

c -------- Integrals needed in the radiative gluino decays:  --------- c
c -------- gluino -> neutralino_i gluon.                     --------- c
c -------- Taken from: Haber/Wyler Nucl.Phys.B323 (1989) 267 --------- c

      complex*16 function NS_iint(mj,mi,m,mc)

      implicit double precision (a-h,k-z)

      complex*16 NS_ccspen,lami,lamj,ctmp1,ctmp2,ctmp3,ctmp4,tmpa,tmpb,
     .           m2s,mc2s

      external NS_ccspen

      eps = 1.D-10

      m2s  = m**2*(1.D0-dcmplx(0.D0,1.D0)*eps)
      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4.D0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4.D0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = NS_ccspen( (mj**2+m2s-mc2s+lamj)/(2.D0*m2s) )
      ctmp2 = NS_ccspen( (mj**2+m2s-mc2s-lamj)/(2.D0*m2s) )
      ctmp3 = NS_ccspen( (mi**2+m2s-mc2s+lami)/(2.D0*m2s) )
      ctmp4 = NS_ccspen( (mi**2+m2s-mc2s-lami)/(2.D0*m2s) )

      if(mc.gt.1.D4.or.m.gt.1.D4) then
         NS_iint = 1.D0/(mc2s-m2s)*(1.D0-
     .                        mc2s/(mc2s-m2s)*cdlog(mc2s/m2s))
      else
         NS_iint = -1.D0/(mj**2-mi**2)*(ctmp1+ctmp2-ctmp3-ctmp4)
      endif

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_jint(mj,mi,m,mc)

      implicit double precision (a-h,k-z)

      complex*16 lami,lamj,tmpa,tmpb,ctmp1,ctmp2,ctmp3,ctmp4,m2s,mc2s,
     .           NS_iint

      external NS_iint

      eps = 1.D-10

      m2s  = m**2*(1.D0-dcmplx(0.D0,1.D0)*eps)
      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4.D0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4.D0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = dcmplx(m2s+mc2s-mj**2+lamj)
      ctmp2 = dcmplx(m2s+mc2s-mi**2+lami)

      ctmp3 = (cdlog(ctmp1/dcmplx(2.D0*m*mc)))**2
      ctmp4 = (cdlog(ctmp2/dcmplx(2.D0*m*mc)))**2

      if(mc.gt.1.D4.or.m.gt.1.D4) then
         NS_jint = 1.D0/(mc2s-m2s)*cdlog(m2s/mc2s) 
     .        - NS_iint(mj,mi,m,mc) 
      else
         NS_jint = 1.D0/(mj**2-mi**2)*( ctmp3 - ctmp4 )  
     .        - NS_iint(mj,mi,m,mc) 
      endif

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_i2int(mj,mi,m,mc)

      implicit double precision (a-h,k-z)

      complex*16 lami,lamj,tmpa,tmpb,m2s,mc2s,ctmp3,ctmp4,ctmp5,ctmp6

      eps = 1.D-10

      m2s  = m**2*(1.D0-dcmplx(0.D0,1.D0)*eps)
      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4.D0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4.D0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp3 = m2s+mc2s-dcmplx(mj**2)-lamj
      ctmp4 = m2s+mc2s-dcmplx(mj**2)+lamj
      ctmp5 = m2s+mc2s-dcmplx(mi**2)-lami
      ctmp6 = m2s+mc2s-dcmplx(mi**2)+lami

      if(mc.gt.1.D4.or.m.gt.1.D4) then
         NS_i2int = -1.D0/(mc2s-m2s)**2*((m2s+mc2s)/2.D0-
     .        m2s*mc2s/(mc2s-m2s)*cdlog(mc2s/m2s))
      else
         NS_i2int = (mc2s-m2s)/(2.D0*mi**2*mj**2)*cdlog(m2s/mc2s) +
     .        1.D0/(mj**2-mi**2)*(
     .        lamj/(2.D0*mj**2)*cdlog(ctmp3/ctmp4) -
     .        lami/(2.D0*mi**2)*cdlog(ctmp5/ctmp6) )
      endif

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_kint(mj,mi,m,mc)

      implicit double precision (a-h,k-z)
      double precision m,mc,mj,mi
      complex*16 NS_iint,NS_jint,NS_i2int,m2s,mc2s

      external NS_iint,NS_jint,NS_i2int

      eps = 1.D-10

      m2s  = m**2*(1.D0-dcmplx(0.D0,1.D0)*eps)
      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      if(mc.gt.1.D4.or.m.gt.1.D4) then
         NS_kint = 1.D0/2.D0*NS_i2int(mj,mi,m,mc)
      else
         NS_kint = -1.D0/(mj**2-mi**2)*(1.D0+m2s*NS_iint(mj,mi,m,mc)+
     .        mc2s*NS_jint(mj,mi,m,mc)-mj**2*NS_i2int(mj,mi,m,mc))
      endif

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_iinthelp(mj,mi,m,mc)

      implicit double precision (a-h,k-z)
      complex*16 NS_ccspen,lami,lamj,ctmp1,ctmp2,ctmp3,ctmp4,tmpa,tmpb,
     .           m2s,mc2s

      external NS_ccspen

      eps = 1.D-10

      m2s  = m**2*(1.D0-dcmplx(0.D0,1.D0)*eps)
      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4.D0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4.D0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = NS_ccspen( (mj**2+m2s-mc2s+lamj)/(2.D0*m2s) )
      ctmp2 = NS_ccspen( (mj**2+m2s-mc2s-lamj)/(2.D0*m2s) )
      ctmp3 = NS_ccspen( (mi**2+m2s-mc2s+lami)/(2.D0*m2s) )
      ctmp4 = NS_ccspen( (mi**2+m2s-mc2s-lami)/(2.D0*m2s) )

      NS_iinthelp = -1.D0/(mj**2-mi**2)*(ctmp1+ctmp2-ctmp3-ctmp4)
      
      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_jint0(mj,mi,mc)

      implicit double precision (a-h,k-z)
      complex*16 NS_iinthelp

      external NS_iinthelp

      NS_jint0 = NS_iinthelp(mj,mi,mc,0.D0) 

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_i2int0(mj,mi,mc)

      implicit double precision (a-h,k-z)
      complex*16 mc2s

      eps = 1.D-10

      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      if(mc.gt.1.D4) then
         NS_i2int0 = -1.D0/2.D0/mc2s
      else 
         NS_i2int0 = -mc2s/mi**2/mj**2*cdlog(mc2s) + 1.D0/(mj**2-mi**2)*
     .        ((mc2s-mi**2)/mi**2*cdlog(mc2s-mi**2) -
     .        (mc2s-mj**2)/mj**2*cdlog(mc2s-mj**2) )
      endif

      return

      end

c -------------------------------------------------------------------- c

      complex*16 function NS_kint0(mj,mi,mc)

      implicit double precision (a-h,k-z)
      complex*16 NS_jint0,NS_i2int0,mc2s

      external NS_jint0,NS_i2int0

      eps = 1.D-10

      mc2s = mc**2*(1.D0-dcmplx(0.D0,1.D0)*eps)

      if(mc.gt.1.D4) then
         NS_kint0 = 1.D0/2.D0*NS_i2int0(mj,mi,mc)
      else
         NS_kint0 = -1.D0/(mj**2-mi**2)*(1.D0+mc2s*NS_jint0(mj,mi,mc)
     .        -mj**2*NS_i2int0(mj,mi,mc))
      endif

      return

      end

c -------------------------------------------------------------------- c
c -------------- Help functions for the QCD corrections -------------- c
c -------------------------------------------------------------------- c

c -------------------------------------------------------------------- c
c ------------------ A.Bartl et al., hep-ph/9710286 ------------------ c
c -------------------------------------------------------------------- c
c -- The vertex corrections -- c

      double precision function NS_gluonvertex(ami,amj,amv,lamv,amuv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_C1_lam,SD_C2_lam,SD_B02
      complex*16 SD_C0_lam
      
      C1Den = - dreal(SD_C0_lam(ami,amj,amv,lamv)) - 
     .    NS_C1_lam(ami,amj,amv,ami,lamv,amj,amuv,lamv)
      C2Den = SD_C2_lam(ami,amj,amv,ami,lamv,amj,amuv,lamv)

      NS_gluonvertex = SD_B02(ami**2,lamv,ami,amuv**2) +
     .     SD_B02(amj**2,lamv,amj,amuv**2) - 2.D0*
     .     (ami**2+amj**2-amv**2)*
     .     (dreal(SD_C0_lam(ami,amj,amv,lamv)) + C1Den + C2Den)

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gluinoZvertex(ami,amj,amv,lamv,amuv,
     .     mgluino,amq,iq3L,eq,sw,thetasq)

      implicit double precision (a-h,k-z)
      double precision iq3L,NS_C1,SD_C2,SD_B02
      complex*16 SD_C03

      C1Den = - dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq)) - 
     .    NS_C1(ami,amv,amj,mgluino,amq,amq,amuv)
      C2Den = SD_C2(ami,amv,amj,mgluino,amq,amq,amuv)

      NS_gluinoZvertex = iq3L*(2.D0*mgluino**2*
     .     dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq))
     .     + ami**2*C1Den + amj**2*C2Den
     .     +(mgluino**2-amq**2)*(C1Den+C2Den) +
     .     SD_B02(amv**2,amq,amq,amuv**2))*dsin(2.D0*thetasq) +
     .     2.D0*mgluino*amq*(iq3L-2.D0*eq*sw**2)*
     .     ( dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq)) 
     .     + C1Den + C2Den )*dcos(2.D0*thetasq)

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gluinoWvertex(ami,amj,amv,lamv,amuv,
     .     mgluino,amq,amqp,thsq,thsqp,i,j)

      implicit double precision (a-h,k-z)
      complex*16 SD_C03
      dimension r(2,2),rp(2,2)
      double precision NS_C1,SD_C2,SD_B02
      
      r(1,1) = dcos(thsq)
      r(1,2) = dsin(thsq)
      r(2,1) = -dsin(thsq)
      r(2,2) = dcos(thsq)

      rp(1,1) = dcos(thsqp)
      rp(1,2) = dsin(thsqp)
      rp(2,1) = -dsin(thsqp)
      rp(2,2) = dcos(thsqp)

      C1Den = - dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp)) - 
     .    NS_C1(ami,amv,amj,mgluino,amq,amqp,amuv)
      C2Den = SD_C2(ami,amv,amj,mgluino,amq,amqp,amuv)

      NS_gluinoWvertex = mgluino*(
     .     dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp)) +
     .     C1Den + C2Den )*
     .     ( amq*r(i,2)*rp(j,1) + amqp*r(i,1)*rp(j,2) )
     .     -(ami**2*C1Den + amj**2*C2Den + mgluino**2*(2.D0*
     .       dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp)) 
     .      + C1Den + C2Den ) 
     .      + SD_B02(amv**2,amq,amqp,amuv**2) )*r(i,1)*rp(j,1)
     .     -amq*amqp*( C1Den + C2Den )*r(i,2)*rp(j,2)

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_wavefuncvertex(amsqi,amsqpj,amq,amqp,
     .     thetasq,thetasqp,vecindex,quarkindex,ii,jj,mgluino,lamv,amuv)

      implicit double precision (a-h,k-z)
      double precision NS_deltaZnngluon,NS_deltaZnngluino,
     .         NS_deltaZnnpgluino,NS_deltaZnnpsquark

      dimension gvqqp(2,2)
      dimension gztt(2,2),gzbb(2,2),gztautau(2,2)
      dimension gwtb(2,2),gwntau(2,2)
      dimension gmst(2),gmsb(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1


      COMMON/NS_coup19/gztt,gzbb,gztautau
      COMMON/NS_coup20/gwtb,gwntau
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1


      external NS_deltaZnngluon,NS_deltaZnngluino,NS_deltaZnnpgluino,
     .         NS_deltaZnnpsquark

      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      if(ii.eq.1) then
         ik = 2
      elseif(ii.eq.2) then
         ik = 1
      endif

      if(jj.eq.1) then
         il = 2
      elseif(jj.eq.2) then
         il = 1
      endif

      if(quarkindex.eq.1.D0) then
         if(vecindex.eq.1.D0) then
            do i=1,2,1
               do j=1,2,1
                  gvqqp(i,j) = 1.D0/2.D0/cw*gztt(i,j)
                  amsqk  = gmst(ik) 
                  amsqpl = gmst(il)
                  amsq1  = gmst(1)
                  amsq2  = gmst(2)
                  amsqp1 = gmst(1)
                  amsqp2 = gmst(2)
               end do
            end do
         elseif(vecindex.eq.2.D0) then
            do i=1,2,1
               do j=1,2,1
                  gvqqp(i,j) = 1.D0/dsqrt(2.D0)*gwtb(i,j)
                  amsqk  = gmst(ik) 
                  amsqpl = gmsb(il)
                  amsq1  = gmst(1)
                  amsq2  = gmst(2)
                  amsqp1 = gmsb(1)
                  amsqp2 = gmsb(2)
               end do
            end do
         endif
      elseif(quarkindex.eq.2.D0) then
         if(vecindex.eq.1.D0) then
            do i=1,2,1
               do j=1,2,1
                  gvqqp(i,j) = 1.D0/2.D0/cw*gzbb(i,j)
                  amsqk  = gmsb(ik) 
                  amsqpl = gmsb(il)
                  amsq1  = gmsb(1)
                  amsq2  = gmsb(2)
                  amsqp1 = gmsb(1)
                  amsqp2 = gmsb(2)
               end do
            end do
         elseif(vecindex.eq.2.D0) then
            do i=1,2,1
               do j=1,2,1
                  gvqqp(i,j) = 1.D0/dsqrt(2.D0)*gwtb(j,i)
                  amsqk  = gmsb(ik) 
                  amsqpl = gmst(il)
                  amsq1  = gmsb(1)
                  amsq2  = gmsb(2)
                  amsqp1 = gmst(1)
                  amsqp2 = gmst(2)
               end do
            end do
         endif
      endif

      NS_wavefuncvertex =1.D0/2.D0*(NS_deltaZnngluon(amsqi,lamv,amuv)+
     .                NS_deltaZnngluon(amsqpj,lamv,amuv))*gvqqp(ii,jj)
     .  +  1.D0/2.D0*(NS_deltaZnngluino(amsqi,mgluino,amq,amuv,dble(ii),
     .                thetasq) +
     .                NS_deltaZnngluino(amsqpj,mgluino,amqp,amuv,
     .                dble(jj),thetasqp) )*gvqqp(ii,jj)
     .    +NS_deltaZnnpgluino(amsqi,amsqk,mgluino,amq,amuv,thetasq)*
     .     gvqqp(ik,jj)
     .    +NS_deltaZnnpgluino(amsqpj,amsqpl,mgluino,amqp,amuv,thetasqp)*
     .     gvqqp(ii,il)
     .    +NS_deltaZnnpsquark(amsq1,amsq2,amsqi,amsqk,amuv,thetasq)*
     .     gvqqp(ik,jj)
     .    +NS_deltaZnnpsquark(amsqp1,amsqp2,amsqpj,amsqpl,amuv,thetasqp)
     .     *gvqqp(ii,il)

c      write(*,*)'NS_wavefuncvertex',NS_wavefuncvertex
      return

      end
c -------------------------------------------------------------------- c

      double precision function NS_quarkmixZ(amsq,thetasq,iq3L,eq,amsq1,
     .     amsq2,amq,mgluino,amuv)

      implicit double precision (a-h,k-z)
      double precision iq3L,NS_A01,SD_B02

      COMMON/NS_weinberg/sw,cw,tw
     
      deltathetasqsq = 1.D0/6.D0*dsin(4.D0*thetasq)/(amsq1**2-amsq2**2)
     .     *( NS_A01(amsq2**2,amuv**2) - NS_A01(amsq1**2,amuv**2) )

      v11 = 4.D0*(iq3L*dcos(thetasq)**2-eq*sw**2)
      v22 = 4.D0*(iq3L*dsin(thetasq)**2-eq*sw**2)

      deltathetasqgl = 1.D0/3.D0*mgluino*amq/iq3L/(amsq1**2-amsq2**2)*
     .     ( SD_B02(amsq2**2,mgluino,amq,amuv**2)*v11 -
     .       SD_B02(amsq1**2,mgluino,amq,amuv**2)*v22 )

      NS_quarkmixZ = -1.D0/cw*iq3L*dcos(2.D0*thetasq)*(
     .     deltathetasqsq + deltathetasqgl )

c      write(*,*)'NS_quarkmixZ',NS_quarkmixZ
      return

      end

c -------------------------------------------------------------------- c


      double precision function NS_quarkmixW(amsq,thetasq,iq3L,eq,amsq1,
     .     amsq2,amq,amsqp,thetasqp,iqp3L,eqp,amsqp1,amsqp2,amqp,ii,jj,
     .     mgluino,amuv)

      implicit double precision (a-h,k-z)
      double precision iq3L,iqp3L,NS_A01,SD_B02
      dimension cijwsq(2,2),cijwsqp(2,2)
      COMMON/NS_weinberg/sw,cw,tw
      

      cijwsq(1,1) = -dsin(thetasq)*dcos(thetasqp)
      cijwsq(1,2) =  dsin(thetasq)*dsin(thetasqp)
      cijwsq(2,1) = -dcos(thetasq)*dcos(thetasqp)
      cijwsq(2,2) =  dcos(thetasq)*dsin(thetasqp)

      cijwsqp(1,1) = -dsin(thetasqp)*dcos(thetasq)
      cijwsqp(1,2) = -dcos(thetasq)*dcos(thetasqp)
      cijwsqp(2,1) =  dsin(thetasq)*dsin(thetasqp)
      cijwsqp(2,2) =  dcos(thetasqp)*dsin(thetasq)

c --------------------------------------

      deltathetasqsq = 1.D0/6.D0*dsin(4.D0*thetasq)/(amsq1**2-amsq2**2)
     .     *( NS_A01(amsq2**2,amuv**2) - NS_A01(amsq1**2,amuv**2) )

      deltathetasqpsq = 1.D0/6.D0*dsin(4.D0*thetasqp)
     .     /(amsqp1**2-amsqp2**2)*
     .     ( NS_A01(amsqp2**2,amuv**2) - NS_A01(amsqp1**2,amuv**2) )

c --------------------------------------

      v11sq = 4.D0*(iq3L*dcos(thetasq)**2-eq*sw**2)
      v22sq = 4.D0*(iq3L*dsin(thetasq)**2-eq*sw**2)

      deltathetasqgl = 1.D0/3.D0*mgluino*amq/iq3L/(amsq1**2-amsq2**2)*
     .     ( SD_B02(amsq2**2,mgluino,amq,amuv**2)*v11sq -
     .       SD_B02(amsq1**2,mgluino,amq,amuv**2)*v22sq )

c --------------------------------------

      v11sqp = 4.D0*(iqp3L*dcos(thetasqp)**2-eqp*sw**2)
      v22sqp = 4.D0*(iqp3L*dsin(thetasqp)**2-eqp*sw**2)

      deltathetasqpgl = 1.D0/3.D0*mgluino*amqp/iqp3L
     .     /(amsqp1**2-amsqp2**2)*
     .     ( SD_B02(amsqp2**2,mgluino,amqp,amuv**2)*v11sqp -
     .       SD_B02(amsqp1**2,mgluino,amqp,amuv**2)*v22sqp )

c --------------------------------------

      NS_quarkmixW = 1.D0/dsqrt(2.D0)*( 
     .     cijwsq(ii,jj)*(deltathetasqsq+deltathetasqgl) +
     .     cijwsqp(ii,jj)*(deltathetasqpsq+deltathetasqpgl) )

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_realgluonem(amsqi,amsqpj,amv,lamv)

      implicit double precision (a-h,k-z)

      complex*16 NS_ccspen,NS_kappa

      external NS_kappa,NS_ccspen

      kap = dreal(NS_kappa(amsqi**2,amsqpj**2,amv**2,0.D0))

      b0 = (amsqi**2-amsqpj**2-amv**2+kap)/2.D0/amsqpj/amv
      b1 = (amsqi**2-amsqpj**2+amv**2-kap)/2.D0/amsqi/amv
      b2 = (amsqi**2+amsqpj**2-amv**2-kap)/2.D0/amsqi/amsqpj

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hint = 1.D0/4.D0/amsqi**2*(kap/2.D0*(amsqi**2+amsqpj**2+amv**2)+
     .     2.D0*amsqi**2*amsqpj**2*lb2 + 
     .     2.D0*amsqi**2*amv**2*lb1 +
     .     2.D0*amsqpj**2*amv**2*lb0 )

      hint0 = 1.D0/4.D0/amsqi**2*(-2.D0*amsqpj**2*lb2-2.D0*amv**2*lb1-
     .     kap)

      hint1 = 1.D0/4.D0/amsqi**2*(-2.D0*amsqi**2*lb2-2.D0*amv**2*lb0-
     .     kap)

      hint00 = 1.D0/4.D0/amsqi**4*(
     .     kap*dlog(kap**2/(lamv*amsqi*amsqpj*amv)) - kap -
     .     (amsqpj**2-amv**2)*lb12 - amsqi**2*lb0 )

      hint11 = 1.D0/4.D0/amsqi**2/amsqpj**2*(
     .     kap*dlog(kap**2/(lamv*amsqi*amsqpj*amv)) - kap -
     .     (amsqi**2-amv**2)*lb02 - amsqpj**2*lb1 )

      hint01 = dreal(1.D0/4.D0/amsqi**2*(
     .     -2.D0*dlog((lamv*amsqi*amsqpj*amv)/kap**2)*lb2 +
     .     2.D0*lb2**2 - lb0**2 - lb1**2 + 
     .     2.D0*NS_ccspen(dcmplx(1.D0-b2**2)) - 
     .     NS_ccspen(dcmplx(1-b0**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )

      NS_realgluonem = 2.D0*hint - kap**2/amv**2*( hint0+hint1+
     .     amsqi**2*hint00+amsqpj**2*hint11+
     .     (amsqi**2+amsqpj**2-amv**2)*hint01 )

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltaZnngluon(amsq,lamv,amuv)

      implicit double precision (a-h,k-z)
      double precision SD_B02,SD_BP02

      NS_deltaZnngluon = 2.D0/3.D0*(
     .     SD_B02(amsq**2,0.D0,amsq,amuv**2) + 
     .     2.D0*amsq**2*SD_BP02(amsq**2,lamv,amsq,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltaZnngluino(amsq,mgluino,amq,
     .     amuv,n,thetasq)

      implicit double precision (a-h,k-z)
      double precision SD_B02,SD_BP02

      NS_deltaZnngluino = -2.D0/3.D0*(
     .     SD_B02(amsq**2,mgluino,amq,amuv**2) + 
     .     (amsq**2-amq**2-mgluino**2)*
     .     SD_BP02(amsq**2,mgluino,amq,amuv**2) -
     .     2.D0*amq*mgluino*(-1.D0)**n*dsin(2.D0*thetasq)*
     .     SD_BP02(amsq**2,mgluino,amq,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltaZnnpgluino(amsq,amsqp,mgluino,
     .     amq,amuv,thetasq)

      implicit double precision (a-h,k-z)
      double precision SD_B02

      NS_deltaZnnpgluino = - 1.D0/(amsq**2-amsqp**2)*4.D0/3.D0*
     .     mgluino*amq*dcos(2.D0*thetasq)*
     .     SD_B02(amsq**2,mgluino,amq,amuv**2)

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltaZnnpsquark(amsq1,amsq2,amsq,
     .     amsqp,amuv,thetasq)

      implicit double precision (a-h,k-z)
      double precision NS_A01

      NS_deltaZnnpsquark = - 1.D0/(amsq**2-amsqp**2)*1.D0/6.D0*
     .     dsin(4.D0*thetasq)*( NS_A01(amsq2**2,amuv**2) -
     .     NS_A01(amsq1**2,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c
c ----- A.Arhrib, A.Djouadi, W.Hollik, C.Juenger, hep-ph/9702426 ----- c
c -------------------------------------------------------------------- c
c -- Virtual corrections for the decays squark_i -> Higgs squark_j' -- c

      double precision function NS_gvirtgl(ami,amhi,amj,lamv,amuv)

      implicit double precision (a-h,k-z)
      double precision SD_B02
      complex*16 SD_C0_lam
      
      amisq = ami**2
      amjsq = amj**2
      amhsq = amhi**2
      amusq = amuv**2

      NS_gvirtgl = SD_B02(amisq,lamv,ami,amusq) + 
     .             SD_B02(amjsq,lamv,amj,amusq) -
     .             SD_B02(amhsq,ami,amj,amusq)  +
     .             2.D0*(amisq+amjsq-amhsq)*
     .             dreal(SD_C0_lam(amj,ami,amhi,lamv))

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gvirtmix(am1,am2,amij,amgl,amq,theq,
     .                                   amuv)

      implicit double precision (a-h,k-z)
      double precision NS_A01,SD_B02

      c2q = dcos(2.D0*theq)
      s2q = dsin(2.D0*theq)

      am1sq  = am1**2
      am2sq  = am2**2
      amijsq = amij**2
      amusq  = amuv**2

      NS_gvirtmix = c2q*s2q*(NS_A01(am2sq,amusq)-
     .           NS_A01(am1sq,amusq))+
     .           4.D0*c2q*amq*amgl*SD_B02(amijsq,amgl,amq,amusq)

      return 

      end

c -------------------------------------------------------------------- c

      double precision function NS_gvirtmixdiv(am1,am2,amij,amgl,amq,
     .                                         theq,amuv)

      implicit double precision (a-h,k-z)
      double precision SD_B02_DIV

      c2q = dcos(2.D0*theq)
      s2q = dsin(2.D0*theq)

      am1sq  = am1**2
      am2sq  = am2**2
      amijsq = amij**2
      amusq  = amuv**2

      NS_gvirtmixdiv = c2q*s2q*(am2sq*dlog(amusq)-am1sq*dlog(amusq) )
     .     +4.D0*c2q*amq*amgl*SD_B02_DIV(amijsq,amgl,amq,amusq)

      return 

      end

c -------------------------------------------------------------------- c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
c -------------------------------------------------------------------- c
C REF:::PRD57,5860
c -------------------------------------------------------------------- c
C**********************************************************************C
C**********************************************************************C
      DOUBLE PRECISION function NS_topneut1719(nh,amuv)
*
      IMPLICIT NONE
*   
      INTEGER nh
      COMPLEX*16 SD_C03
*
      DOUBLE PRECISION sm(2,2,2)
      DOUBLE PRECISION amusq,amuv,aml,amh,amnh,ama,amna,
     .coupphi,squarktopneut,amq,amar,amch,
     .ast1sq,ast2sq,amlsq,amhsq,amnhsq,amasq,amnasq,v1,v2,a1,a2,gluinoex
      DOUBLE PRECISION sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION runmt,runmb,rmtauc
*        
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1

      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_phitoptop/Httr,Attr
*
      call NS_smatrix(sm)
*
      aml = SMASS(1)
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = PMASS(1)
      amar = ama
      amna = PMASS(2)
      amch = CMASS
      amusq  = amuv**2
      amq    = runmt
      ast1sq = ast1**2
      ast2sq = ast2**2
      amlsq  = aml**2
      amhsq  = amh**2
      amnhsq = amnh**2
      amasq  = ama**2
      amnasq = amna**2
*
      v1 = 1.D0/2.D0*(ct-st)
      v2 = -1.D0/2.D0*(ct+st)
      a1 = -v2
      a2 = v1
C------------------EQ.17------------------------
      if(nh.eq.1) then
         coupphi = sm(2,1,1)*Hstopstopr(1,1,1)*sm(1,1,1)*
     .        SD_B02(aml**2,ast1,ast1,amusq) +
     .        sm(2,1,1)*Hstopstopr(1,1,2)*sm(2,1,1)*
     .        SD_B02(aml**2,ast1,ast2,amusq) +
     .        sm(2,2,1)*Hstopstopr(1,1,2)*sm(1,1,1)*
     .        SD_B02(aml**2,ast2,ast1,amusq) +
     .        sm(2,2,1)*Hstopstopr(1,2,2)*sm(2,1,1)*
     .        SD_B02(aml**2,ast2,ast2,amusq) 
      elseif(nh.eq.2) then
         coupphi = sm(2,1,1)*Hstopstopr(2,1,1)*sm(1,1,1)*
     .        SD_B02(amh**2,ast1,ast1,amusq) +
     .        sm(2,1,1)*Hstopstopr(2,1,2)*sm(2,1,1)*
     .        SD_B02(amh**2,ast1,ast2,amusq) +
     .        sm(2,2,1)*Hstopstopr(2,1,2)*sm(1,1,1)*
     .        SD_B02(amh**2,ast2,ast1,amusq) +
     .        sm(2,2,1)*Hstopstopr(2,2,2)*sm(2,1,1)*
     .        SD_B02(amh**2,ast2,ast2,amusq) 
       elseif(nh.eq.3) then
         coupphi = sm(2,1,1)*Hstopstopr(3,1,1)*sm(1,1,1)*
     .        SD_B02(amnh**2,ast1,ast1,amusq) +
     .        sm(2,1,1)*Hstopstopr(3,1,2)*sm(2,1,1)*
     .        SD_B02(amnh**2,ast1,ast2,amusq) +
     .        sm(2,2,1)*Hstopstopr(3,1,2)*sm(1,1,1)*
     .        SD_B02(amnh**2,ast2,ast1,amusq) +
     .        sm(2,2,1)*Hstopstopr(3,2,2)*sm(2,1,1)*
     .        SD_B02(amnh**2,ast2,ast2,amusq) 
       elseif(nh.eq.4) then
         coupphi = sm(2,1,1)*Astopstopr(1,1,1)*sm(1,1,1)*
     .        SD_B02(ama**2,ast1,ast2,amusq) -
     .        sm(2,1,1)*Astopstopr(1,1,2)*sm(2,1,1)*
     .        SD_B02(ama**2,ast1,ast2,amusq)-
     .        sm(2,2,1)*Astopstopr(1,1,2)*sm(1,1,1)*
     .        SD_B02(ama**2,ast2,ast1,amusq)-
     .        sm(2,2,1)*Astopstopr(1,2,2)*sm(2,1,1)*
     .        SD_B02(ama**2,ast2,ast2,amusq)
      elseif(nh.eq.5) then
         coupphi = sm(2,1,1)*Astopstopr(2,1,1)*sm(1,1,1)*
     .        SD_B02(amna**2,ast1,ast2,amusq) -
     .        sm(2,1,1)*Astopstopr(2,1,2)*sm(2,1,1)*
     .        SD_B02(amna**2,ast1,ast2,amusq)-
     .        sm(2,2,1)*Astopstopr(2,1,2)*sm(1,1,1)*
     .        SD_B02(amna**2,ast2,ast1,amusq)-
     .        sm(2,2,1)*Astopstopr(2,2,2)*sm(2,1,1)*
     .        SD_B02(amna**2,ast2,ast2,amusq) 
      endif
      
      coupphi = dsqrt(2.D0)*AMZ**2*coupphi
      
c the result Eq.(17) in the paper

            squarktopneut = coupphi
c --------------------------------------------
c the result Eq.(19) in the paper

      if(nh.eq.1) then
         gluinoex = 4.D0*AMW*Httr(1)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(aml**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(aml**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(ast2sq,amlsq,ast1sq,mgluino,amq,amq)) )
      elseif(nh.eq.2) then
         gluinoex = 4.D0*AMW*Httr(2)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amh**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(ast2sq,amhsq,ast1sq,mgluino,amq,amq)) )
          elseif(nh.eq.3) then
         gluinoex = 4.D0*AMW*Httr(3)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amnh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amnh**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(ast2sq,amnhsq,ast1sq,mgluino,amq,amq)) )
         elseif(nh.eq.4) then
         gluinoex = -4.D0*AMW*Attr(1)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*ama**2
     .         +(v2*a1+a2*v1)*(ast1**2*amq-ast2**2*amq) )*
     .        dreal(SD_C03(ast2sq,amasq,ast1sq,mgluino,amq,amq)) )
      elseif(nh.eq.5) then
         gluinoex = -4.D0*AMW*Attr(2)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*amna**2
     .         +(v2*a1+a2*v1)*(ast1**2*amq-ast2**2*amq) )*
     .        dreal(SD_C03(ast2sq,amnasq,ast1sq,mgluino,amq,amq)) )
      endif

      NS_topneut1719 = squarktopneut + gluinoex

      return
      
      end

c -------------------------------------------------------------------- c
      
      DOUBLE PRECISION function NS_botneut1719(nh,amuv)
      IMPLICIT NONE
*
      INTEGER nh
      COMPLEX*16 SD_C03
*     
      DOUBLE PRECISION sm(2,2,2)
      DOUBLE PRECISION amusq,amuv,aml,amh,amnh,ama,amna,
     .coupphi,squarkbotneut,amq,amar,amch,
     .asb1sq,asb2sq,amlsq,amhsq,amnhsq,amasq,amnasq,v1,v2,a1,a2,gluinoex
      DOUBLE PRECISION sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
***      
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_phibotbot/Hbbr,Abbr
*
      aml = SMASS(1)
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = PMASS(1)
      amar = ama
      amna = PMASS(2)
      amch = CMASS
      call NS_smatrix(sm)
      amusq  = amuv**2
c -----------------------------

      if(nh.eq.1) then
         coupphi = sm(2,1,2)*Hsbotsbotr(1,1,1)*sm(1,1,2)*
     .             SD_B02(aml**2,asb1,asb1,amusq) +
     .             sm(2,1,2)*Hsbotsbotr(1,1,2)*sm(2,1,2)*
     .             SD_B02(aml**2,asb1,asb2,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(1,1,2)*sm(1,1,2)*
     .             SD_B02(aml**2,asb2,asb1,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(1,2,2)*sm(2,1,2)*
     .             SD_B02(aml**2,asb2,asb2,amusq) 
      elseif(nh.eq.2) then
         coupphi = sm(2,1,2)*Hsbotsbotr(2,1,1)*sm(1,1,2)*
     .             SD_B02(amh**2,asb1,asb1,amusq) +
     .             sm(2,1,2)*Hsbotsbotr(2,1,2)*sm(2,1,2)*
     .             SD_B02(amh**2,asb1,asb2,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(2,1,2)*sm(1,1,2)*
     .             SD_B02(amh**2,asb2,asb1,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(2,2,2)*sm(2,1,2)*
     .             SD_B02(amh**2,asb2,asb2,amusq) 
      elseif(nh.eq.3) then
         coupphi = sm(2,1,2)*Hsbotsbotr(3,1,1)*sm(1,1,2)*
     .             SD_B02(amnh**2,asb1,asb1,amusq) +
     .             sm(2,1,2)*Hsbotsbotr(3,1,2)*sm(2,1,2)*
     .             SD_B02(amnh**2,asb1,asb2,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(3,1,2)*sm(1,1,2)*
     .             SD_B02(amnh**2,asb2,asb1,amusq) +
     .             sm(2,2,2)*Hsbotsbotr(3,2,2)*sm(2,1,2)*
     .             SD_B02(amnh**2,asb2,asb2,amusq) 
      elseif(nh.eq.4) then
         coupphi = sm(2,1,2)*Asbotsbotr(1,1,1)*sm(1,1,2)*
     .             SD_B02(ama**2,asb1,asb2,amusq) -
     .              sm(2,1,2)*Asbotsbotr(1,1,2)*sm(2,1,2)*
     .             SD_B02(ama**2,asb1,asb2,amusq)-
     .             sm(2,2,2)*Asbotsbotr(1,1,2)*sm(1,1,2)*
     .             SD_B02(ama**2,asb2,asb1,amusq)-
     .             sm(2,2,2)*Asbotsbotr(1,2,2)*sm(2,1,2)*
     .             SD_B02(ama**2,asb2,asb1,amusq)
      elseif(nh.eq.5) then
         coupphi = sm(2,1,2)*Asbotsbotr(2,1,1)*sm(1,1,2)*
     .             SD_B02(amna**2,asb1,asb2,amusq) -
     .              sm(2,1,2)*Asbotsbotr(2,1,2)*sm(2,1,2)*
     .             SD_B02(amna**2,asb1,asb2,amusq)-
     .             sm(2,2,2)*Asbotsbotr(2,1,2)*sm(1,1,2)*
     .             SD_B02(amna**2,asb2,asb1,amusq)-
     .             sm(2,2,2)*Asbotsbotr(2,2,2)*sm(2,1,2)*
     .             SD_B02(amna**2,asb2,asb1,amusq)

      endif

      coupphi = dsqrt(2.D0)*AMZ**2*coupphi

c the result Eq.(17) in the paper

      squarkbotneut = coupphi
c --------------------------------------------
      amq    = runmb
      asb1sq = asb1**2
      asb2sq = asb2**2
      amlsq  = aml**2
      amhsq  = amh**2
      amnhsq  = amnh**2
      amasq  = ama**2
      amnasq  = amna**2
**
      v1 = 1.D0/2.D0*(cb-sb)
      v2 = -1.D0/2.D0*(cb+sb)
      a1 = -v2
      a2 = v1

c the result Eq.(19) in the paper

      if(nh.eq.1) then
         gluinoex = 4.D0*AMW*Hbbr(1)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(aml**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(aml**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(asb2sq,amlsq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.2) then
         gluinoex = 4.D0*AMW*Hbbr(2)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amh**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(asb2sq,amhsq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.3) then
         gluinoex = 4.D0*AMW*Hbbr(3)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2.D0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amnh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amnh**2-4.D0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2.D0*amq) )*
     .        dreal(SD_C03(asb2sq,amnhsq,asb1sq,mgluino,amq,amq)) )
               elseif(nh.eq.4) then
         gluinoex = -4.D0*AMW*Abbr(1)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*ama**2
     .         +(v2*a1+a2*v1)*(asb1**2*amq-asb2**2*amq) )*
     .        dreal(SD_C03(asb2sq,amasq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.5) then
         gluinoex = -4.D0*AMW*Abbr(2)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*amna**2
     .         +(v2*a1+a2*v1)*(asb1**2*amq-asb2**2*amq) )*
     .        dreal(SD_C03(asb2sq,amnasq,asb1sq,mgluino,amq,amq)) )
      endif

      NS_botneut1719 = squarkbotneut + gluinoex

      return

      end
c -------------------------------------------------------------------- c

      DOUBLE PRECISION function NS_stopsbot1719(amuv,ni,nj)
*
      IMPLICIT NONE 
      COMPLEX*16 SD_C03
*
      DOUBLE PRECISION sm(2,2,2),gctbr(2,2),gmst(2),vq(2),aq(2),vqp(2),
     .          aqp(2),gmsb(2)
      INTEGER ni,nj
      DOUBLE PRECISION amusq,amuv,
     .coupphi,squarkstopsbot,gluinoex,amq,amqp,amchsq,vs,as
      DOUBLE PRECISION sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION ama,aml,amh,amch,amar,amna,amnh
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SD_B02
*
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
**
      aml = SMASS(1)
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = PMASS(1)
      amar = ama
      amna = PMASS(2)
      amch = CMASS
*      
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2
*
      call NS_smatrix(sm)
      amusq = amuv**2
c  --------------

      coupphi = sm(ni,1,1)*gctbr(1,1)*sm(1,nj,2)*
     .     SD_B02(amch**2,ast1,asb1,amusq) +
     .     sm(ni,1,1)*gctbr(1,2)*sm(2,nj,2)*
     .     SD_B02(amch**2,ast1,asb2,amusq) +
     .     sm(ni,2,1)*gctbr(2,1)*sm(1,nj,2)*
     .     SD_B02(amch**2,ast2,asb1,amusq) +
     .     sm(ni,2,1)*gctbr(2,2)*sm(2,nj,2)*
     .     SD_B02(amch**2,ast2,asb2,amusq) 

      coupphi = dsqrt(2.D0)*AMW**2*coupphi

c the result Eq.(17) in the paper

      squarkstopsbot = coupphi

c --------------------------------------------

      amq  =  runmt
      amqp =  runmb

      amchsq = amch**2

      vq(1) = 1.D0/2.D0*(ct-st)
      vq(2) = -1.D0/2.D0*(ct+st)
      aq(1) = -vq(2)
      aq(2) = vq(1)

      vqp(1) = 1.D0/2.D0*(cb-sb)
      vqp(2) = -1.D0/2.D0*(cb+sb)
      aqp(1) = -vqp(2)
      aqp(2) = vqp(1)

      vs = 2.D0*dsqrt(2.D0)*AMW*(chtbrunr+chtbrunl)
      as = 2.D0*dsqrt(2.D0)*AMW*(chtbrunr-chtbrunl)

c the result Eq.(19) in the paper

      gluinoex = (vs*(amq*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj))) 
     .     -as*(amq*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))+
     .     mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))))*
     .     SD_B02(gmst(ni)**2,mgluino,amq,amusq) +
     .     (vs*(amqp*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj))) +
     .     as*(amqp*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj)) -
     .     mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))))*
     .     SD_B02(gmsb(nj)**2,mgluino,amqp,amusq) +
     .     (vs*(amq*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     amqp*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))) -
     .     as*(amq*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))-
     .     amqp*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))))*
     .     SD_B02(amchsq,amq,amqp,amusq) +
     .     (as*mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))*
     .      (amch**2-(amq-amqp)**2) -
     .      vs*mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj))*
     .      (amch**2-(amq+amqp)**2) +
     .      as*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))*(gmsb(nj)**2*amq
     .      -gmst(ni)**2*amqp-(mgluino**2-amq*amqp)*(amq-amqp)) -
     .      vs*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))*(gmsb(nj)**2*amq
     .      +gmst(ni)**2*amqp-(mgluino**2+amq*amqp)*(amq+amqp)) )*
     . dreal(SD_C03(gmst(ni)**2,amchsq,gmsb(nj)**2,mgluino,amq,amqp))


      NS_stopsbot1719 = squarkstopsbot + gluinoex
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION function NS_sbotstop1719(amuv,ni,nj)
*
      IMPLICIT NONE
      INTEGER ni,nj
      COMPLEX*16 SD_C03
*
      DOUBLE PRECISION sm(2,2,2),gctbr(2,2),gmsb(2),vq(2),aq(2),vqp(2),
     .          aqp(2),gmst(2)
      DOUBLE PRECISION amusq,coupphi,amuv,
     .gluinoex,amq,amqp,amchsq,vs,as
      DOUBLE PRECISION squarksbotstop
      DOUBLE PRECISION ama,aml,amh,amch,amar,amna,amnh
      DOUBLE PRECISION sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
**
      aml = SMASS(1)
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = PMASS(1)
      amar = ama
      amna = PMASS(2)
      amch = CMASS
*
      gmsb(1) = asb1
      gmsb(2) = asb2
**
      gmst(1) = ast1
      gmst(2) = ast2
*
      call NS_smatrix(sm)
      amusq = amuv**2
**
      coupphi = sm(nj,1,1)*gctbr(1,1)*sm(1,ni,2)*
     .     SD_B02(amch**2,ast1,asb1,amusq) +
     .     sm(nj,1,1)*gctbr(1,2)*sm(2,ni,2)*
     .     SD_B02(amch**2,ast1,asb2,amusq) +
     .     sm(nj,2,1)*gctbr(2,1)*sm(1,ni,2)*
     .     SD_B02(amch**2,ast2,asb1,amusq) +
     .     sm(nj,2,1)*gctbr(2,2)*sm(2,ni,2)*
     .     SD_B02(amch**2,ast2,asb2,amusq) 

      coupphi = dsqrt(2.D0)*AMW**2*coupphi

c the result Eq.(17) in the paper

      squarksbotstop = coupphi

c --------------------------------------------

      amq  = runmt
      amqp = runmb

      amchsq = amch**2

      vqp(1) = 1.D0/2.D0*(cb-sb)
      vqp(2) = -1.D0/2.D0*(cb+sb)
      aqp(1) = -vqp(2)
      aqp(2) = vqp(1)

      vq(1) = 1.D0/2.D0*(ct-st)
      vq(2) = -1.D0/2.D0*(ct+st)
      aq(1) = -vq(2)
      aq(2) = vq(1)

      vs = 2.D0*dsqrt(2.D0)*AMW*(chtbrunr+chtbrunl)
      as = 2.D0*dsqrt(2.D0)*AMW*(chtbrunr-chtbrunl)

c the result Eq.(19) in the paper

      gluinoex = (vs*(amq*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni))) 
     .     -as*(amq*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))+
     .     mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))))*
     .     SD_B02(gmst(nj)**2,mgluino,amq,amusq) +
     .     (vs*(amqp*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni))) +
     .     as*(amqp*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni)) -
     .     mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))))*
     .     SD_B02(gmsb(ni)**2,mgluino,amqp,amusq) +
     .     (vs*(amq*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     amqp*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))) -
     .     as*(amq*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))-
     .     amqp*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))))*
     .     SD_B02(amchsq,amq,amqp,amusq) +
     .     (as*mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))*
     .      (amch**2-(amq-amqp)**2) -
     .      vs*mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni))*
     .      (amch**2-(amq+amqp)**2) +
     .      as*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))*(gmsb(ni)**2*amq
     .      -gmst(nj)**2*amqp-(mgluino**2-amq*amqp)*(amq-amqp)) -
     .      vs*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))*(gmsb(ni)**2*amq
     .      +gmst(nj)**2*amqp-(mgluino**2+amq*amqp)*(amq+amqp)) )*
     . dreal(SD_C03(gmst(nj)**2,amchsq,gmsb(ni)**2,mgluino,amq,amqp))

      NS_sbotstop1719 = squarksbotstop + gluinoex

      return

      end

c -------------------------------------------------------------------- c

      subroutine NS_smatrix(sm)

      implicit double precision (a-h,k-z)

      dimension sm(2,2,2)
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn

      sm(1,1,1) = ct**2-st**2
      sm(1,2,1) = -2.D0*st*ct
      sm(2,1,1) = sm(1,2,1)
      sm(2,2,1) = -sm(1,1,1)

      sm(1,1,2) = cb**2-sb**2
      sm(1,2,2) = -2.D0*sb*cb
      sm(2,1,2) = sm(1,2,2)
      sm(2,2,2) = -sm(1,1,2)

      end
c -------------------------------------------------------------------- c

c -------------------------------------------------------------------- c
c -- the counterterm for squark2 -> h/H/A squark1 --

      DOUBLE PRECISION function NS_dcounterneut(amsq1,amsq2,amq,theq,
     . mgluino,amuv,amuvdiv,lamv,nq,nh)
*      
      IMPLICIT NONE
      INTEGER nq,nh
      DOUBLE PRECISION lamv,amq,theq,amuv,amuvdiv,amqq,amsq1,amsq2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION dtlttr(3,2,2),dtattr(2,2,2),
     .datlttr(3,2,2),datattr(2,2,2),
     .dthlttr(3,2,2),dthattr(2,2,2)
      DOUBLE PRECISION dblbbr(3,2,2),dbabbr(2,2,2),
     .dablbbr(3,2,2),dababbr(2,2,2),
     .dthlbbr(3,2,2),dthabbr(2,2,2)
      DOUBLE PRECISION gqqr(2,2),dmqqr(2,2),daqr(2,2),dthr(2,2)
      DOUBLE PRECISION mgluino
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION NS_deltaz,NS_deltamqdiv,
     .NS_deltaAq,NS_deltathdiv
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/NS_higgsst1st2deriv/
     .dtlttr,dtattr,datlttr,datattr,dthlttr,dthattr
      COMMON/NS_higgssb1sb2deriv/dblbbr,dbabbr,
     .dablbbr,dababbr,dthlbbr,dthabbr

c --- the running mass ---
      
      if(nq.eq.1) then
         amqq = runmt 
      elseif(nq.eq.2) then
         amqq = runmb 
      endif
C
C      write(*,*)'hi',dtlttr(1,2,1),datlttr(1,2,1),dthlttr(1,2,1)
c --- some definitions ---
C ----- FOR OTHER HIGGS--CORRECTIONS ARE ADDED--------------

      if(nq.eq.1.and.nh.eq.1) then
         gqqr(2,1)  = Hstopstopr(1,1,2)
         dmqqr(2,1) = dtlttr(1,2,1)
         daqr(2,1)  = datlttr(1,2,1)
         dthr(2,1)  = dthlttr(1,2,1)
      elseif(nq.eq.1.and.nh.eq.2) then
         gqqr(2,1)  = Hstopstopr(2,1,2)
         dmqqr(2,1) = dtlttr(2,2,1)
         daqr(2,1)  = datlttr(2,2,1)
         dthr(2,1)  = dthlttr(2,2,1)
C --NMSSM
      elseif(nq.eq.1.and.nh.eq.3) then
         gqqr(2,1)  = Hstopstopr(3,1,2) 
         dmqqr(2,1) = dtlttr(3,2,1)
         daqr(2,1)  = datlttr(3,2,1)
         dthr(2,1)  = dthlttr(3,2,1)
      elseif(nq.eq.1.and.nh.eq.4) then
         gqqr(2,1)  = -Astopstopr(1,1,2)
         dmqqr(2,1) = -dtattr(1,2,1)
         daqr(2,1)  = -datattr(1,2,1)
         dthr(2,1)  = -dthattr(1,2,1)
      elseif(nq.eq.1.and.nh.eq.5) then
         gqqr(2,1)  = -Astopstopr(2,1,2)
         dmqqr(2,1) = -dtattr(2,2,1)
         daqr(2,1)  = -datattr(2,2,1)
         dthr(2,1)  = -dthattr(2,2,1)
C --NMSSM
C! NEGATIVE SIGN IN FRON OF dblbbr(1,2,1) AND dablbbr(1,2,1) HAVE BEEN ADDED
C! FOR CONSISTENCY WITH SDECAY FORMULA. SIMILARLY THERE SHOULD NOT BE ANY
C! NEGATIVE SIGN BEFORE dblbbr(2,2,1) AND dablbbr(2,2,1)

      elseif(nq.eq.2.and.nh.eq.1) then
         gqqr(2,1)  = Hsbotsbotr(1,1,2)
         dmqqr(2,1) = -dblbbr(1,2,1)
         daqr(2,1)  = -dablbbr(1,2,1)
         dthr(2,1)  = dthlbbr(1,2,1)
c       write(*,*)'gqqr(2,1)',gqqr(2,1),dmqqr(2,1),daqr(2,1),dthr(2,1)
      elseif(nq.eq.2.and.nh.eq.2) then
         gqqr(2,1)  = Hsbotsbotr(2,1,2)
         dmqqr(2,1) = dblbbr(2,2,1)
         daqr(2,1)  = dablbbr(2,2,1)
         dthr(2,1)  = dthlbbr(2,2,1)
C --NMSSM
      elseif(nq.eq.2.and.nh.eq.3) then
         gqqr(2,1)  = Hsbotsbotr(3,1,2)
         dmqqr(2,1) = dblbbr(3,2,1)
         daqr(2,1)  = dablbbr(3,2,1)
         dthr(2,1)  = dthlbbr(3,2,1)
      elseif(nq.eq.2.and.nh.eq.4) then
         gqqr(2,1)  = -Asbotsbotr(1,1,2)
         dmqqr(2,1) = -dbabbr(1,2,1)
         daqr(2,1)  = -dababbr(1,2,1)
         dthr(2,1)  = -dthabbr(1,2,1) 
c      write(*,*)'gqqr(2,1)',gqqr(2,1),dmqqr(2,1),daqr(2,1),dthr(2,1),
c     .Asbotsbotr(2,1,2)
      elseif(nq.eq.2.and.nh.eq.5) then
         gqqr(2,1)  = -Asbotsbotr(2,1,2)
         dmqqr(2,1) = -dbabbr(2,2,1)
         daqr(2,1)  = -dababbr(2,2,1)
         dthr(2,1)  = -dthabbr(2,2,1)
c         write(*,*)'gqqr(2,1)',gqqr(2,1),dmqqr(2,1),daqr(2,1),dthr(2,1)
C --NMSSM
      endif

      NS_dcounterneut = gqqr(2,1)/2.D0*(
     .     NS_deltaz(amsq2,mgluino,amq,theq,amuv,lamv,2) 
     .     +NS_deltaz(amsq1,mgluino,amq,theq,amuv,lamv,1) )
     .     +dmqqr(2,1)*amqq*
     .     NS_deltamqdiv(amsq1,amsq2,mgluino,amq,theq,amuvdiv,lamv) +
     .     daqr(2,1)*
     .     NS_deltaAq(amsq2,amsq1,amsq2,mgluino,amq,theq,amuvdiv,lamv,
     .                nq) +
     .   dthr(2,1)*NS_deltathdiv(amsq1,amsq2,mgluino,amq,theq,amuvdiv)

      NS_dcounterneut = -dsqrt(2.D0)*amz**2*NS_dcounterneut

      return

      end
c -------------------------------------------------------------------- c

      double precision function NS_deltaz(amsq,mgluino,amq,theq,amuv,
     .                                 lamv,ni)

      implicit double precision (a-h,k-z)
      integer ni
      DOUBLE PRECISION SD_B02,SD_BP02

      amusq = amuv**2

      NS_deltaz = 2.D0*( 
     .     (amq**2+mgluino**2-amsq**2)*
     .     SD_BP02(amsq**2,amq,mgluino,amusq) - 
     .     SD_B02(amsq**2,amq,mgluino,amusq) 
     .     + SD_B02(amsq**2,lamv,amsq,amusq) +
     .     2.D0*(-1.D0)**ni*dsin(2.D0*theq)*mgluino*amq*
     .     SD_BP02(amsq**2,amq,mgluino,amusq) 
     .     +2.D0*amsq**2*
     .     SD_BP02(amsq**2,lamv,amsq,amusq) )

      return 

      end 

c -------------------------------------------------------------------- c

      double precision function NS_deltaAq(amdec,amsq1,amsq2,mgluino,
     .                                  amq,theq,amuv,lamv,nq)

      implicit double precision (a-h,k-z)
      integer nq
      DOUBLE PRECISION NS_deltamqdiv,NS_deltathdiv,NS_deltamsqdiv
     
c --- the running masses ---

      amqq = amq

c --------------------------

      NS_deltaAq = (amsq1**2-amsq2**2)/(2.D0*amqq)*(
     .     2.D0*dcos(2.D0*theq)*
     .     NS_deltathdiv(amsq1,amsq2,mgluino,amq,theq,amuv)
     .     -dsin(2.D0*theq)*
     .     NS_deltamqdiv(amsq1,amsq2,mgluino,amq,theq,amuv,lamv) )
     .     +dsin(2.D0*theq)*
     .     NS_deltamsqdiv(amsq1,mgluino,amq,theq,amuv,lamv,1,nq)/amqq
     .     -dsin(2.D0*theq)*
     .     NS_deltamsqdiv(amsq2,mgluino,amq,theq,amuv,lamv,2,nq)/amqq

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltathdiv(amsq1,amsq2,mgluino,amq,
     .                                        theq,amuv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV

      amusq = amuv**2

      NS_deltathdiv = 1.D0/(amsq1**2-amsq2**2)*(4.D0*mgluino*amq*
     .     dcos(2.D0*theq)/2.D0*(
     .     SD_B02_DIV(amsq1**2,amq,mgluino,amusq)+
     .     SD_B02_DIV(amsq2**2,amq,mgluino,amusq) ) +
     .     dcos(2.D0*theq)*dsin(2.D0*theq)*(amsq2**2*dlog(amusq)
     .     -amsq1**2*log(amusq) ) )

      return

      end


c -------------------------------------------------------------------- c

      double precision function NS_deltamqdiv(amsq1,amsq2,mgluino,amq,
     .                                        theq,amuv,lamv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV,SD_B02_DIV

      amusq = amuv**2

      NS_deltamqdiv = -(2.D0*NS_B1_DIV(amq**2,amq,lamv,amusq) +
     .     4.D0*SD_B02_DIV(amq**2,amq,lamv,amusq) +
     .     NS_B1_DIV(amq**2,mgluino,amsq1,amusq) +
     .     NS_B1_DIV(amq**2,mgluino,amsq2,amusq) ) +
     .     dsin(2.D0*theq)*mgluino/amq*( 
     .     SD_B02_DIV(amq**2,mgluino,amsq1,amusq) -
     .     SD_B02_DIV(amq**2,mgluino,amsq2,amusq) )

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_deltamsqdiv(amsq,mgluino,amq,theq,
     .                                         amuv,lamv,i,nq)

      implicit double precision (a-h,k-z)
      integer nq
      dimension gmsq(2)
      DOUBLE PRECISION NS_A01_DIV,SD_B02_DIV
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1

      amusq = amuv**2
      ip = 3 - i
      
      if(nq.eq.1) then
         gmsq(1) = ast1 
         gmsq(2) = ast2
      elseif(nq.eq.2) then
         gmsq(1) = asb1
         gmsq(2) = asb2
      endif

      NS_deltamsqdiv = -(amq**2+mgluino**2-amsq**2)*
     .     SD_B02_DIV(amsq**2,amq,mgluino,amusq) - 2.D0*amsq**2*
     .     SD_B02_DIV(amsq**2,lamv,amsq,amusq) -
     .     NS_A01_DIV(mgluino**2,amusq) - NS_A01_DIV(amq**2,amusq) +
     .     1.D0/2.D0*((1.D0+dcos(2.D0*theq)**2)*
     .     NS_A01_DIV(amsq**2,amusq)
     .     + dsin(2.D0*theq)**2*NS_A01_DIV(gmsq(ip)**2,amusq) ) -
     .     2.D0*(-1.D0)**i*dsin(2.D0*theq)*mgluino*amq*
     .     SD_B02_DIV(amsq**2,amq,mgluino,amusq) 

      return
      
      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c


      double precision function NS_realcorr(mphi,msq,msqp,lamv,nh,nq,
     .                                 ni,nj,scala)

      IMPLICIT NONE
      integer nh,nq,ni,nj
      DOUBLE PRECISION gctbr(2,2),lamv
      DOUBLE PRECISION msq,mphi,msqp,
     .GPHI,KAP,B0,B1,B2,lb0,lb1,lb2,scala
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      complex*16 NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/NS_hcsbotstop/gctbr

c --- some definitions ---

      if(nh.eq.6.and.nq.eq.0) then
         gphi = -dsqrt(2.D0)*amw**2*gctbr(ni,nj)
      endif

      if(nh.eq.1.and.nq.eq.1) then
         gphi = -dsqrt(2.D0)*amz**2*Hstopstopr(1,1,2)
      elseif(nh.eq.1.and.nq.eq.2) then
         gphi = -dsqrt(2.D0)*amz**2*Hsbotsbotr(1,1,2)
      elseif(nh.eq.2.and.nq.eq.1) then
         gphi = -dsqrt(2.D0)*amz**2*Hstopstopr(2,1,2)
      elseif(nh.eq.2.and.nq.eq.2) then
         gphi = -dsqrt(2.D0)*amz**2*Hsbotsbotr(2,1,2)
c --NMSSM 
      elseif(nh.eq.4.and.nq.eq.1) then
         gphi = dsqrt(2.D0)*amz**2*Astopstopr(1,1,2)
      elseif(nh.eq.4.and.nq.eq.2) then
         gphi = dsqrt(2.D0)*amz**2*Asbotsbotr(1,1,2)
      elseif(nh.eq.3.and.nq.eq.1) then
         gphi = -dsqrt(2.D0)*amz**2*Hstopstopr(3,1,2)
      elseif(nh.eq.3.and.nq.eq.2) then
         gphi = -dsqrt(2.D0)*amz**2*Hsbotsbotr(3,1,2)
      elseif(nh.eq.5.and.nq.eq.1) then
         gphi = dsqrt(2.D0)*amz**2*Astopstopr(2,1,2)
      elseif(nh.eq.5.and.nq.eq.2) then
         gphi = dsqrt(2.D0)*amz**2*Asbotsbotr(2,1,2)
C --NMSSM
      endif

      kap = dreal(NS_kappa(msq**2,mphi**2,msqp**2,0.D0))

      b0 = (msq**2-mphi**2-msqp**2+kap)/(2.D0*mphi*msqp)

      b1 = (msq**2-mphi**2+msqp**2-kap)/(2.D0*msq*msqp)

      b2 = (mphi**2+msq**2-msqp**2-kap)/(2.D0*mphi*msq)

c -----------------------

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))

      NS_realcorr = 2.D0*gphi/kap*(
     .     (mphi**2-msq**2-msqp**2)*(
     .     -2.D0*dlog((lamv*mphi*msq*msqp)/kap**2)*lb1
     .     +2.D0*lb1**2-lb0**2-lb2**2
     .     +2.D0*dreal(NS_ccspen(dcmplx(1.D0-b1**2)))
     .     -dreal(NS_ccspen(dcmplx(1.D0-b0**2)))
     .     -dreal(NS_ccspen(dcmplx(1.D0-b2**2))) )
     .     +2.D0*kap*dlog((lamv*mphi*msq*msqp)/kap**2)
     .     +4.D0*kap+(2.D0*mphi**2+msq**2+msqp**2)*lb1
     .     +(mphi**2+2.D0*msqp**2)*lb2+(mphi**2+2.D0*msq**2)*lb0 )

      return

      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
c ---------- A.Djouadi, W.Hollik, C.Juenger, hep-ph/9609419 ---------- c
c -------------------------------------------------------------------- c

c --- QCD corrections to the light squark decays --- c

      double precision function NS_ftotqcd(kap,gam)

      implicit double precision (a-h,k-z)

      complex*16 NS_ccspen,funci

      pi = 4.D0*datan(1.D0)

      if(kap*gam.lt.1.D0) then
         funci = NS_ccspen(dcmplx((gam-1.D0)/(gam*kap-1.D0))) -
     .        NS_ccspen(dcmplx(kap*(gam-1.D0)/(gam*kap-1.D0))) -
     .        NS_ccspen(dcmplx((gam+kap-2.D0)/(gam*kap-1.D0))) +
     .        NS_ccspen(dcmplx(kap*(gam+kap-2.D0)/(gam*kap-1.D0)))
      elseif(kap*gam.ge.1.D0) then
         funci = -NS_ccspen(dcmplx((gam*kap-1.D0)/(gam-1.D0))) +
     .        NS_ccspen(dcmplx((gam*kap-1.D0)/(gam+kap-2.D0))) +
     .        NS_ccspen(dcmplx((gam*kap-1.D0)/(kap*(gam-1.D0)))) -
     .        NS_ccspen(dcmplx((gam*kap-1.D0)/(kap*(gam+kap-2.D0)))) -
     .        dlog(kap)*cdlog(dcmplx((gam+kap-2.D0)/(gam-1.D0)))
      endif

      NS_ftotqcd = -1.D0/8.D0*( (4.D0*gam**2-27.D0*gam+25.D0)/
     .(gam-1.D0)
     .     + (3.D0*kap-5.D0)/(kap-1.D0) ) - pi**2/3.D0 
     .     - 2.D0*dreal(NS_ccspen(dcmplx(kap))) 
     .     - 1.D0/2.D0*(gam**2-1.D0)*
     .     dreal(cdlog(dcmplx((gam-1.D0)/gam)))
     .     + (3.D0*gam**2-4.D0*gam+2.D0)/
     .     (4.D0*(1.D0-gam)**2)*dlog(gam) -3.D0/2.D0*dlog(1.D0-kap) +
     .     3.D0/4.D0*(kap**2-4.D0*kap)/(kap-1.D0)**2*dlog(kap)
     .     -dlog(kap)*dlog(1.D0-kap) + dsqrt(kap*gam)*
     .     (1.D0/kap*dlog(1.D0-kap)+1.D0/(1.D0-kap)*(gam*dlog(gam)
     .     -(gam-1.D0)*dreal(cdlog(dcmplx(gam-1.D0))) )
     .     + (kap+gam-2.D0)/(1.D0-kap)**2*dreal(funci) )

      return
    
      end

c -------------------------------------------------------------------- c
C SEE REF PRD55,6975 QCD CORRECTIONS TO SCALAR QURAK DECAYS BY A. DJOUADI
C ET.AL
c --- Heavy squark decays                       --- c
c --- Virtual corrections for the decays        --- c
c --- squark_i -> chargino_j/neutralino_j quark --- c

      DOUBLE PRECISION function NS_gltneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE 
      INTEGER ni,nj,k,idec,I,J
      DOUBLE PRECISION gmst(2),atopr(2,5),btopr(2,5),vt(2),
     . at(2),del(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,Z(5,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1

      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
**        
        Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo
      bet = datan(TANBETA_Z)
      gmst(1) = ast1
      gmst(2) = ast2
      
      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv
      
      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_ftfunctions(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

c --- the running couplings ---



      amqt = amt

      NS_gltneut = 0.D0

      do k=1,2,1
         NS_gltneut = NS_gltneut - 2.D0*(
     .        atopr(k,nj)*( 
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt4ik(ni,nj,k) -
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt5ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt6ik(ni,nj,k) -
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt7ik(ni,nj,k) ) +
     .        btopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt1ik(ni,nj,k) -
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt2ik(ni,nj,k) -
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt3ik(ni,nj,k)) )
      end do

      NS_gltneut = NS_gltneut + btopr(ni,nj)*fnt1(ni,nj) + 
     .     atopr(ni,nj)*fnt2(ni,nj)

      NS_gltneut = NS_gltneut + (-1.D0)**ni*(del(1,ni)*btopr(2,nj)+
     .     del(2,ni)*btopr(1,nj))/(ast1**2-ast2**2)*(
     .     4.D0*amqt*mgluino*dcos(2.D0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     dcos(2.D0*thet)*dsin(2.D0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

      if(ni.eq.1) then
         NS_gltneut = NS_gltneut + 1.D0/2.D0*btopr(1,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dcos(thet)*
     .       NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        runmt/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dsin(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) -
     .        (-dsqrt(2.D0))*sw*(2.D0/3.D0*zp(nj,1)-2.D0/3.D0*sw/cw*
     .        zp(nj,2))*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) 
         
      elseif(ni.eq.2) then
         NS_gltneut = NS_gltneut + 1.D0/2.D0*btopr(2,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*(-dsin(thet))*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        runmt/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) -
     .        (-dsqrt(2.D0))*sw*(2.D0/3.D0*zp(nj,1)-2.D0/3.D0*sw/cw*
     .        zp(nj,2))*(-dsin(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
        
      endif
c      write(*,*)'higi1',NS_gltneut
      NS_gltneut = (-1.D0)*NS_gltneut


      return 

      end

c -------------------------------------------------------------------- c

      double precision function NS_glbneut(ni,nj,amusc,amuscdiv,lamsc)

      IMPLICIT NONE
      integer ni,nj,k,idec,I,J

      DOUBLE PRECISION gmsb(2),abot(2,5),bbot(2,5),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     .NS_A01,NS_delmtdiv
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutsbotbot/abot,bbot
*
      bet=datan(TANBETA_z)
       Do i =1,5
           Do j=1,5
             Z(i,j) = ZZ(i,j)
           enddo
      enddo

C     -------------------------------------
      gmsb(1) = asb1
      gmsb(2) = asb2

      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_fbfunctions(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

c --- the running couplings ---

      amqb = amb

c -------------

      NS_glbneut = 0.D0

      do k=1,2,1
         NS_glbneut = NS_glbneut - 2.D0*(
     .        abot(k,nj)*( 
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb4ik(ni,nj,k) -
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb5ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb6ik(ni,nj,k) -
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb7ik(ni,nj,k) ) +
     .        bbot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb1ik(ni,nj,k) -
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb2ik(ni,nj,k) -
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb3ik(ni,nj,k)) )
      end do

      NS_glbneut = NS_glbneut + bbot(ni,nj)*fnb1(ni,nj) + 
     .     abot(ni,nj)*fnb2(ni,nj)

      NS_glbneut = NS_glbneut + (-1.D0)**ni*(del(1,ni)*bbot(2,nj)+
     .     del(2,ni)*bbot(1,nj))/(asb1**2-asb2**2)*(
     .     4.D0*amqb*mgluino*dcos(2.D0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     dcos(2.D0*theb)*dsin(2.D0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_glbneut = NS_glbneut + 1.D0/2.D0*bbot(1,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dcos(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        runmb/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dsin(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) -
     .        (-dsqrt(2.D0))*sw*(-1.D0/3.D0*zp(nj,1)+1.D0/3.D0*sw/cw*
     .        zp(nj,2))*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) 
      elseif(ni.eq.2) then
         NS_glbneut = NS_glbneut + 1.D0/2.D0*bbot(2,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*(-dsin(theb))*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        runmb/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) -
     .        (-dsqrt(2.D0))*sw*(-1.D0/3.D0*zp(nj,1)+1.D0/3.D0*sw/cw*
     .        zp(nj,2))*(-dsin(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_glbneut = (-1.D0)*NS_glbneut

      return 

      end



c -------------------------------------------------------------------- c
      double precision function NS_grtneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE 
      integer ni,nj,k,idec,I,J
      Double Precision gmst(2),atopr(2,5),btopr(2,5),vt(2),
     .at(2),del(2,2)
      Double Precision fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)

      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
*
      bet=datan(TANBETA_Z)
      Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo
      gmst(1) = ast1
      gmst(2) = ast2

      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_ftfunctions(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

      amqt = amt

      NS_grtneut = 0.D0

      do k=1,2,1
         NS_grtneut = NS_grtneut - 2.D0*(
     .        btopr(k,nj)*( 
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt4ik(ni,nj,k) +
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt5ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt6ik(ni,nj,k) +
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt7ik(ni,nj,k) ) +
     .        atopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt2ik(ni,nj,k) +
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt3ik(ni,nj,k)) )
      end do

      NS_grtneut = NS_grtneut + atopr(ni,nj)*fnt1(ni,nj) + 
     .         btopr(ni,nj)*fnt2(ni,nj)

      NS_grtneut = NS_grtneut + (-1.D0)**ni*(del(1,ni)*atopr(2,nj)+
     .     del(2,ni)*atopr(1,nj))/(ast1**2-ast2**2)*(
     .     4.D0*amqt*mgluino*dcos(2.D0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     dcos(2.D0*thet)*dsin(2.D0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

        if(ni.eq.1) then
         
         NS_grtneut = NS_grtneut + 1.D0/2.D0*atopr(1,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dsin(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt -
     .        runmt/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        dsqrt(2.D0)*sw*(2.D0/3.D0*zp(nj,1)+(1.D0/2.D0-
     .         2.D0/3.D0*sw**2)*1.D0/sw/cw*zp(nj,2))*dsin(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) 

      elseif(ni.eq.2) then
         NS_grtneut = NS_grtneut + 1.D0/2.D0*atopr(2,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*dcos(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt -
     .        runmt/(dsqrt(2.D0)*amw*dsin(bet))*z(nj,4)*(-dsin(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        dsqrt(2.D0)*sw*(2.D0/3.D0*zp(nj,1)+(1.D0/2.D0-
     .         2.D0/3.D0*sw**2)*1.D0/sw/cw*zp(nj,2))*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) 

      endif
c      write(*,*)'higi1',NS_grtneut
      NS_grtneut = (-1.D0)*NS_grtneut
      
      return 

      end

c -------------------------------------------------------------------- c

      double precision function NS_grbneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,idec,I,J
      DOUBLE PRECISION gmsb(2),abot(2,5),bbot(2,5),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqb
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutsbotbot/abot,bbot
*
      bet=datan(TANBETA_Z)
      Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo

      gmsb(1) = asb1
      gmsb(2) = asb2

      idec = ni

      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_fbfunctions(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqb = amb
      NS_grbneut= 0.D0

      do k=1,2,1
        NS_grbneut=NS_grbneut- 2.D0*(
     .        bbot(k,nj)*( 
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb4ik(ni,nj,k) +
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb5ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb6ik(ni,nj,k) +
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb7ik(ni,nj,k) ) +
     .        abot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb2ik(ni,nj,k) +
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb3ik(ni,nj,k)) )
      end do

      NS_grbneut=NS_grbneut+ abot(ni,nj)*fnb1(ni,nj) + 
     .         bbot(ni,nj)*fnb2(ni,nj)

      NS_grbneut=NS_grbneut+ (-1.D0)**ni*(del(1,ni)*abot(2,nj)+
     .     del(2,ni)*abot(1,nj))/(asb1**2-asb2**2)*(
     .     4.D0*amqb*mgluino*dcos(2.D0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     dcos(2.D0*theb)*dsin(2.D0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
        NS_grbneut=NS_grbneut+ 1.D0/2.D0*abot(1,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dsin(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb -
     .        runmb/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        dsqrt(2.D0)*sw*(-1.D0/3.D0*zp(nj,1)+(-1.D0/2.D0+
     .         1.D0/3.D0*sw**2)*1.D0/sw/cw*zp(nj,2))*dsin(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) 
      elseif(ni.eq.2) then
        NS_grbneut=NS_grbneut+ 1.D0/2.D0*abot(2,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) -
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*dcos(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb -
     .        runmb/(dsqrt(2.D0)*amw*dcos(bet))*z(nj,3)*(-dsin(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        dsqrt(2.D0)*sw*(-1.D0/3.D0*zp(nj,1)+(-1.D0/2.D0+
     .         1.D0/3.D0*sw**2)*1.D0/sw/cw*zp(nj,2))*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) 
      endif

      NS_grbneut= (-1.D0)*NS_grbneut

      return 

      end
c -------------------------------------------------------------------- c

      double precision function NS_gltchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      DOUBLE PRECISION  gmst(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1

      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
*
       bet=datan(TANBETA_Z)
       Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
       enddo
*
      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmst(1) = ast1
      gmst(2) = ast2

      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_ftfunctions(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

      amqt = amt
      amqb = amb

      NS_gltchar = 0.D0

      do k=1,2,1
         NS_gltchar = NS_gltchar -2.D0*(
     .        alsbot(k,nj)*( 
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct4ik(ni,nj,k) -
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct5ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct6ik(ni,nj,k) -
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct7ik(ni,nj,k) ) +
     .        aksbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct1ik(ni,nj,k) -
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct2ik(ni,nj,k) -
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct3ik(ni,nj,k)) )
      end do
      
  
     .
      NS_gltchar = NS_gltchar + 
     .         akstor(ni,nj)*fct1(ni,nj) + alstor(ni,nj)*fct2(ni,nj)

      NS_gltchar = NS_gltchar + (-1.D0)**ni*(del(1,ni)*akstor(2,nj)+
     .     del(2,ni)*akstor(1,nj))/(ast1**2-ast2**2)*(
     .     4.D0*amqt*mgluino*dcos(2.D0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     dcos(2.D0*thet)*dsin(2.D0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))
 
         if(ni.eq.1) then
         NS_gltchar = NS_gltchar + 1.D0/2.D0*akstor(1,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1)  +
     .      2.D0*NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv))
     .        - runmb*uu(nj,2)/dsqrt(2.D0)/amw/dcos(bet)*dsin(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      elseif(ni.eq.2) then
         NS_gltchar = NS_gltchar + 1.D0/2.D0*akstor(2,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2)  +
     .      2.D0*NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv))
     .        - runmb*uu(nj,2)/dsqrt(2.D0)/amw/dcos(bet)*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
      endif

      NS_gltchar = (-1.D0)*NS_gltchar
c       write(*,*) "NS_gltchar:",NS_gltchar

      return 

      end

c -------------------------------------------------------------------- c

      double precision function NS_glbchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
      Double Precision gmsb(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      Double Precision alstor(2,2),akstor(2,2),aksbot(2,2),alsbot(2,2)
      Double Precision fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
**
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
*
      bet=datan(TANBETA_Z)
c        
      Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo
    
      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmsb(1) = asb1
      gmsb(2) = asb2

      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_fbfunctions(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqt = amt
      amqb = amb

      NS_glbchar = 0.D0

      do k=1,2,1
         NS_glbchar = NS_glbchar - 2.D0*(
     .        alstor(k,nj)*( 
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb4ik(ni,nj,k) -
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb5ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb6ik(ni,nj,k) -
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb7ik(ni,nj,k) ) +
     .        akstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb1ik(ni,nj,k) -
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb2ik(ni,nj,k) -
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb3ik(ni,nj,k)) )
      end do

      NS_glbchar = NS_glbchar + 
     .         aksbot(ni,nj)*fcb1(ni,nj) + alsbot(ni,nj)*fcb2(ni,nj)

      NS_glbchar = NS_glbchar + (-1.D0)**ni*(del(1,ni)*aksbot(2,nj)+
     .     del(2,ni)*aksbot(1,nj))/(asb1**2-asb2**2)*(
     .     4.D0*amqb*mgluino*dcos(2.D0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     dcos(2.D0*theb)*dsin(2.D0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_glbchar = NS_glbchar + 1.D0/2.D0*aksbot(1,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1)  +
     .      2.D0*NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv))
     .        - runmt*vv(nj,2)/dsqrt(2.D0)/amw/dsin(bet)*dsin(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      elseif(ni.eq.2) then
         NS_glbchar = NS_glbchar + 1.D0/2.D0*aksbot(2,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2)  +
     .      2.D0*NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv))
     .        - runmt*vv(nj,2)/dsqrt(2.D0)/amw/dsin(bet)*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_glbchar = (-1.D0)*NS_glbchar

      return 

      end

c -------------------------------------------------------------------- c
      double precision function NS_grtchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      Double Precision gmst(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      Double Precision alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      Double Precision fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
*     
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
*        
      bet=datan(TANBETA_Z)
        Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo
*
      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmst(1) = ast1
      gmst(2) = ast2

      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_ftfunctions(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)
      amqt = amt
      amqb = amb

      NS_grtchar = 0.D0

      do k=1,2,1
         NS_grtchar = NS_grtchar -2.D0*(
     .        aksbot(k,nj)*( 
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct4ik(ni,nj,k) +
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct5ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct6ik(ni,nj,k) +
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct7ik(ni,nj,k) ) +
     .        alsbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct2ik(ni,nj,k) +
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct3ik(ni,nj,k)) )
      end do

      NS_grtchar = NS_grtchar + 
     .         alstor(ni,nj)*fct1(ni,nj) + akstor(ni,nj)*fct2(ni,nj)

      NS_grtchar = NS_grtchar + (-1.D0)**ni*(del(1,ni)*alstor(2,nj)+
     .     del(2,ni)*alstor(1,nj))/(ast1**2-ast2**2)*(
     .     4.D0*amqt*mgluino*dcos(2.D0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     dcos(2.D0*thet)*dsin(2.D0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

      if(ni.eq.1) then
         NS_grtchar = NS_grtchar + 1.D0/2.D0*alstor(1,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) +
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*vv(nj,2)*dsin(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        vv(nj,1)*dsin(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        runmt*vv(nj,2)/dsqrt(2.D0)/amw/dsin(bet)*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) 
      elseif(ni.eq.2) then
         NS_grtchar = NS_grtchar + 1.D0/2.D0*alstor(2,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) +
     .        1.D0/(dsqrt(2.D0)*amw*dsin(bet))*vv(nj,2)*dcos(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        vv(nj,1)*dcos(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        runmt*vv(nj,2)/dsqrt(2.D0)/amw/dsin(bet)*(-dsin(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
      endif

      NS_grtchar = (-1.D0)*NS_grtchar
      return 

      end

c -------------------------------------------------------------------- c
      double precision function NS_grbchar(ni,nj,amusc,amuscdiv,lamsc)
*       
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      Double Precision gmsb(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      Double Precision alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      Double Precision fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION lamsc,amuvdiv,amusc,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      bet=datan(TANBETA_Z)
        Do i =1,5
           Do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmsb(1) = asb1
      gmsb(2) = asb2

      vt(1) = 1.D0/2.D0*(dcos(thet)-dsin(thet))
      vt(2) = -1.D0/2.D0*(dcos(thet)+dsin(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1.D0/2.D0*(dcos(theb)-dsin(theb))
      vb(2) = -1.D0/2.D0*(dcos(theb)+dsin(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1.D0
      del(1,2) = 0.D0
      del(2,1) = 0.D0
      del(2,2) = 1.D0

      call NS_fbfunctions(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqt = amt
      amqb = amb

      NS_grbchar = 0.D0

      do k=1,2,1
         NS_grbchar = NS_grbchar - 2.D0*(
     .        akstor(k,nj)*( 
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb4ik(ni,nj,k) +
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb5ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb6ik(ni,nj,k) +
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb7ik(ni,nj,k) ) +
     .        alstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb2ik(ni,nj,k) +
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb3ik(ni,nj,k)) )
      end do

      NS_grbchar = NS_grbchar + 
     .         alsbot(ni,nj)*fcb1(ni,nj) + aksbot(ni,nj)*fcb2(ni,nj)

      NS_grbchar = NS_grbchar + (-1.D0)**ni*(del(1,ni)*alsbot(2,nj)+
     .     del(2,ni)*alsbot(1,nj))/(asb1**2-asb2**2)*(
     .     4.D0*amqb*mgluino*dcos(2.D0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     dcos(2.D0*theb)*dsin(2.D0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_grbchar = NS_grbchar + 1.D0/2.D0*alsbot(1,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) +
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*uu(nj,2)*dsin(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        uu(nj,1)*dsin(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        runmb*uu(nj,2)/dsqrt(2.D0)/amw/dcos(bet)*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) 
      elseif(ni.eq.2) then
         NS_grbchar = NS_grbchar + 1.D0/2.D0*alsbot(2,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) +
     .        1.D0/(dsqrt(2.D0)*amw*dcos(bet))*uu(nj,2)*dcos(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        uu(nj,1)*dcos(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        runmb*uu(nj,2)/dsqrt(2.D0)/amw/dcos(bet)*(-dsin(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_grbchar = (-1.D0)*NS_grbchar

      return 

      end

c -------------------------------------------------------------------- c
      subroutine NS_ftfunctions(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,
     .     fnt5ik,fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,
     .     fct5ik,fct6ik,fct7ik)
*
      IMPLICIT NONE
      integer k,idec,I,J
      complex*16 SD_C03,SD_C0_lam
*
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION gmst(2),gmsb(2)
      DOUBLE PRECISION amqt,amqb,amuv,lamv
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_C1_lam,SD_C2_lam,NS_C1,SD_C2
*     
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec   
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      amqt = amt
      amqb = amb

      do i=1,2,1
         do j=1,5,1 
            fnt1(i,j) = SD_B02(gmst(i)**2,lamv,gmst(i),amuv**2) + 
     .        2.D0*amqt**2*dreal(SD_C0_lam(amqt,gmst(i),amneut(j),lamv))
     .         -2.D0*gmst(i)**2*(
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv)-
     .         SD_C2_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv))
     .         +2.D0*amneut(j)**2*
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv) 

            fnt2(i,j) = -2.D0*amqt*xmneut(j)*(
     .         dreal(SD_C0_lam(amqt,gmst(i),amneut(j),lamv)) +
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv) )

            do k=1,2,1
               fnt1ik(i,j,k) = mgluino*xmneut(j)*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .                             gmst(k),mgluino,amqt)) +
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt2ik(i,j,k) = xmneut(j)*(amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt3ik(i,j,k) = xmneut(j)*(-amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt4ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) + amqt*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) )
               fnt5ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) - amqt*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) )
               fnt6ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) + amqt*amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amneut(j)**2*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2) 
               fnt7ik(i,j,k) = gmst(k)**2* 
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) - amqt*amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amneut(j)**2*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2) 
            enddo
         enddo
      enddo

      do i=1,2,1
         do j=1,2,1
            fct1(i,j) = SD_B02(gmst(i)**2,lamv,gmst(i),amuv**2) + 
     .       2.D0*amqb**2*dreal(SD_C0_lam(amqb,gmst(i),amchar(j),lamv)) 
     .         -2.D0*gmst(i)**2*(
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv) -
     .         SD_C2_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv))
     .         +2.D0*amchar(j)**2*
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv)

            fct2(i,j) = -2.D0*amqb*xmchar(j)*(
     .         dreal(SD_C0_lam(amqb,gmst(i),amchar(j),lamv)) +
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv) )

            do k=1,2,1
               fct1ik(i,j,k) = mgluino*xmchar(j)*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct2ik(i,j,k) = xmchar(j)*(amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct3ik(i,j,k) = xmchar(j)*(-amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct4ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) + amqb*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) )
               fct5ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) - amqb*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) )
               fct6ik(i,j,k) = gmsb(k)**2* 
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) + amqt*amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amchar(j)**2*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2)
               fct7ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) - amqt*amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amchar(j)**2*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2) 
            enddo
         enddo
      enddo

      end
C------------------------------------------------------------------------C

      subroutine NS_fbfunctions(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,
     .     fnb5ik,fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,
     .     fcb5ik,fcb6ik,fcb7ik)
*
      IMPLICIT NONE
      integer k,I,J,idec
      complex*16 SD_C03,SD_C0_lam
*
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION gmst(2),gmsb(2)
      DOUBLE PRECISION amuv,lamv
      DOUBLE PRECISION amqt,amqb
      DOUBLE PRECISION SD_B02,NS_C1_lam,SD_C2_lam,NS_C1,SD_C2
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_decindex/idec
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      amqt = amt
      amqb = amb
 
      do i=1,2,1
         do j=1,5,1 
            fnb1(i,j) = SD_B02(gmsb(i)**2,lamv,gmsb(i),amuv**2) + 
     .       2.D0*amqb**2*dreal(SD_C0_lam(amqb,gmsb(i),amneut(j),lamv))
     .         -2.D0*gmsb(i)**2*(
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv)-
     .         SD_C2_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv))
     .         +2.D0*amneut(j)**2*
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv) 

            fnb2(i,j) = -2.D0*amqb*xmneut(j)*(
     .         dreal(SD_C0_lam(amqb,gmsb(i),amneut(j),lamv)) +
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv) )

            do k=1,2,1
               fnb1ik(i,j,k) = mgluino*xmneut(j)*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .                             gmsb(k),mgluino,amqb)) +
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb2ik(i,j,k) = xmneut(j)*(amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb3ik(i,j,k) = xmneut(j)*(-amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb4ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) + amqb*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) )
               fnb5ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) - amqb*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) )
               fnb6ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) + amqb*amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amneut(j)**2*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
               fnb7ik(i,j,k) = gmsb(k)**2* 
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) - amqb*amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amneut(j)**2*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2) 
            enddo
         enddo
      enddo

      do i=1,2,1
         do j=1,2,1
            fcb1(i,j) = SD_B02(gmsb(i)**2,lamv,gmsb(i),amuv**2) + 
     .       2.D0*amqt**2*dreal(SD_C0_lam(amqt,gmsb(i),amchar(j),lamv)) 
     .         -2.D0*gmsb(i)**2*(
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv) -
     .         SD_C2_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv))
     .         +2.D0*amchar(j)**2*
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv)

            fcb2(i,j) = -2.D0*amqt*xmchar(j)*(
     .         dreal(SD_C0_lam(amqt,gmsb(i),amchar(j),lamv)) +
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv) )

            do k=1,2,1
               fcb1ik(i,j,k) = mgluino*xmchar(j)*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb2ik(i,j,k) = xmchar(j)*(amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb3ik(i,j,k) = xmchar(j)*(-amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb4ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) + amqt*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) )
               fcb5ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) - amqt*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) )
               fcb6ik(i,j,k) = gmst(k)**2* 
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) + amqb*amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amchar(j)**2*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
               fcb7ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) - amqb*amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amchar(j)**2*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2) 
            enddo
         enddo
      enddo

      end

c -------------------------------------------------------------------- c
      double precision function NS_delmtdiv(amq,mgluino,ast1,ast2,thet,
     .                                   amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_sigmardiv,NS_sigmaldiv,NS_sigmasdiv

      NS_delmtdiv = 1.D0/2.D0*(
     .     NS_sigmardiv(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmaldiv(amq,mgluino,ast1,ast2,thet,amuv,lamv) ) +
     .     NS_sigmasdiv(amq,mgluino,ast1,ast2,thet,amuv,lamv)

      return 
      
      end 
c -------------------------------------------------------------------- c
      double precision function NS_delztr(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_sigmar,NS_sigmalp,NS_sigmarp,NS_sigmasp

      NS_delztr = -NS_sigmar(amq,mgluino,ast1,ast2,thet,amuv,lamv) -
     .     amq**2*(NS_sigmalp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmarp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     2.D0*NS_sigmasp(amq,mgluino,ast1,ast2,thet,amuv,lamv) )

      return

      end

c -------------------------------------------------------------------- c
      double precision function NS_delztl(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_sigmal,NS_sigmalp,NS_sigmarp,NS_sigmasp

      NS_delztl = -NS_sigmal(amq,mgluino,ast1,ast2,thet,amuv,lamv) -
     .     amq**2*(NS_sigmalp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmarp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     2.D0*NS_sigmasp(amq,mgluino,ast1,ast2,thet,amuv,lamv) )

      return

      end

c -------------------------------------------------------------------- c
      double precision function NS_delzst(amqt,mgluino,amsq,thet,
     .                                 amuv,lamv,ni)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1,SD_BP1,SD_BP02,SD_B02
      integer ni

      NS_delzst = 2.D0*(-2.D0*NS_B1(amsq**2,amsq,lamv,amuv**2)
     .     -2.D0*amsq**2*SD_BP1(amsq**2,amsq,lamv,amuv**2) +
     .     (amqt**2+mgluino**2-amsq**2)*
     .     SD_BP02(amsq**2,amqt,mgluino,amuv**2) -
     .     SD_B02(amsq**2,amqt,mgluino,amuv**2) + 
     .     (-1.D0)**ni*2.D0*dsin(2.D0*thet)*amqt*mgluino*
     .     SD_BP02(amsq**2,amqt,mgluino,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c
      double precision function NS_sigmar(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1

      NS_sigmar = -(2.D0*NS_B1(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0-dcos(2.D0*thet))*NS_B1(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0+dcos(2.D0*thet))*NS_B1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmardiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV

      NS_sigmardiv = -(2.D0*NS_B1_DIV(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0-dcos(2.D0*thet))*NS_B1_DIV(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0+dcos(2.D0*thet))*NS_B1_DIV(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmarp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_BP1

      NS_sigmarp = -(2.D0*SD_BP1(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0-dcos(2.D0*thet))*SD_BP1(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0+dcos(2.D0*thet))*SD_BP1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmal(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1

      NS_sigmal = -(2.D0*NS_B1(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0+dcos(2.D0*thet))*NS_B1(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0-dcos(2.D0*thet))*NS_B1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmaldiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV

      NS_sigmaldiv = -(2.D0*NS_B1_DIV(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0+dcos(2.D0*thet))*NS_B1_DIV(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0-dcos(2.D0*thet))*NS_B1_DIV(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmalp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_BP1

      NS_sigmalp = -(2.D0*SD_BP1(amq**2,amq,lamv,amuv**2)+ 
     .     (1.D0+dcos(2.D0*thet))*SD_BP1(amq**2,mgluino,ast1,amuv**2)
     .    +(1.D0-dcos(2.D0*thet))*SD_BP1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end 

c -------------------------------------------------------------------- c
      double precision function NS_sigmasdiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV

      NS_sigmasdiv = -(4.D0*SD_B02_DIV(amq**2,amq,lamv,amuv**2)+ 
     .     mgluino/amq*dsin(2.D0*thet)*(
     .     SD_B02_DIV(amq**2,mgluino,ast1,amuv**2)
     .    -SD_B02_DIV(amq**2,mgluino,ast2,amuv**2)
     .     ) )

      return

      end 
c -------------------------------------------------------------------- c
      double precision function NS_sigmasp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_BP02

      NS_sigmasp = -(4.D0*SD_BP02(amq**2,amq,lamv,amuv**2)+ 
     .     mgluino/amq*dsin(2.D0*thet)*(
     .     SD_BP02(amq**2,mgluino,ast1,amuv**2)
     .    -SD_BP02(amq**2,mgluino,ast2,amuv**2)
     .     ) )

      return

      end 

c -------------------------------------------------------------------- c
      double precision function NS_delthdiv(amqt,mgluino,ast1,ast2,
     .                                   thet,amuv)
      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV

      NS_delthdiv = 1.D0/(ast1**2-ast2**2)*(4.D0*amqt*mgluino*
     .     dcos(2.D0*thet)*SD_B02_DIV(ast2**2,amqt,mgluino,amuv**2) +
     .     dcos(2.D0*thet)*dsin(2.D0*thet)*
     .     (ast2**2-ast1**2)*dlog(amuv**2) )

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c

      double precision function NS_corrreali(amq,mcharneut,amsti,lamv,
     .     icharneut,isign,ni,nj,idec)
*
      implicit double precision (a-h,k-z)
      dimension atopr(2,5),btopr(2,5),alstor(2,2),akstor(2,2),
     .     abot(2,5),bbot(2,5),alsbot(2,2),aksbot(2,2)
      integer ni,nj
      complex*16 NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2 
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
      External NS_ccspen,NS_kappa
  
      kap = dreal(NS_kappa(amsti**2,amq**2,mcharneut**2,0.D0))

      b0 = (amsti**2-amq**2-mcharneut**2+kap)/(2.D0*amq*mcharneut)
      b1 = (amsti**2-amq**2+mcharneut**2-kap)/(2.D0*amsti*mcharneut)
      b2 = (amsti**2+amq**2-mcharneut**2-kap)/(2.D0*amsti*amq)

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hi00 = 1.D0/4.D0/amsti**4*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (amq**2-mcharneut**2)*lb12 - amsti**2*lb0 )
      hi11 = 1.D0/4.D0/(amq**2*amsti**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (amsti**2-mcharneut**2)*lb02 - amq**2*lb1)
      hi01 = dreal(1.D0/(4.D0*amsti**2)*(-2.D0*
     .     dlog((lamv*amsti*amq*mcharneut)/kap**2)*lb2 + 
     .     2.D0*lb2**2 - lb0**2 - lb1**2 +
     .     2.D0*NS_ccspen(dcmplx(1.D0-b2**2)) - 
     .     NS_ccspen(dcmplx(1-b0**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )
      hi = 1.D0/(4.D0*amsti**2)*(kap/2.D0*
     .     (amsti**2+amq**2+mcharneut**2) + 2.D0*amsti**2*amq**2*
     .     lb2 + 2.D0*amsti**2*mcharneut**2*lb1 + 
     .     2.D0*amq**2*mcharneut**2*lb0)
      hi0 = 1.D0/(4.D0*amsti**2)*(-2.D0*amq**2*lb2
     .     -2.D0*mcharneut**2*lb1-kap)
      hi1 = 1.D0/(4.D0*amsti**2)*(-2.D0*amsti**2*lb2 
     .     -2.D0*mcharneut**2*lb0-kap)
      hi10 = 1.D0/(4.D0*amsti**2)*(amsti**4*lb2-
     .     mcharneut**2*(2.D0*amq**2-2.D0*amsti**2+mcharneut**2)*
     .     lb0 - kap/4.D0*(amq**2-3.D0*amsti**2+
     .     5.D0*mcharneut**2) )

      if(icharneut.eq.1.and.idec.eq.1) then
         cli = - btopr(ni,nj)
         cri = - atopr(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.1) then
         cli = - akstor(ni,nj)
         cri = - alstor(ni,nj)
      elseif(icharneut.eq.1.and.idec.eq.2) then
         cli = - bbot(ni,nj)
         cri = - abot(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.2) then
         cli = - aksbot(ni,nj)
         cri = - alsbot(ni,nj)
      endif

      if(isign.eq.1) then
         epschi = -1.D0
      elseif(isign.eq.0) then
         epschi = 1.D0
      endif
         
      NS_corrreali = 8.D0*cli*cri*amq*mcharneut*epschi*((amsti**2+amq**2
     .     -mcharneut**2)*hi01+amsti**2*hi00+amq**2*hi11+hi0+hi1)
     .     +(cli**2+cri**2)*(2.D0*(amq**2+mcharneut**2-amsti**2)*
     .     (amsti**2*hi00+amq**2*hi11+hi0+hi1) + 2.D0*(amq**4-
     .     (mcharneut**2-amsti**2)**2)*hi01-hi-hi10)

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c
c ----------- for the processes gaugino -> squark + quark ------------ c

      double precision function NS_realicorr(amq,mcharneut,amsti,lamv,
     .     icharneut,isign,ni,nj,idec)

      implicit double precision (a-h,k-z)
      integer ni,nj
*
      DOUBLE PRECISION atopr(2,5),btopr(2,5),alstor(2,2),akstor(2,2),
     .     abot(2,5),bbot(2,5),alsbot(2,2),aksbot(2,2)
     
      complex*16 NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
*
      external NS_kappa

      kap = dreal(NS_kappa(amsti**2,amq**2,mcharneut**2,0.D0))

      b0 = (mcharneut**2-amq**2-amsti**2+kap)/(2.D0*amq*amsti)
      b1 = (mcharneut**2-amq**2+amsti**2-kap)/(2.D0*amsti*mcharneut)
      b2 = (mcharneut**2+amq**2-amsti**2-kap)/(2.D0*mcharneut*amq)

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb01 = dreal(cdlog(dcmplx(b0/b1)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hi11 = 1.D0/4.D0/(amq**2*mcharneut**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (mcharneut**2-amsti**2)*lb02 - amq**2*lb1)
      hi22 = 1.D0/4.D0/(mcharneut**2*amsti**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (mcharneut**2-amq**2)*lb01 - amsti**2*lb2)
      hi21 = dreal(1.D0/(4.D0*mcharneut**2)*(-2.D0*
     .     dlog((lamv*amsti*amq*mcharneut)/kap**2)*lb0 + 
     .     2.D0*lb0**2 - lb2**2 - lb1**2 +
     .     2.D0*NS_ccspen(dcmplx(1.D0-b0**2)) - 
     .     NS_ccspen(dcmplx(1-b2**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )
      hi = 1.D0/(4.D0*mcharneut**2)*(kap/2.D0*
     .     (mcharneut**2+amq**2+amsti**2) + 2.D0*mcharneut**2*amq**2*
     .     lb2 + 2.D0*amsti**2*mcharneut**2*lb1 + 
     .     2.D0*amq**2*amsti**2*lb0)
      hi2 = 1.D0/(4.D0*mcharneut**2)*(-2.D0*mcharneut**2*lb1
     .     -2.D0*amq**2*lb0-kap)
      hi1 = 1.D0/(4.D0*mcharneut**2)*(-2.D0*mcharneut**2*lb2 
     .     -2.D0*amsti**2*lb0-kap)
      hi12 = 1.D0/(4.D0*mcharneut**2)*(amsti**4*lb0-
     .     mcharneut**2*(2.D0*amq**2-2.D0*amsti**2+mcharneut**2)*
     .     lb2 - kap/4.D0*(amq**2-3.D0*amsti**2+
     .     5.D0*mcharneut**2) )

      if(icharneut.eq.1.and.idec.eq.1) then
         cli = - btopr(ni,nj)
         cri = - atopr(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.1) then
         cli = - akstor(ni,nj)
         cri = - alstor(ni,nj)
      elseif(icharneut.eq.1.and.idec.eq.2) then
         cli = - bbot(ni,nj)
         cri = - abot(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.2) then
         cli = - aksbot(ni,nj)
         cri = - alsbot(ni,nj)
      endif

      if(isign.eq.1) then
         epschi = -1.D0
      elseif(isign.eq.0) then
         epschi = 1.D0
      endif
         
      NS_realicorr = 8.D0*cli*cri*amq*mcharneut*epschi*((amsti**2+amq**2
     .     -mcharneut**2)*hi21+amsti**2*hi22+amq**2*hi11+hi2+hi1)
     .     +(cli**2+cri**2)*(2.D0*(amq**2+mcharneut**2-amsti**2)*
     .     (amsti**2*hi22+amq**2*hi11+hi2+hi1) + 2.D0*(amq**4-
     .     (mcharneut**2-amsti**2)**2)*hi21-hi-hi12)

      return

      end

c -------------------------------------------------------------------- c
c ---------- Beenakker, Hoepker and Zerwas, hep-ph/9602378 ----------- c
c -------------------------------------------------------------------- c

c - QCD corrections to the decays light squark -> light quark + gluino c
      double precision function NS_gama(r)
* 
      implicit double precision (a-h,k-z)
      complex*16 NS_ccspen
*
      pi = 4.D0*datan(1.D0)

      NS_gama = 3.D0/(r-1.D0)*dreal(NS_ccspen(dcmplx(1.D0-r))) - 
     .     r/(r-1.D0)*dreal(NS_ccspen(dcmplx(-r))) 
     .     + (5.D0*r-6.D0)/(12.D0*(r-1.D0))*pi**2 +
     .     59.D0/24.D0 + r/(4.D0*(r-1.D0)) + 
     .     ( (3.D0+r)/(2.D0*(r-1.D0))*dlog(r) - 2.D0)*
     .     dlog(dabs(1.D0-r)) + 
     .     ( (r*(5.D0*r-6.D0))/(4.D0*(r-1.D0)**2) - 
     .     r/(r-1.D0)*dlog(1.D0+r) )*dlog(r) 

      return

      end

c -------------------------------------------------------------------- c
      double precision function NS_gamfcap(r)
*
      implicit double precision (a-h,k-z)
      complex*16 NS_ccspen
*
      pi = 4.D0*datan(1.D0)

      NS_gamfcap = -2.D0/(r-1.D0)*dreal(NS_ccspen(dcmplx(1.D0-r))) + 
     .     (2.D0*r)/(r-1.D0)*dreal(NS_ccspen(dcmplx(-r))) + 
     .     (4.D0-3.D0*r)/(6.D0*(r-1.D0))*pi**2 + 5.D0/2.D0 - r/2.D0 + 
     .     ( r- r**2/2.D0 - (r+1.D0)/(r-1.D0)*dlog(r) )*
     .     dlog(dabs(1.D0-r)) + 
     .     ( 2.D0*r/(r-1.D0)*dlog(1.D0+r)-r+r**2/2.D0)*dlog(r) 

      return

      end

c -------------------------------------------------------------------- c
      double precision function NS_gamf(r)

      implicit double precision (a-h,k-z)

      NS_gamf = -3.D0/(4.D0*r) + 
     .     ((r-1.D0)*(r+3.D0))/(4.D0*r**2)*dlog(dabs(1.D0-r))

      return

      end

c -------------------------------------------------------------------- c
c maggie changed with respect to the paper 26/3/03
      double precision function NS_gamrendec(amsq,amst1,amst2,amt,amsb1,
     .     amsb2,amgl,scala)

      implicit double precision (a-h,k-z)

      mur = scala

      NS_gamrendec = 3.D0/4.D0*dlog(mur**2/amsq**2) - 1.D0/4.D0
    
      NS_gamrendec =  NS_gamrendec + 4.D0/12.D0*dlog(mur**2/amsq**2)
     .     + 1.D0/24.D0*dlog(mur**2/amsb1**2) + 1.D0/24.D0*
     .     dlog(mur**2/amsb2**2) 
     .     + 1.D0/24.D0*dlog(mur**2/amst1**2) + 1.D0/24.D0*
     .     dlog(mur**2/amst2**2) + 1.D0/6.D0*dlog(mur**2/amt**2)
     .     + 1.D0/2.D0*dlog(mur**2/amgl**2) 
      return

      end

c end maggie changed

c -------------------------------------------------------------------- c
c ------ Beenakker, Hoepker, Plehn and Zerwas, hep-ph/9610313 -------- c
c -------------------------------------------------------------------- c

c -------- QCD corrections to the decay stop1/2 -> top gluino -------- c
c -------- and sbottom1/2 -> bottom gluino, gluino -> stop1/2 top ---- c
c -------- gluino -> sbottom1/2 bottom ------------------------------- c

      double precision function NS_gamtop1(amst1,amst2,amt,amgl,thet,
     .                                     ival,amuv,lamv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02
      double precision isign

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif 

      sig2t = isign*amt*amgl*dsin(2.D0*thet)

      NS_gamtop1 = 2.D0*NS_A01(amt**2,amuv**2) -2.D0*amt**2 + 
     .   2.D0*NS_A01(amgl**2,amuv**2) - NS_A01(amst1**2,amuv**2)
     .   - NS_A01(amst2**2,amuv**2) +
     .   (amt**2+amst1**2-amgl**2)*SD_B02(amt**2,amgl,amst1,amuv**2)
     .  +(amt**2+amst2**2-amgl**2)*SD_B02(amt**2,amgl,amst2,amuv**2)
     .   -4.D0*amt**2*sig2t*(SD_BP02(amt**2,amgl,amst1,amuv**2) 
     .   - SD_BP02(amt**2,amgl,amst2,amuv**2) ) + 2.D0*amt**2*(
     .   (amgl**2+amt**2-amst1**2)*SD_BP02(amt**2,amgl,amst1,amuv**2)
     .  +(amgl**2+amt**2-amst2**2)*SD_BP02(amt**2,amgl,amst2,amuv**2)
     .  -4.D0*amt**2*SD_BP02(amt**2,lamv,amt,amuv**2) )

      NS_gamtop1 = -1.D0/16.D0/pi**2*NS_gamtop1

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gamtop2(amst1,amst2,amt,amgl,thet,
     .     ival,amuv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02
      double precision isign

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif 

      NS_gamtop2 = NS_A01(amst2**2,amuv**2) - 
     .  NS_A01(amst1**2,amuv**2)
     .  -(amgl**2+amt**2-amst1**2)*SD_B02(amt**2,amgl,amst1,amuv**2)
     .  +(amgl**2+amt**2-amst2**2)*SD_B02(amt**2,amgl,amst2,amuv**2)

      NS_gamtop2 = (isign*dcos(2.D0*thet))**2*(amgl**2+amt**2-amst1**2)*
     .     NS_gamtop2

      NS_gamtop2 = -1.D0/16.D0/pi**2*NS_gamtop2

      return

      end
c -------------------------------------------------------------------- c

      double precision function NS_gamglui1(amst1,amst2,amsq,amt,amgl,
     .     amuv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_BP02

      pi = 4.D0*datan(1.D0)

      NS_gamglui1 = -NS_A01(amsq**2,amuv**2)+(amsq**2+amgl**2)*
     .     SD_B02(amgl**2,amsq,0.D0,amuv**2) + 2.D0*amgl**2*
     .     (amgl**2-amsq**2)*SD_BP02(amgl**2,amsq,0.D0,amuv**2)

      NS_gamglui1 = -1.D0/16.D0/pi**2*NS_gamglui1

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gamglui2(amst1,amst2,amt,thet,amsb1,
     .     amsb2,amb,theb,amgl,ival,amuv)

      implicit double precision (a-h,k-z)
      double precision isign
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif 

      sig2t = isign*amt*amgl*dsin(2.D0*thet)
      sig2b = isign*amb*amgl*dsin(2.D0*theb)

      NS_gamglui2 = 2.D0*NS_A01(amt**2,amuv**2) -
     .     NS_A01(amst1**2,amuv**2) - NS_A01(amst2**2,amuv**2)
     .     + (amst1**2+amgl**2-amt**2)*
     .     SD_B02(amgl**2,amst1,amt,amuv**2) + 
     .     (amst2**2+amgl**2-amt**2)*
     .     SD_B02(amgl**2,amst2,amt,amuv**2) - 
     .     4.D0*amgl**2*sig2t*( SD_BP02(amgl**2,amst1,amt,amuv**2)
     .     - SD_BP02(amgl**2,amst2,amt,amuv**2) )
     .     + 2.D0*amgl**2*( (amgl**2+amt**2-amst1**2)*
     .     SD_BP02(amgl**2,amst1,amt,amuv**2) + 
     .     (amgl**2+amt**2-amst2**2)*
     .     SD_BP02(amgl**2,amst2,amt,amuv**2) )

c maggie changed with respect to the paper 26/3/03
      NS_gamglui2 = NS_gamglui2 +
     .     2.D0*NS_A01(amb**2,amuv**2) -
     .     NS_A01(amsb1**2,amuv**2) - NS_A01(amsb2**2,amuv**2)
     .     + (amsb1**2+amgl**2-amb**2)*
     .     SD_B02(amgl**2,amsb1,amb,amuv**2) + 
     .     (amsb2**2+amgl**2-amb**2)*
     .     SD_B02(amgl**2,amsb2,amb,amuv**2) - 
     .     4.D0*amgl**2*sig2b*( SD_BP02(amgl**2,amsb1,amb,amuv**2)
     .     - SD_BP02(amgl**2,amsb2,amb,amuv**2) )
     .     + 2.D0*amgl**2*( (amgl**2+amb**2-amsb1**2)*
     .     SD_BP02(amgl**2,amsb1,amb,amuv**2) + 
     .     (amgl**2+amb**2-amsb2**2)*
     .     SD_BP02(amgl**2,amsb2,amb,amuv**2) )
c end maggie changed

      NS_gamglui2 = -1.D0/16.D0/pi**2*NS_gamglui2
    
      return

      end


c -------------------------------------------------------------------- c

      double precision function NS_gamglui3(amgl,amuv,lamv)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_BP02

      pi = 4.D0*datan(1.D0)

      NS_gamglui3 =  NS_A01(amgl**2,amuv**2) - amgl**2 - 
     .     4.D0*amgl**4*SD_BP02(amgl**2,lamv,amgl,amuv**2)

      NS_gamglui3 = -1.D0/16.D0/pi**2*NS_gamglui3

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gam11(amst1,amst2,amt,amgl,thet,ival,
     .     amuv,lamv)

      implicit double precision (a-h,k-z)
      double precision isign
      DOUBLE PRECISION SD_B02,SD_BP02

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif 

      sig2t = isign*amt*amgl*dsin(2.D0*thet)

      NS_gam11 = SD_B02(amst1**2,amgl,amt,amuv**2) - 
     .     SD_B02(amst1**2,lamv,amst1,amuv**2) + 
     .     2.D0*sig2t*SD_BP02(amst1**2,amgl,amt,amuv**2) -
     .     (amgl**2+amt**2-amst1**2)*
     .     SD_BP02(amst1**2,amgl,amt,amuv**2) - 2.D0*amst1**2*
     .     SD_BP02(amst1**2,lamv,amst1,amuv**2)

      NS_gam11 = -1.D0/16.D0/pi**2*NS_gam11

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gam12(amst1,amst2,amt,amgl,thet,ival,
     .     amuv,lamv,scala)

      implicit double precision (a-h,k-z)
      double precision isign
      DOUBLE PRECISION NS_A01,SD_B02,NS_A01_DIV,SD_B02_DIV

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif 

      sig2t = isign*amt*amgl*dsin(2.D0*thet)

      NS_gam12 = NS_A01(amst2**2,amuv**2) - NS_A01(amst1**2,amuv**2)
     .     + 4.D0*amt**2*amgl**2/sig2t*
     .     SD_B02(amst1**2,amgl,amt,amuv**2)

      NS_gam12 = NS_gam12 - 
     .     ( NS_A01_DIV(amst2**2,(amuv/scala)**2) 
     .     - NS_A01_DIV(amst1**2,(amuv/scala)**2)
     .     + 4.D0*amt**2*amgl**2/sig2t*
     .       SD_B02_DIV(amst1**2,amgl,amt,(amuv/scala)**2) )

      NS_gam12 = 1.D0/(amst1**2-amst2**2)*sig2t*dcos(2.D0*thet)**2*
     .     NS_gam12

      NS_gam12 = -1.D0/16.D0/pi**2*NS_gam12

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gamvirt(amst1,amst2,amtop,amgl,thet,
     .     ival,amuv,lamv)

      implicit double precision (a-h,k-z)
      double precision isign

      dimension fffunc(3),fafunc(3)

      COMMON/NS_qcdscales/amuvv,lamvv
      COMMON/NS_relmasses/mst1,mst2,mgl,mtop

      amuvv = amuv
      lamvv = lamv
      mst1  = amst1
      mst2  = amst2
      mgl   = amgl
      mtop  = amtop

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif

      sig2t = isign*amtop*amgl*dsin(2.D0*thet)

      call NS_fifafunctions(fffunc,fafunc,fbfunc)

      NS_gamvirt = 64.D0/9.D0*pi*
     .     ( fffunc(1) + sig2t*fffunc(2) + sig2t**2*fffunc(3) ) +
     .     8.D0*pi*
     .     (fafunc(1) + sig2t*fafunc(2) + sig2t**2*fafunc(3)) +
     .     16.D0/3.D0*pi*
     .     (-(amgl**2+mtop**2-amst1**2)+2.D0*sig2t)*fbfunc

      return

      end

c -------------------------------------------------------------------- c
c --- this function is for the stop1/2, sbottom1/2 decays --- 

      double precision function NS_gamreal(amst1,amt,amgl,thet,ival,
     .                                     lamv)

      implicit double precision (a-h,k-z)
      double precision ist1gl,ist1t,iglgl,ist1st1,itt,igl,ist1,it,
     .     itgl,iglst1
      double precision isign

      complex*16 NS_kappa,NS_ccspen

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif

      sig2t = isign*amt*amgl*dsin(2.D0*thet)

      m0 = amst1
      m1 = amt
      m2 = amgl

      kap = dreal(NS_kappa(m0**2,m1**2,m2**2,0.D0))

      b0 = (m0**2-m1**2-m2**2+kap)/(2.D0*m1*m2)
      b1 = (m0**2-m1**2+m2**2-kap)/(2.D0*m0*m2)
      b2 = (m0**2+m1**2-m2**2-kap)/(2.D0*m0*m1)

      ist1gl = dreal( 
     .     1.D0/(4.D0*m0**2)*(-2.D0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b1) + 2.D0*(dlog(b1))**2 - (dlog(b0))**2 - 
     .     (dlog(b2))**2 + 2.D0*NS_ccspen(dcmplx(1.D0-b1**2)) 
     .     - NS_ccspen(dcmplx(1.D0-b0**2)) - 
     .     NS_ccspen(dcmplx(1.D0-b2**2)) ) )

      ist1t = dreal( 
     .     1.D0/(4.D0*m0**2)*(-2.D0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b2) + 2.D0*(dlog(b2))**2 - (dlog(b0))**2 - 
     .     (dlog(b1))**2 + 2.D0*NS_ccspen(dcmplx(1.D0-b2**2)) 
     .     - NS_ccspen(dcmplx(1.D0-b0**2)) - 
     .     NS_ccspen(dcmplx(1.D0-b1**2)) ) )

      iglgl = 1.D0/(4.D0*m2**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m1**2)*dlog(b0/b1)-m2**2*dlog(b2) )

      ist1st1 = 1.D0/(4.D0*m0**4)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m1**2-m2**2)*dlog(b1/b2)-m0**2*dlog(b0))

      itt = 1.D0/(4.D0*m1**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m2**2)*dlog(b0/b2)-m1**2*dlog(b1) )

      igl = 1.D0/(4.D0*m0**2)*(-2.D0*m0**2*dlog(b1)-2.D0*m1**2*dlog(b0)
     .     -kap)

      ist1 = 1.D0/(4.D0*m0**2)*(-2.D0*m1**2*dlog(b2)-
     .     2.D0*m2**2*dlog(b1)-kap)

      it = 1.D0/(4.D0*m0**2)*(-2.D0*m0**2*dlog(b2)-2.D0*m2**2*dlog(b0)
     .     -kap)

      itgl = -1.D0/(4.D0*m0**2)*(-m2**4*dlog(b0)+m0**2*(2.D0*m1**2
     .     -2.D0*m2**2)*dlog(b2) + m0**4*dlog(b2) + 
     .     kap/4.D0*(m1**2+5.D0*m0**2-3.D0*m2**2) )
      
      iglst1 = 1.D0/(4.D0*m0**2)*(m0**4*dlog(b1)-m1**2*(2.D0*m2**2
     .     -2.D0*m0**2+m1**2)*dlog(b0) -kap/4.D0*(m2**2-3.D0*m0**2
     .     +5.D0*m1**2) )

      NS_gamreal = 8.D0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2.D0*sig2t)*(-(amst1**2-amt**2)*ist1gl+amt**2*ist1t
     .     -amgl**2*iglgl-igl) +
     .     32.D0/9.D0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2.D0*sig2t)*(-amst1**2*ist1st1-amt**2*itt-
     .     (amt**2+amst1**2-amgl**2)*ist1t-ist1-it) +
     .     4.D0/3.D0/pi/amst1*(4.D0/3.D0*itgl-3.D0*iglst1)

      return

      end

c -------------------------------------------------------------------- c
c --- this function is for the gluino decays --- 

      double precision function NS_gamrealgl(amst1,amt,amgl,thet,ival,
     .                                       lamv)

      implicit double precision (a-h,k-z)
      double precision ist1gl,ist1t,iglgl,ist1st1,itt,igl,ist1,it,
     .     itgl,iglst1
      double precision isign

      complex*16 NS_kappa,NS_ccspen

      pi = 4.D0*datan(1.D0)

      if(ival.eq.1) then
         isign = 1.D0
      elseif(ival.eq.2) then
         isign = -1.D0
      endif

      sig2t = isign*amt*amgl*dsin(2.D0*thet)

      m0 = amgl
      m1 = amt
      m2 = amst1

      kap = dreal(NS_kappa(m0**2,m1**2,m2**2,0.D0))

      b0 = (m0**2-m1**2-m2**2+kap)/(2.D0*m1*m2)
      b1 = (m0**2-m1**2+m2**2-kap)/(2.D0*m0*m2)
      b2 = (m0**2+m1**2-m2**2-kap)/(2.D0*m0*m1)

      ist1gl = dreal( 
     .     1.D0/(4.D0*m0**2)*(-2.D0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b1) + 2.D0*(dlog(b1))**2 - (dlog(b0))**2 - 
     .     (dlog(b2))**2 + 2.D0*NS_ccspen(dcmplx(1.D0-b1**2)) 
     .     - NS_ccspen(dcmplx(1.D0-b0**2)) - 
     .     NS_ccspen(dcmplx(1.D0-b2**2)) ) )

      ist1t = dreal( 
     .     1.D0/(4.D0*m0**2)*(-2.D0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b0) + 2.D0*(dlog(b0))**2 - (dlog(b1))**2 - 
     .     (dlog(b2))**2 + 2.D0*NS_ccspen(dcmplx(1.D0-b0**2)) 
     .     - NS_ccspen(dcmplx(1.D0-b1**2)) - 
     .     NS_ccspen(dcmplx(1.D0-b2**2)) ) )

      iglgl = 1.D0/(4.D0*m0**4)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m1**2-m2**2)*dlog(b1/b2)-m0**2*dlog(b0) )

      ist1st1 = 1.D0/(4.D0*m2**2*m0**2)*(
     .     kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m1**2)*dlog(b0/b1)-m2**2*dlog(b2))

      itt = 1.D0/(4.D0*m1**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m2**2)*dlog(b0/b2)-m1**2*dlog(b1) )

      igl = 1.D0/(4.D0*m0**2)*(-2.D0*m1**2*dlog(b2)-2.D0*m2**2*dlog(b1)
     .     -kap)

      ist1 = 1.D0/(4.D0*m0**2)*(-2.D0*m0**2*dlog(b1)
     .      -2.D0*m1**2*dlog(b0)-kap)

      it = 1.D0/(4.D0*m0**2)*(-2.D0*m0**2*dlog(b2)-2.D0*m2**2*dlog(b0)
     .     -kap)

      itgl = -1.D0/(4.D0*m0**2)*(-m0**4*dlog(b2)+m2**2*(2.D0*m1**2
     .     -2.D0*m0**2)*dlog(b0) + m2**4*dlog(b0) + 
     .     kap/4.D0*(m1**2+5.D0*m2**2-3.D0*m0**2) )

      iglst1 = 1.D0/(4.D0*m0**2)*(m2**4*dlog(b1)-m1**2*(2.D0*m0**2
     .     -2.D0*m2**2+m1**2)*dlog(b2) -kap/4.D0*(m0**2-3.D0*m2**2
     .     +5.D0*m1**2) )

      NS_gamrealgl = 8.D0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2.D0*sig2t)*(-(amst1**2-amt**2)*ist1gl+amt**2*ist1t
     .     -amgl**2*iglgl-igl) +
     .     32.D0/9.D0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2.D0*sig2t)*(-amst1**2*ist1st1-amt**2*itt-
     .     (amt**2+amst1**2-amgl**2)*ist1t-ist1-it) +
     .     4.D0/3.D0/pi/amst1*(4.D0/3.D0*itgl-3.D0*iglst1)

      return

      end

c -------------------------------------------------------------------- c

      double precision function NS_gamcfdec(amst1,amst2,amt,amsb1,amsb2,
     .     amb,amgl,amsq,amuv,scala)

      implicit double precision (a-h,k-z)

      mur = scala

      NS_gamcfdec = -(-dlog(mur**2/amuv**2) )*3.D0 - 1.D0

      NS_gamcfdec = NS_gamcfdec + 8.D0/3.D0

c maggie changed with respect to the paper 26/3/03
      NS_gamcfdec = NS_gamcfdec + 4.D0*( 4.D0/12.D0*dlog(mur**2/amsq**2)
     .     + 1.D0/24.D0*dlog(mur**2/amsb1**2) + 1.D0/24.D0*
     .     dlog(mur**2/amsb2**2) 
     .     + 1.D0/24.D0*dlog(mur**2/amst1**2) + 1.D0/24.D0*
     .     dlog(mur**2/amst2**2) + 1.D0/6.D0*dlog(mur**2/amt**2)
     .     + 1.D0/2.D0*dlog(mur**2/amgl**2) )
c maggie changed with respect to the paper 26/3/03

      return

      end

c -------------------------------------------------------------------- c

      subroutine NS_fifafunctions(fffunc,fafunc,fbfunc)

      implicit double precision (a-h,k-z)

      complex*16 SD_C03,SD_C0_lam

      dimension fffunc(3),fafunc(3)
      DOUBLE PRECISION SD_B02

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_relmasses/amst1,amst2,amgl,amt

      pi = 4.D0*datan(1.D0)

      fffunc(1) = 2.D0*(amt**2+amgl**2)*
     .   SD_B02(amst1**2,amgl,amt,amuv**2) + 
     .   (amst1**2+amt**2+amgl**2)*SD_B02(amst1**2,lamv,amst1,amuv**2)
     .   + 2.D0*(amgl**2-amst1**2)*SD_B02(amt**2,lamv,amt,amuv**2)
     .   - 2.D0*amt**2*SD_B02(amt**2,amgl,amst2,amuv**2)
     .   - 4.D0*amgl**2*SD_B02(amgl**2,amt,amst1,amuv**2)
     .   + 4.D0*amgl**2*(amst1**2-amgl**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   + 2.D0*amt**2*(amst1**2+amst2**2-2.D0*amt**2)* 
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fffunc(2) = -2.D0*SD_B02(amst1**2,lamv,amst1,amuv**2) 
     .   -2.D0*SD_B02(amt**2,lamv,amt,amuv**2) 
     .   -4.D0*SD_B02(amst1**2,amgl,amt,amuv**2) 
     .   +2.D0*SD_B02(amt**2,amgl,amst1,amuv**2) 
     .   +4.D0*SD_B02(amgl**2,amt,amst1,amuv**2)
     .   +4.D0*(amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   +2.D0*(amst1**2-amst2**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fffunc(3) = 1.D0/(amgl**2*amt**2)*( 2.D0*amt**2*(
     .   SD_B02(amt**2,amgl,amst2,amuv**2) - 
     .   SD_B02(amt**2,amgl,amst1,amuv**2) ) + 
     .   (amgl**2+amt**2-amst1**2)*( (amgl**2+amt**2-amst1**2) - 4.D0*
     .   amt**2 )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) ) 
     .   -((amgl**2+amt**2-amst1**2)*(amgl**2+amt**2-amst2**2) - 
     .   2.D0*amt**2*(amgl**2+amt**2-amst1**2) - 2.D0*amt**2*
     .   (amgl**2+amt**2-amst2**2) )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) ) )

      fafunc(1) = -2.D0*(amgl**2+amt**2-amst1**2) + 
     .   4.D0*(amt**2-amst1**2)*SD_B02(amgl**2,lamv,amgl,amuv**2)
     .   +2.D0*amt**2*( SD_B02(amt**2,amgl,amst2,amuv**2) -
     .   SD_B02(amt**2,amgl,amst1,amuv**2) ) + 4.D0*amgl**2*
     .   SD_B02(amgl**2,amt,amst1,amuv**2)  
     .   + 4.D0*amgl**2*(amgl**2-amst1**2)* 
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   -2.D0*amt**2*(amst1**2+amst2**2-2.D0*amt**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fafunc(2) = 4.D0 - 4.D0*SD_B02(amgl**2,lamv,amgl,amuv**2)
     .   -4.D0*SD_B02(amgl**2,amt,amst1,amuv**2) 
     .   -4.D0*(amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   -2.D0*(amst1**2-amst2**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fafunc(3) = 1.D0/(amgl**2*amt**2)*( 
     .   2.D0*amt**2*( SD_B02(amt**2,amgl,amst1,amuv**2) 
     .   - SD_B02(amt**2,amgl,amst2,amuv**2) ) 
     .   + (amgl**2+amt**2-amst1**2)*(4.D0*amt**2-(amgl**2+amt**2-
     .   amst1**2))*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   + ( (amgl**2+amt**2-amst1**2)*(amgl**2+amt**2-amst2**2) 
     .   -2.D0*amt**2*(amgl**2+amt**2-amst1**2) - 2.D0*amt**2*
     .   (amgl**2+amt**2-amst2**2) )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) ))

      fbfunc = 3.D0*( (amt**2+amst1**2-amgl**2)*
     .   dreal( SD_C0_lam(amst1,amt,amgl,lamv) )
     .   - (amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C0_lam(amt,amgl,amst1,lamv) )
     .   - (amst1**2+amgl**2-amt**2)*
     .   dreal( SD_C0_lam(amgl,amst1,amt,lamv) ) )
     .   - 8.D0/3.D0*(amt**2+amst1**2-amgl**2)*
     .   dreal( SD_C0_lam(amst1,amt,amgl,lamv) ) 

      do i=1,3,1
         fffunc(i) = -1.D0/16.D0/pi**2*fffunc(i)
         fafunc(i) = -1.D0/16.D0/pi**2*fafunc(i)
      end do
         fbfunc    = -1.D0/16.D0/pi**2*fbfunc

      end

c -------------------------------------------------------------------- c
c ---------- The A function for the higher order corrections --------- c

      double precision function NS_A01(s,mu2)

      implicit double precision (a-h,k-z)

      NS_A01 = s*(-dlog(s/mu2)+1.D0)

      return

      end

c -------------------------------------------------------------------- c
c ---------------------- The divergent piece of A01 ------------------ c

      double precision function NS_A01_DIV(s,mu2)

      implicit double precision (a-h,k-z)

      NS_A01_DIV = s*dlog(mu2)

      return

      end

c -------------------------------------------------------------------- c
c -------- The function B1 for the higher order corrections ---------- c

      double precision function NS_B1(s,m1,m2,mu2)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02

      NS_B1 = 1.D0/2.D0/s*( NS_A01(m1**2,mu2)-NS_A01(m2**2,mu2)
     .     +(m2**2-m1**2-s)*SD_B02(s,m1,m2,mu2) )
      
      return 

      end 

c -------------------------------------------------------------------- c
c ----------------------- The divergent piece of B1 ------------------ c

      double precision function NS_B1_DIV(s,m1,m2,mu2)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV

      NS_B1_DIV = 1.D0/2.D0/s*( m1**2*log(mu2)-m2**2*log(mu2)
     .     +(m2**2-m1**2-s)*SD_B02_DIV(s,m1,m2,mu2) )
      
      return 

      end 

c -------------------------------------------------------------------- c
c ------------------- The derivative of B1: dB1/ ds ------------------ c

      double precision function SD_BP1(s,m1,m2,mu2)

      implicit double precision (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02

      SD_BP1 = -1.D0/2.D0/s**2*( NS_A01(m1**2,mu2)-
     .     NS_A01(m2**2,mu2)+(m2**2-m1**2)*SD_B02(s,m1,m2,mu2) )
     .     +1.D0/2.D0/s*(m2**2-m1**2)*SD_BP02(s,m1,m2,mu2)
     .     -1.D0/2.D0*SD_BP02(s,m1,m2,mu2)

      return 

      end 

c -------------------------------------------------------------------- c
c ------- The C function for a small mass lambda: -------------------- c
c ------- C0(msqp**2,msq**2,mphi**2,msqp,lamv,msq) =     ------------- c
c ------- C0(msq**2,mphi**2,msqp**2,lamv,msq,msqp) =     ------------- c
c - 4.2.2003 M. Muehlleitner ----------------------------------------- c
c -------------------------------------------------------------------- c

      function SD_C0_lam(msqp,msq,mphi,lamv)

      implicit real*8 (a-h,o-z)

      real*8 m1,m2,lamv,msq,msqp,mphi

      complex*16 SD_C0_lam,NS_ccspen,dlxs,dlfc,ieps,ys,xs

      ieps = dcmplx(0.d0,1.d-17)
      pi   = 4.D0*datan(1.D0)

      m1 = msqp
      m2 = msq
      s  = mphi**2

      ys = cdsqrt(1.D0 - (4.D0*m1*m2)/(s-(m1-m2)**2+ieps))

      xs = (ys-1.D0)/(ys+1.D0)

      dlxs = cdlog(xs)
      dlfc = cdlog(1.D0-xs**2)

      SD_C0_lam = xs/(m1*m2*(1.D0-xs**2))*( 
     .     dlxs*(-dlog(lamv**2/(m1*m2))-1.D0/2.D0*dlxs+
     .     2.D0*dlfc ) +
     .     NS_ccspen(1.D0-xs*m1/m2) + 
     .     NS_ccspen(1.D0-xs*m2/m1) + NS_ccspen(xs**2) +
     .     1.D0/2.D0*(dlog(m1/m2))**2 - pi**2/6.D0 )

      return

      end

c -------------------------------------------------------------------- c
c --      The B and C functions for the higher order corrections ----- c
c --      Spence function.                                       ----- c
c -- taken from hdecay.f Version 3.0,                            ----- c
c -- authors: A.Djouadi, J.Kalinowski and M.Spira                ----- c
c -------------------------------------------------------------------- c

      function NS_kappa(a,b,c,d)

      real*8 a,b,c,d
      complex*16 NS_kappa,ieps

      ieps = dcmplx(0.d0,1.d-17)

      if(A.eq.B) then 
         NS_KAPPA = cdsqrt((C*(C-4.D0*A))*(1+IEPS*D))
      elseif(B.eq.C) then
         NS_KAPPA = cdsqrt((A*(A-4.D0*B))*(1+IEPS*D))
      elseif(A.eq.C) then
         NS_KAPPA = cdsqrt((B*(B-4.D0*A))*(1+IEPS*D))
      else
         NS_KAPPA = CDSQRT((A**2+B**2+C**2-2*(A*B+A*C+B*C))
     .        * (1+IEPS*D))
      endif

      end

************************************************************************
      FUNCTION SD_C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
      COMPLEX*16 SD_C03,NS_CCSPEN,NS_ETA,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM,CX,CY
C     REAL*8 NS_KAPPA
      COMPLEX*16 NS_KAPPA
c maggie changed 17/2/03
      COMPLEX*16 ALPHA
      EPS = 1.D-8*(P1+P2+P3)
      IM = DCMPLX(0.D0,1.D0)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IEPS = DCMPLX(0.D0,1.D-17)
c     IEPS = DCMPLX(0.D0,1.D-20)
      PI = 4*DATAN(1.D0)
      XX = 0.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ALPHA  = NS_KAPPA(Q10,Q21,Q20,1.D0)
      ALP(0) = NS_KAPPA(Q21,SM1,SM2,DSIGN(1.D0,Q21))
      ALP(1) = NS_KAPPA(Q20,SM2,SM0,DSIGN(1.D0,Q20))
      ALP(2) = NS_KAPPA(Q10,SM0,SM1,DSIGN(1.D0,Q10))
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
      CDUM=0.D0
      DO I=0,2
       DO J=1,2
        CDUM = CDUM+NS_CCSPEN((Y0(I)-1d0)/Y(I,J))
     .         -NS_CCSPEN(Y0(I)/Y(I,J))
        CX = NS_ETA(1d0-X(I,J),1d0/Y(I,J))
        IF(CX.NE.0.D0)THEN
         CDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = NS_ETA(-X(I,J),1d0/Y(I,J))
        IF(CY.NE.0.D0)THEN
         CDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = NS_ETA(-X(I,1),-X(I,2))
       IF(CX.NE.0.D0)THEN
        CDUM = CDUM - CX*CDLOG((1d0-Y0(I))/(-Y0(I)))
       ENDIF
       CY = NS_ETA(Y(I,1),Y(I,2))
       IF(CY.NE.0.D0)THEN
        CDUM = CDUM + CY*CDLOG((1d0-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2d0*PI*IM*CDLOG((1d0-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      SD_C03 = CDUM/ALPHA
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C        SUBROUTINE CALCULATING THE FINITE REAL PART OF THE            C
C          GENERAL MASSIVE TWO POINT FUNCTION                          C
C                                                                      C
C           SD_B02(P.P,M1,M2,MU**2)                                    C
C           SD_BP02(P.P,M1,M2,MU**2)                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c -------------------------------------------------------------------- c

      real*8 function SD_B02(s,m1,m2,mu2)

      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2 

      m12 = m1**2 
      m22 = m2**2 

      if(s.eq.m22) then
         zkappa = cdsqrt(dcmplx(m12*(m12-4.D0*s)))
      elseif(s.eq.m12) then
         zkappa = cdsqrt(dcmplx(m22*(m22-4.D0*s)))
      elseif(m12.eq.m22) then
         zkappa = cdsqrt(dcmplx(s*(s-4.D0*m12)))
      else
         zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     .        -2.D0*(s*m12+s*m22+m12*m22)))
      endif

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            SD_B02=-dlog(m12/mu2) 
         else
            SD_B02=1.D0 - m12/(m12-m22)*dlog(m12/mu2)
     .                 + m22/(m12-m22)*dlog(m22/mu2) 
         endif
      else 
         if ((m12.eq.0.D0).and.(m22.eq.0.D0)) then 
            SD_B02=2.D0 - dlog(s/mu2)
         elseif ((m12.eq.s).and.(m22.eq.0.D0)) then 
            SD_B02=2.D0 - dlog(m12/mu2)
         elseif ((m22.eq.s).and.(m12.eq.0.D0)) then 
            SD_B02=2.D0 - dlog(m22/mu2)
         elseif (m12.eq.0.D0) then
            SD_B02=2.D0 - (s-m22)/s*dlog( dabs(m22-s)/m22 )
     .                 - dlog(m22/mu2)
         elseif (m22.eq.0.D0) then
            SD_B02=2.D0 - (s-m12)/s*dlog( dabs(m12-s)/m12 ) 
     .                 - dlog(m12/mu2)
         else
            x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
            x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
            SD_B02=dreal( 2.D0+ dlog(mu2/m22) 
     .                       + x1*cdlog(1.D0-1.D0/x1) 
     .                       + x2*cdlog(1.D0-1.D0/x2))
         endif
      endif 

      return
      end

c -------------------------------------------------------------------- c

      real*8 function SD_BP02(s,m1,m2,mu2)
      
      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2
      
      m12 = m1**2
      m22 = m2**2 

      if(s.eq.m22) then
         zkappa = cdsqrt(dcmplx(m12*(m12-4.D0*s)))
      elseif(s.eq.m12) then
         zkappa = cdsqrt(dcmplx(m22*(m22-4.D0*s)))
      elseif(m12.eq.m22) then
         zkappa = cdsqrt(dcmplx(s*(s-4.D0*m12)))
      else
         zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     .        -2.D0*(s*m12+s*m22+m12*m22)))
      endif

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            SD_BP02=1.D0/(6.D0*m12)
         else
            SD_BP02=( (m12+m22)/2.D0 
     .        - m12*m22/(m12-m22)*dlog(m12/m22) )/(m12-m22)**2 
         endif
      elseif ((s.eq.m12).and.(m22.eq.0.D0)) then 
         SD_BP02=( -1.D0 + dlog(m12/mu2)/2.D0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0.D0)) then 
         SD_BP02=( -1.D0 + dlog(m22/mu2)/2.D0 )/m22
      elseif (m22.eq.0.D0) then
         if(m12.ge.s) then
            SD_BP02=( -1.D0 - m12/s*dlog((m12-s)/m12) )/s
         elseif(m12.lt.s) then
            SD_BP02=( -1.D0 - m12/s*dlog((-m12+s)/m12) )/s
         endif
      else 
         x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
         x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
         SD_BP02=dreal( -1.D0 + ( x1*(1.D0-x1)*cdlog(1.D0-1.D0/x1)
     .                     - x2*(1.D0-x2)*cdlog(1.D0-1.D0/x2) )  
     .                                                  /(x1-x2) )/s
      endif 

      return
      end

************************************************************************
        FUNCTION NS_ETA(C1,C2)
************************************************************************
*       COMPLEX ETA-FUNKTION                                           *
*----------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       *
************************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 NS_ETA,C1,C2
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            NS_ETA = DCMPLX(0D0,2D0*PI)
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            NS_ETA = DCMPLX(0D0,-2D0*PI)
        ELSE                                                           
            NS_ETA = DCMPLX(0D0)
        END IF                                                         
        END                                  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION NS_CCSPEN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                      C
C----------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 NS_CCSPEN,W,SUM,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)/
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
c      write(*,*) 'z:',z
      Z =Z*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c      write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
        NS_CCSPEN=-CDLOG(1D0-Z)
        RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
        NS_CCSPEN=1.64493406684822643D0
        RETURN
      END IF
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CDLOG(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 2
      DO 1 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2
      SUM=SUM+U*B(K)
 1    CONTINUE
 2    NS_CCSPEN=SUM
      RETURN
10    W=-CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 12

      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12
      SUM=SUM+U*B(K)
11    CONTINUE
12    NS_CCSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 22
      DO 21 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22
      SUM=SUM+U*B(K)
21    CONTINUE
22    NS_CCSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
      RETURN
30    W=CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 32
      DO 31 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32
      SUM=SUM+U*B(K)
31    CONTINUE
32    NS_CCSPEN=SUM+3.28986813369645287D0
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
      END

c ==================================================================== c
c                     The C11 and C12 functions                        c
c                    C_mu = p1_mu*C11 + p2_mu*C12                      c
c  11.2.03 M.Muehlleitner                                              c
c  p1,p2,p3 squared external momenta                                   c
c ==================================================================== c

      real*8 function NS_C1(p1,p2,p3,m1,m2,m3,amuv)

      implicit real*8 (a-h,o-z)
      real*8 m1,m2,m3
      complex*16 SD_C03

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2- (p3**2-p1**2-p2**2)**2/4.D0

      r1 = 1.D0/2.D0*( 
     .     f1*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1.D0/2.D0*( 
     .     f2*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      NS_C1 = 1.D0/den*( p2**2*r1 - (p3**2-p1**2-p2**2)/2.D0*r2 )

      return 

      end

c -------------------------------------------------------------------- c

      real*8 function SD_C2(p1,p2,p3,m1,m2,m3,amuv)

      implicit real*8 (a-h,o-z)
      Double Precision SD_B02
      real*8 m1,m2,m3
      complex*16 SD_C03

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2 - (p3**2-p1**2-p2**2)**2/4.D0

      r1 = 1.D0/2.D0*( 
     .     f1*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1.D0/2.D0*( 
     .     f2*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      SD_C2 = 1.D0/den*( -(p3**2-p1**2-p2**2)/2.D0*r1 + p1**2*r2 )

      return 

      end

c -------------------------------------------------------------------- c

      real*8 function NS_C1_lam(p1,p2,p3,m1,m2,m3,amuv,lamv)

      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION SD_B02
      real*8 m1,m2,m3,lamv
      complex*16 SD_C0_lam

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2- (p3**2-p1**2-p2**2)**2/4.D0

      r1 = 1.D0/2.D0*( f1*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1.D0/2.D0*( f2*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      NS_C1_lam = 1.D0/den*( p2**2*r1 - (p3**2-p1**2-p2**2)/2.D0*r2 )

      return 

      end

c -------------------------------------------------------------------- c

      real*8 function SD_C2_lam(p1,p2,p3,m1,m2,m3,amuv,lamv)

      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION SD_B02
      real*8 m1,m2,m3,lamv
      complex*16 SD_C0_lam

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2 - (p3**2-p1**2-p2**2)**2/4.D0

      r1 = 1.D0/2.D0*( f1*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1.D0/2.D0*( f2*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      SD_C2_lam = 1.D0/den*( -(p3**2-p1**2-p2**2)/2.D0*r1 + p1**2*r2 )

      return 

      end
      
c -------------------------------------------------------------------- c
c            The divergent pieces of the B functions                   c
c -------------------------------------------------------------------- c

      real*8 function SD_B02_DIV(s,m1,m2,mu2)

      implicit none 

      real*8 s,m1,m2,mu2,m12,m22 

      m12 = m1**2 
      m22 = m2**2 

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            SD_B02_DIV=dlog(mu2) 
         else
            SD_B02_DIV= + m12/(m12-m22)*dlog(mu2)
     .               - m22/(m12-m22)*dlog(mu2) 
         endif
      else 
         if ((m12.eq.0.D0).and.(m22.eq.0.D0)) then 
            SD_B02_DIV= dlog(mu2)
         elseif ((m12.eq.s).and.(m22.eq.0.D0)) then 
            SD_B02_DIV= dlog(mu2)
         elseif ((m22.eq.s).and.(m12.eq.0.D0)) then 
            SD_B02_DIV= dlog(mu2)
         elseif (m12.eq.0.D0) then
            SD_B02_DIV= dlog(mu2)
         elseif (m22.eq.0.D0) then
            SD_B02_DIV= dlog(mu2)
         else
            SD_B02_DIV= dlog(mu2) 
         endif
      endif 

      return
      end

c -------------------------------------------------------------------- c

      real*8 function SD_BP02_DIV(s,m1,m2,mu2)
      
      implicit none 

      real*8 s,m1,m2,mu2,m12,m22 
      
      m12 = m1**2
      m22 = m2**2 

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            SD_BP02_DIV=0.D0
         else
            SD_BP02_DIV=0.D0
         endif
      elseif ((s.eq.m12).and.(m22.eq.0.D0)) then 
         SD_BP02_DIV=( - dlog(mu2)/2.D0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0.D0)) then 
         SD_BP02_DIV=( - dlog(mu2)/2.D0 )/m22
      else 
         SD_BP02_DIV= 0.D0
      endif 

      return
      end

c -------------------- Derivatives of couplings ---------------------- c
c -------------------------------------------------------------------- c
c                   H+ - stop1/2 - sbottom1/2 couplings                c
c -------------------------------------------------------------------- c
      subroutine NS_hcsbotstopderiv(gcdmtr,gcdmbr,gcdabr,gcdatr,
     .                      gcdthtr,gcdthbr)
      IMPLICIT NONE
*
      DOUBLE PRECISION chctbdt(2,2),chctbdb(2,2),chctbab(2,2),
     .     chctbat(2,2),
     .     chctbtt(2,2),gcdthtr(2,2),gcdthbr(2,2),chctbbb(2,2),
     .     gcdmtr(2,2),gcdmbr(2,2),gcdabr(2,2),gcdatr(2,2)
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION au,ad,al,amu
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION b,tgbet
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION tmpt,s11t,s12t,s21t,s22t,tmpb,s11b,s12b,
     . s21b,s22b,s11ab,s12ab,s21ab,s22ab,s11at,s12at,s21at,s22at,
     .ctt,stt,cbb,sbb,s11,s12,s21,s22
C
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/NS_trilin_mu/au,ad,al,amu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_tanb/tanbeta
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
**
      b     = datan(tanbeta)
      tgbet = tanbeta

c ---- derivative d/dmtop ----

      tmpt = 1.D0/dsqrt(2.D0)/amw/dsin(b)

      s11t = 1.d0/dsqrt(2.d0)*amw*(2.D0*tmpt*scalt)*dsin(2.d0*b)
      s12t = 0.D0
      s21t = tmpt*dsin(b)*(au/tgbet+amu)
      s22t = dsqrt(2.d0)*amw*scalb*tmpt

      chctbdt(1,1)=(-ct*cb*s11t-st*sb*s22t-ct*sb*s12t-st*cb*s21t)
      chctbdt(1,2)=(ct*sb*s11t-ct*cb*s12t+sb*st*s21t-st*cb*s22t)
      chctbdt(2,1)=(st*cb*s11t+st*sb*s12t-ct*cb*s21t-ct*sb*s22t)
      chctbdt(2,2)=(-st*sb*s11t+st*cb*s12t+ct*sb*s21t-ct*cb*s22t)

      gcdmtr(1,1)=chctbdt(1,1)/amw
      gcdmtr(1,2)=chctbdt(1,2)/amw
      gcdmtr(2,1)=chctbdt(2,1)/amw
      gcdmtr(2,2)=chctbdt(2,2)/amw

c ---- derivative d/dmbottom ----

      tmpb = 1.D0/dsqrt(2.D0)/amw/dcos(b)

      s11b = 1.d0/dsqrt(2.d0)*amw*(2.D0*tmpb*scalb)*dsin(2.d0*b)
      s12b = tmpb*dcos(b)*(ad*tgbet+amu)
      s21b = 0.D0
      s22b = dsqrt(2.d0)*amw*tmpb*scalt

      chctbdb(1,1)=(-ct*cb*s11b-st*sb*s22b-ct*sb*s12b-st*cb*s21b)
      chctbdb(1,2)=(ct*sb*s11b-ct*cb*s12b+sb*st*s21b-st*cb*s22b)
      chctbdb(2,1)=(st*cb*s11b+st*sb*s12b-ct*cb*s21b-ct*sb*s22b)
      chctbdb(2,2)=(-st*sb*s11b+st*cb*s12b+ct*sb*s21b-ct*cb*s22b)

      gcdmbr(1,1)=chctbdb(1,1)/amw
      gcdmbr(1,2)=chctbdb(1,2)/amw
      gcdmbr(2,1)=chctbdb(2,1)/amw
      gcdmbr(2,2)=chctbdb(2,2)/amw

c ---- derivative d/dAb ----

      s11ab = 0.D0
      s12ab = scalb*dcos(b)*tgbet
      s21ab = 0.D0
      s22ab = 0.D0

      chctbab(1,1)=(-ct*cb*s11ab-st*sb*s22ab-ct*sb*s12ab-st*cb*s21ab)
      chctbab(1,2)=(ct*sb*s11ab-ct*cb*s12ab+sb*st*s21ab-st*cb*s22ab)
      chctbab(2,1)=(st*cb*s11ab+st*sb*s12ab-ct*cb*s21ab-ct*sb*s22ab)
      chctbab(2,2)=(-st*sb*s11ab+st*cb*s12ab+ct*sb*s21ab-ct*cb*s22ab)

      gcdabr(1,1)=chctbab(1,1)/amw
      gcdabr(1,2)=chctbab(1,2)/amw
      gcdabr(2,1)=chctbab(2,1)/amw
      gcdabr(2,2)=chctbab(2,2)/amw

c ---- derivative d/dAt ----

      s11at = 0.D0
      s12at = 0.D0
      s21at = scalt*dsin(b)*1.D0/tgbet
      s22at = 0.D0

      chctbat(1,1)=(-ct*cb*s11at-st*sb*s22at-ct*sb*s12at-st*cb*s21at)
      chctbat(1,2)=(ct*sb*s11at-ct*cb*s12at+sb*st*s21at-st*cb*s22at)
      chctbat(2,1)=(st*cb*s11at+st*sb*s12at-ct*cb*s21at-ct*sb*s22at)
      chctbat(2,2)=(-st*sb*s11at+st*cb*s12at+ct*sb*s21at-ct*cb*s22at)

      gcdatr(1,1)=chctbat(1,1)/amw
      gcdatr(1,2)=chctbat(1,2)/amw
      gcdatr(2,1)=chctbat(2,1)/amw
      gcdatr(2,2)=chctbat(2,2)/amw

c ---- derivative d/dtheta_t ----

      ctt = -dsin(thet)
      stt = dcos(thet)

      s11 = 1.d0/dsqrt(2.d0)*amw*(scalb**2+scalt**2-1.D0)*dsin(2.d0*b)
      s12 = scalb*dcos(b)*(ad*tgbet+amu)
      s21 = scalt*dsin(b)*(au/tgbet+amu)
      s22 = dsqrt(2.d0)*amw*scalb*scalt

      chctbtt(1,1)=(-ctt*cb*s11-stt*sb*s22-ctt*sb*s12-stt*cb*s21)
      chctbtt(1,2)=(ctt*sb*s11-ctt*cb*s12+sb*stt*s21-stt*cb*s22)
      chctbtt(2,1)=(stt*cb*s11+stt*sb*s12-ctt*cb*s21-ctt*sb*s22)
      chctbtt(2,2)=(-stt*sb*s11+stt*cb*s12+ctt*sb*s21-ctt*cb*s22)

      gcdthtr(1,1)=chctbtt(1,1)/amw
      gcdthtr(1,2)=chctbtt(1,2)/amw
      gcdthtr(2,1)=chctbtt(2,1)/amw
      gcdthtr(2,2)=chctbtt(2,2)/amw

c ---- derivative d/dtheta_b ----

      cbb = -dsin(theb)
      sbb = dcos(theb)

      s11 = 1.d0/dsqrt(2.d0)*amw*(scalb**2+scalt**2-1.D0)*dsin(2.d0*b)
      s12 = scalb*dcos(b)*(ad*tgbet+amu)
      s21 = scalt*dsin(b)*(au/tgbet+amu)
      s22 = dsqrt(2.d0)*amw*scalb*scalt

      chctbbb(1,1)=(-ct*cbb*s11-st*sbb*s22-ct*sbb*s12-st*cbb*s21)
      chctbbb(1,2)=(ct*sbb*s11-ct*cbb*s12+sbb*st*s21-st*cbb*s22)
      chctbbb(2,1)=(st*cbb*s11+st*sbb*s12-ct*cbb*s21-ct*sbb*s22)
      chctbbb(2,2)=(-st*sbb*s11+st*cbb*s12+ct*sbb*s21-ct*cbb*s22)

      gcdthbr(1,1)=chctbbb(1,1)/amw
      gcdthbr(1,2)=chctbbb(1,2)/amw
      gcdthbr(2,1)=chctbbb(2,1)/amw
      gcdthbr(2,2)=chctbbb(2,2)/amw

      END

c -------------------------------------------------------------------- c
c --------------------------- The counterterms ----------------------- c

      DOUBLE PRECISION function NS_dcounterhc(amsq,amq,theq,ni,amsqp,
     .                   amqp,theqp,nj,mgluino,amuv,amuvdiv,lamv,ic,jc)

      IMPLICIT NONE 
      INTEGER ic,jc,ni,nj
      DOUBLE PRECISION lamv
      DOUBLE PRECISION gctbr(2,2),gcdthtr(2,2),gcdthbr(2,2),gcdmtr(2,2),
     .          gcdmbr(2,2),gcdabr(2,2),gcdatr(2,2)
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION mgluino
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      DOUBLE PRECISION thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION NS_deltaz,NS_deltamqdiv,NS_deltaAq
     .,NS_deltathdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION amtt,ambb,amsq,amq,theq,amuv,amsqp,amqp,
     .theqp,amuvdiv
*
      COMMON/NS_sfmixang/thet,theb,thel,ct,st,cb,sb,cl,sl,
     .cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,MLR,MLL,MNL,
     C     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     C     CST,CSB,CSL,ase2,ase1,asne1
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_hcsbotstop/gctbr
*
      call NS_hcsbotstopderiv(gcdmtr,gcdmbr,gcdabr,gcdatr,gcdthtr,
     .                        gcdthbr)
c --- the running mass ---
      amtt = runmt 
      ambb = runmb 
*
      NS_dcounterhc = 1.D0/2.D0*gctbr(ic,jc)*( 
     .     NS_deltaz(amsq,mgluino,amq,theq,amuv,lamv,ni) +
     .     NS_deltaz(amsqp,mgluino,amqp,theqp,amuv,lamv,nj) ) +
     .     gcdmtr(ic,jc)*amtt*
     .     NS_deltamqdiv(ast1,ast2,mgluino,amtt,thet,amuvdiv,lamv) +
     .     gcdmbr(ic,jc)*ambb*
     .     NS_deltamqdiv(asb1,asb2,mgluino,ambb,theb,amuvdiv,lamv) +
     .     gcdatr(ic,jc)*
     .     NS_deltaAq(amsq,ast1,ast2,mgluino,amtt,thet,amuvdiv,lamv,1) +
     .     gcdabr(ic,jc)*
     .     NS_deltaAq(amsq,asb1,asb2,mgluino,ambb,theb,amuvdiv,lamv,2) +
     .     gcdthtr(ic,jc)*NS_deltathdiv(ast1,ast2,mgluino,amtt,thet,
     .                               amuvdiv) +
     .     gcdthbr(ic,jc)*NS_deltathdiv(asb1,asb2,mgluino,ambb,theb,
     .                               amuvdiv)

      NS_dcounterhc = -dsqrt(2.D0)*amw**2*NS_dcounterhc

      return 

      end
