* This file contains
*
* subroutine twlpyuk(mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,
*     .     vv,DMS,DMP)
* subroutine makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
*     .     Ft,Fb,FA)
* subroutine makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb)
*****************************************************************

      subroutine twlpyuk(mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,
     .     vv,DMS,DMP)

      implicit none

      integer i,j
      real*8 mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,vv,
     .     DMS(3,3),DMP(3,3)
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,t,b,cbe,sbe,ht,hb,pi,k,mu,DMA
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA!1,FA2,FA3

      pi = 4d0*atan(1d0)

      t = mt**2
      b = mb**2
      mu = l*xx

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2-st**2
      c2b = cb**2-sb**2

      Xt = (T1-T2)*s2t/2d0/mt
      Xb = (B1-B2)*s2b/2d0/mb
      At = Xt+mu/tanb
      Ab = Xb+mu*tanb

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      k = 3d0/(16d0*Pi**2)**2


      CALL makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     .     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA)

      DMS(1,1) = .5d0*ht**2*mu**2*s2t**2*F3t
     .    +2d0*hb**2*mb**2*F1b+2d0*hb**2*Ab*mb*s2b*F2b
     .    +.5d0*hb**2*Ab**2*s2b**2*F3b
     .    -2d0*hb*ht*mb*mu*s2t*F4t-ht*hb*mu*Ab*s2t*s2b*F5
     .    +ht**2*tanb*mu*At/(T1-T2)*Ft+hb**2*tanb*mu*Ab/(B1-B2)*Fb

      DMS(1,2) = -ht**2*mu*mt*s2t*F2t-.5d0*ht**2*At*mu*s2t**2*F3t
     .    -hb**2*mu*mb*s2b*F2b-.5d0*hb**2*Ab*mu*s2b**2*F3b
     .    +ht*hb*mb*At*s2t*F4t+hb*ht*mt*Ab*s2b*F4b
     .    +.5d0*ht*hb*s2t*s2b*(At*Ab+mu**2)*F5
     .    +2d0*ht*hb*mt*mb*F6
     .    -ht**2*mu*At/(T1-T2)*Ft-hb**2*mu*Ab/(B1-B2)*Fb

      DMS(2,2) = .5d0*hb**2*mu**2*s2b**2*F3b
     .    +2d0*ht**2*mt**2*F1t+2d0*ht**2*At*mt*s2t*F2t
     .    +.5d0*ht**2*At**2*s2t**2*F3t
     .    -2d0*ht*hb*mt*mu*s2b*F4b-hb*ht*mu*At*s2b*s2t*F5
     .    +ht**2/tanb*mu*At/(T1-T2)*Ft+hb**2/tanb*mu*Ab/(B1-B2)*Fb

      DMS(1,3) = 0d0
*     .     .5d0*ht*s2t**2*l*mu*mt/tanb*F3t
*     .    -ht*l*mt*(At-2d0*mu/tanb)/(T1-T2)*Ft
*     .    -.5d0*hb*s2b**2*l*Ab*mb*tanb*F3b
*     .    -hb*s2b*l*b*tanb*F2b-hb*l*mb*Ab*tanb/(B1-B2)*Fb
*     .    -l*ht*b*s2t*F4t-.5d0*ht*l*mb*s2t*s2b*(Ab-mu*tanb)*F5

      DMS(2,3) = 0d0
*     .    -.5d0*ht*s2t**2*l*At*mt/tanb*F3t-ht*s2t*l*t/tanb*F2t
*     .    -ht*l*mt*At/tanb/(T1-T2)*Ft
*     .    +.5d0*hb*s2b**2*l*mu*mb*tanb*F3b
*     .    -hb*l*mb*(Ab-2d0*mu*tanb)/(B1-B2)*Fb
*     .    -l*hb*t*s2b*F4b-.5d0*hb*l*mt*s2t*s2b*(At-mu/tanb)*F5

      DMS(3,3) = 0d0
*     .     .5d0*l**2*s2t**2*t/tanb**2*F3t
*     .    +l**2*t/tanb*At/mu/(T1-T2)*Ft
*     .    +.5d0*l**2*s2b**2*b*tanb**2*F3b
*     .    +l**2*b*tanb*Ab/mu/(B1-B2)*Fb
*     .    +l**2*mb*mt*s2t*s2b*F5

      DMS(2,1) = DMS(1,2)
      DMS(3,1) = DMS(1,3)
      DMS(3,2) = DMS(2,3)

      DMA = ht**2*mu*At/(T1-T2)*Ft+hb**2*mu*Ab/(B1-B2)*Fb -2d0*ht*hb*FA

      DMP(1,1) = DMA*tanb

      DMP(1,2) = DMA

      DMP(2,2) = DMA/tanb

      DMP(1,3) = 0d0
*     .     l*ht*mt*At/(T1-T2)*Ft
*     .    +l*hb*mb*Ab/(B1-B2)*Fb*tanb
*     .    -2d0*l*hb*mt*FA2

      DMP(2,3) = 0d0
*     .     l*ht*mt*At/(T1-T2)*Ft/tanb
*     .    +l*hb*mb*Ab/(B1-B2)*Fb
*     .    -2d0*l*ht*mb*FA2

      DMP(3,3) = 0d0
*     .     l**2*t*At/mu/(T1-T2)*Ft/tanb
*     .    +l**2*b*Ab/mu/(B1-B2)*Fb*tanb
*     .    -l**2*FA3

      DMP(2,1) = DMP(1,2)
      DMP(3,1) = DMP(1,3)
      DMP(3,2) = DMP(2,3)

      do i=1,3
       do j=1,3
        DMS(i,j) = k*DMS(i,j)
        DMP(i,j) = k*DMP(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      subroutine makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     .     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
     .     Ft,Fb,FA)

      implicit none

      real*8 t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     .     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA!1,FA2,FA3

      real*8 D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      common/listderivtb/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     .     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      real*8 D1b,DB1,DB2,Dc2b,DB1B1,DB2B2,Dbb,Dc2bc2b,DB1b,DB2b,DB1B2,
     .     Dbc2b,DB1c2b,DB2c2b,Dtc2b

      real*8 Xt,Xb,At,Ab


      CALL makederivtb(b,t,A0,B1,B2,T1,T2,s2b,c2b,s2t,c2t,
     .     q,-mu,vv,1d0/tanb)
*     note: the potential was computed with the opposite convention for mu

      D1b = D1t
      DB1 = DT1
      DB2 = DT2
      Dc2b = Dc2t
      DB1B1 = DT1T1
      DB2B2 = DT2T2
      Dbb = Dtt
      Dc2bc2b = Dc2tc2t
      DB1b = DT1t
      DB2b = DT2t
      DB1B2 = DT1T2
      Dbc2b = Dtc2t
      DB1c2b = DT1c2t
      DB2c2b = DT2c2t
      Dtc2b = Dbc2t


      CALL makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,-mu,vv,tanb)

      F1t = Dtt+DT1T1+DT2T2+2d0*(DT1t+DT2t+DT1T2)

      F2t = DT1T1-DT2T2+DT1t-DT2t
     .    -4d0*c2t**2/(T1-T2)*(Dtc2t+DT1c2t+DT2c2t)

      F3t = DT1T1+DT2T2-2d0*DT1T2
     .    -2d0/(T1-T2)*(DT1-DT2)
     .    +16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t+2d0*Dc2t)
     .    -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)

      F4t = DT1b+DT1B1+DT1B2-DT2b-DT2B1-DT2B2
     .    -4d0*c2t**2/(T1-T2)*(DB1c2t+DB2c2t+Dbc2t)

      F1b = Dbb+DB1B1+DB2B2+2d0*(DB1b+DB2b+DB1B2)

      F2b = DB1B1-DB2B2+DB1b-DB2b
     .    -4d0*c2b**2/(B1-B2)*(Dbc2b+DB1c2b+DB2c2b)

      F3b = DB1B1+DB2B2-2d0*DB1B2
     .    -2d0/(B1-B2)*(DB1-DB2)
     .    +16d0*c2b**2/(B1-B2)**2*(c2b**2*Dc2bc2b+2d0*Dc2b)
     .    -8d0*c2b**2/(B1-B2)*(DB1c2b-DB2c2b)

      F4b = DB1t+DT1B1-DT1B2-DB2t+DT2B1-DT2B2
     .    -4d0*c2b**2/(B1-B2)*(DT1c2b+DT2c2b+Dtc2b)

      F5  = DT1B1-DT1B2-DT2B1+DT2B2
     .    +16d0*c2t**2*c2b**2/(T1-T2)/(B1-B2)*Dc2tc2b
     .    -4d0*c2t**2/(T1-T2)*(DB1c2t-DB2c2t)
     .    -4d0*c2b**2/(B1-B2)*(DT1c2b-DT2c2b)

      F6 = Dtb+DT1b+DT2b+DB1t+DB2t
     .    +DT1B1+DT1B2+DT2B1+DT2B2

      Ft = DT1-DT2-4d0*c2t**2/(T1-T2)*Dc2t

      Fb = DB1-DB2-4d0*c2b**2/(B1-B2)*Dc2b

      Xt = (T1-T2)*s2t/2d0/sqrt(t)
      Xb = (B1-B2)*s2b/2d0/sqrt(b)

      At = Xt+mu/tanb
      Ab = Xb+mu*tanb

      FA = Dcptpb/4d0/Sqrt(b)/Sqrt(t)
     .    +4d0*Sqrt(t)*Sqrt(b)*(At*Ab-mu**2)**2
     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
     .    +Sqrt(t)/Sqrt(b)/s2t**2/(T1-T2)**2
     .     *(At**2*Dcpbptt+mu**2/tanb**2*Dcptmptt)
     .    +Sqrt(b)/Sqrt(t)/s2b**2/(B1-B2)**2
     .     *(Ab**2*Dcptptb+mu**2*tanb**2*Dcpbmptb)
     .    +2d0*mu/s2t/s2b/(T1-T2)/(B1-B2)
     .     *(At*tanb*Dspbmptbspbptt+Ab/tanb*Dsptmpttsptptb
     .    -mu*Dsptmpttspbmptb)

*      FA2 = 4d0*Sqrt(t)*Sqrt(b)*(At*Ab-mu**2)*(Ab/tanb+At*tanb-2d0*mu)
*     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
*     .    +Sqrt(t)/Sqrt(b)/tanb/s2t**2/(T1-T2)**2
*     .     *(At*Dcpbptt+mu/tanb*Dcptmptt)
*     .    +Sqrt(b)/Sqrt(t)*tanb/s2b**2/(B1-B2)**2
*     .     *(Ab*Dcptptb+mu*tanb*Dcpbmptb)
*     .    +1d0/s2t/s2b/(T1-T2)/(B1-B2)
*     .     *((At*tanb+mu)*Dspbmptbspbptt+(Ab/tanb+mu)*Dsptmpttsptptb
*     .    -2d0*mu*Dsptmpttspbmptb)
*
*      FA3 = 8d0*t*b*(Ab/tanb+At*tanb-2d0*mu)**2
*     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
*     .    +2d0*t/tanb**2/s2t**2/(T1-T2)**2*(Dcpbptt+Dcptmptt)
*     .    +2d0*b*tanb**2/s2b**2/(B1-B2)**2*(Dcptptb+Dcpbmptb)
*     .    +4d0*Sqrt(t)*Sqrt(b)/s2t/s2b/(T1-T2)/(B1-B2)
*     .     *(Dspbmptbspbptt+Dsptmpttsptptb-Dsptmpttspbmptb)

      end

*
***********************************************************************
*

      subroutine makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     .     q,mu,vv,tanb)

      implicit real*8 (t)
      implicit character (a-s,u-z)

      real*8 t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,mt,mb,Xt,Xb,Yt,Yb,sbe,cbe,mu2,Nc
      real*8 Delt,phi,SLLi2
      external delt,phi,SLLi2

      real*8 D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      common/listderivtb/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     .     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      real*8  Logt,Logb,Logmu2,LogA0,LogT1,LogT2,LogB1,LogB2,
     .     phimu2tT1,phimu2T2T,phiB1tmu2,phiB2tmu2,phiT1bmu2,phiT2bmu2,
     .     phiA0T1T1,phiA0T2T2,phiA0T1T2,phiA0B1B1,phiA0B2B2,
     .     phiA0B1T1,phiA0B2T1,phiA0B1T2,phiA0B2T2,phiA0tt,phiA0bt,
     .     deltmu2tT1,deltmu2T2T,deltB1tmu2,deltB2tmu2,
     .     deltT1bmu2,delT2Tbmu2,deltA0tt,deltA0bt,
     .     deltA0T1T1,deltA0T2T2,deltA0T1T2,deltA0B1B1,deltA0B2B2,
     .     deltA0B1T1,deltA0B2T1,deltA0B1T2,deltA0B2T2,
     .     Li2T1T2,Li2T1B1,Li2T1B2,Li2T2B1,Li2T2B2,Li2bt

      Nc = 3d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      Xt = (T1-T2)*s2t/2d0/mt
      Xb = (B1-B2)*s2b/2d0/mb
      Yt  = Xt-mu/sbe/cbe
      Yb  = Xb-mu/sbe/cbe

      mu2 = mu**2

      Logt = Log(t/q)
      Logb = Log(b/q)
      Logmu2 = Log(mu2/q)
      LogA0 = Log(A0/q)
      LogT1 = Log(T1/q)
      LogT2 = Log(T2/q)
      LogB1 = Log(B1/q)
      LogB2 = Log(B2/q)
      phimu2tT1 = phi(mu2,t,T1)
      phimu2T2T = phi(mu2,t,T2)
      phiB1tmu2 = phi(B1,t,mu2)
      phiB2tmu2 = phi(B2,t,mu2)
      phiT1bmu2 = phi(T1,b,mu2)
      phiT2bmu2 = phi(T2,b,mu2)
      phiA0T1T1 = phi(A0,T1,T1)
      phiA0T1T2 = phi(A0,T1,T2)
      phiA0T2T2 = phi(A0,T2,T2)
      phiA0B1B1 = phi(A0,B1,B1)
      phiA0B2B2 = phi(A0,B2,B2)
      phiA0B1T1 = phi(A0,B1,T1)
      phiA0B1T2 = phi(A0,B1,T2)
      phiA0B2T1 = phi(A0,B2,T1)
      phiA0B2T2 = phi(A0,B2,T2)
      phiA0tt = phi(A0,t,t)
      phiA0bt = phi(A0,b,t)
      deltmu2tT1 = delt(mu2,t,T1)
      deltmu2T2T = delt(mu2,t,T2)
      deltB1tmu2 = delt(B1,t,mu2)
      deltB2tmu2 = delt(B2,t,mu2)
      deltT1bmu2 = delt(T1,b,mu2)
      delT2Tbmu2 = delt(T2,b,mu2)
      deltA0T1T1 = delt(A0,T1,T1)
      deltA0T1T2 = delt(A0,T1,T2)
      deltA0T2T2 = delt(A0,T2,T2)
      deltA0B1B1 = delt(A0,B1,B1)
      deltA0B2B2 = delt(A0,B2,B2)
      deltA0B1T1 = delt(A0,B1,T1)
      deltA0B1T2 = delt(A0,B1,T2)
      deltA0B2T1 = delt(A0,B2,T1)
      deltA0B2T2 = delt(A0,B2,T2)
      deltA0tt = delt(A0,t,t)
      deltA0bt = delt(A0,b,t)
      Li2T1T2 = SLLi2(1d0-T1/T2)
      Li2T1B1 = SLLi2(1d0-T1/B1)
      Li2T2B1 = SLLi2(1d0-T2/B1)
      Li2T1B2 = SLLi2(1d0-T1/B2)
      Li2T2B2 = SLLi2(1d0-T2/B2)
      Li2bt = SLLi2(1d0-b/t)

      tmp1 = 2d0-LogA0-Logb+2d0*Logt
      tmp2 = 2d0-LogB1-Logmu2+2d0*Logt
      tmp3 = 2d0-LogB2-Logmu2+2d0*Logt
      tmp4 = 2d0-Logb-Logmu2+2d0*LogT1
      tmp5 = 2d0-Logb-Logmu2+2d0*LogT2
      tmp6 = 0.25d0*hb**2/c2t-0.25d0*ht**2/c2t
      tmp7 = -(0.25d0*hb**2/c2t)+0.25d0*ht**2/c2t
      tmp8 = c2t**2-0.5d0*((-1d0+Nc)*s2t**2)
      tmp9 = cbe**2*ht**2+hb**2*sbe**2
      tmp10 = cbe**2*hb**2+ht**2*sbe**2
      tmp11 = -(0.5d0*(cbe**2*ht**2*mb*mt*s2b*s2t))-
     .   0.5d0*(hb**2*mb*mt*s2b*s2t*sbe**2)
      tmp12 = 0.5d0*(cbe**2*ht**2*mb*mt*s2b*s2t)+
     .   0.5d0*(hb**2*mb*mt*s2b*s2t*sbe**2)
      tmp13 = -(0.5d0*(cbe**2*hb**2*mb*mt*s2b*s2t))-
     .   0.5d0*(ht**2*mb*mt*s2b*s2t*sbe**2)
      tmp14 = 0.5d0*(cbe**2*hb**2*mb*mt*s2b*s2t)+
     .   0.5d0*(ht**2*mb*mt*s2b*s2t*sbe**2)
      tmp15 = 2d0/t**2-(2d0*(-1d0+Logt))/t**2
      tmp16 = 2d0/t**2-(2d0*Logt)/t**2
      tmp17 = 2d0/T1**2-(2d0*LogT1)/T1**2
      tmp18 = 2d0/T2**2-(2d0*LogT2)/T2**2
      tmp19 = -(0.0625d0*
     .      (B1*(-1d0+LogB1)*(-1d0+LogT1)*T1)/(c2b*c2t))+
     .   0.0625d0*(B2*(-1d0+LogB2)*(-1d0+LogT1)*T1)/(c2b*c2t)+
     .   0.0625d0*(B1*(-1d0+LogB1)*(-1d0+LogT2)*T2)/(c2b*c2t)-
     .   0.0625d0*(B2*(-1d0+LogB2)*(-1d0+LogT2)*T2)/(c2b*c2t)
      tmp20 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Xb)/mt
      tmp21 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)-0.125d0*(mb*Xb)/(c2t*s2b)+
     .   0.125d0*(mt*Xb)/(c2b*s2t)-0.0625d0*Xb**2/(c2b*c2t)
      tmp22 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)+
     .   0.125d0*(mb*Xb)/(c2t*s2b)-0.125d0*(mt*Xb)/(c2b*s2t)+
     .   0.0625d0*Xb**2/(c2b*c2t)
      tmp23 = 0.25d0*Xb/c2b-0.25d0*Xt/c2b
      tmp24 = -(0.25d0*Xb/c2b)+0.25d0*Xt/c2b
      tmp25 = 0.125d0*Xb/c2t**3-0.125d0*Xt/c2t**3
      tmp26 = -(0.125d0*Xb/c2t**3)+0.125d0*Xt/c2t**3
      tmp27 = 0.25d0*Xb/c2t-0.25d0*Xt/c2t
      tmp28 = -(0.25d0*Xb/c2t)+0.25d0*Xt/c2t
      tmp29 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Xt)/mt
      tmp30 = (-2d0*mt*Xt)/s2t-Xt**2
      tmp31 = (2d0*mt*Xt)/s2t-Xt**2
      tmp32 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)+0.125d0*(mb*Xt)/(c2t*s2b)-
     .   0.125d0*(mt*Xt)/(c2b*s2t)-0.0625d0*Xt**2/(c2b*c2t)
      tmp33 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)-
     .   0.125d0*(mb*Xt)/(c2t*s2b)+0.125d0*(mt*Xt)/(c2b*s2t)+
     .   0.0625d0*Xt**2/(c2b*c2t)
      tmp34 = 4d0*t-4d0*mt*s2t*Xt+s2t**2*Xt**2
      tmp35 = 4d0*t+4d0*mt*s2t*Xt+s2t**2*Xt**2
      tmp36 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Yb)/mt
      tmp37 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)-0.125d0*(mb*Yb)/(c2t*s2b)+
     .   0.125d0*(mt*Yb)/(c2b*s2t)-0.0625d0*Yb**2/(c2b*c2t)
      tmp38 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)+
     .   0.125d0*(mb*Yb)/(c2t*s2b)-0.125d0*(mt*Yb)/(c2b*s2t)+
     .   0.0625d0*Yb**2/(c2b*c2t)
      tmp39 = 0.25d0*Yb/c2b-0.25d0*Yt/c2b
      tmp40 = -(0.25d0*Yb/c2b)+0.25d0*Yt/c2b
      tmp41 = 0.125d0*Yb/c2t**3-0.125d0*Yt/c2t**3
      tmp42 = -(0.125d0*Yb/c2t**3)+0.125d0*Yt/c2t**3
      tmp43 = 0.25d0*Yb/c2t-0.25d0*Yt/c2t
      tmp44 = -(0.25d0*Yb/c2t)+0.25d0*Yt/c2t
      tmp45 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Yt)/mt
      tmp46 = (-2d0*mt*Yt)/s2t-Yt**2
      tmp47 = (2d0*mt*Yt)/s2t-Yt**2
      tmp48 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)+0.125d0*(mb*Yt)/(c2t*s2b)-
     .   0.125d0*(mt*Yt)/(c2b*s2t)-0.0625d0*Yt**2/(c2b*c2t)
      tmp49 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)-
     .   0.125d0*(mb*Yt)/(c2t*s2b)+0.125d0*(mt*Yt)/(c2b*s2t)+
     .   0.0625d0*Yt**2/(c2b*c2t)
      tmp50 = 4d0*t-4d0*mt*s2t*Yt+s2t**2*Yt**2
      tmp51 = 4d0*t+4d0*mt*s2t*Yt+s2t**2*Yt**2
      tmp52 = 0.0625d0*(1d0-c2b)/c2t**3-0.0625d0*(1d0+c2b)/c2t**3
      tmp53 = -(0.0625d0*(1d0-c2b)/c2t**3)
     .      +0.0625d0*(1d0+c2b)/c2t**3
      tmp54 = 0.125d0*(1d0-c2b)/c2t-0.125d0*(1d0+c2b)/c2t
      tmp55 = -(0.125d0*(1d0-c2b)/c2t)+0.125d0*(1d0+c2b)/c2t
      tmp56 = 0.25d0*((1d0+c2b)*(1d0-c2t))
     .      +0.25d0*((1d0-c2b)*(1d0+c2t))
      tmp57 = 0.125d0*(1d0-c2t)/c2b-0.125d0*(1d0+c2t)/c2b
      tmp58 = -(0.125d0*(1d0-c2t)/c2b)+0.125d0*(1d0+c2t)/c2b
      tmp59 = 0.25d0*((1d0-c2b)*(1d0-c2t))
     .      +0.25d0*((1d0+c2b)*(1d0+c2t))
      tmp60 = 0.5d0*((1d0+c2b)*hb**2)+0.5d0*((1d0-c2b)*ht**2)
      tmp61 = 0.5d0*((1d0-c2b)*hb**2)+0.5d0*((1d0+c2b)*ht**2)
      tmp62 = 0.5d0*((1d0+c2t)*hb**2)+0.5d0*((1d0-c2t)*ht**2)
      tmp63 = 0.5d0*((1d0-c2t)*hb**2)+0.5d0*((1d0+c2t)*ht**2)
      tmp64 = (A0-b)/b+LogA0-Logb-t/b
      tmp65 = (-LogA0+Logt)*(A0-t)+(-LogA0+Logt)*t
      tmp66 = (A0-b)*(-LogA0+Logb)+(-LogA0-Logb+2d0*Logt)*t
      tmp67 = (-LogB1+Logmu2)*(B1-mu2)+
     .   (-LogB1-Logmu2+2d0*Logt)*t
      tmp68 = (-LogB2+Logmu2)*(B2-mu2)+
     .   (-LogB2-Logmu2+2d0*Logt)*t
      tmp69 = b*(-LogA0+2d0*Logb-Logt)+(LogA0-Logt)*(-A0+t)
      tmp70 = (-LogA0+Logt)*t+(LogA0-Logt)*(-A0+t)
      tmp71 = B1*(2d0*LogB1-Logmu2-Logt)+
     .   (Logmu2-Logt)*(-mu2+t)
      tmp72 = B2*(2d0*LogB2-Logmu2-Logt)+
     .   (Logmu2-Logt)*(-mu2+t)
      tmp73 = Logmu2-Logt-B1/t-(-mu2+t)/t
      tmp74 = Logmu2-Logt-B2/t-(-mu2+t)/t
      tmp75 = 2d0*B1*LogB1-2d0*B2*LogB2-
     .   0.5d0*(deltB1tmu2*phiB1tmu2)/mu2+
     .   0.5d0*(deltB2tmu2*phiB2tmu2)/mu2+
     .   0.5d0*(Logmu2*Logt*(B1-mu2-t))-
     .   0.5d0*(Logmu2*Logt*(B2-mu2-t))+
     .   0.5d0*(LogB1*Logt*(-B1+mu2-t))-
     .   0.5d0*(LogB2*Logt*(-B2+mu2-t))+
     .   0.5d0*(LogB1*Logmu2*(-B1-mu2+t))-
     .   0.5d0*(LogB2*Logmu2*(-B2-mu2+t))-2.5d0*(B1+mu2+t)+
     .   2.5d0*(B2+mu2+t)
      tmp76 = -5d0*b+4d0*b*Logb+Logt**2*(b-t)-5d0*t-
     .   2d0*Li2bt*(-b+t)+Logt*(-2d0*b*Logb+4d0*t)
      tmp77 = -6d0+4d0*Logt+(2d0-2d0*Logt)*Logt+(4d0*t-2d0*Logt*t)/t
      tmp78 = -Logb+Logmu2-(b-mu2)/b-T1/b
      tmp79 = 2d0*(-1d0+LogT1)*T1+2d0*(-1d0+LogT1)**2*T1
      tmp80 = (-LogA0+LogT1)*(A0-T1)+(-LogA0+LogT1)*T1
      tmp81 = (A0-B1)*(-LogA0+LogB1)+
     .   (-LogA0-LogB1+2d0*LogT1)*T1
      tmp82 = (A0-B2)*(-LogA0+LogB2)+
     .   (-LogA0-LogB2+2d0*LogT1)*T1
      tmp83 = (-Logb+Logmu2)*(b-mu2)+
     .   (-Logb-Logmu2+2d0*LogT1)*T1
      tmp84 = (Logmu2-Logt)*(-mu2+t)+
     .   (-Logmu2-Logt+2d0*LogT1)*T1
      tmp85 = B1*(-LogA0+2d0*LogB1-LogT1)+
     .   (LogA0-LogT1)*(-A0+T1)
      tmp86 = B2*(-LogA0+2d0*LogB2-LogT1)+
     .   (LogA0-LogT1)*(-A0+T1)
      tmp87 = (-LogA0+LogT1)*T1+(LogA0-LogT1)*(-A0+T1)
      tmp88 = 2d0*A0*LogA0+2d0*B1*LogB1+2d0*LogT1*T1+
     .   0.5d0*(LogB1*LogT1*(A0-B1-T1))+
     .   0.5d0*(LogA0*LogT1*(-A0+B1-T1))-
     .   0.5d0*(deltA0B1T1*phiA0B1T1)/T1+
     .   0.5d0*(LogA0*LogB1*(-A0-B1+T1))-2.5d0*(A0+B1+T1)
      tmp89 = 2d0*A0*LogA0+2d0*B2*LogB2+2d0*LogT1*T1+
     .   0.5d0*(LogB2*LogT1*(A0-B2-T1))+
     .   0.5d0*(LogA0*LogT1*(-A0+B2-T1))-
     .   0.5d0*(deltA0B2T1*phiA0B2T1)/T1+
     .   0.5d0*(LogA0*LogB2*(-A0-B2+T1))-2.5d0*(A0+B2+T1)
      tmp90 = b*(2d0*Logb-Logmu2-LogT1)+
     .   (Logmu2-LogT1)*(-mu2+T1)
      tmp91 = (-Logmu2+2d0*Logt-LogT1)*t+
     .   (Logmu2-LogT1)*(-mu2+T1)
      tmp92 = 2d0*b*Logb+2d0*Logmu2*mu2+2d0*LogT1*T1-
     .   0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .   0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T1))-2.5d0*(b+mu2+T1)
      tmp93 = 2d0*A0*LogA0-A0*LogA0*LogT1+4d0*LogT1*T1+
     .   0.5d0*(LogT1**2*(A0-2d0*T1))-
     .   0.5d0*(deltA0T1T1*phiA0T1T1)/T1-2.5d0*(A0+2d0*T1)
      tmp94 = -10d0*T1+4d0*LogT1*T1+LogT1*(4d0*T1-2d0*LogT1*T1)
      tmp95 = -6d0+4d0*LogT1+(2d0-2d0*LogT1)*LogT1+
     .   (4d0*T1-2d0*LogT1*T1)/T1
      tmp96 = (-LogA0+2d0*LogT1-LogT2)*T1+
     .   (-LogA0+LogT2)*(A0-T2)
      tmp97 = -Logb+Logmu2-(b-mu2)/b-T2/b
      tmp98 = 2d0*(-1d0+LogT2)*T2+2d0*(-1d0+LogT2)**2*T2
      tmp99 = (-LogA0+LogT2)*(A0-T2)+(-LogA0+LogT2)*T2
      tmp100 = (A0-B1)*(-LogA0+LogB1)+
     .   (-LogA0-LogB1+2d0*LogT2)*T2
      tmp101 = (A0-B2)*(-LogA0+LogB2)+
     .   (-LogA0-LogB2+2d0*LogT2)*T2
      tmp102 = (-Logb+Logmu2)*(b-mu2)+
     .   (-Logb-Logmu2+2d0*LogT2)*T2
      tmp103 = (Logmu2-Logt)*(-mu2+t)+
     .   (-Logmu2-Logt+2d0*LogT2)*T2
      tmp104 = (LogA0-LogT1)*(-A0+T1)+
     .   (-LogA0-LogT1+2d0*LogT2)*T2
      tmp105 = B1*(-LogA0+2d0*LogB1-LogT2)+
     .   (LogA0-LogT2)*(-A0+T2)
      tmp106 = B2*(-LogA0+2d0*LogB2-LogT2)+
     .   (LogA0-LogT2)*(-A0+T2)
      tmp107 = (-LogA0+LogT2)*T2+(LogA0-LogT2)*(-A0+T2)
      tmp108 = 2d0*A0*LogA0+2d0*B1*LogB1+2d0*LogT2*T2+
     .   0.5d0*(LogB1*LogT2*(A0-B1-T2))+
     .   0.5d0*(LogA0*LogT2*(-A0+B1-T2))-
     .   0.5d0*(deltA0B1T2*phiA0B1T2)/T2+
     .   0.5d0*(LogA0*LogB1*(-A0-B1+T2))-2.5d0*(A0+B1+T2)
      tmp109 = 2d0*A0*LogA0+2d0*B2*LogB2+2d0*LogT2*T2+
     .   0.5d0*(LogB2*LogT2*(A0-B2-T2))+
     .   0.5d0*(LogA0*LogT2*(-A0+B2-T2))-
     .   0.5d0*(deltA0B2T2*phiA0B2T2)/T2+
     .   0.5d0*(LogA0*LogB2*(-A0-B2+T2))-2.5d0*(A0+B2+T2)
      tmp110 = b*(2d0*Logb-Logmu2-LogT2)+
     .   (Logmu2-LogT2)*(-mu2+T2)
      tmp111 = (-Logmu2+2d0*Logt-LogT2)*t+
     .   (Logmu2-LogT2)*(-mu2+T2)
      tmp112 = 2d0*b*Logb+2d0*Logmu2*mu2+2d0*LogT2*T2-
     .   0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT2*(b-mu2-T2))+
     .   0.5d0*(Logb*LogT2*(-b+mu2-T2))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T2))-2.5d0*(b+mu2+T2)
      tmp113 = 2d0*LogT1*T1-2d0*LogT2*T2-
     .   0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .   0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .   0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T1))-2.5d0*(b+mu2+T1)-
     .   0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .   0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2)
      tmp114 = 2d0*A0*LogA0-A0*LogA0*LogT2+4d0*LogT2*T2+
     .   0.5d0*(LogT2**2*(A0-2d0*T2))-
     .   0.5d0*(deltA0T2T2*phiA0T2T2)/T2-2.5d0*(A0+2d0*T2)
      tmp115 = -10d0*T2+4d0*LogT2*T2+LogT2*(4d0*T2-2d0*LogT2*T2)
      tmp116 = -6d0+4d0*LogT2+(2d0-2d0*LogT2)*LogT2+
     .   (4d0*T2-2d0*LogT2*T2)/T2
      tmp117 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Xb)/mt
      tmp118 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Xb)/mt
      tmp119 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Xb)/mt
      tmp120 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Xb)-0.5d0*((1d0-c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Xb**2)
      tmp121 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Xb)-0.5d0*((1d0+c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Xb**2)
      tmp122 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xb)/c2t+0.25d0*((1d0-c2b)*mt*Xb)/s2t-
     .   0.125d0*((1d0-c2b)*Xb**2)/c2t
      tmp123 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xb)/c2t-0.25d0*((1d0-c2b)*mt*Xb)/s2t+
     .   0.125d0*((1d0-c2b)*Xb**2)/c2t
      tmp124 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xb)/c2t+0.25d0*((1d0+c2b)*mt*Xb)/s2t-
     .   0.125d0*((1d0+c2b)*Xb**2)/c2t
      tmp125 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xb)/c2t-0.25d0*((1d0+c2b)*mt*Xb)/s2t+
     .   0.125d0*((1d0+c2b)*Xb**2)/c2t
      tmp126 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Xb)+0.5d0*((1d0-c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Xb**2)
      tmp127 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Xb)+0.5d0*((1d0+c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Xb**2)
      tmp128 = 0.5d0*((1d0+c2b)*Xb)+0.5d0*((1d0-c2b)*Xt)
      tmp129 = 0.5d0*((1d0-c2b)*Xb)+0.5d0*((1d0+c2b)*Xt)
      tmp130 = 0.5d0*((1d0+c2t)*Xb)+0.5d0*((1d0-c2t)*Xt)
      tmp131 = 0.5d0*((1d0-c2t)*Xb)+0.5d0*((1d0+c2t)*Xt)
      tmp132 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Xt)/mt
      tmp133 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Xt)/mt
      tmp134 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Xt)/mt
      tmp135 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Xt)+0.5d0*((1d0-c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Xt**2)
      tmp136 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Xt)+0.5d0*((1d0+c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Xt**2)
      tmp137 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xt)/c2t-0.25d0*((1d0-c2b)*mt*Xt)/s2t-
     .   0.125d0*((1d0-c2b)*Xt**2)/c2t
      tmp138 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xt)/c2t+0.25d0*((1d0-c2b)*mt*Xt)/s2t+
     .   0.125d0*((1d0-c2b)*Xt**2)/c2t
      tmp139 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xt)/c2t-0.25d0*((1d0+c2b)*mt*Xt)/s2t-
     .   0.125d0*((1d0+c2b)*Xt**2)/c2t
      tmp140 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xt)/c2t+0.25d0*((1d0+c2b)*mt*Xt)/s2t+
     .   0.125d0*((1d0+c2b)*Xt**2)/c2t
      tmp141 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Xt)-0.5d0*((1d0-c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Xt**2)
      tmp142 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Xt)-0.5d0*((1d0+c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Xt**2)
      tmp143 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Yb)/mt
      tmp144 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Yb)/mt
      tmp145 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Yb)/mt
      tmp146 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Yb)-0.5d0*((1d0-c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Yb**2)
      tmp147 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Yb)-0.5d0*((1d0+c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Yb**2)
      tmp148 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yb)/c2t+0.25d0*((1d0-c2b)*mt*Yb)/s2t-
     .   0.125d0*((1d0-c2b)*Yb**2)/c2t
      tmp149 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yb)/c2t-0.25d0*((1d0-c2b)*mt*Yb)/s2t+
     .   0.125d0*((1d0-c2b)*Yb**2)/c2t
      tmp150 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yb)/c2t+0.25d0*((1d0+c2b)*mt*Yb)/s2t-
     .   0.125d0*((1d0+c2b)*Yb**2)/c2t
      tmp151 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yb)/c2t-0.25d0*((1d0+c2b)*mt*Yb)/s2t+
     .   0.125d0*((1d0+c2b)*Yb**2)/c2t
      tmp152 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Yb)+0.5d0*((1d0-c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Yb**2)
      tmp153 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Yb)+0.5d0*((1d0+c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Yb**2)
      tmp154 = 0.5d0*((1d0+c2b)*Yb)+0.5d0*((1d0-c2b)*Yt)
      tmp155 = 0.5d0*((1d0-c2b)*Yb)+0.5d0*((1d0+c2b)*Yt)
      tmp156 = 0.5d0*((1d0+c2t)*Yb)+0.5d0*((1d0-c2t)*Yt)
      tmp157 = 0.5d0*((1d0-c2t)*Yb)+0.5d0*((1d0+c2t)*Yt)
      tmp158 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Yt)/mt
      tmp159 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Yt)/mt
      tmp160 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Yt)/mt
      tmp161 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Yt)+0.5d0*((1d0-c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Yt**2)
      tmp162 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Yt)+0.5d0*((1d0+c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Yt**2)
      tmp163 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yt)/c2t-0.25d0*((1d0-c2b)*mt*Yt)/s2t-
     .   0.125d0*((1d0-c2b)*Yt**2)/c2t
      tmp164 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yt)/c2t+0.25d0*((1d0-c2b)*mt*Yt)/s2t+
     .   0.125d0*((1d0-c2b)*Yt**2)/c2t
      tmp165 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yt)/c2t-0.25d0*((1d0+c2b)*mt*Yt)/s2t-
     .   0.125d0*((1d0+c2b)*Yt**2)/c2t
      tmp166 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yt)/c2t+0.25d0*((1d0+c2b)*mt*Yt)/s2t+
     .   0.125d0*((1d0+c2b)*Yt**2)/c2t
      tmp167 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Yt)-0.5d0*((1d0-c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Yt**2)
      tmp168 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Yt)-0.5d0*((1d0+c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Yt**2)
      tmp169 = -Li2T1B1-0.5d0*(-LogB1+LogT1)**2
      tmp170 = -Li2T1B2-0.5d0*(-LogB2+LogT1)**2
      tmp171 = -Li2T1T2-0.5d0*(LogT1-LogT2)**2
      tmp172 = -Li2T2B1-0.5d0*(-LogB1+LogT2)**2
      tmp173 = -Li2T2B2-0.5d0*(-LogB2+LogT2)**2
      tmp174 = -A0+(A0-t)**2/b-t
      tmp175 = -A0+(A0-t)**2/t-t
      tmp176 = -A0+(A0-T1)**2/B1-T1
      tmp177 = -A0+(A0-T1)**2/B2-T1
      tmp178 = -A0+(A0-T1)**2/T1-T1
      tmp179 = -A0-T1+(A0-T1)**2/T2
      tmp180 = -A0+(A0-T2)**2/B1-T2
      tmp181 = -A0+(A0-T2)**2/B2-T2
      tmp182 = -A0+(A0-T2)**2/T2-T2
      tmp183 = cbe**2*hb**2*tmp119+ht**2*sbe**2*tmp133-
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp130)/mt)
      tmp184 = cbe**2*hb**2*tmp118+ht**2*sbe**2*tmp29-
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp130)/mt)
      tmp185 = cbe**2*hb**2*tmp127+ht**2*sbe**2*tmp135-
     .   cbe*hb*ht*sbe*(mb*s2t*tmp128-mt*s2b*tmp130+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp186 = cbe**2*hb**2*tmp126+ht**2*sbe**2*tmp136-
     .   cbe*hb*ht*sbe*(mb*s2t*tmp129+mt*s2b*tmp130+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp187 = ht**2*sbe**2*tmp132+cbe**2*hb**2*tmp20-
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp131)/mt)
      tmp188 = cbe**2*hb**2*tmp117+ht**2*sbe**2*tmp134-
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp131)/mt)
      tmp189 = cbe**2*hb**2*tmp121+ht**2*sbe**2*tmp141-
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp128)-mt*s2b*tmp131+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp190 = cbe**2*hb**2*tmp120+ht**2*sbe**2*tmp142-
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp129)+mt*s2b*tmp131+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp191 = hb**2*sbe**2*tmp145+cbe**2*ht**2*tmp159+
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp156)/mt)
      tmp192 = hb**2*sbe**2*tmp144+cbe**2*ht**2*tmp45+
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp156)/mt)
      tmp193 = hb**2*sbe**2*tmp153+cbe**2*ht**2*tmp161+
     .   cbe*hb*ht*sbe*(mb*s2t*tmp154-mt*s2b*tmp156+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp194 = hb**2*sbe**2*tmp152+cbe**2*ht**2*tmp162+
     .   cbe*hb*ht*sbe*(mb*s2t*tmp155+mt*s2b*tmp156+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp195 = cbe**2*ht**2*tmp158+hb**2*sbe**2*tmp36+
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp157)/mt)
      tmp196 = hb**2*sbe**2*tmp143+cbe**2*ht**2*tmp160+
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp157)/mt)
      tmp197 = hb**2*sbe**2*tmp147+cbe**2*ht**2*tmp167+
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp154)-mt*s2b*tmp157+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp198 = hb**2*sbe**2*tmp146+cbe**2*ht**2*tmp168+
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp155)+mt*s2b*tmp157+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp199 = cbe**2*hb**2*tmp125+ht**2*sbe**2*tmp137-
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp27)+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp128)/s2t+
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp200 = cbe**2*hb**2*tmp123+ht**2*sbe**2*tmp139-
     .   cbe*hb*ht*sbe*(mt*s2b*tmp27+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp129)/s2t-
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp201 = cbe**2*hb**2*tmp124+ht**2*sbe**2*tmp138-
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp28)+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp128)/s2t-
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp202 = cbe**2*hb**2*tmp122+ht**2*sbe**2*tmp140-
     .   cbe*hb*ht*sbe*(mt*s2b*tmp28+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp129)/s2t+
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp203 = hb**2*sbe**2*tmp151+cbe**2*ht**2*tmp163+
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp43)+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp154)/s2t+
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp204 = hb**2*sbe**2*tmp149+cbe**2*ht**2*tmp165+
     .   cbe*hb*ht*sbe*(mt*s2b*tmp43+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp155)/s2t-
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp205 = hb**2*sbe**2*tmp150+cbe**2*ht**2*tmp164+
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp44)+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp154)/s2t-
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp206 = hb**2*sbe**2*tmp148+cbe**2*ht**2*tmp166+
     .   cbe*hb*ht*sbe*(mt*s2b*tmp44+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp155)/s2t+
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp207 = (b**2*(Logb-Logt))/((1d0-b/t)**2*t**4)+
     .   b/((1d0-b/t)*t**3)+(2d0*b*(Logb-Logt))/((1d0-b/t)*t**3)
      tmp208 = (-2d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+
     .   (-2d0-2d0*Logb)/t+(2d0*Logt)/t-
     .   (2d0*(Logb-Logt))/((1d0-b/t)*t)+
     .   (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)**2*t**3)+
     .   (2d0*(-b+t))/((1d0-b/t)*t**2)+
     .   (2d0*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)
      tmp209 = -1d0+2d0*Li2bt+4d0*Logb+(-2d0-2d0*Logb)*Logt+
     .   Logt**2-(2d0*(Logb-Logt)*(-b+t))/((1d0-b/t)*t)
      tmp210 = -5d0-2d0*Li2bt+4d0*Logt-Logt**2+
     .   (2d0*Logt*(b-t))/t+
     .   (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)+
     .   (-2d0*b*Logb+4d0*t)/t
      tmp211 = -1d0+4d0*LogB1+(-2d0-2d0*LogB1)*LogT1+
     .   LogT1**2-(2d0*(LogB1-LogT1)*(-B1+T1))/
     .    ((1d0-B1/T1)*T1)+2d0*tmp169
      tmp212 = -1d0+4d0*LogB2+(-2d0-2d0*LogB2)*LogT1+
     .   LogT1**2-(2d0*(LogB2-LogT1)*(-B2+T1))/
     .    ((1d0-B2/T1)*T1)+2d0*tmp170
      tmp213 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(B1-T1))/T1+
     .   (2d0*B1*(LogB1-LogT1)*(-B1+T1))/((1d0-B1/T1)*T1**2)+
     .   (-2d0*B1*LogB1+4d0*T1)/T1-2d0*tmp169
      tmp214 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(B2-T1))/T1+
     .   (2d0*B2*(LogB2-LogT1)*(-B2+T1))/((1d0-B2/T1)*T1**2)+
     .   (-2d0*B2*LogB2+4d0*T1)/T1-2d0*tmp170
      tmp215 = -1d0+4d0*LogB1+(-2d0-2d0*LogB1)*LogT2+
     .   LogT2**2-(2d0*(LogB1-LogT2)*(-B1+T2))/
     .    ((1d0-B1/T2)*T2)+2d0*tmp172
      tmp216 = -1d0+4d0*LogB2+(-2d0-2d0*LogB2)*LogT2+
     .   LogT2**2-(2d0*(LogB2-LogT2)*(-B2+T2))/
     .    ((1d0-B2/T2)*T2)+2d0*tmp173
      tmp217 = -5d0+4d0*LogT2-LogT2**2+
     .   (2d0*LogT2*(B1-T2))/T2+
     .   (2d0*B1*(LogB1-LogT2)*(-B1+T2))/((1d0-B1/T2)*T2**2)+
     .   (-2d0*B1*LogB1+4d0*T2)/T2-2d0*tmp172
      tmp218 = -5d0+4d0*LogT2-LogT2**2+
     .   (2d0*LogT2*(B2-T2))/T2+
     .   (2d0*B2*(LogB2-LogT2)*(-B2+T2))/((1d0-B2/T2)*T2**2)+
     .   (-2d0*B2*LogB2+4d0*T2)/T2-2d0*tmp173
      tmp219 = -1d0+LogT1**2+LogT1*(-2d0-2d0*LogT2)+
     .   4d0*LogT2-(2d0*(-LogT1+LogT2)*(T1-T2))/
     .    (T1*(1d0-T2/T1))+2d0*tmp171
      tmp220 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(-T1+T2))/T1+(4d0*T1-2d0*LogT2*T2)/T1+
     .   (2d0*(-LogT1+LogT2)*(T1-T2)*T2)/(T1**2*(1d0-T2/T1))-
     .   2d0*tmp171
      tmp221 = -0.5d0+2d0*LogT1-
     .   (phiA0T1T2*(-A0+T1-T2))/T2-0.5d0*(LogA0*LogT1)+
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogT1*LogT2)+
     .   0.5d0*(LogT2*(A0-T1-T2))/T1+
     .   0.5d0*(LogA0*(-A0-T1+T2))/T1-
     .   0.5d0*(deltA0T1T2*((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .       (T2*tmp96)/(deltA0T1T2*T1)))/T2
      tmp222 = -0.5d0+2d0*LogT2-
     .   (phiA0B1T2*(-A0-B1+T2))/T2+0.5d0*(LogA0*LogB1)-
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogB1*LogT2)+
     .   0.5d0*(LogB1*(A0-B1-T2))/T2+
     .   0.5d0*(LogA0*(-A0+B1-T2))/T2-
     .   0.5d0*(deltA0B1T2*((B1*phiA0B1T2*(A0+B1-T2))/
     .        (deltA0B1T2*T2)+(B1*tmp100)/(deltA0B1T2*T2)))/B1
      tmp223 = -0.5d0+2d0*LogT2-
     .   (phiA0B2T2*(-A0-B2+T2))/T2+0.5d0*(LogA0*LogB2)-
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogB2*LogT2)+
     .   0.5d0*(LogB2*(A0-B2-T2))/T2+
     .   0.5d0*(LogA0*(-A0+B2-T2))/T2-
     .   0.5d0*(deltA0B2T2*((B2*phiA0B2T2*(A0+B2-T2))/
     .        (deltA0B2T2*T2)+(B2*tmp101)/(deltA0B2T2*T2)))/B2
      tmp224 = -0.5d0+2d0*LogT2-
     .   (phiT2bmu2*(-b-mu2+T2))/mu2+0.5d0*(Logb*Logmu2)-
     .   0.5d0*(Logb*LogT2)-0.5d0*(Logmu2*LogT2)+
     .   0.5d0*(Logmu2*(b-mu2-T2))/T2+
     .   0.5d0*(Logb*(-b+mu2-T2))/T2-
     .   0.5d0*(delT2Tbmu2*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .       (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2
      tmp225 = 0.5d0-2d0*LogT2+
     .   (phiT2bmu2*(-b-mu2+T2))/mu2-0.5d0*(Logb*Logmu2)+
     .   0.5d0*(Logb*LogT2)+0.5d0*(Logmu2*LogT2)-
     .   0.5d0*(Logmu2*(b-mu2-T2))/T2-
     .   0.5d0*(Logb*(-b+mu2-T2))/T2+
     .   0.5d0*(delT2Tbmu2*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .       (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2
      tmp226 = -(phiT2bmu2/delT2Tbmu2)-
     .   (2d0*phiT2bmu2*(b+mu2-T2)*(-b-mu2+T2))/
     .    delT2Tbmu2**2-(mu2*tmp102)/(delT2Tbmu2*T2**2)-
     .   (2d0*mu2*(-b-mu2+T2)*tmp102)/(delT2Tbmu2**2*T2)+
     .   ((b+mu2-T2)*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .      (mu2*tmp102)/(delT2Tbmu2*T2)))/delT2Tbmu2+
     .   (mu2*tmp5)/(delT2Tbmu2*T2)
      tmp227 = phiT2bmu2/delT2Tbmu2-
     .   (2d0*phiT2bmu2*(b-mu2-T2)*(b+mu2-T2))/
     .    delT2Tbmu2**2-(2d0*mu2*(b-mu2-T2)*tmp102)/
     .    (delT2Tbmu2**2*T2)+
     .   ((b+mu2-T2)*((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .      (mu2*tmp110)/(b*delT2Tbmu2)))/delT2Tbmu2+
     .   (mu2*tmp97)/(delT2Tbmu2*T2)
      tmp228 = -0.5d0+2d0*Logt-
     .   (phimu2T2T*(-mu2+t-T2))/T2-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*LogT2)-0.5d0*(Logt*LogT2)+
     .   0.5d0*(LogT2*(mu2-t-T2))/t+
     .   0.5d0*(Logmu2*(-mu2-t+T2))/t-
     .   0.5d0*(deltmu2T2T*((mu2*phimu2T2T*(mu2-t+T2))/
     .        (deltmu2T2T*T2)+(mu2*tmp111)/(deltmu2T2T*t)))/mu2
      tmp229 = -0.5d0+2d0*Logt-(phiA0bt*(-A0-b+t))/t+
     .   0.5d0*(LogA0*Logb)-0.5d0*(LogA0*Logt)-0.5d0*(Logb*Logt)+
     .   0.5d0*(Logb*(A0-b-t))/t+0.5d0*(LogA0*(-A0+b-t))/t-
     .   0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .       (b*tmp66)/(deltA0bt*t)))/b
      tmp230 = -((b*phiA0bt)/(deltA0bt*t))-
     .   (2d0*b*phiA0bt*(A0+b-t)*(-A0-b+t))/
     .    (deltA0bt**2*t)+(b*tmp1)/(deltA0bt*t)-
     .   (b*tmp66)/(deltA0bt*t**2)-
     .   (2d0*b*(-A0-b+t)*tmp66)/(deltA0bt**2*t)+
     .   ((A0+b-t)*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .      (b*tmp66)/(deltA0bt*t)))/deltA0bt
      tmp231 = -0.5d0+2d0*Logt-
     .   (phiB1tmu2*(-B1-mu2+t))/mu2+0.5d0*(LogB1*Logmu2)-
     .   0.5d0*(LogB1*Logt)-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*(B1-mu2-t))/t+
     .   0.5d0*(LogB1*(-B1+mu2-t))/t-
     .   0.5d0*(deltB1tmu2*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .       (mu2*tmp67)/(deltB1tmu2*t)))/mu2
      tmp232 = -(phiB1tmu2/deltB1tmu2)-
     .   (2d0*phiB1tmu2*(B1+mu2-t)*(-B1-mu2+t))/
     .    deltB1tmu2**2+(mu2*tmp2)/(deltB1tmu2*t)-
     .   (mu2*tmp67)/(deltB1tmu2*t**2)-
     .   (2d0*mu2*(-B1-mu2+t)*tmp67)/(deltB1tmu2**2*t)+
     .   ((B1+mu2-t)*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .      (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2
      tmp233 = phiB1tmu2/deltB1tmu2-
     .   (2d0*phiB1tmu2*(-B1-mu2+t)*(-B1+mu2+t))/
     .    deltB1tmu2**2+((-B1+mu2+t)*
     .      ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .      (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2-
     .   (2d0*mu2*(-B1-mu2+t)*tmp71)/(B1*deltB1tmu2**2)+
     .   (mu2*tmp73)/(B1*deltB1tmu2)
      tmp234 = -0.5d0+2d0*Logt-
     .   (phiB2tmu2*(-B2-mu2+t))/mu2+0.5d0*(LogB2*Logmu2)-
     .   0.5d0*(LogB2*Logt)-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*(B2-mu2-t))/t+
     .   0.5d0*(LogB2*(-B2+mu2-t))/t-
     .   0.5d0*(deltB2tmu2*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .       (mu2*tmp68)/(deltB2tmu2*t)))/mu2
      tmp235 = -((phiB1tmu2*(-B1-mu2+t))/mu2)+
     .   (phiB2tmu2*(-B2-mu2+t))/mu2+0.5d0*(LogB1*Logmu2)-
     .   0.5d0*(LogB2*Logmu2)-0.5d0*(LogB1*Logt)+
     .   0.5d0*(LogB2*Logt)+0.5d0*(Logmu2*(B1-mu2-t))/t-
     .   0.5d0*(Logmu2*(B2-mu2-t))/t+
     .   0.5d0*(LogB1*(-B1+mu2-t))/t-
     .   0.5d0*(LogB2*(-B2+mu2-t))/t-
     .   0.5d0*(deltB1tmu2*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .       (mu2*tmp67)/(deltB1tmu2*t)))/mu2+
     .   0.5d0*(deltB2tmu2*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .       (mu2*tmp68)/(deltB2tmu2*t)))/mu2
      tmp236 = -(phiB2tmu2/deltB2tmu2)-
     .   (2d0*phiB2tmu2*(B2+mu2-t)*(-B2-mu2+t))/
     .    deltB2tmu2**2+(mu2*tmp3)/(deltB2tmu2*t)-
     .   (mu2*tmp68)/(deltB2tmu2*t**2)-
     .   (2d0*mu2*(-B2-mu2+t)*tmp68)/(deltB2tmu2**2*t)+
     .   ((B2+mu2-t)*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .      (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2
      tmp237 = phiB2tmu2/deltB2tmu2-
     .   (2d0*phiB2tmu2*(-B2-mu2+t)*(-B2+mu2+t))/
     .    deltB2tmu2**2+((-B2+mu2+t)*
     .      ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .      (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2-
     .   (2d0*mu2*(-B2-mu2+t)*tmp72)/(B2*deltB2tmu2**2)+
     .   (mu2*tmp74)/(B2*deltB2tmu2)
      tmp238 = -0.5d0+2d0*LogT1-
     .   (phiA0B1T1*(-A0-B1+T1))/T1+0.5d0*(LogA0*LogB1)-
     .   0.5d0*(LogA0*LogT1)-0.5d0*(LogB1*LogT1)+
     .   0.5d0*(LogB1*(A0-B1-T1))/T1+
     .   0.5d0*(LogA0*(-A0+B1-T1))/T1-
     .   0.5d0*(deltA0B1T1*((B1*phiA0B1T1*(A0+B1-T1))/
     .        (deltA0B1T1*T1)+(B1*tmp81)/(deltA0B1T1*T1)))/B1
      tmp239 = -0.5d0+2d0*LogT1-
     .   (phiA0B2T1*(-A0-B2+T1))/T1+0.5d0*(LogA0*LogB2)-
     .   0.5d0*(LogA0*LogT1)-0.5d0*(LogB2*LogT1)+
     .   0.5d0*(LogB2*(A0-B2-T1))/T1+
     .   0.5d0*(LogA0*(-A0+B2-T1))/T1-
     .   0.5d0*(deltA0B2T1*((B2*phiA0B2T1*(A0+B2-T1))/
     .        (deltA0B2T1*T1)+(B2*tmp82)/(deltA0B2T1*T1)))/B2
      tmp240 = -0.5d0+2d0*LogT1-
     .   (phiT1bmu2*(-b-mu2+T1))/mu2+0.5d0*(Logb*Logmu2)-
     .   0.5d0*(Logb*LogT1)-0.5d0*(Logmu2*LogT1)+
     .   0.5d0*(Logmu2*(b-mu2-T1))/T1+
     .   0.5d0*(Logb*(-b+mu2-T1))/T1-
     .   0.5d0*(deltT1bmu2*((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .       (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
      tmp241 = -(phiT1bmu2/deltT1bmu2)-
     .   (2d0*phiT1bmu2*(b+mu2-T1)*(-b-mu2+T1))/
     .    deltT1bmu2**2+(mu2*tmp4)/(deltT1bmu2*T1)-
     .   (mu2*tmp83)/(deltT1bmu2*T1**2)-
     .   (2d0*mu2*(-b-mu2+T1)*tmp83)/(deltT1bmu2**2*T1)+
     .   ((b+mu2-T1)*((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .      (mu2*tmp83)/(deltT1bmu2*T1)))/deltT1bmu2
      tmp242 = phiT1bmu2/deltT1bmu2-
     .   (2d0*phiT1bmu2*(b-mu2-T1)*(b+mu2-T1))/
     .    deltT1bmu2**2+(mu2*tmp78)/(deltT1bmu2*T1)-
     .   (2d0*mu2*(b-mu2-T1)*tmp83)/(deltT1bmu2**2*T1)+
     .   ((b+mu2-T1)*((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .      (mu2*tmp90)/(b*deltT1bmu2)))/deltT1bmu2
      tmp243 = -0.5d0+2d0*Logt-
     .   (phimu2tT1*(-mu2+t-T1))/T1-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*LogT1)-0.5d0*(Logt*LogT1)+
     .   0.5d0*(LogT1*(mu2-t-T1))/t+
     .   0.5d0*(Logmu2*(-mu2-t+T1))/t-
     .   0.5d0*(deltmu2tT1*((mu2*phimu2tT1*(mu2-t+T1))/
     .        (deltmu2tT1*T1)+(mu2*tmp91)/(deltmu2tT1*t)))/mu2
      tmp244 = 2d0/T1-Logb/T1-Logmu2/T1-
     .   0.5d0*(Logmu2*(b-mu2-T1))/T1**2-
     .   0.5d0*(Logb*(-b+mu2-T1))/T1**2-
     .   0.5d0*(2d0*phiT1bmu2+deltT1bmu2*tmp241+
     .       4d0*(-b-mu2+T1)*
     .      ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .        (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
      D1t = -(tmp109*tmp195)-tmp108*tmp196-
     .   2d0*hb*ht*mt*mu*s2b*tmp235-(hb*ht*mu*s2b*tmp75)/mt+
     .   tmp61*(-(B1*(-1d0+LogB1))+2d0*B1*LogB1-
     .      B1*(-1d0+LogB1)*(-1d0+Logt)+(-1d0+Logmu2)*mu2+
     .      2d0*Logmu2*mu2+(-1d0+Logmu2)*(-1d0+Logt)*mu2+
     .      2d0*Logt*t-(B1-mu2-t)*tmp231-
     .      0.5d0*(deltB1tmu2*phiB1tmu2)/mu2+
     .      0.5d0*(Logmu2*Logt*(B1-mu2-t))+
     .      0.5d0*(LogB1*Logt*(-B1+mu2-t))+
     .      0.5d0*(LogB1*Logmu2*(-B1-mu2+t))-2.5d0*(B1+mu2+t)
     .      )+tmp60*(-(B2*(-1d0+LogB2))+2d0*B2*LogB2-
     .      B2*(-1d0+LogB2)*(-1d0+Logt)+(-1d0+Logmu2)*mu2+
     .      2d0*Logmu2*mu2+(-1d0+Logmu2)*(-1d0+Logt)*mu2+
     .      2d0*Logt*t-(B2-mu2-t)*tmp234-
     .      0.5d0*(deltB2tmu2*phiB2tmu2)/mu2+
     .      0.5d0*(Logmu2*Logt*(B2-mu2-t))+
     .      0.5d0*(LogB2*Logt*(-B2+mu2-t))+
     .      0.5d0*(LogB2*Logmu2*(-B2-mu2+t))-2.5d0*(B2+mu2+t)
     .      )-0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-
     .      5d0*T1+LogT1*(-2d0*B2*LogB2+4d0*T1)-
     .      2d0*(-B2+T1)*tmp170)*tmp183)-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      tmp184)-0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .      LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      tmp187)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      tmp188)-4d0*cbe*hb*ht*mb*mt*sbe*
     .    (0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .      0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .      0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b)-
     .   (2d0*cbe*hb*ht*mb*sbe*
     .      (-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76))/mt+
     .   tmp10*(b*(-1d0+Logb)+b*(-1d0+Logb)*(-1d0+Logt)+
     .      0.5d0*((b+t)*tmp210)+0.5d0*tmp76)
      D1t = D1t-tmp192*tmp88-tmp191*tmp89+
     .   tmp9*(-(A0*(-1d0+LogA0))+2d0*A0*LogA0+b*(-1d0+Logb)+
     .      2d0*b*Logb-A0*(-1d0+LogA0)*(-1d0+Logt)+
     .      b*(-1d0+Logb)*(-1d0+Logt)+2d0*Logt*t-
     .      (A0-b-t)*tmp229+0.5d0*(Logb*Logt*(A0-b-t))+
     .      0.5d0*(LogA0*Logt*(-A0+b-t))-
     .      0.5d0*(deltA0bt*phiA0bt)/t+
     .      0.5d0*(LogA0*Logb*(-A0-b+t))-2.5d0*(A0+b+t))+
     .   ht**2*(2d0*(-1d0+Logmu2)*mu2+4d0*Logmu2*mu2+
     .      2d0*(-1d0+Logmu2)*(-1d0+Logt)*mu2+4d0*Logt*t-
     .      (-1d0+LogT1)*T1-(-1d0+Logt)*(-1d0+LogT1)*T1+
     .      2d0*LogT1*T1-(-1d0+LogT2)*T2-
     .      (-1d0+Logt)*(-1d0+LogT2)*T2+2d0*LogT2*T2-
     .      (-mu2-t+T2)*tmp228-(-mu2-t+T1)*tmp243+
     .      sbe**2*(-10d0*t+2d0*(-1d0+Logt)*t+2d0*(-1d0+Logt)**2*t+
     .       4d0*Logt*t+Logt*(4d0*t-2d0*Logt*t)+t*tmp77)+
     .      0.5d0*(Logt*LogT1*(mu2-t-T1))+
     .      0.5d0*(Logmu2*LogT1*(-mu2+t-T1))-
     .      0.5d0*(deltmu2tT1*phimu2tT1)/T1+
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T1))-
     .      2.5d0*(mu2+t+T1)+0.5d0*(Logt*LogT2*(mu2-t-T2))+
     .      0.5d0*(Logmu2*LogT2*(-mu2+t-T2))-
     .      0.5d0*(deltmu2T2T*phimu2T2T)/T2+
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T2))-
     .      2.5d0*(mu2+t+T2)+
     .      cbe**2*(-2d0*A0*(-1d0+LogA0)-
     .       2d0*A0*(-1d0+LogA0)*(-1d0+Logt)+2d0*(-1d0+Logt)*t+
     .       2d0*(-1d0+Logt)**2*t+
     .       2d0*(2d0*A0*LogA0-A0*LogA0*Logt+4d0*Logt*t+
     .          0.5d0*(Logt**2*(A0-2d0*t))-
     .          0.5d0*(deltA0tt*phiA0tt)/t-2.5d0*(A0+2d0*t))-
     .       (A0-2d0*t)*(-1d0+4d0*Logt-Logt**2-(A0*LogA0)/t+
     .          (2d0*A0*phiA0tt)/t+(Logt*(A0-2d0*t))/t+
     .          0.5d0*(deltA0tt*phiA0tt)/t**2-
     .          0.5d0*(deltA0tt*
     .            ((A0*phiA0tt)/deltA0tt+
     .              (phiA0tt*tmp175)/deltA0tt+
     .              tmp65/deltA0tt+tmp70/deltA0tt))/t)))-
     .   0.25d0*(ht**2*(cbe**2*tmp114*(4d0-(2d0*s2t*Yt)/mt)+
     .      cbe**2*tmp93*(4d0+(2d0*s2t*Yt)/mt)+
     .      0.5d0*(sbe**2*tmp115*(4d0-(2d0*s2t*Xt)/mt))+
     .      0.5d0*(sbe**2*tmp94*(4d0+(2d0*s2t*Xt)/mt))))
      DT1 = -(tmp194*tmp238)-tmp193*tmp239-
     .   2d0*hb*ht*mb*mu*s2t*tmp240+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      sbe**2*(0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0))+
     .       0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0)*(-1d0+LogT1)))+
     .      0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1))
     .      )+tmp62*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT1)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT1)*mu2-
     .      2d0*LogT1*T1-(-b-mu2+T1)*tmp240+
     .      0.5d0*(deltT1bmu2*phiT1bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))-
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))+2.5d0*(b+mu2+T1))
     .    -0.5d0*(tmp186*tmp213)-0.5d0*(tmp185*tmp214)+
     .   ht**2*((-1d0+LogT2)*T2*tmp8+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T2*tmp8+
     .      cbe**2*(A0*(-1d0+LogA0)*(1d0+0.5d0*(1d0-c2t))+
     .       A0*(-1d0+LogA0)*(-1d0+LogT1)*(1d0+0.5d0*(1d0-c2t)))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2)*
     .       (-1d0+LogT1))+0.25d0*((1d0+Nc)*s2t**2*tmp79))+
     .   ht**2*(-((-1d0+Logmu2)*mu2)-2d0*Logmu2*mu2-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2-(-1d0+Logt)*t-
     .      2d0*Logt*t-(-1d0+Logt)*(-1d0+LogT1)*t-2d0*LogT1*T1-
     .      0.5d0*(Logt*LogT1*(mu2-t-T1))-
     .      0.5d0*(Logmu2*LogT1*(-mu2+t-T1))+
     .      0.5d0*(deltmu2tT1*phimu2tT1)/T1-
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T1))+
     .      2.5d0*(mu2+t+T1)-
     .      (-mu2-t+T1)*
     .       (-0.5d0+2d0*LogT1-(phimu2tT1*(-mu2-t+T1))/T1+
     .       0.5d0*(Logmu2*Logt)-0.5d0*(Logmu2*LogT1)-
     .       0.5d0*(Logt*LogT1)+0.5d0*(Logt*(mu2-t-T1))/T1+
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1-
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1*(mu2+t-T1))/
     .            (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1)
     .             ))/mu2))
      DT1 = DT1-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*tmp220*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*tmp221*Yt**2+
     .      cbe**2*tmp51*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp35*tmp95)))
      DT2 = -(tmp198*tmp222)-tmp197*tmp223-
     .   2d0*hb*ht*mb*mu*s2t*tmp225+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      sbe**2*(0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0))+
     .       0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0)*(-1d0+LogT2)))+
     .      0.25d0*(B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT2))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2))
     .      )+tmp63*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT2)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT2)*mu2-
     .      2d0*LogT2*T2-(-b-mu2+T2)*tmp224+
     .      0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2))
     .    +ht**2*(-((-1d0+Logmu2)*mu2)-2d0*Logmu2*mu2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2-(-1d0+Logt)*t-
     .      2d0*Logt*t-(-1d0+Logt)*(-1d0+LogT2)*t-2d0*LogT2*T2-
     .      0.5d0*(Logt*LogT2*(mu2-t-T2))-
     .      0.5d0*(Logmu2*LogT2*(-mu2+t-T2))+
     .      0.5d0*(deltmu2T2T*phimu2T2T)/T2-
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T2))+
     .      2.5d0*(mu2+t+T2)-
     .      (-mu2-t+T2)*
     .       (-0.5d0+2d0*LogT2-(phimu2T2T*(-mu2-t+T2))/T2+
     .       0.5d0*(Logmu2*Logt)-0.5d0*(Logmu2*LogT2)-
     .       0.5d0*(Logt*LogT2)+0.5d0*(Logt*(mu2-t-T2))/T2+
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2-
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T*(mu2+t-T2))/
     .            (deltmu2T2T*T2)+
     .             (mu2*tmp103)/(deltmu2T2T*T2)))/mu2))-
     .   0.5d0*(tmp190*tmp217)-0.5d0*(tmp189*tmp218)+
     .   ht**2*((-1d0+LogT1)*T1*tmp8+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T1*tmp8+
     .      cbe**2*(A0*(-1d0+LogA0)*(1d0+0.5d0*(1d0+c2t))+
     .       A0*(-1d0+LogA0)*(-1d0+LogT2)*(1d0+0.5d0*(1d0+c2t)))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT2))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2)*
     .       (-1d0+LogT2))+0.25d0*((1d0+Nc)*s2t**2*tmp98))
      DT2 = DT2-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*tmp219*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (-0.5d0+2d0*LogT2-(phiA0T1T2*(-A0-T1+T2))/T2+
     .         0.5d0*(LogA0*LogT1)-0.5d0*(LogA0*LogT2)-
     .         0.5d0*(LogT1*LogT2)+
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2**2+
     .         0.5d0*(LogT1*(A0-T1-T2))/T2+
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2-
     .         0.5d0*(deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*tmp34)+
     .      cbe**2*tmp50*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)))
      Dc2t = (hb*ht*mb*mu*tmp113)/s2t-tmp109*tmp205-
     .   tmp108*tmp206+(b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT2)*T2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2*T2-
     .      (-b-mu2+T2)*tmp112)*tmp7-tmp204*tmp88-
     .   tmp203*tmp89+tmp6*
     .    (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT1)*T1-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2*T1-
     .      (-b-mu2+T1)*tmp92)+
     .   hb**2*(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t+0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*
     .        T1)/c2t+sbe**2*
     .       (0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t)-
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t-0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*
     .        T2)/c2t)+
     .   ht**2*((-1d0+LogT1)*(-1d0+LogT2)*T1*T2*
     .       (1d0+0.5d0*(-1d0+Nc))-
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t-0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*
     .        T1)/c2t+cbe**2*
     .       (-(0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t)+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t)+
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t+0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*
     .        T2)/c2t-0.25d0*
     .       ((1d0+Nc)*((-1d0+LogT1)**2*T1**2+
     .         (-1d0+LogT2)**2*T2**2)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      tmp199)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      tmp200)-0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .      LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      tmp201)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      tmp202)
      Dc2t = Dc2t-0.25d0*
     .    (ht**2*(cbe**2*tmp114*tmp47+cbe**2*tmp46*tmp93+
     .      sbe**2*(-5d0*T1-5d0*T2+4d0*LogT2*T2+
     .         LogT1**2*(-T1+T2)+LogT1*(4d0*T1-2d0*LogT2*T2)-
     .         2d0*(T1-T2)*tmp171)*Xt**2+
     .      2d0*cbe**2*Yt**2*
     .       (2d0*A0*LogA0+2d0*LogT1*T1+2d0*LogT2*T2+
     .         0.5d0*(LogT1*LogT2*(A0-T1-T2))+
     .         0.5d0*(LogA0*LogT2*(-A0+T1-T2))-
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2+
     .         0.5d0*(LogA0*LogT1*(-A0-T1+T2))-
     .         2.5d0*(A0+T1+T2))+0.5d0*(sbe**2*tmp115*tmp31)+
     .      0.5d0*(sbe**2*tmp30*tmp94)))
      DT1T1 = -2d0*hb*ht*mb*mu*s2t*tmp244+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))/T1+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))/T1+
     .      0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0)*sbe**2)/T1)+
     .   ht**2*(((-1d0+LogT2)*T2*tmp8)/T1+
     .      (A0*cbe**2*(-1d0+LogA0)*(1d0+0.5d0*(1d0-c2t)))/T1+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))/T1+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))/T1+
     .      0.25d0*((1d0+Nc)*s2t**2*
     .       (8d0*(-1d0+LogT1)+2d0*(-1d0+LogT1)**2+
     .         (2d0/T1**2-(2d0*(-1d0+LogT1))/T1**2)*T1**2)))-
     .   0.5d0*(((4d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**2)+8d0/T1-
     .      (4d0*LogT1)/T1-
     .      2d0*((B2**2*(LogB2-LogT1))/((1d0-B2/T1)**2*T1**4)+
     .         B2/((1d0-B2/T1)*T1**3)+
     .         (2d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**3))*
     .       (-B2+T1)-(-2d0*B2*LogB2+4d0*T1)/T1**2+
     .      (B2-T1)*tmp17)*tmp185)-
     .   0.5d0*(((4d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**2)+8d0/T1-
     .      (4d0*LogT1)/T1-
     .      2d0*((B1**2*(LogB1-LogT1))/((1d0-B1/T1)**2*T1**4)+
     .         B1/((1d0-B1/T1)*T1**3)+
     .         (2d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**3))*
     .       (-B1+T1)-(-2d0*B1*LogB1+4d0*T1)/T1**2+
     .      (B1-T1)*tmp17)*tmp186)-
     .   tmp194*(2d0/T1-LogA0/T1-LogB1/T1-
     .      0.5d0*(LogB1*(A0-B1-T1))/T1**2-
     .      0.5d0*(LogA0*(-A0+B1-T1))/T1**2-
     .      0.5d0*((2d0*B1*phiA0B1T1)/T1+
     .        4d0*(-A0-B1+T1)*
     .         ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .           (B1*tmp81)/(deltA0B1T1*T1))+
     .        deltA0B1T1*
     .         ((B1*(2d0-LogA0-LogB1+2d0*LogT1))/
     .            (deltA0B1T1*T1)-
     .           (B1*phiA0B1T1)/(deltA0B1T1*T1)-
     .           (2d0*B1*phiA0B1T1*(A0+B1-T1)*(-A0-B1+T1))/
     .            (deltA0B1T1**2*T1)-
     .           (B1*tmp81)/(deltA0B1T1*T1**2)-
     .           (2d0*B1*(-A0-B1+T1)*tmp81)/
     .            (deltA0B1T1**2*T1)+
     .           ((A0+B1-T1)*
     .            ((B1*phiA0B1T1*(A0+B1-T1))/
     .               (deltA0B1T1*T1)+
     .              (B1*tmp81)/(deltA0B1T1*T1)))/deltA0B1T1))/
     .      B1)
      DT1T1 = DT1T1-tmp193*
     .    (2d0/T1-LogA0/T1-LogB2/T1-
     .      0.5d0*(LogB2*(A0-B2-T1))/T1**2-
     .      0.5d0*(LogA0*(-A0+B2-T1))/T1**2-
     .      0.5d0*((2d0*B2*phiA0B2T1)/T1+
     .        4d0*(-A0-B2+T1)*
     .         ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .           (B2*tmp82)/(deltA0B2T1*T1))+
     .        deltA0B2T1*
     .         ((B2*(2d0-LogA0-LogB2+2d0*LogT1))/
     .            (deltA0B2T1*T1)-
     .           (B2*phiA0B2T1)/(deltA0B2T1*T1)-
     .           (2d0*B2*phiA0B2T1*(A0+B2-T1)*(-A0-B2+T1))/
     .            (deltA0B2T1**2*T1)-
     .           (B2*tmp82)/(deltA0B2T1*T1**2)-
     .           (2d0*B2*(-A0-B2+T1)*tmp82)/
     .            (deltA0B2T1**2*T1)+
     .           ((A0+B2-T1)*
     .            ((B2*phiA0B2T1*(A0+B2-T1))/
     .               (deltA0B2T1*T1)+
     .              (B2*tmp82)/(deltA0B2T1*T1)))/deltA0B2T1))/
     .      B2)+tmp62*(-((b*(-1d0+Logb))/T1)-
     .      ((-1d0+Logmu2)*mu2)/T1-(-b-mu2+T1)*tmp244+
     .      2d0*(0.5d0-2d0*LogT1+(phiT1bmu2*(-b-mu2+T1))/mu2-
     .       0.5d0*(Logb*Logmu2)+0.5d0*(Logb*LogT1)+
     .       0.5d0*(Logmu2*LogT1)-
     .       0.5d0*(Logmu2*(b-mu2-T1))/T1-
     .       0.5d0*(Logb*(-b+mu2-T1))/T1+
     .       0.5d0*(deltT1bmu2*
     .           ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .             (mu2*tmp83)/(deltT1bmu2*T1)))/mu2))
      DT1T1 = DT1T1+ht**2*
     .    (-(((-1d0+Logmu2)*mu2)/T1)-((-1d0+Logt)*t)/T1+
     .      2d0*(0.5d0-2d0*LogT1+(phimu2tT1*(-mu2-t+T1))/T1-
     .       0.5d0*(Logmu2*Logt)+0.5d0*(Logmu2*LogT1)+
     .       0.5d0*(Logt*LogT1)-0.5d0*(Logt*(mu2-t-T1))/T1-
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1+
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1*(mu2+t-T1))/
     .            (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1)
     .             ))/mu2)-
     .      (-mu2-t+T1)*
     .       (2d0/T1-Logmu2/T1-Logt/T1-
     .       0.5d0*(Logt*(mu2-t-T1))/T1**2-
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1**2-
     .       0.5d0*((2d0*mu2*phimu2tT1)/T1+
     .           4d0*(-mu2-t+T1)*
     .            ((mu2*phimu2tT1*(mu2+t-T1))/
     .             (deltmu2tT1*T1)+
     .            (mu2*tmp84)/(deltmu2tT1*T1))+
     .           deltmu2tT1*
     .            (((2d0-Logmu2-Logt+2d0*LogT1)*mu2)/
     .             (deltmu2tT1*T1)-
     .            (mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .            (2d0*mu2*phimu2tT1*(mu2+t-T1)*
     .               (-mu2-t+T1))/(deltmu2tT1**2*T1)-
     .            (mu2*tmp84)/(deltmu2tT1*T1**2)-
     .            (2d0*mu2*(-mu2-t+T1)*tmp84)/
     .             (deltmu2tT1**2*T1)+
     .            ((mu2+t-T1)*
     .               ((mu2*phimu2tT1*(mu2+t-T1))/
     .                  (deltmu2tT1*T1)+
     .                 (mu2*tmp84)/(deltmu2tT1*T1)))/deltmu2tT1
     .            ))/mu2))
      DT1T1 = DT1T1-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*
     .       (8d0/T1-(4d0*LogT1)/T1-(4d0*T1-2d0*LogT2*T2)/T1**2+
     .         (4d0*(-LogT1+LogT2)*T2)/(T1**2*(1d0-T2/T1))-
     .         2d0*(T1-T2)*
     .          (((-LogT1+LogT2)*T2**2)/
     .             (T1**4*(1d0-T2/T1)**2)+
     .            T2/(T1**3*(1d0-T2/T1))+
     .            (2d0*(-LogT1+LogT2)*T2)/(T1**3*(1d0-T2/T1)))+
     .         (-T1+T2)*tmp17)*Xt**2+
     .      0.5d0*(sbe**2*(4d0/T1+(2d0*(2d0-2d0*LogT1))/T1-
     .           (2d0*LogT1)/T1-(4d0*T1-2d0*LogT1*T1)/T1**2)*tmp35)
     .       +cbe**2*tmp51*
     .       (-((deltA0T1T1*phiA0T1T1)/T1**3)+
     .         (A0*LogA0)/T1**2+4d0/T1-(4d0*LogT1)/T1+
     .         (-4d0*A0*phiA0T1T1+
     .            deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1**2+
     .           0.5d0*((A0-2d0*T1)*tmp17)-
     .         0.5d0*(-8d0*A0*
     .            ((A0*phiA0T1T1)/deltA0T1T1+
     .              (phiA0T1T1*tmp178)/deltA0T1T1+
     .              tmp80/deltA0T1T1+tmp87/deltA0T1T1)+
     .             deltA0T1T1*
     .            ((4d0*A0**2*phiA0T1T1)/deltA0T1T1**2+
     .              (phiA0T1T1*
     .                 (-1d0-(A0-T1)**2/T1**2-
     .                   (2d0*(A0-T1))/T1))/deltA0T1T1+
     .              (1d0+(A0-T1)/T1)/deltA0T1T1+
     .              (1d0-(-A0+T1)/T1)/deltA0T1T1+
     .              (4d0*A0*phiA0T1T1*tmp178)/deltA0T1T1**2+
     .              (4d0*A0*tmp80)/deltA0T1T1**2+
     .              (4d0*A0*tmp87)/deltA0T1T1**2+
     .              (A0*
     .                 ((A0*phiA0T1T1)/deltA0T1T1+
     .                   (phiA0T1T1*tmp178)/deltA0T1T1+
     .                   tmp80/deltA0T1T1+tmp87/deltA0T1T1))/
     .               deltA0T1T1+
     .              (tmp178*
     .                 ((A0*phiA0T1T1)/deltA0T1T1+
     .                   (phiA0T1T1*tmp178)/deltA0T1T1+
     .                   tmp80/deltA0T1T1+tmp87/deltA0T1T1))/
     .               deltA0T1T1))/T1)+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (2d0/T1-LogA0/T1-LogT2/T1-
     .         0.5d0*(LogT2*(A0-T1-T2))/T1**2-
     .         0.5d0*(LogA0*(-A0-T1+T2))/T1**2-
     .         0.5d0*(2d0*phiA0T1T2+
     .             4d0*(-A0+T1-T2)*
     .            ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .              (T2*tmp96)/(deltA0T1T2*T1))+
     .             deltA0T1T2*
     .            (-(phiA0T1T2/deltA0T1T2)+
     .              ((2d0-LogA0+2d0*LogT1-LogT2)*T2)/
     .               (deltA0T1T2*T1)-
     .              (2d0*phiA0T1T2*(-A0+T1-T2)*
     .                 (A0-T1+T2))/deltA0T1T2**2-
     .              (T2*tmp96)/(deltA0T1T2*T1**2)-
     .              (2d0*(-A0+T1-T2)*T2*tmp96)/
     .               (deltA0T1T2**2*T1)+
     .              ((A0-T1+T2)*
     .                 ((phiA0T1T2*(A0-T1+T2))/
     .                  deltA0T1T2+(T2*tmp96)/(deltA0T1T2*T1)
     .                   ))/deltA0T1T2))/T2)))
      DT2T2 = hb**2*(0.25d0*
     .       (B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))/T2+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))/T2+
     .      0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0)*sbe**2)/T2)+
     .   ht**2*(-(((-1d0+Logmu2)*mu2)/T2)-((-1d0+Logt)*t)/T2+
     .      2d0*(0.5d0-2d0*LogT2+(phimu2T2T*(-mu2-t+T2))/T2-
     .       0.5d0*(Logmu2*Logt)+0.5d0*(Logmu2*LogT2)+
     .       0.5d0*(Logt*LogT2)-0.5d0*(Logt*(mu2-t-T2))/T2-
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2+
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T*(mu2+t-T2))/
     .            (deltmu2T2T*T2)+
     .             (mu2*tmp103)/(deltmu2T2T*T2)))/mu2)-
     .      (-mu2-t+T2)*
     .       (2d0/T2-Logmu2/T2-Logt/T2-
     .       0.5d0*(Logt*(mu2-t-T2))/T2**2-
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2**2-
     .       0.5d0*((2d0*mu2*phimu2T2T)/T2+
     .           4d0*(-mu2-t+T2)*
     .            ((mu2*phimu2T2T*(mu2+t-T2))/
     .             (deltmu2T2T*T2)+
     .            (mu2*tmp103)/(deltmu2T2T*T2))+
     .           deltmu2T2T*
     .            (((2d0-Logmu2-Logt+2d0*LogT2)*mu2)/
     .             (deltmu2T2T*T2)-
     .            (mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .            (2d0*mu2*phimu2T2T*(mu2+t-T2)*
     .               (-mu2-t+T2))/(deltmu2T2T**2*T2)-
     .            (mu2*tmp103)/(deltmu2T2T*T2**2)-
     .            (2d0*mu2*(-mu2-t+T2)*tmp103)/
     .             (deltmu2T2T**2*T2)+
     .            ((mu2+t-T2)*
     .               ((mu2*phimu2T2T*(mu2+t-T2))/
     .                  (deltmu2T2T*T2)+
     .                 (mu2*tmp103)/(deltmu2T2T*T2)))/
     .             deltmu2T2T))/mu2))-
     .   0.5d0*(((4d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**2)+8d0/T2-
     .      (4d0*LogT2)/T2-
     .      2d0*((B2**2*(LogB2-LogT2))/((1d0-B2/T2)**2*T2**4)+
     .         B2/((1d0-B2/T2)*T2**3)+
     .         (2d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**3))*
     .       (-B2+T2)-(-2d0*B2*LogB2+4d0*T2)/T2**2+
     .      (B2-T2)*tmp18)*tmp189)-
     .   0.5d0*(((4d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**2)+8d0/T2-
     .      (4d0*LogT2)/T2-
     .      2d0*((B1**2*(LogB1-LogT2))/((1d0-B1/T2)**2*T2**4)+
     .         B1/((1d0-B1/T2)*T2**3)+
     .         (2d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**3))*
     .       (-B1+T2)-(-2d0*B1*LogB1+4d0*T2)/T2**2+
     .      (B1-T2)*tmp18)*tmp190)
      DT2T2 = DT2T2-tmp198*
     .    (2d0/T2-LogA0/T2-LogB1/T2-
     .      0.5d0*(LogB1*(A0-B1-T2))/T2**2-
     .      0.5d0*(LogA0*(-A0+B1-T2))/T2**2-
     .      0.5d0*((2d0*B1*phiA0B1T2)/T2+
     .        4d0*(-A0-B1+T2)*
     .         ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .           (B1*tmp100)/(deltA0B1T2*T2))+
     .        deltA0B1T2*
     .         ((B1*(2d0-LogA0-LogB1+2d0*LogT2))/
     .            (deltA0B1T2*T2)-
     .           (B1*phiA0B1T2)/(deltA0B1T2*T2)-
     .           (2d0*B1*phiA0B1T2*(A0+B1-T2)*(-A0-B1+T2))/
     .            (deltA0B1T2**2*T2)-
     .           (B1*tmp100)/(deltA0B1T2*T2**2)-
     .           (2d0*B1*(-A0-B1+T2)*tmp100)/
     .            (deltA0B1T2**2*T2)+
     .           ((A0+B1-T2)*
     .            ((B1*phiA0B1T2*(A0+B1-T2))/
     .               (deltA0B1T2*T2)+
     .              (B1*tmp100)/(deltA0B1T2*T2)))/deltA0B1T2))/
     .      B1)-tmp197*
     .    (2d0/T2-LogA0/T2-LogB2/T2-
     .      0.5d0*(LogB2*(A0-B2-T2))/T2**2-
     .      0.5d0*(LogA0*(-A0+B2-T2))/T2**2-
     .      0.5d0*((2d0*B2*phiA0B2T2)/T2+
     .        4d0*(-A0-B2+T2)*
     .         ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .           (B2*tmp101)/(deltA0B2T2*T2))+
     .        deltA0B2T2*
     .         ((B2*(2d0-LogA0-LogB2+2d0*LogT2))/
     .            (deltA0B2T2*T2)-
     .           (B2*phiA0B2T2)/(deltA0B2T2*T2)-
     .           (2d0*B2*phiA0B2T2*(A0+B2-T2)*(-A0-B2+T2))/
     .            (deltA0B2T2**2*T2)-
     .           (B2*tmp101)/(deltA0B2T2*T2**2)-
     .           (2d0*B2*(-A0-B2+T2)*tmp101)/
     .            (deltA0B2T2**2*T2)+
     .           ((A0+B2-T2)*
     .            ((B2*phiA0B2T2*(A0+B2-T2))/
     .               (deltA0B2T2*T2)+
     .              (B2*tmp101)/(deltA0B2T2*T2)))/deltA0B2T2))/
     .      B2)+tmp63*(-((b*(-1d0+Logb))/T2)-
     .      ((-1d0+Logmu2)*mu2)/T2+2d0*tmp225-
     .      (-b-mu2+T2)*
     .       (2d0/T2-Logb/T2-Logmu2/T2-
     .       0.5d0*(Logmu2*(b-mu2-T2))/T2**2-
     .       0.5d0*(Logb*(-b+mu2-T2))/T2**2-
     .       0.5d0*(2d0*phiT2bmu2+
     .           4d0*(-b-mu2+T2)*
     .            ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .            (mu2*tmp102)/(delT2Tbmu2*T2))+
     .           delT2Tbmu2*tmp226)/mu2))-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (-2d0/T2+Logb/T2+Logmu2/T2+
     .      0.5d0*(Logmu2*(b-mu2-T2))/T2**2+
     .      0.5d0*(Logb*(-b+mu2-T2))/T2**2+
     .      0.5d0*(2d0*phiT2bmu2+
     .        4d0*(-b-mu2+T2)*
     .         ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .           (mu2*tmp102)/(delT2Tbmu2*T2))+
     .        delT2Tbmu2*tmp226)/mu2)
      DT2T2 = DT2T2+ht**2*
     .    (((-1d0+LogT1)*T1*tmp8)/T2+
     .      (A0*cbe**2*(-1d0+LogA0)*(1d0+0.5d0*(1d0+c2t)))/T2+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))/T2+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))/T2+
     .      0.25d0*((1d0+Nc)*s2t**2*
     .       (8d0*(-1d0+LogT2)+2d0*(-1d0+LogT2)**2+
     .         (2d0/T2**2-(2d0*(-1d0+LogT2))/T2**2)*T2**2)))
      DT2T2 = DT2T2-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*
     .       (4d0/T2-(2d0*LogT1)/T2+
     .         (4d0*(-LogT1+LogT2))/(T1*(1d0-T2/T1))-
     .         2d0*(T1-T2)*
     .          ((-LogT1+LogT2)/(T1**2*(1d0-T2/T1)**2)+
     .            1/(T1*T2*(1d0-T2/T1))))*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (-((deltA0T1T2*phiA0T1T2)/T2**3)+2d0/T2-
     .         LogA0/T2-LogT1/T2+
     .         (2d0*phiA0T1T2*(-A0-T1+T2)+
     .            deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2**2-
     .         0.5d0*(LogT1*(A0-T1-T2))/T2**2-
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2**2-
     .         0.5d0*(2d0*phiA0T1T2+
     .             4d0*(-A0-T1+T2)*
     .            (tmp104/deltA0T1T2+
     .              (phiA0T1T2*tmp179)/deltA0T1T2)+
     .             deltA0T1T2*
     .            ((2d0-LogA0-LogT1+2d0*LogT2)/deltA0T1T2-
     .              (phiA0T1T2*(A0-T1)**2)/
     .               (deltA0T1T2*T2**2)-
     .              (2d0*(-A0-T1+T2)*tmp104)/deltA0T1T2**2-
     .              (2d0*phiA0T1T2*(-A0-T1+T2)*tmp179)/
     .               deltA0T1T2**2+
     .              (tmp179*
     .                 (tmp104/deltA0T1T2+
     .                   (phiA0T1T2*tmp179)/deltA0T1T2))/
     .               deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*(4d0/T2+(2d0*(2d0-2d0*LogT2))/T2-
     .           (2d0*LogT2)/T2-(4d0*T2-2d0*LogT2*T2)/T2**2)*tmp34)
     .       +cbe**2*tmp50*
     .       (-((deltA0T2T2*phiA0T2T2)/T2**3)+
     .         (A0*LogA0)/T2**2+4d0/T2-(4d0*LogT2)/T2+
     .         (-4d0*A0*phiA0T2T2+
     .            deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2**2+
     .         0.5d0*((A0-2d0*T2)*tmp18)-
     .         0.5d0*(-8d0*A0*
     .            ((A0*phiA0T2T2)/deltA0T2T2+
     .              tmp107/deltA0T2T2+
     .              (phiA0T2T2*tmp182)/deltA0T2T2+
     .              tmp99/deltA0T2T2)+
     .             deltA0T2T2*
     .            ((4d0*A0**2*phiA0T2T2)/deltA0T2T2**2+
     .              (phiA0T2T2*
     .                 (-1d0-(A0-T2)**2/T2**2-
     .                   (2d0*(A0-T2))/T2))/deltA0T2T2+
     .              (1d0+(A0-T2)/T2)/deltA0T2T2+
     .              (1d0-(-A0+T2)/T2)/deltA0T2T2+
     .              (4d0*A0*tmp107)/deltA0T2T2**2+
     .              (4d0*A0*phiA0T2T2*tmp182)/deltA0T2T2**2+
     .              (4d0*A0*tmp99)/deltA0T2T2**2+
     .              (A0*
     .                 ((A0*phiA0T2T2)/deltA0T2T2+
     .                   tmp107/deltA0T2T2+
     .                   (phiA0T2T2*tmp182)/deltA0T2T2+
     .                   tmp99/deltA0T2T2))/deltA0T2T2+
     .              (tmp182*
     .                 ((A0*phiA0T2T2)/deltA0T2T2+
     .                   tmp107/deltA0T2T2+
     .                   (phiA0T2T2*tmp182)/deltA0T2T2+
     .                   tmp99/deltA0T2T2))/deltA0T2T2))/T2)))
      tmp245 = (2d0*(-1d0+Logmu2)*mu2)/t-((-1d0+LogT1)*T1)/t-
     .   ((-1d0+LogT2)*T2)/t+2d0*tmp228+2d0*tmp243+
     .   sbe**2*(8d0*(-1d0+Logt)+2d0*(-1d0+Logt)**2+
     .      t*(4d0/t+(2d0*(2d0-2d0*Logt))/t-(2d0*Logt)/t-
     .       (4d0*t-2d0*Logt*t)/t**2)+t**2*tmp15+2d0*tmp77)-
     .   (-mu2-t+T2)*(2d0/t-Logmu2/t-LogT2/t-
     .      0.5d0*(LogT2*(mu2-t-T2))/t**2-
     .      0.5d0*(Logmu2*(-mu2-t+T2))/t**2-
     .      0.5d0*((2d0*mu2*phimu2T2T)/T2+
     .        4d0*(-mu2+t-T2)*
     .         ((mu2*phimu2T2T*(mu2-t+T2))/(deltmu2T2T*T2)+
     .           (mu2*tmp111)/(deltmu2T2T*t))+
     .        deltmu2T2T*
     .         (((2d0-Logmu2+2d0*Logt-LogT2)*mu2)/
     .            (deltmu2T2T*t)-
     .           (mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .           (2d0*mu2*phimu2T2T*(-mu2+t-T2)*
     .            (mu2-t+T2))/(deltmu2T2T**2*T2)-
     .           (mu2*tmp111)/(deltmu2T2T*t**2)-
     .           (2d0*mu2*(-mu2+t-T2)*tmp111)/
     .            (deltmu2T2T**2*t)+
     .           ((mu2-t+T2)*
     .            ((mu2*phimu2T2T*(mu2-t+T2))/
     .               (deltmu2T2T*T2)+
     .              (mu2*tmp111)/(deltmu2T2T*t)))/deltmu2T2T))/
     .      mu2)+cbe**2*
     .    (8d0*(-1d0+Logt)+2d0*(-1d0+Logt)**2-
     .      (2d0*A0*(-1d0+LogA0))/t+t**2*tmp15+
     .      4d0*(-1d0+4d0*Logt-Logt**2-(A0*LogA0)/t+
     .       (2d0*A0*phiA0tt)/t+(Logt*(A0-2d0*t))/t+
     .       0.5d0*(deltA0tt*phiA0tt)/t**2-
     .       0.5d0*(deltA0tt*
     .           ((A0*phiA0tt)/deltA0tt+
     .             (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .             tmp70/deltA0tt))/t)-
     .      (A0-2d0*t)*(-((deltA0tt*phiA0tt)/t**3)+
     .       (A0*LogA0)/t**2+4d0/t-(4d0*Logt)/t+
     .       (-4d0*A0*phiA0tt+
     .          deltA0tt*
     .           ((A0*phiA0tt)/deltA0tt+
     .             (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .             tmp70/deltA0tt))/t**2+
     .       0.5d0*((A0-2d0*t)*tmp16)-
     .       0.5d0*(-8d0*A0*((A0*phiA0tt)/deltA0tt+
     .            (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .            tmp70/deltA0tt)+
     .           deltA0tt*
     .            ((4d0*A0**2*phiA0tt)/deltA0tt**2+
     .            (phiA0tt*
     .               (-1d0-(A0-t)**2/t**2-(2d0*(A0-t))/t))/
     .             deltA0tt+(1d0+(A0-t)/t)/deltA0tt+
     .            (1d0-(-A0+t)/t)/deltA0tt+
     .            (4d0*A0*phiA0tt*tmp175)/deltA0tt**2+
     .            (4d0*A0*tmp65)/deltA0tt**2+
     .            (4d0*A0*tmp70)/deltA0tt**2+
     .            (A0*((A0*phiA0tt)/deltA0tt+
     .                 (phiA0tt*tmp175)/deltA0tt+
     .                 tmp65/deltA0tt+tmp70/deltA0tt))/
     .             deltA0tt+
     .            (tmp175*
     .               ((A0*phiA0tt)/deltA0tt+
     .                 (phiA0tt*tmp175)/deltA0tt+
     .                 tmp65/deltA0tt+tmp70/deltA0tt))/
     .             deltA0tt))/t))
      tmp245 = tmp245-
     .   (-mu2-t+T1)*(2d0/t-Logmu2/t-LogT1/t-
     .      0.5d0*(LogT1*(mu2-t-T1))/t**2-
     .      0.5d0*(Logmu2*(-mu2-t+T1))/t**2-
     .      0.5d0*((2d0*mu2*phimu2tT1)/T1+
     .        4d0*(-mu2+t-T1)*
     .         ((mu2*phimu2tT1*(mu2-t+T1))/(deltmu2tT1*T1)+
     .           (mu2*tmp91)/(deltmu2tT1*t))+
     .        deltmu2tT1*
     .         (((2d0-Logmu2+2d0*Logt-LogT1)*mu2)/
     .            (deltmu2tT1*t)-
     .           (mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .           (2d0*mu2*phimu2tT1*(-mu2+t-T1)*
     .            (mu2-t+T1))/(deltmu2tT1**2*T1)-
     .           (mu2*tmp91)/(deltmu2tT1*t**2)-
     .           (2d0*mu2*(-mu2+t-T1)*tmp91)/
     .            (deltmu2tT1**2*t)+
     .           ((mu2-t+T1)*
     .            ((mu2*phimu2tT1*(mu2-t+T1))/
     .               (deltmu2tT1*T1)+
     .              (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))/
     .      mu2)
      Dtt = ht**2*tmp245+
     .   tmp10*(-5d0-2d0*Li2bt+4d0*Logt-Logt**2+
     .      (b*(-1d0+Logb))/t+(2d0*Logt*(b-t))/t+
     .      (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)+
     .      (-2d0*b*Logb+4d0*t)/t+
     .      0.5d0*((b+t)*((4d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+
     .         8d0/t-(4d0*Logt)/t-(-2d0*b*Logb+4d0*t)/t**2+
     .         (b-t)*tmp16-2d0*(-b+t)*tmp207)))+
     .   tmp61*(-((B1*(-1d0+LogB1))/t)+((-1d0+Logmu2)*mu2)/t+
     .      2d0*tmp231-(B1-mu2-t)*
     .       (2d0/t-LogB1/t-Logmu2/t-
     .       0.5d0*(Logmu2*(B1-mu2-t))/t**2-
     .       0.5d0*(LogB1*(-B1+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB1tmu2+deltB1tmu2*tmp232+
     .           4d0*(-B1-mu2+t)*
     .            ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2))+
     .   tmp60*(-((B2*(-1d0+LogB2))/t)+((-1d0+Logmu2)*mu2)/t+
     .      2d0*tmp234-(B2-mu2-t)*
     .       (2d0/t-LogB2/t-Logmu2/t-
     .       0.5d0*(Logmu2*(B2-mu2-t))/t**2-
     .       0.5d0*(LogB2*(-B2+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB2tmu2+deltB2tmu2*tmp236+
     .           4d0*(-B2-mu2+t)*
     .            ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2))-
     .   2d0*hb*ht*mu*s2b*(tmp235/mt+
     .      mt*(-(LogB1/t)+LogB2/t-
     .       0.5d0*(Logmu2*(B1-mu2-t))/t**2+
     .       0.5d0*(Logmu2*(B2-mu2-t))/t**2-
     .       0.5d0*(LogB1*(-B1+mu2-t))/t**2+
     .       0.5d0*(LogB2*(-B2+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB1tmu2+deltB1tmu2*tmp232+
     .           4d0*(-B1-mu2+t)*
     .            ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2+
     .       0.5d0*(2d0*phiB2tmu2+deltB2tmu2*tmp236+
     .           4d0*(-B2-mu2+t)*
     .            ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2)-
     .      0.25d0*tmp75/t**1.5d0)
      Dtt = Dtt+tmp9*
     .    (-((A0*(-1d0+LogA0))/t)+(b*(-1d0+Logb))/t+2d0*tmp229-
     .      (A0-b-t)*(2d0/t-LogA0/t-Logb/t-
     .       0.5d0*(Logb*(A0-b-t))/t**2-
     .       0.5d0*(LogA0*(-A0+b-t))/t**2-
     .       0.5d0*((2d0*b*phiA0bt)/t+deltA0bt*tmp230+
     .           4d0*(-A0-b+t)*
     .            ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b))-
     .   4d0*cbe*hb*ht*mb*sbe*
     .    ((0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .       0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .       0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .       0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .             (b*tmp66)/(deltA0bt*t)))/b)/mt+
     .      mt*(-2d0/t+LogA0/t+Logb/t+
     .       0.5d0*(Logb*(A0-b-t))/t**2+
     .       0.5d0*(LogA0*(-A0+b-t))/t**2+
     .       0.5d0*((4d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+8d0/t-
     .          (4d0*Logt)/t-(-2d0*b*Logb+4d0*t)/t**2+
     .          (b-t)*tmp16-2d0*(-b+t)*tmp207)+
     .       0.5d0*((2d0*b*phiA0bt)/t+deltA0bt*tmp230+
     .           4d0*(-A0-b+t)*
     .            ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b)-
     .      0.25d0*(-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .        0.5d0*(Logb*Logt*(A0-b-t))-
     .        0.5d0*(LogA0*Logt*(-A0+b-t))+
     .        0.5d0*(deltA0bt*phiA0bt)/t-
     .        0.5d0*(LogA0*Logb*(-A0-b+t))+
     .        2.5d0*(A0+b+t)+0.5d0*tmp76)/t**1.5d0)+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp130)/t**1.5d0
     .        -0.5d0*(mb*tmp56)/t**1.5d0)-
     .        cbe**2*hb**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .         0.125d0*((1d0+c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .         0.125d0*((1d0-c2b)*s2t*Xt)/t**1.5d0)))
      Dtt = Dtt+tmp89*
     .    (-(cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp156)/t**1.5d0
     .        -0.5d0*(mb*tmp56)/t**1.5d0))-
     .      hb**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .       0.125d0*((1d0+c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .       0.125d0*((1d0-c2b)*s2t*Yt)/t**1.5d0))+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp131)/t**1.5d0
     .        -0.5d0*(mb*tmp59)/t**1.5d0)-
     .        cbe**2*hb**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .         0.125d0*((1d0+c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .         0.125d0*((1d0-c2b)*s2t*Xt)/t**1.5d0)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp130)/t**1.5d0)-
     .         0.5d0*(mb*tmp59)/t**1.5d0)-
     .      cbe**2*hb**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .         0.125d0*((1d0-c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .         0.125d0*((1d0+c2b)*s2t*Xt)/t**1.5d0)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp131)/t**1.5d0)-
     .         0.5d0*(mb*tmp56)/t**1.5d0)-
     .      cbe**2*hb**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .         0.125d0*((1d0-c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .         0.125d0*((1d0+c2b)*s2t*Xt)/t**1.5d0)))-
     .   0.25d0*(ht**2*((cbe**2*s2t*tmp114*Yt)/t**1.5d0-
     .      (cbe**2*s2t*tmp93*Yt)/t**1.5d0+
     .      0.5d0*(s2t*sbe**2*tmp115*Xt)/t**1.5d0-
     .      0.5d0*(s2t*sbe**2*tmp94*Xt)/t**1.5d0))
      Dtt = Dtt+tmp109*
     .    (-(cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp157)/t**1.5d0
     .        -0.5d0*(mb*tmp59)/t**1.5d0))-
     .      hb**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .       0.125d0*((1d0+c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .       0.125d0*((1d0-c2b)*s2t*Yt)/t**1.5d0))+
     .   tmp88*(-(cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp156)/t**1.5d0)-
     .         0.5d0*(mb*tmp59)/t**1.5d0))-
     .      hb**2*sbe**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .       0.125d0*((1d0-c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .       0.125d0*((1d0+c2b)*s2t*Yt)/t**1.5d0))+
     .   tmp108*(-(cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp157)/t**1.5d0)-
     .         0.5d0*(mb*tmp56)/t**1.5d0))-
     .      hb**2*sbe**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .       0.125d0*((1d0-c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .       0.125d0*((1d0+c2b)*s2t*Yt)/t**1.5d0))
      Dc2tc2t = (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT2)*T2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2*T2-
     .      (-b-mu2+T2)*tmp112)*
     .    (0.125d0*hb**2/c2t**3-0.125d0*ht**2/c2t**3)+
     .   (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT1)*T1-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2*T1-
     .      (-b-mu2+T1)*tmp92)*
     .    (-(0.125d0*hb**2/c2t**3)+0.125d0*ht**2/c2t**3)+
     .   ht**2*(0.0625d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t**3+0.0625d0*
     .       (B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t**3+
     .      cbe**2*(0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/
     .         c2t**3-0.125d0*
     .        (A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t**3)-
     .      0.0625d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t**3-0.0625d0*
     .       (B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t**3)+
     .   hb**2*(-(0.0625d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*
     .          T1)/c2t**3)-
     .      0.0625d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/
     .      c2t**3+sbe**2*
     .       (-(0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t**3)+
     .       0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t**3)+
     .      0.0625d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t**3+0.0625d0*
     .       (B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t**3)+
     .   0.5d0*(hb*ht*mb*mu*tmp113)/s2t**3+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe**2*hb**2*
     .         (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .           0.125d0*(mb*mt*s2b)/s2t**3-
     .           0.0625d0*((1d0+c2b)*t)/c2t**3-
     .           0.125d0*(mb*s2b*Xb)/c2t**3+
     .           0.125d0*((1d0+c2b)*mt*Xb)/s2t**3+
     .           0.0625d0*((1d0+c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (-(mt*s2b*tmp25)+2d0*mb*mt*tmp52-
     .         0.125d0*(s2b*(b+t))/s2t**3+
     .         0.25d0*(mb*tmp128)/s2t**3-0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0-c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Xt)/c2t**3+
     .         0.125d0*((1d0-c2b)*mt*Xt)/s2t**3-
     .         0.0625d0*((1d0-c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe**2*hb**2*
     .         (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .           0.125d0*(mb*mt*s2b)/s2t**3+
     .           0.0625d0*((1d0+c2b)*t)/c2t**3+
     .           0.125d0*(mb*s2b*Xb)/c2t**3-
     .           0.125d0*((1d0+c2b)*mt*Xb)/s2t**3-
     .           0.0625d0*((1d0+c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (-(mt*s2b*tmp26)+2d0*mb*mt*tmp53+
     .         0.125d0*(s2b*(b+t))/s2t**3-
     .         0.25d0*(mb*tmp128)/s2t**3+0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0-c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Xt)/c2t**3-
     .         0.125d0*((1d0-c2b)*mt*Xt)/s2t**3+
     .         0.0625d0*((1d0-c2b)*Xt**2)/c2t**3)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe**2*hb**2*
     .         (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .           0.125d0*(mb*mt*s2b)/s2t**3-
     .           0.0625d0*((1d0-c2b)*t)/c2t**3+
     .           0.125d0*(mb*s2b*Xb)/c2t**3+
     .           0.125d0*((1d0-c2b)*mt*Xb)/s2t**3+
     .           0.0625d0*((1d0-c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (mt*s2b*tmp25+2d0*mb*mt*tmp53+
     .         0.125d0*(s2b*(b+t))/s2t**3+
     .         0.25d0*(mb*tmp129)/s2t**3+0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0+c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Xt)/c2t**3+
     .         0.125d0*((1d0+c2b)*mt*Xt)/s2t**3-
     .         0.0625d0*((1d0+c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   tmp109*(-(hb**2*sbe**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0+c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Yb)/c2t**3+
     .         0.125d0*((1d0+c2b)*mt*Yb)/s2t**3+
     .         0.0625d0*((1d0+c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(-(mt*s2b*tmp41)+2d0*mb*mt*tmp52-
     .       0.125d0*(s2b*(b+t))/s2t**3+
     .       0.25d0*(mb*tmp154)/s2t**3-0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .       0.125d0*(mb*mt*s2b)/s2t**3+
     .       0.0625d0*((1d0-c2b)*t)/c2t**3+
     .       0.125d0*(mb*s2b*Yt)/c2t**3+
     .       0.125d0*((1d0-c2b)*mt*Yt)/s2t**3-
     .       0.0625d0*((1d0-c2b)*Yt**2)/c2t**3))-
     .   0.25d0*(ht**2*((cbe**2*mt*tmp114*Yt)/s2t**3-
     .      (cbe**2*mt*tmp93*Yt)/s2t**3+
     .      0.5d0*(mt*sbe**2*tmp115*Xt)/s2t**3-
     .      0.5d0*(mt*sbe**2*tmp94*Xt)/s2t**3))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe**2*hb**2*
     .         (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .           0.125d0*(mb*mt*s2b)/s2t**3+
     .           0.0625d0*((1d0-c2b)*t)/c2t**3-
     .           0.125d0*(mb*s2b*Xb)/c2t**3-
     .           0.125d0*((1d0-c2b)*mt*Xb)/s2t**3-
     .           0.0625d0*((1d0-c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (mt*s2b*tmp26+2d0*mb*mt*tmp52-
     .         0.125d0*(s2b*(b+t))/s2t**3-
     .         0.25d0*(mb*tmp129)/s2t**3-0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0+c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Xt)/c2t**3-
     .         0.125d0*((1d0+c2b)*mt*Xt)/s2t**3+
     .         0.0625d0*((1d0+c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   tmp89*(-(hb**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0+c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Yb)/c2t**3-
     .         0.125d0*((1d0+c2b)*mt*Yb)/s2t**3-
     .         0.0625d0*((1d0+c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(-(mt*s2b*tmp42)+2d0*mb*mt*tmp53+
     .       0.125d0*(s2b*(b+t))/s2t**3-
     .       0.25d0*(mb*tmp154)/s2t**3+0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .       0.125d0*(mb*mt*s2b)/s2t**3-
     .       0.0625d0*((1d0-c2b)*t)/c2t**3-
     .       0.125d0*(mb*s2b*Yt)/c2t**3-
     .       0.125d0*((1d0-c2b)*mt*Yt)/s2t**3+
     .       0.0625d0*((1d0-c2b)*Yt**2)/c2t**3))+
     .   tmp108*(-(hb**2*sbe**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0-c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Yb)/c2t**3+
     .         0.125d0*((1d0-c2b)*mt*Yb)/s2t**3+
     .         0.0625d0*((1d0-c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(mt*s2b*tmp41+2d0*mb*mt*tmp53+
     .       0.125d0*(s2b*(b+t))/s2t**3+
     .       0.25d0*(mb*tmp155)/s2t**3+0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .       0.125d0*(mb*mt*s2b)/s2t**3+
     .       0.0625d0*((1d0+c2b)*t)/c2t**3-
     .       0.125d0*(mb*s2b*Yt)/c2t**3+
     .       0.125d0*((1d0+c2b)*mt*Yt)/s2t**3-
     .       0.0625d0*((1d0+c2b)*Yt**2)/c2t**3))+
     .   tmp88*(-(hb**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0-c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Yb)/c2t**3-
     .         0.125d0*((1d0-c2b)*mt*Yb)/s2t**3-
     .         0.0625d0*((1d0-c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(mt*s2b*tmp42+2d0*mb*mt*tmp52-
     .       0.125d0*(s2b*(b+t))/s2t**3-
     .       0.25d0*(mb*tmp155)/s2t**3-0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .       0.125d0*(mb*mt*s2b)/s2t**3-
     .       0.0625d0*((1d0+c2b)*t)/c2t**3+
     .       0.125d0*(mb*s2b*Yt)/c2t**3-
     .       0.125d0*((1d0+c2b)*mt*Yt)/s2t**3+
     .       0.0625d0*((1d0+c2b)*Yt**2)/c2t**3))
      DT1t = -(tmp192*tmp238)-tmp191*tmp239-
     .   0.5d0*(tmp184*tmp213)-0.5d0*(tmp183*tmp214)+
     .   ht**2*(1d0-3d0*Logt+Logmu2*Logt-
     .      (-1d0+Logt)*(-1d0+LogT1)+LogT1-Logmu2*LogT1+
     .      (phimu2tT1*(-mu2+t-T1))/T1-
     .      (phimu2tT1*(-mu2-t+T1))/T1-
     .      0.5d0*(LogT1*(mu2-t-T1))/t+
     .      0.5d0*(Logt*(mu2-t-T1))/T1+
     .      0.5d0*(Logmu2*(-mu2+t-T1))/T1-
     .      0.5d0*(Logmu2*(-mu2-t+T1))/t-
     .      0.5d0*(deltmu2tT1*
     .        ((mu2*phimu2tT1*(mu2+t-T1))/(deltmu2tT1*T1)+
     .          (mu2*tmp84)/(deltmu2tT1*T1)))/mu2+
     .      0.5d0*(deltmu2tT1*
     .        ((mu2*phimu2tT1*(mu2-t+T1))/(deltmu2tT1*T1)+
     .          (mu2*tmp91)/(deltmu2tT1*t)))/mu2-
     .      (-mu2-t+T1)*
     .       (phimu2tT1/T1-
     .       ((-mu2+t-T1)*
     .          ((mu2*phimu2tT1*(mu2+t-T1))/
     .             (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1))
     .          )/mu2-((-mu2-t+T1)*
     .          ((mu2*phimu2tT1*(mu2-t+T1))/
     .             (deltmu2tT1*T1)+(mu2*tmp91)/(deltmu2tT1*t)))
     .         /mu2+0.5d0*Logmu2/t-0.5d0*LogT1/t+
     .       0.5d0*Logmu2/T1-0.5d0*Logt/T1+
     .       0.5d0*(mu2-t-T1)/(t*T1)-
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .             (2d0*mu2*phimu2tT1*(-mu2+t-T1)*
     .              (mu2+t-T1))/(deltmu2tT1**2*T1)+
     .             (mu2*(Logmu2-Logt-(-mu2+t)/t-T1/t))/
     .            (deltmu2tT1*T1)-
     .             (2d0*mu2*(-mu2+t-T1)*tmp84)/
     .            (deltmu2tT1**2*T1)+
     .             ((mu2+t-T1)*
     .              ((mu2*phimu2tT1*(mu2-t+T1))/
     .                 (deltmu2tT1*T1)+
     .                (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))
     .          /mu2))-
     .   0.25d0*(ht**2*(cbe**2*(4d0+(2d0*s2t*Yt)/mt)*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp95*(4d0+(2d0*s2t*Xt)/mt))))
      DT2t = -(tmp196*tmp222)-tmp195*tmp223+
     .   ht**2*(1d0-3d0*Logt+Logmu2*Logt-
     .      (-1d0+Logt)*(-1d0+LogT2)+LogT2-Logmu2*LogT2+
     .      (phimu2T2T*(-mu2+t-T2))/T2-
     .      (phimu2T2T*(-mu2-t+T2))/T2-
     .      0.5d0*(LogT2*(mu2-t-T2))/t+
     .      0.5d0*(Logt*(mu2-t-T2))/T2+
     .      0.5d0*(Logmu2*(-mu2+t-T2))/T2-
     .      0.5d0*(Logmu2*(-mu2-t+T2))/t-
     .      0.5d0*(deltmu2T2T*
     .        ((mu2*phimu2T2T*(mu2+t-T2))/(deltmu2T2T*T2)+
     .          (mu2*tmp103)/(deltmu2T2T*T2)))/mu2+
     .      0.5d0*(deltmu2T2T*
     .        ((mu2*phimu2T2T*(mu2-t+T2))/(deltmu2T2T*T2)+
     .          (mu2*tmp111)/(deltmu2T2T*t)))/mu2-
     .      (-mu2-t+T2)*
     .       (phimu2T2T/T2-
     .       ((-mu2+t-T2)*
     .          ((mu2*phimu2T2T*(mu2+t-T2))/
     .             (deltmu2T2T*T2)+(mu2*tmp103)/(deltmu2T2T*T2)
     .            ))/mu2-
     .       ((-mu2-t+T2)*
     .          ((mu2*phimu2T2T*(mu2-t+T2))/
     .             (deltmu2T2T*T2)+(mu2*tmp111)/(deltmu2T2T*t))
     .          )/mu2+0.5d0*Logmu2/t-0.5d0*LogT2/t+
     .       0.5d0*Logmu2/T2-0.5d0*Logt/T2+
     .       0.5d0*(mu2-t-T2)/(t*T2)-
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .             (2d0*mu2*phimu2T2T*(-mu2+t-T2)*
     .              (mu2+t-T2))/(deltmu2T2T**2*T2)+
     .             (mu2*(Logmu2-Logt-(-mu2+t)/t-T2/t))/
     .            (deltmu2T2T*T2)-
     .             (2d0*mu2*(-mu2+t-T2)*tmp103)/
     .            (deltmu2T2T**2*T2)+
     .             ((mu2+t-T2)*
     .              ((mu2*phimu2T2T*(mu2-t+T2))/
     .                 (deltmu2T2T*T2)+
     .                (mu2*tmp111)/(deltmu2T2T*t)))/deltmu2T2T)
     .           )/mu2))-0.5d0*(tmp188*tmp217)-
     .   0.5d0*(tmp187*tmp218)-
     .   0.25d0*(ht**2*(cbe**2*(4d0-(2d0*s2t*Yt)/mt)*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*(4d0-(2d0*s2t*Xt)/mt))))
      DT1T2 = ht**2*(c2t**2+(-1d0+LogT1)*tmp8+
     .      (-1d0+LogT2)*tmp8+(-1d0+LogT1)*(-1d0+LogT2)*tmp8-
     .      0.5d0*((-1d0+Nc)*s2t**2))-
     .   0.25d0*(ht**2*((1d0+c2t**2)*sbe**2*
     .       ((2d0*LogT1)/T1+(-2d0-2d0*LogT2)/T1+
     .         (2d0*(-LogT1+LogT2)*(T1-T2)*T2)/
     .          (T1**3*(1d0-T2/T1)**2)-
     .         (2d0*(-LogT1+LogT2))/(T1*(1d0-T2/T1))+
     .         (2d0*(T1-T2))/(T1**2*(1d0-T2/T1))+
     .         (2d0*(-LogT1+LogT2)*(T1-T2))/
     .          (T1**2*(1d0-T2/T1))-
     .         (2d0*(-LogT1+LogT2)*T2)/(T1**2*(1d0-T2/T1)))*Xt**2
     .       +2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       ((phiA0T1T2*(-A0+T1-T2))/T2**2+phiA0T1T2/T2-
     .         ((-A0+T1-T2)*
     .            (tmp104/deltA0T1T2+
     .            (phiA0T1T2*tmp179)/deltA0T1T2))/T2-
     .         ((-A0-T1+T2)*
     .            ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .            (T2*tmp96)/(deltA0T1T2*T1)))/T2+
     .         0.5d0*LogA0/T1-0.5d0*LogT2/T1+0.5d0*LogA0/T2-
     .         0.5d0*LogT1/T2+0.5d0*(A0-T1-T2)/(T1*T2)+
     .         0.5d0*(deltA0T1T2*
     .             ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .             (T2*tmp96)/(deltA0T1T2*T1)))/T2**2-
     .         0.5d0*(deltA0T1T2*
     .             (phiA0T1T2/deltA0T1T2+
     .             ((LogA0-LogT2-T1/T2+(A0-T2)/T2)*T2)/
     .              (deltA0T1T2*T1)-
     .             (2d0*phiA0T1T2*(-A0-T1+T2)*
     .                (A0-T1+T2))/deltA0T1T2**2+
     .             ((A0-T1+T2)*
     .                (tmp104/deltA0T1T2+
     .                  (phiA0T1T2*tmp179)/deltA0T1T2))/
     .              deltA0T1T2+tmp96/(deltA0T1T2*T1)-
     .             (2d0*T2*(-A0-T1+T2)*tmp96)/
     .              (deltA0T1T2**2*T1)))/T2)))
      Dtc2t = -(0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .        LogT1**2*(B2-T1)-5d0*T1+
     .        LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe*hb*ht*sbe*
     .           ((mb*tmp54)/mt+0.25d0*s2b/s2t-
     .             0.5d0*(s2b*tmp27)/mt))+
     .        cbe**2*hb**2*
     .         (-(0.125d0*(1d0+c2b)/c2t)+
     .           0.125d0*(mb*s2b)/(mt*s2t)-
     .           0.125d0*((1d0+c2b)*Xb)/(mt*s2t))+
     .        ht**2*sbe**2*
     .         (0.125d0*(1d0-c2b)/c2t+0.125d0*(mb*s2b)/(mt*s2t)-
     .           0.125d0*((1d0-c2b)*Xt)/(mt*s2t)))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp55)/mt-0.25d0*s2b/s2t-
     .           0.5d0*(s2b*tmp28)/mt))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0+c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0-c2b)*Xt)/(mt*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp55)/mt-0.25d0*s2b/s2t+
     .           0.5d0*(s2b*tmp27)/mt))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mb*s2b)/(mt*s2t)-
     .         0.125d0*((1d0-c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mb*s2b)/(mt*s2t)-
     .         0.125d0*((1d0+c2b)*Xt)/(mt*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp54)/mt+0.25d0*s2b/s2t+
     .           0.5d0*(s2b*tmp28)/mt))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0-c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0+c2b)*Xt)/(mt*s2t))))
      Dtc2t = Dtc2t-tmp89*
     .    (cbe*hb*ht*sbe*((mb*tmp54)/mt+0.25d0*s2b/s2t-
     .       0.5d0*(s2b*tmp43)/mt)+
     .      hb**2*sbe**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0+c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0-c2b)*Yt)/(mt*s2t)))-
     .   tmp109*(cbe*hb*ht*sbe*
     .       ((mb*tmp55)/mt-0.25d0*s2b/s2t-0.5d0*(s2b*tmp44)/mt)+
     .      hb**2*sbe**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0+c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0-c2b)/c2t)-
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0-c2b)*Yt)/(mt*s2t)))-
     .   tmp88*(cbe*hb*ht*sbe*
     .       ((mb*tmp55)/mt-0.25d0*s2b/s2t+0.5d0*(s2b*tmp43)/mt)+
     .      hb**2*sbe**2*(-(0.125d0*(1d0-c2b)/c2t)-
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0-c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0+c2b)*Yt)/(mt*s2t)))-
     .   tmp108*(cbe*hb*ht*sbe*
     .       ((mb*tmp54)/mt+0.25d0*s2b/s2t+0.5d0*(s2b*tmp44)/mt)+
     .      hb**2*sbe**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0-c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0+c2b)*Yt)/(mt*s2t)))-
     .   0.25d0*(ht**2*((cbe**2*tmp114*Yt)/(mt*s2t)-
     .      (cbe**2*tmp93*Yt)/(mt*s2t)+
     .      0.5d0*(sbe**2*tmp115*Xt)/(mt*s2t)-
     .      0.5d0*(sbe**2*tmp94*Xt)/(mt*s2t)))
      DT1c2t = -(tmp204*tmp238)-tmp203*tmp239+
     .   (hb*ht*mb*mu*tmp240)/s2t+
     .   hb**2*(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1))/c2t+
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2))/c2t+
     .      sbe**2*(0.25d0*(A0*(-1d0+LogA0))/c2t+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1))/c2t)+
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1))/c2t+
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1))/c2t)+
     .   tmp6*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT1)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT1)*mu2-
     .      2d0*LogT1*T1-(-b-mu2+T1)*tmp240+
     .      0.5d0*(deltT1bmu2*phiT1bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))-
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))+2.5d0*(b+mu2+T1))
     .    -0.5d0*(tmp200*tmp213)-0.5d0*(tmp199*tmp214)+
     .   ht**2*(-(0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1))/c2t)-
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2))/c2t+
     .      cbe**2*(-(0.25d0*(A0*(-1d0+LogA0))/c2t)-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1))/c2t)-
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1))/c2t-
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1))/c2t+
     .      (-1d0+LogT2)*T2*(1d0+0.5d0*(-1d0+Nc))+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T2*(1d0+0.5d0*(-1d0+Nc))-
     .      0.25d0*((1d0+Nc)*tmp79))-
     .   0.25d0*(ht**2*(sbe**2*tmp220*Xt**2+
     .      2d0*cbe**2*tmp221*Yt**2+
     .      cbe**2*tmp46*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp30*tmp95)))
      DT2c2t = -(tmp206*tmp222)-tmp205*tmp223+
     .   (hb*ht*mb*mu*tmp225)/s2t+
     .   hb**2*(-(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1))/c2t)-
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2))/c2t+
     .      sbe**2*(-(0.25d0*(A0*(-1d0+LogA0))/c2t)-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2))/c2t)-
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2))/c2t-
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2))/c2t)+
     .   tmp7*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT2)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT2)*mu2-
     .      2d0*LogT2*T2-(-b-mu2+T2)*tmp224+
     .      0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2))
     .    -0.5d0*(tmp202*tmp217)-0.5d0*(tmp201*tmp218)+
     .   ht**2*(0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1))/c2t+
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2))/c2t+
     .      cbe**2*(0.25d0*(A0*(-1d0+LogA0))/c2t+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2))/c2t)+
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2))/c2t+
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2))/c2t+
     .      (-1d0+LogT1)*T1*(1d0+0.5d0*(-1d0+Nc))+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T1*(1d0+0.5d0*(-1d0+Nc))-
     .      0.25d0*((1d0+Nc)*tmp98))-
     .   0.25d0*(ht**2*(sbe**2*tmp219*Xt**2+
     .      2d0*cbe**2*Yt**2*
     .       (-0.5d0+2d0*LogT2-(phiA0T1T2*(-A0-T1+T2))/T2+
     .         0.5d0*(LogA0*LogT1)-0.5d0*(LogA0*LogT2)-
     .         0.5d0*(LogT1*LogT2)+
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2**2+
     .         0.5d0*(LogT1*(A0-T1-T2))/T2+
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2-
     .         0.5d0*(deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*tmp31)+
     .      cbe**2*tmp47*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)))
      Dtb = tmp10*(-1d0+Logb+(-1d0+Logb)*(-1d0+Logt)+
     .      Logt+0.5d0*((b+t)*tmp208)
     .      +0.5d0*tmp209+0.5d0*tmp210)-
     .     tmp108*(-(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))-
     .   tmp109*(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt))-
     .   (2d0*cbe*hb*ht*mt*sbe*
     .      (0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .      0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .      0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .      0.5d0*(deltA0bt*
     .          ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b))/mb-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)))
      Dtb = Dtb-tmp89*
     .    (-(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))-
     .   tmp88*(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt))-
     .   (2d0*cbe*hb*ht*mb*sbe*
     .      (0.5d0-2d0*Logb+(phiA0bt*(-A0+b-t))/t+
     .      0.5d0*(LogA0*Logb)-0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logt*(A0-b-t))/b-
     .      0.5d0*(deltA0bt*phiA0bt)/(b*t)-
     .      0.5d0*(LogA0*(-A0-b+t))/b+0.5d0*tmp209+
     .      0.5d0*(deltA0bt*
     .          ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .            tmp69/deltA0bt))/b))/mt-
     .   4d0*cbe*hb*ht*mb*mt*sbe*
     .    (-(phiA0bt/t)-(phiA0bt*(-A0-b+t))/(b*t)+
     .      ((-A0+b-t)*
     .       ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .         (b*tmp66)/(deltA0bt*t)))/b+
     .      ((-A0-b+t)*
     .       ((b*phiA0bt*tmp174)/(deltA0bt*t)+tmp69/deltA0bt))/
     .       b-0.5d0*LogA0/b+0.5d0*Logt/b-0.5d0*LogA0/t+
     .      0.5d0*Logb/t-0.5d0*(A0-b-t)/(b*t)+0.5d0*tmp208-
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b**2+
     .      0.5d0*(deltA0bt*((b*phiA0bt)/(deltA0bt*t)-
     .          (2d0*b*phiA0bt*(-A0+b-t)*(A0+b-t))/
     .           (deltA0bt**2*t)+(b*tmp64)/(deltA0bt*t)+
     .          tmp66/(deltA0bt*t)-
     .          (2d0*b*(-A0+b-t)*tmp66)/(deltA0bt**2*t)+
     .          ((A0+b-t)*
     .             ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .             tmp69/deltA0bt))/deltA0bt))/b)-
     .   (cbe*hb*ht*sbe*(-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76))/(mb*mt)
      Dtb = Dtb+tmp9*
     .    (-2d0+3d0*Logb+(-1d0+Logb)*(-1d0+Logt)+3d0*Logt-
     .      Logb*Logt-(phiA0bt*(-A0+b-t))/t-
     .      (phiA0bt*(-A0-b+t))/t+
     .      0.5d0*(Logt*(A0-b-t))/b+
     .      0.5d0*(deltA0bt*phiA0bt)/(b*t)+
     .      0.5d0*(Logb*(A0-b-t))/t+
     .      0.5d0*(LogA0*(-A0+b-t))/t+
     .      0.5d0*(LogA0*(-A0-b+t))/b-
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b-
     .      0.5d0*(deltA0bt*((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .          tmp69/deltA0bt))/b-
     .      (A0-b-t)*(phiA0bt/t+
     .       (phiA0bt*(-A0-b+t))/(b*t)-
     .       ((-A0+b-t)*
     .          ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b-
     .       ((-A0-b+t)*
     .          ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .            tmp69/deltA0bt))/b+0.5d0*LogA0/b-
     .       0.5d0*Logt/b+0.5d0*LogA0/t-0.5d0*Logb/t+
     .       0.5d0*(A0-b-t)/(b*t)+
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .             (b*tmp66)/(deltA0bt*t)))/b**2-
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt)/(deltA0bt*t)-
     .             (2d0*b*phiA0bt*(-A0+b-t)*(A0+b-t))/
     .            (deltA0bt**2*t)+(b*tmp64)/(deltA0bt*t)+
     .             tmp66/(deltA0bt*t)-
     .             (2d0*b*(-A0+b-t)*tmp66)/(deltA0bt**2*t)+
     .             ((A0+b-t)*
     .              ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .                tmp69/deltA0bt))/deltA0bt))/b))
      DT1b = -((hb*ht*mu*s2t*tmp240)/mb)-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (phiT1bmu2/mu2-
     .      ((b-mu2-T1)*
     .       ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .         (mu2*tmp83)/(deltT1bmu2*T1)))/mu2-
     .      ((-b-mu2+T1)*
     .       ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .         (mu2*tmp90)/(b*deltT1bmu2)))/mu2+0.5d0*Logmu2/b-
     .      0.5d0*LogT1/b-0.5d0*Logb/T1+0.5d0*Logmu2/T1+
     .      0.5d0*(-b+mu2-T1)/(b*T1)-
     .      0.5d0*(deltT1bmu2*tmp242)/mu2)+
     .   tmp62*(1d0-3d0*Logb+Logb*Logmu2-
     .      (-1d0+Logb)*(-1d0+LogT1)+LogT1-Logmu2*LogT1+
     .      (phiT1bmu2*(b-mu2-T1))/mu2-
     .      (phiT1bmu2*(-b-mu2+T1))/mu2-
     .      0.5d0*(LogT1*(-b+mu2-T1))/b+
     .      0.5d0*(Logmu2*(b-mu2-T1))/T1+
     .      0.5d0*(Logb*(-b+mu2-T1))/T1-
     .      0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .      (-b-mu2+T1)*
     .       (phiT1bmu2/mu2-
     .       ((b-mu2-T1)*
     .          ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .            (mu2*tmp83)/(deltT1bmu2*T1)))/mu2-
     .       ((-b-mu2+T1)*
     .          ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .            (mu2*tmp90)/(b*deltT1bmu2)))/mu2+
     .       0.5d0*Logmu2/b-0.5d0*LogT1/b-0.5d0*Logb/T1+
     .       0.5d0*Logmu2/T1+0.5d0*(-b+mu2-T1)/(b*T1)-
     .       0.5d0*(deltT1bmu2*tmp242)/mu2)-
     .      0.5d0*(deltT1bmu2*
     .        ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .          (mu2*tmp83)/(deltT1bmu2*T1)))/mu2+
     .      0.5d0*(deltT1bmu2*
     .        ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .          (mu2*tmp90)/(b*deltT1bmu2)))/mu2)-
     .   0.5d0*(tmp214*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp56)/mb-0.5d0*(s2b*s2t)+
     .           0.5d0*(s2t*tmp128)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0-c2b)*(1d0+c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Xt)/mb)
     .      ))
      DT1b = DT1b-tmp239*
     .    (cbe*hb*ht*sbe*((mt*tmp56)/mb-0.5d0*(s2b*s2t)+
     .       0.5d0*(s2t*tmp154)/mb)+
     .      hb**2*sbe**2*(0.25d0*((1d0-c2b)*(1d0+c2t))-
     .       0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-0.25d0*(mt*s2b*s2t)/mb-
     .       0.25d0*((1d0-c2t)*s2b*Yt)/mb))-
     .   tmp238*(cbe*hb*ht*sbe*
     .       ((mt*tmp59)/mb+0.5d0*(s2b*s2t)+0.5d0*(s2t*tmp155)/mb)+
     .      hb**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+0.25d0*(mt*s2b*s2t)/mb+
     .       0.25d0*((1d0+c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .       0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Yt)/mb))-
     .     0.5d0*(tmp213*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp59)/mb+0.5d0*(s2b*s2t)+
     .           0.5d0*(s2t*tmp129)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0-c2b)*(1d0-c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Xt)/mb)
     .      ))
      DT2b = -((hb*ht*mu*s2t*tmp225)/mb)+
     .   tmp63*(1d0-3d0*Logb+Logb*Logmu2-
     .      (-1d0+Logb)*(-1d0+LogT2)+LogT2-Logmu2*LogT2+
     .      (phiT2bmu2*(b-mu2-T2))/mu2-
     .      (phiT2bmu2*(-b-mu2+T2))/mu2-
     .      0.5d0*(LogT2*(-b+mu2-T2))/b+
     .      0.5d0*(Logmu2*(b-mu2-T2))/T2+
     .      0.5d0*(Logb*(-b+mu2-T2))/T2-
     .      0.5d0*(Logmu2*(-b-mu2+T2))/b-
     .      0.5d0*(delT2Tbmu2*
     .        ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .          (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2+
     .      0.5d0*(delT2Tbmu2*
     .        ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .          (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      (-b-mu2+T2)*
     .       (phiT2bmu2/mu2-
     .       ((b-mu2-T2)*
     .          ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .            (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2-
     .       ((-b-mu2+T2)*
     .          ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .            (mu2*tmp110)/(b*delT2Tbmu2)))/mu2+
     .       0.5d0*Logmu2/b-0.5d0*LogT2/b-0.5d0*Logb/T2+
     .       0.5d0*Logmu2/T2+0.5d0*(-b+mu2-T2)/(b*T2)-
     .       0.5d0*(delT2Tbmu2*tmp227)/mu2))-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (-(phiT2bmu2/mu2)+
     .      ((b-mu2-T2)*
     .       ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .         (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2+
     .      ((-b-mu2+T2)*
     .       ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .         (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      0.5d0*Logmu2/b+0.5d0*LogT2/b+0.5d0*Logb/T2-
     .      0.5d0*Logmu2/T2-0.5d0*(-b+mu2-T2)/(b*T2)+
     .      0.5d0*(delT2Tbmu2*tmp227)/mu2)-
     .   0.5d0*(tmp218*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp59)/mb+0.5d0*(s2b*s2t)-
     .           0.5d0*(s2t*tmp128)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0-c2b)*(1d0-c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Xt)/mb)
     .      ))
      DT2b = DT2b-tmp223*
     .    (cbe*hb*ht*sbe*((mt*tmp59)/mb+0.5d0*(s2b*s2t)-
     .       0.5d0*(s2t*tmp154)/mb)+
     .      hb**2*sbe**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .       0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+0.25d0*(mt*s2b*s2t)/mb-
     .       0.25d0*((1d0+c2t)*s2b*Yt)/mb))-
     .   tmp222*(cbe*hb*ht*sbe*
     .       ((mt*tmp56)/mb-0.5d0*(s2b*s2t)-0.5d0*(s2t*tmp155)/mb)+
     .      hb**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-0.25d0*(mt*s2b*s2t)/mb+
     .       0.25d0*((1d0-c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*(0.25d0*((1d0-c2b)*(1d0+c2t))-
     .       0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Yt)/mb))-
     .     0.5d0*(tmp217*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp56)/mb-0.5d0*(s2b*s2t)-
     .           0.5d0*(s2t*tmp129)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0-c2b)*(1d0+c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Xt)/mb)
     .      ))
      DB1t = -(tmp196*(-0.5d0+2d0*LogB1-
     .      (phiA0B1T2*(-A0+B1-T2))/T2-0.5d0*(LogA0*LogB1)+
     .      0.5d0*(LogA0*LogT2)-0.5d0*(LogB1*LogT2)+
     .      0.5d0*(LogT2*(A0-B1-T2))/B1+
     .      0.5d0*(deltA0B1T2*phiA0B1T2)/(B1*T2)+
     .      0.5d0*(LogA0*(-A0-B1+T2))/B1-
     .      0.5d0*(deltA0B1T2*
     .          (tmp105/deltA0B1T2+
     .            (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1))-
     .   0.5d0*(tmp184*tmp211)-0.5d0*(tmp188*tmp215)-
     .   2d0*hb*ht*mt*mu*s2b*
     .    (phiB1tmu2/mu2-
     .      ((B1-mu2-t)*
     .       ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .         (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .      ((-B1-mu2+t)*
     .       ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .         (mu2*tmp71)/(B1*deltB1tmu2)))/mu2+
     .      0.5d0*Logmu2/B1-0.5d0*Logt/B1-0.5d0*LogB1/t+
     .      0.5d0*Logmu2/t+0.5d0*(-B1+mu2-t)/(B1*t)-
     .      0.5d0*(deltB1tmu2*tmp233)/mu2)-
     .   (hb*ht*mu*s2b*(-0.5d0+2d0*LogB1-
     .      (phiB1tmu2*(B1-mu2-t))/mu2-
     .      0.5d0*(LogB1*Logmu2)-0.5d0*(LogB1*Logt)+
     .      0.5d0*(Logmu2*Logt)+0.5d0*(Logt*(-B1+mu2-t))/B1+
     .      0.5d0*(Logmu2*(-B1-mu2+t))/B1-
     .      0.5d0*(deltB1tmu2*
     .          ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .            (mu2*tmp71)/(B1*deltB1tmu2)))/mu2))/mt+
     .   tmp61*(1d0+LogB1-LogB1*Logmu2-
     .      (-1d0+LogB1)*(-1d0+Logt)-3d0*Logt+Logmu2*Logt-
     .      (phiB1tmu2*(B1-mu2-t))/mu2+
     .      (phiB1tmu2*(-B1-mu2+t))/mu2+
     .      0.5d0*(Logt*(-B1+mu2-t))/B1-
     .      0.5d0*(Logmu2*(B1-mu2-t))/t-
     .      0.5d0*(LogB1*(-B1+mu2-t))/t+
     .      0.5d0*(Logmu2*(-B1-mu2+t))/B1-
     .      (B1-mu2-t)*
     .       (phiB1tmu2/mu2-
     .       ((B1-mu2-t)*
     .          ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .       ((-B1-mu2+t)*
     .          ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .            (mu2*tmp71)/(B1*deltB1tmu2)))/mu2+
     .       0.5d0*Logmu2/B1-0.5d0*Logt/B1-0.5d0*LogB1/t+
     .       0.5d0*Logmu2/t+0.5d0*(-B1+mu2-t)/(B1*t)-
     .       0.5d0*(deltB1tmu2*tmp233)/mu2)+
     .      0.5d0*(deltB1tmu2*
     .        ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .          (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .      0.5d0*(deltB1tmu2*
     .        ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .          (mu2*tmp71)/(B1*deltB1tmu2)))/mu2)
      DB1t = DB1t-tmp192*
     .    (-0.5d0+2d0*LogB1-(phiA0B1T1*(-A0+B1-T1))/T1-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB1*LogT1)+0.5d0*(LogT1*(A0-B1-T1))/B1+
     .      0.5d0*(deltA0B1T1*phiA0B1T1)/(B1*T1)+
     .      0.5d0*(LogA0*(-A0-B1+T1))/B1-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .          tmp85/deltA0B1T1))/B1)
      DB2t = -(tmp195*(-0.5d0+2d0*LogB2-
     .      (phiA0B2T2*(-A0+B2-T2))/T2-0.5d0*(LogA0*LogB2)+
     .      0.5d0*(LogA0*LogT2)-0.5d0*(LogB2*LogT2)+
     .      0.5d0*(LogT2*(A0-B2-T2))/B2+
     .      0.5d0*(deltA0B2T2*phiA0B2T2)/(B2*T2)+
     .      0.5d0*(LogA0*(-A0-B2+T2))/B2-
     .      0.5d0*(deltA0B2T2*
     .          (tmp106/deltA0B2T2+
     .            (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2))-
     .   0.5d0*(tmp183*tmp212)-0.5d0*(tmp187*tmp216)-
     .   2d0*hb*ht*mt*mu*s2b*
     .    (-(phiB2tmu2/mu2)+
     .      ((B2-mu2-t)*
     .       ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .         (mu2*tmp68)/(deltB2tmu2*t)))/mu2+
     .      ((-B2-mu2+t)*
     .       ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .         (mu2*tmp72)/(B2*deltB2tmu2)))/mu2-
     .      0.5d0*Logmu2/B2+0.5d0*Logt/B2+0.5d0*LogB2/t-
     .      0.5d0*Logmu2/t-0.5d0*(-B2+mu2-t)/(B2*t)+
     .      0.5d0*(deltB2tmu2*tmp237)/mu2)+
     .   tmp60*(1d0+LogB2-LogB2*Logmu2-
     .      (-1d0+LogB2)*(-1d0+Logt)-3d0*Logt+Logmu2*Logt-
     .      (phiB2tmu2*(B2-mu2-t))/mu2+
     .      (phiB2tmu2*(-B2-mu2+t))/mu2+
     .      0.5d0*(Logt*(-B2+mu2-t))/B2-
     .      0.5d0*(Logmu2*(B2-mu2-t))/t-
     .      0.5d0*(LogB2*(-B2+mu2-t))/t+
     .      0.5d0*(Logmu2*(-B2-mu2+t))/B2-
     .      (B2-mu2-t)*
     .       (phiB2tmu2/mu2-
     .       ((B2-mu2-t)*
     .          ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2-
     .       ((-B2-mu2+t)*
     .          ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .            (mu2*tmp72)/(B2*deltB2tmu2)))/mu2+
     .       0.5d0*Logmu2/B2-0.5d0*Logt/B2-0.5d0*LogB2/t+
     .       0.5d0*Logmu2/t+0.5d0*(-B2+mu2-t)/(B2*t)-
     .       0.5d0*(deltB2tmu2*tmp237)/mu2)+
     .      0.5d0*(deltB2tmu2*
     .        ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .          (mu2*tmp68)/(deltB2tmu2*t)))/mu2-
     .      0.5d0*(deltB2tmu2*
     .        ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .          (mu2*tmp72)/(B2*deltB2tmu2)))/mu2)-
     .   (hb*ht*mu*s2b*(0.5d0-2d0*LogB2+
     .      (phiB2tmu2*(B2-mu2-t))/mu2+
     .      0.5d0*(LogB2*Logmu2)+0.5d0*(LogB2*Logt)-
     .      0.5d0*(Logmu2*Logt)-0.5d0*(Logt*(-B2+mu2-t))/B2-
     .      0.5d0*(Logmu2*(-B2-mu2+t))/B2+
     .      0.5d0*(deltB2tmu2*
     .          ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .            (mu2*tmp72)/(B2*deltB2tmu2)))/mu2))/mt
      DB2t = DB2t-tmp191*
     .    (-0.5d0+2d0*LogB2-(phiA0B2T1*(-A0+B2-T1))/T1-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB2*LogT1)+0.5d0*(LogT1*(A0-B2-T1))/B2+
     .      0.5d0*(deltA0B2T1*phiA0B2T1)/(B2*T1)+
     .      0.5d0*(LogA0*(-A0-B2+T1))/B2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .          tmp86/deltA0B2T1))/B2)
      DT1B1 = ht**2*(0.25d0*((1d0+c2b)*(1d0-c2t))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT1)))+
     .     hb**2*(0.25d0*((1d0-c2b)*(1d0+c2t))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT1)))-
     .     0.5d0*(((-2d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**2)+
     .      (-2d0-2d0*LogB1)/T1+(2d0*LogT1)/T1-
     .      (2d0*(LogB1-LogT1))/((1d0-B1/T1)*T1)+
     .      (2d0*B1*(LogB1-LogT1)*(-B1+T1))/
     .       ((1d0-B1/T1)**2*T1**3)+
     .      (2d0*(-B1+T1))/((1d0-B1/T1)*T1**2)+
     .      (2d0*(LogB1-LogT1)*(-B1+T1))/((1d0-B1/T1)*T1**2))*
     .      tmp186)-tmp194*
     .    (phiA0B1T1/T1+(phiA0B1T1*(-A0-B1+T1))/(B1*T1)-
     .      ((-A0+B1-T1)*
     .       ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .         (B1*tmp81)/(deltA0B1T1*T1)))/B1-
     .      ((-A0-B1+T1)*
     .       ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .         tmp85/deltA0B1T1))/B1+0.5d0*LogA0/B1-
     .      0.5d0*LogT1/B1+0.5d0*LogA0/T1-0.5d0*LogB1/T1+
     .      0.5d0*(A0-B1-T1)/(B1*T1)+
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .          (B1*tmp81)/(deltA0B1T1*T1)))/B1**2-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1)/(deltA0B1T1*T1)-
     .          (2d0*B1*phiA0B1T1*(-A0+B1-T1)*(A0+B1-T1))/
     .           (deltA0B1T1**2*T1)+
     .          (B1*((A0-B1)/B1+LogA0-LogB1-T1/B1))/
     .           (deltA0B1T1*T1)+tmp81/(deltA0B1T1*T1)-
     .          (2d0*B1*(-A0+B1-T1)*tmp81)/
     .           (deltA0B1T1**2*T1)+
     .          ((A0+B1-T1)*
     .             ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .             tmp85/deltA0B1T1))/deltA0B1T1))/B1)
      DT2B1 = hb**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT2)))+
     .     ht**2*(0.25d0*((1d0+c2b)*(1d0+c2t))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT2)))-
     .     tmp198*(phiA0B1T2/T2+
     .      (phiA0B1T2*(-A0-B1+T2))/(B1*T2)-
     .      ((-A0+B1-T2)*
     .       ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .         (B1*tmp100)/(deltA0B1T2*T2)))/B1-
     .      ((-A0-B1+T2)*
     .       (tmp105/deltA0B1T2+
     .         (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1+
     .      0.5d0*LogA0/B1-0.5d0*LogT2/B1+0.5d0*LogA0/T2-
     .      0.5d0*LogB1/T2+0.5d0*(A0-B1-T2)/(B1*T2)+
     .      0.5d0*(deltA0B1T2*
     .        ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .          (B1*tmp100)/(deltA0B1T2*T2)))/B1**2-
     .      0.5d0*(deltA0B1T2*
     .        ((B1*phiA0B1T2)/(deltA0B1T2*T2)-
     .          (2d0*B1*phiA0B1T2*(-A0+B1-T2)*(A0+B1-T2))/
     .           (deltA0B1T2**2*T2)+
     .          (B1*((A0-B1)/B1+LogA0-LogB1-T2/B1))/
     .           (deltA0B1T2*T2)+tmp100/(deltA0B1T2*T2)-
     .          (2d0*B1*(-A0+B1-T2)*tmp100)/
     .           (deltA0B1T2**2*T2)+
     .          ((A0+B1-T2)*
     .             (tmp105/deltA0B1T2+
     .             (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/
     .           deltA0B1T2))/B1)-
     .   0.5d0*(((-2d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**2)+
     .      (-2d0-2d0*LogB1)/T2+(2d0*LogT2)/T2-
     .      (2d0*(LogB1-LogT2))/((1d0-B1/T2)*T2)+
     .      (2d0*B1*(LogB1-LogT2)*(-B1+T2))/
     .       ((1d0-B1/T2)**2*T2**3)+
     .      (2d0*(-B1+T2))/((1d0-B1/T2)*T2**2)+
     .      (2d0*(LogB1-LogT2)*(-B1+T2))/((1d0-B1/T2)*T2**2))*
     .      tmp190)
      DT1B2 = ht**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT1)))+
     .     hb**2*(0.25d0*((1d0+c2b)*(1d0+c2t))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1)))-
     .     0.5d0*(((-2d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**2)+
     .      (-2d0-2d0*LogB2)/T1+(2d0*LogT1)/T1-
     .      (2d0*(LogB2-LogT1))/((1d0-B2/T1)*T1)+
     .      (2d0*B2*(LogB2-LogT1)*(-B2+T1))/
     .       ((1d0-B2/T1)**2*T1**3)+
     .      (2d0*(-B2+T1))/((1d0-B2/T1)*T1**2)+
     .      (2d0*(LogB2-LogT1)*(-B2+T1))/((1d0-B2/T1)*T1**2))*
     .      tmp185)-tmp193*
     .    (phiA0B2T1/T1+(phiA0B2T1*(-A0-B2+T1))/(B2*T1)-
     .      ((-A0+B2-T1)*
     .       ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .         (B2*tmp82)/(deltA0B2T1*T1)))/B2-
     .      ((-A0-B2+T1)*
     .       ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .         tmp86/deltA0B2T1))/B2+0.5d0*LogA0/B2-
     .      0.5d0*LogT1/B2+0.5d0*LogA0/T1-0.5d0*LogB2/T1+
     .      0.5d0*(A0-B2-T1)/(B2*T1)+
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .          (B2*tmp82)/(deltA0B2T1*T1)))/B2**2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1)/(deltA0B2T1*T1)-
     .          (2d0*B2*phiA0B2T1*(-A0+B2-T1)*(A0+B2-T1))/
     .           (deltA0B2T1**2*T1)+
     .          (B2*((A0-B2)/B2+LogA0-LogB2-T1/B2))/
     .           (deltA0B2T1*T1)+tmp82/(deltA0B2T1*T1)-
     .          (2d0*B2*(-A0+B2-T1)*tmp82)/
     .           (deltA0B2T1**2*T1)+
     .          ((A0+B2-T1)*
     .             ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .             tmp86/deltA0B2T1))/deltA0B2T1))/B2)
      DT2B2 = hb**2*(0.25d0*((1d0+c2b)*(1d0-c2t))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2)))+
     .     ht**2*(0.25d0*((1d0-c2b)*(1d0+c2t))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT2)))-
     .     tmp197*(phiA0B2T2/T2+
     .      (phiA0B2T2*(-A0-B2+T2))/(B2*T2)-
     .      ((-A0+B2-T2)*
     .       ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .         (B2*tmp101)/(deltA0B2T2*T2)))/B2-
     .      ((-A0-B2+T2)*
     .       (tmp106/deltA0B2T2+
     .         (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2+
     .      0.5d0*LogA0/B2-0.5d0*LogT2/B2+0.5d0*LogA0/T2-
     .      0.5d0*LogB2/T2+0.5d0*(A0-B2-T2)/(B2*T2)+
     .      0.5d0*(deltA0B2T2*
     .        ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .          (B2*tmp101)/(deltA0B2T2*T2)))/B2**2-
     .      0.5d0*(deltA0B2T2*
     .        ((B2*phiA0B2T2)/(deltA0B2T2*T2)-
     .          (2d0*B2*phiA0B2T2*(-A0+B2-T2)*(A0+B2-T2))/
     .           (deltA0B2T2**2*T2)+
     .          (B2*((A0-B2)/B2+LogA0-LogB2-T2/B2))/
     .           (deltA0B2T2*T2)+tmp101/(deltA0B2T2*T2)-
     .          (2d0*B2*(-A0+B2-T2)*tmp101)/
     .           (deltA0B2T2**2*T2)+
     .          ((A0+B2-T2)*
     .             (tmp106/deltA0B2T2+
     .             (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/
     .           deltA0B2T2))/B2)-
     .   0.5d0*(((-2d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**2)+
     .      (-2d0-2d0*LogB2)/T2+(2d0*LogT2)/T2-
     .      (2d0*(LogB2-LogT2))/((1d0-B2/T2)*T2)+
     .      (2d0*B2*(LogB2-LogT2)*(-B2+T2))/
     .       ((1d0-B2/T2)**2*T2**3)+
     .      (2d0*(-B2+T2))/((1d0-B2/T2)*T2**2)+
     .      (2d0*(LogB2-LogT2)*(-B2+T2))/((1d0-B2/T2)*T2**2))*
     .      tmp189)
      Dbc2t = tmp7*(2d0*b*Logb+(-1d0+Logmu2)*mu2+
     .      (-1d0+Logb)*(-1d0+Logmu2)*mu2+2d0*Logmu2*mu2-
     .      (-1d0+LogT2)*T2-(-1d0+Logb)*(-1d0+LogT2)*T2+
     .      2d0*LogT2*T2-0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))+
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))+
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))-
     .      2.5d0*(b+mu2+T2)-
     .      (-b-mu2+T2)*
     .       (-0.5d0+2d0*Logb-(phiT2bmu2*(b-mu2-T2))/mu2-
     .       0.5d0*(Logb*Logmu2)-0.5d0*(Logb*LogT2)+
     .       0.5d0*(Logmu2*LogT2)+
     .       0.5d0*(LogT2*(-b+mu2-T2))/b+
     .       0.5d0*(Logmu2*(-b-mu2+T2))/b-
     .       0.5d0*(delT2Tbmu2*
     .           ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .             (mu2*tmp110)/(b*delT2Tbmu2)))/mu2))+
     .   0.5d0*(hb*ht*mu*tmp113)/(mb*s2t)+
     .   tmp6*(2d0*b*Logb+(-1d0+Logmu2)*mu2+
     .      (-1d0+Logb)*(-1d0+Logmu2)*mu2+2d0*Logmu2*mu2-
     .      (-1d0+LogT1)*T1-(-1d0+Logb)*(-1d0+LogT1)*T1+
     .      2d0*LogT1*T1-0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))-
     .      2.5d0*(b+mu2+T1)-
     .      (-b-mu2+T1)*
     .       (-0.5d0+2d0*Logb-(phiT1bmu2*(b-mu2-T1))/mu2-
     .       0.5d0*(Logb*Logmu2)-0.5d0*(Logb*LogT1)+
     .       0.5d0*(Logmu2*LogT1)+
     .       0.5d0*(LogT1*(-b+mu2-T1))/b+
     .       0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .       0.5d0*(deltT1bmu2*
     .           ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .             (mu2*tmp90)/(b*deltT1bmu2)))/mu2))+
     .   (hb*ht*mb*mu*(-((phiT1bmu2*(b-mu2-T1))/mu2)+
     .      (phiT2bmu2*(b-mu2-T2))/mu2-0.5d0*(Logb*LogT1)+
     .      0.5d0*(Logmu2*LogT1)+0.5d0*(Logb*LogT2)-
     .      0.5d0*(Logmu2*LogT2)+0.5d0*(LogT1*(-b+mu2-T1))/b+
     .      0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .      0.5d0*(LogT2*(-b+mu2-T2))/b-
     .      0.5d0*(Logmu2*(-b-mu2+T2))/b+
     .      0.5d0*(delT2Tbmu2*
     .          ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .            (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      0.5d0*(deltT1bmu2*
     .          ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .            (mu2*tmp90)/(b*deltT1bmu2)))/mu2))/s2t
      Dbc2t = Dbc2t-0.5d0*
     .    ((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp55)/mb-0.25d0*s2b/s2t-
     .           0.25d0*tmp129/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp55)/mb-0.25d0*s2b/s2t+
     .           0.25d0*tmp128/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp54)/mb+0.25d0*s2b/s2t+
     .           0.25d0*tmp129/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp54)/mb+0.25d0*s2b/s2t-
     .           0.25d0*tmp128/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xt)/(c2t*mb))))
      Dbc2t = Dbc2t-tmp88*
     .    (cbe*hb*ht*sbe*((mt*tmp55)/mb-0.25d0*s2b/s2t-
     .       0.25d0*tmp155/(mb*s2t))+
     .      hb**2*sbe**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yt)/(c2t*mb)))-
     .   tmp109*(cbe*hb*ht*sbe*
     .       ((mt*tmp55)/mb-0.25d0*s2b/s2t+0.25d0*tmp154/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)+
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mt*s2b)/(mb*s2t)-0.125d0*(s2b*Yt)/(c2t*mb)))-
     .     tmp108*(cbe*hb*ht*sbe*
     .       ((mt*tmp54)/mb+0.25d0*s2b/s2t+0.25d0*tmp155/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yt)/(c2t*mb)))-
     .     tmp89*(cbe*hb*ht*sbe*
     .       ((mt*tmp54)/mb+0.25d0*s2b/s2t-0.25d0*tmp154/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yt)/(c2t*mb)))
      DB1c2t = hb**2*(0.125d0*((1d0-c2b)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/c2t-
     .      0.125d0*((1d0-c2b)*(-1d0+LogT2)*T2)/c2t-
     .      0.125d0*((1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/c2t)+
     .   ht**2*(-(0.125d0*((1d0+c2b)*(-1d0+LogT1)*T1)/c2t)-
     .      0.125d0*((1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogT2)*T2)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/c2t)-
     .   tmp206*(-0.5d0+2d0*LogB1-(phiA0B1T2*(-A0+B1-T2))/T2-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT2)-
     .      0.5d0*(LogB1*LogT2)+0.5d0*(LogT2*(A0-B1-T2))/B1+
     .      0.5d0*(deltA0B1T2*phiA0B1T2)/(B1*T2)+
     .      0.5d0*(LogA0*(-A0-B1+T2))/B1-
     .      0.5d0*(deltA0B1T2*
     .        (tmp105/deltA0B1T2+
     .          (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1)-
     .   0.5d0*(tmp200*tmp211)-0.5d0*(tmp202*tmp215)-
     .   tmp204*(-0.5d0+2d0*LogB1-(phiA0B1T1*(-A0+B1-T1))/T1-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB1*LogT1)+0.5d0*(LogT1*(A0-B1-T1))/B1+
     .      0.5d0*(deltA0B1T1*phiA0B1T1)/(B1*T1)+
     .      0.5d0*(LogA0*(-A0-B1+T1))/B1-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .          tmp85/deltA0B1T1))/B1)
      DB2c2t = ht**2*(-(0.125d0*
     .       ((1d0-c2b)*(-1d0+LogT1)*T1)/c2t)-
     .      0.125d0*((1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogT2)*T2)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t)+
     .   hb**2*(0.125d0*((1d0+c2b)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t-
     .      0.125d0*((1d0+c2b)*(-1d0+LogT2)*T2)/c2t-
     .      0.125d0*((1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t)-
     .   tmp205*(-0.5d0+2d0*LogB2-(phiA0B2T2*(-A0+B2-T2))/T2-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT2)-
     .      0.5d0*(LogB2*LogT2)+0.5d0*(LogT2*(A0-B2-T2))/B2+
     .      0.5d0*(deltA0B2T2*phiA0B2T2)/(B2*T2)+
     .      0.5d0*(LogA0*(-A0-B2+T2))/B2-
     .      0.5d0*(deltA0B2T2*
     .        (tmp106/deltA0B2T2+
     .          (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2)-
     .   0.5d0*(tmp199*tmp212)-0.5d0*(tmp201*tmp216)-
     .   tmp203*(-0.5d0+2d0*LogB2-(phiA0B2T1*(-A0+B2-T1))/T1-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB2*LogT1)+0.5d0*(LogT1*(A0-B2-T1))/B2+
     .      0.5d0*(deltA0B2T1*phiA0B2T1)/(B2*T1)+
     .      0.5d0*(LogA0*(-A0-B2+T1))/B2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .          tmp86/deltA0B2T1))/B2)
      DT1c2b = ht**2*(0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1))/c2b-
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2))/c2b+
     .      0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT1))/c2b-
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT1))/c2b)+
     .   hb**2*(-(0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1))/c2b)+
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2))/c2b-
     .      0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT1))/c2b+
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1))/c2b)-
     .   0.5d0*(tmp214*(cbe**2*hb**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b+
     .         0.25d0*((1d0+c2t)*mb*Xb)/s2b+
     .         0.25d0*(mt*s2t*Xb)/c2b+0.125d0*((1d0+c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (mb*s2t*tmp23+2d0*mb*mt*tmp57+
     .         0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp130)/s2b+
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0+c2t)*t)/c2b+
     .         0.25d0*((1d0-c2t)*mb*Xt)/s2b-
     .         0.25d0*(mt*s2t*Xt)/c2b-0.125d0*((1d0-c2t)*Xt**2)/c2b
     .         )))-0.5d0*
     .    (tmp213*(cbe**2*hb**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0-c2t)*t)/c2b-
     .         0.25d0*((1d0+c2t)*mb*Xb)/s2b-
     .         0.25d0*(mt*s2t*Xb)/c2b-0.125d0*((1d0+c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (mb*s2t*tmp24+2d0*mb*mt*tmp58-
     .         0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp130)/s2b-
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b-
     .         0.25d0*((1d0-c2t)*mb*Xt)/s2b+
     .         0.25d0*(mt*s2t*Xt)/c2b+0.125d0*((1d0-c2t)*Xt**2)/c2b
     .         )))
      DT1c2b = DT1c2b-
     .   tmp239*(hb**2*sbe**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+0.25d0*(mb*mt*s2t)/s2b+
     .       0.125d0*((1d0-c2t)*t)/c2b+
     .       0.25d0*((1d0+c2t)*mb*Yb)/s2b+0.25d0*(mt*s2t*Yb)/c2b+
     .       0.125d0*((1d0+c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(mb*s2t*tmp39+2d0*mb*mt*tmp57+
     .       0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp156)/s2b+
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(0.125d0*(b*(1d0-c2t))/c2b+
     .       0.25d0*(mb*mt*s2t)/s2b-0.125d0*((1d0+c2t)*t)/c2b+
     .       0.25d0*((1d0-c2t)*mb*Yt)/s2b-0.25d0*(mt*s2t*Yt)/c2b-
     .       0.125d0*((1d0-c2t)*Yt**2)/c2b))-
     .   tmp238*(hb**2*sbe**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .       0.125d0*((1d0-c2t)*t)/c2b-
     .       0.25d0*((1d0+c2t)*mb*Yb)/s2b-0.25d0*(mt*s2t*Yb)/c2b-
     .       0.125d0*((1d0+c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(mb*s2t*tmp40+2d0*mb*mt*tmp58-
     .       0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp156)/s2b-
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(-(0.125d0*(b*(1d0-c2t))/c2b)-
     .       0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b-
     .       0.25d0*((1d0-c2t)*mb*Yt)/s2b+0.25d0*(mt*s2t*Yt)/c2b+
     .       0.125d0*((1d0-c2t)*Yt**2)/c2b))
      DT2c2b = hb**2*(-(0.125d0*
     .       (B1*(1d0-c2t)*(-1d0+LogB1))/c2b)+
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2))/c2b-
     .      0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT2))/c2b+
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2))/c2b)+
     .   ht**2*(0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1))/c2b-
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2))/c2b+
     .      0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT2))/c2b-
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT2))/c2b)-
     .   0.5d0*(tmp218*(cbe**2*hb**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b+
     .         0.25d0*((1d0-c2t)*mb*Xb)/s2b-
     .         0.25d0*(mt*s2t*Xb)/c2b+0.125d0*((1d0-c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (-(mb*s2t*tmp23)+2d0*mb*mt*tmp58-
     .         0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp131)/s2b-
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0-c2t)*t)/c2b+
     .         0.25d0*((1d0+c2t)*mb*Xt)/s2b+
     .         0.25d0*(mt*s2t*Xt)/c2b-0.125d0*((1d0+c2t)*Xt**2)/c2b
     .         )))-0.5d0*
     .    (tmp217*(cbe**2*hb**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0+c2t)*t)/c2b-
     .         0.25d0*((1d0-c2t)*mb*Xb)/s2b+
     .         0.25d0*(mt*s2t*Xb)/c2b-0.125d0*((1d0-c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (-(mb*s2t*tmp24)+2d0*mb*mt*tmp57+
     .         0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp131)/s2b+
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b-
     .         0.25d0*((1d0+c2t)*mb*Xt)/s2b-
     .         0.25d0*(mt*s2t*Xt)/c2b+0.125d0*((1d0+c2t)*Xt**2)/c2b
     .         )))
      DT2c2b = DT2c2b-
     .   tmp223*(hb**2*sbe**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-0.25d0*(mb*mt*s2t)/s2b+
     .       0.125d0*((1d0+c2t)*t)/c2b+
     .       0.25d0*((1d0-c2t)*mb*Yb)/s2b-0.25d0*(mt*s2t*Yb)/c2b+
     .       0.125d0*((1d0-c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(-(mb*s2t*tmp39)+2d0*mb*mt*tmp58-
     .       0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp157)/s2b-
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(0.125d0*(b*(1d0+c2t))/c2b-
     .       0.25d0*(mb*mt*s2t)/s2b-0.125d0*((1d0-c2t)*t)/c2b+
     .       0.25d0*((1d0+c2t)*mb*Yt)/s2b+0.25d0*(mt*s2t*Yt)/c2b-
     .       0.125d0*((1d0+c2t)*Yt**2)/c2b))-
     .   tmp222*(hb**2*sbe**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .       0.125d0*((1d0+c2t)*t)/c2b-
     .       0.25d0*((1d0-c2t)*mb*Yb)/s2b+0.25d0*(mt*s2t*Yb)/c2b-
     .       0.125d0*((1d0-c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(-(mb*s2t*tmp40)+2d0*mb*mt*tmp57+
     .       0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp157)/s2b+
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(-(0.125d0*(b*(1d0+c2t))/c2b)+
     .       0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b-
     .       0.25d0*((1d0+c2t)*mb*Yt)/s2b-0.25d0*(mt*s2t*Yt)/c2b+
     .       0.125d0*((1d0+c2t)*Yt**2)/c2b))
      Dc2tc2b = hb**2*tmp19+ht**2*tmp19-
     .   tmp89*(hb**2*sbe**2*tmp38+cbe**2*ht**2*tmp49+
     .      cbe*hb*ht*sbe*(-(0.25d0*(mb*mt)/(c2b*c2t))-
     .       0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp39)/s2t+
     .       0.5d0*(mt*tmp43)/s2b-0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   tmp108*(hb**2*sbe**2*tmp38+cbe**2*ht**2*tmp49+
     .      cbe*hb*ht*sbe*(-(0.25d0*(mb*mt)/(c2b*c2t))-
     .       0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp40)/s2t-
     .       0.5d0*(mt*tmp44)/s2b-0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (cbe**2*hb**2*tmp22+ht**2*sbe**2*tmp33-
     .      cbe*hb*ht*sbe*
     .       (-(0.25d0*(mb*mt)/(c2b*c2t))-
     .         0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp23)/s2t+
     .         0.5d0*(mt*tmp27)/s2b-0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (cbe**2*hb**2*tmp22+ht**2*sbe**2*tmp33-
     .      cbe*hb*ht*sbe*
     .       (-(0.25d0*(mb*mt)/(c2b*c2t))-
     .         0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp24)/s2t-
     .         0.5d0*(mt*tmp28)/s2b-0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (cbe**2*hb**2*tmp21+ht**2*sbe**2*tmp32-
     .      cbe*hb*ht*sbe*
     .       (0.25d0*(mb*mt)/(c2b*c2t)+0.125d0*(b+t)/(s2b*s2t)-
     .         0.5d0*(mb*tmp24)/s2t-0.5d0*(mt*tmp27)/s2b+
     .         0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (cbe**2*hb**2*tmp21+ht**2*sbe**2*tmp32-
     .      cbe*hb*ht*sbe*
     .       (0.25d0*(mb*mt)/(c2b*c2t)+0.125d0*(b+t)/(s2b*s2t)+
     .         0.5d0*(mb*tmp23)/s2t+0.5d0*(mt*tmp28)/s2b+
     .         0.125d0*(Xb*Xt)/(s2b*s2t))))
      Dc2tc2b = Dc2tc2b-
     .   tmp88*(hb**2*sbe**2*tmp37+cbe**2*ht**2*tmp48+
     .      cbe*hb*ht*sbe*(0.25d0*(mb*mt)/(c2b*c2t)+
     .       0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp40)/s2t-
     .       0.5d0*(mt*tmp43)/s2b+0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   tmp109*(hb**2*sbe**2*tmp37+cbe**2*ht**2*tmp48+
     .      cbe*hb*ht*sbe*(0.25d0*(mb*mt)/(c2b*c2t)+
     .       0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp39)/s2t+
     .       0.5d0*(mt*tmp44)/s2b+0.125d0*(Yb*Yt)/(s2b*s2t)))
      Dcptpb = -2d0*cbe*hb*ht*mb*mt*sbe*tmp108*tmp56+
     .   cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*tmp56
     .    +cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*tmp56
     .    -2d0*cbe*hb*ht*mb*mt*sbe*tmp109*tmp59+
     .   cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*tmp59
     .    +cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*tmp59
     .    -2d0*cbe*hb*ht*mb*mt*sbe*tmp59*tmp88-
     .   2d0*cbe*hb*ht*mb*mt*sbe*tmp56*tmp89-
     .   4d0*cbe*hb*ht*mb*mt*sbe*
     .    (-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76)
      Dcpttptb = 0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*Xb*
     .      Xt)-0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*Xb*
     .      Xt)-0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*Xb*
     .      Xt)+0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*Xb*
     .      Xt)+0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp108*Yb*Yt)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp109*Yb*Yt)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp88*Yb*Yt)+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp89*Yb*Yt)
      Dcpbptt = -2d0*hb*ht*mb*mu*s2t*tmp113-
     .   cbe*hb*ht*sbe*tmp89*(mb*s2t*tmp154-0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp108*
     .    (-(mb*s2t*tmp155)-0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp109*
     .    (-(mb*s2t*tmp154)+0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp88*(mb*s2t*tmp155+0.5d0*(b*s2b*s2t))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (mb*s2t*tmp128-0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(mb*s2t*tmp129)-0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(mb*s2t*tmp128)+0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (mb*s2t*tmp129+0.5d0*(b*s2b*s2t)))
      Dcptptb = -2d0*hb*ht*mt*mu*s2b*tmp75-
     .   cbe*hb*ht*sbe*tmp89*
     .    (-(mt*s2b*tmp156)-0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp108*(mt*s2b*tmp157-0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp88*(mt*s2b*tmp156+0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp109*
     .    (-(mt*s2b*tmp157)+0.5d0*(s2b*s2t*t))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(mt*s2b*tmp130)-0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (mt*s2b*tmp131-0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (mt*s2b*tmp130+0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(mt*s2b*tmp131)+0.5d0*(s2b*s2t*t)))
      Dcptmptt = -(tmp109*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      hb**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2b)*mt*s2t*Yt))))-
     .     tmp89*(-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0+c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0-c2b)*mt*s2t*Yt)))-
     .   tmp108*(-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))-
     .       0.5d0*((1d0-c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))-
     .       0.5d0*((1d0+c2b)*mt*s2t*Yt)))-
     .   0.25d0*(ht**2*(-2d0*mt*s2t*sbe**2*tmp115*Xt+
     .      2d0*mt*s2t*sbe**2*tmp94*Xt-
     .      4d0*cbe**2*mt*s2t*tmp114*Yt+4d0*cbe**2*mt*s2t*tmp93*Yt)
     .      )-0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-
     .      5d0*T2+LogT2*(-2d0*B2*LogB2+4d0*T2)-
     .      2d0*(-B2+T2)*tmp173)*
     .      (-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2b)*mt*s2t*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2b)*mt*s2t*Xt))))-
     .     0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0+c2b)*mt*s2t*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0-c2b)*mt*s2t*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-
     .      5d0*T2+LogT2*(-2d0*B1*LogB1+4d0*T2)-
     .      2d0*(-B1+T2)*tmp172)*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2b)*mt*s2t*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2b)*mt*s2t*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-
     .      5d0*T1+LogT1*(-2d0*B1*LogB1+4d0*T1)-
     .      2d0*(-B1+T1)*tmp169)*
     .      (-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0-c2b)*mt*s2t*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0+c2b)*mt*s2t*Xt))))
      Dcptmptt = Dcptmptt-
     .   tmp88*(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0-c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0+c2b)*mt*s2t*Yt)))
      Dcpbmptb = -(tmp89*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      hb**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2t)*mb*s2b*Yb))
     .       +cbe**2*ht**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2t)*mb*s2b*Yt))
     .      ))-0.25d0*(hb**2*
     .      (2d0*cbe**2*(-10d0*B1+4d0*B1*LogB1+
     .         LogB1*(4d0*B1-2d0*B1*LogB1))*mb*s2b*Xb-
     .      2d0*cbe**2*(-10d0*B2+4d0*B2*LogB2+
     .         LogB2*(4d0*B2-2d0*B2*LogB2))*mb*s2b*Xb+
     .      4d0*mb*s2b*sbe**2*Yb*
     .       (2d0*A0*LogA0+4d0*B1*LogB1-A0*LogA0*LogB1-
     .         2.5d0*(A0+2d0*B1)+0.5d0*((A0-2d0*B1)*LogB1**2)-
     .         0.5d0*(deltA0B1B1*phiA0B1B1)/B1)-
     .      4d0*mb*s2b*sbe**2*Yb*
     .       (2d0*A0*LogA0+4d0*B2*LogB2-A0*LogA0*LogB2-
     .         2.5d0*(A0+2d0*B2)+0.5d0*((A0-2d0*B2)*LogB2**2)-
     .         0.5d0*(deltA0B2B2*phiA0B2B2)/B2)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2t)*mb*s2b*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2t)*mb*s2b*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-
     .      5d0*T1+LogT1*(-2d0*B1*LogB1+4d0*T1)-
     .      2d0*(-B1+T1)*tmp169)*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0+c2t)*mb*s2b*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0-c2t)*mb*s2b*Xt))))-
     .     0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2t)*mb*s2b*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2t)*mb*s2b*Xt))))-
     .     0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0-c2t)*mb*s2b*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0+c2t)*mb*s2b*Xt))
     .      ))
      Dcpbmptb = Dcpbmptb-
     .   tmp88*(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0+c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0-c2t)*mb*s2b*Yt)))-
     .   tmp109*(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)-
     .       0.5d0*((1d0-c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)-
     .       0.5d0*((1d0+c2t)*mb*s2b*Yt)))-
     .   tmp108*(-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0-c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0+c2t)*mb*s2b*Yt)))
      Dspbmptbspbptt =
     .  -(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp108))+
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp109)-
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))+
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170))+
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172))-
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173))+
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp88)-
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp89)
      Dsptmpttsptptb =
     .  -(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp108))+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp109)-
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))+
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170))+
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172))-
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173))+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp88)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp89)
      Dsptmpttspbmptb =
     .  -(tmp108*tmp11)-tmp109*tmp12-tmp12*tmp88-
     .   tmp11*tmp89-0.5d0*
     .    (tmp14*(-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))-
     .   0.5d0*(tmp13*(-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-
     .      5d0*T1+LogT1*(-2d0*B2*LogB2+4d0*T1)-
     .      2d0*(-B2+T1)*tmp170))-
     .   0.5d0*(tmp13*(-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-
     .      5d0*T2+LogT2*(-2d0*B1*LogB1+4d0*T2)-
     .      2d0*(-B1+T2)*tmp172))-
     .   0.5d0*(tmp14*(-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-
     .      5d0*T2+LogT2*(-2d0*B2*LogB2+4d0*T2)-
     .      2d0*(-B2+T2)*tmp173))

      end
